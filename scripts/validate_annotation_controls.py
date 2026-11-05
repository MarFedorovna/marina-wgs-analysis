#!/usr/bin/env python3
"""Validate offline annotation using public ClinVar variants and synthetic genotypes."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import uuid

from workflow_common import Runtime, artifact_manifest, canonical_hash, promote_exclusive, sha256
import annotation_stages
import alphamissense_adapter

TERMS = ['stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant',
         'splice_acceptor_variant', 'splice_donor_variant', 'missense_variant', 'synonymous_variant']
CLINVAR_MC = {'stop_gained': ('SO:0001587', 'nonsense'),
    'stop_lost': ('SO:0001578', 'stop_lost'), 'start_lost': ('SO:0002012', 'SO:0001582', 'start_lost', 'initiator_codon_variant'),
    'frameshift_variant': ('SO:0001589',), 'splice_acceptor_variant': ('SO:0001574',),
    'splice_donor_variant': ('SO:0001575',), 'missense_variant': ('SO:0001583',),
    'synonymous_variant': ('SO:0001819',)}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def select_public_controls(path):
    categories, genes, fillers = {}, {'MYH7': [], 'MYBPC3': []}, []
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            chrom, pos, vid, ref, alt = fields[:5]
            if (chrom not in {str(i) for i in range(1, 23)} or not re.fullmatch('[ACGT]+', ref + alt)
                or max(len(ref), len(alt)) > 40 or ',' in alt or vid == '.'):
                continue
            info = dict((entry.split('=', 1) + [''])[:2] for entry in fields[7].split(';'))
            record = dict(chrom='chr' + chrom, pos=int(pos), id=vid, ref=ref, alt=alt,
                          source_gene=info.get('GENEINFO', ''), source_MC=info.get('MC', ''))
            if len(fillers) < 40:
                fillers.append(record)
            for term in TERMS:
                if any(alias in info.get('MC', '') for alias in CLINVAR_MC[term]):
                    categories.setdefault(term, record)
            if re.search(r'(^|[|/])(?:Pathogenic|Likely_pathogenic)(?=$|[|/])', info.get('CLNSIG', ''), re.I):
                categories.setdefault('external_path_assertion', record)
            for gene in genes:
                if re.search(r'(^|[|])' + gene + ':', info.get('GENEINFO', '')) and len(genes[gene]) < 5:
                    genes[gene].append(record)
            if len(categories) == len(TERMS) + 1 and all(len(value) == 5 for value in genes.values()):
                break
    need(all(term in categories for term in TERMS) and all(genes.values()),
         'Public ClinVar control classes/target genes not found; observed categories=' + ','.join(sorted(categories)))
    selected = {}
    for record in [*categories.values(), *genes['MYH7'], *genes['MYBPC3'], *fillers]:
        key = record['chrom'], record['pos'], record['ref'], record['alt']
        selected.setdefault(key, record)
        if len(selected) == 20:
            break
    need(len(selected) == 20, 'Need twenty distinct public control alleles.')
    return sorted(selected.values(), key=lambda r: (int(r['chrom'][3:]), r['pos'], r['ref'], r['alt']))


def write_public_vcf(path, records, contigs):
    with path.open('x') as stream:
        stream.write('##fileformat=VCFv4.2\n##source=PUBLIC_ClinVar_with_SYNTHETIC_control_genotypes_not_patients\n')
        for chrom in sorted({record['chrom'] for record in records}, key=lambda name: int(name[3:])):
            need(chrom in contigs, 'Public control chromosome absent from the analysis reference.')
            stream.write('##contig=<ID=' + chrom + ',length=' + str(contigs[chrom]) + '>\n')
        stream.write('##INFO=<ID=CONTROL_ID,Number=1,Type=String,Description="Original public ClinVar variation ID">\n')
        for field, number, kind in [('GT', '1', 'String'), ('FT', '1', 'String'), ('DP', '1', 'Integer'), ('AD', 'R', 'Integer'), ('GQ', '1', 'Integer')]:
            stream.write(f'##FORMAT=<ID={field},Number={number},Type={kind},Description="Synthetic public-control field">\n')
        stream.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPUBLIC_ALT\tPUBLIC_REF\tPUBLIC_UNKNOWN\tPUBLIC_FILTERED\n')
        for record in records:
            stream.write('\t'.join(map(str, [record['chrom'], record['pos'], record['id'], record['ref'], record['alt'],
                100, 'PASS', 'CONTROL_ID=' + record['id'], 'GT:FT:DP:AD:GQ',
                '0/1:PASS:20:10,10:99', '0/0:PASS:20:20,0:99', './.:NoCall:0:0,0:0', '0/1:LowGQ:2:1,1:1'])) + '\n')


def check_existing(target, identity):
    receipt = json.loads((target / 'public-control-receipt.json').read_text())
    need(receipt.get('status') == 'PASS' and receipt.get('identity') == identity, 'Control receipt belongs to another source/config/code.')
    for name, item in receipt['artifacts'].items():
        path = target / name
        need(path.resolve(strict=True).is_relative_to(target.resolve()) and path.stat().st_size == item['bytes']
             and sha256(path) == item['sha256'], 'A completed public-control artifact changed.')
    return target / 'public-control-receipt.json'


def validate(cfg, output_root, execute=False, package_sha=None):
    resources_path = Path(cfg['annotation']['resource_manifest'])
    alpha_path = alphamissense_adapter.manifest_location(cfg)
    identity = canonical_hash(dict(reference_sha256=cfg['reference']['sha256'],
        resource_manifest_sha256=sha256(resources_path), alpha_manifest_sha256=sha256(alpha_path),
        code=package_sha or {name: sha256(Path(__file__).parent / name) for name in
            ['annotation_stages.py', 'alphamissense_adapter.py', 'validate_annotation_controls.py', 'workflow_common.py']},
        tools=cfg['tools']))
    target = output_root / ('controls-' + identity[:20])
    if not execute:
        print(json.dumps({'status': 'PLAN', 'target': str(target), 'public_only': True}))
        return target / 'public-control-receipt.json'
    if target.exists():
        return check_existing(target, identity)
    output_root.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + '.partial-' + uuid.uuid4().hex[:12])
    partial.mkdir(mode=0o700)
    resources = json.loads(resources_path.read_text())
    need(resources.get('status') == 'PASS', 'Public resources must be ready before controls.')
    clinvar = annotation_stages.resource_path(resources_path, resources['resources']['clinvar'])
    records = select_public_controls(clinvar)
    with Path(cfg['reference']['fasta'] + '.fai').open() as stream:
        contigs = {row[0]: int(row[1]) for row in csv.reader(stream, delimiter='\t')}
    write_public_vcf(partial / 'public-controls.vcf', records, contigs)
    (partial / 'selection.json').write_text(json.dumps({'source': 'ClinVar public snapshot', 'patient_data': False, 'records': records}, indent=2))
    annotation = dict(resource_manifest=str(resources_path), alphamissense_manifest=str(alpha_path),
        assembly='GRCh38', cache_version=115, input_vcf=str(partial / 'public-controls.normalized.vcf'),
        forks=2, memory_gib=8, make_plots=False, export_docx=False, pathway_mode='disabled')
    control_cfg = dict(schema_version=2, data_root=cfg['data_root'], output_root=str(partial / 'workflow'),
        reference=cfg['reference'], tools=cfg['tools'], annotation=annotation,
        budget=cfg['budget'], resources={'max_hours': 1, 'min_free_gib': 100},
        synthetic_validation=True, analysis_scope={'public_controls_only': True, 'clinical_diagnosis': False})
    runtime = Runtime(control_cfg, execute=True, workdir=partial, stage='annotate')
    runtime.run([cfg['tools'].get('bcftools', 'bcftools'), 'norm', '-f', cfg['reference']['fasta'], '-m', '-any',
        '-c', 'e', '-Ov', '-o', str(partial / 'public-controls.normalized.vcf'), str(partial / 'public-controls.vcf')], logname='normalize-public-controls')
    ann = runtime.stage_path('annotate')
    ann.mkdir(parents=True)
    annotation_stages.run_stage('annotate', control_cfg, ann, runtime.child(ann, 'annotate'))
    ids, observed_genes, sift, polyphen = {}, set(), 0, 0
    with gzip.open(ann / 'transcript_samples.tsv.gz', 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            ids.setdefault(row['variant_id'], set()).update(re.findall(r'\d+', row['clinvar_ids']))
            observed_genes.add(row['gene'])
            sift += int(bool(row['SIFT']))
            polyphen += int(bool(row['PolyPhen']))
    expected = {}
    with (partial / 'public-controls.normalized.vcf').open() as stream:
        for line in stream:
            if not line.startswith('#'):
                f = line.rstrip().split('\t')
                expected[annotation_stages.variant_key(f[0], f[1], f[3], f[4])] = f[2]
    need(len(expected) == 20 and all(value in ids.get(key, set()) for key, value in expected.items()),
         'Not all twenty normalized public alleles joined their exact ClinVar variation IDs; patient annotation blocked.')
    need({'MYH7', 'MYBPC3'}.issubset(observed_genes), 'Public HCM target-gene consequence controls missing.')
    need(sift > 0 and polyphen > 0, 'SIFT/PolyPhen cache fields did not populate the public missense controls.')
    counters = json.loads((ann / 'consequence_counts.json').read_text())
    need(counters['samples']['PUBLIC_ALT']['qc']['unique_normalized_alleles'] == 20
         and counters['samples']['PUBLIC_REF']['raw']['unique_normalized_alleles'] == 0
         and counters['samples']['PUBLIC_FILTERED']['qc']['unique_normalized_alleles'] == 0,
         'Public control genotype/quality/unique-allele counts failed.')
    with gzip.open(ann / 'intersections.tsv.gz', 'rt') as stream:
        matrix = list(csv.DictReader(stream, delimiter='\t'))
    need(all(row['verified_unique_to_one_sample'] == '0' and 'PUBLIC_UNKNOWN' in row['unknown_or_filtered_samples']
             and row['verified_reference_samples'] == 'PUBLIC_REF' for row in matrix), 'Unknown/no-call was incorrectly treated as reference.')
    report_dir = runtime.stage_path('report')
    report_dir.mkdir(parents=True)
    annotation_stages.run_stage('report', control_cfg, report_dir, runtime.child(report_dir, 'report'))
    for name, limit in [('top10_PUBLIC_ALT.tsv', 10), ('top20_PUBLIC_ALT.tsv', 20), ('cohort_top20.tsv', 20)]:
        with (report_dir / name).open() as stream:
            need(len(list(csv.DictReader(stream, delimiter='\t'))) <= limit, 'Public ranking test exceeded its requested size.')
    receipt = dict(schema_version=1, status='PASS', identity=identity, public_variants=20,
        ClinVar_exact_ID_matches=20, target_genes=['MYH7', 'MYBPC3'], SIFT_fields=sift, PolyPhen_fields=polyphen,
        consequence_flags=counters['cohort']['raw']['nonexclusive_flags'], unknown_not_reference=True,
        patient_data_used=False, genotypes='synthetic controls', VEP_execution=True,
        resource_manifest_sha256=sha256(resources_path), alpha_manifest_sha256=sha256(alpha_path),
        reference_sha256=cfg['reference']['sha256'], package_sha256=package_sha,
        artifacts=artifact_manifest(partial), clinical_diagnosis=False)
    (partial / 'public-control-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    promote_exclusive(partial, target)
    print('PUBLIC_ANNOTATION_CONTROLS_PASS ' + str(target / 'public-control-receipt.json'), flush=True)
    return target / 'public-control-receipt.json'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, default=Path('/srv/marina/resources/annotation-validation'))
    parser.add_argument('--package-sha')
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    validate(json.loads(args.config.read_text()), args.output_root, args.execute, args.package_sha)
