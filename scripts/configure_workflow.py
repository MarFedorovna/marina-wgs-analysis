#!/usr/bin/env python3
"""Build a private downstream plan from existing verified alignments and offline resources."""
import argparse
import csv
import io
import json
import os
from pathlib import Path
import re
import tempfile

from workflow_common import canonical_hash, require, sha256

GATK_IMAGE = 'broadinstitute/gatk@sha256:e7996ba655225c1cde0a1faec6a113e217758310af2cf99b00d61dae8ec6e9f2'
VEP_IMAGE = 'ensemblorg/ensembl-vep@sha256:ed6b660e278458109afc00d5e730aee240a3bf77e817796386033e1a18963152'


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def tsv_read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def tsv_bytes(rows, fields):
    buffer = io.StringIO(newline='')
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    return buffer.getvalue().encode()


def exclusive(path, data):
    path = Path(path)
    require(not path.is_symlink(), 'Generated configuration cannot be a symlink.')
    if path.exists():
        require(path.read_bytes() == data, 'Existing workflow plan differs; use another --label rather than overwrite.')
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.workflow-config-', dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(data); stream.flush(); os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def configure(root, resource_manifest, label='v1', execute=False, callsets_manifest=None, allow_pending_resources=False):
    root = Path(root).resolve(strict=True)
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,60}', label), 'Invalid workflow label.')
    private = root / 'private'
    pilot = load(private / 'pilot-run.json')
    remainder = load(private / 'remaining-align.json')
    samples_file = private / 'samples.tsv'
    samples = tsv_read(samples_file)
    ids = [row['sample_id'] for row in samples]
    require(ids and len(set(ids)) == len(ids), 'Invalid canonical cohort.')
    selected = tsv_read(pilot['samples_manifest'])
    require(len(selected) == 1 and selected[0] in samples, 'Pilot does not match the canonical cohort.')
    first = selected[0]['sample_id']
    reference_root = Path(pilot['reference_root']).resolve(strict=True)
    reference = (reference_root / pilot['reference']['fasta']).resolve(strict=True)
    require(sha256(reference) == pilot['reference']['sha256'], 'Verified reference changed.')
    require(remainder['image_digest'] == pilot['image_digest']
            and remainder['reference']['sha256'] == pilot['reference']['sha256'], 'Alignment runs use different software/reference.')
    alignments = []
    for sample in samples:
        sid = sample['sample_id']
        require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}', sid), 'Invalid canonical identifier.')
        directory = Path(pilot['output_root'] if sid == first else remainder['output_root']) / sid / 'align'
        receipt_path = directory / 'receipt.json'
        require(receipt_path.is_file(), 'An alignment is not completed yet; finish stage 2 before configuring downstream.')
        receipt = load(receipt_path)
        require(receipt.get('sample') == sid and receipt['config']['image_digest'] == pilot['image_digest']
                and receipt['config']['reference']['sha256'] == pilot['reference']['sha256'], 'Alignment provenance mismatch.')
        for name, expected in receipt['outputs'].items():
            path = (directory / name).resolve(strict=True)
            require(path.is_relative_to(directory.resolve()) and sha256(path) == expected,
                    'An alignment output differs from its completed receipt.')
        alignments.append({'sample_id': sid, 'bam': str(directory / 'marked.bam'), 'receipt': str(receipt_path)})
    resource_manifest = Path(resource_manifest).resolve()
    if resource_manifest.is_file():
        resources = load(resource_manifest)
        require(resources.get('status') == 'PASS' and resources.get('assembly') == 'GRCh38', 'Offline resource preparation is not complete.')
        require(resources['vep']['image'] == VEP_IMAGE, 'Unexpected VEP image; use a separately reviewed plan.')
    else:
        require(allow_pending_resources, 'Resources are not ready. --allow-pending-resources permits calling/filtering only; annotation still requires verified resources.')
    final_private = private / ('workflow-' + label)
    alignment_path = final_private / 'alignments.tsv'
    interval_path = reference_root / 'analysis-plans' / ('autosomes-' + label + '.interval_list')
    lengths = dict((row[0], int(row[1])) for row in (line.split('\t') for line in Path(str(reference) + '.fai').read_text().splitlines()))
    headers = [line for line in reference.with_suffix('.dict').read_text().splitlines() if line.startswith('@')]
    contigs = ['chr' + str(i) for i in range(1, 23)]
    require(all(contig in lengths for contig in contigs), 'The selected reference lacks expected autosomes.')
    intervals = ('\n'.join(headers) + '\n' + ''.join(
        f'{contig}\t1\t{lengths[contig]}\t+\t{contig}_diploid\n' for contig in contigs)).encode()
    callsets = [{'name': 'autosomes', 'intervals': str(interval_path),
                 'intervals_sha256': __import__('hashlib').sha256(intervals).hexdigest(), 'ploidy': 2,
                 'samples': ids, 'scope': 'autosomes', 'excluded_contigs': ['chrX', 'chrY', 'chrM', 'ALT/decoy']}]
    if callsets_manifest:
        callsets = load(callsets_manifest)
        require(isinstance(callsets, list) and callsets, 'Explicit callsets must be a nonempty JSON list.')
    scope = {'label': 'Approved germline SNV/indel intervals',
             'default': 'chr1-chr22 diploid autosomes' if not callsets_manifest else 'explicit supplied callsets',
             'not_assessed': ['CNV', 'SV', 'repeat expansions', 'mtDNA heteroplasmy', 'clinical diagnosis'],
             'sex_chromosomes': 'not_assessed_without_explicit_interval_and_per_sample_ploidy_policy',
             'coverage': 'Strongly nonuniform coverage requires per-region quality; absence of a call is not absence of variation.'}
    cfg = {'schema_version': 2, 'data_root': str(root), 'samples_manifest': str(samples_file),
           'alignment_manifest': str(alignment_path), 'output_root': str(root / 'outputs' / ('workflow-' + label)),
           'scratch_root': str(root / 'scratch' / ('workflow-' + label)),
           'reference': {'fasta': str(reference), 'sha256': pilot['reference']['sha256'],
                         'assembly_label': pilot['reference']['assembly_label']},
           'tools': {'parabricks_image': pilot['image_digest'], 'gatk_image': GATK_IMAGE,
                     'vep_image': VEP_IMAGE, 'bcftools': 'bcftools'},
           'callsets': callsets, 'scope': scope, 'analysis_scope': scope,
           'resources': {'max_hours': 6, 'min_free_gib': 100, 'gpu_device': 0,
                         'call_cpus': 24, 'call_memory_gib': 160, 'joint_cpus': 8,
                         'joint_memory_gib': 48, 'joint_java_heap_gib': 32},
           'budget': {'hourly_usd': 2.39, 'planned_vm_usd': 40, 'user_budget_usd': 50},
           'filters': {'min_dp': 10, 'min_gq': 20},
           'annotation': {'resource_manifest': str(resource_manifest), 'assembly': 'GRCh38', 'cache_version': 115,
                          'forks': 4, 'memory_gib': 24, 'fourth_predictor': 'CADD_unconfirmed',
                          'disease_context': 'hypertrophic cardiomyopathy research cohort',
                          'pathway_mode': 'disabled', 'top_sizes': [10, 20]}}
    server_state = private / 'server-budget.json'
    if server_state.is_file():
        cfg['budget'].update(load(server_state))
    config_path = private / ('workflow-' + label + '.json')
    files = {alignment_path: tsv_bytes(alignments, ['sample_id', 'bam', 'receipt']),
             config_path: (json.dumps(cfg, indent=2, ensure_ascii=False) + '\n').encode()}
    if not callsets_manifest:
        files[interval_path] = intervals
    for path, data in files.items():
        require(not path.is_symlink() and (not path.exists() or path.read_bytes() == data),
                'An existing workflow plan differs; use another label.')
    if execute:
        for path, data in files.items():
            exclusive(path, data)
    result = {'status': 'PREPARED' if execute else 'DRY_RUN', 'samples': len(ids),
              'config': str(config_path), 'analysis_started': False, 'scope': scope,
              'configuration_sha256': canonical_hash(cfg)}
    print(json.dumps(result, ensure_ascii=False))
    return cfg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/srv/marina'))
    parser.add_argument('--resource-manifest', type=Path, required=True)
    parser.add_argument('--label', default='v1')
    parser.add_argument('--callsets-manifest', type=Path)
    parser.add_argument('--execute', action='store_true', help='Create private plan only, without calling variants.')
    parser.add_argument('--allow-pending-resources', action='store_true', help='Allow stages3-5 before databases finish; annotation still checks the final verified manifest.')
    args = parser.parse_args()
    os.umask(0o077)
    configure(args.root, args.resource_manifest, args.label, args.execute, args.callsets_manifest, args.allow_pending_resources)


if __name__ == '__main__':
    main()
