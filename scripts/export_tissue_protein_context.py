#!/usr/bin/env python3
"""Build public gene context, then join it locally into a separate gated sidecar."""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import re


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def dump(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write('\n')


def identity(value):
    value = str(value or '')
    if value.isdigit():
        value = 'HGNC:' + value
    return value if re.fullmatch(r'HGNC:\d+', value) else None


def version_match(transcript, mapping):
    if not re.fullmatch(r'ENST\d+\.\d+', transcript or ''):
        return 'VERSIONED_TRANSCRIPT_REQUIRED', []
    if transcript in mapping:
        return 'EXACT_TRANSCRIPT_VERSION_MATCH', mapping[transcript]
    base = transcript.split('.')[0]
    if any(key.split('.')[0] == base for key in mapping):
        return 'TRANSCRIPT_VERSION_MISMATCH', []
    return 'TRANSCRIPT_NOT_IN_SOURCE', []


def checked_source(path, receipt):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError('Missing or symlink source: ' + str(path))
    if path.stat().st_size != receipt['bytes'] or sha(path) != receipt['sha256']:
        raise ValueError('Source receipt mismatch: ' + str(path))
    return load(path)


def build(cfg, output):
    panel_path = Path(cfg['panel_manifest'])
    if sha(panel_path) != cfg['panel_sha256']:
        raise ValueError('Public panel SHA differs')
    panel = load(panel_path)
    genes = panel['genes']
    if panel.get('no_patient_selection') is not True or len(genes) != 39:
        raise ValueError('Only fixed public HGNC39 panel accepted')
    contexts = {gene['hgnc_id']: dict(hgnc_id=gene['hgnc_id'], gene=gene['symbol'],
                gtex_gene=[], gtex_transcripts={}, pext_transcripts={}, uniprot=[]) for gene in genes}
    if len(contexts) != 39 or len({g['symbol'] for g in genes}) != 39:
        raise ValueError('Duplicate public panel identity')
    source_root = Path(cfg['root'])
    manifests = {}
    verified = {}
    for profile in ('gtex', 'uniprot', 'pext'):
        manifest_path = source_root / profile / 'resource-manifest.json'
        manifest = load(manifest_path)
        if not manifest['status'].startswith('AVAILABLE_PUBLIC_'):
            raise ValueError('Public source is not complete: ' + profile)
        manifests[profile] = dict(path=str(manifest_path), sha256=sha(manifest_path), status=manifest['status'])
        # Raw source bytes, not derived files, are verified against fetch receipts.
        for receipt in manifest['source_receipts']:
            if profile == 'pext':
                symbol = receipt['query']['variables']['symbol']
                if symbol not in {g['symbol'] for g in genes}:
                    raise ValueError('Unapproved public query identity')
                path = source_root / profile / (symbol + '.json')
            else:
                receipt_files = list((source_root / profile).glob('*.json.receipt.json'))
                matches = [p for p in receipt_files if load(p).get('url') == receipt['url']]
                if len(matches) != 1:
                    raise ValueError('Ambiguous raw source receipt')
                path = Path(str(matches[0])[:-len('.receipt.json')])
            verified[str(path)] = dict(bytes=receipt['bytes'], sha256=receipt['sha256'])
            checked_source(path, receipt)
    gtex = load(source_root / 'gtex/resource-manifest.json')
    if (gtex['dataset'], gtex['gencode'], gtex['genome']) != ('gtex_v11', 'v47', 'GRCh38/hg38'):
        raise ValueError('GTEx version differs from reviewed contract')
    gtex_ids = {}
    for row in gtex['gene_reference_records']:
        context = contexts[row['hgnc_id']]
        if context['gene'] != row['gene'] or row['source_gencode_id'] in gtex_ids:
            raise ValueError('Conflicting exact GTEx gene identity')
        gtex_ids[row['source_gencode_id']] = context
        context['gtex_gencode_id'] = row['source_gencode_id']
    for endpoint, key in [('medianGeneExpression', 'gtex_gene'), ('medianTranscriptExpression', 'gtex_transcripts')]:
        for row in load(source_root / 'gtex' / (endpoint + '.json'))['data']:
            context = gtex_ids[row['gencodeId']]
            if row['geneSymbol'] != context['gene'] or row['unit'] != 'TPM' or row['datasetId'] != 'gtex_v11':
                raise ValueError('GTEx row identity/unit differs')
            entry = {k: row[k] for k in ('tissueSiteDetailId', 'ontologyId', 'median', 'unit', 'datasetId')}
            if key == 'gtex_gene':
                context[key].append(entry)
            else:
                transcript = row['transcriptId']
                if not re.fullmatch(r'ENST\d+\.\d+', transcript):
                    raise ValueError('GTEx transcript version missing')
                context[key].setdefault(transcript, []).append(entry)
    raw_proteins = load(source_root / 'uniprot/uniprot-reviewed-human-hcm39.json')['results']
    for protein in raw_proteins:
        if protein['organism']['taxonId'] != 9606 or protein['entryType'] != 'UniProtKB reviewed (Swiss-Prot)':
            raise ValueError('UniProt source is not reviewed human')
        for crossref in protein.get('uniProtKBCrossReferences', []):
            if crossref['database'] == 'HGNC' and crossref['id'] in contexts:
                contexts[crossref['id']]['uniprot'].append(dict(accession=protein['primaryAccession'],
                    release=load(source_root / 'uniprot/resource-manifest.json')['release'],
                    sequence_length=protein['sequence']['length'], sequence_sha256=hashlib.sha256(protein['sequence']['value'].encode()).hexdigest(),
                    features=len(protein.get('features', [])), entry_audit=protein.get('entryAudit'),
                    raw_source=str(source_root / 'uniprot/uniprot-reviewed-human-hcm39.json'),
                    scope='GENE_PROTEIN_CONTEXT_ONLY; exact sequence/isoform mapping required for allele-specific domains'))
    for gene in genes:
        context = contexts[gene['hgnc_id']]
        path = source_root / 'pext' / (gene['symbol'] + '.json')
        value = load(path)['data']['gene']
        context['pext_raw_source'] = str(path)
        if value is None:
            context['pext_status'] = 'MISSING_NOT_ZERO_OR_BENIGN'
            continue
        if identity(value['hgnc_id']) != gene['hgnc_id'] or value['symbol'] != gene['symbol'] or value['reference_genome'] != 'GRCh38':
            raise ValueError('pext public gene identity differs')
        for transcript in value.get('transcripts', []):
            key = transcript['transcript_id'] + '.' + str(transcript['transcript_version'])
            if key in context['pext_transcripts']:
                raise ValueError('Duplicate pext source transcript')
            context['pext_transcripts'][key] = [r for r in transcript.get('gtex_tissue_expression', [])
                                                if r['tissue'] in ('heart_left_ventricle', 'heart_atrial_appendage')]
        pext = value.get('pext') or {}
        context['pext_region_count'] = len(pext.get('regions') or [])
        context['pext_flags'] = pext.get('flags')
        context['pext_status'] = 'GENE_CONTEXT_AVAILABLE_ALLELE_COORDINATE_CONTRACT_UNVERIFIED' if context['pext_region_count'] else 'MISSING_NOT_ZERO_OR_BENIGN'
    output = Path(output)
    output.mkdir(mode=0o700, parents=True, exist_ok=False)
    index = dict(schema_version=1, status='PUBLIC_GENE_CONTEXT_READY_SEPARATE_SIDECAR_ONLY',
        public_panel_sha256=cfg['panel_sha256'], source_manifests=manifests, verified_raw_sources=verified,
        versions=dict(gtex='v11/GENCODE47/GRCh38', pext='gnomAD browser GRCh38; GTExv10 documented; API dataset release not returned',
                      uniprot=load(source_root / 'uniprot/resource-manifest.json')['release'], workflow='GENCODE49; no cross-version exon coordinate mapping'),
        genes=contexts, patient_input_used=False, allele_pext_assignment_enabled=False, allele_protein_domain_assignment_enabled=False)
    dump(output / 'context-index.json', index)
    with (output / 'gene-context.tsv').open('x', encoding='utf-8', newline='') as stream:
        fields = ['hgnc_id', 'gene', 'gtex_gene_context_json', 'gtex_versioned_transcripts', 'pext_status', 'pext_versioned_transcripts', 'uniprot_context_json']
        writer = csv.DictWriter(stream, fields, delimiter='\t'); writer.writeheader()
        for row in contexts.values():
            writer.writerow(dict(hgnc_id=row['hgnc_id'], gene=row['gene'], gtex_gene_context_json=json.dumps(row['gtex_gene']),
                gtex_versioned_transcripts=len(row['gtex_transcripts']), pext_status=row['pext_status'],
                pext_versioned_transcripts=len(row['pext_transcripts']), uniprot_context_json=json.dumps(row['uniprot'])))
    manifest = dict(status=index['status'], genes=len(contexts), source_manifests=manifests,
        files={p.name: dict(bytes=p.stat().st_size, sha256=sha(p)) for p in output.iterdir() if p.is_file()},
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), original_pipeline_changed=False)
    dump(output / 'manifest.json', manifest)
    return manifest


EXTRA = ['reference_context_hgnc', 'reference_context_scope', 'gtex_v11_gene_heart_medians_json',
         'gtex_v11_transcript_match_status', 'gtex_v11_transcript_heart_medians_json',
         'pext_gtex_v10_status', 'pext_source_transcript_match_status', 'pext_transcript_heart_expression_json',
         'pext_allele_assignment_status', 'uniprot_gene_protein_context_json', 'protein_domain_allele_assignment_status',
         'reference_context_source_row', 'reference_context_index_sha256']


def supplement(row, index, source_row):
    hgnc = identity(row.get('hgnc_id'))
    context = index['genes'].get(hgnc)
    if context is None:
        return None
    # Symbol is an additional consistency check, never a fallback for absent HGNC.
    if row.get('gene') and row['gene'] != context['gene']:
        raise ValueError('Annotation HGNC/symbol conflict at row ' + str(source_row))
    transcript = row.get('transcript', '')
    gtex_status, gtex_values = version_match(transcript, context['gtex_transcripts'])
    pext_status, pext_values = version_match(transcript, context['pext_transcripts'])
    return dict(reference_context_hgnc=hgnc, reference_context_scope='PUBLIC_REFERENCE_GENE_CONTEXT_NOT_INDIVIDUAL_EXPRESSION_OR_CLINICAL_CLASS',
        gtex_v11_gene_heart_medians_json=json.dumps(context['gtex_gene']), gtex_v11_transcript_match_status=gtex_status,
        gtex_v11_transcript_heart_medians_json=json.dumps(gtex_values), pext_gtex_v10_status=context['pext_status'],
        pext_source_transcript_match_status=pext_status, pext_transcript_heart_expression_json=json.dumps(pext_values),
        pext_allele_assignment_status='NOT_ASSIGNED_COORDINATE_TRANSCRIPT_CONTRACT_REQUIRES_REVIEW',
        uniprot_gene_protein_context_json=json.dumps(context['uniprot']),
        protein_domain_allele_assignment_status='NOT_ASSIGNED_EXACT_SEQUENCE_ISOFORM_MAPPING_REQUIRED', reference_context_source_row=source_row)


def enrich(workflow, context_path, output, deadline_utc='2026-11-06T03:46:48.046025+00:00'):
    workflow, context_path, output = Path(workflow).resolve(), Path(context_path).resolve(), Path(output).resolve()
    if output == workflow or workflow in output.parents:
        raise ValueError('Sidecar must be outside the immutable workflow root')
    receipt_paths = [workflow / 'stages' / stage / 'workflow-receipt.json' for stage in ('06_annotate', '07_report')]
    receipts = [load(p) for p in receipt_paths]
    receipt_hashes = {str(p): sha(p) for p in receipt_paths}
    deadline = min(dt.datetime.fromisoformat(deadline_utc), dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=8))
    if dt.datetime.now(dt.timezone.utc) >= deadline:
        raise ValueError('Absolute budget deadline reached')
    if any(r.get('status') != 'PASS' for r in receipts):
        raise ValueError('NOT_READY: both stage 6 and stage 7 must PASS')
    source = workflow / 'stages/06_annotate/transcript_samples.tsv.gz'
    expected = receipts[0].get('artifacts', {}).get(source.name, {})
    source_hash = sha(source)
    if source_hash != expected.get('sha256') or source.stat().st_size != expected.get('bytes'):
        raise ValueError('Completed source artifact differs from stage6 receipt')
    context_manifest = load(context_path.parent / 'manifest.json')
    checked_source(context_path, context_manifest['files'][context_path.name])
    index = load(context_path)
    index_hash = sha(context_path)
    if index.get('status') != 'PUBLIC_GENE_CONTEXT_READY_SEPARATE_SIDECAR_ONLY' or index.get('patient_input_used') is not False:
        raise ValueError('Unreviewed context index status')
    # Recheck raw-source pins before using normalized context, without any network.
    for path, receipt in index['verified_raw_sources'].items():
        checked_source(path, receipt)
    output.mkdir(mode=0o700, parents=True, exist_ok=False)
    counters = dict(input_rows=0, public39_matched_rows=0, samples={}, gtex_transcript_status={})
    target = output / 'public39-reference-context.tsv.gz'
    with gzip.open(source, 'rt', encoding='utf-8', newline='') as inp, gzip.open(target, 'wt', encoding='utf-8', newline='') as out:
        reader = csv.DictReader(inp, delimiter='\t')
        if not {'sample', 'hgnc_id', 'transcript', 'GT', 'FT'} <= set(reader.fieldnames or []):
            raise ValueError('Required per-sample annotation columns absent')
        if set(EXTRA) & set(reader.fieldnames):
            raise ValueError('Source already contains supplemental columns')
        writer = csv.DictWriter(out, reader.fieldnames + EXTRA, delimiter='\t'); writer.writeheader()
        for number, row in enumerate(reader, 1):
            if number % 10000 == 1 and dt.datetime.now(dt.timezone.utc) >= deadline:
                raise ValueError('Absolute budget deadline reached during sidecar export')
            counters['input_rows'] += 1
            extra = supplement(row, index, number)
            if extra is None:
                continue
            extra['reference_context_index_sha256'] = index_hash
            writer.writerow(dict(row, **extra))
            counters['public39_matched_rows'] += 1
            sample = row['sample']; counters['samples'][sample] = counters['samples'].get(sample, 0) + 1
            state = extra['gtex_v11_transcript_match_status']
            counters['gtex_transcript_status'][state] = counters['gtex_transcript_status'].get(state, 0) + 1
    # A changed source during export invalidates the sidecar, never alters original outputs.
    if sha(source) != source_hash or sha(context_path) != index_hash or any(sha(p) != receipt_hashes[str(p)] for p in receipt_paths):
        raise ValueError('Source changed while sidecar was being exported')
    manifest = dict(status='PASS_SEPARATE_RESEARCH_REFERENCE_CONTEXT_EXPORT', source_artifact=str(source),
        source_sha256=source_hash, source_bytes=source.stat().st_size, context_index_sha256=sha(context_path),
        stage_receipts={str(p): dict(sha256=receipt_hashes[str(p)]) for p in receipt_paths}, counters=counters,
        artifact=dict(path=str(target), bytes=target.stat().st_size, sha256=sha(target)),
        scope='Matched public HGNC39 rows only; each original per-sample row and GT/FT/quality retained; no clinical burden or mechanism assignment',
        original_pipeline_changed=False, network_requests=0, created_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    dump(output / 'manifest.json', manifest)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('build'); p.add_argument('--config', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('enrich'); p.add_argument('--workflow-root', type=Path, required=True); p.add_argument('--context-index', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = build(load(args.config), args.output) if args.command == 'build' else enrich(args.workflow_root, args.context_index, args.output)
    print(json.dumps(result, indent=2))
