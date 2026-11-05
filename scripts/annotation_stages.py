#!/usr/bin/env python3
"""Offline transcript annotation and transparent research evidence reports.

Public interface: run_stage('annotate'|'report', cfg, workdir, runtime).
The workflow owns stage atomics, input fingerprints and completed-stage receipts.
No variant, genotype or phenotype is submitted to an external service.
"""
import csv
from contextlib import ExitStack
import gzip
import hashlib
import html
import importlib.util
import json
import math
from pathlib import Path
import re
import sqlite3
from urllib.parse import unquote


PLOF = {'stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant'}
FUNCTIONAL = PLOF | {'start_lost', 'stop_lost', 'missense_variant', 'inframe_insertion',
                     'inframe_deletion', 'protein_altering_variant', 'splice_region_variant'}
ACMG = ['PVS1', *['PS' + str(i) for i in range(1, 5)], *['PM' + str(i) for i in range(1, 7)],
        *['PP' + str(i) for i in range(1, 6)], 'BA1', *['BS' + str(i) for i in range(1, 5)],
        *['BP' + str(i) for i in range(1, 8)]]
MANUAL_SOURCES = ['HGMD', 'OMIM', 'dbNSFP', 'GeneReviews', 'LOVD', 'Franklin', 'VarSome',
                  'VariantValidator', 'DECIPHER', 'DGV', 'MITOMAP', 'UniProt', 'HGNC',
                  'HPO', 'UCSC', 'NCBI', 'dbSNP', 'pathways', 'ExAC', 'RuSeq', 'All_of_Us',
                  'Orphanet', 'SAMS', 'REVEL', 'AlphaMissense', 'SpliceAI', 'Varity',
                  'MT_unresolved_label', 'MUT_Assessor_unresolved_label', 'LRG', 'Mutalyzer',
                  'COSMIC', 'MyCancerGenome', 'DoCM', 'PharmGKB', 'DGIdb', 'OncoKB', 'JAX_CKB',
                  'pLI_LOEUF', 'pHaplo_pTriplo', 'ACMG_secondary_findings_consent_dependent']
FLAG_TERMS = {'stop_gained': {'stop_gained'}, 'stop_lost': {'stop_lost'}, 'start_lost': {'start_lost'},
    'frameshift': {'frameshift_variant'}, 'splice_donor_acceptor': {'splice_donor_variant', 'splice_acceptor_variant'},
    'missense': {'missense_variant'}, 'synonymous': {'synonymous_variant'}, 'plof_candidate': PLOF}
PRIMARY_ORDER = ['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant',
                 'start_lost', 'stop_lost', 'inframe_insertion', 'inframe_deletion', 'missense_variant',
                 'protein_altering_variant', 'splice_region_variant', 'synonymous_variant']
FIELDS = [
    'sample', 'variant_id', 'chrom', 'pos', 'ref', 'alt', 'variant_type', 'GT', 'FT', 'DP', 'AD', 'GQ', 'QUAL', 'AD_derived_VAF', 'sample_state',
    'site_filter', 'phase_set', 'gene', 'gene_id', 'hgnc_id', 'transcript', 'transcript_source', 'protein_id', 'amino_acids', 'codons', 'strand', 'cds_position', 'consequence', 'impact',
    'exon', 'intron', 'hgvs_g', 'hgvs_g_status', 'hgvs_c', 'hgvs_p', 'protein_position', 'protein_length',
    'protein_region_heuristic', 'nmd_assessment', 'mechanism_status', 'mane_select', 'canonical', 'biotype',
    'stop_gain', 'stop_loss', 'start_lost', 'frameshift', 'splice', 'synonymous', 'missense',
    'plof_candidate', 'SIFT', 'PolyPhen', 'FATHMM', 'CADD_PHRED',
    'predictors_available', 'predictors_damaging', 'research_vote_3of4',
    'fourth_predictor_identity', 'af_max_observed', 'af_status', 'af_fields_json',
    'ancestry_af_json', 'gnomad_version', 'clinvar_significance', 'clinvar_review',
    'clinvar_conflicts', 'clinvar_conditions', 'clinvar_ids', 'clinvar_allele_id', 'clinvar_disease_ids', 'clinvar_presence', 'clinvar_snapshot',
    'clinvar_assertion_dates_status', 'clinvar_cache_historical', 'clingen_evidence_json', 'gene_context_json', 'local_evidence_json', 'vep_all_csq_json',
    'hpo_exact_term_matches', 'resource_versions_json', 'candidate_eligible',
    'candidate_reason', 'external_pathogenic_assertion', 'quality_rank',
    'gene_validity_rank', 'gene_validity_any_context', 'gene_validity_context_status', 'phenotype_context_rank', 'target_gene_context_rank', 'clinvar_review_rank',
    'consequence_rank', 'population_rank', 'ranking_rule', 'classification_status']
FIELDS += ['AlphaMissense_score', 'AlphaMissense_model_label', 'AlphaMissense_match_status',
           'AlphaMissense_source_transcript', 'AlphaMissense_source_protein_variant',
           'AlphaMissense_source_uniprot', 'AlphaMissense_source_dataset',
           'AlphaMissense_version', 'AlphaMissense_research_rank', 'AlphaMissense_allele_predictions_json']
_ALPHA = None


def alpha_module():
    global _ALPHA
    if _ALPHA is None:
        path = Path(__file__).resolve().parent / 'alphamissense_adapter.py'
        spec = importlib.util.spec_from_file_location('_marina_alphamissense_adapter', path)
        _ALPHA = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_ALPHA)
    return _ALPHA


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def jdump(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False)


def text_open(path, mode='rt'):
    return gzip.open(path, mode, encoding='utf-8', newline='') if str(path).endswith('.gz') else Path(path).open(mode.replace('t', ''), encoding='utf-8', newline='')


def write_json(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False)
        stream.write('\n')


def writer(stack, path, fields):
    stream = stack.enter_context(text_open(path, 'xt'))
    result = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
    result.writeheader()
    return result


def sample_state(gt, ft, site_filter='PASS'):
    """A missing/failed FT is never silently treated as a confidently reference sample."""
    alleles = re.split(r'[/|]', gt or '.')
    if any(a == '.' for a in alleles):
        return 'unknown_no_call'
    need(all(a in ('0', '1') for a in alleles), 'Expected decomposed biallelic GT with alleles 0/1 only.')
    alternate = '1' in alleles
    if ft in ('', '.', None):
        return 'unknown_quality_alt' if alternate else 'unknown_quality_reference'
    if ft != 'PASS' or site_filter != 'PASS':
        return 'filtered_alt' if alternate else 'filtered_reference'
    return 'alternate' if alternate else 'reference'


def protein_region(value):
    """Research-only thirds of transcript protein length, not a functional-domain/NMD claim."""
    try:
        position, length = value.split('/')
        first = int(position.split('-')[0])
        total = int(length)
        if total <= 0 or first <= 0 or first > total:
            return '', '', 'unknown'
        fraction = first / total
        return position, length, 'beginning' if fraction <= 1 / 3 else 'middle' if fraction <= 2 / 3 else 'end'
    except (ValueError, AttributeError):
        return value or '', '', 'unknown'


def predictor_vote(sift='', polyphen='', fathmm='', cadd='', cadd_threshold=20.0, fourth='CADD_unconfirmed'):
    results = []
    for value, positive, negative in [
        (sift, ('deleterious', 'd'), ('tolerated', 't')),
        (polyphen, ('probably_damaging', 'possibly_damaging', 'd', 'p'), ('benign', 'b')),
        (fathmm, ('damaging', 'd'), ('tolerated', 't', 'neutral', 'n'))]:
        label = (value or '').lower().split('(')[0].strip()
        results.append(True if label in positive else False if label in negative else None)
    try:
        cadd_vote = float(cadd) >= cadd_threshold if cadd not in ('', '.', None) else None
    except (ValueError, TypeError):
        cadd_vote = None
    results.append(cadd_vote if fourth.startswith('CADD') else None)
    available = sum(v is not None for v in results)
    damaging = sum(v is True for v in results)
    vote = 'yes_research_only' if damaging >= 3 else 'no' if available == 4 else 'incomplete_not_negative'
    return available, damaging, vote


def normalize_chrom(chrom):
    value = chrom[3:] if chrom.startswith('chr') else chrom
    return 'MT' if value == 'M' else value


def variant_key(chrom, pos, ref, alt):
    return ':'.join((normalize_chrom(chrom), str(int(pos)), ref.upper(), alt.upper()))


def sequence_allele(ref, alt):
    return bool(re.fullmatch('[ACGTacgt]+', ref) and re.fullmatch('[ACGTacgt]+', alt))


def record_state(record, data):
    state = sample_state(data.get('GT', '.'), data.get('FT', '.'), record['filter'])
    return state if sequence_allele(record['ref'], record['alt']) else 'not_sequence_allele_placeholder'


def vcf_samples(path):
    with text_open(path) as stream:
        for line in stream:
            if line.startswith('#CHROM'):
                return line.rstrip('\r\n').split('\t')[9:]
    raise ValueError('VCF sample header missing.')


def partition_sequence_input(source, destination):
    count = 0
    with text_open(source) as inp, text_open(destination, 'xt') as out:
        for line in inp:
            if line.startswith('#'):
                out.write(line)
            else:
                fields = line.rstrip('\r\n').split('\t')
                need(',' not in fields[4], 'Expected decomposed biallelic input before annotation.')
                if sequence_allele(fields[3], fields[4]):
                    out.write(line)
                    count += 1
    return count


def merge_vep_preserving_all(source, vep_output, destination):
    """Bypass *, symbolic/breakend/nonvariant records and retain the original soft labels.

    Require a one-to-one ordered result for real alleles. VEP skipping a real input
    (e.g. reference mismatch) is a failed stage, not a silently negative finding.
    """
    with text_open(vep_output) as annotated, text_open(source) as original, destination.open('x') as out:
        annotated_header, first = [], None
        for line in annotated:
            if line.startswith('#'):
                annotated_header.append(line)
            else:
                first = line
                break
        need(any(line.startswith('#CHROM') for line in annotated_header), 'VEP output lacks sample header.')
        for line in annotated_header:
            if line.startswith('#CHROM'):
                out.write('##INFO=<ID=ANNOTATION_STATUS,Number=1,Type=String,Description="Nonsequence input retained without small-variant annotation">\n')
            out.write(line)
        for line in original:
            if line.startswith('#'):
                continue
            raw = line.rstrip('\r\n').split('\t')
            if sequence_allele(raw[3], raw[4]):
                need(first is not None, 'VEP omitted a real input allele; incomplete annotation refused.')
                fields = first.rstrip('\r\n').split('\t')
                need(fields[:2] == raw[:2] and fields[3:5] == raw[3:5] and fields[8:] == raw[8:],
                     'VEP changed input allele order or sample GT/FT; refusing merge.')
                out.write(first)
                first = next(annotated, None)
            else:
                raw[7] = ('' if raw[7] == '.' else raw[7] + ';') + 'ANNOTATION_STATUS=not_sequence_allele'
                out.write('\t'.join(raw) + '\n')
        need(first is None, 'VEP produced an unexpected extra allele.')


def resource_path(manifest_path, item):
    path = Path(item['path'])
    if path.is_absolute():
        need(item.get('local_adapter') is True, 'Absolute resource paths require an explicit local adapter.')
        return path.resolve(strict=True)
    root = Path(manifest_path).parent.resolve(strict=True)
    result = (root / path).resolve(strict=True)
    need(result.is_relative_to(root), 'Resource path escapes its prepared bundle.')
    return result


def load_resources(path, verify_tree=True, local_adapters=None):
    path = Path(path)
    manifest = json.loads(path.read_text(encoding='utf-8-sig'))
    need(manifest.get('schema_version') == 1 and manifest.get('status') == 'PASS'
         and manifest.get('assembly') == 'GRCh38', 'Prepared GRCh38 resource manifest is required.')
    for adapter in local_adapters or []:
        name = adapter.get('name', '')
        need(re.fullmatch(r'[A-Za-z0-9_.-]+', name) and adapter.get('local_adapter') is True
             and adapter.get('status') == 'available' and Path(adapter.get('path', '')).is_absolute()
             and adapter.get('license_status') and adapter.get('source_url'),
             'Local adapter requires name, absolute path, availability, provenance and explicit license status.')
        need(name not in {'clinvar', 'clinvar_index', 'vep_cache', 'clingen_validity', 'clingen_dosage',
                          'sift', 'polyphen2', 'mane_cache', 'gnomad', 'hgnc', 'mane', 'hpo_ontology'},
             'Attach additional local evidence under a new name; do not replace locked core annotation resources.')
        format_name = adapter.get('format')
        need(format_name == 'variant_evidence_tsv' or (name, format_name) in {
            ('gene_evidence', 'gene_evidence_tsv'), ('hpo_gene', 'hpo_gene_tsv'), ('pathways', 'pathway_tsv')},
            'Unsupported local adapter schema/name; refusing to label an unparsed file as loaded.')
        prior = manifest['resources'].get(name, {})
        manifest['resources'][name] = dict(adapter, previous_status=prior.get('status'))
    for name, item in manifest['resources'].items():
        if item.get('status') == 'available':
            need(item.get('version') and re.fullmatch('[0-9a-f]{64}', item.get('sha256', '')),
                 'Available resources require version and SHA256: ' + name)
            source = resource_path(path, item)
            need(source.stat().st_size == item['bytes'] and sha256(source) == item['sha256'],
                 'Annotation resource differs from its preparation receipt: ' + name)
    if verify_tree:
        for item in manifest.get('artifacts', []):
            source = resource_path(path, item)
            need(source.stat().st_size == item['bytes'] and sha256(source) == item['sha256'],
                 'Prepared cache/resource artifact changed.')
    if all(manifest['resources'].get(key, {}).get('status') == 'available' for key in ('hpo_gene', 'hpo_ontology')):
        previous = manifest['resources'].get('HPO', {}).get('status')
        manifest['resources']['HPO'] = dict(status='available_gene_associations_and_ontology; patient_phenotype_separate',
            version=manifest['resources']['hpo_gene']['version'], previous_status=previous,
            components={key: {'version': manifest['resources'][key]['version'], 'sha256': manifest['resources'][key]['sha256']}
                        for key in ('hpo_gene', 'hpo_ontology')})
    for name in MANUAL_SOURCES:
        if not any(key.casefold() == name.casefold() for key in manifest['resources']):
            manifest['resources'][name] = dict(status='unavailable_manual_or_local_adapter_required', version=None)
    return manifest


def requirements_coverage(config, manifest):
    """Account for every audited source without claiming every database is installed."""
    path = Path(config.get('evidence_registry') or Path(__file__).resolve().parents[1] / 'config/evidence-sources.example.json')
    if not path.is_file():
        return dict(status='requirements_registry_not_supplied', sources=[])
    registry = json.loads(path.read_text())
    aliases = {'clingen_gene_validity': ['clingen_validity'], 'hpo': ['hpo_gene', 'hpo_ontology'],
        'ncbi_refseq_gene_mane': ['mane'], 'ensembl': ['vep_cache'], 'gnomad': ['vep_cache'],
        'sift': ['vep_cache'], 'polyphen2': ['vep_cache'], 'dbsnp': ['vep_cache']}
    rows = []
    for requirement in registry['resources']:
        resource_id = requirement['id']
        keys = aliases.get(resource_id, [resource_id])
        matches = [(key, manifest['resources'].get(key, {})) for key in keys]
        available = [key for key, item in matches if item.get('status') in ('available', 'bundled_in_verified_cache', 'available_local_supplement')]
        status = 'available_local' if len(available) == len(keys) else 'partially_available_local' if available else 'unavailable_or_manual'
        if requirement.get('scope') in ('somatic_only', 'somatic_pharmacogenomics', 'mitochondrial', 'cnv_sv') and not available:
            status = 'separate_branch_not_analysed'
        rows.append(dict(resource_id=resource_id, role=requirement.get('role'), scope=requirement.get('scope'),
            status=status, actual_resource_keys=available,
            versions={key: manifest['resources'][key].get('version') for key in available},
            source_url=requirement.get('official_url'), access=requirement.get('access'),
            license_review_status=requirement.get('license_review_status'),
            limitation=('cache_AF_only; no AF=0/AN/coverage/constraint assumed' if resource_id == 'gnomad' else
                        'rs identifiers are not a clinical classification' if resource_id == 'dbsnp' else
                        'summary metadata/transcript labels; no complete clinical disease narrative claimed' if resource_id == 'ncbi_refseq_gene_mane' else None)))
    return dict(status='audited_registry_mapped_to_actual_resources', registry_sha256=sha256(path),
                source_count=len(rows), sources=rows)


def gene_evidence(manifest, manifest_path):
    genes = {}
    for name in ('clingen_validity', 'clingen_dosage', 'gene_evidence', 'hgnc', 'mane'):
        item = manifest['resources'].get(name, {})
        if item.get('status') != 'available':
            continue
        path = resource_path(manifest_path, item)
        with text_open(path) as stream:
            if name.startswith('clingen_'):
                rows = csv.reader(stream)
                header = next((r for r in rows if r and r[0] == 'GENE SYMBOL'), None)
                need(header is not None, 'ClinGen CSV header changed.')
                records = (dict(zip(header, row)) for row in rows if row and not row[0].startswith('+'))
            else:
                records = csv.DictReader(stream, delimiter='\t')
            for row in records:
                gene = row.get('GENE SYMBOL') or row.get('gene') or row.get('symbol')
                if gene:
                    if name == 'hgnc':
                        row = {key: row.get(key, '') for key in ('hgnc_id', 'symbol', 'name', 'status', 'alias_symbol', 'prev_symbol', 'mane_select', 'refseq_accession', 'ensembl_gene_id', 'uniprot_ids')}
                    genes.setdefault(gene, []).append(dict(source=name, resource_version=item['version'], evidence=row))
    return genes


def make_adapter_database(manifest, manifest_path, destination):
    """Local normalized exports only; no dbNSFP/HGMD licensing is implied by this adapter."""
    connection = sqlite3.connect(destination)
    connection.execute('CREATE TABLE evidence (variant TEXT, transcript TEXT, evidence TEXT)')
    for name, item in manifest['resources'].items():
        if item.get('status') != 'available' or item.get('format') != 'variant_evidence_tsv':
            continue
        with text_open(resource_path(manifest_path, item)) as stream:
            rows = csv.DictReader(stream, delimiter='\t')
            need({'chrom', 'pos', 'ref', 'alt'}.issubset(rows.fieldnames or []), 'Local adapter lacks exact-allele keys.')
            for row in rows:
                need(',' not in row['alt'] and row['alt'] not in ('.', '*'), 'Adapter must be normalized and biallelic.')
                evidence = dict(row, resource=name, resource_version=item['version'], resource_sha256=item['sha256'])
                connection.execute('INSERT INTO evidence VALUES (?,?,?)',
                    (variant_key(row['chrom'], row['pos'], row['ref'], row['alt']), row.get('transcript') or '*', jdump(evidence)))
    connection.execute('CREATE INDEX exact_variant_transcript ON evidence (variant,transcript)')
    connection.commit()
    return connection


def local_evidence(connection, key, transcript):
    return [json.loads(row[0]) for row in connection.execute(
        'SELECT evidence FROM evidence WHERE variant=? AND transcript IN (?,?)', (key, transcript, '*'))]


def hpo_map(manifest, manifest_path):
    item = manifest['resources'].get('hpo_gene', {})
    mapping = {}
    if item.get('status') == 'available':
        with text_open(resource_path(manifest_path, item)) as stream:
            rows = csv.DictReader(stream, delimiter='\t')
            need('hpo_id' in (rows.fieldnames or []) and ('gene' in rows.fieldnames or 'gene_symbol' in rows.fieldnames), 'HPO adapter needs gene/gene_symbol and hpo_id.')
            for row in rows:
                mapping.setdefault(row.get('gene') or row['gene_symbol'], set()).add(row['hpo_id'])
    return mapping


def parse_vcf(path):
    """Yield all records including soft-filtered sites; no genotype is discarded here."""
    with text_open(path) as stream:
        csq_fields, samples = [], None
        previous, seen_at_position, completed_contigs = None, set(), set()
        for line in stream:
            if line.startswith('##INFO=<ID=CSQ,'):
                match = re.search(r'Format: ([^"\n]+)', line)
                need(match is not None, 'Cannot parse VEP CSQ schema.')
                csq_fields = match.group(1).rstrip('>').split('|')
            elif line.startswith('#CHROM'):
                samples = line.rstrip('\n\r').split('\t')[9:]
                need(samples and len(set(samples)) == len(samples), 'Unique VCF sample names are required.')
                need(all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}', s) for s in samples), 'Unsafe report sample name.')
            elif not line.startswith('#'):
                need(samples is not None, 'VCF sample header missing.')
                fields = line.rstrip('\n\r').split('\t')
                need(len(fields) == 9 + len(samples), 'VCF column count mismatch.')
                chrom, pos, vid, ref, alt, qual, filt, info, fmt = fields[:9]
                coordinate = chrom, int(pos)
                if previous != coordinate:
                    if previous is not None:
                        need(chrom != previous[0] or int(pos) >= previous[1], 'Annotation VCF must be coordinate sorted.')
                        if chrom != previous[0]:
                            completed_contigs.add(previous[0])
                            need(chrom not in completed_contigs, 'Annotation VCF contig order repeats a completed block.')
                    previous, seen_at_position = coordinate, set()
                need((ref, alt) not in seen_at_position, 'Duplicate normalized allele in annotation VCF.')
                seen_at_position.add((ref, alt))
                need(',' not in alt, 'Annotation requires decomposed biallelic records; placeholders are preserved separately.')
                infos = dict((entry.split('=', 1) + [''])[:2] for entry in info.split(';'))
                formats = fmt.split(':')
                need('GT' in formats, 'Per-sample GT is required; missing FT remains unknown quality.')
                genotypes = dict(zip(samples, (dict(zip(formats, value.split(':'))) for value in fields[9:])))
                consequences = []
                for entry in infos.get('CSQ', '').split(','):
                    if entry:
                        values = entry.split('|')
                        need(len(values) == len(csq_fields), 'CSQ value/schema length differs.')
                        row = dict(zip(csq_fields, map(unquote, values)))
                        need(row.get('ALLELE_NUM', '1') == '1', 'CSQ allele mapping differs from decomposed ALT.')
                        consequences.append(row)
                yield dict(chrom=chrom, pos=pos, ref=ref, alt=alt, qual=qual, filter=filt,
                           samples=genotypes, csq=consequences or [{}])


def annotation_row(record, csq, sample, data, context):
    config, manifest, genes, hpogenes, db = context[:5]
    alpha = context[5] if len(context) > 5 else None
    key = variant_key(record['chrom'], record['pos'], record['ref'], record['alt'])
    local = local_evidence(db, key, csq.get('Feature', ''))
    # Multiple conflicting adapter assertions stay visible; ambiguous predictor values are not voted.
    def one(field):
        values = {entry[field] for entry in local if entry.get(field) not in ('', '.', None)}
        return next(iter(values)) if len(values) == 1 else ''
    sift, polyphen = csq.get('SIFT', ''), csq.get('PolyPhen', '')
    fathmm, cadd = one('FATHMM'), one('CADD_PHRED')
    fourth = config.get('fourth_predictor', 'CADD_unconfirmed')
    available, damaging, vote = predictor_vote(sift, polyphen, fathmm, cadd,
        config.get('cadd_research_threshold', 20), fourth)
    afs, ancestry = {}, {}
    for field, value in csq.items():
        if field.startswith(('gnomADe_', 'gnomADg_')) and field.endswith('AF') and value not in ('', '.'):
            try:
                number = float(value)
                need(math.isfinite(number) and 0 <= number <= 1, 'Invalid cached population AF.')
                afs['cache:' + field] = number
                if field not in ('gnomADe_AF', 'gnomADg_AF'):
                    ancestry['cache:' + field] = number
            except ValueError:
                raise ValueError('Malformed cached population AF; no rarity assigned.') from None
    for evidence in local:
        for field, value in evidence.items():
            if field.startswith('gnomAD_') and field.endswith('AF') and value not in ('', '.', None):
                try:
                    number = float(value)
                    need(0 <= number <= 1, 'Invalid local population AF.')
                    afs[evidence['resource'] + ':' + field] = number
                    if field != 'gnomAD_AF':
                        ancestry[evidence['resource'] + ':' + field] = number
                except ValueError:
                    if value not in ('', '.'):
                        raise
    maximum = max(afs.values()) if afs else None
    rare = maximum is not None and maximum <= config.get('research_max_af', 0.01)
    state = record_state(record, data)
    consequences = set(csq.get('Consequence', '').split('&'))
    gene = csq.get('SYMBOL', '')
    evidence = genes.get(gene, [])
    strong = [e for e in evidence if e['source'] == 'clingen_validity' and e['evidence'].get('CLASSIFICATION') in ('Definitive', 'Strong')]
    supplied_disease = config.get('disease_context_by_sample', {}).get(sample, config.get('disease_context', {}))
    disease = supplied_disease if isinstance(supplied_disease, dict) else {}
    def matches_context(entry):
        row = entry['evidence']
        selected = (row.get('DISEASE ID (MONDO)') in disease.get('mondo_ids', [])
                    or row.get('DISEASE LABEL', '').lower() in {label.lower() for label in disease.get('disease_labels', [])})
        return selected and (not disease.get('inheritance') or row.get('MOI') in disease['inheritance'])
    validity = any(matches_context(e) for e in strong)
    hpo_terms = config.get('hpo_terms_by_sample', {}).get(sample, config.get('hpo_terms', []))
    hpo = sorted(set(hpo_terms) & hpogenes.get(gene, set()))
    significance = csq.get('ClinVar_CLNSIG', '')
    review = csq.get('ClinVar_CLNREVSTAT', '')
    conflict = csq.get('ClinVar_CLNSIGCONF', '')
    external_path = bool(re.search(r'(^|[&|/])(?:Pathogenic|Likely_pathogenic)(?:_low_penetrance)?(?=$|[&|/])', significance, re.I)
                         or re.search(r'(?:^|[&|/])(?:Pathogenic|Likely_pathogenic)\(', conflict, re.I))
    review_rank = (4 if 'practice_guideline' in review else 3 if 'reviewed_by_expert_panel' in review else
                   2 if 'multiple_submitters' in review and 'no_conflicts' in review else 1 if 'criteria_provided' in review else 0)
    if conflict or 'conflicting' in significance.lower() or 'conflicting' in review.lower():
        review_rank = 0
    functional = bool(consequences & FUNCTIONAL)
    eligible = state == 'alternate' and (external_path or (functional and (rare or maximum is None)))
    reason = ('external_ClinVar_assertion_for_manual_review' if external_path else
              'functional_and_observed_AF_below_research_threshold' if rare else
              'functional_with_population_data_missing_requires_review' if maximum is None and functional else 'not_research_candidate')
    if state != 'alternate':
        reason = 'not_quality_passing_alternate:' + state
    position, length, region = protein_region(csq.get('Protein_position', ''))
    try:
        depths = [int(value) for value in data.get('AD', '').split(',')]
        derived_vaf = depths[1] / sum(depths) if len(depths) == 2 and all(n >= 0 for n in depths) and sum(depths) else ''
    except ValueError:
        derived_vaf = ''
    variant_type = ('spanning_deletion_placeholder' if record['alt'] == '*' else 'nonsequence_placeholder_or_symbolic' if not sequence_allele(record['ref'], record['alt']) else
                    'SNV' if len(record['ref']) == len(record['alt']) == 1 else
                    'insertion' if len(record['alt']) > len(record['ref']) else
                    'deletion' if len(record['ref']) > len(record['alt']) else 'MNV_or_complex')
    result = dict(sample=sample, variant_id=key, chrom=record['chrom'], pos=record['pos'],
        ref=record['ref'], alt=record['alt'], variant_type=variant_type, GT=data.get('GT', '.'), FT=data.get('FT', '.'),
        DP=data.get('DP', '.'), AD=data.get('AD', '.'), GQ=data.get('GQ', '.'), QUAL=record.get('qual', '.'), AD_derived_VAF=derived_vaf,
        sample_state=state, site_filter=record['filter'], phase_set=data.get('PS', ''),
        gene=gene, gene_id=csq.get('Gene', ''), hgnc_id=csq.get('HGNC_ID', ''), transcript=csq.get('Feature', ''),
        transcript_source=csq.get('SOURCE') or 'Ensembl115_cache', protein_id=csq.get('ENSP', ''),
        amino_acids=csq.get('Amino_acids', ''), codons=csq.get('Codons', ''), strand=csq.get('STRAND', ''), cds_position=csq.get('CDS_position', ''),
        consequence=csq.get('Consequence', ''), impact=csq.get('IMPACT', ''), exon=csq.get('EXON', ''),
        intron=csq.get('INTRON', ''), hgvs_g=csq.get('HGVSg') or one('HGVSg'),
        hgvs_g_status='VEP_or_local_adapter_generated; reference_accession_and_nomenclature_require_review',
        hgvs_c=csq.get('HGVSc', ''), hgvs_p=csq.get('HGVSp', ''),
        protein_position=position, protein_length=length, protein_region_heuristic=region,
        nmd_assessment='not_assessed; manual_transcript_exon_and_mechanism_review_required',
        mechanism_status='not_assessed; no_automatic_LoF_GoF_or_PVS1',
        mane_select=csq.get('MANE_SELECT', ''), canonical=csq.get('CANONICAL', ''), biotype=csq.get('BIOTYPE', ''),
        stop_gain=int('stop_gained' in consequences), stop_loss=int('stop_lost' in consequences),
        start_lost=int('start_lost' in consequences), frameshift=int('frameshift_variant' in consequences),
        splice=int(any(term.startswith('splice_') for term in consequences)),
        synonymous=int('synonymous_variant' in consequences), missense=int('missense_variant' in consequences),
        plof_candidate=int(bool(consequences & PLOF)), SIFT=sift, PolyPhen=polyphen, FATHMM=fathmm,
        CADD_PHRED=cadd, predictors_available=available, predictors_damaging=damaging,
        research_vote_3of4=vote, fourth_predictor_identity=fourth,
        af_max_observed='' if maximum is None else maximum,
        af_status='missing_not_zero' if maximum is None else 'observed_below_research_threshold' if rare else 'observed_above_research_threshold',
        af_fields_json=jdump(afs), ancestry_af_json=jdump(ancestry),
        gnomad_version=manifest['resources'].get('gnomad', {}).get('version', 'cache info.txt; inspect resource report'),
        clinvar_significance=significance, clinvar_review=review, clinvar_conflicts=conflict,
        clinvar_conditions=csq.get('ClinVar_CLNDN', ''), clinvar_ids=csq.get('ClinVar', ''),
        clinvar_allele_id=csq.get('ClinVar_ALLELEID', ''), clinvar_disease_ids=csq.get('ClinVar_CLNDISDB', ''),
        clinvar_presence='exact_allele_summary_match' if csq.get('ClinVar') else 'no_exact_summary_match' if manifest['resources'].get('clinvar', {}).get('status') == 'available' else 'resource_unavailable',
        clinvar_snapshot=manifest['resources'].get('clinvar', {}).get('version', 'unavailable'),
        clinvar_assertion_dates_status='not_provided_by_summary_VCF; inspect dated assertion export',
        clinvar_cache_historical=csq.get('CLIN_SIG', ''),
        clingen_evidence_json=jdump([e for e in evidence if e['source'].startswith('clingen_')]),
        gene_context_json=jdump([e for e in evidence if not e['source'].startswith('clingen_')]),
        local_evidence_json=jdump(local), vep_all_csq_json=jdump(csq),
        hpo_exact_term_matches='&'.join(hpo), resource_versions_json=jdump({n: i.get('version') for n, i in manifest['resources'].items()}),
        candidate_eligible=int(eligible), candidate_reason=reason, external_pathogenic_assertion=int(external_path),
        quality_rank=int(state == 'alternate'), gene_validity_rank=int(validity), gene_validity_any_context=int(bool(strong)),
        gene_validity_context_status='explicit_context_match' if validity else 'no_matching_explicit_context' if disease else 'descriptive_context_without_controlled_disease_ids; agnostic_evidence_only' if supplied_disease else 'disease_context_not_supplied; agnostic_evidence_only',
        phenotype_context_rank=int(bool(hpo)),
        target_gene_context_rank=int(gene in config.get('target_genes', [])),
        clinvar_review_rank=review_rank if external_path else 0,
        consequence_rank=3 if consequences & PLOF else 2 if functional else 0,
        population_rank=2 if rare else 0 if maximum is not None else 1,
        ranking_rule='lexicographic_v1; phenotype,target_gene_context,review,genevalidity,consequence,population,predictorvotes,AlphaMissense_compatible_research_rank,MANE; identifiers break ties',
        classification_status='manual_ACMG_review_required; no_automatic_pathogenicity_or_GoF_or_PVS1')
    if alpha:
        result.update(alpha.match(key, csq, config.get('alphamissense_research_threshold', .564)))
    else:
        result.update(AlphaMissense_match_status='not_prepared_in_this_context; not_benign', AlphaMissense_research_rank=0)
    return result


def make_tables(annotated, workdir, config, manifest, manifest_path, alpha_inputs=None):
    dbpath = workdir / 'local-adapter-index.sqlite'
    db = make_adapter_database(manifest, manifest_path, dbpath)
    alpha = alpha_module().Lookup(db, alpha_inputs[1], alpha_inputs[0]) if alpha_inputs else None
    context = config, manifest, gene_evidence(manifest, manifest_path), hpo_map(manifest, manifest_path), db, alpha
    names = vcf_samples(annotated)
    counts = dict(records=0, transcript_sample_rows=0, missing_csq=0, nonsequence_records_preserved=0, samples=names,
                  alphamissense_transcript_sample_statuses={}, alphamissense_matched_source_rows=alpha.source_rows if alpha else 0)
    def empty_counter():
        return dict(unique_normalized_alleles=0, nonexclusive_flags={flag: 0 for flag in FLAG_TERMS}, exclusive_primary_consequence={}, variant_type_histogram={})
    counters = dict(schema_version=1, cohort={level: empty_counter() for level in ('raw', 'qc', 'candidate')}, samples={},
        definitions={'raw': 'observed ALT allele in GT, including partial/filtered calls; not a quality claim',
                     'qc': 'GT fully called, sample FT PASS and site FILTER PASS, ALT present',
                     'candidate': 'at least one transcript meets documented research candidate rule',
                     'count_unit': 'unique normalized CHROM:POS:REF:ALT per sample or cohort; transcripts unioned',
                     'primary_order': PRIMARY_ORDER, 'flags': 'nonexclusive; do not sum flags as variant total'})
    for sample in names:
        counters['samples'][sample] = {level: empty_counter() for level in ('raw', 'qc', 'candidate')}
    def increment(counter, terms, record):
        counter['unique_normalized_alleles'] += 1
        for flag, relevant in FLAG_TERMS.items():
            counter['nonexclusive_flags'][flag] += int(bool(terms & relevant))
        primary = next((term for term in PRIMARY_ORDER if term in terms), sorted(terms - {''})[0] if terms - {''} else 'unannotated')
        histogram = counter['exclusive_primary_consequence']
        histogram[primary] = histogram.get(primary, 0) + 1
        ref, alt = record['ref'], record['alt']
        kind = 'SNV' if len(ref) == len(alt) == 1 else 'insertion' if len(alt) > len(ref) else 'deletion' if len(ref) > len(alt) else 'MNV_or_complex'
        counter['variant_type_histogram'][kind] = counter['variant_type_histogram'].get(kind, 0) + 1
    try:
        with ExitStack() as stack:
            long = writer(stack, workdir / 'transcript_samples.tsv.gz', FIELDS)
            intersections = writer(stack, workdir / 'intersections.tsv.gz',
                ['variant_id', 'carrier_count', 'quality_passing_carriers', 'verified_reference_count', 'verified_reference_samples', 'unknown_or_filtered_count', 'unknown_or_filtered_samples',
                 'sample_states_json', 'verified_unique_to_one_sample', 'absence_scope'])
            for record in parse_vcf(annotated):
                counts['records'] += 1
                states = {s: record_state(record, d) for s, d in record['samples'].items()}
                carriers = sorted(s for s, state in states.items() if state == 'alternate')
                references = sorted(s for s, state in states.items() if state == 'reference')
                unknown = sorted(s for s, state in states.items() if state not in ('alternate', 'reference'))
                intersections.writerow(dict(variant_id=variant_key(record['chrom'], record['pos'], record['ref'], record['alt']),
                    quality_passing_carriers='&'.join(carriers), verified_reference_samples='&'.join(references),
                    carrier_count=len(carriers), verified_reference_count=len(references), unknown_or_filtered_count=len(unknown),
                    unknown_or_filtered_samples='&'.join(unknown), sample_states_json=jdump(states),
                    verified_unique_to_one_sample=int(len(carriers) == 1 and not unknown),
                    absence_scope='only_explicit_quality_passing_reference_GT; missing_records_not_reference'))
                counts['missing_csq'] += int(record['csq'] == [{}])
                terms = set().union(*(set(csq.get('Consequence', '').split('&')) for csq in record['csq']))
                candidates_for_record = set()
                for csq in record['csq']:
                    for sample, data in record['samples'].items():
                        row = annotation_row(record, csq, sample, data, context)
                        long.writerow(row)
                        if row['candidate_eligible']:
                            candidates_for_record.add(sample)
                        counts['transcript_sample_rows'] += 1
                        status = row['AlphaMissense_match_status']
                        counts['alphamissense_transcript_sample_statuses'][status] = counts['alphamissense_transcript_sample_statuses'].get(status, 0) + 1
                real = sequence_allele(record['ref'], record['alt'])
                counts['nonsequence_records_preserved'] += int(not real)
                memberships = {'raw': {s for s, data in record['samples'].items() if real and '1' in re.split(r'[/|]', data.get('GT', '.'))},
                               'qc': set(carriers), 'candidate': candidates_for_record}
                for sample in record['samples']:
                    counters['samples'].setdefault(sample, {level: empty_counter() for level in memberships})
                for level, members in memberships.items():
                    if members:
                        increment(counters['cohort'][level], terms, record)
                    for sample in members:
                        increment(counters['samples'][sample][level], terms, record)
    finally:
        db.close()
        # Index is a transient implementation detail; preserve only portable evidence tables.
        dbpath.unlink()
    write_json(workdir / 'consequence_counts.json', counters)
    return counts


def rank_key(row):
    descending = ['phenotype_context_rank', 'target_gene_context_rank', 'clinvar_review_rank', 'gene_validity_rank',
                  'consequence_rank', 'population_rank', 'predictors_damaging', 'AlphaMissense_research_rank']
    return (*[-int(row.get(field) or 0) for field in descending], -int(bool(row.get('mane_select'))),
            row['variant_id'], row.get('gene', ''), row.get('transcript', ''), row['sample'])


def annotate(cfg, workdir, runtime):
    config = cfg['annotation']
    input_vcf = Path(config.get('input_vcf') or runtime.stage_path('filter') / 'cohort.filtered.vcf.gz')
    manifest_path = Path(config['resource_manifest'])
    if not runtime.execute:
        return {'outputs': [], 'metadata': dict(status='DRY_RUN', input=str(input_vcf), resource_manifest=str(manifest_path), offline=True)}
    need(config.get('assembly', 'GRCh38') == 'GRCh38', 'Annotation assembly must match GRCh38.')
    manifest = load_resources(manifest_path, local_adapters=config.get('local_adapters'))
    if config.get('require_alphamissense', True):
        alpha_module().load_manifest(alpha_module().manifest_location(cfg))
    image = cfg['tools']['vep_image']
    need(image == manifest['vep']['image'] and str(config.get('cache_version', 115)) == str(manifest['vep']['cache_version']),
         'VEP image/cache differ from the locked prepared resources.')
    fasta = Path(cfg['reference']['fasta']).resolve(strict=True)
    need(sha256(fasta) == cfg['reference']['sha256'], 'Annotation reference differs from the analysis reference.')
    source_root = manifest_path.parent.resolve()
    need((source_root / 'cache/homo_sapiens/115_GRCh38/info.txt').is_file(), 'Prepared VEP cache missing.')
    output = workdir / 'annotated.vcf.gz'
    sequence_input = workdir / 'sequence_input.vcf.gz'
    sequence_count = partition_sequence_input(input_vcf, sequence_input)
    mounts = [dict(source=str(workdir.resolve()), target='/input', readonly=True),
              dict(source=str(source_root), target='/resources', readonly=True),
              dict(source=str(fasta.parent), target='/reference', readonly=True),
              dict(source=str(workdir.resolve()), target='/output', readonly=False)]
    arguments = ['vep', '--offline', '--cache', '--dir_cache', '/resources/cache', '--cache_version', '115',
        '--species', 'homo_sapiens', '--assembly', 'GRCh38', '--input_file', '/input/' + sequence_input.name,
        '--output_file', '/output/sequence.annotated.vcf.gz', '--format', 'vcf', '--vcf', '--compress_output', 'bgzip',
        '--fasta', '/reference/' + fasta.name, '--check_ref', '--hgvs', '--hgvsg', '--symbol', '--numbers',
        '--canonical', '--mane', '--transcript_version', '--protein', '--biotype', '--total_length',
        '--allele_number', '--sift', 'b', '--polyphen', 'b', '--check_existing', '--af_gnomade', '--af_gnomadg',
        '--domains', '--uniprot', '--tsl', '--appris', '--variant_class', '--pubmed', '--no_stats',
        '--fork', str(config.get('forks', 4))]
    clinvar = manifest['resources'].get('clinvar', {})
    if clinvar.get('status') == 'available':
        path = resource_path(manifest_path, clinvar)
        need(path.with_name(path.name + '.tbi').is_file(), 'ClinVar index missing.')
        arguments += ['--custom', 'file=/resources/' + path.relative_to(source_root).as_posix() +
            ',short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNSIGCONF%CLNDN%CLNDISDB%ALLELEID']
    argv = runtime.docker(image, arguments, mounts=mounts, gpu=False, cpus=config.get('forks', 4), memory_gib=config.get('memory_gib', 24))
    # docker --network none is mandatory, regardless of the runtime helper's defaults.
    if '--network' not in argv and not any(str(arg).startswith('--network=') for arg in argv):
        position = argv.index('run') + 1
        argv[position:position] = ['--network', 'none']
    else:
        need('--network=none' in argv or ('--network' in argv and argv[argv.index('--network') + 1] == 'none'), 'Annotation requires Docker network none.')
    if sequence_count:
        runtime.run(argv, logname='offline-vep', timeout=config.get('timeout_seconds', 21600))
    else:
        # Empty sequence partition still has the complete sample/header contract.
        import shutil
        shutil.copyfile(sequence_input, workdir / 'sequence.annotated.vcf.gz')
    merged = workdir / 'annotated.all-records.vcf'
    merge_vep_preserving_all(input_vcf, workdir / 'sequence.annotated.vcf.gz', merged)
    bcftools = cfg['tools'].get('bcftools', 'bcftools')
    runtime.run([bcftools, 'view', '-Oz', '-o', str(output), str(merged)], logname='compress-all-annotation-records')
    runtime.run([bcftools, 'index', '-t', str(output)], logname='index-annotated-vcf')
    # Remove only this stage's transient partitions after successful validated merge.
    for temporary in (sequence_input, workdir / 'sequence.annotated.vcf.gz', merged):
        temporary.unlink()
    need(output.is_file() and output.stat().st_size > 0, 'VEP did not produce annotated VCF.')
    alpha_inputs = alpha_module().prepare_matches(cfg, workdir, runtime, output) if config.get('require_alphamissense', True) else None
    counts = make_tables(output, workdir, config, manifest, manifest_path, alpha_inputs)
    if alpha_inputs:
        for paths in alpha_inputs[1].values():
            paths['regions'].unlink()
            paths['matches'].unlink()
    alpha_availability = (dict(status='available_local', manifest_sha256=sha256(alpha_inputs[2]),
        dataset_version=alpha_inputs[0]['dataset_version'], source_reference=alpha_inputs[0]['source_reference'],
        transcript_annotation=alpha_inputs[0]['transcript_annotation'], license=alpha_inputs[0]['license'],
        citation=alpha_inputs[0]['citation'], exact_transcript_version_and_protein_substitution_required=True)
        if alpha_inputs else dict(status='disabled_explicitly_in_this_context'))
    if alpha_inputs:
        manifest['resources']['alphamissense'] = dict(status='available_local_supplement',
            version=alpha_inputs[0]['dataset_version'], manifest_sha256=sha256(alpha_inputs[2]),
            source=alpha_inputs[0]['citation'], license=alpha_inputs[0]['license'],
            scope='exact_allele_transcript_version_and_protein_change; model_labels_not_ACMG')
        if manifest['resources'].get('AlphaMissense', {}).get('status', '').startswith('unavailable'):
            del manifest['resources']['AlphaMissense']
    write_json(workdir / 'resource_availability.json', dict(manifest_sha256=sha256(manifest_path),
        vep=manifest['vep'], cache_info=manifest['cache_info'], resources=manifest['resources'],
        alphamissense=alpha_availability,
        audited_source_coverage=requirements_coverage(config, manifest),
        no_patient_variant_remote_queries=True))
    write_json(workdir / 'annotation_counts.json', counts)
    return {'outputs': ['annotated.vcf.gz', 'annotated.vcf.gz.tbi', 'transcript_samples.tsv.gz', 'intersections.tsv.gz',
                         'resource_availability.json', 'annotation_counts.json', 'consequence_counts.json'],
            'metadata': dict(counts, offline=True, manual_classification_required=True)}


def local_pathways(candidate_genes, manifest, manifest_path, config):
    need(config.get('pathway_mode', 'disabled') in ('disabled', 'descriptive', 'enrichment'), 'Unknown explicit pathway mode.')
    if config.get('pathway_mode', 'disabled') == 'disabled':
        return dict(status='disabled_primary_transcript_withdrew_enrichment', results=[])
    item = manifest['resources'].get('pathways', {})
    if item.get('status') != 'available':
        return dict(status='unavailable_local_pathway_adapter', results=[])
    pathways = {}
    with text_open(resource_path(manifest_path, item)) as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        need({'gene', 'pathway'}.issubset(reader.fieldnames or []), 'Pathway adapter needs gene/pathway columns.')
        for row in reader:
            pathways.setdefault(row['pathway'], set()).add(row['gene'])
    background_path = config.get('pathway_background')
    mode = config.get('pathway_mode', 'descriptive')
    rows = [dict(pathway=name, observed_genes=sorted(genes & candidate_genes), count=len(genes & candidate_genes))
            for name, genes in sorted(pathways.items()) if genes & candidate_genes]
    if mode != 'enrichment':
        return dict(status='descriptive_only', source_version=item['version'], results=rows)
    if not background_path:
        return dict(status='missing_background', results=rows)
    background = {line.strip() for line in Path(background_path).read_text().splitlines() if line.strip() and not line.startswith('#')}
    if not background or not candidate_genes.issubset(background):
        return dict(status='invalid_or_incomplete_explicit_background', results=rows)
    N, n = len(background), len(candidate_genes)
    tests = []
    for name, genes in sorted(pathways.items()):
        K, k = len(genes & background), len(genes & candidate_genes)
        if K == 0:
            continue
        p = sum(math.exp(math.lgamma(K + 1) - math.lgamma(i + 1) - math.lgamma(K - i + 1)
            + math.lgamma(N - K + 1) - math.lgamma(n - i + 1) - math.lgamma(N - K - n + i + 1)
            - math.lgamma(N + 1) + math.lgamma(n + 1) + math.lgamma(N - n + 1))
            for i in range(max(k, n - (N - K)), min(K, n) + 1))
        tests.append(dict(pathway=name, p_value=min(1, p), observed_genes=sorted(genes & candidate_genes),
                          foreground=n, background=N, pathway_background=K, count=k))
    ordered = sorted(tests, key=lambda r: r['p_value'])
    adjusted = 1.0
    for i in range(len(ordered) - 1, -1, -1):
        adjusted = min(adjusted, ordered[i]['p_value'] * len(ordered) / (i + 1))
        ordered[i]['BH_q_value'] = adjusted
    return dict(status='research_only_explicit_background', background_sha256=sha256(background_path),
                source_version=item['version'], tests=len(tests), results=ordered)


def report(cfg, workdir, runtime):
    source = runtime.stage_path('annotate')
    if not runtime.execute:
        return {'outputs': [], 'metadata': dict(status='DRY_RUN', input=str(source), main_report_section='7.2')}
    availability = json.loads((source / 'resource_availability.json').read_text())
    config = cfg['annotation']
    consequence_counts = json.loads((source / 'consequence_counts.json').read_text())
    outputs, best, counts, sample_names = [], {}, {}, set(consequence_counts['samples'])
    gene_carriers, previous_variant, counted_gene_sample = {}, None, set()
    with ExitStack() as stack, text_open(source / 'transcript_samples.tsv.gz') as stream:
        supplementary, candidates, urgent = {}, {}, {}
        master = writer(stack, workdir / 'master_evidence.tsv.gz', FIELDS)
        outputs.append('master_evidence.tsv.gz')
        def ensure_sample(sample):
            need(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}', sample), 'Unsafe sample filename.')
            if sample not in supplementary:
                for name, collection in [('Supplementary_', supplementary), ('main_candidates_', candidates), ('urgent_manual_review_', urgent)]:
                    filename = name + sample + '.tsv.gz'
                    collection[sample] = writer(stack, workdir / filename, FIELDS)
                    outputs.append(filename)
        for sample in sorted(sample_names):
            ensure_sample(sample)
        for row in csv.DictReader(stream, delimiter='\t'):
            sample = row['sample']
            sample_names.add(sample)
            need(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}', sample), 'Unsafe sample filename.')
            master.writerow(row)
            if row['variant_id'] != previous_variant:
                previous_variant, counted_gene_sample = row['variant_id'], set()
            gene_sample = row['gene'], sample
            if row['gene'] and row['sample_state'] == 'alternate' and gene_sample not in counted_gene_sample:
                counted_gene_sample.add(gene_sample)
                carrier_counts = gene_carriers.setdefault(row['gene'], {})
                carrier_counts[sample] = carrier_counts.get(sample, 0) + 1
            counts[sample] = counts.get(sample, 0) + 1
            ensure_sample(sample)
            observed_alt = '1' in re.split(r'[/|]', row['GT'])
            if observed_alt:
                supplementary[sample].writerow(row)
            if row['external_pathogenic_assertion'] == '1' and observed_alt and row['sample_state'] != 'alternate':
                urgent[sample].writerow(row)
            if row['candidate_eligible'] == '1':
                candidates[sample].writerow(row)
                key = sample, row['variant_id']
                if key not in best or rank_key(row) < rank_key(best[key]):
                    best[key] = row
    all_ranked = sorted(best.values(), key=rank_key)
    with ExitStack() as stack:
        fields = ['rank', *FIELDS]
        for sample in sorted(sample_names):
            ranked = [r for r in all_ranked if r['sample'] == sample]
            for limit in (10, 20):
                filename = 'top' + str(limit) + '_' + sample + '.tsv'
                out = writer(stack, workdir / filename, fields)
                for index, row in enumerate(ranked[:limit], 1):
                    out.writerow(dict(row, rank=index))
                outputs.append(filename)
        cohort = {}
        for row in all_ranked:
            if row['variant_id'] not in cohort:
                cohort[row['variant_id']] = row
        cohort_ranked = sorted(cohort.values(), key=rank_key)
        for limit in (10, 20):
            filename = 'cohort_top' + str(limit) + '.tsv'
            out = writer(stack, workdir / filename, ['rank', 'carrier_samples', *FIELDS])
            for index, row in enumerate(cohort_ranked[:limit], 1):
                carriers = sorted({r['sample'] for r in all_ranked if r['variant_id'] == row['variant_id']})
                out.writerow(dict(row, rank=index, carrier_samples='&'.join(carriers)))
            outputs.append(filename)
        worksheet = writer(stack, workdir / 'ACMG_manual_worksheet.tsv.gz',
            ['sample', 'variant_id', 'gene', 'criterion', 'status', 'strength', 'evidence_source', 'assessor', 'date', 'rationale'])
        for row in all_ranked:
            for criterion in ACMG:
                worksheet.writerow(dict(sample=row['sample'], variant_id=row['variant_id'], gene=row['gene'], criterion=criterion,
                    status='not_assessed_manual_review_required'))
        outputs.append('ACMG_manual_worksheet.tsv.gz')
        gene_table = writer(stack, workdir / 'gene_intersections.tsv.gz',
            ['gene', 'quality_passing_carrier_samples', 'carrier_count', 'distinct_variant_counts_by_sample_json', 'absence_status'])
        for gene, carriers in sorted(gene_carriers.items()):
            gene_table.writerow(dict(gene=gene, quality_passing_carrier_samples='&'.join(sorted(carriers)),
                carrier_count=len(carriers), distinct_variant_counts_by_sample_json=jdump(carriers),
                absence_status='gene_negative_status_not_established; no_gene_wide_callability_claim'))
        outputs.append('gene_intersections.tsv.gz')
    # Stream-copy intersections; unknown/no-call remains distinct from true reference.
    with (source / 'intersections.tsv.gz').open('rb') as src, (workdir / 'intersections.tsv.gz').open('xb') as dst:
        import shutil
        shutil.copyfileobj(src, dst)
    outputs.append('intersections.tsv.gz')
    write_json(workdir / 'consequence_counts.json', consequence_counts)
    outputs.append('consequence_counts.json')
    manifest_path = Path(config['resource_manifest'])
    manifest = load_resources(manifest_path, verify_tree=False, local_adapters=config.get('local_adapters'))
    pathway = local_pathways({r['gene'] for r in all_ranked if r['gene']}, manifest, manifest_path, config)
    write_json(workdir / 'pathway_context.json', pathway)
    write_json(workdir / 'resource_availability.json', availability)
    outputs += ['pathway_context.json', 'resource_availability.json']
    limitations = [
        'Research prioritization is not a diagnosis or an ACMG classification; every criterion requires manual assessment and sources.',
        'stop_gain/frameshift/splice are pLoF candidates, not automatic LoF, NMD, PVS1, gain-of-function or pathogenicity.',
        'Protein beginning/middle/end divides transcript protein length into thirds; it is not a functional-domain or NMD rule.',
        'SIFT/PolyPhen/FATHMM/CADD vote is a correlated research heuristic, not four independent ACMG evidence items. Fourth predictor identity remains explicitly unconfirmed unless configured.',
        'AlphaMissense is a separate model prediction, not an ACMG class or an independent clinical evidence item. Its GENCODE v32 transcript versions must match the current transcript and protein substitution; mismatches retain source context without assigning a score to all transcripts.',
        'Missing predictors and absent ancestry AF are unavailable, not benign predictions or AF=0. Ancestry is not inferred from the patient.',
        'ClinVar values are external summary assertions; review/conflict/snapshot fields remain separate. Summary VCF does not supply all SCV-level assertions, dates or supporting evidence.',
        'ClinGen gene validity/dosage belongs to a gene-disease/inheritance context; dosage scores do not establish a specific small-variant mechanism.',
        'HPO matches, if supplied locally, are exact-term descriptive matches; absence of phenotype is not evidence against a gene.',
        'Only explicit quality-passing reference GT supports within-cohort absence. Unknown/no-call/filtered samples are not reference and cannot establish a private or de novo variant.',
        'No familial relationships, phase-based compound heterozygosity or de novo status is inferred.',
        'This small-variant workflow does not establish SV/CNV, repeats, mosaicism, mitochondrial or sex-chromosome completeness; evaluated intervals/ploidy must be read from the call plan.',
        'Coverage/callability metrics must be reviewed separately; VCF absence does not establish adequate coverage or a negative clinical result.',
        'Enrichment was withdrawn in the primary transcript and pathways are disabled by default. Explicitly requested local descriptions are not enrichment; inference requires an assay-appropriate background and reports all tested pathways with BH adjustment.',
        'Top10/top20 are deterministic review queues with no padding. Low-quality alternate ClinVar P/LP assertions appear separately in urgent manual review, never as validated candidates.'
    ]
    methods = dict(schema_version=1, section='7.2', analysis_scope=cfg.get('analysis_scope', 'Not supplied; inspect call intervals and ploidy'),
        reference=cfg.get('reference', 'not_supplied_in_this_report_context'), tools=cfg.get('tools', {}),
        callsets=cfg.get('callsets', []), disease_context=config.get('disease_context', 'not_supplied'),
        ranking=dict(rule='lexicographic_v1', order=['phenotype_exact_term_match', 'explicit_target_gene_context', 'external_ClinVar_review_rank_if_P_LP',
            'ClinGen_Strong_Definitive_gene_context', 'consequence_research_rank', 'observed_population_AF_rank',
            'research_damaging_votes', 'AlphaMissense_compatible_research_rank', 'MANE_select', 'variant_gene_transcript_sample_lexical_ties'],
            max_af=config.get('research_max_af', 0.01), cadd_threshold=config.get('cadd_research_threshold', 20),
            alphamissense_research_threshold=config.get('alphamissense_research_threshold', .564),
            candidate_rule='quality_passing_ALT AND (external_ClinVar_P_LP OR (functional_consequence AND (observed_AF<=threshold OR AF_missing_for_manual_review)))',
            missing_AF='eligible_for_manual_population_review; never labelled rare', no_padding=True),
        limitations=limitations, resources=availability, pathway_status=pathway['status'],
        phenotype_status='local_terms_supplied' if config.get('hpo_terms') or config.get('hpo_terms_by_sample') else 'not_supplied',
        target_gene_list_status='explicit_list_supplied' if config.get('target_genes') else 'not_supplied; no_previous_gene_list_inferred',
        target_gene_list_provenance=config.get('target_gene_list_provenance', 'not_supplied; target membership alone is not gene-disease evidence'),
        target_gene_list_sha256=hashlib.sha256(jdump(config.get('target_genes', [])).encode()).hexdigest(),
        disease_conclusion=None,
        manual_criteria_sources={'ACMG_AMP_2015': 'https://pubmed.ncbi.nlm.nih.gov/25741868/',
            'ClinGen_PVS1_2018': 'https://pubmed.ncbi.nlm.nih.gov/30192042/',
            'current_gene_specific_specifications': 'https://cspec.clinicalgenome.org/'},
        deprecated_criteria_warning='PP5/BP6 are legacy worksheet labels, not automatically applicable. Use current ClinGen/gene-specific specifications.')
    filter_path = runtime.stage_path('filter') / 'filtering-policy.json'
    methods['filtering'] = json.loads(filter_path.read_text()) if filter_path.is_file() else dict(status='filter_policy_not_present_in_this_report_context')
    methods['alphamissense'] = availability.get('alphamissense', dict(status='not_prepared_in_this_report_context'))
    annotation_counts_path = source / 'annotation_counts.json'
    if annotation_counts_path.is_file():
        annotations = json.loads(annotation_counts_path.read_text())
        alpha_coverage = dict(unit='transcript_sample_rows; not unique allele counts',
            source_rows=annotations.get('alphamissense_matched_source_rows'),
            matching_statuses=annotations.get('alphamissense_transcript_sample_statuses', {}))
        write_json(workdir / 'alphamissense_coverage.json', alpha_coverage)
        outputs.append('alphamissense_coverage.json')
        methods['alphamissense_coverage'] = alpha_coverage
    methods['stage_receipts'] = {}
    for stage in ('call', 'joint', 'filter', 'annotate'):
        receipt_path = runtime.stage_path(stage) / 'workflow-receipt.json'
        methods['stage_receipts'][stage] = (dict(path=str(receipt_path), sha256=sha256(receipt_path))
            if receipt_path.is_file() else dict(status='not_present_in_this_report_context'))
    outputs.append('methods_and_limitations.json')
    plot_images = []
    report_python = cfg.get('tools', {}).get('report_python') or str(Path(cfg.get('data_root', '/srv/marina')) / 'tools/report-venv/bin/python')
    if config.get('make_plots', True):
        plot_script = Path(__file__).resolve().parent / 'plot_report.py'
        runtime.run([str(report_python), str(plot_script), '--counts', str(workdir / 'consequence_counts.json'),
                     '--output-dir', str(workdir / 'figures')], logname='count-histograms', timeout=300)
        plot_manifest = json.loads((workdir / 'figures/plot_manifest.json').read_text())
        outputs.append('figures/plot_manifest.json')
        for plot in plot_manifest['plots']:
            outputs.extend(['figures/' + plot['svg'], 'figures/' + plot['png']])
            plot_images.append('figures/' + plot['svg'])
        methods['plotting'] = dict(matplotlib_version=plot_manifest['matplotlib_version'],
            count_file_sha256=plot_manifest['count_file_sha256'], manifest='figures/plot_manifest.json')
    write_json(workdir / 'methods_and_limitations.json', methods)
    if config.get('export_docx', True):
        outputs.append('research_report.docx')
    links = ''.join('<li><a href="' + html.escape(p, quote=True) + '">' + html.escape(p) + '</a></li>' for p in outputs)
    body = '<!doctype html><meta charset="utf-8"><title>Evidence report</title><h1>7.2 Review priorities and full evidence</h1>'
    body += '<p>Research review queues; no automatic disease conclusion or ACMG classification. Tables are not padded.</p><ul>' + links + '</ul>'
    body += '<h2>Methods and limitations</h2><ul>' + ''.join('<li>' + html.escape(s) + '</li>' for s in limitations) + '</ul>'
    if plot_images:
        body += '<h2>Histograms: raw, QC and research candidates</h2><p>Counts use unique normalized alleles. Flags overlap; primary consequences and variant types are separate exclusive histograms.</p>'
        body += ''.join('<figure><img style="max-width:100%" src="' + html.escape(path, quote=True) + '" alt="Unique allele histogram"><figcaption>' + html.escape(path) + '</figcaption></figure>' for path in plot_images)
    (workdir / 'report.html').write_text(body, encoding='utf-8')
    outputs.append('report.html')
    if config.get('export_docx', True):
        docx_command = [str(report_python), str(Path(__file__).resolve().parent / 'export_report_docx.py'),
            '--source', str(workdir), '--output', str(workdir / 'research_report.docx')]
        if cfg.get('synthetic_validation') is True:
            docx_command.append('--synthetic')
        runtime.run(docx_command, logname='export-research-report-docx', timeout=300)
        need((workdir / 'research_report.docx').is_file(), 'DOCX exporter did not produce its declared output.')
    return {'outputs': outputs, 'metadata': dict(samples=len(sample_names), eligible_sample_variants=len(best),
        cohort_candidate_variants=len(cohort), report_section='7.2', automatic_classification=False)}


def run_stage(stage, cfg, workdir, runtime):
    workdir = Path(workdir)
    need(stage in ('annotate', 'report'), 'Unknown annotation/report stage.')
    if runtime.execute:
        workdir.mkdir(parents=True, exist_ok=True)
    return annotate(cfg, workdir, runtime) if stage == 'annotate' else report(cfg, workdir, runtime)
