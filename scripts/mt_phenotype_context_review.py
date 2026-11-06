"""Local descriptive mitochondrial phenotype context; no network or clinical scoring."""
import argparse
import collections
import csv
import datetime
import gzip
import json
import math
from pathlib import Path

NEARBY = {'upstream_gene_variant', 'downstream_gene_variant'}

def read_tsv(path):
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')

def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

def value_or_none(value):
    return None if value in ('', '.', None) else value

def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) and 0 <= result <= 1 else None
    except (TypeError, ValueError):
        return None

def original_signature(row):
    return (row['alt'], row['allele_number'], row['GT'], row['AD'], row['AF'], row['HL'], row['DP'], row['FILTER'], row['source_FORMAT'], row['source_sample'])

def table_signature(record):
    return (record['ALT_vector'], str(record['ALT_index_1']), record['GT'], record['AD_vector'], record['AF_vector'], record['HL_raw'], record['DP_raw'], record['FILTER'], record['FORMAT'], record['raw_sample'])

def matched_fraction(record, field):
    raw = record[field]
    values = raw.split(',')
    alternatives = record['ALT_vector'].split(',')
    index = record['ALT_index_1'] - 1
    if raw in ('', '.'):
        return {'raw': raw, 'exact_allele_raw': None, 'status': 'MISSING_NOT_ZERO'}
    if len(values) != len(alternatives) or not 0 <= index < len(values) or number(values[index]) is None:
        return {'raw': raw, 'exact_allele_raw': None, 'status': 'AMBIGUOUS_OR_INVALID_VECTOR_NOT_FLATTENED'}
    return {'raw': raw, 'exact_allele_raw': values[index], 'status': 'EXACT_ALT_INDEX_NUMBER_A'}

def build_review(results_root, hpo_genes):
    hpo_term = 'HP:0001639'
    RESULTS = Path(results_root).resolve()
    FUNCTIONAL = RESULTS / 'Reports/mitochondrial-functional-views-v1'
    CONSEQUENCE = RESULTS / 'Reports/mitochondrial-consequence-annotation-v1'
    manifest = load(FUNCTIONAL / 'manifest.json')
    source_manifest = load(CONSEQUENCE / 'manifest.json')
    receipt = load(FUNCTIONAL / 'workflow-receipt.json')
    assert manifest['status'] == 'PASS_MT_FUNCTIONAL_VIEWS' and receipt['status'] == 'PASS'
    assert source_manifest['status'] == 'PASS_MT_CONSEQUENCE_ANNOTATION'
    samples = manifest['samples']
    pass_path = FUNCTIONAL / 'tables/cohort/site_PASS/all.functional.tsv.gz'
    source_rows = list(read_tsv(pass_path))
    assert len({r['allele_id'] for r in source_rows}) == len(source_rows), 'Duplicate exact allele rows'
    pass_by_key = {(r['chrom'], r['pos'], r['ref'], r['alt']): r for r in source_rows}
    hpo_path = Path(hpo_genes).resolve()
    with hpo_path.open(encoding='utf-8', newline='') as handle:
        hpo_rows = [r for r in csv.DictReader(handle, delimiter='\t') if r['hpo_id'] == hpo_term]
    hpo_genes = {r['gene_symbol'] for r in hpo_rows}
    if not hpo_genes:
        raise ValueError('Exact phenotype term has no gene associations in source')
    hpo_by_gene = collections.defaultdict(list)
    for row in hpo_rows:
        hpo_by_gene[row['gene_symbol']].append(row)
    effects = collections.defaultdict(dict)
    record_signatures = collections.defaultdict(set)
    master_counts = {}
    master_paths = {}
    for sample in samples:
        path = CONSEQUENCE / sample / 'master.from-vcf.tsv.gz'
        master_paths[sample] = str(path)
        counts = collections.Counter()
        for row in read_tsv(path):
            counts['transcript_rows_all_filters'] += 1
            key = (row['chrom'], row['pos'], row['ref'], row['matched_alt'])
            if key not in pass_by_key:
                continue
            counts['transcript_rows_for_selected_PASS_alleles'] += 1
            record_signatures[(key, sample)].add(original_signature(row))
            consequences = sorted(set(row['consequence'].split('&')) - NEARBY - {''})
            if not consequences:
                counts['nearby_only_rows_excluded_from_gene_HCM_match'] += 1
                continue
            csq = json.loads(row['csq_json'])
            effect = {'gene': row['gene'], 'transcript': row['transcript'], 'direct_consequences': consequences,
                      'hgvs_c_raw': row['hgvs_c'], 'hgvs_p_raw': row['hgvs_p'],
                      'biotype': csq.get('BIOTYPE'), 'hgnc_id': csq.get('HGNC_ID'),
                      'hpo_exact_gene_match': row['gene'] in hpo_genes,
                      'existing_ClinVar_source_fields': {k: v for k, v in csq.items() if k.startswith('ClinVar')},
                      'samples_with_this_transcript_effect': [], 'site_PASS_samples_with_this_effect': []}
            # The two sample lists are excluded from the effect identity.
            identity = json.dumps({k: v for k, v in effect.items() if k not in ('samples_with_this_transcript_effect', 'site_PASS_samples_with_this_effect')}, sort_keys=True)
            target = effects[key].setdefault(identity, effect)
            if sample not in target['samples_with_this_transcript_effect']:
                target['samples_with_this_transcript_effect'].append(sample)
            if row['FILTER'] == 'PASS' and sample not in target['site_PASS_samples_with_this_effect']:
                target['site_PASS_samples_with_this_effect'].append(sample)
        master_counts[sample] = dict(counts)
    alleles = []
    signature_checks = 0
    for key, source in sorted(pass_by_key.items(), key=lambda item: (int(item[0][1]), item[0][2], item[0][3])):
        context = json.loads(source['sample_source_context_json'])
        sample_values = {}
        for sample in samples:
            observed = context.get(sample, {'state': 'NOT_OBSERVED_MUTECT2_CALL_RECORD', 'original_records': []})
            originals = observed.get('original_records', [])
            if originals:
                assert {table_signature(r) for r in originals} == record_signatures[(key, sample)], (key, sample)
                signature_checks += 1
            else:
                assert not record_signatures[(key, sample)], (key, sample)
            sample_values[sample] = {'state': observed['state'], 'literal_site_PASS': sample in source['site_PASS_samples'].split('&'),
                                     'original_records': originals,
                                     'exact_ALT_fractions': [{'AF': matched_fraction(r, 'AF_vector'), 'HL': matched_fraction(r, 'HL_raw')} for r in originals]}
        direct = list(effects[key].values())
        direct_genes = sorted({r['gene'] for r in direct if r['gene']})
        direct_consequences = sorted({c for r in direct for c in r['direct_consequences']})
        all_genes = json.loads(source['genes_json'])
        nearby_only_genes = sorted(set(all_genes) - set(direct_genes))
        frequencies = {}
        for sample, strings in json.loads(source['frequency_contexts_json']).items():
            for text in strings:
                context_row = json.loads(text)
                assert (context_row['CHROM'], context_row['POS'], context_row['REF']) == key[:3]
                for dataset, status_key, raw_key in [('gnomAD_MT', 'gnomAD_MT_match_status', 'gnomAD_MT_reference_json'), ('HelixMTdb', 'Helix_match_status', 'Helix_reference_json')]:
                    raw = context_row[raw_key]
                    status = context_row[status_key]
                    identity = json.dumps([dataset, status, raw])
                    ref = json.loads(raw) if raw else None
                    if ref is not None:
                        assert (str(ref['pos']), ref['ref'], ref['alt']) == key[1:]
                    entry = frequencies.setdefault(identity, {'dataset': dataset, 'match_status': status, 'reference_raw_json': raw,
                                                               'source_AF_fields_raw': {k: ref.get(k) for k in ('AF_hom', 'AF_het')} if ref else None,
                                                               'samples_with_this_context': []})
                    if sample not in entry['samples_with_this_context']:
                        entry['samples_with_this_context'].append(sample)
        descriptive_common = []
        for entry in frequencies.values():
            for field, raw in (entry['source_AF_fields_raw'] or {}).items():
                value = number(raw)
                if value is not None and value > 0.01:
                    descriptive_common.append({'dataset': entry['dataset'], 'field': field, 'raw_value': raw})
        alleles.append({'allele_id': source['allele_id'], 'chrom': key[0], 'pos': int(key[1]), 'ref': key[2], 'alt': key[3],
                        'observed_call_count': int(source['observed_call_count']), 'observed_call_samples': source['observed_call_samples'].split('&'),
                        'literal_site_PASS_count': int(source['site_PASS_count']), 'literal_site_PASS_samples': source['site_PASS_samples'].split('&'),
                        'direct_consequences': direct_consequences, 'direct_genes': direct_genes, 'direct_transcript_effects': direct,
                        'excluded_nearby_only_genes': nearby_only_genes,
                        'exact_HPO_HP0001639_direct_genes': sorted(set(direct_genes) & hpo_genes),
                        'excluded_nearby_only_HPO_genes': sorted(set(nearby_only_genes) & hpo_genes),
                        'sample_original_AF_HL_and_call_context': sample_values,
                        'existing_frequency_contexts': list(frequencies.values()),
                        'descriptive_source_AF_above_0_01': descriptive_common,
                        'frequency_missing_is_not_zero': True,
                        'clinical_classification_created': False})
    assert all(a['direct_consequences'] for a in alleles)
    consequences = collections.Counter(c for a in alleles for c in a['direct_consequences'])
    missense = [a for a in alleles if 'missense_variant' in a['direct_consequences']]
    common = [a for a in alleles if a['descriptive_source_AF_above_0_01']]
    no_frequency = [a for a in alleles if not any(x['source_AF_fields_raw'] for x in a['existing_frequency_contexts'])]
    all4 = [a['allele_id'] for a in alleles if a['literal_site_PASS_count'] == 4]
    result = {'schema_version': 1, 'status': 'PASS_LOCAL_DESCRIPTIVE_MT_HCM_AUDIT',
              'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'source_provenance': {'functional_manifest_path': str(FUNCTIONAL / 'manifest.json'),
                                    'functional_reviewed_profile_sha256': manifest['reviewed_profile_sha256'],
                                    'selected_table': {'path': str(pass_path), **manifest['artifacts']['tables/cohort/site_PASS/all.functional.tsv.gz']},
                                    'consequence_manifest_path': str(CONSEQUENCE / 'manifest.json'),
                                    'consequence_reviewed_profile': source_manifest['reviewed_profile'],
                                    'consequence_masters': master_paths,
                                    'HPO_file': {'local_path': str(hpo_path), 'bytes': hpo_path.stat().st_size}},
              'HPO_contract': {'term': hpo_term, 'name': hpo_rows[0]['hpo_name'], 'exact_term_rows': len(hpo_rows),
                               'unique_gene_symbols': len(hpo_genes), 'direct_gene_exact_symbol_matching_only': True,
                               'nearest_upstream_downstream_genes_excluded': True,
                               'gene_disease_context_is_not_variant_causality': True,
                               'matched_MT_gene_source_rows': {g: hpo_by_gene[g] for g in sorted({g for a in alleles for g in a['exact_HPO_HP0001639_direct_genes']})}},
              'summary': {'unique_literal_PASS_alleles': len(alleles), 'literal_PASS_sample_allele_observations': sum(a['literal_site_PASS_count'] for a in alleles),
                          'PASS_alleles_per_sample': {s: sum(a['sample_original_AF_HL_and_call_context'][s]['literal_site_PASS'] for a in alleles) for s in samples},
                          'distinct_direct_consequences_allele_counts': dict(consequences),
                          'direct_HPO_gene_matched_alleles': sum(bool(a['exact_HPO_HP0001639_direct_genes']) for a in alleles),
                          'nearest_only_HPO_match_without_direct_match_alleles': sum(bool(a['excluded_nearby_only_HPO_genes']) and not a['exact_HPO_HP0001639_direct_genes'] for a in alleles),
                          'descriptive_any_existing_population_AF_above_0_01_count': len(common),
                          'no_existing_exact_population_AF_context_count': len(no_frequency),
                          'all_four_literal_PASS_exact_alleles': all4,
                          'master_original_AF_HL_GT_AD_FILTER_record_vectors_crosschecked': signature_checks,
                          'all_original_AF_HL_vectors_unchanged': True, 'no_new_clinical_classification': True},
              'top_functional_observations_descriptive_not_causal': [{'allele_id': a['allele_id'], 'direct_genes': a['direct_genes'],
                   'direct_HPO_genes': a['exact_HPO_HP0001639_direct_genes'], 'literal_PASS_samples': a['literal_site_PASS_samples'],
                   'protein_changes_raw': sorted({e['hgvs_p_raw'] for e in a['direct_transcript_effects']}),
                   'descriptive_population_AF_above_0_01': a['descriptive_source_AF_above_0_01'],
                   'sample_fractions': {s: c['exact_ALT_fractions'] for s, c in a['sample_original_AF_HL_and_call_context'].items() if c['original_records']}}
                  for a in missense],
              'alleles': alleles, 'source_master_read_counts': master_counts,
              'limits': ['PASS denotes source site-filter status, not disease causality.',
                         'Only direct overlapping transcript effects supply HPO gene context; nearest upstream/downstream genes supply no HCM boost.',
                         'HP:0001639 is broad gene-level disease annotation; this audit performs no variant pathogenicity or individual penetrance inference.',
                         'AF_hom/AF_het above 0.01 is descriptive source frequency, not an ACMG threshold or benign classification.',
                         'Missing HL remains missing. AF is retained as measured caller allele fraction and is not replaced by diploid GT or a tissue-independent heteroplasmy estimate.',
                         'Multiallelic source population no-join statuses are preserved; no external requests, indel normalization, or missing-to-zero conversion.',
                         'No input was modified; existing frequency/reference hashes were retained without rehashing source payloads.']}
    return result


def write_outputs(result, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError('New output directory required; source files never overwritten')
    output.mkdir(parents=True)
    (output / 'review.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    fields = ['allele_id', 'chrom', 'pos', 'ref', 'alt', 'literal_site_PASS_count', 'literal_site_PASS_samples',
              'direct_genes', 'direct_consequences', 'exact_HPO_HP0001639_direct_genes', 'excluded_nearby_only_genes',
              'excluded_nearby_only_HPO_genes', 'sample_original_AF_HL_and_call_context', 'existing_frequency_contexts',
              'descriptive_source_AF_above_0_01', 'direct_transcript_effects']
    with (output / 'alleles.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for allele in result['alleles']:
            writer.writerow({k: json.dumps(allele[k], ensure_ascii=False, sort_keys=True) if isinstance(allele[k], (dict, list)) else allele[k] for k in fields})
    with (output / 'sample-allele-fractions.tsv').open('w', encoding='utf-8', newline='') as handle:
        fields = ['sample', 'allele_id', 'state', 'literal_site_PASS', 'ALT_vector', 'ALT_index_1', 'GT', 'AD_vector', 'DP_raw', 'FILTER', 'AF_raw', 'AF_exact_ALT_raw', 'AF_status', 'HL_raw', 'HL_exact_ALT_raw', 'HL_status']
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for allele in result['alleles']:
            for sample, context in allele['sample_original_AF_HL_and_call_context'].items():
                records = context['original_records']
                if not records:
                    writer.writerow({'sample': sample, 'allele_id': allele['allele_id'], 'state': context['state'], 'literal_site_PASS': False, 'AF_status': 'NOT_OBSERVED_NOT_REFERENCE', 'HL_status': 'NOT_OBSERVED_NOT_REFERENCE'})
                for record, fractions in zip(records, context['exact_ALT_fractions']):
                    row = {'sample': sample, 'allele_id': allele['allele_id'], 'state': context['state'], 'literal_site_PASS': context['literal_site_PASS']}
                    row.update({k: record[k] for k in ['ALT_vector', 'ALT_index_1', 'GT', 'AD_vector', 'DP_raw', 'FILTER']})
                    for label in ['AF', 'HL']:
                        row.update({label + '_raw': fractions[label]['raw'], label + '_exact_ALT_raw': fractions[label]['exact_allele_raw'], label + '_status': fractions[label]['status']})
                    writer.writerow(row)
    (output / 'summary.json').write_text(json.dumps(result['summary'], indent=2) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', required=True)
    parser.add_argument('--hpo-genes', required=True, help='Local genes_to_phenotype TSV')
    parser.add_argument('--output', required=True, help='New, absent output directory')
    args = parser.parse_args()
    result = build_review(args.results_root, args.hpo_genes)
    write_outputs(result, args.output)
    print(json.dumps({'status': result['status'], 'summary': result['summary'], 'output': str(Path(args.output).resolve())}))


if __name__ == '__main__':
    main()
