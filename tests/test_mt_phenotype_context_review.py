import csv
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('mt_context', Path(__file__).resolve().parents[1] / 'scripts/mt_phenotype_context_review.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def json_file(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj), encoding='utf-8')


def gzip_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


class ContextTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        functional = self.root / 'Reports/mitochondrial-functional-views-v1'
        consequence = self.root / 'Reports/mitochondrial-consequence-annotation-v1'
        table = functional / 'tables/cohort/site_PASS/all.functional.tsv.gz'
        record = {'ALT_vector': 'G,T', 'ALT_index_1': 2, 'GT': '0/1/2', 'AD_vector': '0,53,1', 'AF_vector': '0.945,0.037', 'HL_raw': '.', 'DP_raw': '54', 'FILTER': 'PASS', 'FORMAT': 'GT:AD:AF:DP', 'raw_sample': '0/1/2:0,53,1:0.945,0.037:54'}
        frequency = {'CHROM': 'chrM', 'POS': '100', 'REF': 'A', 'ALT': 'G,T', 'gnomAD_MT_match_status': 'MULTIALLELIC_OR_SYMBOLIC_NO_FLATTEN_NO_JOIN', 'gnomAD_MT_reference_json': '', 'Helix_match_status': 'MULTIALLELIC_OR_SYMBOLIC_NO_FLATTEN_NO_JOIN', 'Helix_reference_json': ''}
        row = {'allele_id': 'chrM:100:A:T', 'chrom': 'chrM', 'pos': '100', 'ref': 'A', 'alt': 'T', 'observed_call_count': '1', 'observed_call_samples': 'PUBLIC_SYNTHETIC', 'site_PASS_count': '1', 'site_PASS_samples': 'PUBLIC_SYNTHETIC', 'sample_source_context_json': json.dumps({'PUBLIC_SYNTHETIC': {'state': 'OBSERVED_MUTECT2_CALL_RECORD', 'original_records': [record]}}), 'genes_json': json.dumps(['MT-ND4L', 'MT-ND5']), 'frequency_contexts_json': json.dumps({'PUBLIC_SYNTHETIC': [json.dumps(frequency)]})}
        gzip_tsv(table, [row])
        master = {'sample': 'PUBLIC_SYNTHETIC', 'chrom': 'chrM', 'pos': '100', 'ref': 'A', 'alt': 'G,T', 'matched_alt': 'T', 'allele_number': '2', 'GT': record['GT'], 'AD': record['AD_vector'], 'AF': record['AF_vector'], 'HL': '.', 'DP': '54', 'FILTER': 'PASS', 'source_FORMAT': record['FORMAT'], 'source_sample': record['raw_sample'], 'gene': 'MT-ND5', 'transcript': 'PUBLIC_NEIGHBOR.1', 'consequence': 'upstream_gene_variant', 'hgvs_c': '', 'hgvs_p': '', 'csq_json': '{}'}
        direct = {**master, 'gene': 'MT-ND4L', 'transcript': 'PUBLIC_DIRECT.1', 'consequence': 'missense_variant', 'hgvs_p': 'PUBLIC_PROTEIN:p.Met27Ile'}
        gzip_tsv(consequence / 'PUBLIC_SYNTHETIC/master.from-vcf.tsv.gz', [master, direct])
        json_file(functional / 'manifest.json', {'status': 'PASS_MT_FUNCTIONAL_VIEWS', 'samples': ['PUBLIC_SYNTHETIC'], 'reviewed_profile_sha256': '0' * 64, 'artifacts': {'tables/cohort/site_PASS/all.functional.tsv.gz': {'bytes': table.stat().st_size, 'sha256': '0' * 64}}})
        json_file(functional / 'workflow-receipt.json', {'status': 'PASS'})
        json_file(consequence / 'manifest.json', {'status': 'PASS_MT_CONSEQUENCE_ANNOTATION', 'reviewed_profile': {'sha256': '0' * 64}})
        self.hpo = self.root / 'public-hpo.tsv'
        self.hpo.write_text('ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n1\tMT-ND5\tHP:0001639\tHypertrophic cardiomyopathy\t-\tPUBLIC:1\n2\tMT-ND4L\tHP:0000001\tOther term\t-\tPUBLIC:2\n', encoding='utf-8')

    def review(self):
        return MODULE.build_review(self.root, self.hpo)

    def test_multiallelic_exact_T_fraction_not_maximum(self):
        result = self.review()
        context = result['alleles'][0]['sample_original_AF_HL_and_call_context']['PUBLIC_SYNTHETIC']
        self.assertEqual(context['exact_ALT_fractions'][0]['AF']['exact_allele_raw'], '0.037')
        self.assertEqual(context['original_records'][0]['AF_vector'], '0.945,0.037')
        self.assertEqual(context['exact_ALT_fractions'][0]['HL']['status'], 'MISSING_NOT_ZERO')

    def test_nearby_only_HPO_gene_is_excluded(self):
        allele = self.review()['alleles'][0]
        self.assertEqual(allele['direct_genes'], ['MT-ND4L'])
        self.assertEqual(allele['exact_HPO_HP0001639_direct_genes'], [])
        self.assertEqual(allele['excluded_nearby_only_HPO_genes'], ['MT-ND5'])

    def test_missing_population_context_is_not_zero(self):
        result = self.review()
        allele = result['alleles'][0]
        self.assertEqual(result['summary']['no_existing_exact_population_AF_context_count'], 1)
        self.assertTrue(all(x['source_AF_fields_raw'] is None for x in allele['existing_frequency_contexts']))
        self.assertEqual(allele['descriptive_source_AF_above_0_01'], [])

    def test_export_full_rows_and_no_overwrite(self):
        result = self.review()
        output = self.root / 'new-review'
        MODULE.write_outputs(result, output)
        with (output / 'sample-allele-fractions.tsv').open(encoding='utf-8', newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(rows[0]['AF_exact_ALT_raw'], '0.037')
        self.assertEqual(len((output / 'alleles.tsv').read_text().splitlines()), 2)
        self.assertEqual(json.loads((output / 'review.json').read_text())['summary'], result['summary'])
        with self.assertRaises(FileExistsError):
            MODULE.write_outputs(result, output)

    def test_bad_AF_vector_is_explicitly_unavailable(self):
        value = MODULE.matched_fraction({'AF_vector': '0.9', 'ALT_vector': 'G,T', 'ALT_index_1': 2}, 'AF_vector')
        self.assertIsNone(value['exact_allele_raw'])
        self.assertEqual(value['status'], 'AMBIGUOUS_OR_INVALID_VECTOR_NOT_FLATTENED')


if __name__ == '__main__':
    unittest.main()
