"""Synthetic-only checks: no patient data, network or annotation container required."""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import sqlite3
import tarfile
import io
import unittest
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / 'scripts/annotation_stages.py'
spec = importlib.util.spec_from_file_location('annotation_stages', MODULE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
common_spec = importlib.util.spec_from_file_location('workflow_common', MODULE.parent / 'workflow_common.py')
common = importlib.util.module_from_spec(common_spec)
common_spec.loader.exec_module(common)
prep_spec = importlib.util.spec_from_file_location('prepare_annotation_resources', MODULE.parent / 'prepare_annotation_resources.py')
prep = importlib.util.module_from_spec(prep_spec)
prep_spec.loader.exec_module(prep)


class AnnotationTests(unittest.TestCase):
    def test_missing_and_filtered_genotypes_are_not_reference(self):
        self.assertEqual(m.sample_state('0/0', 'PASS'), 'reference')
        self.assertEqual(m.sample_state('0/1', 'PASS'), 'alternate')
        self.assertEqual(m.sample_state('0/.', 'PASS'), 'unknown_no_call')
        self.assertEqual(m.sample_state('0/0', '.'), 'unknown_quality_reference')
        self.assertEqual(m.sample_state('0/1', 'LowGQ'), 'filtered_alt')
        self.assertEqual(m.sample_state('0/0', 'PASS', 'LowQUAL'), 'filtered_reference')
        with self.assertRaises(ValueError):
            m.sample_state('0/2', 'PASS')

    def test_predictors_missing_is_not_a_negative_vote(self):
        self.assertEqual(m.predictor_vote('deleterious(0.01)', 'probably_damaging(0.99)', '', ''),
                         (2, 2, 'incomplete_not_negative'))
        self.assertEqual(m.predictor_vote('deleterious(0.01)', 'probably_damaging(0.99)', 'D', '25'),
                         (4, 4, 'yes_research_only'))
        self.assertEqual(m.predictor_vote('tolerated(0.8)', 'benign(0.01)', 'T', '1'), (4, 0, 'no'))

    def test_protein_region_is_transcript_relative_not_a_mechanism(self):
        self.assertEqual(m.protein_region('10/100'), ('10', '100', 'beginning'))
        self.assertEqual(m.protein_region('50/100'), ('50', '100', 'middle'))
        self.assertEqual(m.protein_region('90/100'), ('90', '100', 'end'))
        self.assertEqual(m.protein_region('10')[2], 'unknown')
        self.assertEqual(m.protein_region('101/100')[2], 'unknown')

    def test_all_records_and_unknowns_survive_report_top_tables_not_padded(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            annotation = root / 'stages/06_annotate'
            report = root / 'stages/07_report'
            annotation.mkdir(parents=True)
            report.mkdir(parents=True)
            manifest = dict(schema_version=1, status='PASS', assembly='GRCh38', resources={},
                            vep={'version': 'synthetic'}, cache_info='synthetic fixture only', artifacts=[])
            manifest_path = root / 'resource-manifest.json'
            manifest_path.write_text(json.dumps(manifest))
            csq_fields = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature',
                          'Protein_position', 'SIFT', 'PolyPhen', 'ClinVar_CLNSIG', 'ClinVar_CLNREVSTAT', 'ALLELE_NUM']
            csq = '|'.join(['A', 'stop_gained', 'HIGH', 'GENE_A', 'ENSG_TEST', 'ENST_TEST',
                            '10/100', '', '', 'Pathogenic', 'reviewed_by_expert_panel', '1'])
            header = ['##fileformat=VCFv4.2',
                      '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: ' + '|'.join(csq_fields) + '">',
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B\tSAMPLE_C\tSAMPLE_D']
            vcf = annotation / 'annotated.vcf.gz'
            with gzip.open(vcf, 'wt') as stream:
                stream.write('\n'.join(header) + '\n')
                stream.write('chr1\t10\t.\tC\tA\t80\tPASS\tCSQ=' + csq + ',' + csq.replace('ENST_TEST', 'ENST_OTHER').replace('stop_gained', 'missense_variant') + '\tGT:FT:DP:AD:GQ\t0/1:PASS:20:8,12:99\t0/0:PASS:20:20,0:99\t./.:PASS:0:0,0:0\t0/1:LowGQ:10:5,5:2\n')
            counts = m.make_tables(vcf, annotation, {}, manifest, manifest_path)
            self.assertEqual(counts['transcript_sample_rows'], 8)
            counters = json.loads((annotation / 'consequence_counts.json').read_text())
            self.assertEqual(counters['cohort']['raw']['unique_normalized_alleles'], 1)
            self.assertEqual(counters['cohort']['raw']['nonexclusive_flags']['stop_gained'], 1)
            self.assertEqual(counters['cohort']['raw']['nonexclusive_flags']['missense'], 1)
            self.assertEqual(counters['samples']['SAMPLE_D']['raw']['unique_normalized_alleles'], 1)
            self.assertEqual(counters['samples']['SAMPLE_D']['qc']['unique_normalized_alleles'], 0)
            self.assertEqual(counters['cohort']['candidate']['exclusive_primary_consequence'], {'stop_gained': 1})
            with gzip.open(annotation / 'intersections.tsv.gz', 'rt') as stream:
                row = next(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual(row['quality_passing_carriers'], 'SAMPLE_A')
            self.assertEqual(row['verified_reference_samples'], 'SAMPLE_B')
            self.assertEqual(row['unknown_or_filtered_samples'], 'SAMPLE_C&SAMPLE_D')
            self.assertEqual(row['verified_unique_to_one_sample'], '0')
            availability = dict(resources={}, cache_info='synthetic', vep={'version': 'synthetic'})
            (annotation / 'resource_availability.json').write_text(json.dumps(availability))
            cfg = {'output_root': str(root), 'annotation': {'resource_manifest': str(manifest_path), 'make_plots': False, 'export_docx': False}}
            runtime = common.Runtime(cfg, execute=True, workdir=report, stage='report')
            result = m.run_stage('report', cfg, report, runtime)
            self.assertIn('report.html', result['outputs'])
            self.assertEqual(result['metadata']['cohort_candidate_variants'], 1)
            with (report / 'top20_SAMPLE_A.tsv').open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream, delimiter='\t'))), 1)
            with (report / 'top10_SAMPLE_D.tsv').open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream, delimiter='\t'))), 0)
            with gzip.open(report / 'urgent_manual_review_SAMPLE_D.tsv.gz', 'rt') as stream:
                urgent = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual(len(urgent), 2)
            self.assertEqual(urgent[0]['sample_state'], 'filtered_alt')
            with gzip.open(report / 'master_evidence.tsv.gz', 'rt') as stream:
                master = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual(len(master), 8)
            self.assertEqual(master[0]['DP'], '20')
            self.assertEqual(master[0]['GQ'], '99')
            self.assertEqual(master[0]['AD_derived_VAF'], '0.6')
            self.assertEqual(master[0]['QUAL'], '80')
            self.assertTrue(all(r['classification_status'].startswith('manual_ACMG') for r in master))
            methods = json.loads((report / 'methods_and_limitations.json').read_text())
            self.assertIsNone(methods['disease_conclusion'])

    def test_pathway_inference_requires_a_valid_explicit_background(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = root / 'pathways.tsv'
            data.write_text('gene\tpathway\nGENE_A\tSET_A\nGENE_B\tSET_A\n')
            manifest = {'resources': {'pathways': {'status': 'available', 'path': data.name, 'version': 'test'}}}
            receipt = root / 'manifest.json'
            self.assertEqual(m.local_pathways({'GENE_A'}, manifest, receipt, {})['status'], 'disabled_primary_transcript_withdrew_enrichment')
            self.assertEqual(m.local_pathways({'GENE_A'}, manifest, receipt, {'pathway_mode': 'descriptive'})['status'], 'descriptive_only')
            self.assertEqual(m.local_pathways({'GENE_A'}, manifest, receipt, {'pathway_mode': 'enrichment'})['status'], 'missing_background')
            background = root / 'background.txt'
            background.write_text('GENE_A\nGENE_B\nGENE_C\n')
            result = m.local_pathways({'GENE_A'}, manifest, receipt, {'pathway_mode': 'enrichment', 'pathway_background': str(background)})
            self.assertEqual(result['status'], 'research_only_explicit_background')
            self.assertGreaterEqual(result['results'][0]['p_value'], 0)
            self.assertLessEqual(result['results'][0]['BH_q_value'], 1)

    def test_star_split_preserved_but_not_counted_as_SNV_or_candidate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'raw.vcf'
            header = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\n'
            source.write_text(header + 'chr1\t10\t.\tC\tA\t50\tPASS\t.\tGT:FT\t0/1:PASS\n'
                              'chr1\t10\t.\tC\t*\t50\tPASS\t.\tGT:FT\t0/1:PASS\n')
            sequence = root / 'sequence.vcf.gz'
            self.assertEqual(m.partition_sequence_input(source, sequence), 1)
            merged = root / 'merged.vcf'
            m.merge_vep_preserving_all(source, sequence, merged)
            records = list(m.parse_vcf(merged))
            self.assertEqual(len(records), 2)
            self.assertEqual(records[1]['alt'], '*')
            manifest = dict(resources={})
            m.make_tables(merged, root, {}, manifest, root / 'manifest.json')
            counters = json.loads((root / 'consequence_counts.json').read_text())
            self.assertEqual(counters['cohort']['raw']['unique_normalized_alleles'], 1)
            self.assertEqual(counters['cohort']['raw']['variant_type_histogram'], {'SNV': 1})
            with gzip.open(root / 'transcript_samples.tsv.gz', 'rt') as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual(rows[1]['variant_type'], 'spanning_deletion_placeholder')
            self.assertEqual(rows[1]['candidate_eligible'], '0')
            self.assertEqual(rows[1]['sample_state'], 'not_sequence_allele_placeholder')

    def test_missing_FT_is_unknown_quality_and_duplicates_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            vcf = root / 'test.vcf'
            header = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\n'
            record = 'chr1\t10\t.\tC\tA\t50\tPASS\t.\tGT\t0/0\n'
            vcf.write_text(header + record)
            row = next(m.parse_vcf(vcf))
            self.assertEqual(m.record_state(row, row['samples']['SAMPLE_A']), 'unknown_quality_reference')
            vcf.write_text(header + record + record)
            with self.assertRaisesRegex(ValueError, 'Duplicate normalized'):
                list(m.parse_vcf(vcf))

    def test_gene_validity_is_not_assigned_to_an_unmatched_disease(self):
        with tempfile.TemporaryDirectory() as temporary:
            connection = sqlite3.connect(':memory:')
            connection.execute('CREATE TABLE evidence (variant TEXT, transcript TEXT, evidence TEXT)')
            genes = {'GENE_A': [{'source': 'clingen_validity', 'evidence': {
                'CLASSIFICATION': 'Definitive', 'DISEASE ID (MONDO)': 'MONDO:TEST1', 'DISEASE LABEL': 'Other disease', 'MOI': 'AR'}}]}
            record = dict(chrom='chr1', pos='10', ref='C', alt='A', filter='PASS')
            csq = dict(SYMBOL='GENE_A', Consequence='missense_variant')
            manifest = {'resources': {}}
            context = ({'disease_context': {'mondo_ids': ['MONDO:TEST2']}}, manifest, genes, {}, connection)
            row = m.annotation_row(record, csq, 'SAMPLE_A', {'GT': '0/1', 'FT': 'PASS'}, context)
            self.assertEqual(row['gene_validity_any_context'], 1)
            self.assertEqual(row['gene_validity_rank'], 0)
            context[0]['disease_context']['mondo_ids'] = ['MONDO:TEST1']
            row = m.annotation_row(record, csq, 'SAMPLE_A', {'GT': '0/1', 'FT': 'PASS'}, context)
            self.assertEqual(row['gene_validity_rank'], 1)
            context[0]['disease_context'] = 'hypertrophic cardiomyopathy research cohort'
            row = m.annotation_row(record, csq, 'SAMPLE_A', {'GT': '0/1', 'FT': 'PASS'}, context)
            self.assertEqual(row['gene_validity_rank'], 0)
            self.assertIn('descriptive_context_without_controlled', row['gene_validity_context_status'])
            connection.close()

    def test_local_adapter_hash_and_schema_are_required(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'local.tsv'
            source.write_text('chrom\tpos\tref\talt\tFATHMM\n1\t10\tC\tA\tD\n')
            manifest_path = root / 'manifest.json'
            manifest_path.write_text(json.dumps({'schema_version': 1, 'status': 'PASS', 'assembly': 'GRCh38', 'resources': {}, 'artifacts': []}))
            item = dict(name='FATHMM', local_adapter=True, status='available', path=str(source), bytes=source.stat().st_size,
                sha256=m.sha256(source), version='synthetic', source_url='https://example.org/synthetic-fixture',
                license_status='synthetic_test_data', format='variant_evidence_tsv')
            result = m.load_resources(manifest_path, local_adapters=[item])
            self.assertEqual(result['resources']['FATHMM']['status'], 'available')
            with self.assertRaisesRegex(ValueError, 'replace locked core'):
                m.load_resources(manifest_path, local_adapters=[dict(item, name='clinvar')])
            source.write_text('CORRUPTED')
            with self.assertRaisesRegex(ValueError, 'differs'):
                m.load_resources(manifest_path, local_adapters=[item])

    def test_malicious_cache_tar_stops_before_extraction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / 'bad.tar.gz'
            with tarfile.open(archive, 'w:gz') as tar:
                member = tarfile.TarInfo('homo_sapiens/../../escape')
                member.size = 1
                tar.addfile(member, io.BytesIO(b'x'))
            with self.assertRaisesRegex(ValueError, 'Unsafe cache'):
                prep.extract_cache(archive, root / 'cache', limit=1024, reserve=0)
            self.assertFalse((root / 'cache').exists())

    def test_resource_base_and_v2_reuse_validated_hardlinks_without_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / 'synthetic-cache.tar.gz'
            cache_info = b'gnomAD 4.1\nsift 6.2.1\npolyphen 2.2.3\nMANE 1.4\n'
            with tarfile.open(archive, 'w:gz') as tar:
                member = tarfile.TarInfo('homo_sapiens/115_GRCh38/info.txt')
                member.size = len(cache_info)
                tar.addfile(member, io.BytesIO(cache_info))
            cfg = dict(schema_version=1, bundle_id='synthetic-base', assembly='GRCh38',
                vep=dict(image='ensemblorg/ensembl-vep@sha256:' + '0' * 64, version='115.2', cache_version=115),
                resources=[dict(name='vep_cache', version='synthetic115', enabled=True, destination='synthetic-cache.tar.gz',
                    url='https://ftp.ensembl.org/synthetic-not-requested', bytes=archive.stat().st_size, sha256=prep.sha256(archive))])
            base = root / 'base.json'
            base.write_text(json.dumps(cfg))
            final = root / 'final.json'
            final.write_text(json.dumps(dict(cfg, bundle_id='synthetic-v2', base_config=base.name, resources=[])))
            free = type('Disk', (), {'free': 1000 * 2 ** 30})()
            with patch.object(prep, 'command', return_value='VEP 115.2'), patch.object(prep.shutil, 'disk_usage', return_value=free):
                result = prep.prepare(final, root / 'resources', True, False, archive, prep.sha256(archive))
                again = prep.prepare(final, root / 'resources', True, False, archive, prep.sha256(archive))
            self.assertEqual(result, again)
            original = root / 'resources/synthetic-base/cache/homo_sapiens/115_GRCh38/info.txt'
            reused = root / 'resources/synthetic-v2/cache/homo_sapiens/115_GRCh38/info.txt'
            self.assertEqual(original.stat().st_ino, reused.stat().st_ino)
            self.assertEqual(json.loads((result / 'resource-manifest.json').read_text())['status'], 'PASS')
            self.assertEqual(original.read_bytes(), cache_info)


if __name__ == '__main__':
    unittest.main()
