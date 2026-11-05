import csv
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('context', Path(__file__).resolve().parents[1] / 'scripts/export_tissue_protein_context.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def index():
    return dict(status='PUBLIC_GENE_CONTEXT_READY_SEPARATE_SIDECAR_ONLY', patient_input_used=False, verified_raw_sources={},
        genes={'HGNC:7551': dict(gene='MYBPC3', gtex_gene=[dict(median=5, unit='TPM')],
            gtex_transcripts={'ENST000001.2': [dict(median=3, unit='TPM')]},
            pext_transcripts={'ENST000001.1': [dict(tissue='heart_left_ventricle', value=2)]},
            pext_status='GENE_CONTEXT_AVAILABLE_ALLELE_COORDINATE_CONTRACT_UNVERIFIED', uniprot=[dict(accession='Q14896')])})


class ContextTests(unittest.TestCase):
    def test_exact_transcript_version_never_borrow_other_version(self):
        row = dict(hgnc_id='7551', gene='MYBPC3', transcript='ENST000001.3')
        result = m.supplement(row, index(), 1)
        self.assertEqual(result['gtex_v11_transcript_match_status'], 'TRANSCRIPT_VERSION_MISMATCH')
        self.assertEqual(json.loads(result['gtex_v11_transcript_heart_medians_json']), [])
        row['transcript'] = 'ENST000001.2'
        self.assertEqual(json.loads(m.supplement(row, index(), 1)['gtex_v11_transcript_heart_medians_json'])[0]['median'], 3)
        row['transcript'] = 'ENST000001'
        self.assertEqual(m.supplement(row, index(), 1)['gtex_v11_transcript_match_status'], 'VERSIONED_TRANSCRIPT_REQUIRED')

    def test_symbol_is_not_identity_fallback_and_conflict_rejected(self):
        self.assertIsNone(m.supplement(dict(gene='MYBPC3', hgnc_id=''), index(), 1))
        with self.assertRaises(ValueError):
            m.supplement(dict(gene='MYH7', hgnc_id='7551'), index(), 1)

    def test_no_allele_or_individual_expression_claim(self):
        result = m.supplement(dict(hgnc_id='7551', transcript='ENST000001.2'), index(), 1)
        self.assertTrue(result['pext_allele_assignment_status'].startswith('NOT_ASSIGNED'))
        self.assertTrue(result['protein_domain_allele_assignment_status'].startswith('NOT_ASSIGNED'))
        self.assertIn('NOT_INDIVIDUAL', result['reference_context_scope'])

    def test_stage_gates_hash_and_per_sample_preservation(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); workflow = root / 'workflow'
            stages = workflow / 'stages'
            (stages / '06_annotate').mkdir(parents=True); (stages / '07_report').mkdir()
            source = stages / '06_annotate/transcript_samples.tsv.gz'
            with gzip.open(source, 'wt', newline='') as stream:
                writer = csv.DictWriter(stream, ['sample', 'variant_id', 'gene', 'hgnc_id', 'transcript', 'GT', 'FT', 'GQ'], delimiter='\t'); writer.writeheader()
                for sample, gt, ft in [('A', '0/1', 'PASS'), ('B', './.', 'LOW_GQ')]:
                    writer.writerow(dict(sample=sample, variant_id='11:1:T:G', gene='MYBPC3', hgnc_id='7551', transcript='ENST000001.2', GT=gt, FT=ft, GQ='12'))
            receipts = [stages / p / 'workflow-receipt.json' for p in ('06_annotate', '07_report')]
            receipts[0].write_text(json.dumps(dict(status='PASS', artifacts={source.name: dict(bytes=source.stat().st_size, sha256=m.sha(source))})))
            receipts[1].write_text(json.dumps(dict(status='WAITING')))
            context = root / 'context'; context.mkdir(); path = context / 'context-index.json'; path.write_text(json.dumps(index()))
            (context / 'manifest.json').write_text(json.dumps(dict(files={path.name: dict(bytes=path.stat().st_size, sha256=m.sha(path))})))
            with self.assertRaisesRegex(ValueError, 'NOT_READY'):
                m.enrich(workflow, path, root / 'blocked', '2099-01-01T00:00:00+00:00')
            self.assertFalse((root / 'blocked').exists())
            receipts[1].write_text(json.dumps(dict(status='PASS')))
            before = m.sha(source)
            result = m.enrich(workflow, path, root / 'sidecar', '2099-01-01T00:00:00+00:00')
            self.assertEqual(result['counters']['samples'], {'A': 1, 'B': 1})
            self.assertEqual(m.sha(source), before)
            with gzip.open(root / 'sidecar/public39-reference-context.tsv.gz', 'rt') as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual([(r['sample'], r['GT'], r['FT'], r['GQ']) for r in rows], [('A', '0/1', 'PASS', '12'), ('B', './.', 'LOW_GQ', '12')])
            with self.assertRaises(ValueError):
                m.enrich(workflow, path, workflow / 'forbidden', '2099-01-01T00:00:00+00:00')
            source.write_bytes(source.read_bytes() + b'X')
            with self.assertRaisesRegex(ValueError, 'differs'):
                m.enrich(workflow, path, root / 'tampered', '2099-01-01T00:00:00+00:00')

    def test_raw_source_hash_tamper_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'source.json'; path.write_text('{}')
            receipt = dict(bytes=2, sha256=m.sha(path)); path.write_text('[]')
            with self.assertRaises(ValueError): m.checked_source(path, receipt)


if __name__ == '__main__': unittest.main()
