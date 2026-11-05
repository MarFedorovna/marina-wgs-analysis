import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/export_consequence_tables.py'
SPEC = importlib.util.spec_from_file_location('consequence_export', SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ConsequenceExportTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'synthetic.tsv.gz'
        self.rows = []
        def add(variant, sample, gt, consequence, **extra):
            row = dict(sample=sample, variant_id=variant, chrom='1', pos=variant[1:], ref='A', alt='T',
                       GT=gt, FT='PASS', site_filter='PASS', consequence=consequence,
                       transcript='ENST_SYNTHETIC', clinvar_significance='UNASSESSED',
                       candidate_eligible='false', mechanism_status='not_inferred',
                       AlphaMissense_score='', local_evidence_json='{"synthetic":true}')
            row.update(extra)
            self.rows.append(row)
        add('v1', 'A', '0/1', 'stop_gained', FT='LowDP', transcript='TX1')
        add('v1', 'A', '0/1', 'stop_gained', FT='LowDP', transcript='TX2')
        add('v1', 'A', '0/1', 'missense_variant', FT='LowDP', transcript='TX3')
        add('v1', 'B', './.', 'stop_gained')
        add('v2', 'A', '1/.', 'stop_lost&novel_unknown_term', site_filter='LowQUAL')
        add('v2', 'B', '0/0', 'stop_lost')
        add('v3', 'A', './.', 'intron_variant')
        add('v3', 'B', './.', 'intron_variant')
        add('v4', '=sample()', '0/2', '5_prime_UTR_variant', FT='.')
        with gzip.open(self.source, 'wt', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(self.rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(self.rows)
        self.original = self.source.read_bytes()

    def run_export(self):
        self.out = self.root / 'export'
        self.manifest = MODULE.export(self.source, self.out, synthetic=True)
        return self.manifest

    def table(self, category, sample='cohort'):
        item = next(x for x in self.manifest['counts'] if x['category'] == category and x['scope_id'] == sample)
        with gzip.open(self.out / item['path'], 'rt', encoding='utf-8', newline='') as stream:
            return list(csv.DictReader(stream, delimiter='\t'))

    def test_overlap_dedup_quality_and_all_source_evidence(self):
        self.run_export()
        self.assertEqual(len(self.table('stop_gained', 'A')), 1)
        self.assertEqual(len(self.table('missense_variant', 'A')), 1)
        self.assertEqual(len(self.table('pLoF_candidates', 'A')), 1)
        self.assertEqual(self.table('stop_gained', 'A')[0]['source_transcript_rows'], '3')
        self.assertEqual(json.loads(self.table('stop_gained', 'A')[0]['quality_levels_json']), ['FILTERED'])
        self.assertEqual(json.loads(self.table('stop_gained', 'A')[0]['source_FT_json']), ['LowDP'])
        self.assertEqual(self.source.read_bytes(), self.original)
        with gzip.open(self.out / 'evidence/transcript_effects.tsv.gz', 'rt', encoding='utf-8', newline='') as stream:
            raw = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual(len(raw), len(self.rows))
        for actual, expected in zip(raw, self.rows):
            self.assertEqual({key: actual[key] for key in expected}, expected)
        for name, artifact in self.manifest['artifacts'].items():
            self.assertEqual(MODULE.sha256(self.out / name), artifact['sha256'])

    def test_partial_unknown_dynamic_zero_and_no_automatic_gof(self):
        self.run_export()
        self.assertEqual(self.table('stop_gained', 'B'), [])
        self.assertEqual(self.table('intron_variant', 'A'), [])
        self.assertEqual(len(self.table('intron_variant')), 1)
        unknown = self.table('intron_variant')[0]
        self.assertEqual(json.loads(unknown['observed_alt_samples_json']), [])
        self.assertEqual(json.loads(unknown['unknown_gt_samples_json']), ['A', 'B'])
        partial = self.table('stop_lost', 'A')[0]
        self.assertEqual(json.loads(partial['partial_alt_samples_json']), ['A'])
        self.assertEqual(json.loads(partial['unknown_gt_samples_json']), ['A'])
        self.assertEqual(json.loads(partial['mechanism_review_statuses_json']), ['UNASSESSED'])
        self.assertEqual(json.loads(partial['source_mechanism_status_json']), ['not_inferred'])
        self.assertEqual(self.table('pLoF_candidates', 'A')[0]['variant_id'], 'v1')
        self.assertEqual(len(self.table('novel_unknown_term', 'A')), 1)
        self.assertEqual(len(self.table('5_prime_UTR_variant', '=sample()')), 1)
        self.assertEqual(self.table('start_lost', 'A'), [])
        self.assertIn("'=sample()", (self.out / 'index.tsv').read_text(encoding='utf-8'))
        self.assertIn('SYNTHETIC FIXTURE', (self.out / 'index.html').read_text(encoding='utf-8'))

    def test_existing_output_and_nonpass_receipt_are_rejected(self):
        self.run_export()
        with self.assertRaisesRegex(ValueError, 'refusing overwrite'):
            MODULE.export(self.source, self.out)
        receipt = self.root / 'receipt.json'
        receipt.write_text(json.dumps({'status': 'RUNNING'}))
        with self.assertRaisesRegex(ValueError, 'PASS'):
            MODULE.export(self.source, self.root / 'new', receipt=receipt)
        self.assertFalse((self.root / 'new').exists())

    def test_genotype_token_semantics(self):
        self.assertEqual(MODULE.gt_flags('./.'), (False, True, False, False, False))
        self.assertEqual(MODULE.gt_flags('1|.'), (True, True, False, True, False))
        self.assertEqual(MODULE.gt_flags('0/2'), (True, False, False, False, False))
        self.assertEqual(MODULE.gt_flags('0'), (False, False, True, False, False))
        self.assertEqual(MODULE.gt_flags('bad'), (False, True, False, False, True))

    @unittest.skipUnless(os.name == 'posix' and shutil.which('bash'), 'Unix wrapper tested on server')
    def test_wrapper_gate_and_receipt_checksum(self):
        wrapper = SCRIPT.with_name('run_consequence_export.sh')
        workflow = self.root / 'workflow'
        stage6 = workflow / 'stages/06_annotate'
        stage7 = workflow / 'stages/07_report'
        stage6.mkdir(parents=True)
        stage7.mkdir(parents=True)
        destination = self.root / 'sidecar'
        command = ['bash', str(wrapper), '--workflow-root', str(workflow), '--output', str(destination), '--synthetic']
        not_ready = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(not_ready.returncode, 0)
        self.assertIn('NOT_READY', not_ready.stderr)
        self.assertFalse(destination.exists())
        master = stage6 / 'transcript_samples.tsv.gz'
        master.write_bytes(self.original)
        stage6_receipt = {'status': 'PASS', 'identity': {'synthetic': True},
                          'artifacts': {master.name: {'sha256': 'wrong'}}}
        (stage6 / 'workflow-receipt.json').write_text(json.dumps(stage6_receipt))
        (stage7 / 'workflow-receipt.json').write_text(json.dumps({'status': 'PASS', 'synthetic': True}))
        mismatch = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(mismatch.returncode, 0)
        self.assertIn('SHA256', mismatch.stderr)
        self.assertFalse(destination.exists())
        stage6_receipt['artifacts'][master.name]['sha256'] = MODULE.sha256(master)
        (stage6 / 'workflow-receipt.json').write_text(json.dumps(stage6_receipt))
        complete = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(complete.returncode, 0, complete.stderr)
        manifest = json.loads((destination / 'manifest.json').read_text())
        self.assertEqual(manifest['status'], 'PASS')
        self.assertTrue(manifest['synthetic'])
        self.assertIn('continuation-receipts.json', manifest['artifacts'])
        self.assertEqual(master.read_bytes(), self.original)


if __name__ == '__main__':
    unittest.main()
