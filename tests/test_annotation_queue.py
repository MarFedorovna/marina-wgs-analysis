"""Continuation must reject altered accepted artifacts and sealed source files."""
import datetime as dt
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/annotation_queue.py'
spec = importlib.util.spec_from_file_location('annotation_queue', SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class QueueTests(unittest.TestCase):
    def test_budget_is_anchored_to_actual_creation_and_cannot_increase_allowance(self):
        end = m.budget_deadline('2026-11-05T11:02:37Z', 2.39, 40)
        self.assertEqual(end.isoformat(), '2026-11-06T03:46:48.046025+00:00')
        with self.assertRaises(ValueError):
            m.budget_deadline('2026-11-05T11:02:37Z', 2.39, 41)

    def test_accepted_stage_requires_actual_artifact_hash_not_only_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact = root / 'synthetic.vcf'
            artifact.write_text('PUBLIC_SYNTHETIC\n')
            receipt = {'status': 'PASS', 'artifacts': {'synthetic.vcf':
                {'bytes': artifact.stat().st_size, 'sha256': m.digest(artifact)}}}
            (root / 'workflow-receipt.json').write_text(json.dumps(receipt))
            self.assertTrue(m.verify_stage(root))
            artifact.write_text('PUBLIC_ALTERED__\n')
            with self.assertRaises(ValueError):
                m.verify_stage(root)

    def test_sealed_file_is_idempotent_but_changed_source_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / 'code.py'
            m.exclusive(target, b'print("synthetic")\n')
            m.exclusive(target, b'print("synthetic")\n')
            with self.assertRaises(ValueError):
                m.exclusive(target, b'print("changed")\n')
            self.assertEqual(target.read_bytes(), b'print("synthetic")\n')


if __name__ == '__main__':
    unittest.main()
