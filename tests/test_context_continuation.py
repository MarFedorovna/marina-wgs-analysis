import datetime as dt
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('continuation', Path(__file__).resolve().parents[1] / 'scripts/continue_tissue_protein_context.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class ContinuationTests(unittest.TestCase):
    def fixture(self, root):
        context = root / 'context'; context.mkdir()
        (context / 'context-index.json').write_text('{}'); (context / 'manifest.json').write_text('{}')
        scripts = Path(m.__file__).resolve().parent
        rule = dict(schema_version=1, review_status='ROOT_REVIEWED_EXACT_LOCAL_ONLY_SIDECAR_RULE', reviewed_utc='2026-11-05T19:00:00Z',
            code_sha256={str(p): m.sha(p) for p in [scripts / 'continue_tissue_protein_context.py', scripts / 'export_tissue_protein_context.py']},
            context_index=str(context / 'context-index.json'), context_index_sha256=m.sha(context / 'context-index.json'),
            context_manifest_sha256=m.sha(context / 'manifest.json'), workflow_root=str(root / 'workflow'), output=str(root / 'sidecar'),
            state_directory=str(root / 'state'), deadline_utc=m.ABSOLUTE_DEADLINE.isoformat(), interval_seconds=60, max_hours=8,
            no_network=True, no_primary_workflow_mutation=True)
        path = root / 'rule.json'; path.write_text(json.dumps(rule))
        return rule, path

    def test_review_pins_and_unknown_fields(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); rule, path = self.fixture(root)
            now = m.ABSOLUTE_DEADLINE - dt.timedelta(hours=1)
            m.validate(rule, path, m.sha(path), now)
            with self.assertRaisesRegex(ValueError, 'rule SHA'): m.validate(rule, path, '0' * 64, now)
            changed = dict(rule, unexpected=True)
            with self.assertRaisesRegex(ValueError, 'Unknown'): m.validate(changed, path, m.sha(path), now)
            rule['code_sha256'][str(Path(m.__file__).resolve())] = '0' * 64
            with self.assertRaisesRegex(ValueError, 'code digest'): m.validate(rule, path, m.sha(path), now)

    def test_deadline_and_state_gates(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); rule, path = self.fixture(root)
            with self.assertRaisesRegex(ValueError, 'deadline'): m.validate(rule, path, m.sha(path), m.ABSOLUTE_DEADLINE)
            rule['deadline_utc'] = (m.ABSOLUTE_DEADLINE + dt.timedelta(seconds=1)).isoformat()
            with self.assertRaisesRegex(ValueError, 'deadline'): m.validate(rule, path, m.sha(path), m.ABSOLUTE_DEADLINE - dt.timedelta(hours=1))
            self.assertFalse(m.readiness(root / 'workflow'))
            for stage in ['06_annotate', '07_report']:
                folder = root / 'workflow/stages' / stage; folder.mkdir(parents=True)
                (folder / 'workflow-receipt.json').write_text(json.dumps(dict(status='PASS')))
            self.assertTrue(m.readiness(root / 'workflow'))
            (root / 'workflow/stages/07_report/workflow-receipt.json').write_text(json.dumps(dict(status='FAIL')))
            with self.assertRaises(ValueError): m.readiness(root / 'workflow')


if __name__ == '__main__': unittest.main()
