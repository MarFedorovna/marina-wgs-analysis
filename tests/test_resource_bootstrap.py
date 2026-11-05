"""No network: duplicate prevention, read-only plans, budget anchors and range preservation."""
import contextlib
import datetime as dt
import fcntl
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/prepare_all_resources.py'
spec = importlib.util.spec_from_file_location('resource_bootstrap', SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class BootstrapTests(unittest.TestCase):
    def seed(self, root):
        (root / 'private').mkdir()
        created = (dt.datetime.now(dt.timezone.utc) - dt.timedelta(hours=1)).isoformat()
        state = dict(identity=dict(created_utc=created, hourly_usd=2.39, planned_vm_usd=40),
            status='INSTALLER_RUNNING', current_pid=123456)
        (root / 'private/databases-supervisor-status.json').write_text(json.dumps(state))
        return created

    def test_dry_run_is_read_only(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.seed(root)
            before = {p.relative_to(root).as_posix(): p.read_bytes() for p in root.rglob('*') if p.is_file()}
            with patch.object(m, 'check_package', return_value=None), patch.object(m, 'inventory', return_value=[]), \
                    patch.object(m, 'lock_busy', return_value=False), patch.object(m.subprocess, 'Popen', side_effect=AssertionError('No jobs in plan')), \
                    contextlib.redirect_stdout(io.StringIO()):
                m.prepare(root, execute=False)
            after = {p.relative_to(root).as_posix(): p.read_bytes() for p in root.rglob('*') if p.is_file()}
            self.assertEqual(before, after)

    def test_repeated_execute_reuses_active_supervisor_and_partial_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.seed(root)
            partial = root / 'owned.agent-partial'; partial.write_bytes(b'PUBLIC VERIFIED RANGE')
            state = root / 'owned.agent-parts.json'; state.write_text('{"completed":{"0":"synthetic"}}')
            originals = partial.read_bytes(), state.read_bytes()
            rows = [dict(resource='alphamissense', mode='auto', status='ready_receipt')]
            with patch.object(m, 'check_package', return_value=None), patch.object(m, 'inventory', return_value=rows), \
                    patch.object(m, 'lock_busy', return_value=True), patch.object(m.subprocess, 'Popen', side_effect=AssertionError('Duplicate job')), \
                    patch.object(m, 'run_bounded', side_effect=AssertionError('No preparation needed')), contextlib.redirect_stdout(io.StringIO()):
                for _ in range(2):
                    self.assertEqual(m.prepare(root, execute=True)['status'], 'ALREADY_RUNNING_REUSING_READY_SMALLS_AND_ALPHA')
            self.assertEqual(originals, (partial.read_bytes(), state.read_bytes()))

    def test_supervisor_lock_is_detected_without_writing(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'lock'; path.write_text('unchanged')
            with path.open('r+') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                self.assertTrue(m.lock_busy(path))
            self.assertFalse(m.lock_busy(path))
            self.assertEqual(path.read_text(), 'unchanged')

    def test_creation_anchor_and_existing_rate_cannot_be_reset(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); original = self.seed(root)
            self.assertEqual(m.budget(root)[0], original)
            with self.assertRaises(ValueError):
                m.budget(root, dt.datetime.now(dt.timezone.utc).isoformat())
            path = root / 'private/databases-supervisor-status.json'
            state = json.loads(path.read_text()); state['identity']['hourly_usd'] = 3
            path.write_text(json.dumps(state))
            with self.assertRaises(ValueError):
                m.budget(root)


if __name__ == '__main__':
    unittest.main()
