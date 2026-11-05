"""Small offline fixtures: current receipts, sparse sizes, PID reuse and bounds."""
import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('monitor', Path(__file__).resolve().parents[1] / 'scripts/monitor_pipeline.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class MonitorTests(unittest.TestCase):
    def test_sparse_size_does_not_count_as_verified(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'sparse'
            with path.open('wb') as stream:
                stream.truncate(1000)
            item = m.range_progress({'identity': {'bytes': 1000, 'chunk_bytes': 100},
                'completed': {'0': 'a' * 64, '9': 'b' * 64, '99': 'c' * 64, '1': 'invalid'}}, path)
            self.assertEqual(item['logical_bytes_not_progress'], 1000)
            self.assertEqual(item['verified_bytes'], 200)
            self.assertEqual(item['percent'], 20)

    def test_current_v2_not_historical_failure_and_read_only(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); private = root / 'private'; private.mkdir()
            receipts = {'mitochondria-recovery-v2-queue-status.json': {'status': 'RUNNING', 'pid': 90},
                'mitochondria-recovery-queue-status.json': {'status': 'FAIL'},
                'mitochondria-recovery-status.json': {'status': 'RUNNING', 'completed': 3, 'policy': {'version': 2}},
                'databases-supervisor-status.json': {'status': 'INSTALLER_RUNNING', 'attempts_started': 4, 'current_pid': 91}}
            for name, data in receipts.items():
                (private / name).write_text(json.dumps(data))
            before = {p.name: p.read_bytes() for p in private.iterdir()}
            result = m.snapshot(root, proc=root / 'no-proc')
            self.assertEqual(result['mt_v2']['status'], 'RUNNING')
            self.assertEqual(result['mt_v2']['progress']['completed'], 3)
            self.assertEqual(result['databases']['attempts_started'], 4)
            self.assertFalse(result['databases']['installer']['verified_active'])
            self.assertEqual(before, {p.name: p.read_bytes() for p in private.iterdir()})

    def test_pid_reuse_and_zombies_are_not_current(self):
        with tempfile.TemporaryDirectory() as temp:
            proc = Path(temp); folder = proc / '123'; folder.mkdir()
            (folder / 'stat').write_text('123 (python3) S 1 2')
            (folder / 'cmdline').write_bytes(b'python3\0unrelated.py\0')
            self.assertFalse(m.process(123, 'installer.py', proc)['verified_active'])
            (folder / 'cmdline').write_bytes(b'python3\0installer.py\0')
            self.assertTrue(m.process(123, 'installer.py', proc)['verified_active'])
            (folder / 'stat').write_text('123 (python3) Z 1 2')
            self.assertFalse(m.process(123, 'installer.py', proc)['verified_active'])

    def test_budget_stops_observer_without_sleep_or_actuation(self):
        args = argparse.Namespace(max_hours=8, interval=60, root=Path('/missing'), local=False, watch=True, report_dir=None)
        calls = []
        with contextlib.redirect_stdout(io.StringIO()):
            m.run(args, collect=lambda *_: {'databases': {'budget_deadline_reached': True}},
                  monotonic=lambda: 0, sleep=lambda _: calls.append('sleep'))
        self.assertEqual(calls, [])
        args.max_hours = 9
        with self.assertRaises(ValueError):
            m.run(args)

    def test_local_mode_never_reads_private_config(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); (root / '.local').mkdir()
            (root / '.local/native-sftp-direct-config.json').write_text('SECRET')
            original = m.load
            def allowed(path):
                self.assertEqual(path.name, 'native-sftp-direct-status.json')
                return original(path)
            with patch.object(m, 'load', side_effect=allowed), patch.object(m, 'discover', side_effect=AssertionError('no proc')):
                self.assertNotIn('SECRET', json.dumps(m.snapshot(root, local=True)))

    def test_historical_log_errors_excluded(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'db.log'
            old = b'Traceback: checksum mismatch\n'
            path.write_bytes(old + b'connection refused; retry continues\n')
            result = m.current_attempt_errors(path, len(old))
            self.assertEqual(result['signals'], ['connection_refused'])

    def test_watch_has_actual_elapsed_time_bound(self):
        args = argparse.Namespace(max_hours=0.001, interval=60, root=Path('/missing'), local=False, watch=True, report_dir=None)
        clock = [0.0]
        def tick(seconds):
            clock[0] += seconds
        with contextlib.redirect_stdout(io.StringIO()):
            m.run(args, collect=lambda *_: {}, monotonic=lambda: clock[0], sleep=tick)
        self.assertLessEqual(clock[0], 3.6)


if __name__ == '__main__':
    unittest.main()
