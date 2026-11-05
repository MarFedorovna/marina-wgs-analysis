"""Synthetic budget/process identity checks; no network or provider operations."""
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
from types import SimpleNamespace

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/supervise_databases.py'
spec = importlib.util.spec_from_file_location('database_supervisor', SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class SupervisorTests(unittest.TestCase):
    def test_allowance_is_anchored_to_creation_and_cannot_be_increased(self):
        created = m.utc('2026-11-05T11:02:37Z')
        self.assertAlmostEqual((m.deadline(created, 2.39, 40) - created).total_seconds(), 40 / 2.39 * 3600, places=5)
        with self.assertRaises(ValueError):
            m.deadline(created, 2.39, 41)
        with self.assertRaises(ValueError):
            m.utc('2026-11-05T11:02:37')

    def test_only_transient_failures_are_retryable(self):
        self.assertTrue(m.transient_failure('ValueError: Direct download time limit reached; owned partial is resumable.'))
        self.assertTrue(m.transient_failure('TimeoutError: read timed out'))
        self.assertTrue(m.transient_failure('urllib.error.URLError: <urlopen error [Errno 111] Connection refused>'))
        self.assertFalse(m.transient_failure('ValueError: Resource SHA256 mismatch'))
        self.assertFalse(m.transient_failure('TimeoutError followed by cache metadata version unexpected'))
        self.assertFalse(m.transient_failure('ConnectionRefusedError followed by object identity changed'))
        self.assertFalse(m.transient_failure('HTTPError: HTTP Error 403: Forbidden'))

    def test_ready_requires_every_artifact_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            recipe = root / 'recipe.json'
            recipe.write_text('{"synthetic":true}')
            source = root / 'public.txt'
            source.write_text('synthetic resource')
            manifest = root / 'resource-manifest.json'
            manifest.write_text(json.dumps(dict(schema_version=1, status='PASS', assembly='GRCh38',
                config_sha256=m.digest(recipe), artifacts=[dict(path=source.name, bytes=source.stat().st_size, sha256=m.digest(source))])))
            self.assertTrue(m.verified_ready(manifest, recipe))
            source.write_text('corrupted')
            with self.assertRaises(ValueError):
                m.verified_ready(manifest, recipe)

    def test_network_backoff_is_capped_and_progress_resets_the_streak(self):
        self.assertEqual(m.retry_backoff_seconds(0, False), (1, 60))
        self.assertEqual(m.retry_backoff_seconds(1, False), (2, 120))
        self.assertEqual(m.retry_backoff_seconds(9, False), (10, 900))
        self.assertEqual(m.retry_backoff_seconds(9, True), (0, 60))

    @unittest.skipUnless(sys.platform.startswith('linux'), 'Linux owned installer integration')
    def test_deadline_policy_resumes_after_two_but_integrity_failure_stays_terminal(self):
        for prior_error, allowed in [('ConnectionRefusedError: Connection refused', True),
                                     ('Resource SHA256 mismatch', False)]:
            with self.subTest(prior_error=prior_error), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                for folder in ['private', 'logs', 'publicstages', 'scripts', 'config']:
                    (root / folder).mkdir()
                recipe = root / 'config/recipe.json'
                recipe.write_text('{"bundle_id":"synthetic"}')
                installer = root / 'publicstages/01_prepare_databases.sh'
                installer.write_text('exec python3 "' + str(root / 'scripts/prepare_annotation_resources.py') + '"\n')
                script = root / 'scripts/prepare_annotation_resources.py'
                script.write_text('''import json,hashlib,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
counter=root/'attempt-counter'
n=int(counter.read_text())+1 if counter.exists() else 1
counter.write_text(str(n))
if n==1:
 print('ConnectionRefusedError: Connection refused',flush=True);sys.exit(1)
target=root/'resources/annotation/synthetic';target.mkdir(parents=True)
artifact=target/'public.txt';artifact.write_text('PUBLIC SYNTHETIC RESOURCE')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
receipt=dict(status='PASS',schema_version=1,assembly='GRCh38',config_sha256=sha(root/'config/recipe.json'),artifacts=[dict(path='public.txt',bytes=artifact.stat().st_size,sha256=sha(artifact))])
(target/'resource-manifest.json').write_text(json.dumps(receipt))
''')
                (root / 'logs/annotation-databases.log').write_text(prior_error)
                state_path = root / 'private/databases-supervisor-status.json'
                state_path.write_text(json.dumps(dict(attempts_started=2, attempt_log_offset=0,
                    status='MAX_ATTEMPTS_REACHED_PARTIALS_PRESERVED')))
                args = SimpleNamespace(root=root, installer=installer, config=recipe,
                    created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), hourly_usd=2.39,
                    planned_usd=40, max_attempts=2, workers=1, retry_until_deadline=True, adopt_pid=None)
                with patch.object(m, 'retry_backoff_seconds', return_value=(1, 0)):
                    if allowed:
                        self.assertEqual(m.run(args), 0)
                        result = json.loads(state_path.read_text())
                        self.assertEqual(result['status'], 'READY_VALIDATED')
                        self.assertEqual(result['attempts_started'], 4)
                    else:
                        with self.assertRaises(ValueError):
                            m.run(args)
                        self.assertFalse((root / 'attempt-counter').exists())
                        # A second manual invocation must not bypass that terminal failure
                        # merely because the first invocation recorded SUPERVISOR_FAILURE.
                        with self.assertRaises(ValueError):
                            m.run(args)
                        self.assertFalse((root / 'attempt-counter').exists())

    @unittest.skipUnless(sys.platform.startswith('linux'), 'Linux process-group ownership checks')
    def test_stops_only_exact_owned_installer_process_group(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'publicstages').mkdir()
            (root / 'scripts').mkdir()
            installer = root / 'publicstages/01_prepare_databases.sh'
            installer.write_text('#!/bin/sh\n')
            script = root / 'scripts/prepare_annotation_resources.py'
            script.write_text('import time\ntime.sleep(60)\n')
            intended = subprocess.Popen([sys.executable, str(script)], start_new_session=True)
            unrelated = subprocess.Popen([sys.executable, '-c', 'import time;time.sleep(60)'], start_new_session=True)
            try:
                time.sleep(.1)
                token = m.process_token(intended.pid, installer)
                self.assertTrue(m.still_owned(token, installer))
                with self.assertRaises(ValueError):
                    m.process_token(unrelated.pid, installer)
                m.stop_owned(token, installer)
                intended.wait(timeout=3)
                self.assertIsNone(unrelated.poll())
            finally:
                if intended.poll() is None:
                    intended.terminate()
                if unrelated.poll() is None:
                    unrelated.terminate()
                intended.wait(timeout=3)
                unrelated.wait(timeout=3)


if __name__ == '__main__':
    unittest.main()
