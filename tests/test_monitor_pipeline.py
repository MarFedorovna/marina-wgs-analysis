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
    def direct_fixture(self, root, reviewed=True, fork_child=False, version='v2'):
        private = root / 'private'; scripts = private / 'candidate/scripts'; scripts.mkdir(parents=True)
        script = scripts / ('direct_fully_annotated_vcf_' + version + '.py'); script.write_text('# synthetic observer fixture')
        state = private / ('direct-core33-state-' + version); state.mkdir()
        (state / 'status.json').write_text(json.dumps(dict(status='DIRECT_FROZEN_ROWS_INLINE_CONTEXT_RUNNING',
            pid=501, workers=24, expected_rows=100000, batches=2)))
        rule_path = private / 'candidate/rule.json'
        rule_path.write_text(json.dumps(dict(root_reviewed=reviewed, status='ROOT_REVIEWED_DIRECT_CORE33_V2_RULE',
            driver={'path': str(script), 'sha256': 'a' * 64}, state_directory=str(state),
            output=str(root / 'outputs/direct-core-context-v1'), workers=24, worker_cpu_ids=list(range(2, 26)))))
        proc = root / 'proc'; proc.mkdir()
        command = ['python3', str(script), 'run', '--rule', str(rule_path), '--reviewed-rule-sha256', 'b' * 64, '--execute']
        for pid in ([501, 502] if fork_child else [501]):
            folder = proc / str(pid); folder.mkdir(); (folder / 'stat').write_text(str(pid) + ' (python3) S 1 2')
            (folder / 'cmdline').write_bytes(b'\0'.join(a.encode() for a in command) + b'\0')
        batch = state / 'work-one/result-00000'; batch.mkdir(parents=True)
        (batch / 'batch-proof.json').write_text(json.dumps(dict(id=0, first_i=1, last_i=50000, row_count=50000, record_count=1000)))
        old = private / 'direct-core33-state-v1'; old.mkdir()
        (old / 'ROOT-superseded.json').write_text(json.dumps(dict(status='ROOT_SUPERSEDED', new_pid=501)))
        return proc, state

    def resume_fixture(self, root):
        proc, state = self.direct_fixture(root, fork_child=True, version='v3')
        status = json.loads((state / 'status.json').read_text())
        status.update(reused_rows=50000, reused_batches=1)
        (state / 'status.json').write_text(json.dumps(status))
        audit = dict(reused_rows=50000, reused_batches=1, source_origin_mapping=[dict(reuse=True,
            old_job=dict(first_i=50001,last_i=100000,expected_rows=50000,records=1000),new_batch_ids=[1])])
        (state / 'work-one/resume-audit.json').write_text(json.dumps(audit))
        return proc, state

    def test_resume_v3_includes_cached_rows_without_reopening_cached_bodies(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, _ = self.resume_fixture(root)
            original = Path.open
            def metadata_only(path, *args, **kwargs):
                self.assertNotIn(path.suffix, ('.gz','.sqlite','.bam'))
                return original(path,*args,**kwargs)
            with patch.object(Path,'open',metadata_only):
                result = m.direct_core_controller(root/'private',proc)
            self.assertEqual(result['active_pids'],[501])
            self.assertEqual(result['progress']['completed_rows'],100000)
            self.assertEqual(result['progress']['reused_rows_observed'],50000)
            self.assertEqual(result['progress']['row_progress_percent'],100)

    def test_resume_v3_missing_or_conflicting_cache_metadata_never_claims_percentage(self):
        for missing in (True,False):
            with self.subTest(missing=missing),tempfile.TemporaryDirectory() as temp:
                root=Path(temp);proc,state=self.resume_fixture(root)
                path=state/'work-one/resume-audit.json'
                if missing:path.unlink()
                else:
                    audit=json.loads(path.read_text());audit['reused_rows']=49999
                    path.write_text(json.dumps(audit))
                progress=m.direct_core_controller(root/'private',proc)['progress']
                self.assertNotIn('row_progress_percent',progress)
                self.assertIn(progress['status'],('RESUME_CACHE_METADATA_NOT_OBSERVED','INCONSISTENT_RESUME_CACHE_METADATA'))

    def test_resume_v3_assembly_rows_are_generation_progress_not_final_delivery(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);proc,state=self.resume_fixture(root)
            (state/'status.json').write_text(json.dumps(dict(pid=501,status='ASSEMBLING_AND_INDEXING_DIRECT_CORE33_VCF',rows=100000)))
            result=m.snapshot(root,proc=proc)
            self.assertEqual(result['deliverables']['direct_core_context']['controller']['progress']['row_progress_percent'],100)
            self.assertNotEqual(result['deliverables']['fully_annotated_vcf']['status'],'COMPLETE_RECEIPT_OBSERVED')

    def test_completed_spliceai_receipt_does_not_display_stale_partial_progress_as_final_count(self):
        value=dict(spliceai=dict(status='COMPLETE_RECEIPT_OBSERVED',progress=dict(variants_completed=9,variants_estimated=10)))
        output=m.human(value)
        self.assertIn('SpliceAI: COMPLETE_RECEIPT_OBSERVED',output)
        self.assertNotIn('вариантов 9/10',output)

    def test_direct24_current_authority_ignores_fork_worker_argv_and_reports_real_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, _ = self.direct_fixture(root, fork_child=True)
            result = m.snapshot(root, proc=proc)
            self.assertEqual(result['annotation']['status'], 'DIRECT_CORE33_RUNNING')
            direct = result['deliverables']['direct_core_context']
            self.assertEqual(direct['controller']['active_pids'], [501])
            self.assertEqual(direct['controller']['workers'], 24)
            progress = direct['controller']['progress']
            self.assertEqual(progress['completed_rows'], 50000)
            self.assertEqual(progress['row_progress_percent'], 50)
            self.assertFalse(progress['payloads_read'])
            self.assertEqual(direct['historical_v1_superseded']['status'], 'ROOT_SUPERSEDED')

    def test_unreviewed_direct_rule_and_orphan_pool_worker_are_not_current_authority(self):
        for reviewed in (False, True):
            with self.subTest(reviewed=reviewed), tempfile.TemporaryDirectory() as temp:
                root = Path(temp); proc, _ = self.direct_fixture(root, reviewed=reviewed, fork_child=True)
                if reviewed: (proc / '501/stat').write_text('501 (python3) Z 1 2')
                result = m.direct_core_controller(root / 'private', proc)
                self.assertEqual(result['active_pids'], [])

    def test_overlapping_batch_metadata_never_becomes_a_progress_percentage(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, state = self.direct_fixture(root)
            second = state / 'work-one/result-00001'; second.mkdir()
            (second / 'batch-proof.json').write_text(json.dumps(dict(id=1, first_i=40001, last_i=90000, row_count=50000, record_count=1000)))
            progress = m.direct_core_controller(root / 'private', proc)['progress']
            self.assertEqual(progress['status'], 'OVERLAPPING_BATCH_METADATA')
            self.assertNotIn('row_progress_percent', progress)

    def test_intermediate_core33_pass_does_not_claim_final_spliceai_delivery(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, _ = self.direct_fixture(root)
            output = root / 'outputs/direct-core-context-v1'; output.mkdir(parents=True)
            (output / 'manifest.json').write_text(json.dumps(dict(status='INTERMEDIATE_CORE33_PASS', final_delivery_complete=False,
                original_records_fields_preserved=True, packet_decode_verified=True, direct_core_prefix_verified=True)))
            (output / 'workflow-receipt.json').write_text(json.dumps(dict(status='PASS', stage='direct-core-context-v1', identity='c' * 64,
                artifacts={name: {'bytes': 1, 'sha256': 'd' * 64} for name in ('cohort.core-context.vcf.gz', 'cohort.core-context.vcf.gz.tbi', 'core-manifest.json')})))
            result = m.snapshot(root, proc=proc)['deliverables']
            self.assertEqual(result['direct_core_context']['status'], 'COMPLETE_RECEIPT_OBSERVED')
            self.assertNotEqual(result['fully_annotated_vcf']['status'], 'COMPLETE_RECEIPT_OBSERVED')
            (proc/'501/stat').write_text('501 (python3) Z 1 2')
            complete=m.snapshot(root,proc=proc)
            self.assertEqual(complete['annotation']['status'],'DIRECT_CORE33_COMPLETE_WAITING_FINAL')

    def test_direct_batch_observer_never_opens_sqlite_or_packet_files(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, state = self.direct_fixture(root)
            (state / 'work-one/result-00000/body.vcf.gz').write_bytes(b'PRIVATE_PAYLOAD')
            (state / 'work-one/snapshot.sqlite').write_bytes(b'PRIVATE_DATABASE')
            original = Path.open
            def allowed(path, *args, **kwargs):
                self.assertNotIn(path.suffix, ('.gz', '.sqlite', '.bam'))
                return original(path, *args, **kwargs)
            with patch.object(Path, 'open', allowed):
                result = m.snapshot(root, proc=proc)
            self.assertTrue(result['no_hashing'])
            self.assertNotIn('PRIVATE_PAYLOAD', json.dumps(result))

    def views_fixture(self, root, *, active=True):
        state = root / 'private/packet-fastviews-continuation-v5'
        state.mkdir(parents=True)
        (state / 'status.json').write_text(json.dumps(dict(
            status='WAITING_DIRECT_FULLY_ANNOTATED_VCF_PASS_NO_PRIMARY_ACTUATION',
            pid=222, rule_sha256='a' * 64)))
        old = root / 'private/cohort-sample-views-continuation-v2'; old.mkdir()
        (old / 'status.json').write_text(json.dumps(dict(status='STOPPED_OLD_VIEWS_RULE', pid=111)))
        proc = root / 'proc'; proc.mkdir()
        if active:
            folder = proc / '222'; folder.mkdir()
            (folder / 'stat').write_text('222 (python3) S 1 2')
            (folder / 'cmdline').write_bytes(b'python3\0/private/continue_packet_views_fast_v2.py\0--execute\0')
        output = root / 'outputs/cohort-sample-views-v1'; output.mkdir(parents=True)
        return proc, output

    def test_vcf_first_waiter_supersedes_stopped_old_views_controller(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, _ = self.views_fixture(root)
            result = m.snapshot(root, proc=proc)['deliverables']['cohort_and_sample_views']
            self.assertEqual(result['status'], 'WAITING_FOR_COMPLETED_INPUTS')
            self.assertTrue(result['controller']['state_pid_matches_active'])
            self.assertEqual(result['historical_views_v2']['status'], 'STOPPED_OLD_VIEWS_RULE')

    def test_guarded_v2_1_views_controller_is_current_same_authority(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, _ = self.views_fixture(root)
            (proc / '222/cmdline').write_bytes(b'python3\0/private/continue_packet_views_fast_v2_1.py\0--execute\0')
            result = m.snapshot(root, proc=proc)['deliverables']['cohort_and_sample_views']
            self.assertEqual(result['status'], 'WAITING_FOR_COMPLETED_INPUTS')
            self.assertTrue(result['controller']['state_pid_matches_active'])
            self.assertEqual(result['controller']['active_pids'], [222])

    def test_v2_2_controller_process_registration_release_is_observed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, _ = self.views_fixture(root)
            (proc / '222/cmdline').write_bytes(b'python3\0/private/continue_packet_views_fast_v2_2.py\0--execute\0')
            result = m.snapshot(root, proc=proc)['deliverables']['cohort_and_sample_views']
            self.assertEqual(result['status'], 'WAITING_FOR_COMPLETED_INPUTS')
            self.assertTrue(result['controller']['state_pid_matches_active'])
            self.assertEqual(result['controller']['active_pids'], [222])

    def test_v2_3_bounded_queue_controller_release_is_observed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, _ = self.views_fixture(root)
            (proc / '222/cmdline').write_bytes(b'python3\0/private/continue_packet_views_fast_v2_3.py\0--execute\0')
            result = m.snapshot(root, proc=proc)['deliverables']['cohort_and_sample_views']
            self.assertEqual(result['status'], 'WAITING_FOR_COMPLETED_INPUTS')
            self.assertTrue(result['controller']['state_pid_matches_active'])
            self.assertEqual(result['controller']['active_pids'], [222])

    def test_old_views_pass_is_not_vcf_first_completion_and_dead_waiter_is_not_active(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, output = self.views_fixture(root, active=False)
            (output / 'manifest.json').write_text(json.dumps(dict(status='PASS', report_input_kind='COMBINED_ANNOTATION_MASTER_V1')))
            result = m.snapshot(root, proc=proc)['deliverables']['cohort_and_sample_views']
            self.assertEqual(result['status'], 'HISTORICAL_RESULT_NOT_VCF_FIRST')
            self.assertEqual(result['controller']['status'], 'CONTROLLER_NOT_RUNNING_LAST_STATUS_OBSERVED')
            self.assertTrue(result['historical_result_receipt_observed'])

    def test_vcf_first_views_completion_requires_sample_scope_and_vcf_index_pair(self):
        variants = [(True, True, True), (False, True, False), (True, False, False)]
        for sample_scope, has_index, complete in variants:
            with self.subTest(sample_scope=sample_scope, has_index=has_index), tempfile.TemporaryDirectory() as temp:
                root = Path(temp); proc, output = self.views_fixture(root, active=False)
                artifacts = {'vcf/cohort.fully-annotated.vcf.gz': {'bytes': 20, 'sha256': 'a' * 64}}
                if has_index:artifacts['vcf/cohort.fully-annotated.vcf.gz.tbi'] = {'bytes': 10, 'sha256': 'b' * 64}
                (output / 'manifest.json').write_text(json.dumps(dict(status='PASS',
                    report_input_kind='FULLY_ANNOTATED_VCF_PACKETS_V1', report_producer_kind='FAST_DIRECT_PACKET_VIEWS_V2',
                    public_header_bundles_verified=True, no_legacy_receipt_created=True, sample_packet_scope_verified=sample_scope,
                    artifacts=artifacts)))
                result = m.snapshot(root, proc=proc)['deliverables']['cohort_and_sample_views']
                self.assertEqual(result['status'] == 'COMPLETE_RECEIPT_OBSERVED', complete)

    def test_vcf_first_monitor_never_opens_vcf_bam_or_packet_payloads(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); proc, output = self.views_fixture(root)
            (output / 'cohort.vcf.gz').write_bytes(b'PATIENT_PACKET_MUST_NOT_BE_READ')
            (output / 'sample.bam').write_bytes(b'PATIENT_BAM_MUST_NOT_BE_READ')
            original = Path.open
            def permitted(path, *args, **kwargs):
                self.assertNotIn(path.suffix, ('.gz', '.bam'))
                return original(path, *args, **kwargs)
            with patch.object(Path, 'open', permitted):
                result = m.snapshot(root, proc=proc)
            self.assertNotIn('PATIENT_PACKET', json.dumps(result))
            self.assertTrue(result['no_hashing'])

    def test_fully_annotated_vcf_controller_waiting_and_embedding_are_observed(self):
        for phase, expected in [('WAITING_CORE33_AND_SPLICEAI_PASS', 'WAITING_FOR_COMPLETED_INPUTS'),
                                ('APPENDING_ACTUAL_SPLICEAI_ALLELE_INFO_CORE_PACKETS_UNCHANGED', 'CONTROLLER_RUNNING_RESULT_NOT_COMPLETE')]:
            with self.subTest(phase=phase), tempfile.TemporaryDirectory() as temp:
                root = Path(temp); state = root / 'private/spliceai-allele-overlay-state-v1'; state.mkdir(parents=True)
                (state / 'status.json').write_text(json.dumps(dict(status=phase, pid=223, rule_sha256='c' * 64)))
                proc = root / 'proc'; folder = proc / '223'; folder.mkdir(parents=True)
                (folder / 'stat').write_text('223 (python3) S 1 2')
                (folder / 'cmdline').write_bytes(b'python3\0/private/continue_spliceai_allele_overlay.py\0run\0--execute\0')
                result = m.snapshot(root, proc=proc)['deliverables']['fully_annotated_vcf']
                self.assertEqual(result['status'], expected)

    def test_fully_annotated_vcf_completion_requires_own_receipt_and_index_artifact(self):
        for receipt_kind, expected in [('missing', False), ('no_index', False), ('packet_unverified', False), ('complete', True)]:
            with self.subTest(receipt_kind=receipt_kind), tempfile.TemporaryDirectory() as temp:
                root = Path(temp); output = root / 'outputs/fully-annotated-vcf-v1'; output.mkdir(parents=True)
                (output / 'manifest.json').write_text(json.dumps(dict(status='PASS_FULLY_ANNOTATED_VCF_DIRECT', final_delivery_complete=True,
                    original_records_fields_preserved=True, packet_decode_verified=(receipt_kind != 'packet_unverified'),
                    direct_core_prefix_verified=True, core33_packets_byte_identical=True, original_headers_preserved=True,
                    projection={'projection_contract': 'ALLELE_INFO_V1_ROW_PROJECTION_ON_READ', 'projection_fields': ['field'] * 13,
                                'physical_core_packets_unchanged': True})))
                artifacts = {'cohort.fully-annotated.vcf.gz': {'bytes': 20, 'sha256': 'a' * 64}}
                if receipt_kind in ('complete', 'packet_unverified'):
                    artifacts['cohort.fully-annotated.vcf.gz.tbi'] = {'bytes': 10, 'sha256': 'b' * 64}
                if receipt_kind != 'missing':
                    (output / 'workflow-receipt.json').write_text(json.dumps(dict(status='PASS',
                        stage='fully-annotated-vcf-v1', identity='d' * 64, artifacts=artifacts)))
                result = m.snapshot(root, proc=root / 'proc')['deliverables']['fully_annotated_vcf']
                self.assertEqual(result['status'] == 'COMPLETE_RECEIPT_OBSERVED', expected)

    def parallel_fixture(self, root, mode='run', execute=True):
        private = root / 'private'; private.mkdir()
        candidate = private / 'annotation-parallel-candidate'; (candidate / 'scripts').mkdir(parents=True)
        producer = candidate / 'scripts/parallel_annotation.py'; producer.write_text('# synthetic')
        rule = candidate / 'rule.json'
        value = dict(schema_version=1, status='CANDIDATE_USER_LAUNCH_ONLY', producer_path=str(producer),
                     producer_sha256='b' * 64, output_root=str(root / 'outputs/workflow-v1'),
                     workers=5, cpus_per_worker=4, memory_gib_per_worker=16,
                     deadline_utc='2030-01-01T00:00:00+00:00')
        rule.write_text(json.dumps(value))
        proc = root / 'proc'; folder = proc / '123'; folder.mkdir(parents=True)
        (folder / 'stat').write_text('123 (python3) S 1 2')
        args = ['python3', str(producer), mode, '--rule', str(rule), '--reviewed-rule-sha256', 'a' * 64]
        if execute:args.append('--execute')
        (folder / 'cmdline').write_bytes(b'\0'.join(a.encode() for a in args) + b'\0')
        queue = private / 'annotation-public-hcm-v2-review'; queue.mkdir()
        (queue / 'queue-status.json').write_text(json.dumps({'status': 'FAILED_OR_BUDGET_STOPPED'}))
        return private, proc, rule

    def test_live_parallel_primary_supersedes_cancelled_serial_queue(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); private, proc, _ = self.parallel_fixture(root)
            result = m.annotation_status(private, proc)
            self.assertEqual(result['status'], 'PARALLEL_ANNOTATION_RUNNING')
            self.assertEqual(result['source'], 'active_parallel_primary_runtime')
            self.assertEqual(result['historical_reference_v2_queue_status'], 'FAILED_OR_BUDGET_STOPPED')
            self.assertEqual(result['parallel_processes'][0]['pid'], 123)
            self.assertFalse(result['parallel_processes'][0]['rule_pins_rehashed'])
            self.assertNotIn('argv', json.dumps(result))

    def test_proof_dryrun_bad_rule_and_zombie_do_not_claim_parallel_running(self):
        for mode, execute, change in [('prove', True, ''), ('run', False, ''), ('run', True, 'rule'), ('run', True, 'zombie')]:
            with self.subTest(mode=mode, execute=execute, change=change), tempfile.TemporaryDirectory() as temp:
                root = Path(temp); private, proc, rule = self.parallel_fixture(root, mode, execute)
                if change == 'rule':
                    value = json.loads(rule.read_text()); value['output_root'] = '/unrelated'; rule.write_text(json.dumps(value))
                if change == 'zombie':(proc / '123/stat').write_text('123 (python3) Z 1 2')
                self.assertNotEqual(m.annotation_status(private, proc).get('source'), 'active_parallel_primary_runtime')

    def test_completed_primary_report_receipts_override_old_queue_without_live_pid(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); private, proc, _ = self.parallel_fixture(root)
            (proc / '123/stat').write_text('123 (python3) Z 1 2')
            for name in ['06_annotate', '07_report']:
                path = root / 'outputs/workflow-v1/stages' / name; path.mkdir(parents=True)
                (path / 'workflow-receipt.json').write_text(json.dumps(dict(status='PASS', stage=name,
                    identity='c' * 64, artifacts={'example.tsv': {'bytes': 10, 'sha256': 'd' * 64}})))
            result = m.annotation_status(private, proc)
            self.assertEqual(result['source'], 'primary_stage_07_receipt')
            self.assertEqual(result['status'], 'RESEARCH_REPORT_COMPLETE_RECEIPT_OBSERVED')
            self.assertEqual(result['parallel_processes'], [])
            self.assertIn('not rehashed', result['completion_basis'])

    def test_missing_outputs_are_not_reported_as_configured_or_complete(self):
        with tempfile.TemporaryDirectory() as temp:
            result = m.snapshot(Path(temp), proc=Path(temp) / 'proc')
            self.assertEqual(result['deliverables']['combined_master']['status'], 'NOT_CREATED_OR_WAITING')
            self.assertIn('ещё не созданы', m.human(result))

    def test_human_is_opt_in_and_default_stdout_remains_json(self):
        args = argparse.Namespace(max_hours=1, interval=60, root=Path('/missing'), local=False, watch=False, report_dir=None)
        for human in [False, True]:
            args.human = human
            with contextlib.redirect_stdout(io.StringIO()) as output:
                m.run(args, collect=lambda *_: {'annotation': {'status': 'PARALLEL_ANNOTATION_RUNNING'}})
            if human:self.assertIn('выполняется параллельно', output.getvalue())
            else:self.assertEqual(json.loads(output.getvalue())['annotation']['status'], 'PARALLEL_ANNOTATION_RUNNING')

    def test_resource_recovery_is_distinct_from_historical_installer_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); private = root / 'private'; private.mkdir()
            (private / 'databases-supervisor-status.json').write_text(json.dumps({'status': 'NONRETRYABLE_FAILURE_PARTIALS_PRESERVED'}))
            folder = private / 'annotation-resource-recovery-v1'; folder.mkdir()
            (folder / 'status.json').write_text(json.dumps({'status': 'VERIFYING_ALL_EXTRACTED_FILES_AGAINST_ARCHIVE', 'pid': 999, 'verified_files': 250}))
            result = m.snapshot(root, proc=root / 'no-proc')['databases']
            self.assertEqual(result['historical_installer_status'], 'NONRETRYABLE_FAILURE_PARTIALS_PRESERVED')
            self.assertEqual(result['status'], 'VERIFYING_ALL_EXTRACTED_FILES_AGAINST_ARCHIVE')
            self.assertEqual(result['recovery']['verified_files'], 250)
            self.assertFalse(result['recovery']['process']['verified_active'])

    def test_reviewed_controller_supersedes_stopped_historical_queue(self):
        with tempfile.TemporaryDirectory() as temp:
            private = Path(temp)
            (private / 'annotation-queue-status.json').write_text(json.dumps({'status': 'STOPPED_VM_ALLOWANCE_OR_CANCEL'}))
            folder = private / 'annotation-reference-controller-v1'; folder.mkdir()
            (folder / 'controller-status.json').write_text(json.dumps({'status': 'WAITING_ORIGINAL_VEP_CORE_PASS_NO_APPROVAL_NO_QUEUE', 'rule_sha256': 'a' * 64, 'pid': 123}))
            result = m.annotation_status(private, private / 'no-proc')
            self.assertEqual(result['source'], 'reviewed_reference_v2_controller')
            self.assertEqual(result['status'], 'WAITING_ORIGINAL_VEP_CORE_PASS_NO_APPROVAL_NO_QUEUE')
            self.assertEqual(result['historical_queue_status'], 'STOPPED_VM_ALLOWANCE_OR_CANCEL')
            self.assertFalse(result['controller']['verified_active'])
            queue = private / 'annotation-public-hcm-v2-review'; queue.mkdir()
            (queue / 'queue-status.json').write_text(json.dumps({'status': 'PUBLIC_CONTROLS_PASS_PATIENT_ANNOTATION_RUNNING'}))
            result = m.annotation_status(private, private / 'no-proc')
            self.assertEqual(result['source'], 'reviewed_reference_v2_queue')
            self.assertEqual(result['status'], 'PUBLIC_CONTROLS_PASS_PATIENT_ANNOTATION_RUNNING')

    def test_invalid_controller_receipt_does_not_hide_original_queue(self):
        with tempfile.TemporaryDirectory() as temp:
            private = Path(temp)
            (private / 'annotation-queue-status.json').write_text(json.dumps({'status': 'WAITING_VEP_RESOURCES'}))
            folder = private / 'annotation-reference-controller-v1'; folder.mkdir()
            (folder / 'controller-status.json').write_text('{"status":"COMPLETE","rule_sha256":"invalid"}')
            self.assertEqual(m.annotation_status(private, private / 'no-proc')['source'], 'original_queue')

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
