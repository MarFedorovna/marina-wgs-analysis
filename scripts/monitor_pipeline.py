#!/usr/bin/env python3
"""Bounded, read-only pipeline observer. Never starts/stops jobs or hashes data."""
import argparse
import datetime as dt
import html
import json
import math
from pathlib import Path
import re
import time

MAX_JSON = 2 * 1024 * 1024


def load(path):
    try:
        if path.is_symlink() or path.stat().st_size > MAX_JSON:
            return {'read_error': 'symlink_or_oversized_receipt'}
        return json.loads(path.read_text(encoding='utf-8'))
    except FileNotFoundError:
        return {}
    except (OSError, ValueError):
        return {'read_error': 'unreadable_or_changing_receipt'}


def fields(value, *names):
    return {key: value[key] for key in names if key in value}


def parallel_processes(private, proc):
    """Observe an executing private primary wrapper, never a proof or dry run.

    Small rule fields bind process purpose/output. Pins are reported, not rehashed;
    this observer does not replace the launcher's exact-rule verification.
    """
    result = []
    for folder in proc.glob('[0-9]*'):
        try:
            args = [v.decode() for v in (folder / 'cmdline').read_bytes().split(b'\0') if v]
            scripts = [v for v in args if Path(v).name == 'parallel_annotation.py']
            if len(scripts) != 1 or 'run' not in args or '--execute' not in args:
                continue
            if args.count('--rule') != 1 or args.count('--reviewed-rule-sha256') != 1:
                continue
            rule_path = Path(args[args.index('--rule') + 1]).resolve(strict=True)
            producer = Path(scripts[0]).resolve(strict=True)
            if not all(p.is_relative_to(private.resolve()) for p in (rule_path, producer)):
                continue
            rule = load(rule_path)
            reviewed_sha = args[args.index('--reviewed-rule-sha256') + 1]
            if (rule.get('schema_version') != 1 or rule.get('status') != 'CANDIDATE_USER_LAUNCH_ONLY'
                or Path(rule.get('producer_path', '')).resolve() != producer
                or Path(rule.get('output_root', '')).resolve() != (private.parent / 'outputs/workflow-v1').resolve()
                or not re.fullmatch(r'[0-9a-f]{64}', reviewed_sha)
                or not re.fullmatch(r'[0-9a-f]{64}', str(rule.get('producer_sha256', '')))
                or rule.get('workers') != 5 or rule.get('cpus_per_worker') != 4
                or rule.get('memory_gib_per_worker') != 16):
                continue
            item = process(folder.name, scripts[0], proc)
            if item['verified_active']:
                result.append(dict(pid=item['pid'], verified_active=True, workers_max=5,
                    cpus_max=20, memory_gib_max=80, reviewed_rule_sha256=reviewed_sha,
                    deadline_utc=rule.get('deadline_utc'), rule_pins_rehashed=False))
        except (OSError, ValueError, UnicodeError, IndexError):
            continue
    return sorted(result, key=lambda item: item['pid'])


def completed_stage(root, name):
    receipt = load(root / 'outputs/workflow-v1/stages' / name / 'workflow-receipt.json')
    return (receipt.get('status') == 'PASS' and receipt.get('stage') == name
            and bool(re.fullmatch(r'[0-9a-f]{64}', str(receipt.get('identity', ''))))
            and isinstance(receipt.get('artifacts'), dict) and bool(receipt['artifacts']))


def annotation_status(private, proc):
    historical = load(private / 'annotation-queue-status.json')
    controller = load(private / 'annotation-reference-controller-v1/controller-status.json')
    current = load(private / 'annotation-public-hcm-v2-review/queue-status.json')
    valid_controller = isinstance(controller.get('rule_sha256'), str) and bool(
        re.fullmatch(r'[0-9a-f]{64}', controller['rule_sha256'])) and bool(controller.get('status'))
    if valid_controller:
        selected = current if current.get('status') else controller
        result = fields(selected, 'status', 'updated_utc', 'child_pid', 'read_error')
        result['source'] = 'reviewed_reference_v2_queue' if current.get('status') else 'reviewed_reference_v2_controller'
        result['controller'] = process(controller.get('pid'), 'continue_annotation_reference_v2.py', proc)
        result['controller_status'] = controller['status']
        result['historical_queue_status'] = historical.get('status')
    else:
        result = fields(historical, 'status', 'updated_utc', 'child_pid', 'read_error')
        result['source'] = 'original_queue'
    result['queue_pids'] = discover('annotation_queue.py run', proc)
    result['public_controls_then_stages_6_7'] = 'CONSULT_CURRENT_QUEUE_RECEIPT'
    parallel = parallel_processes(private, proc)
    stage6 = completed_stage(private.parent, '06_annotate')
    stage7 = completed_stage(private.parent, '07_report')
    if parallel or stage6 or stage7:
        result['historical_reference_v2_queue_status'] = current.get('status')
        result['historical_original_queue_status'] = historical.get('status')
        result['parallel_processes'] = parallel
        result['stage06_receipt_observed_PASS'] = stage6
        result['stage07_receipt_observed_PASS'] = stage7
        result['completion_basis'] = 'small primary receipts observed; artifacts not rehashed'
        result['public_controls_then_stages_6_7'] = 'PRIMARY_STAGE_RECEIPTS_AND_CURRENT_PRODUCER'
        if stage6 and stage7:
            result.update(source='primary_stage_07_receipt', status='RESEARCH_REPORT_COMPLETE_RECEIPT_OBSERVED')
        elif stage6:
            result.update(source='primary_stage_06_receipt', status=(
                'ANNOTATION_PASS_REPORT_RUNNING' if parallel else 'ANNOTATION_PASS_REPORT_PENDING_OR_NOT_OBSERVED'))
        elif parallel:
            result.update(source='active_parallel_primary_runtime', status='PARALLEL_ANNOTATION_RUNNING')
        result.pop('child_pid', None)
    return result


def process(pid, expected, proc=Path('/proc')):
    """Check PID identity instead of treating an old log/PID as a live job."""
    try:
        pid = int(pid)
        if pid <= 0:
            raise ValueError()
        folder = proc / str(pid)
        stat = (folder / 'stat').read_text()
        state = stat[stat.rfind(')') + 2:].split()[0]
        command = (folder / 'cmdline').read_bytes().replace(b'\0', b' ').decode(errors='replace')
        match = expected in command
        return dict(pid=pid, alive=state not in ('Z', 'X'), identity_matches=match,
                    verified_active=state not in ('Z', 'X') and match, state=state)
    except (OSError, ValueError, IndexError, TypeError):
        return dict(pid=pid, alive=False, identity_matches=False, verified_active=False)


def discover(expected, proc=Path('/proc')):
    # Read-only proc inventory; never prints command lines or private arguments.
    result = []
    for folder in proc.glob('[0-9]*'):
        item = process(folder.name, expected, proc)
        if item['verified_active']:
            result.append(item['pid'])
    return sorted(result)


def executing_controllers(script, proc):
    """Exact script basename and execution flag; never expose private argv."""
    result = []
    allowed = {script} if isinstance(script, str) else set(script)
    for folder in proc.glob('[0-9]*'):
        try:
            args = [a.decode() for a in (folder / 'cmdline').read_bytes().split(b'\0') if a]
            matching = [Path(a).name for a in args if Path(a).name in allowed]
            if '--execute' in args and len(matching) == 1:
                item = process(folder.name, matching[0], proc)
                if item['verified_active']:
                    result.append(item['pid'])
        except (OSError, UnicodeError):
            continue
    return sorted(result)


def controller_status(path, script, proc):
    state = load(path)
    pids = executing_controllers(script, proc)
    bound = (state.get('pid') in pids and bool(re.fullmatch(r'[0-9a-f]{64}', str(state.get('rule_sha256', '')))))
    # Fork workers inherit the parent's command line. Only the PID that owns
    # this state file is the controller, not every process with that argv.
    if bound:
        pids = [state['pid']]
    result = dict(active_pids=pids, state_pid_matches_active=bound,
                  last_observed=fields(state, 'status', 'updated_utc', 'deadline_utc', 'error_type'),
                  rule_pins_rehashed=False)
    if bound:
        result['status'] = state.get('status', 'ACTIVE_CONTROLLER_STATUS_NOT_OBSERVED')
    elif pids:
        result['status'] = 'ACTIVE_CONTROLLER_STATE_NOT_BOUND'
    elif state:
        result['status'] = 'CONTROLLER_NOT_RUNNING_LAST_STATUS_OBSERVED'
    else:
        result['status'] = 'CONTROLLER_NOT_OBSERVED'
    return result


def batch_progress(state_directory, state, *, resume=False):
    works = [p for p in state_directory.glob('work-*') if p.is_dir() and not p.is_symlink()]
    result = dict(basis='completed batch row-count receipts only; not elapsed-work or ETA proof',
                  payloads_read=False, artifacts_rehashed=False)
    if len(works) != 1:
        result['status'] = 'WORK_DIRECTORY_NOT_UNAMBIGUOUS'
        return result
    receipts, seen = [], set()
    if resume:
        audit = load(works[0] / 'resume-audit.json')
        origins = audit.get('source_origin_mapping')
        if not isinstance(origins, list):
            return dict(result, status='RESUME_CACHE_METADATA_NOT_OBSERVED')
        for item in origins:
            if item.get('reuse') is not True: continue
            job, ids = item.get('old_job', {}), item.get('new_batch_ids', [])
            if (len(ids) != 1 or type(ids[0]) is not int or not all(type(job.get(k)) is int
                for k in ('first_i','last_i','expected_rows','records')) or job['expected_rows'] <= 0
                or job['last_i'] - job['first_i'] + 1 != job['expected_rows'] or ids[0] in seen):
                return dict(result, status='INCONSISTENT_RESUME_CACHE_METADATA')
            seen.add(ids[0]); receipts.append(dict(id=ids[0], first_i=job['first_i'], last_i=job['last_i'],
                row_count=job['expected_rows'], record_count=job['records']))
        reused_rows = sum(v['row_count'] for v in receipts)
        if (len(receipts) != audit.get('reused_batches') or reused_rows != audit.get('reused_rows')
            or (state.get('reused_rows') is not None and reused_rows != state['reused_rows'])
            or (state.get('reused_batches') is not None and len(receipts) != state['reused_batches'])):
            return dict(result, status='INCONSISTENT_RESUME_CACHE_METADATA')
        result.update(reused_batches_observed=len(receipts), reused_rows_observed=reused_rows,
            cache_basis='reviewed producer resume audit metadata observed; cached bodies not reopened or rehashed')
    for path in works[0].glob('result-*/batch-proof.json'):
        value = load(path)
        names = ('id', 'first_i', 'last_i', 'row_count', 'record_count')
        if not all(type(value.get(k)) is int for k in names):
            result['status'] = 'INCOMPLETE_OR_CHANGING_BATCH_RECEIPT'
            continue
        if value['id'] in seen or value['row_count'] <= 0 or value['last_i'] - value['first_i'] + 1 != value['row_count']:
            return dict(result, status='INCONSISTENT_BATCH_METADATA')
        seen.add(value['id']); receipts.append(value)
    ordered = sorted(receipts, key=lambda v: v['first_i'])
    if any(a['last_i'] >= b['first_i'] for a, b in zip(ordered, ordered[1:])):
        return dict(result, status='OVERLAPPING_BATCH_METADATA')
    rows = sum(v['row_count'] for v in receipts); total = state.get('expected_rows', state.get('rows'))
    expected_batches = state.get('batches')
    if resume and expected_batches is None and type(audit.get('regenerated_tasks')) is int:
        expected_batches = audit['reused_batches'] + audit['regenerated_tasks']
    result.update(status='COMPLETED_BATCH_METADATA_OBSERVED', completed_batches=len(receipts),
        completed_rows=rows, completed_records=sum(v['record_count'] for v in receipts),
        expected_rows=total, expected_batches=expected_batches)
    if type(total) is int and total > 0 and rows <= total:
        result['row_progress_percent'] = round(100 * rows / total, 2)
    return result


def direct_core_controller(private, proc):
    profiles = []
    allowed = {'direct_fully_annotated_vcf.py', 'direct_fully_annotated_vcf_v2.py', 'direct_fully_annotated_vcf_v3.py'}
    for folder in proc.glob('[0-9]*'):
        try:
            args = [a.decode() for a in (folder / 'cmdline').read_bytes().split(b'\0') if a]
            scripts = [a for a in args if Path(a).name in allowed]
            if len(scripts) != 1 or 'run' not in args or '--execute' not in args or args.count('--rule') != 1:
                continue
            path = Path(args[args.index('--rule') + 1]).resolve(strict=True)
            if not path.is_relative_to(private.resolve()) or args.count('--reviewed-rule-sha256') != 1:
                continue
            rule = load(path); reviewed = args[args.index('--reviewed-rule-sha256') + 1]
            state_directory = Path(rule.get('state_directory', '')).resolve()
            if (rule.get('root_reviewed') is not True or rule.get('status') not in
                ('ROOT_REVIEWED_DIRECT_CORE33_RULE', 'ROOT_REVIEWED_DIRECT_CORE33_V2_RULE')
                or not re.fullmatch('[0-9a-f]{64}', reviewed)
                or Path(rule.get('driver', {}).get('path', '')).resolve() != Path(scripts[0]).resolve()
                or not state_directory.is_relative_to(private.resolve())
                or Path(rule.get('output', '')).resolve() != (private.parent / 'outputs/direct-core-context-v1').resolve()):
                continue
            state = load(state_directory / 'status.json')
            # This also excludes Pool workers with inherited identical argv.
            if state.get('pid') != int(folder.name) or not process(folder.name, scripts[0], proc)['verified_active']:
                continue
            profile = dict(status=state.get('status', 'DIRECT_STATE_NOT_OBSERVED'), active_pids=[int(folder.name)],
                state_pid_matches_active=True, reviewed_rule_sha256_observed=reviewed, rule_pins_rehashed=False,
                workers=rule.get('workers'), worker_cpu_ids=rule.get('worker_cpu_ids'),
                last_observed=fields(state, 'status', 'updated_utc', 'deadline_utc', 'expected_rows', 'batches', 'workers'),
                progress=batch_progress(state_directory, state, resume=Path(scripts[0]).name == 'direct_fully_annotated_vcf_v3.py'))
            profiles.append(profile)
        except (OSError, ValueError, UnicodeError, IndexError, TypeError):
            continue
    if len(profiles) == 1:
        return profiles[0]
    if len(profiles) > 1:
        return dict(status='MULTIPLE_DIRECT_CONTROLLERS_OBSERVED', active_pids=[p['active_pids'][0] for p in profiles],
                    state_pid_matches_active=False, rule_pins_rehashed=False)
    return dict(status='DIRECT_CONTROLLER_NOT_OBSERVED', active_pids=[], state_pid_matches_active=False,
                rule_pins_rehashed=False)


def controlled_delivery(path, expected, controller, required_fields=None, required_artifacts=(),
                        workflow_stage=None, workflow_artifacts=(), mismatch_is_historical=False):
    item = deliverable(path, expected)
    manifest = load(path)
    matches = all(manifest.get(k) == v for k, v in (required_fields or {}).items())
    artifacts = manifest.get('artifacts', {})
    matches = matches and isinstance(artifacts, dict) and set(required_artifacts).issubset(artifacts)
    if workflow_stage:
        receipt = load(path.parent / 'workflow-receipt.json')
        artifacts = receipt.get('artifacts', {})
        ready = (receipt.get('status') == 'PASS' and receipt.get('stage') == workflow_stage
                 and bool(re.fullmatch(r'[0-9a-f]{64}', str(receipt.get('identity', ''))))
                 and isinstance(artifacts, dict) and set(workflow_artifacts).issubset(artifacts))
        item['workflow_receipt_observed_PASS'] = ready
        matches = matches and ready
    if item['status'] == 'COMPLETE_RECEIPT_OBSERVED' and not matches:
        item['status'] = 'HISTORICAL_RESULT_NOT_VCF_FIRST' if mismatch_is_historical else 'NO_MATCHING_PASS_RECEIPT'
        if mismatch_is_historical:
            item['historical_result_receipt_observed'] = True
    item['controller'] = controller
    if item['status'] != 'COMPLETE_RECEIPT_OBSERVED' and controller['active_pids']:
        state = controller['status']
        if state.startswith('WAITING_'):
            item['status'] = 'WAITING_FOR_COMPLETED_INPUTS'
        elif state.startswith(('FAIL', 'REFUSED', 'STOPPED')):
            item['status'] = 'CONTROLLER_REPORTED_FAILURE_OR_STOP'
        else:
            item['status'] = 'CONTROLLER_RUNNING_RESULT_NOT_COMPLETE'
    return item


def range_progress(receipt, partial):
    identity = receipt.get('identity', {})
    total, chunk = identity.get('bytes'), identity.get('chunk_bytes')
    result = {'verified_bytes': 0, 'verified_ranges': 0, 'total_bytes': total,
              'basis': 'installer SHA range receipts; no rehash; sparse size is not progress'}
    if isinstance(total, int) and isinstance(chunk, int) and total > 0 and chunk > 0:
        completed = receipt.get('completed', {})
        if isinstance(completed, dict):
            for key, digest in completed.items():
                if not str(key).isdigit() or not isinstance(digest, str) or not re.fullmatch(r'[0-9a-fA-F]{64}', digest):
                    continue
                index = int(key)
                if index * chunk < total:
                    result['verified_ranges'] += 1
                    result['verified_bytes'] += min(chunk, total - index * chunk)
        result['percent'] = round(100 * result['verified_bytes'] / total, 2)
    try:
        info = partial.stat()
        result['logical_bytes_not_progress'] = info.st_size
        result['allocated_bytes_not_progress'] = getattr(info, 'st_blocks', 0) * 512
    except OSError:
        pass
    if receipt.get('read_error'):
        result['read_error'] = receipt['read_error']
    return result


def current_attempt_errors(path, offset):
    """Bounded tail after supervisor's attempt offset; historical failures excluded."""
    try:
        if path.is_symlink() or not isinstance(offset, int) or offset < 0:
            return {'status': 'NO_CURRENT_ATTEMPT_OFFSET'}
        size = path.stat().st_size
        if offset > size:
            return {'status': 'LOG_ROTATED_OR_OFFSET_UNAVAILABLE'}
        with path.open('rb') as stream:
            stream.seek(max(offset, size - 16384))
            tail = stream.read(16384).decode(errors='replace').lower()
        labels = {'connection_refused': 'connection refused', 'timeout': 'timed out',
                  'checksum_error': 'checksum mismatch', 'traceback': 'traceback',
                  'permission_denied': 'permission denied'}
        return {'status': 'BOUNDED_CURRENT_ATTEMPT_TAIL', 'tail_bytes_max': 16384,
                'signals': [name for name, text in labels.items() if text in tail],
                'note': 'Signals are context, not job failure; use live process and current receipt.'}
    except OSError:
        return {'status': 'LOG_UNAVAILABLE'}


def deliverable(path, expected):
    value = load(path)
    if value.get('read_error'):
        status = 'UNREADABLE_COMPLETION_RECEIPT'
    elif value.get('status') == expected:
        status = 'COMPLETE_RECEIPT_OBSERVED'
    elif value:
        status = 'NO_MATCHING_PASS_RECEIPT'
    else:
        status = 'NOT_CREATED_OR_WAITING'
    return dict(status=status, receipt_status=value.get('status'), artifacts_rehashed=False)


def snapshot(root, local=False, proc=Path('/proc')):
    now = dt.datetime.now(dt.timezone.utc)
    result = dict(utc=now.isoformat(), mode='LOCAL_RECEIPTS_ONLY' if local else 'SERVER_READ_ONLY',
                  no_process_actuation=True, no_patient_analysis=True, no_hashing=True)
    if local:
        # Deliberately never read native-sftp-direct-config.json, SSH keys, BAMs or VCFs.
        transfer = load(root / '.local/native-sftp-direct-status.json')
        result['local_transfer'] = fields(transfer, 'status', 'updated_utc', 'completed_files',
                                         'received_bytes', 'total_bytes', 'bytes_downloaded', 'elapsed_seconds',
                                         'Status', 'UpdatedUtc', 'FileIndex', 'FileCount', 'CompletedBytes',
                                         'ExpectedBytes', 'MiBPerSecond', 'PatientHashesComputed')
        return result
    private = root / 'private'
    db = load(private / 'databases-supervisor-status.json')
    result['databases'] = fields(db, 'status', 'attempts_started', 'updated_utc', 'deadline_utc',
                                'retry_policy', 'next_retry_delay_seconds', 'read_error')
    result['databases']['installer'] = process(db.get('current_pid'), 'prepare_annotation_resources.py', proc)
    result['databases']['supervisor_pids'] = discover('supervise_databases.py', proc)
    result['databases']['current_attempt_log'] = current_attempt_errors(root / 'logs/annotation-databases.log', db.get('attempt_log_offset'))
    recovery = load(private / 'annotation-resource-recovery-v1/status.json')
    if recovery.get('status') and isinstance(recovery.get('pid'), int):
        result['databases']['historical_installer_status'] = db.get('status')
        result['databases']['status'] = recovery['status']
        result['databases']['status_authority'] = 'reviewed_resource_recovery_receipt'
        result['databases']['recovery'] = fields(recovery, 'status', 'updated_utc',
            'verified_files', 'verified_expanded_bytes', 'error_type', 'deadline_utc')
        result['databases']['recovery']['process'] = process(recovery['pid'], 'recover_annotation_extracted_stage.py', proc)
    deadline = db.get('deadline_utc')
    try:
        seconds = (dt.datetime.fromisoformat(deadline.replace('Z', '+00:00')) - now).total_seconds()
        result['databases']['budget_deadline_reached'] = seconds <= 0
        result['databases']['seconds_until_original_deadline'] = seconds
    except (TypeError, ValueError, AttributeError):
        result['databases']['budget_deadline_reached'] = None
    base = root / 'resources/.downloads/homo_sapiens_vep_115_GRCh38.tar.gz'
    result['vep_cache'] = range_progress(load(base.with_name(base.name + '.agent-parts.json')),
                                         base.with_name(base.name + '.agent-partial'))
    mt = load(private / 'mitochondria-recovery-v2-queue-status.json')
    # The generic progress file is shared with the v2 worker: bind by policy version.
    progress = load(private / 'mitochondria-recovery-status.json')
    if progress.get('policy', {}).get('version') != 2:
        progress = {'status': 'UNBOUND_PROGRESS_NOT_USED'}
    result['mt_v2'] = fields(mt, 'status', 'updated_utc', 'returncode', 'samples', 'read_error')
    result['mt_v2']['controller'] = process(mt.get('pid'), 'run-mito-recovery-v2-queue.py', proc)
    result['mt_v2']['progress'] = fields(progress, 'status', 'sample_index', 'samples', 'completed')
    result['mt_v2']['worker'] = process(progress.get('pid'), 'recover_mitochondria.py', proc)
    audit = load(private / 'mitochondria-recovery-v2-technical-audit.json')
    result['mt_v2']['final_audit'] = fields(audit, 'status', 'samples', 'completed_samples', 'read_error') or {'status': 'NOT_YET_PRESENT'}
    waiter = load(private / 'mitochondria-recovery-v2-audit-status.json')
    result['mt_v2']['audit_waiter'] = fields(waiter, 'status', 'updated_utc')
    result['mt_v2']['audit_waiter']['process'] = process(waiter.get('pid'), 'run-mito-recovery-v2-audit.py', proc)
    result['coverage'] = fields(load(private / 'gene-coverage-status.json'), 'status', 'completed_samples', 'updated_utc', 'read_error')
    result['annotation'] = annotation_status(private, proc)
    result['annotation']['role'] = 'BASE_VEP_AND_STAGE07_CONTROL_REPORT_NOT_FINAL_VCF_FIRST_DELIVERY'
    result['deliverables'] = dict(combined_master=deliverable(
        root / 'outputs/combined-annotation-master-v1/manifest.json', 'PASS_COMBINED_ANNOTATION_MASTER'))
    direct = direct_core_controller(private, proc)
    result['deliverables']['direct_core_context'] = controlled_delivery(
        root / 'outputs/direct-core-context-v1/manifest.json', 'INTERMEDIATE_CORE33_PASS', direct,
        {'final_delivery_complete': False, 'original_records_fields_preserved': True,
         'packet_decode_verified': True, 'direct_core_prefix_verified': True},
        workflow_stage='direct-core-context-v1',
        workflow_artifacts=['cohort.core-context.vcf.gz', 'cohort.core-context.vcf.gz.tbi', 'core-manifest.json'])
    result['deliverables']['direct_core_context']['role'] = 'INTERMEDIATE_CORE33_NOT_FINAL_SPLICEAI_DELIVERY'
    result['deliverables']['direct_core_context']['historical_v1_superseded'] = fields(
        load(private / 'direct-core33-state-v1/ROOT-superseded.json'), 'status', 'new_pid', 'updated_utc')
    model = controller_status(private / 'spliceai-reference-context-state-v1/status.json',
                              'continue_spliceai_reference_context.py', proc)
    result['spliceai'] = controlled_delivery(root / 'outputs/spliceai-reference-context-v1/manifest.json',
        'PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT', model, workflow_stage='spliceai-reference-context-v1',
        workflow_artifacts=['predictions.tsv.gz'])
    result['spliceai']['progress'] = fields(load(private / 'spliceai-reference-context-state-v1/progress.json'),
        'phase', 'updated_utc', 'variants_completed', 'variants_estimated', 'inferences_this_run', 'elapsed_seconds')
    result['deliverables']['fully_annotated_vcf'] = controlled_delivery(
        root / 'outputs/fully-annotated-vcf-v1/manifest.json', 'PASS_FULLY_ANNOTATED_VCF_DIRECT',
        controller_status(private / 'spliceai-allele-overlay-state-v1/status.json', 'continue_spliceai_allele_overlay.py', proc),
        {'final_delivery_complete': True, 'original_records_fields_preserved': True, 'packet_decode_verified': True,
         'direct_core_prefix_verified': True, 'core33_packets_byte_identical': True, 'original_headers_preserved': True},
        workflow_stage='fully-annotated-vcf-v1', workflow_artifacts=['cohort.fully-annotated.vcf.gz', 'cohort.fully-annotated.vcf.gz.tbi'])
    result['deliverables']['fully_annotated_vcf']['progress'] = fields(
        load(private / 'spliceai-allele-overlay-state-v1/progress.json'), 'phase', 'updated_utc',
        'records', 'records_total', 'records_verified', 'matched_model_alleles')
    final_manifest = load(root / 'outputs/fully-annotated-vcf-v1/manifest.json')
    projection = final_manifest.get('projection', {})
    result['deliverables']['fully_annotated_vcf']['projection_contract_observed'] = projection.get('projection_contract')
    if result['deliverables']['fully_annotated_vcf']['status'] == 'COMPLETE_RECEIPT_OBSERVED':
        if (projection.get('projection_contract') != 'ALLELE_INFO_V1_ROW_PROJECTION_ON_READ'
            or not isinstance(projection.get('projection_fields'), list) or len(projection['projection_fields']) != 13
            or projection.get('physical_core_packets_unchanged') is not True):
            result['deliverables']['fully_annotated_vcf']['status'] = 'NO_MATCHING_PASS_RECEIPT'
    result['deliverables']['cohort_and_sample_views'] = controlled_delivery(
        root / 'outputs/cohort-sample-views-v1/manifest.json', 'PASS',
        controller_status(private / 'packet-fastviews-continuation-v5/status.json',
            ('continue_packet_views_fast_v2.py', 'continue_packet_views_fast_v2_1.py', 'continue_packet_views_fast_v2_2.py', 'continue_packet_views_fast_v2_3.py'), proc),
        {'report_input_kind': 'FULLY_ANNOTATED_VCF_PACKETS_V1', 'report_producer_kind': 'FAST_DIRECT_PACKET_VIEWS_V2',
         'sample_packet_scope_verified': True, 'public_header_bundles_verified': True, 'no_legacy_receipt_created': True},
        ['vcf/cohort.fully-annotated.vcf.gz', 'vcf/cohort.fully-annotated.vcf.gz.tbi'], mismatch_is_historical=True)
    result['deliverables']['cohort_and_sample_views']['historical_views_v3'] = fields(
        load(private / 'cohort-sample-views-continuation-v3/status.json'), 'status', 'updated_utc')
    result['deliverables']['cohort_and_sample_views']['historical_views_v2'] = fields(
        load(private / 'cohort-sample-views-continuation-v2/status.json'), 'status', 'updated_utc')
    result['deliverables']['cohort_and_sample_views']['historical_fastviews_v2'] = fields(
        load(private / 'packet-fastviews-continuation-v2/status.json'), 'status', 'updated_utc')
    if direct['active_pids']:
        result['annotation']['historical_base_vep_status'] = result['annotation'].get('status')
        result['annotation'].update(status='DIRECT_CORE33_RUNNING', source='root_reviewed_direct_core33', direct_controller=direct,
                                    role='CURRENT_DIRECT_CORE33_PRODUCER_FINAL_OVERLAY_AND_VIEWS_SEPARATE')
    elif result['deliverables']['direct_core_context']['status'] == 'COMPLETE_RECEIPT_OBSERVED':
        result['annotation'].update(status='DIRECT_CORE33_COMPLETE_WAITING_FINAL',source='direct_core33_receipt')
    if result['deliverables']['fully_annotated_vcf']['status'] == 'COMPLETE_RECEIPT_OBSERVED':
        result['annotation'].update(status='FULLY_ANNOTATED_VCF_COMPLETE_WAITING_VIEWS',source='final_vcf_receipt')
    if result['deliverables']['cohort_and_sample_views']['status'] == 'COMPLETE_RECEIPT_OBSERVED':
        result['annotation'].update(status='VCF_FIRST_RESULTS_COMPLETE',source='vcf_first_views_receipt')
    result['deliverables']['combined_master']['role'] = 'INTERMEDIATE_REFERENCE_JOIN_NOT_PRIMARY_RESULT'
    alpha = load(private / 'alphamissense-control-v1/public-controls-receipt.json')
    result['annotation']['separate_alphamissense_controls'] = fields(alpha, 'status', 'controls', 'tested')
    result['xy'] = dict(status='PLAN_ONLY_NO_EXECUTOR', clinical_sex_assignment=False,
                        somatic_analysis=False, ploidy='unknown_all_four_until_review',
                        par_alignment_contract='unverified_until_separate_review')
    result['limitations'] = ['No old MT v1 failure is used as current v2 status.',
        'Receipts are reported, not a fresh technical/clinical audit.',
        'Observed primary PASS receipts do not establish X/Y completion or clinical interpretation.']
    return result


def report(result, directory):
    # Only these dedicated observer files may be written. Existing streams stay untouched.
    if directory.name != 'monitor-reports' or directory.is_symlink():
        raise ValueError('Report directory must be a dedicated non-symlink monitor-reports directory.')
    directory.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(result, ensure_ascii=False, indent=2) + '\n'
    for name, content in [('status.json', payload), ('status.html', '<!doctype html><meta charset="utf-8"><title>Pipeline status</title><pre>' + html.escape(payload) + '</pre>')]:
        path = directory / name
        if path.is_symlink() or path.with_suffix(path.suffix + '.tmp').is_symlink():
            raise ValueError('Refusing report symlink.')
        temporary = path.with_suffix(path.suffix + '.tmp')
        temporary.write_text(content, encoding='utf-8')
        temporary.replace(path)


def human(result):
    if result.get('mode') == 'LOCAL_RECEIPTS_ONLY':
        item = result.get('local_transfer', {})
        return 'Локальная копия: ' + str(item.get('status', item.get('Status', 'нет квитанции')))
    annotation = result.get('annotation', {})
    state = annotation.get('status', 'нет текущей квитанции')
    names = {'PARALLEL_ANNOTATION_RUNNING': 'выполняется параллельно',
             'DIRECT_CORE33_RUNNING': 'VEP готов; формируются пакеты строк и33 дополнительных поля',
             'DIRECT_CORE33_COMPLETE_WAITING_FINAL': 'VEP и33 поля готовы; итоговый VCF ещё ожидается',
             'FULLY_ANNOTATED_VCF_COMPLETE_WAITING_VIEWS': 'итоговый VCF готов; таблицы ещё ожидаются',
             'VCF_FIRST_RESULTS_COMPLETE': 'итоговый VCF и таблицы готовы по квитанциям',
             'ANNOTATION_PASS_REPORT_RUNNING': 'готова; формируется отчёт',
             'ANNOTATION_PASS_REPORT_PENDING_OR_NOT_OBSERVED': 'готова; завершение отчёта ещё не подтверждено',
             'RESEARCH_REPORT_COMPLETE_RECEIPT_OBSERVED': 'базовая аннотация и контрольный отчёт этапа 07 готовы по квитанциям'}
    lines = ['Основная VEP-аннотация: ' + names.get(state, state)]
    active = annotation.get('parallel_processes', [])
    if active:
        lines.append('Основной процесс: PID ' + ', '.join(str(p['pid']) for p in active)
                     + '; до 5 задач, 20 CPU, 80 GiB.')
    if 'historical_reference_v2_queue_status' in annotation:
        lines.append('Прежняя очередь: ' + str(annotation['historical_reference_v2_queue_status'])
                     + ' — историческая запись, не статус текущего процесса.')
    deadline = next((p.get('deadline_utc') for p in active if p.get('deadline_utc')), None)
    if deadline:
        lines.append('Предельное время текущего плана (UTC): ' + deadline)
    direct = annotation.get('direct_controller', {})
    if direct.get('active_pids'):
        lines.append('Текущий прямой производитель: PID ' + ', '.join(map(str, direct['active_pids'])) + '; процессов ' + str(direct.get('workers')) + '.')
        progress = direct.get('progress', {})
        if 'completed_rows' in progress:
            lines.append('Завершённые блоки: ' + str(progress['completed_batches']) + '/' + str(progress.get('expected_batches'))
                         + '; строк ' + str(progress['completed_rows']) + '/' + str(progress.get('expected_rows')) + '.')
        if direct.get('last_observed', {}).get('deadline_utc'):
            lines.append('Исходный предел времени (UTC): ' + direct['last_observed']['deadline_utc'])
    splice = result.get('spliceai', {})
    if splice:
        p = splice.get('progress', {})
        show_progress = splice.get('status') != 'COMPLETE_RECEIPT_OBSERVED' and 'variants_completed' in p
        lines.append('SpliceAI: ' + str(splice.get('status')) + ('; вариантов ' + str(p['variants_completed']) + '/' + str(p.get('variants_estimated')) if show_progress else '') + '.')
    for key, label in [('direct_core_context', 'Промежуточный VCF с33 полями'),
                       ('fully_annotated_vcf', 'Итоговый VCF с фактическим SpliceAI'),
                       ('cohort_and_sample_views', 'Таблицы и четыре набора образцов из этого VCF')]:
        item = result.get('deliverables', {}).get(key, {})
        value = item.get('status', 'NOT_CREATED_OR_WAITING')
        names = {'COMPLETE_RECEIPT_OBSERVED': 'готовы по квитанции',
                 'NOT_CREATED_OR_WAITING': 'ещё не созданы; завершение не подтверждено',
                 'WAITING_FOR_COMPLETED_INPUTS': 'контроллер ожидает завершённые входы',
                 'CONTROLLER_RUNNING_RESULT_NOT_COMPLETE': 'контроллер работает; результат ещё не завершён',
                 'CONTROLLER_REPORTED_FAILURE_OR_STOP': 'контроллер сообщил об остановке или ошибке',
                 'HISTORICAL_RESULT_NOT_VCF_FIRST': 'наблюдается прежняя выдача; новая выдача из VCF ещё не подтверждена',
                 'NO_MATCHING_PASS_RECEIPT': 'квитанция есть, успешное завершение не подтверждено',
                 'UNREADABLE_COMPLETION_RECEIPT': 'квитанция временно недоступна'}
        lines.append(label + ': ' + names.get(value, value) + '.')
    lines.append('Монитор только читает квитанции и PID; файлы результатов заново не проверяет.')
    return '\n'.join(lines)


def run(args, collect=snapshot, monotonic=time.monotonic, sleep=time.sleep):
    if not math.isfinite(args.max_hours) or not 0 < args.max_hours <= 8:
        raise ValueError('--max-hours must be >0 and <=8')
    if not math.isfinite(args.interval) or not 1 <= args.interval <= 3600:
        raise ValueError('--interval must be between 1 and 3600 seconds')
    end = monotonic() + args.max_hours * 3600
    while True:
        result = collect(args.root, args.local)
        print(human(result) if getattr(args, 'human', False) else json.dumps(result, ensure_ascii=False), flush=True)
        if args.report_dir:
            report(result, args.report_dir)
        remaining = end - monotonic()
        budget_remaining = result.get('databases', {}).get('seconds_until_original_deadline')
        if isinstance(budget_remaining, (int, float)):
            remaining = min(remaining, max(0, budget_remaining))
        if not args.watch or remaining <= 0 or result.get('databases', {}).get('budget_deadline_reached'):
            return
        # Short sleeps allow Ctrl-C; observing does not extend VM budget or stop billing.
        for _ in range(math.ceil(min(args.interval, remaining))):
            sleep(min(1, max(0, end - monotonic()), max(0, remaining)))
            remaining -= 1
            if remaining <= 0:
                return
            if monotonic() >= end:
                return


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--once', action='store_true', help='default: one read-only snapshot')
    mode.add_argument('--watch', action='store_true', help='bounded observation; no restart/download/kill')
    parser.add_argument('--root', type=Path, default=Path('/srv/marina'))
    parser.add_argument('--local', action='store_true', help='only allowlisted local transfer status; no secrets or proc')
    parser.add_argument('--interval', type=float, default=60)
    parser.add_argument('--max-hours', type=float, default=8)
    parser.add_argument('--human', action='store_true', help='short Russian status; default stdout stays JSON')
    parser.add_argument('--report-dir', type=Path, help='optional dedicated directory named monitor-reports')
    try:
        run(parser.parse_args())
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
