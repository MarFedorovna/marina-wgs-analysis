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
    alpha = load(private / 'alphamissense-control-v1/public-controls-receipt.json')
    result['annotation']['separate_alphamissense_controls'] = fields(alpha, 'status', 'controls', 'tested')
    result['xy'] = dict(status='PLAN_ONLY_NO_EXECUTOR', clinical_sex_assignment=False,
                        somatic_analysis=False, ploidy='unknown_all_four_until_review',
                        par_alignment_contract='unverified_until_separate_review')
    result['limitations'] = ['No old MT v1 failure is used as current v2 status.',
        'Receipts are reported, not a fresh technical/clinical audit.',
        'No completion claim for annotation, X/Y, full genome or interpretation.']
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


def run(args, collect=snapshot, monotonic=time.monotonic, sleep=time.sleep):
    if not math.isfinite(args.max_hours) or not 0 < args.max_hours <= 8:
        raise ValueError('--max-hours must be >0 and <=8')
    if not math.isfinite(args.interval) or not 1 <= args.interval <= 3600:
        raise ValueError('--interval must be between 1 and 3600 seconds')
    end = monotonic() + args.max_hours * 3600
    while True:
        result = collect(args.root, args.local)
        print(json.dumps(result, ensure_ascii=False), flush=True)
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
    parser.add_argument('--report-dir', type=Path, help='optional dedicated directory named monitor-reports')
    try:
        run(parser.parse_args())
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
