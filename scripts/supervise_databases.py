#!/usr/bin/env python3
"""Bounded server-only database installation/resume within the existing VM allowance.

This stops only a verified owned installer group. It does not stop VM billing.
"""
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time
import uuid

STOP = False


def need(ok, message):
    if not ok:
        raise ValueError(message)


def utc(value):
    result = dt.datetime.fromisoformat(value.replace('Z', '+00:00'))
    need(result.tzinfo is not None, 'VM creation time requires an explicit timezone.')
    return result.astimezone(dt.timezone.utc)


def deadline(created, rate, cap):
    need(0 < rate and 0 < cap <= 40, 'The supervisor cannot increase the planned $40 VM allowance.')
    return created + dt.timedelta(hours=cap / rate)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name('.' + path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with temporary.open('x') as stream:
            json.dump(value, stream, indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def process_token(pid, installer):
    """Pin UID, start time, process-group identity and exact allowed source path."""
    proc = Path('/proc') / str(pid)
    if not proc.exists():
        return None
    stat = (proc / 'stat').read_text().rsplit(') ', 1)[1].split()
    if stat[0] == 'Z':
        return None
    command = (proc / 'cmdline').read_bytes().split(b'\0')
    expected = {os.fsencode(installer.resolve()),
                os.fsencode(installer.resolve().parents[1] / 'scripts/prepare_annotation_resources.py')}
    need(proc.stat().st_uid == os.getuid() and os.getpgid(pid) == pid and any(value in expected for value in command),
         'Refusing to adopt/stop a process outside the exact owned database installer group.')
    return dict(pid=pid, start_ticks=stat[19], uid=os.getuid())


def still_owned(token, installer):
    if not token:
        return False
    try:
        current = process_token(token['pid'], installer)
    except (FileNotFoundError, ProcessLookupError, ValueError):
        return False
    return current == token


def stop_owned(token, installer):
    if still_owned(token, installer):
        os.killpg(token['pid'], signal.SIGTERM)
        stop = time.monotonic() + 15
        while still_owned(token, installer) and time.monotonic() < stop:
            time.sleep(.2)
        if still_owned(token, installer):
            os.killpg(token['pid'], signal.SIGKILL)


def transient_failure(text):
    text = text.lower()
    fatal = ('checksum mismatch', 'sha256 mismatch', 'byte count mismatch', 'object identity changed',
             'unsafe ', 'refusing', 'snapshot header/date changed', 'metadata version unexpected',
             'version does not match', 'insufficient', 'no-replace', 'httperror: http error 403',
             'httperror: http error 404', 'httperror: http error 412', 'httperror: http error 416')
    if any(word in text for word in fatal):
        return False
    return any(word in text for word in ('download time limit reached', 'timed out', 'timeouterror',
        'connection reset', 'connection aborted', 'connection refused', 'connectionrefusederror',
        'remote end closed', 'temporary failure',
        'http error 429', 'http error 500', 'http error 502', 'http error 503', 'http error 504'))


def retry_backoff_seconds(previous_failures, progress_made):
    """Provider refusals without progress back off; the absolute VM deadline still wins."""
    failures = 0 if progress_made else previous_failures + 1
    return failures, min(900, 60 * 2 ** min(4, max(0, failures - 1)))


def verified_ready(manifest_path, recipe_path):
    if not manifest_path.is_file():
        return False
    manifest = json.loads(manifest_path.read_text())
    need(manifest.get('status') == 'PASS' and manifest.get('schema_version') == 1
         and manifest.get('assembly') == 'GRCh38' and manifest.get('config_sha256') == digest(recipe_path),
         'Completed annotation resource manifest differs from the pinned recipe.')
    need(manifest.get('artifacts'), 'A READY resource manifest cannot have an empty artifact list.')
    root = manifest_path.parent.resolve()
    for item in manifest['artifacts']:
        relative = Path(item['path'])
        need(not relative.is_absolute() and '..' not in relative.parts, 'Unsafe completed resource path.')
        path = root / relative
        need(not path.is_symlink() and path.resolve(strict=True).is_relative_to(root)
             and path.stat().st_size == item['bytes'] and digest(path) == item['sha256'],
             'Completed resource artifact checksum mismatch; not READY.')
    return True


def run(args):
    global STOP
    os.umask(0o077)
    root = args.root.resolve(strict=True)
    installer, recipe = args.installer.resolve(strict=True), args.config.resolve(strict=True)
    cfg = json.loads(recipe.read_text())
    need(re.fullmatch(r'[A-Za-z0-9_.-]+', cfg['bundle_id']), 'Invalid resource bundle ID.')
    ready = root / 'resources/annotation' / cfg['bundle_id'] / 'resource-manifest.json'
    settings = root / 'private/server-budget.json'
    budget = json.loads(settings.read_text()) if settings.is_file() else {}
    created_text = args.created_utc or budget.get('vm_created_utc')
    need(created_text, 'Supply the actual VM creation timestamp; it is never reset to the current time.')
    created = utc(created_text)
    now = dt.datetime.now(dt.timezone.utc)
    need(created <= now + dt.timedelta(seconds=5), 'VM creation timestamp cannot be in the future.')
    end = deadline(created, args.hourly_usd, args.planned_usd)
    monotonic_end = time.monotonic() + max(0, (end - now).total_seconds())
    need(1 <= args.max_attempts <= 2 and 1 <= args.workers <= 32,
         'Fixed-attempt mode allows at most two attempts; at most 32 connections are allowed.')
    state_path = root / 'private/databases-supervisor-status.json'
    lockpath = root / 'private/.databases-supervisor.lock'
    identity = dict(created_utc=created.isoformat(), hourly_usd=args.hourly_usd, planned_vm_usd=args.planned_usd,
        installer=str(installer), installer_sha256=digest(installer), recipe_sha256=digest(recipe))
    logpath = root / 'logs/annotation-databases.log'
    logpath.parent.mkdir(parents=True, exist_ok=True)
    token, child = None, None
    with lockpath.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = json.loads(state_path.read_text()) if state_path.exists() else dict(attempts_started=0)
        need(not state_path.is_symlink() and (not state.get('identity') or state['identity'] == identity),
             'Existing supervisor state belongs to a different VM allowance or recipe; refusing reset.')
        def remaining():
            return min((end - dt.datetime.now(dt.timezone.utc)).total_seconds(), monotonic_end - time.monotonic())
        def save(status, **extra):
            state.update(schema_version=1, status=status, identity=identity, deadline_utc=end.isoformat(),
                updated_utc=dt.datetime.now(dt.timezone.utc).isoformat(), current_pid=token['pid'] if token else None,
                maximum_attempts=None if args.retry_until_deadline else args.max_attempts,
                retry_policy='transient_only_until_original_absolute_deadline' if args.retry_until_deadline else 'fixed_attempt_limit',
                ready_manifest=str(ready), **extra)
            atomic_json(state_path, state)
            print(json.dumps({'status': status, 'attempts': state['attempts_started'],
                'remaining_budget_hours': round(max(0, remaining()) / 3600, 3), 'pid': state['current_pid']}), flush=True)
        if args.adopt_pid:
            token = process_token(args.adopt_pid, installer)
            need(token is not None, 'The requested installer is not running; no adoption or budget reset performed.')
            if token:
                state['attempts_started'] = max(1, state['attempts_started'])
                state['attempt_log_offset'] = logpath.stat().st_size if logpath.exists() else 0
                save('SUPERVISING_EXISTING_INSTALLER')
        elif state.get('current_pid'):
            token = process_token(state['current_pid'], installer)
            if token:
                save('SUPERVISING_EXISTING_INSTALLER')
        try:
            while True:
                if STOP:
                    stop_owned(token, installer)
                    save('CANCELLED_OWN_INSTALLER_STOPPED')
                    return 130
                if remaining() <= 0:
                    stop_owned(token, installer)
                    save('VM_ALLOWANCE_DEADLINE_OWN_INSTALLER_STOPPED', billing_continues_until_VM_termination=True)
                    return 3
                if not token:
                    if verified_ready(ready, recipe):
                        save('READY_VALIDATED', resource_manifest_sha256=digest(ready))
                        return 0
                    if not args.retry_until_deadline and state['attempts_started'] >= args.max_attempts:
                        save('MAX_ATTEMPTS_REACHED_PARTIALS_PRESERVED')
                        return 4
                    # A manual supervisor restart must not bypass an integrity/access failure.
                    if state['attempts_started']:
                        prior = ''
                        if logpath.is_file():
                            with logpath.open('rb') as log:
                                log.seek(max(state.get('attempt_log_offset', 0), logpath.stat().st_size - 200000))
                                prior = log.read().decode(errors='replace')
                        need(transient_failure(prior), 'Refusing restart after a nontransient installer failure.')
                    env = dict(os.environ, MARINA_DOWNLOAD_WORKERS=str(args.workers), MARINA_DATA_ROOT=str(root))
                    with logpath.open('a') as log:
                        log.write('\nBOUNDED_SUPERVISOR_ATTEMPT ' + str(state['attempts_started'] + 1) + '\n')
                        log.flush()
                        state['attempt_log_offset'] = log.tell()
                        child = subprocess.Popen(['bash', str(installer)], stdin=subprocess.DEVNULL,
                            stdout=log, stderr=subprocess.STDOUT, start_new_session=True, env=env)
                    state['attempts_started'] += 1
                    for _ in range(40):
                        token = process_token(child.pid, installer)
                        if token or child.poll() is not None:
                            break
                        time.sleep(.05)
                    save('INSTALLER_RUNNING')
                # For a child started here, Popen's actual exit/reap is authoritative.
                # /proc cmdline can be momentarily empty while bash execs Python;
                # that transition must not be mistaken for an installer failure.
                while ((child is not None and child.poll() is None) or
                       (token and still_owned(token, installer))) and not STOP and remaining() > 0:
                    time.sleep(min(5, max(0, remaining())))
                if STOP or remaining() <= 0:
                    continue
                token = None
                if child:
                    child.poll()
                if verified_ready(ready, recipe):
                    save('READY_VALIDATED', resource_manifest_sha256=digest(ready))
                    return 0
                text = ''
                if logpath.is_file():
                    with logpath.open('rb') as log:
                        log.seek(max(state.get('attempt_log_offset', 0), logpath.stat().st_size - 200000))
                        text = log.read().decode(errors='replace')
                if not transient_failure(text):
                    save('NONRETRYABLE_FAILURE_PARTIALS_PRESERVED')
                    return 5
                failures, delay = retry_backoff_seconds(state.get('consecutive_no_progress_failures', 0), 'CACHE_PROGRESS ' in text)
                save('TRANSIENT_DOWNLOAD_FAILURE_RETRY_ELIGIBLE', consecutive_no_progress_failures=failures,
                     next_retry_delay_seconds=delay)
                until = time.monotonic() + min(delay, max(0, remaining()))
                while time.monotonic() < until and not STOP and remaining() > 0:
                    time.sleep(min(5, until - time.monotonic()))
        except BaseException:
            stop_owned(token, installer)
            save('SUPERVISOR_FAILURE_OWN_INSTALLER_STOPPED')
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    code = Path(__file__).resolve().parents[1]
    parser.add_argument('--root', type=Path, default=Path('/srv/marina'))
    parser.add_argument('--installer', type=Path, default=code / 'publicstages/01_prepare_databases.sh')
    parser.add_argument('--config', type=Path, default=code / 'config/annotation-resources.example.json')
    parser.add_argument('--created-utc', help='Actual VM creation time; alternatively read private/server-budget.json.')
    parser.add_argument('--hourly-usd', type=float, default=2.39)
    parser.add_argument('--planned-usd', type=float, default=40)
    parser.add_argument('--max-attempts', type=int, default=2)
    parser.add_argument('--retry-until-deadline', action='store_true',
        help='Retry only transient network failures with capped exponential backoff until the original absolute VM allowance deadline; do not reset attempts or time.')
    parser.add_argument('--workers', type=int, default=32)
    parser.add_argument('--adopt-pid', type=int, help='Adopt only an exact owned, separate installer process group.')
    args = parser.parse_args()
    def cancel(signum, frame):
        global STOP
        STOP = True
    signal.signal(signal.SIGTERM, cancel)
    signal.signal(signal.SIGINT, cancel)
    raise SystemExit(run(args))


if __name__ == '__main__':
    main()
