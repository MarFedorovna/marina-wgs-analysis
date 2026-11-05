#!/usr/bin/env python3
"""Bounded local watcher for completed, explicitly planned browser downloads.

It accesses no browser authentication and starts no cloud transfers. Expected
filenames, sizes and SHA256 come from a private trusted native-download plan.
Only import_verified_download.py publishes files after full SHA256 verification.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from write_local_index import atomic_json, size_record_valid


def local_path(value):
    if re.match(r'^[A-Za-z]:[\\/]', value):
        value = subprocess.check_output(['wslpath', '-u', value], text=True).strip()
    return Path(value).expanduser().resolve()


def verified(destination, entry):
    state_file = destination / 'verified-objects.json'
    if not state_file.is_file():
        return False
    state = json.loads(state_file.read_text(encoding='utf-8-sig'))
    obj = destination / 'objects' / entry['Sha256']
    return obj.is_file() and not obj.is_symlink() and obj.stat().st_size == entry['SizeBytes'] and state.get(entry['Sha256']) == {
        'size': obj.stat().st_size, 'mtime_ns': obj.stat().st_mtime_ns}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True)
    parser.add_argument('--download-directory', required=True)
    parser.add_argument('--destination', required=True)
    parser.add_argument('--status-file', help='Optional ignored aggregate status mirror')
    parser.add_argument('--allow-onedrive', action='store_true')
    parser.add_argument('--verification', choices=('sha256', 'size'), default='sha256')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--once', action='store_true')
    parser.add_argument('--interval', type=int, default=60)
    parser.add_argument('--max-hours', type=float, default=8)
    parser.add_argument('--deadline-utc', help='Absolute local-download stop time, ISO8601 UTC')
    parser.add_argument('--resume', action='store_true', help='Explicitly resume after a saved user-requested pause')
    args = parser.parse_args()
    if args.interval < 60 or not 0 < args.max_hours <= 8:
        raise ValueError('Interval >=60 seconds and max-hours <=8 required')
    plan, downloads, destination = map(local_path, (args.plan, args.download_directory, args.destination))
    if not args.execute:
        entries = json.loads(plan.read_text(encoding='utf-8-sig'))['Files']
        print(f'DRY RUN: {len(entries)} planned files; no authentication, transfers or imports started')
        return
    import fcntl
    os.umask(0o077)
    if any((p / '.git').exists() for p in (destination, *destination.parents)):
        raise ValueError('Backup destination cannot be inside Git')
    if 'onedrive' in str(destination).lower() and not args.allow_onedrive:
        raise ValueError('Use --allow-onedrive for an explicitly user-selected destination')
    destination.mkdir(parents=True, exist_ok=True)
    pause_flag = destination / 'Settings' / 'LOCAL_DOWNLOADS_PAUSED.json'
    if pause_flag.exists():
        if not args.resume:
            print('PAUSED by user request; use the saved explicit resume command', flush=True)
            return
        pause_flag.unlink()
    lock = (destination / '.native-watch.lock').open('a')
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    deadline, failed, observed = time.monotonic() + args.max_hours * 3600, {}, {}
    wall_deadline = time.time() + args.max_hours * 3600
    if args.deadline_utc:
        absolute = dt.datetime.fromisoformat(args.deadline_utc.replace('Z', '+00:00'))
        if absolute.tzinfo is None:
            raise ValueError('Deadline must contain a UTC offset')
        deadline = min(deadline, time.monotonic() + max(0, absolute.timestamp() - time.time()))
        wall_deadline = min(wall_deadline, absolute.timestamp())
    try:
        while time.monotonic() < deadline and time.time() < wall_deadline:
            if pause_flag.exists():
                break
            entries = json.loads(plan.read_text(encoding='utf-8-sig'))['Files']
            imported, errors = 0, 0
            for entry in entries:
                if time.monotonic() >= deadline or time.time() >= wall_deadline:
                    break
                if verified(destination, entry) or (args.verification == 'size' and size_record_valid(destination, entry)):
                    continue
                name = entry['DownloadFileName']
                if Path(name).name != name or any(c in name for c in ('/', '\\', ':')) or any(ord(c) < 32 for c in name):
                    raise ValueError('Plan download name must be a basename')
                source = downloads / name
                if not source.is_file() or source.stat().st_size != entry['SizeBytes']:
                    continue
                identity = (source.stat().st_size, source.stat().st_mtime_ns)
                if failed.get(entry['RelativePath']) == identity:
                    continue  # Bad checksum is visible; retry only if the source changes.
                if args.verification == 'size':
                    if Path(str(source) + '.crdownload').exists() or Path(str(source) + '.part').exists():
                        continue
                    previous = observed.get(entry['RelativePath'])
                    if previous is None or previous[0] != identity:
                        observed[entry['RelativePath']] = (identity, time.monotonic())
                        continue
                    if time.monotonic() - previous[1] < 60:
                        continue
                command = ['python3', str(Path(__file__).with_name('import_verified_download.py')),
                           '--source', str(source), '--plan', str(plan), '--relative-path', entry['RelativePath'],
                           '--destination', str(destination), '--verification', args.verification, '--execute']
                if args.verification == 'size':
                    command += ['--completed-download', '--completion-evidence', 'filesystem_stable']
                if args.allow_onedrive:
                    command.append('--allow-onedrive')
                result = subprocess.run(command, text=True, capture_output=True,
                                        timeout=max(1, deadline - time.monotonic()))
                print(result.stdout, end='', flush=True)
                if result.returncode:
                    print(result.stderr, end='', file=sys.stderr, flush=True)
                    errors += 1
                    if 'Stop the owned backup controller' not in result.stderr:
                        failed[entry['RelativePath']] = identity
                else:
                    imported += 1
                    if args.status_file and (destination / 'backup-status.json').is_file():
                        target = local_path(args.status_file)
                        target.parent.mkdir(parents=True, exist_ok=True)
                        atomic_json(target, json.loads((destination / 'backup-status.json').read_text(encoding='utf-8-sig')))
            sha_count = sum(verified(destination, entry) for entry in entries)
            size_count = sum(size_record_valid(destination, entry) and not verified(destination, entry) for entry in entries)
            count = sha_count + (size_count if args.verification == 'size' else 0)
            stopping = args.once or time.monotonic() + args.interval >= deadline
            atomic_json(destination / 'native-watch-status.json', {
                'Status': 'STOPPED' if stopping else 'WATCHING', 'CompletedPlannedFiles': count, 'PlannedFiles': len(entries),
                'VerificationMode': args.verification, 'SHA256VerifiedPlannedFiles': sha_count, 'SizeOnlyPlannedFiles': size_count,
                'ImportedThisPoll': imported, 'ErrorsThisPoll': errors, 'FailedUnchangedSources': len(failed),
                'CompletedSnapshot': False, 'Note': 'Native file imports only; no browser automation or full-snapshot completion claim'})
            print(f'Native import watcher: {count}/{len(entries)} completed; SHA256 {sha_count}, size-only {size_count}; errors {errors}', flush=True)
            if stopping:
                break
            time.sleep(min(args.interval, max(0, deadline - time.monotonic()), max(0, wall_deadline-time.time())))
    finally:
        if (time.monotonic() >= deadline or time.time() >= wall_deadline) and not args.once:
            pause_flag.parent.mkdir(parents=True, exist_ok=True)
            atomic_json(pause_flag, {'Status': 'PAUSED', 'Reason': 'User-requested local-download deadline',
                        'PausedUtc': dt.datetime.now(dt.timezone.utc).isoformat(),
                        'PartialsPreserved': True, 'ServerWorkUnaffected': True, 'RequiresUserResume': True})
            atomic_json(destination / 'native-watch-status.json', {'Status': 'PAUSED', 'VerificationMode': args.verification,
                        'Reason': 'User-requested local-download deadline', 'CompletedSnapshot': False,
                        'PartialsPreserved': True, 'ServerWorkUnaffected': True})
            print('PAUSED at user deadline; partials, plan and INDEX preserved; server work unaffected', flush=True)
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()


if __name__ == '__main__':
    main()
