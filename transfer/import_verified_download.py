#!/usr/bin/env python3
"""Import a completed native download, with explicit SHA256 or size-only checks.

No cookies, tokens, browser state or HTTP authentication are accessed. The plan
must contain Files[{RelativePath,SizeBytes,Sha256}] from completed stage receipts
and SSH-verified remote inventory, never hashes calculated from an untrusted
download. --verification=size never hashes file bytes and uses a separate store.
By default this only checks the source; --execute publishes it.
Stop the owned backup controller before import; its partial files remain resumable.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import time
from write_local_index import atomic_json, publish_verified, publish_size_only, size_key, size_record_valid, write_index


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def path(value):
    if re.match(r'^[A-Za-z]:[\\/]', value):
        value = subprocess.check_output(['wslpath', '-u', value], text=True).strip()
    return Path(value).expanduser().resolve()


def sha256(file):
    result = hashlib.sha256()
    with file.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, help='Completed native download; never .crdownload/.part')
    parser.add_argument('--plan', required=True, help='Private trusted SSH-generated expected hashes/sizes')
    parser.add_argument('--relative-path', required=True, help='Exact RelativePath from the plan')
    parser.add_argument('--destination', required=True, help='Private backup root outside Git')
    parser.add_argument('--allow-onedrive', action='store_true', help='Explicit user-selected OneDrive destination')
    parser.add_argument('--verification', choices=('sha256', 'size'), default='sha256')
    parser.add_argument('--completed-download', action='store_true', help='Native completion or final-filename/stability evidence obtained')
    parser.add_argument('--completion-evidence', choices=('native_receipt', 'filesystem_stable', 'ssh_completed_response'), default='native_receipt')
    parser.add_argument('--stable-seconds', type=int, default=60)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    source, plan_path, destination = path(args.source), path(args.plan), path(args.destination)
    require(source.is_file() and source.suffix not in ('.crdownload', '.part', '.tmp'), 'Use a completed download')
    require(args.allow_onedrive or 'onedrive' not in str(destination).lower(), 'Explicit --allow-onedrive needed for this user-selected destination')
    require(not any((p / '.git').exists() for p in (destination, *destination.parents)), 'Destination cannot be in Git')
    plan = json.loads(plan_path.read_text(encoding='utf-8-sig'))
    selected = [entry for entry in plan['Files'] if entry['RelativePath'] == args.relative_path]
    require(len(selected) == 1, 'Exact unique trusted plan entry required')
    entry = selected[0]
    rel = PurePosixPath(entry['RelativePath'])
    require(rel.parts and not rel.is_absolute() and '..' not in rel.parts and '\\' not in entry['RelativePath'], 'Unsafe plan path')
    require(rel.parts[0] in ('outputs', 'metadata'), 'Only completed results and private code/config metadata permitted')
    require(not any(x in entry['RelativePath'].lower() for x in ('inputs/raw/', '/cache/', '/references/')), 'Public/raw resources excluded')
    require(re.fullmatch('[0-9a-f]{64}', entry['Sha256']) and isinstance(entry['SizeBytes'], int), 'Invalid expected identity')
    require(source.stat().st_size == entry['SizeBytes'], 'Download size differs from trusted expected size')
    if args.verification == 'size':
        require(args.completed_download, 'Completed-download evidence required for size-only import')
        require(not Path(str(source) + '.crdownload').exists() and not Path(str(source) + '.part').exists(),
                'Named partial download still exists')
        require(args.stable_seconds >= 0 and time.time() - source.stat().st_mtime >= args.stable_seconds,
                'Completed filename must have stable mtime before size-only import')
    fd = None
    if args.execute:
        import fcntl
        os.umask(0o077)
        destination.mkdir(parents=True, exist_ok=True, mode=0o700)
        fd = os.open(destination / '.backup.lock', os.O_RDWR | os.O_CREAT, 0o600)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(fd)
            raise RuntimeError('Stop the owned backup controller before import; partials remain resumable')
    try:
        start, before = time.monotonic(), source.stat()
        if args.verification == 'size':
            # User-selected mode: never read file bytes for a checksum here.
            if not args.execute:
                print(f'SIZE ONLY: expected and actual {entry["SizeBytes"]} bytes; SHA256 not checked; dry run', flush=True)
                return
            key = size_key(entry)
            folder = destination / 'size-only' / key
            folder.mkdir(parents=True, exist_ok=True, mode=0o700)
            stored = folder / rel.name
            if stored.exists():
                require(not stored.is_symlink() and stored.stat().st_size == entry['SizeBytes'],
                        'Differing size-only destination; refusing overwrite')
                require(size_record_valid(destination, entry) or os.path.samefile(stored, source),
                        'Changed existing size-only file; refusing overwrite')
            else:
                require(shutil.disk_usage(destination).free >= 15 * 2**30, '15 GiB reserve required')
                if source.stat().st_dev == folder.stat().st_dev:
                    os.link(source, stored)
                else:
                    require(shutil.disk_usage(destination).free >= entry['SizeBytes'] + 15 * 2**30,
                            'Insufficient space for cross-volume size-only copy')
                    part = folder / (rel.name + '.size-copy.part')
                    shutil.copyfile(source, part)
                    require(part.stat().st_size == entry['SizeBytes'], 'Copied size differs')
                    os.link(part, stored)
                    part.unlink()
            after = source.stat()
            require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'Source changed during import')
            publish_size_only(destination, entry, stored, args.completion_evidence)
            count = write_index(destination)
            verified_map_path = destination / 'verified-files.json'
            verified_map = json.loads(verified_map_path.read_text(encoding='utf-8-sig')) if verified_map_path.exists() else {}
            size_map = json.loads((destination / 'size-only-files.json').read_text(encoding='utf-8-sig'))
            sha_bams = {v['Sha256'] for v in verified_map.values() if v['LocalPath'].endswith('.bam')}
            size_bams = {v['SourceKey'] for v in size_map.values() if v['LocalPath'].endswith('.bam')}
            status = {'Status': 'PARTIAL_BACKUP', 'UpdatedUtc': dt.datetime.now(dt.timezone.utc).isoformat(),
                      'VerificationModeForNewFiles': 'size_only', 'CompletedFileLinks': count,
                      'SHA256VerifiedFileLinks': len(verified_map), 'SizeOnlyFileLinks': len(size_map),
                      'SHA256VerifiedBAMs': len(sha_bams), 'SizeOnlyBAMs': len(size_bams),
                      'CompletedBAMs': len(sha_bams) + len(size_bams), 'CompletedSnapshot': False,
                      'FullSnapshotSHA256Pass': False, 'Note': 'New files checked by completion and size only; prior SHA256 history preserved'}
            atomic_json(destination / 'backup-status.json', status)
            atomic_json(destination / 'native-import-status.json', {
                'Status': 'SIZE_ONLY', 'SourceKey': key, 'ExpectedSizeBytes': entry['SizeBytes'],
                'ActualSizeBytes': stored.stat().st_size, 'ExpectedServerSha256': entry['Sha256'],
                'ActualSha256': None, 'CompletionEvidence': args.completion_evidence})
            print(f'SIZE ONLY imported: {entry["SizeBytes"]} bytes; SHA256 not checked; INDEX {count} files', flush=True)
            return
        actual = sha256(source)
        after = source.stat()
        require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'Download changed during checksum')
        require(actual == entry['Sha256'], 'Downloaded SHA256 differs; no promotion performed')
        print(f'PASS source SHA256; {entry["SizeBytes"]} bytes; {time.monotonic()-start:.1f} seconds', flush=True)
        if not args.execute:
            print('Dry run: no files, links or backup state changed.', flush=True)
            return
        objects, staging = destination / 'objects', destination / '.staging'
        objects.mkdir(exist_ok=True, mode=0o700)
        staging.mkdir(exist_ok=True, mode=0o700)
        obj = objects / entry['Sha256']
        if obj.exists():
            require(not obj.is_symlink() and obj.stat().st_size == entry['SizeBytes'] and sha256(obj) == entry['Sha256'],
                    'Differing final object; refusing overwrite')
        else:
            require(shutil.disk_usage(destination).free >= 15 * 2**30, '15 GiB reserve required')
            if source.stat().st_dev == objects.stat().st_dev:
                os.link(source, obj)  # Same-volume native download: no second full BAM copy.
            else:
                require(shutil.disk_usage(destination).free >= entry['SizeBytes'] + 15 * 2**30,
                        'Insufficient space for cross-volume copy plus reserve')
                part = staging / (entry['Sha256'] + '.native-copy.part')
                shutil.copyfile(source, part)
                require(sha256(part) == entry['Sha256'], 'Cross-volume copy checksum failed')
                os.link(part, obj)
                part.unlink()
        state_path = destination / 'verified-objects.json'
        state = json.loads(state_path.read_text()) if state_path.exists() else {}
        state[entry['Sha256']] = {'size': obj.stat().st_size, 'mtime_ns': obj.stat().st_mtime_ns}
        atomic_json(state_path, state)
        map_path = destination / 'verified-files.json'
        mapping = json.loads(map_path.read_text(encoding='utf-8')) if map_path.exists() else {}
        publish_verified(destination, entry, mapping)
        count = write_index(destination)
        atomic_json(destination / 'native-import-status.json', {'Status': 'PASS', 'VerifiedBytes': entry['SizeBytes'],
                    'VerifiedFilesInIndex': count, 'Sha256': entry['Sha256'],
                    'Note': 'One file imported; full snapshot PASS requires all planned objects'})
        atomic_json(destination / 'backup-status.json', {
            'Status': 'PARTIAL_BACKUP', 'UpdatedUtc': dt.datetime.now(dt.timezone.utc).isoformat(),
            'VerifiedObjects': len(state), 'VerifiedFileLinks': count,
            'VerifiedBAMs': len({v['Sha256'] for v in mapping.values() if v['LocalPath'].endswith('.bam')}),
            'CompletedSnapshot': False,
            'Note': 'Native download imported after SHA256; full snapshot requires every planned object'})
        print(f'PASS imported verified object; INDEX contains {count} files. Source download preserved.', flush=True)
    finally:
        if fd is not None:
            fcntl.flock(fd, fcntl.LOCK_UN)
            os.close(fd)


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, ValueError, OSError) as exc:
        print('IMPORT FAILED: ' + str(exc), file=sys.stderr)
        sys.exit(1)
