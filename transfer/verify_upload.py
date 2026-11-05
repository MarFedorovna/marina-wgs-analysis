#!/usr/bin/env python3
"""Verify SFTP staging files and atomically promote complete files inside one root."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def contained(root, relative):
    rel = PurePosixPath(relative)
    if not relative or rel.is_absolute() or '..' in rel.parts or '\\' in relative:
        raise ValueError('Unsafe relative path')
    if any(ord(char) < 32 for char in relative):
        raise ValueError('Control characters in relative path')
    candidate = root.joinpath(*rel.parts)
    # Resolve existing parents as well, refusing symlinks out of the transfer root.
    candidate.resolve(strict=False).relative_to(root)
    if candidate.is_symlink():
        raise ValueError('Destination/staging symlinks are not allowed')
    return candidate


def matches(path, entry):
    return path.is_file() and path.stat().st_size == entry['SizeBytes'] and digest(path) == entry['Sha256']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--plan-id', required=True)
    parser.add_argument('--mode', choices=('plan', 'finalize'), required=True)
    args = parser.parse_args()
    if not re.fullmatch(r'[0-9a-f]{64}', args.plan_id):
        raise ValueError('Invalid plan identifier')
    root = Path(args.root).resolve(strict=True)
    manifest_path = Path(args.manifest).resolve(strict=True)
    manifest_path.relative_to(root)
    payload = manifest_path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != args.plan_id:
        raise ValueError('Transferred manifest checksum mismatch')
    manifest = json.loads(payload.decode('utf-8'))
    if manifest.get('SchemaVersion') != 1 or not isinstance(manifest.get('Files'), list):
        raise ValueError('Unsupported manifest')
    transfer_dir = contained(root, '.transfer')
    transfer_dir.mkdir(mode=0o700, exist_ok=True)
    lock_path = contained(root, '.transfer/verify.lock')
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        pending, resume, ready, failures, complete = [], [], [], [], 0
        seen = set()
        for index, entry in enumerate(manifest['Files']):
            relative = entry['RelativePath']
            if relative in seen:
                raise ValueError('Duplicate destination path')
            seen.add(relative)
            if not re.fullmatch(r'[0-9a-f]{64}', entry['Sha256']):
                raise ValueError('Invalid checksum')
            if not isinstance(entry['SizeBytes'], int) or entry['SizeBytes'] < 0:
                raise ValueError('Invalid file size')
            final = contained(root, relative)
            staging = contained(root, f'.upload-staging/{args.plan_id}/{index:04d}.part')
            if final.exists():
                if not matches(final, entry):
                    raise ValueError(f'Existing destination differs; refusing overwrite at file index {index}')
                complete += 1
                continue
            if args.mode == 'plan':
                if staging.exists():
                    if not staging.is_file():
                        raise ValueError(f'Staging path is not a regular file at index {index}')
                    if matches(staging, entry):
                        ready.append(index)
                        continue
                    if staging.stat().st_size < entry['SizeBytes']:
                        resume.append(index)
                pending.append(index)
                continue
            if not matches(staging, entry):
                failures.append(index)
                continue
            final.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
            # Every file is promoted only after size and SHA256 verification.
            os.chmod(staging, 0o600)
            staging.rename(final)
            complete += 1
        print(json.dumps({'Complete': complete, 'Pending': pending, 'Resume': resume, 'Ready': ready, 'Failed': failures}))
        return 2 if failures else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as error:
        print(f'Transfer verification failed: {error}', file=sys.stderr)
        sys.exit(1)
