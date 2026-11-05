#!/usr/bin/env python3
"""Resume verified data uploads from Windows files through WSL OpenSSH and rsync."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import shutil
import subprocess
import time


def windows_source(value):
    if re.match(r'^[A-Za-z]:[\\/]', value):
        value = subprocess.check_output(['wslpath', '-u', value], text=True).strip()
    return Path(value).resolve(strict=True)


def sha256(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--host', required=True, help='Confirmed SSH configuration alias')
    parser.add_argument('--ssh-config', required=True)
    parser.add_argument('--known-hosts', required=True)
    parser.add_argument('--remote-root', required=True)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--rehash-local', action='store_true')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.host):
        raise ValueError('Use a confirmed SSH alias; shell expressions are not accepted')
    root = PurePosixPath(args.remote_root)
    if not root.is_absolute() or '..' in root.parts or len(root.parts) < 3:
        raise ValueError('Use a dedicated absolute project directory')
    if any(ord(char) < 32 for char in args.remote_root):
        raise ValueError('Control characters in remote root')
    root = str(root)
    manifest_path = Path(args.manifest).resolve(strict=True)
    source_manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    if source_manifest.get('SchemaVersion') != 1 or not source_manifest.get('Files'):
        raise ValueError('Unsupported or empty manifest')
    sources, entries, seen = [], [], set()
    for entry in source_manifest['Files']:
        source = windows_source(entry['SourcePath'])
        rel = PurePosixPath(entry['RelativePath'])
        if rel.is_absolute() or '..' in rel.parts or '\\' in entry['RelativePath']:
            raise ValueError('Unsafe remote relative path')
        if entry['RelativePath'] in seen:
            raise ValueError('Duplicate remote relative path')
        seen.add(entry['RelativePath'])
        if not re.fullmatch(r'[0-9a-f]{64}', entry['Sha256']):
            raise ValueError('Invalid SHA256')
        if not source.is_file() or source.stat().st_size != entry['SizeBytes']:
            raise ValueError('A source changed size or is not a regular file')
        if args.rehash_local and sha256(source) != entry['Sha256']:
            raise ValueError('Source checksum changed; review and rebuild the manifest')
        sources.append(source)
        entries.append({key: entry[key] for key in ('RelativePath', 'SizeBytes', 'Sha256')})
    total = sum(entry['SizeBytes'] for entry in entries)
    print(f'Plan: {len(entries)} files, {total / 2**30:.3f} GiB, verified staging under {args.host}:{root}', flush=True)
    if not args.execute:
        print('Dry run: no SSH connection or transfer was made.', flush=True)
        return
    for tool in ('ssh', 'rsync'):
        if not shutil.which(tool):
            raise RuntimeError(f'{tool} is required inside WSL')
    ssh_config = Path(args.ssh_config).expanduser().resolve(strict=True)
    known_hosts = Path(args.known_hosts).expanduser().resolve(strict=True)
    ssh = ['ssh', '-F', str(ssh_config), '-o', 'BatchMode=yes', '-o', 'StrictHostKeyChecking=yes',
           '-o', f'UserKnownHostsFile={known_hosts}', '-o', 'ConnectTimeout=30', '-o', 'ServerAliveInterval=30',
           '-o', 'ServerAliveCountMax=4']

    def remote(command):
        return subprocess.check_output(ssh + [args.host, command], text=True).strip()

    def upload(source, destination, append=False):
        command = ['rsync', '--partial', '--protect-args', '--chmod=F600', '--info=progress2', '--human-readable']
        if append:
            command.append('--append-verify')
        command += ['-e', shlex.join(ssh), '--', str(source), f'{args.host}:{destination}']
        subprocess.run(command, check=True)

    state = manifest_path.parent / 'transfer-state-wsl'
    state.mkdir(parents=True, mode=0o700, exist_ok=True)
    payload = json.dumps({'SchemaVersion': 1, 'Files': entries}, ensure_ascii=False, separators=(',', ':')).encode('utf-8')
    plan_id = hashlib.sha256(payload).hexdigest()
    local_payload = state / f'{plan_id}.json'
    local_payload.write_bytes(payload)
    helper = Path(__file__).with_name('verify_upload.py').resolve(strict=True)
    helper_id = sha256(helper)
    helper_remote = f'{root}/.transfer/verify-{helper_id}.py'
    manifest_remote = f'{root}/.transfer/{plan_id}.json'
    staging_root = f'{root}/.upload-staging/{plan_id}'
    remote('umask 077; command -v python3 >/dev/null && mkdir -p -- ' + shlex.quote(f'{root}/.transfer') + ' ' + shlex.quote(staging_root))
    upload(helper, helper_remote)
    upload(local_payload, manifest_remote)
    helper_command = shlex.join(['python3', helper_remote, '--root', root, '--manifest', manifest_remote, '--plan-id', plan_id])
    plan = json.loads(remote(helper_command + ' --mode plan'))
    pending = plan['Pending']
    remaining = sum(entries[index]['SizeBytes'] for index in pending)
    free = int(remote('df -B1 --output=avail -- ' + shlex.quote(root) + ' | tail -n 1'))
    if free < remaining + 10 * 2**30:
        raise RuntimeError('Insufficient disk for pending uploads plus 10 GiB reserve')
    started = time.monotonic()
    sent = 0
    print(f'Already verified: {plan["Complete"]}; pending uploads: {len(pending)}; staged ready: {len(plan["Ready"])}', flush=True)
    for sequence, index in enumerate(pending, 1):
        print(f'Uploading file index {index + 1}/{len(entries)}; pending sequence {sequence}/{len(pending)}', flush=True)
        file_started = time.monotonic()
        upload(sources[index], f'{staging_root}/{index:04d}.part', append=index in plan['Resume'])
        elapsed = max(time.monotonic() - file_started, 0.001)
        sent += entries[index]['SizeBytes']
        print(f'Uploaded file index {index + 1}; source size {entries[index]["SizeBytes"] / 2**20:.1f} MiB; elapsed {elapsed:.1f}s; effective rate {entries[index]["SizeBytes"] / 2**20 / elapsed:.2f} MiB/s', flush=True)
    print('Verifying SHA256 and promoting completed files; originals remain unchanged.', flush=True)
    summary = json.loads(remote(helper_command + ' --mode finalize'))
    if summary['Failed'] or summary['Complete'] != len(entries):
        raise RuntimeError('Not all remote files passed SHA256 verification')
    report = {'SchemaVersion': 1, 'VerifiedFiles': summary['Complete'], 'TotalBytes': total,
              'ElapsedSeconds': round(time.monotonic() - started, 2), 'UploadedSourceBytes': sent,
              'HostAlias': args.host, 'RemoteRoot': root, 'PlanId': plan_id}
    (state / 'verified-transfer-summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
