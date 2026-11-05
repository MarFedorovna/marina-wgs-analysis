#!/usr/bin/env python3
"""Back up completed stages through WSL SSH/rsync, with SHA256 and immutable snapshots.

Default is a local dry-run. --execute snapshots small remote configuration/log/code
files and pulls only receipt-completed result artifacts. Inputs/raw, public caches,
references, credentials, symlinks and partial output directories are never selected.
Snapshots contain hard links to a content-addressed object store, avoiding duplicate
BAM storage. Keep the whole destination (objects + snapshots), not one tree alone.
"""
import argparse
import base64
import ctypes
import errno
import datetime as dt
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import signal
import shlex
import shutil
import subprocess
import sys
import time


# This helper runs only with --execute. It never modifies pipeline inputs/results.
# Small mutable metadata is frozen separately; completed BAMs are read in place.
REMOTE_SNAPSHOT = r'''
import datetime as dt, hashlib, json, os, pathlib, shutil, subprocess, sys, uuid
os.umask(0o077)
cfg=json.loads(sys.argv[1]); root=pathlib.Path(cfg['root']).resolve(strict=True)
if not root.is_dir(): raise RuntimeError('Project root missing')
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def contained(p, base):
 if p.is_symlink() or base not in p.resolve(strict=True).parents:raise RuntimeError('Unsafe source')
 return p
def safe_rel(s):
 p=pathlib.PurePosixPath(s)
 if not p.parts or p.is_absolute() or '..' in p.parts or '\\' in s or any(ord(c)<32 for c in s):raise RuntimeError('Unsafe relative path')
 return p
files={}; completed=0; warnings=[]
out=root/'outputs'
if out.is_dir():
 for rp in sorted(out.rglob('*receipt.json')):
  if rp.name not in ('receipt.json','workflow-receipt.json'):continue
  rel=rp.relative_to(out)
  if any(x.startswith('.') or 'partial' in x.lower() for x in rel.parts):continue
  contained(rp,root); receipt=json.loads(rp.read_text())
  if rp.name=='workflow-receipt.json':
   if receipt.get('status')!='PASS':continue
   artifacts=receipt.get('artifacts',{})
  else:
   if not receipt.get('completed_utc') or not receipt.get('fingerprint'):continue
   artifacts=receipt.get('outputs',{})
  if not artifacts:continue
  for name, value in artifacts.items():
   safe_rel(name); p=contained(rp.parent/name,root)
   if not p.is_file():raise RuntimeError('Receipt artifact is not a regular file')
   sha=value.get('sha256') if isinstance(value,dict) else value
   expected_size=value.get('bytes',value.get('size_bytes')) if isinstance(value,dict) else None
   if not isinstance(sha,str) or len(sha)!=64 or any(c not in '0123456789abcdef' for c in sha):raise RuntimeError('Invalid receipt digest')
   size=p.stat().st_size
   if expected_size is not None and size!=expected_size:raise RuntimeError('Completed artifact size changed')
   key=str(p.relative_to(root))
   if key in files and files[key]['Sha256']!=sha:raise RuntimeError('Conflicting completed receipts')
   files[key]={'Source':str(p),'SizeBytes':size,'Sha256':sha,'Receipt':str(rp.relative_to(root))}
  files[str(rp.relative_to(root))]={'Source':str(rp),'SizeBytes':rp.stat().st_size,'Sha256':digest(rp),'Receipt':str(rp.relative_to(root))}
  completed+=1
scratch=root/'.backup-snapshots'/('pending-'+uuid.uuid4().hex)
scratch.mkdir(parents=True,mode=0o700)
exclude_dirs={'.git','.local','.ssh','.aws','.azure','.config','__pycache__','node_modules','venv','.venv','cache','caches','vep','homo_sapiens','private'}
exclude_suffix={'.pem','.key','.p12','.pfx','.fastq','.fq','.bam','.cram','.fasta','.fa','.bwt','.sqlite','.db'}
def metadata(base, label, locks_only=False):
 if not base.exists():return
 for current, dirs, names in os.walk(base,followlinks=False):
  dirs[:]=[d for d in dirs if d not in exclude_dirs and not d.startswith('.') and not (pathlib.Path(current)/d).is_symlink()]
  for name in names:
   p=pathlib.Path(current)/name
   if p.is_symlink() or not p.is_file():continue
   if name.startswith('.') or p.suffix.lower() in exclude_suffix or any(x in name.lower() for x in ('credential','private-key','token','secret')):continue
   if locks_only and not (p.suffix in ('.json','.md','.yaml','.yml','.toml','.lock') and any(x in name.lower() for x in ('lock','receipt','version','source','manifest','provenance','ready'))):continue
   if p.stat().st_size>64*1024*1024:
    warnings.append('Oversized metadata skipped; manual review required');continue
   target=scratch/'metadata'/label/p.relative_to(base);target.parent.mkdir(parents=True,exist_ok=True,mode=0o700)
   before=p.stat()
   try:shutil.copyfile(p,target)
   except PermissionError:
    # Bootstrap may leave root-owned logs. Read only the already selected,
    # bounded metadata file; output remains owned by the ordinary SSH user.
    if root not in p.resolve(strict=True).parents:raise
    with target.open('wb') as stream:subprocess.run(['sudo','-n','cat','--',str(p)],stdout=stream,check=True)
   after=p.stat()
   if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
    target.unlink();warnings.append('Changing metadata skipped; next watch may capture it');continue
   target.chmod(0o600);key=str(target.relative_to(scratch))
   files[key]={'Source':str(target),'SizeBytes':target.stat().st_size,'Sha256':digest(target),'Receipt':None}
metadata(root/'private','private')
metadata(root/'logs','logs')
metadata(root/'resources','resource-locks',True)
metadata(root/'references','reference-locks',True)
if cfg.get('code'):
 code=pathlib.Path(cfg['code']).resolve(strict=True)
 if not code.is_dir() or len(code.parts)<4:raise RuntimeError('Use dedicated runtime code directory')
 metadata(code,'code')
manifest={'SchemaVersion':1,'CreatedUtc':dt.datetime.now(dt.timezone.utc).isoformat(),'CompletedStages':completed,'Files':[dict(RelativePath=k,**v) for k,v in sorted(files.items())],'Warnings':sorted(set(warnings))}
payload=json.dumps(manifest,sort_keys=True,separators=(',',':')).encode();sid=hashlib.sha256(payload).hexdigest()
(scratch/'manifest.json').write_bytes(payload);(scratch/'manifest.json').chmod(0o600)
print(json.dumps({'SnapshotId':sid,'Manifest':manifest}))
'''


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def local_path(value):
    if re.match(r'^[A-Za-z]:[\\/]', value):
        value = subprocess.check_output(['wslpath', '-u', value], text=True).strip()
    return Path(value).expanduser().resolve()


def relative(value):
    p = PurePosixPath(value)
    require(value and p.parts and not p.is_absolute() and '..' not in p.parts and '\\' not in value
            and not any(ord(c) < 32 for c in value), 'Unsafe snapshot path')
    return p


def atomic_json(path, value):
    temp = path.with_name(path.name + '.tmp')
    temp.write_bytes(canonical(value))
    os.replace(temp, path)


def status(args, destination, **fields):
    payload = {'UpdatedUtc': dt.datetime.now(dt.timezone.utc).isoformat(), **fields}
    atomic_json(destination / 'backup-status.json', payload)
    if args.status_file:
        target = local_path(args.status_file)
        target.parent.mkdir(parents=True, exist_ok=True)
        atomic_json(target, payload)


def promote_exclusive(source, target):
    libc = ctypes.CDLL(None, use_errno=True)
    require(hasattr(libc, 'renameat2'), 'Linux atomic no-replace renameat2 required')
    fn = libc.renameat2
    fn.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    fn.restype = ctypes.c_int
    if fn(-100, os.fsencode(source), -100, os.fsencode(target), 1):
        error = ctypes.get_errno()
        if error in (errno.EINVAL, errno.ENOSYS, errno.EOPNOTSUPP) and re.match(r'^/mnt/[a-z]/', str(source)):
            # DrvFS does not implement renameat2 flags. Directory.Move uses the
            # native Windows no-replace move; never substitute os.rename here.
            require(source.parent.parent == target.parent.parent, 'Promotion must remain within the backup root')
            paths = [subprocess.check_output(['wslpath', '-w', str(p)], text=True).strip()
                     for p in (source, target)]
            quoted = ["'" + p.replace("'", "''") + "'" for p in paths]
            ps = "$ErrorActionPreference='Stop'; [System.IO.Directory]::Move(" + ','.join(quoted) + ')'
            encoded = base64.b64encode(ps.encode('utf-16le')).decode('ascii')
            result = subprocess.run(['powershell.exe', '-NoProfile', '-NonInteractive',
                                     '-EncodedCommand', encoded], capture_output=True)
            if result.returncode:
                raise OSError(errno.EEXIST if target.exists() else errno.EIO,
                              'Native Windows no-replace snapshot promotion failed')
            return
        raise OSError(error, os.strerror(error), str(target))


def bounded_transfer(command, deadline):
    require(time.monotonic() < deadline, 'Backup deadline reached before next transfer')
    process = subprocess.Popen(command, start_new_session=True)
    try:
        code = process.wait(timeout=max(1, deadline - time.monotonic()))
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGINT)
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        raise RuntimeError('Backup time limit reached; verified objects and resumable partials remain')
    if code:
        raise subprocess.CalledProcessError(code, command)


def run_snapshot(args, ssh, destination):
    from write_local_index import publish_verified, write_index
    print('Discovering completed stages and freezing small metadata...', flush=True)
    command = shlex.join(['python3', '-', json.dumps({'root': args.remote_root, 'code': args.remote_code})])
    response = subprocess.run(ssh + [args.host, command], input=REMOTE_SNAPSHOT, text=True,
                              capture_output=True, timeout=max(1, args.deadline - time.monotonic()))
    if response.returncode:
        print(response.stderr, file=sys.stderr, flush=True)
        response.check_returncode()
    result = json.loads(response.stdout)
    manifest, sid = result['Manifest'], result['SnapshotId']
    require(digest_bytes(canonical(manifest)) == sid, 'Remote snapshot manifest digest mismatch')
    entries = manifest['Files']
    require(entries and manifest['CompletedStages'] > 0, 'No completed result stages found')
    objects, staging = destination / 'objects', destination / '.staging'
    snapshots = destination / 'snapshots'
    for p in (objects, staging, snapshots):
        p.mkdir(parents=True, exist_ok=True, mode=0o700)
    state_path = destination / 'verified-objects.json'
    state = json.loads(state_path.read_text()) if state_path.exists() else {}
    map_path = destination / 'verified-files.json'
    mapping = json.loads(map_path.read_text(encoding='utf-8')) if map_path.exists() else {}
    required = {}
    for entry in entries:
        relative(entry['RelativePath'])
        sha = entry['Sha256']
        require(re.fullmatch('[0-9a-f]{64}', sha) and isinstance(entry['SizeBytes'], int)
                and entry['SizeBytes'] >= 0, 'Invalid snapshot entry')
        source = PurePosixPath(entry['Source'])
        require(source.is_absolute() and not '..' in source.parts
                and str(source).startswith(args.remote_root.rstrip('/') + '/'), 'Source escaped remote project')
        obj = objects / sha
        if obj.exists():
            require(not obj.is_symlink() and obj.stat().st_size == entry['SizeBytes'], 'Existing object has wrong size')
            identity = {'size': obj.stat().st_size, 'mtime_ns': obj.stat().st_mtime_ns}
            if state.get(sha) != identity:
                require(digest(obj) == sha, 'Existing object checksum differs; refusing overwrite')
                state[sha] = identity
            publish_verified(destination, entry, mapping)
        else:
            if sha in required:
                require(required[sha]['SizeBytes'] == entry['SizeBytes'], 'Conflicting content address sizes')
            required[sha] = entry
    atomic_json(state_path, state)
    write_index(destination)
    required = dict(sorted(required.items(), key=lambda item:
        (not item[1]['RelativePath'].endswith('.bam'), -item[1]['SizeBytes'])))
    needed = sum(e['SizeBytes'] for e in required.values())
    free = shutil.disk_usage(destination).free
    reserve = int(args.reserve_gib * 2**30)
    require(free >= needed + reserve, 'Insufficient local free space including reserve; no data pulled')
    print(f'Plan: {manifest["CompletedStages"]} completed stages; {len(entries)} files; '
          f'{len(required)} new objects / {needed / 2**30:.3f} GiB; free {free / 2**30:.1f} GiB', flush=True)
    started, received = time.monotonic(), 0
    status(args, destination, Status='PULLING', CompletedStages=manifest['CompletedStages'],
           Files=len(entries), NewObjects=len(required), VerifiedNewObjects=0, NewBytes=needed, ReceivedBytes=0)
    for index, (sha, entry) in enumerate(required.items(), 1):
        part = staging / (sha + '.part')
        require(not part.is_symlink(), 'Symlink in local staging')
        partial_bytes = part.stat().st_size if part.is_file() else 0
        require(shutil.disk_usage(destination).free >= max(0, entry['SizeBytes'] - partial_bytes) + reserve,
                'Local disk reserve changed during backup; completed objects remain safe')
        cmd = ['rsync', '--partial', '--protect-args', '--timeout=300', '--info=progress2', '--human-readable']
        if part.exists():
            cmd.append('--append-verify')
        cmd += ['-e', shlex.join(ssh), '--', args.host + ':' + entry['Source'], str(part)]
        print(f'Pulling object {index}/{len(required)} ({entry["SizeBytes"] / 2**20:.1f} MiB)', flush=True)
        bounded_transfer(cmd, args.deadline)
        require(part.stat().st_size == entry['SizeBytes'] and digest(part) == sha,
                'Downloaded checksum failed; object remains unverified staging')
        obj = objects / sha
        try:
            os.link(part, obj)  # Atomic no-replace publication on the same filesystem.
        except FileExistsError:
            require(digest(obj) == sha, 'Differing final object exists; refusing overwrite')
        part.unlink()
        state[sha] = {'size': obj.stat().st_size, 'mtime_ns': obj.stat().st_mtime_ns}
        atomic_json(state_path, state)
        for selected in entries:
            if selected['Sha256'] == sha:
                publish_verified(destination, selected, mapping)
        write_index(destination)
        received += entry['SizeBytes']
        status(args, destination, Status='PULLING', CompletedStages=manifest['CompletedStages'],
               Files=len(entries), NewObjects=len(required), VerifiedNewObjects=index,
               NewBytes=needed, ReceivedBytes=received,
               MeanMiBPerSecond=round(received / 2**20 / max(time.monotonic()-started, .001), 2))
        print(f'Verified {index}/{len(required)} new objects; aggregate '
              f'{received / 2**20 / max(time.monotonic()-started, .001):.2f} MiB/s', flush=True)
    pending = staging / ('snapshot-' + sid)
    pending.mkdir(exist_ok=True, mode=0o700)
    for entry in entries:
        rel = relative(entry['RelativePath'])
        link = pending / 'files' / Path(*rel.parts)
        link.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        obj = objects / entry['Sha256']
        if link.exists():
            require(os.path.samefile(obj, link), 'Differing snapshot link; refusing overwrite')
        else:
            os.link(obj, link)  # No second full BAM copy.
    (pending / 'manifest.json').write_bytes(canonical(manifest))
    atomic_json(pending / 'backup-receipt.json', {
        'SchemaVersion': 1, 'Status': 'PASS', 'SnapshotId': sid,
        'VerifiedFiles': len(entries), 'TotalBytes': sum(e['SizeBytes'] for e in entries),
        'CompletedStages': manifest['CompletedStages'], 'NewBytes': received,
        'CompletedUtc': dt.datetime.now(dt.timezone.utc).isoformat()})
    final = snapshots / sid
    require(not final.exists(), 'Snapshot already exists; refusing replacement')
    promote_exclusive(pending, final)
    atomic_json(destination / 'latest.json', {'SnapshotId': sid, 'Path': str(final), 'Status': 'PASS'})
    print(f'PASS: {len(entries)} files verified; snapshot {sid}', flush=True)
    status(args, destination, Status='PASS', CompletedStages=manifest['CompletedStages'],
           VerifiedFiles=len(entries), TotalBytes=sum(e['SizeBytes'] for e in entries),
           NewBytes=received, SnapshotId=sid)
    for warning in manifest['Warnings']:
        print('Metadata notice: ' + warning, flush=True)


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', required=True, help='Confirmed SSH config alias')
    parser.add_argument('--ssh-config', required=True)
    parser.add_argument('--known-hosts', required=True)
    parser.add_argument('--remote-root', required=True)
    parser.add_argument('--remote-code', help='Optional dedicated code-only runtime directory')
    parser.add_argument('--ssh-hostname', help='Optional confirmed direct host; requires existing host-key alias')
    parser.add_argument('--ssh-port', type=int, default=22)
    parser.add_argument('--host-key-alias', help='Existing pinned known-hosts identity for the same server')
    parser.add_argument('--destination', required=True, help='Local private backup outside Git')
    parser.add_argument('--allow-onedrive', action='store_true', help='Use an explicitly user-selected OneDrive destination')
    parser.add_argument('--status-file', help='Optional private aggregate status JSON outside Git')
    parser.add_argument('--reserve-gib', type=float, default=15)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--watch', action='store_true')
    parser.add_argument('--interval', type=int, default=60)
    parser.add_argument('--max-hours', type=float, default=8)
    args = parser.parse_args()
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.host), 'Invalid SSH alias')
    root = PurePosixPath(args.remote_root)
    require(root.is_absolute() and '..' not in root.parts and len(root.parts) >= 3
            and not any(ord(c) < 32 for c in args.remote_root), 'Use dedicated absolute remote project root')
    require(args.interval >= 60 and 0 < args.max_hours <= 8 and args.reserve_gib >= 15,
            'Watch interval >=60 seconds, max-hours <=8, reserve >=15 GiB required')
    require(not args.ssh_hostname or (re.fullmatch(r'[A-Za-z0-9_.:-]+', args.ssh_hostname)
            and args.host_key_alias and not any(ord(c) < 33 for c in args.host_key_alias)
            and 1 <= args.ssh_port <= 65535), 'Direct transport requires confirmed host and existing pinned key alias')
    destination = local_path(args.destination)
    require(args.allow_onedrive or 'onedrive' not in str(destination).lower(), 'Explicit --allow-onedrive needed for user-selected OneDrive destination')
    for p in (destination, *destination.parents):
        require(not (p / '.git').exists(), 'Backup destination cannot be inside a Git checkout')
    print('Backup policy: completed stages only; no raw reads/public caches/references; 15+ GiB reserve.', flush=True)
    if not args.execute:
        print('DRY RUN: no SSH, remote snapshots, directories or downloads created.', flush=True)
        return
    import fcntl
    os.umask(0o077)
    for name in ('ssh', 'rsync'):
        require(shutil.which(name), name + ' required inside WSL; script installs nothing')
    config, hosts = local_path(args.ssh_config), local_path(args.known_hosts)
    require(config.is_file() and hosts.is_file(), 'Existing SSH config and pinned known-hosts required')
    ssh = ['ssh', '-T', '-F', str(config), '-o', 'BatchMode=yes', '-o', 'StrictHostKeyChecking=yes',
           '-o', 'UserKnownHostsFile=' + str(hosts), '-o', 'ConnectTimeout=30',
           '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=4',
           '-o', 'ControlMaster=no', '-o', 'ControlPath=none', '-o', 'IPQoS=none']
    if args.ssh_hostname:
        ssh += ['-o', 'HostName=' + args.ssh_hostname, '-o', 'Port=' + str(args.ssh_port),
                '-o', 'HostKeyAlias=' + args.host_key_alias]
    destination.mkdir(parents=True, exist_ok=True, mode=0o700)
    lock = destination / '.backup.lock'
    fd = os.open(lock, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        os.close(fd)
        raise RuntimeError('A backup process already holds the destination lock')
    os.ftruncate(fd, 0)
    os.write(fd, str(os.getpid()).encode())
    deadline = time.monotonic() + args.max_hours * 3600
    args.deadline = deadline
    try:
        while True:
            for attempt in range(1, 4):
                try:
                    run_snapshot(args, ssh, destination)
                    break
                except subprocess.CalledProcessError as exc:
                    transient = exc.returncode in (255, 10, 12, 30, 35)
                    if not args.watch or not transient or attempt == 3 or time.monotonic() + 60 >= deadline:
                        status(args, destination, Status='FAIL', Reason='SSH/rsync failed', ExitCode=exc.returncode)
                        raise
                    status(args, destination, Status='RETRY_WAIT', Attempt=attempt, ExitCode=exc.returncode)
                    print(f'Transient SSH/rsync exit {exc.returncode}; retry {attempt}/3 after 60 seconds', flush=True)
                    time.sleep(60)
            if not args.watch or time.monotonic() + args.interval >= deadline:
                break
            time.sleep(args.interval)
    except Exception as exc:
        status(args, destination, Status='FAIL', Reason=type(exc).__name__,
               Detail='Review the private progress log; completed snapshots remain intact')
        raise
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, ValueError, OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        print('BACKUP FAILED: ' + str(exc), file=sys.stderr)
        sys.exit(1)
