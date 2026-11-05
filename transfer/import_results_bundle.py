#!/usr/bin/env python3
"""Import a completed size-checked TAR against its separate pinned-SSH manifest.

Only regular numbered members explicitly enumerated by the trusted manifest are
accepted. Does not call extractall, compute patient checksums or replace files.
"""
import argparse
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from import_verified_download import path
from write_local_index import size_record_valid, size_key, publish_size_only, publish_verified, atomic_json, write_index
from watch_native_imports import verified


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--destination', required=True)
    ap.add_argument('--allow-onedrive', action='store_true')
    ap.add_argument('--completed-download', action='store_true')
    ap.add_argument('--stable-seconds', type=int, default=60)
    ap.add_argument('--wait-max-hours', type=float, default=0, help='Bounded local-only wait for the planned completed native file (max8h)')
    ap.add_argument('--execute', action='store_true')
    args = ap.parse_args()
    source, manifest, destination = map(path, (args.source, args.manifest, args.destination))
    if any((p/'.git').exists() for p in (destination,*destination.parents)):
        raise ValueError('Destination must be outside Git')
    if 'onedrive' in str(destination).lower() and not args.allow_onedrive:
        raise ValueError('Explicit OneDrive destination authorization flag required')
    plan = json.loads(manifest.read_text())
    if not 0 <= args.wait_max_hours <= 8:
        raise ValueError('Bounded wait must be between0 and8hours')
    if args.wait_max_hours:
        deadline=time.monotonic()+args.wait_max_hours*3600
        while not (source.is_file() and source.stat().st_size==plan['SizeBytes'] and time.time()-source.stat().st_mtime>=args.stable_seconds and not any(Path(str(source)+s).exists() for s in ('.part','.crdownload','.tmp'))):
            if (destination/'Settings/LOCAL_DOWNLOADS_PAUSED.json').exists():
                print('PAUSED by user request; native files preserved');return
            if time.monotonic()>=deadline:
                raise TimeoutError('No completed bundle before local wait deadline; partials preserved')
            time.sleep(min(30,max(0,deadline-time.monotonic())))
    if not args.completed_download or source.suffix != '.tar' or source.stat().st_size != plan['SizeBytes']:
        raise ValueError('Completed final TAR matching expected SSH size required')
    if args.stable_seconds < 0 or time.time()-source.stat().st_mtime < args.stable_seconds:
        raise ValueError('Final archive must have stable modification time')
    if any(Path(str(source)+s).exists() for s in ('.part', '.crdownload', '.tmp')):
        raise ValueError('Partial sibling exists')
    expected = {}
    for entry in plan['Files']:
        member, rel = entry['ArchiveMember'], PurePosixPath(entry['RelativePath'])
        if not re.fullmatch(r'files/[0-9]{5}', member) or member in expected:
            raise ValueError('Unsafe or duplicate archive member')
        if not rel.parts or rel.is_absolute() or '..' in rel.parts or '\\' in str(rel) or rel.parts[0] not in ('outputs', 'metadata'):
            raise ValueError('Unsafe intended relative path')
        expected[member] = entry
    with tarfile.open(source, 'r:') as archive:
        members = archive.getmembers()
        if len(members) != len(expected) or {m.name for m in members} != set(expected):
            raise ValueError('Archive membership differs from trusted SSH manifest')
        if any(not m.isfile() or m.size != expected[m.name]['SizeBytes'] for m in members):
            raise ValueError('Only expected regular files of exact size allowed')
        if not args.execute:
            print(f'DRY RUN SIZE-ONLY: {len(members)} safe expected members')
            return
        if (destination/'Settings/LOCAL_DOWNLOADS_PAUSED.json').exists():
            print('PAUSED by user request; no extraction started');return
        if shutil.disk_usage(destination).free < plan['PayloadBytes'] + 15*2**30:
            raise ValueError('Need payload bytes plus15GiB reserve')
        os.umask(0o077)
        import fcntl
        lock=(destination/'.bundle-import.lock').open('a')
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        controller_lock=(destination/'.backup.lock').open('a')
        fcntl.flock(controller_lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        state_path=destination/'size-only-files.json'
        records=json.loads(state_path.read_text(encoding='utf-8-sig')) if state_path.exists() else {}
        sha_state=json.loads((destination/'verified-objects.json').read_text()) if (destination/'verified-objects.json').exists() else {}
        sha_mapping=json.loads((destination/'verified-files.json').read_text()) if (destination/'verified-files.json').exists() else {}
        stage = destination / '.staging'
        stage.mkdir(exist_ok=True)
        trusted = stage / ('bundle-plan-' + source.stem + '.json')
        # Keep exact plan as local provenance; no overwriting different metadata.
        payload = json.dumps(plan, sort_keys=True)
        if trusted.exists() and trusted.read_text() != payload:
            raise ValueError('Differing existing bundle plan')
        if not trusted.exists():
            trusted.write_text(payload)
        imported=0
        def priority(member):
            name=PurePosixPath(expected[member.name]['RelativePath']).name
            return (0 if name.endswith('.bai') else 1 if name.startswith('cohort.') else 2,member.name)
        for member in sorted(members,key=priority):
            if (destination/'Settings/LOCAL_DOWNLOADS_PAUSED.json').exists():
                atomic_json(state_path,records);write_index(destination)
                print('PAUSED by user request; completed imports and source archive retained');return
            entry = expected[member.name]
            obj=destination/'objects'/entry['Sha256']
            if obj.is_file() and not obj.is_symlink() and sha_state.get(entry['Sha256'])=={'size':obj.stat().st_size,'mtime_ns':obj.stat().st_mtime_ns} and obj.stat().st_size==entry['SizeBytes']:
                mapkey=entry['Sha256']+':'+entry['RelativePath']
                if mapkey not in sha_mapping or not (destination/sha_mapping[mapkey]['LocalPath']).is_file():
                    publish_verified(destination,entry,sha_mapping)
                continue
            existing=records.get(size_key(entry))
            if existing and (destination/existing['StoragePath']).is_file() and (destination/existing['LocalPath']).is_file() and (destination/existing['StoragePath']).stat().st_size==entry['SizeBytes'] and (destination/existing['StoragePath']).stat().st_mtime_ns==existing['MtimeNs']:
                continue
            fd,tempname=tempfile.mkstemp(prefix=source.stem+'-'+member.name.split('/')[1]+'-',dir=stage)
            temp=Path(tempname)
            with archive.extractfile(member) as src, os.fdopen(fd,'wb') as out:
                shutil.copyfileobj(src, out)
            if temp.stat().st_size != entry['SizeBytes']:
                raise ValueError('Extracted member size mismatch')
            folder=destination/'size-only'/size_key(entry);folder.mkdir(exist_ok=True,parents=True)
            stored=folder/PurePosixPath(entry['RelativePath']).name
            if stored.exists():
                if stored.is_symlink() or stored.stat().st_size!=entry['SizeBytes']:
                    raise ValueError('Differing existing size-only file; refusing overwrite')
            else:
                os.link(temp,stored)
            publish_size_only(destination,entry,stored,'native_receipt_completed_bundle_exact_member_size',records=records,persist=False)
            imported+=1
            if imported%25==0:
                atomic_json(state_path,records)
                print(f'SIZE ONLY bundle: {imported} newly imported files',flush=True)
            temp.unlink()  # Store/friendly hardlinks remain; source archive remains untouched.
        atomic_json(state_path,records)
        count=write_index(destination)
        atomic_json(destination/'backup-status.json',{'Status':'PARTIAL_BACKUP','VerificationMode':'size_only','CompletedFileLinks':count,'CompletedSnapshot':False,'FullSnapshotSHA256Pass':False,'ActualSha256':None,'BundleImportedFiles':len(expected),'Note':'New bundle checked by native completion and exact member sizes; prior SHA history retained'})
    print(json.dumps({'Status': 'IMPORTED', 'Files': len(expected), 'Verification': 'size_only', 'ActualSha256': None}))


if __name__ == '__main__':
    main()
