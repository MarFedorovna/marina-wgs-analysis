#!/usr/bin/env python3
"""Save an explicitly reviewed completed VEP component without claiming stage06 PASS."""
import argparse
import datetime as dt
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import sys


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def run(plan, auth_path=None, execute=False):
    required = {'schema_version', 'status', 'scientific_status', 'deadline_utc', 'destination',
                'state_directory', 'pause_file', 'source_provenance', 'files'}
    if set(plan) != required or plan['schema_version'] != 1 or plan['status'] != 'ROOT_REVIEWED_COMPLETED_COMPONENT_DOWNLOAD':
        raise ValueError('An exact reviewed completed-component plan is required')
    if plan['scientific_status'] != 'COMPLETED_VEP_COMPONENT_NOT_FINAL_ALL_SOURCES_OR_REPORT':
        raise ValueError('This download must not imply final all-source annotation or legacy stage completion')
    if len(plan['files']) != 2 or sorted(f['DownloadFileName'].endswith('.tbi') for f in plan['files']) != [False, True]:
        raise ValueError('Exactly one VEP VCF and its index are required')
    for entry in plan['files']:
        if (Path(entry['DownloadFileName']).name != entry['DownloadFileName'] or
                any(c in entry['DownloadFileName'] for c in '/\\:') or
                not entry['DownloadFileName'].endswith(('.vcf.gz', '.vcf.gz.tbi')) or
                type(entry['SizeBytes']) is not int or not 0 < entry['SizeBytes'] < 2 * 1024**3 or
                re.fullmatch('[a-f0-9]{64}', entry['Sha256']) is None):
            raise ValueError('Invalid completed VEP component file')
    source = next(f for f in plan['files'] if not f['DownloadFileName'].endswith('.tbi'))
    index = next(f for f in plan['files'] if f['DownloadFileName'].endswith('.tbi'))
    if index['ServerPath'] != source['ServerPath'] + '.tbi' or index['DownloadFileName'] != source['DownloadFileName'] + '.tbi':
        raise ValueError('The VCF and index must be the same completed component pair')
    if not execute:
        return dict(status='PLAN_ONLY_NO_AUTH_OR_NETWORK', files=len(plan['files']), bytes=sum(e['SizeBytes'] for e in plan['files']))
    if os.name != 'nt':
        raise ValueError('Use the existing native Windows SFTP runtime')
    deadline = dt.datetime.fromisoformat(plan['deadline_utc'])
    if deadline.tzinfo is None or deadline <= dt.datetime.now(dt.timezone.utc):
        raise ValueError('The reviewed download deadline has expired')
    destination, state = Path(plan['destination']), Path(plan['state_directory'])
    for folder in [destination, state]:
        if any(p.is_symlink() for p in [folder, *folder.parents]) or any((p / '.git').exists() for p in [folder, *folder.parents]):
            raise ValueError('Private download paths must be outside Git and symlinks')
        folder.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(destination).free < 15 * 1024**3 + sum(e['SizeBytes'] for e in plan['files']):
        raise ValueError('15 GiB of local free space must remain')
    import msvcrt
    lock = (state / 'component.lock').open('a+b'); lock.write(b'0'); lock.flush(); lock.seek(0)
    msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
    spec = importlib.util.spec_from_file_location('reviewed_component_transfer', Path(__file__).with_name('sync_completed_small_results.py'))
    helper = importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
    from native_sftp_download import connect, blocks
    import paramiko
    if paramiko.__version__ != '5.1.0':
        raise ValueError('Use the existing reviewed Paramiko runtime')
    def status(label, **extra):
        helper.atomic(state / 'status.json', dict(status=label, updated_utc=helper.now().isoformat(), deadline_utc=deadline.isoformat(),
                      patient_hashes_computed=False, final_all_source_annotation=False, **extra))
    client = sftp = None
    try:
        client, sftp = connect(load(auth_path)); sftp.get_channel().settimeout(60)
        for entry in sorted(plan['files'], key=lambda e:e['SizeBytes']):
            status('DOWNLOADING_COMPLETED_VEP_COMPONENT', filename=entry['DownloadFileName'])
            helper.download_one(sftp, blocks, entry, destination, deadline, Path(plan['pause_file']),
                progress=lambda done, total:status('DOWNLOADING_COMPLETED_VEP_COMPONENT', filename=entry['DownloadFileName'], received_bytes=done, expected_bytes=total))
        metadata = destination / 'cohort.vep115.annotated.component.json'
        helper.atomic(metadata, dict(status='COMPLETED_SIZE_ONLY', scientific_status=plan['scientific_status'],
                      source_provenance=plan['source_provenance'], files=plan['files'], patient_hashes_computed=False,
                      final_all_source_annotation=False, completed_utc=helper.now().isoformat()))
        status('PASS_COMPLETED_VEP_COMPONENT_COPIED_SIZE_ONLY', files=2)
        return dict(status='PASS_COMPLETED_VEP_COMPONENT_COPIED_SIZE_ONLY', files=2)
    except BaseException as error:
        status('FAILED_PARTIALS_RETAINED', error_type=type(error).__name__)
        raise
    finally:
        if sftp is not None:sftp.close()
        if client is not None:client.close()
        lock.seek(0); msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1); lock.close()


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--auth-config',type=Path);parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    try:print(json.dumps(run(load(args.plan),args.auth_config,args.execute)))
    except BaseException as error:
        print('FAILED: '+type(error).__name__+'; completed files and resumable partials retained; no private connection values logged',file=sys.stderr)
        sys.exit(1)
