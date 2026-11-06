#!/usr/bin/env python3
"""Receipt-bound final cohort VCF/TBI download; existing native trust, SIZE ONLY."""
import argparse
import base64
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import sys

END=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
REMOTE='/srv/marina/outputs/fully-annotated-vcf-v1/'
NAME='cohort.fully-annotated.vcf.gz'
RESERVE=15*1024**3


def need(ok,message):
    if not ok:raise ValueError(message)


def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def captured(item,filename):
    need(set(item)=={'path','bytes','sha256','content_base64'} and item['path']==REMOTE+filename,'Exact captured source metadata required')
    raw=base64.b64decode(item['content_base64'],validate=True)
    need(len(raw)<=32*1024**2 and len(raw)==item['bytes'] and hashlib.sha256(raw).hexdigest()==item['sha256'],'Captured metadata SHA/size differs')
    return json.loads(raw)


def validate(plan):
    keys={'schema_version','status','scientific_status','deadline_utc','destination','state_directory','pause_file','source_provenance','files'}
    need(set(plan)==keys and plan['schema_version']==1 and plan['status']=='ROOT_REVIEWED_FINAL_ANNOTATED_VCF_DOWNLOAD','Exact reviewed final-pair plan required')
    need(plan['scientific_status']=='PASS_FULLY_ANNOTATED_VCF_DIRECT','Final all-source scientific status required')
    deadline=dt.datetime.fromisoformat(plan['deadline_utc'])
    need(deadline.tzinfo is not None and deadline<=END,'Original absolute deadline cannot increase')
    provenance=plan['source_provenance'];need(set(provenance)=={'manifest','receipt'},'Exactly final manifest and receipt required')
    manifest=captured(provenance['manifest'],'manifest.json');receipt=captured(provenance['receipt'],'workflow-receipt.json')
    need(manifest.get('schema_version')==1 and manifest.get('status')=='PASS_FULLY_ANNOTATED_VCF_DIRECT' and manifest.get('final_delivery_complete') is True,'Captured producer is not final PASS')
    need(all(manifest.get(k) is True for k in ('original_records_fields_preserved','packet_decode_verified','direct_core_prefix_verified','core33_packets_byte_identical')),'Final source preservation proofs required')
    need(receipt.get('schema_version')==1 and receipt.get('status')=='PASS' and receipt.get('stage')=='fully-annotated-vcf-v1' and re.fullmatch('[0-9a-f]{64}',receipt.get('identity','')) is not None,'Authentic final completion receipt required')
    need(receipt.get('artifacts',{}).get('manifest.json')=={k:provenance['manifest'][k] for k in ('bytes','sha256')},'Final manifest not bound to completion receipt')
    need(len(plan['files'])==2 and {e['DownloadFileName'] for e in plan['files']}=={NAME,NAME+'.tbi'},'Exactly canonical final VCF/TBI pair required')
    for entry in plan['files']:
        name=entry['DownloadFileName'];relative='outputs/fully-annotated-vcf-v1/'+name
        need(entry['RelativePath']==relative and entry['ServerPath']==REMOTE+name and entry['SourceReceipt']=='outputs/fully-annotated-vcf-v1/workflow-receipt.json','Wrong completed pair source path')
        need(type(entry['SizeBytes']) is int and 0<entry['SizeBytes']<=32*1024**3 and re.fullmatch('[0-9a-f]{64}',entry['Sha256']) is not None,'Invalid final pair expected identity')
        expected=manifest['index' if name.endswith('.tbi') else 'enriched_vcf']
        need(expected==dict(path=entry['ServerPath'],bytes=entry['SizeBytes'],sha256=entry['Sha256']),'Pair differs from final manifest')
        need(receipt.get('artifacts',{}).get(name)==dict(bytes=entry['SizeBytes'],sha256=entry['Sha256']),'Pair differs from final completion receipt')
    return deadline


def remaining_bytes(plan,destination):
    remaining=0
    for entry in plan['files']:
        target=destination/entry['DownloadFileName'];partial=target.with_name(target.name+'.part');offset=0
        need(not target.is_symlink() and not partial.is_symlink(),'Symlink download target refused')
        if target.exists():
            evidence=target.with_name(target.name+'.completed.json')
            need(target.is_file() and target.stat().st_size==entry['SizeBytes'] and evidence.is_file(),'Unrelated final file; refusing overwrite')
            offset=entry['SizeBytes']
        elif partial.exists():
            evidence=partial.with_name(partial.name+'.source.json');need(evidence.is_file(),'Unbound partial retained; no overwrite')
            identity=load(evidence)
            need(all(identity.get(k)==v for k,v in dict(RelativePath=entry['RelativePath'],RemotePath=entry['ServerPath'],ExpectedBytes=entry['SizeBytes'],ExpectedServerSha256=entry['Sha256'],SourceReceipt=entry['SourceReceipt']).items()),'Partial belongs to another source')
            offset=partial.stat().st_size;need(offset<=entry['SizeBytes'],'Partial exceeds expected size')
        remaining+=entry['SizeBytes']-offset
    return remaining


def check_space(plan,destination,free_bytes=None):
    pending=remaining_bytes(plan,destination)
    free=shutil.disk_usage(destination).free if free_bytes is None else free_bytes
    need(free>=RESERVE+pending,'Insufficient space above unchanged 15GiB reserve')
    return pending


def run(plan,auth_path=None,execute=False):
    deadline=validate(plan)
    if not execute:return dict(status='PLAN_ONLY_NO_AUTH_KEY_NETWORK_OR_FILES',files=2,bytes=sum(e['SizeBytes'] for e in plan['files']),scientific_status=plan['scientific_status'])
    need(os.name=='nt','Use the existing native Windows SFTP runtime')
    need(dt.datetime.now(dt.timezone.utc)<deadline,'Original download deadline exhausted')
    destination,state=Path(plan['destination']),Path(plan['state_directory'])
    for folder in (destination,state):
        need(folder.is_absolute() and not any(p.is_symlink() or (p/'.git').exists() for p in (folder,*folder.parents)),'Private absolute paths outside Git/symlinks required')
    need(destination.is_dir(),'Existing user-selected VCF directory required')
    check_space(plan,destination)
    need(auth_path is not None,'Existing auth configuration reference required')
    state.mkdir(parents=True,exist_ok=True)
    import msvcrt
    lock=(state/'final-pair.lock').open('a+b');lock.write(b'0');lock.flush();lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
    spec=importlib.util.spec_from_file_location('reviewed_final_pair_transfer',Path(__file__).with_name('sync_completed_small_results.py'));helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    client=sftp=None
    def status(label,**extra):helper.atomic(state/'status.json',dict(status=label,deadline_utc=deadline.isoformat(),updated_utc=helper.now().isoformat(),patient_hashes_computed=False,scientific_status=plan['scientific_status'],**extra))
    try:
        import paramiko
        need(paramiko.__version__=='5.1.0','Existing reviewed Paramiko5.1.0 runtime required')
        from native_sftp_download import connect,blocks
        client,sftp=connect(load(auth_path));sftp.get_channel().settimeout(min(60,max(1,(deadline-helper.now()).total_seconds())))
        for entry in sorted(plan['files'],key=lambda e:e['SizeBytes']):
            helper.download_one(sftp,blocks,entry,destination,deadline,Path(plan['pause_file']),reserve=RESERVE,
                progress=lambda done,total:status('DOWNLOADING_FINAL_ANNOTATED_PAIR_SIZE_ONLY',filename=entry['DownloadFileName'],received_bytes=done,expected_bytes=total))
        helper.atomic(destination/'cohort.fully-annotated.download.json',dict(status='COMPLETED_SIZE_ONLY',scientific_status=plan['scientific_status'],source_provenance=plan['source_provenance'],files=plan['files'],patient_hashes_computed=False,completed_utc=helper.now().isoformat()))
        status('PASS_FINAL_ANNOTATED_PAIR_COPIED_SIZE_ONLY',files=2)
        return dict(status='PASS_FINAL_ANNOTATED_PAIR_COPIED_SIZE_ONLY',files=2)
    except BaseException as error:
        status('FAILED_PARTIALS_RETAINED_NO_AUTOMATIC_RESTART',error_type=type(error).__name__);raise
    finally:
        if sftp is not None:sftp.close()
        if client is not None:client.close()
        lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--auth-config',type=Path);parser.add_argument('--execute',action='store_true');args=parser.parse_args()
    try:print(json.dumps(run(load(args.plan),args.auth_config,args.execute)))
    except BaseException as error:
        print('FAILED: '+type(error).__name__+'; resumable files retained; private connection values not logged',file=sys.stderr);sys.exit(1)
