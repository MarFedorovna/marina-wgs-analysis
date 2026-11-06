#!/usr/bin/env python3
"""Independent native SFTP result sync, reusing reviewed trust and SIZE ONLY import.

Default is local inventory only: no auth config, key, network, directories or
importer are opened. Existing BAM writers and their files are never controlled.
"""
import argparse
import base64
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import re
import shutil
import subprocess
import sys
import time

DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
SIDE_ROOTS={'consequence-tables-v1','tissue-protein-context-v1','population-reference-context-v1',
            'combined-annotation-master-v1','cohort-sample-views-v1','mitochondrial-reference-context-v1'}
REQUIRED_ROOTS={'cohort-sample-views-v1'}
SUFFIXES=('.vcf.gz','.vcf.gz.tbi','.vcf.gz.csi','.tsv','.csv','.tsv.gz','.csv.gz','.txt','.json','.log','.md','.html','.pdf','.docx','.xlsx','.png','.jpg','.jpeg','.svg')


class SyncError(ValueError):pass


def require(value,message):
    if not value:raise SyncError(message)


def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def now():return dt.datetime.now(dt.timezone.utc)


def atomic(path,value):
    path=Path(path);temp=path.with_name(path.name+'.tmp')
    require(not path.is_symlink() and not temp.is_symlink(),'Unsafe sync status path')
    temp.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8');temp.replace(path)


def selected(relative):
    rel=PurePosixPath(relative);parts=rel.parts;name=rel.name.lower()
    require(parts and not rel.is_absolute() and '..' not in parts and rel.as_posix()==relative and '\\' not in relative,'Unsafe planned result path')
    if len(parts)<3 or parts[0]!='outputs' or name.endswith(('.bam','.bai')) or '.g.vcf' in name:return False
    if any(p.lower() in {'resources','databases','cache','raw','inputs','references'} or p.startswith('.') for p in parts):return False
    if not name.endswith(SUFFIXES):return False
    if parts[1] in SIDE_ROOTS:return True
    if parts[1]=='workflow-v1':
        return len(parts)>=5 and parts[2]=='stages' and parts[3] in {'04_joint','05_filter','06_annotate','07_report'}
    if parts[1]=='gene-coverage-v1':return True
    if parts[1]=='mitochondria-recovery-v2':
        return len(parts)==6 and parts[-2]=='calls' and name in {
            'mitochondria.filtered.vcf.gz','mitochondria.filtered.vcf.gz.tbi','mitochondria.pass.vcf.gz','mitochondria.pass.vcf.gz.tbi',
            'mitochondria.raw.vcf.gz','mitochondria.raw.vcf.gz.tbi','research-report.json','numt-assessment.json','contamination.json','heteroplasmy.tsv','circular-coverage.tsv'}
    return False


def receipt_complete(relative,value):
    parts=PurePosixPath(relative).parts
    if len(parts)<3 or parts[0]!='outputs':return False
    root=parts[1]
    if root=='tissue-protein-context-v1':return value.get('status')=='PASS_SEPARATE_RESEARCH_REFERENCE_CONTEXT_EXPORT' and value.get('original_pipeline_changed') is False and value.get('network_requests')==0
    if root in {'cohort-sample-views-v1','consequence-tables-v1'}:return value.get('status')=='PASS' and value.get('schema_version')==1 and value.get('synthetic') is False
    if root in {'population-reference-context-v1','combined-annotation-master-v1','mitochondrial-reference-context-v1'}:return value.get('status')=='PASS' and value.get('schema_version')==1 and value.get('stage')==root and bool(value.get('identity'))
    return bool(value.get('status')=='PASS' and value.get('identity') or value.get('completed_utc') and value.get('fingerprint'))


def inventory(plan,profile):
    require(plan.get('SchemaVersion')==1 and plan.get('PlanGeneratorVersion',0)>=2,'Reviewed native sidecar plan v4 required')
    receipts={}
    for receipt in plan.get('Receipts',[]):
        raw=base64.b64decode(receipt['ContentBase64'],validate=True)
        require(len(raw)==receipt['SizeBytes'] and hashlib.sha256(raw).hexdigest()==receipt['Sha256'],'Source receipt metadata integrity differs')
        require(receipt['RelativePath'] not in receipts,'Duplicate source receipt metadata')
        receipts[receipt['RelativePath']]=(receipt,json.loads(raw))
    files=[];excluded=0;seen=set();names=set();total=0
    remote_root=PurePosixPath(profile['remote_root'])
    require(remote_root.is_absolute() and '..' not in remote_root.parts,'Dedicated server result root required')
    for entry in plan.get('Files',[]):
        relative=entry['RelativePath']
        if not selected(relative):excluded+=1;continue
        require(relative not in seen and entry['DownloadFileName'] not in names,'Duplicate planned small result/name')
        seen.add(relative);names.add(entry['DownloadFileName'])
        require(Path(entry['DownloadFileName']).name==entry['DownloadFileName'] and not any(c in entry['DownloadFileName'] for c in '/\\:') and not any(ord(c)<32 for c in entry['DownloadFileName']),'Unsafe planned download basename')
        require(type(entry['SizeBytes']) is int and entry['SizeBytes']>=0 and re.fullmatch('[0-9a-f]{64}',entry['Sha256']),'Invalid expected result identity')
        require(entry['ServerPath']==str(remote_root/relative),'Server result path differs from trusted relative identity')
        require(entry['SourceReceipt'] in receipts,'Selected result lacks captured completion receipt')
        receipt,source=receipts[entry['SourceReceipt']]
        require(receipt_complete(entry['SourceReceipt'],source),'Selected result source is not a recognized completed producer')
        source_path=PurePosixPath(entry['SourceReceipt'])
        if relative==entry['SourceReceipt']:
            expected=dict(bytes=receipt['SizeBytes'],sha256=receipt['Sha256'])
        elif source.get('status')=='PASS_SEPARATE_RESEARCH_REFERENCE_CONTEXT_EXPORT':
            expected=source.get('artifact',{})
            require(expected.get('path')==entry['ServerPath'],'Tissue artifact path differs')
        else:
            require(source.get('status')=='PASS' or (source.get('completed_utc') and source.get('fingerprint')),'Result producer is not complete')
            artifact_name=PurePosixPath(relative).relative_to(source_path.parent).as_posix()
            expected=source.get('artifacts',source.get('outputs',{})).get(artifact_name,{})
            if isinstance(expected,str):expected=dict(sha256=expected,bytes=entry['SizeBytes'])
        require(expected.get('sha256')==entry['Sha256'] and expected.get('bytes',expected.get('size_bytes'))==entry['SizeBytes'],'Planned result differs from captured producer receipt')
        total+=entry['SizeBytes'];files.append(entry)
    require(total<=profile['max_total_bytes'] and all(e['SizeBytes']<=profile['max_file_bytes'] for e in files),'Selected result byte cap exceeded; nothing silently skipped')
    completed_roots={PurePosixPath(path).parts[1] for path,item in receipts.items() if receipt_complete(path,item[1])}
    missing=sorted(set(profile['required_completed_roots'])-completed_roots)
    return sorted(files,key=lambda e:(e['SizeBytes'],e['RelativePath'])),dict(selected_files=len(files),selected_bytes=total,excluded_files=excluded,
        missing_completed_roots=missing,individual_annotated_vcf_artifacts=sum('/cohort-sample-views-v1/vcf/' in e['RelativePath'] for e in files),patient_file_hashes_computed=False)


def validate_profile(profile):
    keys={'schema_version','plan','auth_config','download_directory','destination','state_directory','remote_root','wsl_distribution',
          'max_hours','interval_seconds','deadline_utc','max_file_bytes','max_total_bytes','min_free_bytes','required_completed_roots','allow_onedrive','pause_file'}
    require(set(profile)==keys and profile['schema_version']==1,'Unknown/missing small-result profile fields')
    require(profile['wsl_distribution']=='Ubuntu-22.04','Reviewed existing WSL distribution required')
    require(0<profile['max_hours']<=8 and profile['interval_seconds']>=60,'Bounded sync runtime required')
    require(set(profile['required_completed_roots'])==REQUIRED_ROOTS,'Exact final VCF-first cohort views completion root required')
    require(0<profile['max_file_bytes']<=profile['max_total_bytes']<=128*1024**3 and profile['min_free_bytes']>=15*1024**3,'Bounded result sizes and 15GiB reserve required')
    due=dt.datetime.fromisoformat(profile['deadline_utc'])
    require(due.tzinfo is not None and due<=DEADLINE,'Original absolute VM deadline cannot increase; UTC offset required')
    for name in ('download_directory','destination','state_directory'):
        path=Path(profile[name]).resolve()
        require(not any((p/'.git').exists() for p in (path,*path.parents)),'Patient sync data/state must remain outside Git')
        require(not any(p.is_symlink() for p in (path,*path.parents)),'Symlink sync destination refused')
        require(profile['allow_onedrive'] or 'onedrive' not in str(path).lower(),'Explicit user-selected OneDrive destination required')


def wsl_path(path):
    path=PureWindowsPath(str(path));require(re.fullmatch('[A-Za-z]:',path.drive) and path.is_absolute(),'Native absolute drive path required')
    return '/mnt/'+path.drive[0].lower()+'/'+'/'.join(path.parts[1:])


def imported(destination,entry):
    # Existing index semantics: metadata/stat only; no new patient-file checksums.
    spec=importlib.util.spec_from_file_location('small_sync_existing_index',Path(__file__).with_name('write_local_index.py'))
    index=importlib.util.module_from_spec(spec);spec.loader.exec_module(index)
    if index.size_record_valid(destination,entry):return True
    states=destination/'verified-objects.json';obj=destination/'objects'/entry['Sha256']
    if states.is_file() and obj.is_file() and not obj.is_symlink() and obj.stat().st_size==entry['SizeBytes']:
        mapping=destination/'verified-files.json';record=load(mapping).get(entry['Sha256']+':'+entry['RelativePath']) if mapping.is_file() else None
        target=destination/record['LocalPath'] if record else None
        return bool(record and target.is_file() and not target.is_symlink() and os.path.samefile(obj,target)
            and load(states).get(entry['Sha256'])==dict(size=obj.stat().st_size,mtime_ns=obj.stat().st_mtime_ns))
    return False


def download_one(sftp,blocks,entry,directory,deadline,pause,*,reserve=15*1024**3,progress=None):
    require(now()<deadline,'Original sync deadline reached; no source access')
    require(not pause.exists(),'User pause marker present; no source access')
    final=directory/entry['DownloadFileName'];part=final.with_name(final.name+'.part');provenance=part.with_name(part.name+'.source.json')
    before=sftp.stat(entry['ServerPath'])
    require(before.st_size==entry['SizeBytes'],'Remote completed result size differs')
    identity=dict(RelativePath=entry['RelativePath'],RemotePath=entry['ServerPath'],ExpectedBytes=entry['SizeBytes'],
        ExpectedServerSha256=entry['Sha256'],RemoteMtime=before.st_mtime,SourceReceipt=entry['SourceReceipt'],Verification='size_only')
    completion=final.with_name(final.name+'.completed.json')
    if final.exists():
        require(not final.is_symlink() and final.stat().st_size==entry['SizeBytes'] and completion.is_file() and load(completion).get('SourceIdentity')==identity and load(completion).get('LocalMtimeNs')==final.stat().st_mtime_ns,'Unrelated/changed final file; refusing overwrite')
        return final
    require(not part.is_symlink() and not provenance.is_symlink(),'Unsafe partial result')
    if part.exists():require(provenance.is_file() and load(provenance)==identity,'Unrelated partial source identity; retained')
    else:
        require(not provenance.exists(),'Orphan partial metadata; retained for review')
        part.touch(exist_ok=False);atomic(provenance,identity)
    offset=part.stat().st_size;require(offset<=entry['SizeBytes'],'Partial result exceeds expected bytes')
    require(shutil.disk_usage(directory).free>=entry['SizeBytes']-offset+reserve,'Local result reserve required')
    last_progress=time.monotonic()
    if progress:progress(offset,entry['SizeBytes'])
    with sftp.open(entry['ServerPath'],'rb') as remote,part.open('ab') as out:
        while offset<entry['SizeBytes']:
            require(now()<deadline,'Original sync deadline reached; partial retained')
            require(not pause.exists(),'User pause marker present; partial retained')
            amount=min(8*1024**2,entry['SizeBytes']-offset);received=0
            require(shutil.disk_usage(directory).free>=reserve+amount,'Local disk reserve reached; partial retained')
            for block in blocks(remote,offset,amount):
                require(now()<deadline,'Original sync deadline reached; partial retained')
                require(not pause.exists(),'User pause marker present; partial retained')
                require(len(block)==min(32768,amount-received),'Incomplete contiguous SFTP chunk; valid prefix retained')
                out.write(block);received+=len(block)
                if progress and time.monotonic()-last_progress>=5:
                    out.flush();progress(offset+received,entry['SizeBytes']);last_progress=time.monotonic()
            require(received==amount,'Incomplete SFTP range; partial retained')
            out.flush();offset+=received
        os.fsync(out.fileno())
    require(now()<deadline,'Original sync deadline reached before final source verification; partial retained')
    after=sftp.stat(entry['ServerPath'])
    require((before.st_size,before.st_mtime)==(after.st_size,after.st_mtime) and part.stat().st_size==entry['SizeBytes'],'Remote result changed or incomplete')
    require(not final.exists(),'Final result appeared; partial retained')
    os.rename(part,final)  # Native Windows rename refuses replacement.
    atomic(completion,dict(Status='COMPLETED_SIZE_ONLY',SourceIdentity=identity,LocalMtimeNs=final.stat().st_mtime_ns,ActualSha256=None,CompletedUtc=now().isoformat()))
    return final


def import_one(profile,source,entry,snapshot,deadline):
    command=['wsl.exe','-d',profile['wsl_distribution'],'-u','root','--','python3',wsl_path(Path(__file__).with_name('import_verified_download.py').resolve()),
        '--source',wsl_path(source),'--plan',wsl_path(snapshot),'--relative-path',entry['RelativePath'],'--destination',wsl_path(Path(profile['destination']).resolve()),
        '--verification','size','--completed-download','--completion-evidence','ssh_completed_response','--stable-seconds','0','--execute']
    if profile['allow_onedrive']:command.append('--allow-onedrive')
    result=subprocess.run(command,text=True,capture_output=True,timeout=min(120,max(1,(deadline-now()).total_seconds())))
    if result.returncode and 'Stop the owned backup controller' in result.stderr:return False
    require(result.returncode==0,'Existing SIZE ONLY importer failed; completed download retained')
    return True


def run(profile,execute=False,watch=False):
    validate_profile(profile)
    plan=load(profile['plan']);files,summary=inventory(plan,profile)
    if not execute:return dict(status='DRY_RUN_NO_AUTH_CONFIG_KEY_NETWORK_OR_IMPORT',**summary)
    require(os.name=='nt','Execute with the existing native Windows Python/Paramiko environment')
    import paramiko
    require(paramiko.__version__=='5.1.0','Existing reviewed Paramiko5.1.0 required')
    deadline=min(DEADLINE,dt.datetime.fromisoformat(profile['deadline_utc']),now()+dt.timedelta(hours=profile['max_hours']))
    require(now()<deadline,'Original sync deadline already reached')
    destination=Path(profile['destination']).resolve();directory=Path(profile['download_directory']).resolve();state=Path(profile['state_directory']).resolve()
    for path in (destination,directory,state):path.mkdir(mode=0o700,parents=True,exist_ok=True)
    import msvcrt
    lock=(state/'small-results.lock').open('a+b');lock.write(b'0');lock.flush();lock.seek(0)
    try:msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
    except OSError:lock.close();raise ValueError('This small-result sync already owns its lock')
    def status(value,**extra):
        doc=dict(status=value,updated_utc=now().isoformat(),deadline_utc=deadline.isoformat(),patient_file_hashes_computed=False,existing_BAM_writer_changed=False,**extra)
        atomic(state/'status.json',doc);return doc
    try:
        # Read existing endpoint/pin/key-path configuration in memory only.
        # No private key, auth values or alternate credential config is written.
        auth=load(profile['auth_config'])
        require(directory!=Path(auth['download_directory']).resolve(),'Dedicated small-result download subdirectory required')
        from native_sftp_download import connect,blocks
        while now()<deadline:
            if Path(profile['pause_file']).exists():return status('PAUSED_USER_MARKER_PARTIALS_RETAINED')
            plan_raw=Path(profile['plan']).read_bytes();plan=json.loads(plan_raw);files,summary=inventory(plan,profile)
            pending=[entry for entry in files if not imported(destination,entry)]
            if not pending and not summary['missing_completed_roots']:return status('PASS_ALL_CURRENT_REQUIRED_COMPLETED_RESULTS_IMPORTED_SIZE_ONLY',**summary)
            if pending:
                snapshot=state/('source-plan-'+hashlib.sha256(plan_raw).hexdigest()+'.json')
                if not snapshot.exists():
                    snapshot.write_text(json.dumps(dict(SchemaVersion=1,Files=files,Receipts=plan['Receipts'],SourcePlanMetadataSha256=hashlib.sha256(plan_raw).hexdigest()),indent=2)+'\n')
                client,sftp=connect(auth)
                sftp.get_channel().settimeout(min(60,max(1,(deadline-now()).total_seconds())))
                try:
                    for number,entry in enumerate(pending,1):
                        require(now()<deadline,'Original sync deadline reached; completed files/partials retained')
                        status('SYNCING_COMPLETED_SMALL_RESULTS_SIZE_ONLY',selected_files=summary['selected_files'],pending_files=len(pending),file_index=number)
                        source=download_one(sftp,blocks,entry,directory,deadline,Path(profile['pause_file']),reserve=profile['min_free_bytes'],
                            progress=lambda received,expected:status('DOWNLOADING_SMALL_RESULT_SIZE_ONLY',file_index=number,pending_files=len(pending),received_bytes=received,expected_bytes=expected))
                        import_one(profile,source,entry,snapshot,deadline)
                finally:sftp.close();client.close()
            if not watch:return status('ONE_PASS_FINISHED_FUTURE_OUTPUTS_MAY_BE_PENDING',**summary)
            status('WAITING_COMPLETED_OUTPUTS_OR_IMPORT_LOCK_SIZE_ONLY',**summary)
            time.sleep(min(profile['interval_seconds'],max(0,(deadline-now()).total_seconds())))
        return status('STOPPED_ORIGINAL_DEADLINE_PARTIALS_RETAINED')
    except BaseException as error:
        status('FAILED_PARTIALS_RETAINED_NO_AUTOMATIC_RESTART',error_type=type(error).__name__,reason=str(error) if isinstance(error,SyncError) else 'Underlying transport/import error; auth configuration is not logged')
        raise
    finally:
        lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--profile',type=Path,required=True);parser.add_argument('--execute',action='store_true');parser.add_argument('--watch',action='store_true');args=parser.parse_args()
    try:print(json.dumps(run(load(args.profile),args.execute,args.watch),indent=2))
    except BaseException as error:
        # Endpoint/key paths/credentials from third-party exceptions never go to stdout.
        print('FAILED: '+(str(error) if isinstance(error,SyncError) else type(error).__name__)+'; completed files and resumable partials retained; existing BAM writer unchanged',file=sys.stderr);sys.exit(1)
