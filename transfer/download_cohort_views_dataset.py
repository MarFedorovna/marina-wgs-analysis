#!/usr/bin/env python3
"""Closed-receipt dataset transfer: sixteen clients, in-place resumable chunks."""
import argparse
import base64
from concurrent.futures import ThreadPoolExecutor,as_completed
import ctypes
import datetime as dt
import hashlib
import json
import os
from pathlib import Path,PurePosixPath
import queue
import re
import shutil
import sys
import threading
import time
import unicodedata

END=dt.datetime.fromisoformat('2026-11-06T06:45:00+00:00')
REMOTE='/srv/marina/outputs/cohort-sample-views-v1/'
RESERVE=15*1024**3
CHUNK=64*1024**2
WINDOW=8*1024**2


def need(ok,message):
    if not ok:raise ValueError(message)


def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def now():return dt.datetime.now(dt.timezone.utc)
def key(name):return hashlib.sha256(name.encode()).hexdigest()[:32]


def atomic(path,value):
    path=Path(path);temporary=path.with_name(path.name+'.tmp')
    need(not path.is_symlink() and not temporary.is_symlink(),'Metadata symlink refused')
    with temporary.open('w',encoding='utf-8') as stream:json.dump(value,stream,ensure_ascii=False,indent=2);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)


def capture(item,name):
    need(set(item)=={'path','bytes','sha256','content_base64'} and item['path']==REMOTE+name,'Exact captured dataset metadata required')
    raw=base64.b64decode(item['content_base64'],validate=True)
    need(len(raw)<=32*1024**2 and len(raw)==item['bytes'] and hashlib.sha256(raw).hexdigest()==item['sha256'],'Captured metadata identity differs')
    return json.loads(raw)


def safe_relative(name):
    p=PurePosixPath(name)
    need(p.parts and not p.is_absolute() and '..' not in p.parts and '\\' not in name and p.as_posix()==name,'Unsafe artifact relative path')
    need(not any(c in ':<>"|?*' or unicodedata.category(c).startswith('C') for c in name),'Windows drive/invalid/control characters refused')
    need(all(x==x.rstrip(' .') and re.fullmatch(r'(CON|PRN|AUX|NUL|COM[1-9¹²³]|LPT[1-9¹²³])',x.split('.')[0].upper()) is None for x in p.parts),'Windows reserved/ambiguous name refused')
    need(not any(x.lower() in {'resources','databases','references','cache','inputs','raw'} for x in p.parts) and '.g.vcf' not in p.name.lower(),'Database/raw/gVCF artifact refused')
    need(p.name.lower().endswith(('.tsv','.csv','.tsv.gz','.csv.gz','.json','.html','.txt','.md','.log','.vcf.gz','.vcf.gz.tbi','.png','.jpg','.jpeg','.svg','.pdf')),'Unrecognized result artifact type')
    need(not p.name.lower().endswith(('.bam','.bai')),'BAM duplicates excluded')
    return p


def validate(plan):
    need(set(plan)=={'schema_version','status','deadline_utc','destination','state_directory','pause_file','source_provenance','source_aliases'} and plan['schema_version']==1 and plan['status']=='ROOT_REVIEWED_COHORT_VIEWS_DATASET_DOWNLOAD','Exact reviewed dataset plan required')
    deadline=dt.datetime.fromisoformat(plan['deadline_utc']);need(deadline.tzinfo is not None and deadline<=END,'Transfer deadline cannot exceed06:45UTC')
    provenance=plan['source_provenance'];need(set(provenance)=={'manifest','receipt'},'Exactly two captured producer records required')
    m=capture(provenance['manifest'],'manifest.json');r=capture(provenance['receipt'],'workflow-receipt.json')
    need(m.get('schema_version')==1 and m.get('status')=='PASS' and m.get('synthetic') is False and m.get('report_input_kind')=='FULLY_ANNOTATED_VCF_PACKETS_V1' and m.get('report_producer_kind')=='FAST_DIRECT_PACKET_VIEWS_V2','Final VCF-first real dataset PASS required')
    need(all(m.get(f) is True for f in ('sample_packet_scope_verified','public_header_bundles_verified','no_external_requests','no_legacy_receipt_created')) and m.get('primary_workflow_changed') is False and m.get('automatic_ACMG') is False,'Dataset preservation/privacy contract required')
    need(r.get('schema_version')==1 and r.get('status')=='PASS' and r.get('stage')=='cohort-sample-views-v1' and r.get('synthetic') is False and re.fullmatch('[0-9a-f]{64}',r.get('identity','')) is not None,'Exact dataset completion receipt required')
    expected=dict(m['artifacts'],**{'manifest.json':{k:provenance['manifest'][k] for k in ('bytes','sha256')}})
    need(r.get('artifacts')==expected and r.get('producer_code_sha256')==m.get('producer_code_sha256'),'Closed manifest/receipt output set differs')
    files=dict(expected,**{'workflow-receipt.json':{k:provenance['receipt'][k] for k in ('bytes','sha256')}})
    need(0<len(files)<=2048,'Bounded closed artifact list required')
    for name,item in files.items():
        safe_relative(name);need(set(item)=={'bytes','sha256'} and type(item['bytes']) is int and 0<=item['bytes']<=64*1024**3 and re.fullmatch('[0-9a-f]{64}',item['sha256']) is not None,'Invalid expected artifact identity')
    need(sum(v['bytes'] for v in files.values())<=128*1024**3,'Bounded dataset byte ceiling exceeded')
    aliases=plan['source_aliases'];need(isinstance(aliases,list) and len(aliases)<=2 and len({a['artifact'] for a in aliases})==len(aliases),'Only reviewed main VCF/TBI aliases permitted')
    for a in aliases:
        need(set(a)=={'artifact','local_path','completion_receipt'} and a['artifact'] in ('vcf/cohort.fully-annotated.vcf.gz','vcf/cohort.fully-annotated.vcf.gz.tbi') and a['artifact'] in files,'Unrecognized source alias')
    return deadline,files


def private_path(path):
    p=Path(path);need(p.is_absolute() and not any(x.is_symlink() or (x/'.git').exists() for x in (p,*p.parents)),'Private absolute path outside Git/symlinks required');return p


def allocation_bytes(path):
    if os.name!='nt':return path.stat().st_blocks*512
    kernel=ctypes.WinDLL('kernel32',use_last_error=True);high=ctypes.c_ulong();kernel.GetCompressedFileSizeW.restype=ctypes.c_ulong
    ctypes.set_last_error(0)
    low=kernel.GetCompressedFileSizeW(str(path),ctypes.byref(high))
    need(low!=0xffffffff or ctypes.get_last_error()==0,'Cannot inspect existing allocated bytes')
    return (high.value<<32)|low


def preallocate(path,size):
    need(not path.exists() and not path.is_symlink(),'Refusing unrelated owned-part overwrite')
    with path.open('x+b') as stream:
        if size and os.name=='nt':
            import msvcrt
            kernel=ctypes.WinDLL('kernel32',use_last_error=True);allocation=ctypes.c_longlong(size)
            need(kernel.SetFileInformationByHandle(ctypes.c_void_p(msvcrt.get_osfhandle(stream.fileno())),5,ctypes.byref(allocation),ctypes.sizeof(allocation)),'Physical in-place allocation failed')
        elif size and hasattr(os,'posix_fallocate'):os.posix_fallocate(stream.fileno(),0,size)
        stream.truncate(size);stream.flush();os.fsync(stream.fileno())


def local_paths(destination,state,name):
    safe_relative(name);target=destination/Path(name);private_path(target)
    need(target.resolve(strict=False).is_relative_to(destination.resolve()),'Artifact escapes local result root')
    return target,target.with_name(target.name+'.part'),state/'files'/(key(name)+'.chunks.json'),state/'files'/(key(name)+'.completed.json')


def completed(target,note,expected):
    if not target.exists():return False
    need(target.is_file() and not target.is_symlink() and note.is_file(),'Unrelated existing final result; no overwrite')
    value=load(note);need(value['expected']==expected and target.stat().st_size==expected['bytes'] and target.stat().st_mtime_ns==value['mtime_ns'],'Existing completed result changed')
    return True


def alias_file(alias,expected,target,note):
    source=private_path(alias['local_path']);record=load(alias['completion_receipt']);identity=record.get('SourceIdentity',{})
    need(record.get('Status')=='COMPLETED_SIZE_ONLY' and identity.get('ExpectedBytes')==expected['bytes'] and identity.get('ExpectedServerSha256')==expected['sha256'] and source.is_file() and source.stat().st_size==expected['bytes'] and source.stat().st_mtime_ns==record.get('LocalMtimeNs'),'Reviewed source alias identity/size/mtime differs')
    target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists():need(os.path.samefile(target,source),'Existing alias destination differs; refusing overwrite')
    else:os.link(source,target)
    atomic(note,dict(status='ROOT_REVIEWED_SOURCE_ALIAS_SIZE_ONLY',expected=expected,mtime_ns=target.stat().st_mtime_ns,actual_sha256=None))


def required_additional(files,destination,state,aliases):
    total=0
    for name,expected in files.items():
        target,part,meta,note=local_paths(destination,state,name)
        if name in aliases or completed(target,note,expected):continue
        if part.exists():
            need(meta.is_file() and not part.is_symlink(),'Unbound existing part retained')
            value=load(meta);need(value['expected']==expected and value['artifact']==name and value['file_id']==part.stat().st_ino and part.stat().st_size==expected['bytes'],'Owned part identity differs')
            total+=max(0,expected['bytes']-allocation_bytes(part))
        else:total+=expected['bytes']
    return total


def write_all(output,data):
    view=memoryview(data);offset=0
    while offset<len(view):
        written=output.write(view[offset:]);need(type(written) is int and 0<written<=len(view)-offset,'Short/zero local write; chunk not committed');offset+=written


def transfer_chunk(sftp,blocks,task,stop,deadline,pause,progress):
    name,expected,target,part,meta,note,source,lock,index,offset,length=task
    need(index>=0 and offset==index*CHUNK and 0<length==min(CHUNK,expected['bytes']-offset),'Chunk bounds differ from fixed source partition')
    need(not stop.is_set() and now()<deadline and not pause.exists(),'Deadline/pause/peer stop; chunks retained')
    stat=sftp.stat(REMOTE+name);need((stat.st_size,stat.st_mtime)==(source['bytes'],source['remote_mtime']),'Remote immutable result changed')
    received=0
    with sftp.open(REMOTE+name,'rb') as remote,part.open('r+b',buffering=0) as output:
        output.seek(offset)
        while received<length:
            need(not stop.is_set() and now()<deadline and not pause.exists(),'Chunk stop; uncommitted chunk will resume from boundary')
            need(shutil.disk_usage(part.parent).free>=RESERVE,'15GiB runtime reserve reached')
            amount=min(WINDOW,length-received);got=0
            for data in blocks(remote,offset+received,amount):
                need(not stop.is_set() and now()<deadline,'Deadline/peer stop')
                need(len(data)==min(32768,amount-got),'Short contiguous SFTP chunk');write_all(output,data);got+=len(data)
            need(got==amount,'Incomplete SFTP chunk');received+=got
        os.fsync(output.fileno())
    stat=sftp.stat(REMOTE+name);need((stat.st_size,stat.st_mtime)==(source['bytes'],source['remote_mtime']),'Source changed after chunk')
    with lock:
        value=load(meta);need(value['source']==source and value['expected']==expected,'Chunk source ledger changed')
        need(index not in value['completed'],'Duplicate completed chunk write refused');value['completed'].append(index);value['completed'].sort();atomic(meta,value)
        count=(expected['bytes']+CHUNK-1)//CHUNK
        if value['completed']==list(range(count)):
            final=sftp.stat(REMOTE+name)
            need((final.st_size,final.st_mtime)==(source['bytes'],source['remote_mtime']) and part.stat().st_size==expected['bytes'] and not target.exists() and now()<deadline,'Final file immutable/nooverwrite/deadline gate failed')
            os.rename(part,target)
            atomic(note,dict(status='COMPLETED_DATASET_FILE_SIZE_ONLY',expected=expected,source=source,mtime_ns=target.stat().st_mtime_ns,actual_sha256=None))
    progress(length)


def run(plan,auth=None,execute=False,small_only=False):
    deadline,all_files=validate(plan);files={n:v for n,v in all_files.items() if not small_only or v['bytes']<=20*1024**2};aliases={a['artifact']:a for a in plan['source_aliases']}
    if not execute:return dict(status='DRY_NO_AUTH_KEY_NETWORK_OR_FILES',selected_files=len(files),selected_bytes=sum(v['bytes'] for v in files.values()),pending_large_files=len(all_files)-len(files),whole_dataset_complete=False)
    need(os.name=='nt' and now()<deadline,'Existing Windows runtime and current transfer window required')
    destination=private_path(plan['destination']);state=private_path(plan['state_directory']);pause=Path(plan['pause_file'])
    need(destination.parent.is_dir() and state.parent.is_dir(),'Existing user-selected parent folders required')
    destination.mkdir(exist_ok=True);state.mkdir(exist_ok=True);(state/'files').mkdir(exist_ok=True)
    import msvcrt
    own=(state/'dataset.lock').open('a+b');own.write(b'0');own.flush();own.seek(0);msvcrt.locking(own.fileno(),msvcrt.LK_NBLCK,1)
    stop=threading.Event();clients=[];mutex=threading.Lock();tasks=queue.Queue();locks={};done_bytes=[0];last=[0.];pool=None
    def status(label,force=False,**extra):
        with mutex:
            if not force and time.monotonic()-last[0]<2:return
            last[0]=time.monotonic();value=dict(status=label,updated_utc=now().isoformat(),deadline_utc=deadline.isoformat(),committed_selected_bytes=done_bytes[0],selected_bytes=sum(v['bytes'] for v in files.values()),selected_files=len(files),pending_large_files=len(all_files)-len(files),whole_dataset_complete=False,patient_hashes_computed=False);value.update(extra);atomic(state/'status.json',value)
    try:
        additional=required_additional(files,destination,state,aliases);need(shutil.disk_usage(destination).free>=RESERVE+additional,'Whole selected dataset does not fit above15GiB; no files omitted')
        for name,a in aliases.items():
            if name not in files:continue
            target,part,meta,note=local_paths(destination,state,name);alias_file(a,files[name],target,note)
        import paramiko
        need(paramiko.__version__=='5.1.0','Existing Paramiko5.1.0 required');from native_sftp_download import connect,blocks
        config=load(auth)
        def connection(slot):
            stop.wait(.2*slot)
            for attempt in range(3):
                need(not stop.is_set() and now()+dt.timedelta(seconds=50)<deadline,'Connection time/peer stop')
                try:
                    client,sftp=connect(config)
                    with mutex:clients.append((client,sftp))
                    sftp.get_channel().settimeout(15);return client,sftp
                except Exception:
                    if attempt==2:raise
                    stop.wait(.5*(attempt+1))
        control,control_sftp=connection(0);pending=[]
        for number,(name,expected) in enumerate(sorted(files.items(),key=lambda pair:(not pair[0].endswith('top20.tsv'),pair[1]['bytes'],pair[0])),1):
            target,part,meta,note=local_paths(destination,state,name)
            if completed(target,note,expected):done_bytes[0]+=expected['bytes'];continue
            target.parent.mkdir(parents=True,exist_ok=True);remote=control_sftp.stat(REMOTE+name);need(remote.st_size==expected['bytes'],'Captured source size changed')
            source=dict(artifact=name,bytes=expected['bytes'],sha256=expected['sha256'],remote_mtime=remote.st_mtime,manifest_sha256=plan['source_provenance']['manifest']['sha256'])
            if meta.exists():
                value=load(meta);need(value['source']==source and value['expected']==expected and value['chunk_bytes']==CHUNK and value['file_id']==part.stat().st_ino and part.stat().st_size==expected['bytes'],'Existing owned chunk ledger/source differs')
            else:
                need(not part.exists(),'Orphan part retained without overwrite');status('PREPARING_PHYSICALLY_PREALLOCATED_RESULT',force=True,artifact=name,file_index=number);preallocate(part,expected['bytes']);value=dict(schema_version=1,artifact=name,expected=expected,source=source,chunk_bytes=CHUNK,file_id=part.stat().st_ino,completed=[]);atomic(meta,value)
            count=(expected['bytes']+CHUNK-1)//CHUNK;need(len(value['completed'])==len(set(value['completed'])) and all(type(i) is int and 0<=i<count for i in value['completed']),'Invalid durable chunk bitmap')
            lock=locks.setdefault(name,threading.Lock());pending.append((name,expected,target,part,meta,note,source,lock,count))
            for i in range(count):
                length=min(CHUNK,expected['bytes']-i*CHUNK)
                if i in value['completed']:done_bytes[0]+=length
                else:tasks.put((name,expected,target,part,meta,note,source,lock,i,i*CHUNK,length))
        def progress(length):
            with mutex:done_bytes[0]+=length
            status('DOWNLOADING_CLOSED_DATASET_CHUNKS_SIZE_ONLY')
        def worker(slot,existing=None):
            client=sftp=None
            try:
                client,sftp=existing or connection(slot)
                while not stop.is_set():
                    try:task=tasks.get_nowait()
                    except queue.Empty:return
                    transfer_chunk(sftp,blocks,task,stop,deadline,pause,progress)
            finally:
                if existing is None and client is not None:
                    try:sftp.close();client.close()
                    except Exception:pass
        status('DOWNLOADING_CLOSED_DATASET_CHUNKS_SIZE_ONLY',force=True)
        pool=ThreadPoolExecutor(max_workers=16);futures=[pool.submit(worker,i,(control,control_sftp) if i==0 else None) for i in range(16)]
        for future in as_completed(futures,timeout=max(1,(deadline-now()).total_seconds())):future.result()
        for name,expected,target,part,meta,note,source,lock,count in pending:
            if completed(target,note,expected):continue
            value=load(meta);need(value['completed']==list(range(count)) and part.stat().st_size==expected['bytes'],'Not every required chunk is committed')
            remote=control_sftp.stat(REMOTE+name);need((remote.st_size,remote.st_mtime)==(source['bytes'],source['remote_mtime']) and not target.exists() and now()<deadline,'Final immutable source/time/nooverwrite gate failed')
            os.rename(part,target);atomic(note,dict(status='COMPLETED_DATASET_FILE_SIZE_ONLY',expected=expected,source=source,mtime_ns=target.stat().st_mtime_ns,actual_sha256=None))
        status('SMALL_PHASE_COMPLETE_LARGE_FILES_PENDING' if small_only and len(files)<len(all_files) else 'PASS_ALL_CLOSED_DATASET_FILES_SIZE_ONLY',force=True,whole_dataset_complete=not small_only or len(files)==len(all_files))
        return dict(status='SMALL_PHASE_COMPLETE_LARGE_FILES_PENDING' if small_only and len(files)<len(all_files) else 'PASS_ALL_CLOSED_DATASET_FILES_SIZE_ONLY',files=len(files),whole_dataset_complete=not small_only or len(files)==len(all_files),patient_hashes_computed=False)
    except BaseException as error:status('FAILED_OWN_PARTS_AND_BITMAPS_RETAINED',force=True,error_type=type(error).__name__);raise
    finally:
        stop.set()
        for client,sftp in clients:
            try:sftp.close();client.close()
            except Exception:pass
        if pool is not None:pool.shutdown(wait=True,cancel_futures=True)
        own.seek(0);msvcrt.locking(own.fileno(),msvcrt.LK_UNLCK,1);own.close()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--auth-config',type=Path);p.add_argument('--small-only',action='store_true');p.add_argument('--execute',action='store_true');a=p.parse_args()
    try:print(json.dumps(run(load(a.plan),a.auth_config,a.execute,a.small_only)))
    except BaseException as error:print('FAILED: '+type(error).__name__+'; owned parts/bitmaps retained; no private connection values logged',file=sys.stderr);sys.exit(1)
