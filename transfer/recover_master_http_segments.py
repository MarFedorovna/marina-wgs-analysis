#!/usr/bin/env python3
"""Reviewed local HTTP ZIP recovery of only originally missing master chunks."""
import argparse,datetime as dt,hashlib,importlib.util,json,os
from pathlib import Path
import re,shutil,sys,time,zipfile
END=dt.datetime.fromisoformat('2026-11-06T07:10:00+00:00')
CHUNK=64*1024**2
SUB=32*1024**2
ARTIFACT='masters/cohort.master.tsv.gz'

def need(ok,msg):
    if not ok:raise ValueError(msg)
def now():return dt.datetime.now(dt.timezone.utc)
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pinned(p):
    f=Path(p['path']);need(f.is_file() and not f.is_symlink() and f.stat().st_size==p['bytes'] and f.stat().st_size<16*1024**2 and sha(f)==p['sha256'],'Small input metadata pin changed');return f

def modules():
    f=Path(__file__).with_name('download_cohort_views_dataset_v2.py');need(sha(f)=='d6ad2bc5349e609e71ef67420284ade6dc6c25b78e4290e83e18df134c27b800','Frozen physical allocation helper changed')
    spec=importlib.util.spec_from_file_location('http_physical',f);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m.base_module(),m

def dead(pid):
    import ctypes
    k=ctypes.WinDLL('kernel32',use_last_error=True);op=k.OpenProcess;op.argtypes=[ctypes.c_uint32,ctypes.c_int,ctypes.c_uint32];op.restype=ctypes.c_void_p
    handle=op(0x1000,False,pid)
    if not handle:need(ctypes.get_last_error()==87,'Cannot verify original worker exit');return True
    try:
        code=ctypes.c_uint32();need(k.GetExitCodeProcess(ctypes.c_void_p(handle),ctypes.byref(code)),'Cannot inspect original worker');return code.value!=259
    finally:k.CloseHandle(ctypes.c_void_p(handle))

def spec_for(snapshot,index,suboffset):
    return dict(schema_version=1,status='CLOSED_MISSING_MASTER_SEGMENT',source=snapshot['source'],expected=snapshot['expected'],file_id=snapshot['file_id'],index=index,suboffset=suboffset,offset=index*CHUNK+suboffset,bytes=min(SUB,snapshot['expected']['bytes']-(index*CHUNK+suboffset)),initial_missing_indices=snapshot['missing'])

def segment_contract(value,snapshot):
    need(set(value)=={'schema_version','status','source','expected','file_id','index','suboffset','offset','bytes','initial_missing_indices','payload_crc32'},'Exact segment metadata schema required')
    need(type(value['index']) is int and value['index'] in snapshot['missing'] and value['suboffset'] in (0,SUB),'Segment not in original missing chunks')
    expected=spec_for(snapshot,value['index'],value['suboffset']);need(expected['bytes']>0 and all(value[k]==v for k,v in expected.items()),'Segment source/offset/length differs from approved snapshot');need(re.fullmatch('[0-9a-f]{8}',value['payload_crc32']),'Payload CRC declaration required');return expected

def verify_local(ctx):
    b=ctx['base'];snap=ctx['snapshot'];p=ctx['part'];current=load(ctx['meta'])
    need(now()<ctx['deadline'],'New HTTP recovery deadline reached')
    need(p.is_file() and not p.is_symlink() and p.stat().st_ino==snap['file_id'] and p.stat().st_size==snap['expected']['bytes'],'Unowned/changed master part refused')
    need(ctx['physical'](p)>=snap['expected']['bytes'],'Master is not fully physically allocated')
    for key in ['schema_version','artifact','expected','source','chunk_bytes','file_id']:need(current[key]==snap[key],'Original source ledger pins changed')
    initial=set(snap['completed']);done=set(current['completed']);count=(snap['expected']['bytes']+CHUNK-1)//CHUNK
    need(current['completed']==sorted(done) and initial<=done<=set(range(count)),'Existing64MiB bitmap invalid')
    progress=ctx['progress'];ready=set()
    for index in snap['missing']:
        subs=[s for s in (0,SUB) if index*CHUNK+s<snap['expected']['bytes']]
        if all(f'{index}:{s}' in progress['segments'] for s in subs):ready.add(index)
    need(done-initial<=ready,'Bitmap claims completion without own durable HTTP segment receipts')
    return current,ready

def apply_segment(ctx,path,declaration):
    current,ready=verify_local(ctx);before=path.stat();need(before.st_size==declaration['bytes'],'Closed segment ZIP size differs')
    with zipfile.ZipFile(path) as z:
        need(len(z.infolist())==2 and set(z.namelist())=={'block.bin','segment.json'},'ZIP must contain exactly block.bin and segment.json')
        need(all(x.compress_type==zipfile.ZIP_STORED and not x.is_dir() and not x.flag_bits&1 for x in z.infolist()),'Only unencrypted stored ZIP segments accepted')
        need(z.getinfo('segment.json').file_size<=16384,'Segment metadata cap');value=json.loads(z.read('segment.json'));expected=segment_contract(value,ctx['snapshot']);info=z.getinfo('block.bin');need(info.file_size==expected['bytes'] and f'{info.CRC:08x}'==value['payload_crc32'],'ZIP payload size/CRC differs')
        need(all(value[k]==declaration[k] for k in ['index','suboffset','offset','payload_crc32']) and value['bytes']==declaration['payload_bytes'],'Segment differs from closed HTTP manifest')
        need(z.testzip() is None,'ZIP CRC failed before any master write');key=f"{value['index']}:{value['suboffset']}"
        if key in ctx['progress']['segments']:return False
        need(value['index'] not in current['completed'],'Never rewrite an existing completed64MiB chunk')
        need(shutil.disk_usage(ctx['part'].parent).free>=5*1024**3,'Runtime5GiB reserve reached')
        with z.open('block.bin') as inp,ctx['part'].open('r+b',buffering=0) as out:
            out.seek(value['offset']);written=0
            while written<value['bytes']:
                need(now()<ctx['deadline'],'HTTP recovery deadline; uncommitted segment retained');data=inp.read(min(8*1024**2,value['bytes']-written));need(data,'Truncated ZIP block');ctx['base'].write_all(out,data);written+=len(data)
            need(not inp.read(1),'Unexpected extra payload');os.fsync(out.fileno())
    need(path.stat().st_size==before.st_size and path.stat().st_mtime_ns==before.st_mtime_ns,'ZIP changed during write; segment not committed')
    ctx['progress']['segments'][key]=dict(index=value['index'],suboffset=value['suboffset'],offset=value['offset'],bytes=value['bytes'],payload_crc32=value['payload_crc32'],ZIP_CRC_verified=True,completed_utc=now().isoformat());ctx['base'].atomic(ctx['progress_path'],ctx['progress'])
    current,ready=verify_local(ctx);current['completed']=sorted(set(current['completed'])|ready);ctx['base'].atomic(ctx['meta'],current);return True

def validate_other_files(base,plan,files,skip_master=True):
    destination=Path(plan['destination']);state=Path(plan['state_directory']);count=0
    for name,expected in files.items():
        if name==ARTIFACT and skip_master:continue
        target,part,meta,note=base.local_paths(destination,state,name);need(base.completed(target,note,expected),'Another closed dataset result or receipt is missing/changed');n=load(note);need(n.get('status') in ['COMPLETED_DATASET_FILE_SIZE_ONLY','ROOT_REVIEWED_SOURCE_ALIAS_SIZE_ONLY'],'Unexpected existing completion status');count+=1
    return count

def run(profile_path,profile_sha,execute=False,watch=False):
    need(os.name=='nt','Existing Windows local receiver required');need(sha(profile_path)==profile_sha,'Exact ROOT-reviewed HTTP profile SHA required');p=load(profile_path)
    need(p['schema_version']==1 and p['status']=='ROOT_REVIEWED_MASTER_HTTP_RECOVERY_V1' and p['producer_sha256']==sha(__file__),'Reviewed new HTTP recovery profile required');end=dt.datetime.fromisoformat(p['deadline_utc']);need(now()<end<=END,'New explicit recovery deadline invalid');need(2<=p['interval_seconds']<=30,'Bounded local watch interval required')
    base,physical=modules();plan=load(pinned(p['original_plan']));snapshot=load(pinned(p['captured_ledger']));http=load(pinned(p['http_manifest']));_,files=base.validate(plan)
    need(snapshot['artifact']==ARTIFACT and snapshot['expected']==files[ARTIFACT] and snapshot['chunk_bytes']==CHUNK,'Captured master authority differs')
    all_chunks=set(range((snapshot['expected']['bytes']+CHUNK-1)//CHUNK));need(snapshot['missing']==sorted(all_chunks-set(snapshot['completed'])) and not(set(snapshot['completed'])&set(snapshot['missing'])),'Original missing set differs')
    need(http.get('source_size_mtime_unchanged') is True and http.get('source_payload_rehashed') is False and http.get('zip_crc_verified') is True and http.get('source_path')==base.REMOTE+ARTIFACT and http['schema_version']==1 and http['status']=='PASS_CLOSED_MISSING_MASTER_HTTP_SEGMENTS' and all(http[k]==snapshot[k] for k in ['source','expected','file_id']) and http['initial_missing_indices']==snapshot['missing'],'Closed HTTP source authority differs')
    segment_map={};pairs=set()
    for s in http['segments']:
        need(re.fullmatch(r'chunk-[0-9]{4}-sub-0[01]\.zip',s['name']) and s['name'] not in segment_map,'Unsafe/duplicate closed segment filename')
        spec=spec_for(snapshot,s['index'],s['suboffset']);need(s['name']==f"chunk-{s['index']:04d}-sub-{s['suboffset']//SUB:02d}.zip",'Segment filename/index mismatch');need(s['index'] in snapshot['missing'] and s['suboffset'] in (0,SUB) and spec['bytes']>0 and s['offset']==spec['offset'] and s['payload_bytes']==spec['bytes'] and re.fullmatch('[0-9a-f]{8}',s['payload_crc32']) and type(s['bytes']) is int and 0<s['bytes']<=SUB+65536,'HTTP segment declaration invalid');segment_map[s['name']]=s;pairs.add((s['index'],s['suboffset']))
    expected_pairs={(i,s) for i in snapshot['missing'] for s in (0,SUB) if i*CHUNK+s<snapshot['expected']['bytes']};need(pairs==expected_pairs and len(pairs)==len(segment_map),'HTTP segment set not exactly all original missing bytes')
    need(isinstance(p['old_pids'],list) and 1<=len(p['old_pids'])<=4 and len(set(p['old_pids']))==len(p['old_pids']) and all(type(pid) is int and pid>0 and dead(pid) for pid in p['old_pids']),'Old worker PID is alive/reused; no concurrent writer allowed')
    destination=Path(plan['destination']);original_state=Path(plan['state_directory']);target,part,meta,note=base.local_paths(destination,original_state,ARTIFACT);need(str(meta)==snapshot['ledger_path'],'Original ledger path differs')
    state=Path(p['state_directory']).resolve();downloads=Path(p['download_directory']).resolve();need(downloads.is_dir() and not downloads.is_symlink() and state!=original_state.resolve(),'Dedicated HTTP state/download directory required');identity=dict(profile_sha256=profile_sha,ledger_sha256=p['captured_ledger']['sha256'],http_manifest_sha256=p['http_manifest']['sha256']);progress_path=state/'segment-progress.json'
    progress=load(progress_path) if progress_path.exists() else dict(schema_version=1,status='HTTP_MISSING_SEGMENTS_PROGRESS',identity=identity,segments={});need(progress['identity']==identity,'Existing HTTP recovery source/profile changed')
    ctx=dict(base=base,snapshot=snapshot,part=part,meta=meta,physical=physical.physical_allocation,deadline=end,progress=progress,progress_path=progress_path)
    if target.exists():need(base.completed(target,note,snapshot['expected']),'Unrelated final master; no overwrite');return dict(status='ALREADY_FINAL_SIZE_ONLY')
    verify_local(ctx);validate_other_files(base,plan,files)
    found=sum((downloads/name).is_file() and (downloads/name).stat().st_size==s['bytes'] for name,s in segment_map.items())
    if not execute:return dict(status='DRY_HTTP_RECOVERY_NO_WRITES_OR_NETWORK',segments=len(segment_map),available_segment_ZIPs=found,initial_missing_chunks=len(snapshot['missing']),unchanged_original_plan=True,new_deadline_utc=end.isoformat())
    import msvcrt
    own=(original_state/'dataset.lock').open('r+b');own.seek(0);msvcrt.locking(own.fileno(),msvcrt.LK_NBLCK,1)
    state.mkdir(parents=True,exist_ok=True)
    for name,source in [('original-failure-status.json',original_state/'status.json'),('original-master-ledger.json',meta)]:
        f=state/name
        if not f.exists():f.write_bytes(source.read_bytes())
    def status(s,**extra):base.atomic(state/'status.json',dict(status=s,new_deadline_utc=end.isoformat(),profile_sha256=profile_sha,completed32MiBsegments=len(progress['segments']),total_segments=len(segment_map),updated_utc=now().isoformat(),patient_file_hashes_computed=False,**extra))
    try:
        while now()<end:
            for name,declaration in segment_map.items():
                path=downloads/name;key=f"{declaration['index']}:{declaration['suboffset']}"
                if key in progress['segments']:continue
                if path.is_file() and path.stat().st_size==declaration['bytes']:
                    need(not path.is_symlink(),'Segment symlink refused');apply_segment(ctx,path,declaration);status('RECOVERING_ONLY_ORIGINAL_MISSING_MASTER_BYTES')
            current,ready=verify_local(ctx)
            if current['completed']==sorted(all_chunks):
                validate_other_files(base,plan,files);need(not target.exists() and part.stat().st_ino==snapshot['file_id'],'Final nooverwrite/identity gate failed');os.rename(part,target);base.atomic(note,dict(status='COMPLETED_DATASET_FILE_SIZE_ONLY',expected=snapshot['expected'],source=snapshot['source'],mtime_ns=target.stat().st_mtime_ns,actual_sha256=None,http_recovery_profile=identity));validate_other_files(base,plan,files,False)
                final=dict(status='COMPLETED_CLOSED_DATASET_SIZE_ONLY_HTTP_RECOVERY',whole_dataset_complete=True,selected_files=len(files),selected_bytes=sum(v['bytes'] for v in files.values()),committed_selected_bytes=sum(v['bytes'] for v in files.values()),updated_utc=now().isoformat(),patient_hashes_computed=False,original_plan_unchanged=True,new_recovery_deadline_utc=end.isoformat(),recovery_producer_sha256=sha(__file__),recovery_profile=identity);base.atomic(state/'LOCAL-HTTP-RECOVERY-RECEIPT.json',final);base.atomic(original_state/'status.json',final);status(final['status'],whole_dataset_complete=True);return final
            status('WAITING_REVIEWED_HTTP_SEGMENT_ZIPS')
            if not watch:return dict(status='ONE_PASS_PARTIAL_SEGMENTS_PRESERVED',completed_segments=len(progress['segments']))
            time.sleep(min(p['interval_seconds'],max(0,(end-now()).total_seconds())))
        status('STOPPED_NEW_EXPLICIT_DEADLINE_PARTIALS_RETAINED');return dict(status='STOPPED_NEW_EXPLICIT_DEADLINE_PARTIALS_RETAINED')
    except BaseException as e:status('FAILED_OWN_PROGRESS_RETAINED_NO_RESTART',error_type=type(e).__name__);raise
    finally:own.seek(0);msvcrt.locking(own.fileno(),msvcrt.LK_UNLCK,1);own.close()

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--profile',type=Path,required=True);a.add_argument('--profile-sha256',required=True);a.add_argument('--execute',action='store_true');a.add_argument('--watch',action='store_true');v=a.parse_args()
    try:print(json.dumps(run(v.profile,v.profile_sha256,v.execute,v.watch)))
    except BaseException as e:print(json.dumps(dict(status='FAILED_HTTP_RECOVERY_INPUTS_PRESERVED',error_type=type(e).__name__)),file=sys.stderr);sys.exit(1)
