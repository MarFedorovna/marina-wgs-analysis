#!/usr/bin/env python3
"""Sixteen bounded receipt-bound SFTP ranges; preserve the existing valid native prefix."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import ctypes
import datetime as dt
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import threading
import time
import download_final_annotated_vcf_v2 as base


def need(ok,message):
    if not ok:raise ValueError(message)


def partition(prefix,total,connections=16):
    need(0<=prefix<total,'A strict incomplete prefix is required')
    need(1<=connections<=16,'At most sixteen connections permitted')
    remaining=total-prefix;result=[];offset=prefix
    for i in range(connections):
        length=remaining//connections+int(i<remaining%connections)
        need(length>0,'Enough remaining bytes for all ranges required')
        result.append(dict(index=i,offset=offset,length=length));offset+=length
    need(offset==total,'Range partition must be contiguous and complete')
    return result


def check_space(directory,remaining,free=None):
    available=shutil.disk_usage(directory).free if free is None else free
    need(available>=base.RESERVE+2*remaining,'15GiB reserve plus two copies of remaining bytes required')


def transfer_range(sftp,blocks,entry,identity,path,helper,deadline,stop,progress,pause=None):
    metadata=path.with_name(path.name+'.source.json');complete=path.with_name(path.name+'.completed.json')
    need(not path.is_symlink() and not metadata.is_symlink() and not complete.is_symlink(),'Range symlink refused')
    if metadata.exists():need(helper.load(metadata)==identity,'Range source/offset/length/mtime identity changed')
    else:
        need(not path.exists(),'Unbound range retained; refusing overwrite');helper.atomic(metadata,identity);path.touch(exist_ok=False)
    received=path.stat().st_size;need(received<=identity['length'],'Range exceeds expected length')
    if complete.exists():
        receipt=helper.load(complete)
        need(receipt.get('identity')==identity and receipt.get('bytes')==received==identity['length'] and receipt.get('local_mtime_ns')==path.stat().st_mtime_ns,'Completed range changed; retained without overwrite')
    progress(received)
    with sftp.open(entry['ServerPath'],'rb') as remote,path.open('ab') as out:
        while received<identity['length']:
            need(not stop.is_set() and helper.now()<deadline,'Range deadline/peer failure; partial retained')
            need(pause is None or not pause.exists(),'User pause marker; range retained')
            amount=min(8*1024**2,identity['length']-received);got=0
            need(shutil.disk_usage(path.parent).free>=base.RESERVE+amount,'15GiB runtime range reserve required')
            for block in blocks(remote,identity['offset']+received,amount):
                need(not stop.is_set() and helper.now()<deadline,'Range stopped; partial retained')
                need(len(block)==min(32768,amount-got),'Incomplete contiguous range response')
                out.write(block);got+=len(block);progress(received+got)
            need(got==amount,'Short range response');out.flush();received+=got
        os.fsync(out.fileno())
    after=sftp.stat(entry['ServerPath'])
    need((after.st_size,after.st_mtime)==(identity['ExpectedBytes'],identity['RemoteMtime']),'Remote immutable source changed')
    helper.atomic(complete,dict(status='COMPLETED_RANGE_SIZE_ONLY',identity=identity,bytes=received,local_mtime_ns=path.stat().st_mtime_ns,actual_sha256=None))
    return complete


def verify_completed_ranges(layout,paths,total,helper):
    expected=layout['prefix_bytes']
    need(len(paths)==len(layout['ranges'])==16,'Exact original sixteen-range layout required')
    for desc,path in zip(layout['ranges'],paths):
        receipt=helper.load(path.with_name(path.name+'.completed.json'))
        need(desc['offset']==expected and path.stat().st_size==desc['length'] and receipt.get('status')=='COMPLETED_RANGE_SIZE_ONLY' and receipt.get('identity')==dict(layout['source_identity'],**desc) and receipt.get('local_mtime_ns')==path.stat().st_mtime_ns,'All exact completed contiguous range receipts required')
        expected+=desc['length']
    need(expected==total,'Range coverage has a gap or overlap')


def assembly_marker(layout,total):
    return dict(schema_version=1,status='ASSEMBLY_AUTHORIZED_AFTER_ALL_RANGES_COMPLETE',layout=layout,total_bytes=total)


def verify_assembled_tail(prefix,layout,paths,helper,deadline):
    current=prefix.stat().st_size
    need(layout['prefix_bytes']<=current<=layout['source_identity']['ExpectedBytes'],'Prefix outside original assembly bounds')
    left=current-layout['prefix_bytes']
    with prefix.open('rb') as existing:
        existing.seek(layout['prefix_bytes'])
        for path in paths:
            with path.open('rb') as stream:
                remaining=min(left,path.stat().st_size)
                while remaining:
                    need(helper.now()<deadline,'Resume verification deadline; files retained')
                    block=stream.read(min(8*1024**2,remaining))
                    need(block and existing.read(len(block))==block,'Existing appended prefix differs from completed source ranges')
                    left-=len(block);remaining-=len(block)
            if not left:break
    need(left==0,'Appended prefix extends beyond completed ranges')


def assemble(prefix,layout,paths,total,helper,deadline,marker=None):
    marker=marker or prefix.with_name(prefix.name+'.assembly-phase.json')
    need(not marker.is_symlink(),'Assembly marker symlink refused')
    verify_completed_ranges(layout,paths,total,helper)
    if marker.exists():
        need(helper.load(marker)==assembly_marker(layout,total),'Assembly marker source/layout changed')
        verify_assembled_tail(prefix,layout,paths,helper,deadline)
    else:
        need(prefix.stat().st_size==layout['prefix_bytes'],'Prefix changed during download without authorized assembly marker')
        helper.atomic(marker,assembly_marker(layout,total))
    current=prefix.stat().st_size
    with prefix.open('ab') as out:
        for desc,path in zip(layout['ranges'],paths):
            if current>=desc['offset']+desc['length']:continue
            offset=max(0,current-desc['offset'])
            with path.open('rb') as stream:
                stream.seek(offset)
                for block in iter(lambda:stream.read(8*1024**2),b''):
                    need(helper.now()<deadline,'Assembly deadline; ranges and prefix retained');need(shutil.disk_usage(prefix.parent).free>=base.RESERVE+len(block),'15GiB assembly reserve required');out.write(block);current+=len(block)
                    out.flush();os.fsync(out.fileno())
    need(prefix.stat().st_size==total,'Assembled source size differs')


def pid_alive(pid):
    need(pid>0,'Exact previous writer PID required')
    kernel=ctypes.WinDLL('kernel32',use_last_error=True);kernel.OpenProcess.restype=ctypes.c_void_p
    handle=kernel.OpenProcess(0x1000,False,pid)
    if not handle:
        need(ctypes.get_last_error()==87,'Previous writer identity could not be safely checked');return False
    try:
        code=ctypes.c_ulong();need(kernel.GetExitCodeProcess(ctypes.c_void_p(handle),ctypes.byref(code)),'Cannot check previous writer');return code.value==259
    finally:kernel.CloseHandle(ctypes.c_void_p(handle))


def run(plan,auth=None,execute=False,previous_writer_pid=None):
    deadline=base.validate(plan)
    if not execute:return dict(status='PLAN_ONLY_SIXTEEN_RANGES_NO_AUTH_OR_NETWORK',connections=16,deadline_utc=deadline.isoformat())
    need(os.name=='nt' and dt.datetime.now(dt.timezone.utc)<deadline,'Existing Windows runtime and original deadline required')
    need(previous_writer_pid is not None and not pid_alive(previous_writer_pid),'Previous owned sequential writer must already be stopped by ROOT')
    destination=Path(plan['destination']);state=Path(plan['state_directory']);ranges=state/'parallel-ranges-v1'
    need(destination.is_dir() and not any(p.is_symlink() or (p/'.git').exists() for p in (destination,*destination.parents,state,*state.parents)),'Private existing destination outside Git/symlinks required')
    entry=next(e for e in plan['files'] if e['DownloadFileName']==base.NAME);index=next(e for e in plan['files'] if e['DownloadFileName'].endswith('.tbi'))
    final=destination/base.NAME;prefix=destination/(base.NAME+'.part');need(not final.exists() and prefix.is_file() and not prefix.is_symlink(),'Existing incomplete native prefix required; final cannot be overwritten')
    check_space(destination,entry['SizeBytes']-prefix.stat().st_size)
    need(not ranges.is_symlink(),'Range state directory symlink refused');ranges.mkdir(parents=True,exist_ok=True)
    spec=importlib.util.spec_from_file_location('range_existing_native_helper',Path(__file__).with_name('sync_completed_small_results.py'));helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    import msvcrt
    locks=[]
    for path in (state/'final-pair.lock',ranges/'ranges.lock'):
        lock=path.open('a+b');lock.write(b'0');lock.flush();lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1);locks.append(lock)
    stop=threading.Event();clients=[];mutex=threading.Lock();progresses={i:0 for i in range(16)};last=[0.];pool=None;layout=None
    def status(label,force=False,**extra):
        with mutex:
            if not force and time.monotonic()-last[0]<2:return
            last[0]=time.monotonic();value=dict(status=label,deadline_utc=deadline.isoformat(),received_bytes=layout['prefix_bytes']+sum(progresses.values()) if layout else prefix.stat().st_size,expected_bytes=entry['SizeBytes'],range_bytes=dict(progresses),patient_hashes_computed=False);value.update(extra);helper.atomic(ranges/'status.json',value)
    try:
        import paramiko
        need(paramiko.__version__=='5.1.0','Existing reviewed Paramiko runtime required')
        from native_sftp_download import connect,blocks
        config=helper.load(auth)
        def bounded_connect():
            for attempt in range(3):
                need(not stop.is_set() and helper.now()+dt.timedelta(seconds=50)<deadline,'Connection deadline/peer stop')
                try:return connect(config)
                except Exception:
                    if attempt==2:raise
                    stop.wait(.5*(attempt+1))
            raise ValueError('Bounded connection attempts exhausted')
        client,sftp=bounded_connect();clients.append((client,sftp));sftp.get_channel().settimeout(15)
        before=sftp.stat(entry['ServerPath']);need(before.st_size==entry['SizeBytes'],'Remote size differs from captured receipt')
        identity=dict(RelativePath=entry['RelativePath'],RemotePath=entry['ServerPath'],ExpectedBytes=entry['SizeBytes'],ExpectedServerSha256=entry['Sha256'],RemoteMtime=before.st_mtime,SourceReceipt=entry['SourceReceipt'],Verification='size_only')
        need(helper.load(prefix.with_name(prefix.name+'.source.json'))==identity,'Existing native prefix identity differs')
        layout_path=ranges/'ranges-layout.json'
        need(not layout_path.is_symlink(),'Range layout symlink refused')
        if layout_path.exists():
            layout=helper.load(layout_path);need(layout['source_identity']==identity,'Range source identity changed; retained for review')
            marker=ranges/'assembly-phase.json'
            if prefix.stat().st_size!=layout['prefix_bytes']:
                need(marker.is_file() and helper.load(marker)==assembly_marker(layout,entry['SizeBytes']),'Prefix changed during download without assembly marker')
                verify_completed_ranges(layout,[ranges/('range-%02d.part'%i) for i in range(16)],entry['SizeBytes'],helper)
        else:
            layout=dict(source_identity=identity,prefix_bytes=prefix.stat().st_size,ranges=partition(prefix.stat().st_size,entry['SizeBytes']));helper.atomic(layout_path,layout)
        paths=[ranges/('range-%02d.part'%i) for i in range(16)]
        def worker(desc,path,existing=None):
            c=s=None
            try:
                if existing is None:stop.wait(.2*desc['index'])
                c,s=existing if existing else bounded_connect()
                if existing is None:
                    with mutex:clients.append((c,s))
                need(not stop.is_set(),'Peer failed before range start');s.get_channel().settimeout(min(15,max(1,(deadline-helper.now()).total_seconds())))
                actual=s.stat(entry['ServerPath']);need((actual.st_size,actual.st_mtime)==(before.st_size,before.st_mtime),'Source changed before range')
                if desc['index']==0:helper.download_one(s,blocks,index,destination,deadline,Path(plan['pause_file']),reserve=base.RESERVE)
                def progress(value):
                    with mutex:progresses[desc['index']]=value
                    status('DOWNLOADING_FINAL_VCF_SIXTEEN_RANGES_SIZE_ONLY')
                return transfer_range(s,blocks,entry,dict(identity,**desc),path,helper,deadline,stop,progress,Path(plan['pause_file']))
            finally:
                if existing is None and c is not None:
                    try:s.close();c.close()
                    except Exception:pass
        pool=ThreadPoolExecutor(max_workers=16);futures=[pool.submit(worker,desc,path,(client,sftp) if i==0 else None) for i,(desc,path) in enumerate(zip(layout['ranges'],paths))]
        for future in as_completed(futures,timeout=max(1,(deadline-helper.now()).total_seconds())):future.result()
        after=sftp.stat(entry['ServerPath']);need((after.st_size,after.st_mtime)==(before.st_size,before.st_mtime),'Source changed before assembly')
        status('ASSEMBLING_COMPLETED_ORDERED_RANGES',force=True);assemble(prefix,layout,paths,entry['SizeBytes'],helper,deadline,marker=ranges/'assembly-phase.json')
        after=sftp.stat(entry['ServerPath']);need((after.st_size,after.st_mtime)==(before.st_size,before.st_mtime) and not final.exists(),'Source changed or final appeared')
        need(helper.now()<deadline,'Original deadline before final rename');os.rename(prefix,final)
        helper.atomic(final.with_name(final.name+'.completed.json'),dict(Status='COMPLETED_SIZE_ONLY',SourceIdentity=identity,LocalMtimeNs=final.stat().st_mtime_ns,ActualSha256=None,CompletedUtc=helper.now().isoformat()))
        status('PASS_FINAL_VCF_SIXTEEN_RANGES_SIZE_ONLY',force=True,range_files_retained=True,received_bytes=entry['SizeBytes'])
        return dict(status='PASS_FINAL_VCF_SIXTEEN_RANGES_SIZE_ONLY',files=2,range_files_retained=True)
    except BaseException as error:
        status('FAILED_ALL_PREFIX_AND_RANGES_RETAINED',force=True,error_type=type(error).__name__);raise
    finally:
        stop.set()
        for c,s in clients:
            try:s.close();c.close()
            except Exception:pass
        if pool is not None:pool.shutdown(wait=True,cancel_futures=True)
        for lock in reversed(locks):lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--auth-config',type=Path);p.add_argument('--previous-writer-pid',type=int);p.add_argument('--execute',action='store_true');a=p.parse_args()
    try:print(json.dumps(run(base.load(a.plan),a.auth_config,a.execute,a.previous_writer_pid)))
    except BaseException as error:print('FAILED: '+type(error).__name__+'; prefix/ranges retained; no connection values logged',file=sys.stderr);sys.exit(1)
