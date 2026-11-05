#!/usr/bin/env python3
"""Resume a size-only local BAM assembly against a separately trusted SSH manifest.

Never hashes BAM/chunk content, replaces a completed file, changes credentials,
or stops another downloader. Browser partials and native SFTP partials are refused.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import time

RESERVE_BYTES=15*2**30


def require(ok,message):
    if not ok:raise ValueError(message)


def basename(value):
    require(isinstance(value,str) and value and Path(value).name==value and
            not any(c in value for c in '/\\:') and not any(ord(c)<32 for c in value),'Unsafe basename')
    return value


def atomic_json(path,value):
    temp=path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
    os.replace(temp,path)


def acquire_lock(path):
    stream=path.open('a+b')
    if stream.tell()==0:
        stream.write(b'0');stream.flush()
    stream.seek(0)
    try:
        if os.name=='nt':
            import msvcrt
            msvcrt.locking(stream.fileno(),msvcrt.LK_NBLCK,1)
        else:
            import fcntl
            fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except Exception:
        stream.close();raise
    return stream


def validate_manifest(plan):
    require(plan.get('schema_version')==1 and plan.get('verification')=='size_only' and
            plan.get('content_hashes_computed') is False,'Expected trusted size-only manifest')
    basename(plan['download_filename'])
    require(isinstance(plan['size_bytes'],int) and plan['size_bytes']>0,'Invalid full size')
    require(isinstance(plan.get('export_id'),str) and len(plan['export_id'])==24 and
            all(c in '0123456789abcdef' for c in plan['export_id']),'Invalid export metadata identity')
    chunks=plan['chunks']
    require(chunks and len(chunks)==plan['chunk_count'],'Wrong chunk count')
    offset=0;names=set()
    for index,chunk in enumerate(chunks):
        name=basename(chunk['filename'])
        require(name not in names,'Duplicate chunk basename')
        require(chunk['index']==index and chunk['offset']==offset,'Noncontiguous chunk order/offset')
        require(isinstance(chunk['size_bytes'],int) and 0<chunk['size_bytes']<=plan['chunk_bytes'],'Invalid chunk size')
        if index<len(chunks)-1:require(chunk['size_bytes']==plan['chunk_bytes'],'Short nonfinal chunk')
        names.add(name);offset+=chunk['size_bytes']
    require(offset==plan['size_bytes'],'Chunk total differs from source size')
    require(plan['source_identity']['size_bytes']==plan['size_bytes'],'Source snapshot size differs')
    return plan


def check_download(chunk,downloads,stable_seconds):
    target=downloads/chunk['filename']
    if not target.exists():return None
    require(not target.is_symlink() and target.is_file(),'Unsafe chunk source')
    require(target.stat().st_size==chunk['size_bytes'],'Downloaded chunk has wrong size')
    if any(Path(str(target)+suffix).exists() for suffix in ('.part','.crdownload','.tmp')):return None
    before=target.stat()
    if time.time()-before.st_mtime<stable_seconds:return None
    return target,(before.st_size,before.st_mtime_ns)


def assemble(plan,downloads,destination,execute=False,stable_seconds=2):
    plan=validate_manifest(plan)
    downloads=Path(downloads).resolve(strict=True)
    destination=Path(destination).resolve(strict=True)
    require(stable_seconds>=0,'Negative stable interval')
    require(not any((p/'.git').exists() for p in (destination,*destination.parents)),'BAM destination cannot be Git')
    final=destination/plan['download_filename']
    if final.exists():
        require(not final.is_symlink() and final.is_file() and final.stat().st_size==plan['size_bytes'],'Existing final differs; never overwrite')
        return dict(status='EXISTING_FINAL_SIZE_ONLY',final_name=final.name,size_bytes=final.stat().st_size,content_hashes_computed=False)
    native_partial=final.with_name(final.name+'.part')
    require(not native_partial.exists(),'Native SFTP partial exists; do not interfere with its writer')
    ready=[]
    for chunk in plan['chunks']:
        if check_download(chunk,downloads,stable_seconds):ready.append(chunk['index'])
    if not execute:
        return dict(status='DRY_RUN',ready_indices=ready,chunk_count=plan['chunk_count'],content_hashes_computed=False)
    os.umask(0o077)
    partial=final.with_name(final.name+'.chunk-assembly.part')
    journal=partial.with_name(partial.name+'.json')
    lock=partial.with_name(partial.name+'.lock')
    lock_stream=acquire_lock(lock)
    try:
        state_identity={key:plan[key] for key in ('export_id','source_identity','size_bytes','download_filename','chunk_bytes','chunks')}
        if partial.exists():
            require(not partial.is_symlink() and journal.is_file(),'Owned assembly provenance required')
            state=json.loads(journal.read_text(encoding='utf-8-sig'))
            require(state.get('identity')==state_identity,'Unrelated assembly provenance')
        else:
            require(not journal.exists(),'Stale assembly journal without data')
            require(shutil.disk_usage(destination).free>=plan['size_bytes']+RESERVE_BYTES,'PC15GiB reserve required')
            partial.touch(exist_ok=False)
            state={'identity':state_identity,'completed_chunks':0,'committed_bytes':0,'verification':'size_only','content_hashes_computed':False}
            atomic_json(journal,state)
        completed=state['completed_chunks'];committed=state['committed_bytes']
        require(isinstance(completed,int) and 0<=completed<=len(plan['chunks']),'Invalid assembly checkpoint')
        require(committed==sum(c['size_bytes'] for c in plan['chunks'][:completed]),'Checkpoint offset differs')
        actual=partial.stat().st_size
        require(actual>=committed,'Owned assembly shorter than committed checkpoint')
        # Recover an interrupted append only for this exact owned manifest.
        if actual>committed:
            require(completed<len(plan['chunks']) and actual<=committed+plan['chunks'][completed]['size_bytes'],'Unexpected assembly tail')
            with partial.open('r+b') as stream:stream.truncate(committed)
        require(shutil.disk_usage(destination).free>=plan['size_bytes']-committed+RESERVE_BYTES,'PC15GiB reserve required')
        with partial.open('ab') as out:
            for chunk in plan['chunks'][completed:]:
                require(not native_partial.exists(),'Native writer started; separate assembly retained')
                checked=check_download(chunk,downloads,stable_seconds)
                if checked is None:break
                source,before=checked
                with source.open('rb') as stream:
                    shutil.copyfileobj(stream,out,length=4*2**20)
                after=source.stat()
                require(before==(after.st_size,after.st_mtime_ns),'Chunk changed while assembling')
                out.flush();os.fsync(out.fileno())
                state['completed_chunks']+=1;state['committed_bytes']+=chunk['size_bytes']
                require(partial.stat().st_size==state['committed_bytes'],'Assembly length differs')
                atomic_json(journal,state)
        if state['completed_chunks']<plan['chunk_count']:
            return dict(status='WAITING_CHUNKS',completed_chunks=state['completed_chunks'],chunk_count=plan['chunk_count'],assembled_bytes=state['committed_bytes'],next_chunk=plan['chunks'][state['completed_chunks']]['filename'],content_hashes_computed=False)
        require(partial.stat().st_size==plan['size_bytes'],'Final assembly size differs')
        require(not native_partial.exists(),'Native writer active; no promotion')
        os.link(partial,final)  # Atomic no-replace on NTFS/POSIX; same destination volume.
        partial.unlink()
        result=dict(status='COMPLETED_SIZE_ONLY',Name=final.name,SizeBytes=plan['size_bytes'],
                    export_id=plan['export_id'],verification='size_only',PatientHashesComputed=False,
                    source_historical_sha256_preserved_not_recomputed=plan.get('historical_source_sha256_preserved_not_recomputed'),
                    CompletedUtc=dt.datetime.now(dt.timezone.utc).isoformat())
        atomic_json(final.with_name(final.name+'.completed.json'),result)
        return result
    finally:
        lock_stream.close()  # OS releases the lock after normal exit or a crash.


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True,help='Private trusted manifest delivered through pinned SSH')
    parser.add_argument('--download-directory',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    parser.add_argument('--stable-seconds',type=int,default=2)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    plan=json.loads(args.manifest.read_text(encoding='utf-8-sig'))
    print(json.dumps(assemble(plan,args.download_directory,args.destination,args.execute,args.stable_seconds)))


if __name__=='__main__':main()
