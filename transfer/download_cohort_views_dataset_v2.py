#!/usr/bin/env python3
"""Resume only an entirely physically allocated, closed v1 dataset. Default dry."""
import argparse
import ast
import ctypes
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys

BASE_SHA='d73774c1b20062858c0097d81aa6c6bf740d9a3ad7f8324eaa2aaccfe2a2731f'
FLOOR=5*1024**3

def need(ok,message):
    if not ok:raise ValueError(message)

def base_module():
    path=Path(__file__).with_name('download_cohort_views_dataset.py')
    need(hashlib.sha256(path.read_bytes()).hexdigest()==BASE_SHA,'Reviewed v1 source identity mismatch')
    spec=importlib.util.spec_from_file_location('reviewed_dataset_v1',path);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
    return m

def physical_allocation(path):
    need(os.name=='nt','Windows FileStandardInfo allocation check required')
    import msvcrt
    class Info(ctypes.Structure):
        _fields_=[('AllocationSize',ctypes.c_longlong),('EndOfFile',ctypes.c_longlong),('NumberOfLinks',ctypes.c_uint32),('DeletePending',ctypes.c_ubyte),('Directory',ctypes.c_ubyte)]
    k=ctypes.WinDLL('kernel32',use_last_error=True)
    query=k.GetFileInformationByHandleEx;query.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.c_void_p,ctypes.c_uint32];query.restype=ctypes.c_int
    attrs=k.GetFileAttributesW;attrs.argtypes=[ctypes.c_wchar_p];attrs.restype=ctypes.c_uint32
    attribute=attrs(str(path));need(attribute!=0xffffffff and not attribute&(0x200|0x800),'Sparse or compressed part refused for physical-only resume')
    info=Info()
    with Path(path).open('rb') as f:need(query(msvcrt.get_osfhandle(f.fileno()),1,ctypes.byref(info),ctypes.sizeof(info)),'FileStandardInfo physical allocation query failed')
    need(not info.Directory and not info.DeletePending and info.EndOfFile==Path(path).stat().st_size,'FileStandardInfo logical identity mismatch')
    return info.AllocationSize

def source_ledger(value,name,expected,manifest_sha):
    need(value.get('artifact')==name and value.get('bytes')==expected['bytes'] and value.get('sha256')==expected['sha256'] and value.get('manifest_sha256')==manifest_sha and type(value.get('remote_mtime')) is int and value['remote_mtime']>0,'Captured source ledger identity mismatch')

def preflight(base,plan,allocation=physical_allocation):
    deadline,files=base.validate(plan);need(base.now()<deadline,'Resume transfer deadline expired')
    destination=base.private_path(plan['destination']);state=base.private_path(plan['state_directory'])
    need(destination.is_dir() and state.is_dir() and (state/'files').is_dir(),'Resume requires existing destination and ledger directories')
    aliases={a['artifact']:a for a in plan['source_aliases']};completed=parts=allocated=committed=0
    manifest_sha=plan['source_provenance']['manifest']['sha256']
    for name,expected in files.items():
        target,part,meta,note=base.local_paths(destination,state,name)
        if target.exists():
            need(base.completed(target,note,expected),'Completed result identity mismatch')
            need(not note.is_symlink(),'Completed receipt symlink refused')
            record=base.load(note)
            if name in aliases:
                a=aliases[name];src=base.private_path(a['local_path']);original=base.load(a['completion_receipt']);identity=original.get('SourceIdentity',{})
                need(record.get('status')=='ROOT_REVIEWED_SOURCE_ALIAS_SIZE_ONLY' and src.is_file() and os.path.samefile(src,target) and original.get('Status')=='COMPLETED_SIZE_ONLY' and identity.get('ExpectedBytes')==expected['bytes'] and identity.get('ExpectedServerSha256')==expected['sha256'] and src.stat().st_mtime_ns==original.get('LocalMtimeNs'),'Existing main pair alias identity mismatch')
            else:
                need(record.get('status')=='COMPLETED_DATASET_FILE_SIZE_ONLY','Completed output receipt status mismatch')
                source_ledger(record.get('source',{}),name,expected,manifest_sha)
            actual=allocation(target);need(actual>=expected['bytes'],'Completed file physical allocation below required size')
            allocated+=actual;completed+=1;committed+=expected['bytes'];continue
        need(name not in aliases,'Existing main VCF alias is required; no new hardlink or download')
        need(part.is_file() and not part.is_symlink() and meta.is_file() and not meta.is_symlink(),'Every unfinished file requires an existing owned preallocated part')
        value=base.load(meta)
        need(value.get('schema_version')==1 and value.get('artifact')==name and value.get('expected')==expected and value.get('chunk_bytes')==base.CHUNK and value.get('file_id')==part.stat().st_ino and part.stat().st_size==expected['bytes'],'Owned preallocated part/ledger identity mismatch')
        source_ledger(value.get('source',{}),name,expected,manifest_sha)
        count=(expected['bytes']+base.CHUNK-1)//base.CHUNK;bitmap=value.get('completed')
        need(isinstance(bitmap,list) and bitmap==sorted(set(bitmap)) and all(type(i) is int and 0<=i<count for i in bitmap),'Owned durable chunk bitmap invalid')
        actual=allocation(part);need(actual>=expected['bytes'],'Owned part physical allocation below required size')
        allocated+=actual;parts+=1;committed+=sum(min(base.CHUNK,expected['bytes']-i*base.CHUNK) for i in bitmap)
    need(shutil.disk_usage(destination).free>=FLOOR,'5GiB preallocated-only runtime floor reached')
    return dict(status='PASS_RESUME_ONLY_ALL_FILES_PHYSICALLY_ALLOCATED_NO_AUTH',files=len(files),completed_files=completed,owned_parts=parts,expected_bytes=sum(x['bytes'] for x in files.values()),committed_bytes=committed,physical_allocated_bytes=allocated,additional_allocation_bytes=0,runtime_floor_bytes=FLOOR,deadline_utc=deadline.isoformat(),patient_hashes_computed=False)

def literal_reasons(base):
    texts=[Path(base.__file__).read_text(encoding='utf-8'),Path(__file__).read_text(encoding='utf-8')];reasons=set()
    for text in texts:
        for node in ast.walk(ast.parse(text)):
            if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='need' and len(node.args)>1 and isinstance(node.args[1],ast.Constant) and isinstance(node.args[1].value,str):reasons.add(node.args[1].value)
    return reasons

def safe_reason(error,base):
    if not isinstance(error,ValueError):return 'UNCLASSIFIED_'+type(error).__name__
    message=str(error)
    if message=='15GiB runtime reserve reached':return '5GiB preallocated-only runtime floor reached'
    return message if message in literal_reasons(base) else 'UNCLASSIFIED_ValueError'

def run(plan,auth=None,execute=False):
    base=base_module();proof=preflight(base,plan)
    if not execute:return proof
    def no_allocation(*_):raise ValueError('Resume-only forbids any new preallocation')
    def resume_additional(files,destination,state,aliases):
        preflight(base,plan);return 0
    base.RESERVE=FLOOR;base.allocation_bytes=physical_allocation;base.preallocate=no_allocation;base.required_additional=resume_additional
    try:return base.run(plan,auth,True,False)
    except BaseException as error:
        status=Path(plan['state_directory'])/'status.json'
        if status.is_file():
            record=base.load(status)
            if record.get('status')=='FAILED_OWN_PARTS_AND_BITMAPS_RETAINED':record['safe_error_reason']=safe_reason(error,base);record['runtime_floor_bytes']=FLOOR;record['resume_only_all_preallocated']=True;base.atomic(status,record)
        raise

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--auth-config',type=Path);p.add_argument('--execute',action='store_true');a=p.parse_args()
    try:print(json.dumps(run(json.loads(a.plan.read_text(encoding='utf-8-sig')),a.auth_config,a.execute)))
    except BaseException as error:
        try:reason=safe_reason(error,base_module())
        except Exception:reason='UNCLASSIFIED_'+type(error).__name__
        print(json.dumps(dict(status='FAILED_RESUME_ONLY_PARTS_RETAINED',error_type=type(error).__name__,safe_error_reason=reason)),file=sys.stderr);sys.exit(1)
if __name__=='__main__':main()
