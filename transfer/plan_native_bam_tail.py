#!/usr/bin/env python3
"""Read-only native SFTP prefix snapshot and optional tail handoff plan.

Never pauses/kills a writer, modifies a partial, reads key material or exports
patient bytes. The output is private and cannot authorize an automatic handoff.
"""
import argparse
import ctypes
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import time


def writer_state(pid):
    if os.name != 'nt':
        return 'ACTIVE' if Path('/proc',str(pid)).exists() else 'NOT_OBSERVED'
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.OpenProcess.restype=ctypes.c_void_p
    handle=kernel.OpenProcess(0x1000,False,pid)
    if not handle:
        return 'UNKNOWN_ACCESS_DENIED' if ctypes.get_last_error()==5 else 'NOT_OBSERVED'
    try:
        code=ctypes.c_ulong()
        if not kernel.GetExitCodeProcess(ctypes.c_void_p(handle),ctypes.byref(code)):
            return 'UNKNOWN_QUERY_FAILED'
        return 'ACTIVE' if code.value==259 else 'EXITED'
    finally:
        kernel.CloseHandle(ctypes.c_void_p(handle))


def snapshot(path):
    s=path.stat();return {'size_bytes':s.st_size,'mtime_ns':s.st_mtime_ns}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-config',type=Path,required=True)
    parser.add_argument('--writer-pid',type=int,required=True)
    parser.add_argument('--entry-index',type=int,default=0)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--stable-seconds',type=int,default=3)
    args=parser.parse_args()
    if not 1<=args.stable_seconds<=10:raise ValueError('Bounded1..10second snapshot required')
    config=json.loads(args.native_config.read_text(encoding='utf-8-sig'))
    entry=config['files'][args.entry_index]
    folder=Path(config['download_directory']).resolve(strict=True)
    final=folder/entry['name'];partial=final.with_name(final.name+'.part')
    provenance=partial.with_name(partial.name+'.source.json')
    if not partial.is_file() or partial.is_symlink():raise ValueError('Regular existing native partial required')
    native=json.loads(provenance.read_text(encoding='utf-8-sig'))
    expected={'RemotePath':entry['remote_path'],'ExpectedBytes':entry['size_bytes'],
              'PinnedHostKey':config['expected_fingerprint'],'Verification':'size_only'}
    identity_matches=all(native.get(k)==v for k,v in expected.items())
    status=json.loads(Path(config['status_file']).read_text(encoding='utf-8-sig'))
    before=snapshot(partial);state=writer_state(args.writer_pid)
    time.sleep(args.stable_seconds)
    after=snapshot(partial);stable=before==after
    if after['size_bytes']>entry['size_bytes']:raise ValueError('Native partial exceeds full expected size')
    flushed=status.get('CompletedBytes') if status.get('ExpectedBytes')==entry['size_bytes'] else None
    can_freeze=identity_matches and stable and state in ('EXITED','NOT_OBSERVED')
    prefix=after['size_bytes'] if can_freeze else None
    plan={'schema_version':1,'status':'ROOT_REVIEW_REQUIRED_QUIESCED_SNAPSHOT' if can_freeze else 'NATIVE_WRITER_ACTIVE_OR_UNSTABLE_NO_HANDOFF',
          'read_only':True,'writer_pid':args.writer_pid,'writer_state':state,
          'source_provenance_matches_config':identity_matches,'source_remote_path':entry['remote_path'],
          'expected_full_bytes':entry['size_bytes'],'native_partial_path':str(partial),
          'source_state_path':str(provenance),'expected_remote_mtime_seconds':native.get('RemoteMtime'),
          'partial_before':before,'partial_after':after,'size_and_mtime_stable':stable,
          'last_status_flushed_bytes':flushed,'freeze_prefix_bytes':prefix,
          'tail_start_offset_if_frozen':prefix,'tail_bytes_if_frozen':entry['size_bytes']-prefix if prefix is not None else None,
          'pc_free_bytes':shutil.disk_usage(folder).free,'reserve_bytes':15*2**30,
          'new_patient_hashes':False,'actual_handoff_authorized_by_plan':False,
          'actual_tail_chunks_exported':False,'native_writer_modified':False,
          'original_absolute_vm_deadline':'2026-11-06T03:46:48.046025+00:00',
          'next_gates':[
              'Do not stop the running native writer for this planning step',
              'ROOT first accepts second-BAM completion and decides whether tail acceleration is still useful',
              'ROOT separately requests existing pause mechanism and verifies all matching native writers exited',
              'Repeat stable prefix snapshot; confirm remote size/mtime through existing pinned SSH',
              'Prepare a new composite manifest: local prefix [0,prefix_bytes) plus remote tail [prefix_bytes,full_bytes), with no gaps/overlaps',
              'Do not feed a tail-only manifest to the current v1 full-chunk assembler',
              'Implement/review prefix-aware seeding and checkpointing before any adoption; preserve original native partial/provenance',
              'Reserve space for tail chunks plus complete assembly plus15GiB; never overwrite a final or resume two writers on one file'],
          'created_utc':dt.datetime.now(dt.timezone.utc).isoformat()}
    args.output.write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:plan[key] for key in ['status','writer_state','source_provenance_matches_config','size_and_mtime_stable','partial_after','last_status_flushed_bytes','pc_free_bytes','actual_tail_chunks_exported']}))


if __name__=='__main__':main()
