#!/usr/bin/env python3
"""Four-MiB total read-only probe: one MiB serial, three MiB over pinned connections."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import threading
import time
from native_sftp_download import blocks,connect,require,save

def run(cfg,connections):
 require(connections in (4,8),'Only four or eight connections allowed')
 entry=cfg['files'][0];size=entry['size_bytes'];require(size>=4*2**20,'Probe source too small')
 client,sftp=connect(cfg)
 try:
  before=sftp.stat(entry['remote_path']);require(before.st_size==size,'Expected source size mismatch')
  started=time.monotonic();read=0
  with sftp.open(entry['remote_path'],'rb') as remote:
   for block in blocks(remote,0,2**20):read+=len(block)
  elapsed=max(1e-9,time.monotonic()-started);require(read==2**20,'Incomplete baseline probe')
 finally:sftp.close();client.close()
 baseline=dict(bytes=read,seconds=elapsed,mib_per_second=read/2**20/elapsed)
 barrier=threading.Barrier(connections);beginnings=[];endings=[];chunk=3*2**20//connections
 def worker(index):
  client,sftp=connect(cfg)
  try:
   actual=sftp.stat(entry['remote_path']);require((actual.st_size,actual.st_mtime)==(before.st_size,before.st_mtime),'Source identity changed')
   with sftp.open(entry['remote_path'],'rb') as remote:
    barrier.wait(timeout=45);start=time.monotonic();received=0
    for block in blocks(remote,2**20+index*chunk,chunk):received+=len(block)
    stop=time.monotonic();require(received==chunk,'Incomplete parallel probe')
   after=sftp.stat(entry['remote_path']);require((after.st_size,after.st_mtime)==(before.st_size,before.st_mtime),'Source identity changed')
   return received,start,stop
  except BaseException:
   barrier.abort();raise
  finally:sftp.close();client.close()
 started=time.monotonic()
 with ThreadPoolExecutor(max_workers=connections) as executor:results=list(executor.map(worker,range(connections)))
 total=sum(x[0] for x in results);elapsed=max(1e-9,max(x[2] for x in results)-min(x[1] for x in results))
 require(total+baseline['bytes']==4*2**20,'Four-MiB total ceiling mismatch')
 aggregate=total/2**20/elapsed
 return dict(status='PASS_READ_ONLY_FOUR_MIB_TOTAL_PROBE',connections=connections,baseline=baseline,parallel=dict(bytes=total,seconds=elapsed,mib_per_second=aggregate,including_connection_seconds=time.monotonic()-started),speedup_vs_measured_baseline=aggregate/baseline['mib_per_second'],total_source_bytes_read=4*2**20,pinned_full_host_key=True,patient_hashes_computed=False)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--connections',type=int,choices=[4,8],default=4);p.add_argument('--receipt',required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
 if not a.execute:print('DRY_RUN_NO_KEY_OR_NETWORK_READ')
 else:
  value=run(json.loads(Path(a.config).read_text(encoding='utf-8-sig')),a.connections);save(Path(a.receipt),value);print(json.dumps(value),flush=True)
