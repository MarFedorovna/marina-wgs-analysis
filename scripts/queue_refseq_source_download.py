#!/usr/bin/env python3
"""Wait for primary VEP PASS, then bounded source-only full RefSeq-cache download."""
import datetime as dt
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
import urllib.request
from download_external_reference_source import safe,save,load,identity,validate_range,sha256

ROOT=Path('/srv/marina/resources/refseq-vep115-source-only-v1')
PRIMARY=Path('/srv/marina/resources/annotation/vep115.2-cache115-clinvar20261025-clingen20261105-publicextras-v2/resource-manifest.json')
BASE='https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/'
NAME='homo_sapiens_refseq_vep_115_GRCh38.tar.gz'
EXPECTED=24971548035
DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
def now():return dt.datetime.now(dt.timezone.utc)
def space_requirement(remaining):return 150*1024**3+128*1024**3+remaining+10*1024**3
def valid_primary(manifest):return manifest.get('status')=='PASS' and manifest.get('assembly')=='GRCh38' and manifest.get('vep',{}).get('cache_version')==115
def main():
 import fcntl
 os.umask(0o077);safe(ROOT);ROOT.mkdir(exist_ok=True,parents=True)
 with (ROOT/'.queue.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  def state(label,**extra):
   value=dict(status=label,pid=os.getpid(),updated_utc=now().isoformat(),deadline_utc=DEADLINE.isoformat(),primary_manifest=str(PRIMARY),extracted=False,clinical_pipeline_enabled=False)
   value.update(extra);save(ROOT/'status.json',value)
  target=safe(ROOT/NAME);partial=safe(ROOT/(NAME+'.partial'))
  receipt=load(ROOT/'resource-receipt.json')
  if receipt:
   if target.stat().st_size!=EXPECTED or sha256(target)!=receipt['sha256']:raise ValueError('Existing RefSeq source mismatch')
   return
  if target.exists():raise ValueError('Unreceipted source cannot be overwritten')
  while now()<DEADLINE:
   if PRIMARY.exists():
    safe(PRIMARY);manifest=load(PRIMARY)
    if not valid_primary(manifest):state('FAIL_PRIMARY_MANIFEST_INCOMPATIBLE');return
    break
   state('WAITING_FOR_PRIMARY_VEP',available_gib=round(shutil.disk_usage(ROOT).free/1024**3,2),reserve_analysis_gib=150,conservative_expanded_refseq_allowance_gib=128)
   time.sleep(min(30,max(0,(DEADLINE-now()).total_seconds())))
  else:state('STOPPED_ORIGINAL_BUDGET_DEADLINE');return
  started=now();deadline=min(DEADLINE,started+dt.timedelta(hours=2))
  for attempt in range(1,5):
   try:
    if now()>=deadline:raise TimeoutError('Original deadline/download runtime cap')
    state('PROBING_OFFICIAL_REFSEQ_SOURCE',attempt=attempt)
    with urllib.request.urlopen(BASE+'CHECKSUMS',timeout=30) as r:checks=r.read(1_000_001)
    if len(checks)>1_000_000:raise ValueError('Checksum metadata too large')
    rows=[x.split() for x in checks.decode().splitlines() if x.split()[-1:]==[NAME]]
    if len(rows)!=1 or len(rows[0])!=3:raise ValueError('No unique publisher checksum')
    publisher_sum,publisher_blocks=map(int,rows[0][:2])
    if publisher_sum!=476 or publisher_blocks!=math.ceil(EXPECTED/1024):raise ValueError('Publisher snapshot changed')
    safe(ROOT/'CHECKSUMS.verbatim').write_bytes(checks)
    with urllib.request.urlopen(urllib.request.Request(BASE+NAME,method='HEAD'),timeout=30) as r:
     if r.url!=BASE+NAME:raise ValueError('Unexpected redirect')
     current=identity(r.headers,r.url)
    if current['bytes']!=EXPECTED:raise ValueError('Publisher archive byte size changed')
    previous=load(ROOT/'source-identity.json')
    if previous and current!=previous:raise ValueError('Publisher source identity changed')
    save(ROOT/'source-identity.json',current)
    offset=partial.stat().st_size if partial.exists() else 0
    if offset>EXPECTED:raise ValueError('Oversized partial')
    free=shutil.disk_usage(ROOT).free
    if free<space_requirement(EXPECTED-offset):
     state('WAITING_FOR_DISK_RESERVE',available_bytes=free,required_bytes=space_requirement(EXPECTED-offset),reserve_analysis_gib=150,conservative_expanded_refseq_allowance_gib=128)
     # No competing writer or deletion is started. A later explicit retry can recheck space.
     return
    if offset<EXPECTED:
     headers={'Accept-Encoding':'identity'}
     if offset:headers.update({'Range':f'bytes={offset}-','If-Range':current['etag'] or current['last_modified'] or ''})
     with urllib.request.urlopen(urllib.request.Request(BASE+NAME,headers=headers),timeout=30) as r:
      if r.url!=BASE+NAME or r.headers.get('Content-Encoding','identity')!='identity':raise ValueError('Unexpected source response')
      if offset:
       if r.status!=206:raise ValueError('Safe resume unavailable')
       validate_range(r.headers.get('Content-Range'),offset,EXPECTED)
      elif r.status!=200:raise ValueError('Unexpected HTTP status')
      with partial.open('ab' if offset else 'xb') as stream:
       tick=0
       while True:
        if now()>=deadline:raise TimeoutError('Original deadline/download runtime cap')
        block=r.read(8*1024**2)
        if not block:break
        offset+=len(block)
        if offset>EXPECTED:raise ValueError('Byte ceiling exceeded')
        stream.write(block)
        if time.monotonic()-tick>15:
         if shutil.disk_usage(ROOT).free<space_requirement(EXPECTED-offset):raise ValueError('Disk reserve lost: partial retained')
         state('DOWNLOADING_REFSEQ_SOURCE_ONLY',attempt=attempt,received_bytes=offset,expected_bytes=EXPECTED);tick=time.monotonic()
       stream.flush();os.fsync(stream.fileno())
    if offset!=EXPECTED:raise OSError('Incomplete source retained')
    state('VERIFYING_PUBLISHER_REFSEQ_CHECKSUM',received_bytes=offset)
    left=min(600,int((DEADLINE-now()).total_seconds()))
    if left<=0:raise TimeoutError('Original deadline')
    got=subprocess.check_output(['sum',str(partial)],text=True,timeout=left).split()
    if [int(x) for x in got[:2]]!=[publisher_sum,publisher_blocks]:raise ValueError('Publisher BSD checksum mismatch')
    digest=sha256(partial)
    result=dict(status='DOWNLOADED_REFSEQ_SOURCE_ONLY_NOT_INSTALLED',url=BASE+NAME,bytes=offset,sha256=digest,publisher_checksum=dict(algorithm='BSD_sum_noncryptographic',sum=publisher_sum,blocks=publisher_blocks,verified=True),primary_manifest_sha256=sha256(PRIMARY),transport_identity=current,extracted=False,cache_version=115,assembly='GRCh38',clinical_pipeline_enabled=False,requires_next_steps=['validate archive expanded size and integrity','extract into isolated refseq cache, never overwrite Ensembl cache','VEP115 --refseq matching public controls'],retrieved_utc=now().isoformat())
    os.replace(partial,target);save(ROOT/'resource-receipt.json',result);state(result.pop('status'),**result);return
   except ValueError as exc:state('FAIL_REFSEQ_REVIEW_REQUIRED',error=str(exc));return
   except Exception as exc:
    state('REFSEQ_NETWORK_RETRY',attempt=attempt,error=str(exc)[:240])
    if attempt<4 and now()<deadline:time.sleep(min(15*attempt,max(0,(deadline-now()).total_seconds())))
  state('STOPPED_BOUNDED_RETRIES_OR_DEADLINE')
if __name__=='__main__':main()
