#!/usr/bin/env python3
"""Bounded gnomAD v4.1.1 public HCM39 CDS range extraction, never patient queries."""
import datetime as dt
import gzip
import hashlib
import http.server
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import threading
import time
import urllib.request
from medinter_population_mt_sources import safe,save,load,now,HARD_DEADLINE,crc32c
from medinter_population_mt_common import panel,PublicFasta,population_fields,population_info_projection

ROOT=Path('/srv/marina/resources/medinter-population-hcm39-v2')
PRIOR_ROOT=Path('/srv/marina/resources/medinter-population-hcm39-v1')
CAP=1900000000
TILE=262144
class Budget:
 def __init__(self,root):
  self.root=root;self.lock=threading.RLock();seen=set();self.used=0
  for base in [root,PRIOR_ROOT]:
   for p in list(base.rglob('*.tile'))+list(base.rglob('*.tbi')):
    stat=p.stat();key=(stat.st_dev,stat.st_ino)
    if key not in seen:self.used+=stat.st_size;seen.add(key)
  self.deadline=min(HARD_DEADLINE,now()+dt.timedelta(hours=2))
 def take(self,count):
  if now()>=self.deadline or self.used+count>CAP:raise ValueError('Original deadline or actual HTTP byte cap')
  if shutil.disk_usage(self.root).free<150*1024**3+count:raise ValueError('Keep 150GiB analysis reserve')
  self.used+=count
class Source:
 def __init__(self,url,root,budget):
  self.url=url;self.root=safe(root);root.mkdir(exist_ok=True);self.budget=budget;self.stop=threading.Event();self.failure=None
  with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=30) as r:
   if r.url!=url:raise ValueError('Unreviewed source redirect')
   self.size=int(r.headers['Content-Length']);self.identity=dict(url=url,bytes=self.size,etag=r.headers.get('ETag'),last_modified=r.headers.get('Last-Modified'))
  prior=load(root/'identity.json')
  inherited=load(PRIOR_ROOT/root.name/'identity.json')
  if inherited and inherited!=self.identity:raise ValueError('Prior public cache identity changed')
  if prior and prior!=self.identity:raise ValueError('Source changed; preserve cached tiles')
  save(root/'identity.json',self.identity)
 def tile(self,start):
  with self.budget.lock:
   if self.stop.is_set():raise ValueError('Reader stopped')
   left=start//TILE*TILE;right=min(left+TILE,self.size)-1;path=safe(self.root/f'{left:015d}.tile')
   reused=PRIOR_ROOT/self.root.name/path.name
   if path.exists() or reused.exists():
    data=(path if path.exists() else safe(reused)).read_bytes()
    if len(data)!=right-left+1:raise ValueError('Cached public source tile truncated')
   else:
    self.budget.take(right-left+1)
    headers={'Range':f'bytes={left}-{right}','Accept-Encoding':'identity','If-Match':self.identity['etag']}
    with urllib.request.urlopen(urllib.request.Request(self.url,headers=headers),timeout=30) as r:
     if r.url!=self.url or r.status!=206 or r.headers.get('Content-Range')!=f'bytes {left}-{right}/{self.size}' or r.headers.get('ETag')!=self.identity['etag']:raise ValueError('Unexpected source range identity')
     data=r.read(right-left+2)
    if len(data)!=right-left+1:raise ValueError('Truncated source range')
    temporary=safe(path.with_suffix('.tmp'));temporary.write_bytes(data);os.replace(temporary,path)
   return right,data[start-left:]
 def index(self,path):
  url=self.url+'.tbi'
  with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=30) as r:
   expected=int(r.headers['Content-Length']);identity=dict(url=url,bytes=expected,etag=r.headers.get('ETag'),crc32c=r.headers.get('x-goog-hash'))
  prior=load(path.with_suffix('.identity.json'))
  if prior and prior!=identity:raise ValueError('Publisher index changed')
  save(path.with_suffix('.identity.json'),identity)
  inherited=PRIOR_ROOT/self.root.name/path.name
  if not path.exists() and inherited.exists():
   if load(inherited.with_suffix('.identity.json'))!=identity:raise ValueError('Prior index identity mismatch')
   os.link(safe(inherited),path)
  if not path.exists():
   with self.budget.lock:self.budget.take(expected)
   with urllib.request.urlopen(url,timeout=30) as r:body=r.read(expected+1)
   if len(body)!=expected:raise ValueError('Index size mismatch')
   safe(path).write_bytes(body)
  if path.stat().st_size!=expected:raise ValueError('Cached index size mismatch')
  match=re.search(r'(?:^|,)crc32c=([^,]+)',identity['crc32c'] or '')
  if match and crc32c(path)!=match[1]:raise ValueError('Publisher index CRC32C mismatch')
  return identity
def proxy(source):
 class Handler(http.server.BaseHTTPRequestHandler):
  def log_message(self,*args):pass
  def do_HEAD(self):
   self.send_response(200);self.send_header('Content-Length',str(source.size));self.send_header('Accept-Ranges','bytes');self.end_headers()
  def do_GET(self):
   try:
    request=self.headers.get('Range');match=re.fullmatch(r'bytes=(\d+)-(\d*)',request or 'bytes=0-')
    if not match:raise ValueError('Only one byte range')
    start=int(match[1]);end=min(int(match[2]) if match[2] else source.size-1,source.size-1)
    if not 0<=start<=end<source.size:raise ValueError('Bad source bounds')
    self.send_response(206 if request else 200);self.send_header('Content-Length',str(end-start+1));self.send_header('Accept-Ranges','bytes')
    if request:self.send_header('Content-Range',f'bytes {start}-{end}/{source.size}')
    self.end_headers()
    while start<=end and not source.stop.is_set():
     right,data=source.tile(start);data=data[:end-start+1];self.wfile.write(data);self.wfile.flush();start+=len(data)
   except (BrokenPipeError,ConnectionResetError):pass
   except Exception as exc:source.failure=str(exc)
 server=http.server.ThreadingHTTPServer(('127.0.0.1',0),Handler);threading.Thread(target=server.serve_forever,daemon=True).start();return server,f'http://127.0.0.1:{server.server_port}/public.vcf.bgz'
def main():
 import sys
 if '--execute' not in sys.argv:print('DRY_RUN_NO_NETWORK_OR_SOURCE_READ');return
 import fcntl
 import pysam
 os.umask(0o077);safe(ROOT);ROOT.mkdir(parents=True,exist_ok=True)
 with safe(ROOT/'.extract.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  manifest,intervals=panel();budget=Budget(ROOT);fasta=PublicFasta();db=sqlite3.connect(ROOT/'public-population.sqlite')
  db.executescript('CREATE TABLE IF NOT EXISTS variants(source TEXT,chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,filter TEXT,AC TEXT,AN INTEGER,AF TEXT,nhomalt TEXT,info TEXT,PRIMARY KEY(source,chrom,pos,ref,alt));CREATE TABLE IF NOT EXISTS completed(source TEXT,chrom TEXT,start INTEGER,end INTEGER,PRIMARY KEY(source,chrom,start,end));')
  migration=ROOT/'prior-public-source-migration.json'
  if not migration.exists() and (PRIOR_ROOT/'public-population.sqlite').exists():
   old=sqlite3.connect(f'file:{PRIOR_ROOT}/public-population.sqlite?mode=ro',uri=True)
   old.execute('BEGIN')  # One consistent public-source snapshot while prior job is alive.
   with db:
    for row in old.execute('SELECT * FROM variants'):
     values=list(row);values[-1]=population_info_projection(values[-1]);db.execute('INSERT OR IGNORE INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?)',values)
    for row in old.execute('SELECT * FROM completed'):db.execute('INSERT OR IGNORE INTO completed VALUES (?,?,?,?)',row)
   old.close();save(migration,dict(status='PRIOR_PUBLIC_SOURCE_INTERVALS_REUSED',preserved_prior_directory=str(PRIOR_ROOT),only_omitted_INFO_fields=['vep','CSQ'],frequency_ancestry_quality_fields_retained=True,patient_data=False))
  if '--stage-migration-only' in sys.argv:
   db.close();fasta.close();print('STAGED_PUBLIC_SOURCE_MIGRATION_NO_NETWORK');return
  def state(label,**extra):save(ROOT/'status.json',dict(status=label,pid=os.getpid(),actual_http_bytes=budget.used,max_http_bytes=CAP,updated_utc=now().isoformat(),patient_inputs=False,annotation_binding_updated=False,**extra))
  try:
   for dataset in ['genomes','exomes']:
    for chrom in sorted({c for c,s,e in intervals},key=lambda c:(c=='chrX',int(c[3:]) if c[3:].isdigit() else 99)):
     key=f'{dataset}-{chrom}';source_dir=safe(ROOT/key);url=f'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1.1/vcf/{dataset}/gnomad.{dataset}.v4.1.1.sites.{chrom}.vcf.bgz'
     state('FETCHING_PUBLIC_HCM39_GNOMAD_V4_1_1',source=key)
     source=Source(url,source_dir,budget);index=source_dir/'public.vcf.bgz.tbi';index_identity=source.index(index);server,local=proxy(source)
     try:
      tabix=pysam.TabixFile(local,index=str(index));headers=list(tabix.header)
      if not all(any(f'ID={field},' in h for h in headers) for field in ['AC','AN','AF','nhomalt']):raise ValueError('gnomAD source schema missing required fields')
      (source_dir/'headers.vcf.txt').write_text('\n'.join(headers)+'\n')
      for c,start,end in intervals:
       if c!=chrom:continue
       if db.execute('SELECT 1 FROM completed WHERE source=? AND chrom=? AND start=? AND end=?',(key,c,start,end)).fetchone():continue
       if now()>=budget.deadline:raise TimeoutError('Original budget/runtime deadline')
       with db:
        for line in tabix.fetch(c,start,end):
         fields=population_fields(line.split('\t'))
         if not start<fields['pos']<=end:continue
         if len(fields['ref'])<10000 and fasta.sequence(c,fields['pos'],len(fields['ref']))!=fields['ref']:raise ValueError('Public GRCh38 REF mismatch')
         values=(key,fields['chrom'],fields['pos'],fields['ref'],fields['alt'],fields['filter'],fields['AC'],fields['AN'],fields['AF'],fields['nhomalt'],population_info_projection(fields['info']))
         db.execute('INSERT OR IGNORE INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?)',values)
        if source.failure:raise ValueError(source.failure)
        db.execute('INSERT INTO completed VALUES (?,?,?,?)',(key,c,start,end))
       if (ROOT/'public-population.sqlite').stat().st_size>1400000000:raise ValueError('Public source output size cap')
       done=db.execute('SELECT COUNT(*) FROM completed').fetchone()[0]
       state('FETCHING_PUBLIC_HCM39_GNOMAD_V4_1_1',source=key,completed_intervals=done,total_intervals=2*len(intervals),rows=db.execute('SELECT COUNT(*) FROM variants').fetchone()[0])
      tabix.close();save(source_dir/'source-receipt.json',dict(status='COMPLETE_PUBLIC_CDS_SOURCE_RANGES',transport_identity=source.identity,index_identity=index_identity,fields=['AC','AN','AF','nhomalt','FILTER','all_ancestry_INFO_preserved'],full_source_publisher_checksum_verified=False,no_patient_inputs=True))
     finally:source.stop.set();server.shutdown();server.server_close()
   rows=db.execute('SELECT COUNT(*) FROM variants').fetchone()[0];counts=dict(db.execute('SELECT source,COUNT(*) FROM variants GROUP BY source'))
   result=dict(status='READY_PUBLIC_HCM39_GNOMAD_V4_1_1_SIDECAR_NOT_BOUND',rows=rows,counts=counts,scope='Fixed public39-gene GENCODE49 CDS intervals; variants selected by POS in CDS; no genome-wide/noncoding claim',panel_sha256=hashlib.sha256((Path('/srv/marina/resources/public-hcm-cds-v2')/'panel-manifest.json').read_bytes()).hexdigest(),actual_http_bytes=budget.used,only_omitted_INFO_fields=['vep','CSQ'],frequency_ancestry_quality_fields_retained=True,source_range_cache_resume=True,sqlite_transaction_interval_resume=True,reference_alignment_checked=True,missing_record_is_not_assumed_AC_zero=True,all_sites_AN_separate_release='4.1 exomes only; genomes AN not downloaded',annotation_binding_updated=False,patient_inputs=False,completed_utc=now().isoformat())
   save(ROOT/'resource-receipt.json',result);state(result.pop('status'),rows=rows)
  except Exception as exc:state('PARTIAL_PUBLIC_REFERENCE_RETAINED_NOT_BOUND',error=str(exc)[:250]);raise
  finally:db.close();fasta.close()
if __name__=='__main__':main()
