#!/usr/bin/env python3
"""Allowlisted modest public reference downloads, never a patient query or stage binding."""
import argparse
import base64
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time
import urllib.error
import urllib.parse
import urllib.request

ALLOWED_PREFIXES=('https://helix-research-public.s3.amazonaws.com/mito/HelixMTdb_20200327.tsv','https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz','https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/tsv/exomes/gnomad.exomes.v4.1.allele_number_all_sites.tsv.bgz','https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/')
HARD_DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
def now():return dt.datetime.now(dt.timezone.utc)
def safe(path):
 if any(p.is_symlink() for p in (path,*path.parents)):raise ValueError('Symlink path rejected')
 return path
def load(path):return json.loads(path.read_text()) if path.is_file() else {}
def save(path,value):
 temp=safe(path.with_name(path.name+'.tmp'));safe(path);temp.write_text(json.dumps(value,indent=2)+'\n');os.replace(temp,path)
def allowed(url):
 parsed=urllib.parse.urlsplit(url)
 return parsed.scheme=='https' and not parsed.query and not parsed.fragment and any(url.startswith(prefix) for prefix in ALLOWED_PREFIXES) and '..' not in parsed.path.split('/')
def range_ok(header,offset,total):return header==f'bytes {offset}-{total-1}/{total}'
def crc32c(path):
 table=[]
 for value in range(256):
  for _ in range(8):value=(value>>1)^0x82f63b78 if value&1 else value>>1
  table.append(value)
 crc=0xffffffff
 with path.open('rb') as stream:
  for block in iter(lambda:stream.read(1024**2),b''):
   for value in block:crc=table[(crc^value)&255]^(crc>>8)
 return base64.b64encode((crc^0xffffffff).to_bytes(4,'big')).decode()
def check_format(path,kind):
 opener=gzip.open if kind in ['tsv_gzip','vcf_bgzf','bgzf_tbi'] else open
 if kind=='bgzf_tbi':
  with opener(path,'rb') as stream:
   if stream.read(4)!=b'TBI\x01':raise ValueError('Unexpected tabix magic')
  return
 lines=0;first='';tabular_header=False
 with opener(path,'rt',encoding='utf-8-sig') as stream:
  for line in stream:
   if not first:first=line
   if lines<100 and '\t' in line:tabular_header=True
   lines+=1
   if lines%10000==0 and now()>=HARD_DEADLINE:raise TimeoutError('Original budget deadline during format check')
 if lines<2 or (kind.startswith('vcf') and not first.startswith('##fileformat=VCF')):raise ValueError('Unexpected reference file format')
 if kind.startswith('tsv') and not tabular_header:raise ValueError('Expected tab-delimited source header within first 100 lines')
def download(item,root,deadline):
 name=item['id'];url=item['url'];expected=item['expected_bytes'];path=safe(root/(name+'__'+Path(urllib.parse.urlsplit(url).path).name));partial=safe(path.with_name(path.name+'.part'));receipt=root/(name+'-receipt.json');status=root/(name+'-status.json')
 def state(label,**kwargs):save(status,dict(status=label,pid=os.getpid(),updated_utc=now().isoformat(),resource=name,**kwargs))
 prior=load(receipt)
 if prior:
  if not path.is_file() or path.stat().st_size!=expected:raise ValueError('Existing source receipt/file mismatch')
  return prior
 if path.exists():raise ValueError('Unreceipted source retained, no overwrite')
 publisher_md5=item.get('publisher_md5')
 if item.get('publisher_md5_url'):
  if not allowed(item['publisher_md5_url']):raise ValueError('Unapproved metadata URL')
  with urllib.request.urlopen(item['publisher_md5_url'],timeout=30) as r:body=r.read(4097)
  if len(body)>4096:raise ValueError('Checksum metadata ceiling')
  matches=re.findall(rb'\b[0-9a-fA-F]{32}\b',body)
  if len(matches)!=1 or Path(urllib.parse.urlsplit(url).path).name.encode() not in body:raise ValueError('Publisher checksum metadata mismatch')
  publisher_md5=matches[0].decode().lower();safe(root/(name+'-publisher.md5.verbatim')).write_bytes(body)
 for attempt in range(1,5):
  try:
   if now()>=deadline:raise TimeoutError('Original budget/runtime deadline')
   with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=30) as r:
    if r.url!=url or int(r.headers['Content-Length'])!=expected:raise ValueError('Publisher source identity/size changed')
    identity=dict(url=url,bytes=expected,etag=r.headers.get('ETag'),last_modified=r.headers.get('Last-Modified'),publisher_hash_header=r.headers.get('x-goog-hash'))
   old=load(root/(name+'-identity.json'))
   if old and old!=identity:raise ValueError('Source identity changed; preserve partial')
   save(root/(name+'-identity.json'),identity)
   offset=partial.stat().st_size if partial.exists() else 0
   if offset>expected:raise ValueError('Oversized partial retained')
   if shutil.disk_usage(root).free<150*1024**3+(expected-offset):raise ValueError('Keep 150GiB VM analysis reserve')
   if offset<expected:
    headers={'Accept-Encoding':'identity'}
    if offset:headers.update(Range=f'bytes={offset}-',**{'If-Range':identity['etag'] or identity['last_modified'] or ''})
    with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=30) as r:
     if r.url!=url or r.headers.get('Content-Encoding','identity')!='identity':raise ValueError('Unexpected download response')
     if offset and (r.status!=206 or not range_ok(r.headers.get('Content-Range'),offset,expected)):raise ValueError('Safe resume unavailable')
     if not offset and r.status!=200:raise ValueError('Unexpected initial HTTP status')
     with partial.open('ab' if offset else 'xb') as stream:
      tick=0
      while True:
       if now()>=deadline:raise TimeoutError('Original budget/runtime deadline')
       block=r.read(1024**2)
       if not block:break
       offset+=len(block)
       if offset>expected:raise ValueError('Source byte ceiling exceeded')
       stream.write(block)
       if time.monotonic()-tick>5:state('DOWNLOADING_PUBLIC_REFERENCE_SOURCE_ONLY',received_bytes=offset,expected_bytes=expected,attempt=attempt);tick=time.monotonic()
      stream.flush();os.fsync(stream.fileno())
   if offset!=expected:raise OSError('Truncated source retained')
   state('VERIFYING_PUBLIC_SOURCE_CHECKSUM_AND_FORMAT',received_bytes=offset)
   sha=hashlib.sha256();md5=hashlib.md5()
   with partial.open('rb') as stream:
    for block in iter(lambda:stream.read(1024**2),b''):sha.update(block);md5.update(block)
   if publisher_md5 and publisher_md5!=md5.hexdigest():raise ValueError('Publisher MD5 mismatch')
   crc_header=re.search(r'(?:^|,)crc32c=([^,]+)',identity.get('publisher_hash_header') or '')
   if crc_header and crc32c(partial)!=crc_header[1]:raise ValueError('Publisher CRC32C mismatch')
   check_format(partial,item['format'])
   result=dict(status='DOWNLOADED_VERIFIED_PUBLIC_SOURCE_ONLY_NOT_BOUND',source=item,bytes=expected,sha256=sha.hexdigest(),publisher_md5_verified=bool(publisher_md5),publisher_crc32c_verified=bool(crc_header),transport_identity=identity,patient_inputs=False,annotation_binding_updated=False,downloaded_utc=now().isoformat())
   os.replace(partial,path);save(receipt,result);state(result.pop('status'),bytes=expected,sha256=sha.hexdigest());return result
  except ValueError as exc:state('FAIL_SOURCE_REVIEW_REQUIRED',error=str(exc));raise
  except urllib.error.HTTPError as exc:
   state('SOURCE_HTTP_ERROR',http_status=exc.code,attempt=attempt)
   if exc.code not in [408,429] and exc.code<500:raise
  except Exception as exc:state('SOURCE_NETWORK_RETRY',attempt=attempt,error_class=type(exc).__name__)
  if attempt<4 and now()<deadline:time.sleep(min(attempt*5,max(0,(deadline-now()).total_seconds())))
 raise RuntimeError('Bounded public-source retries exhausted')
def main():
 import fcntl
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--group',choices=['mt','clinvar','an'],required=True);p.add_argument('--execute',action='store_true');a=p.parse_args();cfg=load(a.config)
 if cfg['patient_inputs_allowed'] is not False or cfg['annotation_binding_updated'] is not False:raise ValueError('Source-only policy required')
 if sum(x['expected_bytes'] for x in cfg['sources'])>cfg['max_total_source_bytes'] or cfg['max_total_source_bytes']>1300000000:raise ValueError('Modest source ceiling required')
 if any(not allowed(x['url']) for x in cfg['sources']):raise ValueError('Unapproved public source')
 root=safe(Path(cfg['source_root']))
 if root.name!='medinter-population-mt-sources-v1' or not root.is_absolute():raise ValueError('Wrong isolated resource directory')
 if not a.execute:print(json.dumps(dict(status='DRY_RUN_NO_NETWORK',expected_total_bytes=sum(x['expected_bytes'] for x in cfg['sources'])),indent=2));return
 os.umask(0o077);root.mkdir(parents=True,exist_ok=True);deadline=min(HARD_DEADLINE,dt.datetime.fromisoformat(cfg['deadline_utc']),now()+dt.timedelta(minutes=45))
 with safe(root/('.'+a.group+'.lock')).open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  for item in cfg['sources']:
   if item['group']==a.group:download(item,root,deadline)
if __name__=='__main__':main()
