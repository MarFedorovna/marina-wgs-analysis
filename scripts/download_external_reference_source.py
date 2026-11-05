#!/usr/bin/env python3
"""Allowlisted source-only downloads. Never enables clinical annotations."""
import argparse
import datetime as dt
import hashlib
import gzip
import json
import os
from pathlib import Path
import re
import tarfile
import time
import urllib.request
import zipfile

ALLOW = {
 'genereviews': ('https://ftp.ncbi.nlm.nih.gov/pub/litarch/ca/84/gene_NBK1116.tar.gz', 'gene_NBK1116.tar.gz', 640489521, None),
 'spliceai': ('https://files.pythonhosted.org/packages/d6/2b/9dbf72fdd948cd606c21826cc3735a5beea52633dab72d95d9936a9454d4/spliceai-1.3.1-py2.py3-none-any.whl', 'spliceai-1.3.1-py2.py3-none-any.whl', 16676295, '63c633b8de6803d4ffb613e6ef62f1a896fec9fc13245c0f000ff7adea279776'),
}
FINAL_DEADLINE = dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
def now(): return dt.datetime.now(dt.timezone.utc)
def safe(path):
 for item in [path, *path.parents]:
  if item.is_symlink(): raise ValueError('Symlink path rejected')
 return path
def save(path, obj):
 safe(path); tmp=safe(path.with_name(path.name+'.tmp'))
 tmp.write_text(json.dumps(obj,indent=2)+'\n'); os.replace(tmp,path)
def load(path): return json.loads(path.read_text()) if path.exists() else {}
def sha256(path):
 digest=hashlib.sha256()
 with path.open('rb') as stream:
  for block in iter(lambda:stream.read(1024*1024),b''): digest.update(block)
 return digest.hexdigest()
def identity(headers, url):
 return dict(url=url,bytes=int(headers['Content-Length']),last_modified=headers.get('Last-Modified'),etag=headers.get('ETag'))
def validate_range(value, offset, expected):
 match=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',value or '')
 if not match or tuple(map(int,match.groups())) != (offset,expected-1,expected): raise ValueError('Resume Content-Range mismatch')
def validate_archive(path, kind):
 if kind=='spliceai':
  with zipfile.ZipFile(path) as archive:
   if archive.testzip(): raise ValueError('Wheel CRC failure')
   members=[x.filename for x in archive.infolist()]
   if len([x for x in members if x.startswith('spliceai/models/') and x.endswith('.h5')])!=5: raise ValueError('Expected five model weights')
 else:
  with tarfile.open(path,'r:gz') as archive:
   members=[]
   for entry in archive:
    if entry.issym() or entry.islnk() or entry.isdev(): raise ValueError('Unexpected archive link/device')
    members.append(entry.name)
  with gzip.open(path,'rb') as stream:
   for block in iter(lambda:stream.read(1024*1024),b''): pass
 if not members or any(Path(x).is_absolute() or '..' in Path(x).parts for x in members): raise ValueError('Unsafe archive path')
 return dict(member_count=len(members),model_files=[x for x in members if x.endswith('.h5')])
def run(kind, root):
 import fcntl
 os.umask(0o077); safe(root)
 if not root.is_absolute() or root.name != kind+'-source-only-v1': raise ValueError('Wrong isolated root')
 root.mkdir(parents=True,exist_ok=True)
 with safe(root/'.download.lock').open('a') as lock:
  try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError: return dict(status='ALREADY_RUNNING')
  url,name,expected,publisher_sha=ALLOW[kind]
  target=safe(root/name); part=safe(root/(name+'.partial')); status=root/'status.json'; receipt=root/'resource-receipt.json'
  if receipt.exists():
   prior=load(receipt)
   if target.is_file() and target.stat().st_size==expected and sha256(target)==prior['sha256']: return prior
   raise ValueError('Existing receipt/file mismatch')
  if target.exists(): raise ValueError('Unreceipted file will not be overwritten')
  deadline=min(now()+dt.timedelta(minutes=30),FINAL_DEADLINE)
  metadata_url='https://ftp.ncbi.nlm.nih.gov/pub/litarch/ca/84/' if kind=='genereviews' else 'https://pypi.org/pypi/spliceai/1.3.1/json'
  with urllib.request.urlopen(metadata_url,timeout=30) as r: metadata=r.read(2_000_001)
  if len(metadata)>2_000_000: raise ValueError('Metadata size limit')
  safe(root/'publisher-metadata.verbatim').write_bytes(metadata)
  if kind=='spliceai':
   release=json.loads(metadata); item=next(x for x in release['urls'] if x['filename']==name)
   if item['url']!=url or item['size']!=expected or item['digests']['sha256']!=publisher_sha: raise ValueError('PyPI metadata mismatch')
  else:
   if name.encode() not in metadata: raise ValueError('Archive absent from official directory')
  for attempt in range(1,5):
   try:
    if now()>=deadline: raise TimeoutError('Original budget deadline/runtime limit')
    with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=30) as r:
     if r.url!=url: raise ValueError('Unexpected redirect')
     current=identity(r.headers,url)
    if current['bytes']!=expected: raise ValueError('Publisher size changed')
    old=load(root/'source-identity.json')
    if old and old!=current: raise ValueError('Publisher identity changed')
    save(root/'source-identity.json',current)
    offset=part.stat().st_size if part.exists() else 0
    if offset>expected: raise ValueError('Oversized partial')
    if offset<expected:
     headers={'Accept-Encoding':'identity'}
     if offset: headers.update({'Range':f'bytes={offset}-','If-Range':current['etag'] or current['last_modified'] or ''})
     with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=30) as r:
      if r.url!=url: raise ValueError('Unexpected redirect')
      if offset:
       if r.status!=206: raise ValueError('Server refused safe resume')
       validate_range(r.headers.get('Content-Range'),offset,expected)
      elif r.status!=200: raise ValueError('Unexpected download status')
      if r.headers.get('Content-Encoding','identity')!='identity': raise ValueError('Encoded response rejected')
      with part.open('ab' if offset else 'xb') as out:
       tick=0
       while True:
        if now()>=deadline: raise TimeoutError('Original budget deadline/runtime limit')
        block=r.read(1024*1024)
        if not block: break
        offset+=len(block)
        if offset>expected: raise ValueError('Source byte ceiling exceeded')
        out.write(block)
        if time.monotonic()-tick>=5:
         save(status,dict(status='DOWNLOADING_SOURCE_ONLY_NOT_ENABLED',pid=os.getpid(),received_bytes=offset,expected_bytes=expected,attempt=attempt,updated_utc=now().isoformat()));tick=time.monotonic()
       out.flush(); os.fsync(out.fileno())
    if offset!=expected: raise OSError('Incomplete source retained for safe resume')
    sha=hashlib.sha256()
    with part.open('rb') as source:
     for block in iter(lambda:source.read(1024*1024),b''): sha.update(block)
    if publisher_sha and sha.hexdigest()!=publisher_sha: raise ValueError('Publisher SHA256 mismatch')
    inventory=validate_archive(part,kind)
    result=dict(status='DOWNLOADED_SOURCE_ONLY_NOT_ENABLED',resource=kind,url=url,bytes=offset,sha256=sha.hexdigest(),publisher_checksum={'algorithm':'sha256','value':publisher_sha} if publisher_sha else {'status':'NOT_PUBLISHED_IN_INSPECTED_OFFICIAL_DIRECTORY','metadata_url':metadata_url},transport_identity=current,archive_inventory=inventory,extracted=False,patient_scores_enabled=False,clinical_use_enabled=False,purpose_and_terms_review_required=True,retrieved_utc=now().isoformat())
    os.replace(part,target);save(receipt,result);save(status,result);return result
   except ValueError as exc:
    save(status,dict(status='FAIL_REVIEW_REQUIRED',error=str(exc),updated_utc=now().isoformat()));raise
   except Exception as exc:
    save(status,dict(status='NETWORK_RETRY_SOURCE_ONLY',pid=os.getpid(),attempt=attempt,error=str(exc)[:240],updated_utc=now().isoformat()))
    if attempt<4 and now()<deadline: time.sleep(min(attempt*5,max(0,(deadline-now()).total_seconds())))
  save(status,dict(status='FAIL_BOUNDED_RETRIES_OR_DEADLINE',updated_utc=now().isoformat()))
  raise RuntimeError('Source download did not complete')
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--resource',choices=ALLOW,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
 print(json.dumps(run(a.resource,a.root) if a.execute else dict(status='DRY_RUN_SOURCE_ONLY',resource=a.resource,root=str(a.root)),indent=2))
