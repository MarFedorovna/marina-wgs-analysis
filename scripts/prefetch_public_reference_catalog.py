#!/usr/bin/env python3
"""Fixed public gene/identifier catalog prefetch; no arbitrary query or patient inputs."""
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.error
import urllib.request
from download_external_reference_source import safe,save

ROOT=Path('/srv/marina/resources/reference-api-public-catalogs-v1')
SOURCES={
 'lovd_myh7_gene':('https://databases.lovd.nl/shared/api/rest.php/genes/MYH7',1_000_000),
 'omim_mim2gene':('https://omim.org/static/omim/data/mim2gene.txt',10_000_000),
}
def main():
 import fcntl
 os.umask(0o077);safe(ROOT);ROOT.mkdir(parents=True,exist_ok=True)
 with (ROOT/'.prefetch.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  for name,(url,ceiling) in SOURCES.items():
   path=safe(ROOT/(name+'.verbatim'));receipt=ROOT/(name+'-receipt.json')
   if receipt.exists() and path.exists():continue
   if path.exists():raise ValueError('Unreceipted catalog retained; no overwrite')
   for attempt in range(1,4):
    try:
     with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'MarinaPublicReferencePreparation/1.0','Accept':'text/plain, application/json, application/xml'}),timeout=30) as r:
      if r.url!=url:raise ValueError('Unexpected redirect')
      body=r.read(ceiling+1);kind=r.headers.get('Content-Type','')
     if len(body)>ceiling or not body:raise ValueError('Public catalog byte ceiling/empty response')
     text=body.decode('utf-8-sig')
     if '<html' in text[:1000].lower() or '<!doctype html' in text[:1000].lower():raise ValueError('HTML/challenge response, not catalog')
     if name.startswith('lovd') and 'MYH7' not in text:raise ValueError('Expected public gene identifier absent')
     if name=='omim_mim2gene' and not ('MIM' in text[:3000] and 'Gene' in text[:3000]):raise ValueError('Unexpected identifier-map header')
     path.write_bytes(body)
     coordinate_review={}
     if name.startswith('lovd'):
      import xml.etree.ElementTree as ET
      entry=ET.fromstring(text);ns={'a':'http://www.w3.org/2005/Atom'}
      content=entry.findtext('a:content',default='',namespaces=ns)
      fields=dict(line.strip().split(':',1) for line in content.splitlines() if ':' in line)
      coordinate_review=dict(symbol=fields.get('symbol'),refseq_genomic=fields.get('refseq_genomic'),refseq_mrna=fields.get('refseq_mrna'),position_start=fields.get('position_start'),position_end=fields.get('position_end'),source_updated=entry.findtext('a:updated',namespaces=ns),grch38_variant_mapping_accepted=False)
     save(receipt,dict(status='DOWNLOADED_PUBLIC_GENE_OR_IDENTIFIER_CATALOG_NOT_VARIANT_EVIDENCE',source=url,bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),content_type=kind,patient_inputs=False,patient_data_requested=False,variant_pathogenicity_assertions=False,coordinate_review=coordinate_review,license_review_required=True,retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat()));break
    except urllib.error.HTTPError as exc:
     save(receipt,dict(status='PUBLIC_SOURCE_ACCESS_ERROR',http_status=exc.code,source=url,attempt=attempt,patient_inputs=False))
     if exc.code not in [408,429] and exc.code<500:break
    except Exception as exc:save(receipt,dict(status='PUBLIC_SOURCE_ERROR',error=str(exc)[:240],source=url,attempt=attempt,patient_inputs=False))
    if attempt<3:time.sleep(attempt*3)
   time.sleep(1)
if __name__=='__main__':main()
