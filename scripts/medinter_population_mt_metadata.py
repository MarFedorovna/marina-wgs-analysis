#!/usr/bin/env python3
"""Retain public provider notices/reference metadata only; never retry an access denial."""
import hashlib
import json
import urllib.error
import urllib.request
from pathlib import Path
from medinter_population_mt_sources import now,safe,save,HARD_DEADLINE

ROOT=Path('/srv/marina/resources/medinter-population-mt-metadata-v1')
SOURCES={
 'gnomad_terms':('https://raw.githubusercontent.com/broadinstitute/gnomad-browser/main/browser/about/policies/terms.md',100000),
 'gnomad_download_manifest':('https://raw.githubusercontent.com/broadinstitute/gnomad-browser/main/browser/src/DataPage/GnomadV4Downloads.tsx',200000),
 'helix_download_page':('https://www.helix.com/mitochondrial-variant-database',1000000),
 'rcrs_NC_012920_1':('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text',50000),
 'mitomap_resources_access':('https://www.mitomap.org/MITOMAP/Resources',1000000),
}
def main():
 safe(ROOT);ROOT.mkdir(parents=True,exist_ok=True)
 for name,(url,cap) in SOURCES.items():
  receipt=ROOT/(name+'-receipt.json');target=safe(ROOT/(name+'.verbatim'))
  if receipt.exists():continue
  if now()>=HARD_DEADLINE:raise TimeoutError('Original budget deadline')
  try:
   with urllib.request.urlopen(url,timeout=30) as response:
    body=response.read(cap+1);actual_url=response.url
   if len(body)>cap:raise ValueError('Metadata byte cap')
   if name=='rcrs_NC_012920_1' and not body.startswith(b'>NC_012920.1'):raise ValueError('Expected public rCRS accession')
   target.write_bytes(body);save(receipt,dict(status='OFFICIAL_PUBLIC_METADATA_RETAINED_VERBATIM',url=url,actual_url=actual_url,bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),retrieved_utc=now().isoformat(),patient_inputs=False))
  except urllib.error.HTTPError as error:save(receipt,dict(status='HTTP_ACCESS_DENIED_NO_BYPASS' if error.code in [401,403] else 'HTTP_METADATA_ERROR',http_status=error.code,url=url,patient_inputs=False,checked_utc=now().isoformat()))
if __name__=='__main__':main()
