#!/usr/bin/env python3
"""Retain official upstream notices verbatim, with immutable Git commit provenance."""
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import urllib.request
from download_external_reference_source import safe,save

ROOT=Path('/srv/marina/resources/reference-provider-notices-v1')
def fetch(url,cap=2_000_000):
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'PublicReferenceSourcePreparation/1.0'}),timeout=30) as response:content=response.read(cap+1)
 if len(content)>cap:raise ValueError('Notice size ceiling')
 return content
def main():
 safe(ROOT);ROOT.mkdir(parents=True,exist_ok=True)
 if (ROOT/'notice-receipt.json').exists():return
 commit_url='https://api.github.com/repos/Illumina/SpliceAI/commits/master'
 metadata=fetch(commit_url);commit=json.loads(metadata)['sha']
 if not re.fullmatch('[0-9a-f]{40}',commit):raise ValueError('Unexpected upstream commit')
 (ROOT/'spliceai-commit-metadata.verbatim.json').write_bytes(metadata)
 records=[]
 for name in ['LICENSE','COPYRIGHT','NOTICE','README.md']:
  url=f'https://raw.githubusercontent.com/Illumina/SpliceAI/{commit}/{name}';body=fetch(url);path=safe(ROOT/('spliceai-'+name+'.verbatim'));path.write_bytes(body)
  records.append(dict(source=url,file=path.name,bytes=len(body),sha256=hashlib.sha256(body).hexdigest()))
 save(ROOT/'notice-receipt.json',dict(status='OFFICIAL_UPSTREAM_NOTICES_RETAINED_VERBATIM',upstream_commit=commit,scope='Current archived official repository; historical PyPI wheel metadata remains separately preserved',records=records,retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
if __name__=='__main__':main()
