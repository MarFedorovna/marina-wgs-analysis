#!/usr/bin/env python3
"""Public39 gene coverage context only, coordinated sequential gnomAD API requests."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
import urllib.error
import urllib.request
from medinter_population_mt_sources import safe,save,load,now,HARD_DEADLINE
from medinter_population_mt_common import panel

ROOT=Path('/srv/marina/resources/medinter-population-gene-coverage-v1')
ENDPOINT='https://gnomad.broadinstitute.org/api/'
QUERY='''query PublicGeneCoverage($symbol:String!) { gene(gene_symbol:$symbol,reference_genome:GRCh38) {
reference_genome gene_id gene_version symbol hgnc_id chrom start stop
coverage(dataset:gnomad_r4) { exome { pos mean median over_10 over_20 over_30 } genome { pos mean median over_10 over_20 over_30 } }
} }'''
def validate(gene,expected):
 if gene is None or gene.get('hgnc_id')!=expected['hgnc_id'] or gene.get('reference_genome')!='GRCh38':raise ValueError('Public gene identity/build mismatch')
 coverage=gene['coverage'];counts={}
 for cohort in ['exome','genome']:
  rows=coverage[cohort];counts[cohort]=len(rows)
  for row in rows:
   if not isinstance(row['pos'],int) or row['pos']<1:raise ValueError('Coverage source position invalid')
   for name in ['mean','median','over_10','over_20','over_30']:
    value=row.get(name)
    if value is not None and (not math.isfinite(value) or value<0 or (name.startswith('over_') and value>1)):raise ValueError('Coverage value invalid; never coerce missing to zero')
 return counts
def main():
 import fcntl
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');a=p.parse_args()
 if not a.execute:print('DRY_RUN_NO_NETWORK');return
 os.umask(0o077);safe(ROOT);ROOT.mkdir(parents=True,exist_ok=True)
 with safe(ROOT/'.query.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  manifest,_=panel();receipts=[];missing=[];total=0;last=0
  for gene in manifest['genes']:
   symbol=gene['symbol'];receipt_path=ROOT/(symbol+'-receipt.json');path=safe(ROOT/(symbol+'.json'))
   if symbol.startswith('MT-'):
    missing.append(dict(gene=gene,reason='NUCLEAR_COVERAGE_DATASET_NOT_APPLICABLE_TO_MITOCHONDRIAL_GENE'));continue
   if now()>=HARD_DEADLINE:raise TimeoutError('Original budget deadline')
   receipt=load(receipt_path)
   if receipt:
    if not path.is_file() or path.stat().st_size!=receipt['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest()!=receipt['sha256']:raise ValueError('Existing public coverage source mismatch')
    receipts.append(receipt);total+=receipt['bytes'];continue
   body=json.dumps({'query':QUERY,'variables':{'symbol':symbol}}).encode()
   for attempt in range(1,4):
    delay=7-(time.monotonic()-last)
    if delay>0:time.sleep(delay)
    last=time.monotonic()
    try:
     request=urllib.request.Request(ENDPOINT,data=body,headers={'Content-Type':'application/json','Accept':'application/json'})
     with urllib.request.urlopen(request,timeout=45) as response:
      if response.url!=ENDPOINT:raise ValueError('Unexpected API redirect')
      raw=response.read(5_000_001)
     if len(raw)>5_000_000 or total+len(raw)>100_000_000:raise ValueError('Public coverage download byte cap')
     data=json.loads(raw)
     if data.get('errors'):raise ValueError('GraphQL schema/service error; raw source not promoted')
     counts=validate(data['data']['gene'],gene)
     path.write_bytes(raw);receipt=dict(status='PUBLIC_GENE_BINNED_COVERAGE_CONTEXT_RECEIVED',query={'query':QUERY,'variables':{'symbol':symbol}},endpoint=ENDPOINT,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),bins=counts,source_mapping='Official current code: gnomad_r4 -> gnomad_v4_exome_coverage and gnomad_v3_genome_coverage',retrieved_utc=now().isoformat(),patient_inputs=False)
     save(receipt_path,receipt);receipts.append(receipt);total+=len(raw);break
    except urllib.error.HTTPError as error:
     save(ROOT/'status.json',dict(status='PUBLIC_COVERAGE_API_HTTP_ERROR',http_status=error.code,gene=symbol,attempt=attempt,patient_inputs=False))
     if error.code not in [408,429] and error.code<500:raise
     if attempt==3:raise
     time.sleep(30)
   save(ROOT/'status.json',dict(status='FETCHING_PUBLIC_GENE_BINNED_COVERAGE',pid=os.getpid(),completed_genes=len(receipts),download_bytes=total,updated_utc=now().isoformat(),patient_inputs=False))
  result=dict(status='READY_PUBLIC39_GENE_BINNED_COVERAGE_CONTEXT_NOT_EXACT_BASE_AN_NOT_BOUND',downloaded_genes=len(receipts),download_bytes=total,not_applicable=missing,source_mapping='gnomad_r4 API uses v4 exome coverage and v3 genome coverage',source_mapping_url='https://raw.githubusercontent.com/broadinstitute/gnomad-browser/main/graphql-api/src/queries/coverage-queries.ts',binning_note='API aggregates across padded merged exons, roughly 500 bins; mean/median are bin averages, fractions are rounded up to two decimals by source. No exact individual locus depth/callability claim.',API_does_not_return_exact_release_stamp=True,cohort_counts_source='VCFv4.1.1 separately; all-sites AN v4.1 exomes separately',patient_inputs=False,annotation_binding_updated=False,automatic_benign_or_pathogenicity=False,completed_utc=now().isoformat())
  save(ROOT/'resource-receipt.json',result);save(ROOT/'status.json',result)
if __name__=='__main__':main()
