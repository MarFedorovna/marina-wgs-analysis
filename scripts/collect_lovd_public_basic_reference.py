#!/usr/bin/env python3
"""Bounded fixed-public-gene LOVD retrieval; never accepts patient variants or VCFs."""
import argparse
from collections import Counter
import datetime as dt
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

BASE='https://databases.lovd.nl/shared/api/rest.php/'
CONTROL=BASE+'variants/IVD/0000000563?format=application/json'
PANEL_SHA='09c518f2659bd2443ee024d65dd83070922864304ed92f705271e250c227fedb'
END=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
def now():return dt.datetime.now(dt.timezone.utc)
def safe(path):
 path=Path(path)
 if any(p.is_symlink() for p in [path,*path.parents]):raise ValueError('Symlink source/output refused')
 return path
def sha(path):return hashlib.sha256(safe(path).read_bytes()).hexdigest()
def pin(path):
 path=safe(Path(path).absolute());return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
def load(path):return json.loads(safe(path).read_text(encoding='utf-8-sig'))
def dump(path,value):
 with safe(path).open('x',encoding='utf-8') as f:json.dump(value,f,ensure_ascii=False,indent=2);f.write('\n')
def url(kind,symbol):
 if kind not in ['genes','variants'] or not re.fullmatch('[A-Z0-9][A-Z0-9-]{0,30}',symbol):raise ValueError('Fixed gene endpoint only')
 return BASE+kind+'/'+symbol+'?format=application/json'
class NoRedirect(urllib.request.HTTPRedirectHandler):
 def redirect_request(self,*args,**kwargs):return None
def fetch(target,limit,deadline):
 if not (target==CONTROL or re.fullmatch(re.escape(BASE)+r'(genes|variants)/[A-Z0-9][A-Z0-9-]{0,30}\?format=application/json',target)):raise ValueError('Only official basic public gene retrieval endpoints allowed')
 if now()>=deadline:raise TimeoutError('Original compute deadline reached')
 request=urllib.request.Request(target,headers={'User-Agent':'MarinaWGS-PublicGeneReference/1.0 (noncommercial research; fixed public panel)','Accept':'application/json','Accept-Encoding':'identity'})
 try:
  with urllib.request.build_opener(NoRedirect).open(request,timeout=min(30,max(1,(deadline-now()).total_seconds()))) as response:
   length=response.headers.get('Content-Length')
   if length and int(length)>limit:raise ValueError('Public source response exceeds reviewed byte cap')
   body=response.read(limit+1)
   if len(body)>limit:raise ValueError('Public source response exceeds reviewed byte cap')
   return body,dict(url=target,status=response.status,content_type=response.headers.get('Content-Type',''),bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),etag=response.headers.get('ETag'),last_modified=response.headers.get('Last-Modified'),retrieved_utc=now().isoformat(),patient_query=False)
 except urllib.error.HTTPError as error:
  return None,dict(url=target,status=error.code,access_status='HTTP_UNAVAILABLE_NOT_AN_EMPTY_VARIANT_SET',retry_not_before_utc=(now()+dt.timedelta(hours=4)).isoformat(),patient_query=False)
 except (urllib.error.URLError,TimeoutError) as error:
  return None,dict(url=target,status=None,access_status='NETWORK_UNAVAILABLE_NOT_AN_EMPTY_VARIANT_SET',error_type=type(error).__name__,patient_query=False)
def rows(body,kind,symbol):
 try:value=json.loads(body.decode('utf-8-sig'))
 except Exception as error:raise ValueError('Basic endpoint did not return JSON') from error
 # Official single-ID example returns an object; gene collections return arrays.
 if isinstance(value,dict):value=[value]
 if not isinstance(value,list) or any(not isinstance(row,dict) for row in value):raise ValueError('Unexpected LOVD basic API schema')
 for row in value:
  if row.get('symbol')!=symbol:raise ValueError('Public response gene identity differs')
  if any(key.lower() in ['individuals','patients','screenings','phenotypes'] for key in row):raise ValueError('Unexpected case-level payload; refuse basic-reference claim')
  if kind=='genes' and not all(key in row for key in ['symbol','refseq_build']):raise ValueError('Gene response must declare reference build')
  if kind=='variants' and not all(key in row for key in ['id','position_genomic','Variant/DNA','Variant/DBID']):raise ValueError('Variant response lacks basic identity/description fields')
 return value
def mapping_status(row,gene_metadata):
 explicit=row.get('variant_genomic')
 if isinstance(explicit,dict) and 'hg38' in explicit:
  if re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y):g\.[0-9]+[ACGT]>[ACGT]',str(explicit['hg38'])):return 'GRCH38_SNV_CANDIDATE_PENDING_LOCAL_FASTA_VALIDATION'
  return 'UNMAPPED_EXPLICIT_GRCH38_COMPLEX_DESCRIPTION'
 build=gene_metadata.get('refseq_build') if gene_metadata else None
 if build not in ['hg38','GRCh38']:return 'UNMAPPED_SOURCE_BUILD_'+(str(build) if build else 'UNSPECIFIED')
 # A position without a DNA change is not REF/ALT. Transcript HGVS is never
 # guessed into genome coordinates and GRCh37 is never relabelled GRCh38.
 if not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y):g\.[0-9]+[ACGT]>[ACGT]',str(row.get('position_genomic',''))):return 'UNMAPPED_NO_EXPLICIT_GRCH38_REF_ALT'
 return 'UNMAPPED_PENDING_LOCAL_REFERENCE_VALIDATION'
def prepare(panel_path,profile_path,output):
 panel=load(panel_path)
 if sha(panel_path)!=PANEL_SHA or panel.get('no_patient_selection') is not True or len(panel.get('genes',[]))!=39:raise ValueError('Exact pre-existing fixed public39 panel required')
 symbols=[r['symbol'] for r in panel['genes']]
 if len(set(symbols))!=39:raise ValueError('Unique fixed public genes required')
 value=dict(schema_version=1,review_status='PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING',purpose='noncommercial_research',scope='LOVD_BASIC_PUBLIC_GENE_REFERENCE_ONLY_NOT_CASE_JSON',panel=pin(panel_path),symbols=symbols,base_url=BASE,output=str(output),requests_per_second=1,max_requests=78,max_response_bytes=16*1024**2,max_total_bytes=256*1024**2,max_minutes=10,deadline_utc=END.isoformat(),patient_variants_allowed=False,case_level_json_allowed=False,no_source_vcf_mutation=True,source_reuse_policy='Basic retrieval under official API fair-use; preserve source and per-submitter rights; raw case-level/licensed exports excluded',source_terms_url='https://www.lovd.nl/3.0/docs/LOVD_manual_3.0.pdf',code_sha256=sha(Path(__file__)))
 dump(profile_path,value);return dict(status=value['review_status'],profile=str(profile_path),sha256=sha(profile_path),request_count=78,collector_started=False)
def collect(path,digest,execute=False):
 profile=load(path)
 if not execute:return dict(status='PLAN_ONLY_NO_FETCH_NO_PATIENT_QUERY',requests=[url(kind,gene) for gene in profile['symbols'] for kind in ['genes','variants']])
 if not digest or sha(path)!=digest or profile.get('review_status')!='ROOT_REVIEWED_PUBLIC_LOVD_BASIC_REFERENCE_PROFILE':raise ValueError('Exact ROOT-reviewed public reference profile required')
 if profile['code_sha256']!=sha(Path(__file__)) or profile['panel']['sha256']!=PANEL_SHA or pin(profile['panel']['path'])!=profile['panel'] or len(profile['symbols'])!=39 or profile['symbols']!=[r['symbol'] for r in load(profile['panel']['path'])['genes']]:raise ValueError('Reviewed code and fixed public panel changed')
 if profile.get('patient_variants_allowed') is not False or profile.get('case_level_json_allowed') is not False or profile.get('no_source_vcf_mutation') is not True or profile['base_url']!=BASE or not 0<profile['requests_per_second']<=4 or profile['max_requests']!=78 or not 0<profile['max_response_bytes']<=16*1024**2 or not 0<profile['max_total_bytes']<=256*1024**2 or not 0<profile['max_minutes']<=10:raise ValueError('Reviewed public-only source and request/byte limits required')
 end=min(END,dt.datetime.fromisoformat(profile['deadline_utc']),now()+dt.timedelta(minutes=profile['max_minutes']))
 if end<=now():raise TimeoutError('Original compute deadline reached; do not extend')
 output=safe(profile['output'])
 if output.exists():raise ValueError('Fresh source snapshot required; no overwrite/resume across releases')
 output.mkdir(mode=0o700,parents=True);receipts=[];status_counts=Counter();total=0;request_count=0
 for symbol in profile['symbols']:
  gene_metadata=None
  for kind in ['genes','variants']:
   if now()>=end:raise TimeoutError('Reviewed source collection deadline reached; incomplete snapshot retained')
   if shutil.disk_usage(output).free<150*1024**3:raise OSError('Reserve150GiB guard; source collection stopped')
   start=time.monotonic();body,receipt=fetch(url(kind,symbol),min(profile['max_response_bytes'],profile['max_total_bytes']-total),end);request_count+=1;receipt.update(gene=symbol,kind=kind)
   if body is not None:
    total+=len(body);target=output/(symbol+'.'+kind+'.json');target.write_bytes(body);receipt['path']=str(target)
    try:parsed=rows(body,kind,symbol)
    except ValueError:
     receipt['application_status']='SCHEMA_REJECTED_NOT_APPLIED';dump(output/(symbol+'.'+kind+'.receipt.json'),receipt);raise
    receipt['record_count']=len(parsed)
    if kind=='genes':
     if len(parsed)!=1:raise ValueError('Expected one public gene metadata record')
     gene_metadata=parsed[0];receipt['source_build']=gene_metadata['refseq_build']
    else:
     counts=Counter(mapping_status(row,gene_metadata) for row in parsed);receipt['mapping_status_counts']=dict(counts);status_counts.update(counts)
    receipt['application_status']='REFERENCE_RETRIEVED_NOT_APPLIED'
   receipts.append(receipt);dump(output/(symbol+'.'+kind+'.receipt.json'),receipt)
   wait=1/profile['requests_per_second']-(time.monotonic()-start)
   if wait>0:time.sleep(min(wait,max(0,(end-now()).total_seconds())))
 manifest=dict(schema_version=1,status='PUBLIC_LOVD_BASIC_REFERENCE_SNAPSHOT_NOT_APPLIED',coverage='BASIC_API_RESPONSES_ONLY_COMPLETENESS_NOT_ESTABLISHED',profile=pin(path),code_sha256=profile['code_sha256'],public_panel=profile['panel'],requests=request_count,bytes=total,source_receipts=receipts,mapping_status_counts=dict(status_counts),matched_patient_variants=0,patient_variant_queries=0,source_variants_queried_only_by_fixed_public_gene=True,case_level_json_retrieved=False,existing_vcf_changed=False,clinical_classification=False,source_rights='Per-source/per-submitter rights retained; no blanket data license inferred from software/manual license',created_utc=now().isoformat());dump(output/'manifest.json',manifest);return manifest
def control(output):
 output=safe(output)
 if output.exists():raise ValueError('Fresh public control directory required')
 output.mkdir(mode=0o700,parents=True);body,receipt=fetch(CONTROL,1024**2,min(END,now()+dt.timedelta(minutes=1)))
 if body is not None:
  (output/'IVD-0000000563.verbatim.json').write_bytes(body);parsed=rows(body,'variants','IVD');receipt.update(schema_fields=sorted({key for row in parsed for key in row}),rows=len(parsed),mapping_status_counts=dict(Counter(mapping_status(row,None) for row in parsed)),application_status='PUBLIC_CONTROL_ONLY_NOT_APPLIED',manual_example=True)
 dump(output/'receipt.json',receipt);return receipt
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True);a=sub.add_parser('prepare');a.add_argument('--panel',type=Path,required=True);a.add_argument('--profile',type=Path,required=True);a.add_argument('--output',type=Path,required=True);a=sub.add_parser('control');a.add_argument('--output',type=Path,required=True);a=sub.add_parser('collect');a.add_argument('--profile',type=Path,required=True);a.add_argument('--reviewed-profile-sha256');a.add_argument('--execute',action='store_true');a=p.parse_args();result=prepare(a.panel,a.profile,a.output) if a.command=='prepare' else control(a.output) if a.command=='control' else collect(a.profile,a.reviewed_profile_sha256,a.execute);print(json.dumps(result,indent=2))
