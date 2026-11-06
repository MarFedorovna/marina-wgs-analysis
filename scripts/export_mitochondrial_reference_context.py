#!/usr/bin/env python3
"""Independent MT-only wrapper of preserved, reviewed local reference export."""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import signal
from types import SimpleNamespace
import uuid

LEGACY_SHA='4fc954c6cacce92a199fa477fcb75a3c662ecb1fcd0b0dce0ef65b6e13301832'
path=Path(__file__).resolve().parent/'export_population_mt_reference_context.py'
if path.is_symlink() or hashlib.sha256(path.read_bytes()).hexdigest()!=LEGACY_SHA:raise ValueError('Exact preserved MT exporter source required')
spec=importlib.util.spec_from_file_location('preserved_mt_reference_export',path);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
END=old.ABSOLUTE_DEADLINE
MT_SHA='f069164654316043d8a96f47ad5b29413d7a28b0b9ba320648206110e42739b0'
def atomic(path,value):
 path=old.safe(path);temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(value,indent=2)+'\n');os.replace(temporary,path)
def prepare(population_rule,rule_path,output,state):
 if Path(rule_path).exists() or Path(output).exists():raise ValueError('Fresh candidate and independent output required')
 source=old.load(population_rule)
 if source.get('review_status')!='ROOT_REVIEWED_EXACT_LOCAL_ONLY_POPULATION_RULE' or source.get('no_network') is not True:raise ValueError('Actual ROOT-reviewed original MT source profile required')
 ready=old.load(old.verify(source['resource_ready_manifest']))
 if source['resource_ready_manifest']['sha256']!=old.EXPECTED_READY_SHA:raise ValueError('Exact public resource release required')
 database=next(p for p in ready['artifacts'] if Path(p['path']).name=='public-mt-population-v2.sqlite');database={k:database[k] for k in ['path','bytes','sha256']}
 if database['sha256']!=MT_SHA:raise ValueError('Exact verified MT reference database required')
 old.verify(database);old.verify(source['mt_audit'])
 if old.load(source['mt_audit']['path']).get('status')!='PASS_RESEARCH_WITH_LIMITATIONS':raise ValueError('Actual MT-v2 technical audit required')
 entries=[]
 for value in source['mt_sources']:
  entry=dict(value);receipt=old.load(old.verify(entry['receipt']))
  if receipt.get('status')!='PASS' or not receipt.get('identity'):raise ValueError('Actual successful MT-v2 calls required')
  for item in [entry['vcf'],*entry['qc']]:old.verify(item)
  for item in [entry['vcf'],*entry['qc']]:
   relative=Path(item['path']).relative_to(Path(entry['receipt']['path']).parent).as_posix();expected=receipt['artifacts'].get(relative,{})
   if any(expected.get(k)!=item[k] for k in ['bytes','sha256']):raise ValueError('Original MT call/QC artifact differs from actual source receipt')
  index=[]
  for suffix in ['.tbi','.csi']:
   p=Path(entry['vcf']['path']+suffix)
   if p.exists() and p.name in receipt['artifacts']:
    item=old.pin(p);expected=receipt['artifacts'][p.name]
    if any(item[k]!=expected[k] for k in ['bytes','sha256']):raise ValueError('MT source index receipt differs')
    index.append(item)
  if not index:raise ValueError('Receipt-bound original MT VCF index required')
  entry['indices']=index;entries.append(entry)
 if len(entries)!=4 or len({p['sample'] for p in entries})!=4:raise ValueError('Exactly four original MT sources required')
 scripts=Path(__file__).resolve().parent;code={str(p):old.sha(p) for p in [Path(__file__).resolve(),scripts/'export_population_mt_reference_context.py']};rule=dict(schema_version=1,review_status='PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING',source_population_rule=old.pin(population_rule),source_resource_manifest=source['resource_ready_manifest'],source_database=database,mt_audit=source['mt_audit'],mt_sources=entries,expected_samples=source['expected_samples'],code_sha256=code,output=str(output),state_directory=str(state),deadline_utc=END.isoformat(),max_hours=1,cpu_limit=1,memory_limit_bytes=1024**3,max_output_bytes=512*1024**2,min_free_bytes=150*1024**3,no_network=True,no_nuclear_dependency=True,no_primary_or_mt_source_mutation=True)
 Path(rule_path).parent.mkdir(parents=True,exist_ok=True);old.dump(rule_path,rule);return dict(status=rule['review_status'],rule=str(rule_path),rule_sha256=old.sha(rule_path),source_count=4,controller_started=False)
def validate(rule,path,digest,current=None):
 if not digest or old.sha(path)!=digest or rule.get('review_status')!='ROOT_REVIEWED_EXACT_INDEPENDENT_MT_REFERENCE_RULE':raise ValueError('Exact ROOT-reviewed independent MT rule required')
 if rule.get('schema_version')!=1 or rule.get('no_network') is not True or rule.get('no_nuclear_dependency') is not True or rule.get('no_primary_or_mt_source_mutation') is not True or rule['cpu_limit']!=1 or rule['memory_limit_bytes']!=1024**3 or not 0<rule['max_output_bytes']<=512*1024**2 or rule['min_free_bytes']<150*1024**3 or not 0<rule['max_hours']<=1:raise ValueError('Bounded independent local-only MT contract required')
 deadline=dt.datetime.fromisoformat(rule['deadline_utc']);current=current or old.now()
 if deadline>END or deadline<=current:raise ValueError('Original VM deadline cannot extend')
 if Path(rule['output'])!=Path('/srv/marina/outputs/mitochondrial-reference-context-v1') or Path(rule['state_directory'])!=Path('/srv/marina/private/mitochondrial-reference-context-state-v1') or Path(rule['output']).exists():raise ValueError('Separate new fixed MT output/private state required')
 expected={str(Path(__file__).resolve()),str(Path(__file__).resolve().parent/'export_population_mt_reference_context.py')}
 if set(rule['code_sha256'])!=expected:raise ValueError('Only exact wrapper/preserved exporter pair allowed')
 for p,d in rule['code_sha256'].items():
  if old.sha(p)!=d:raise ValueError('Reviewed MT code changed')
 for item in [rule['source_population_rule'],rule['source_resource_manifest'],rule['source_database'],rule['mt_audit']]:old.verify(item)
 if rule['source_database']['sha256']!=MT_SHA or rule['source_resource_manifest']['sha256']!=old.EXPECTED_READY_SHA:raise ValueError('Unreviewed reference release')
 if len(rule['mt_sources'])!=4 or sorted(p['sample'] for p in rule['mt_sources'])!=sorted(rule['expected_samples']):raise ValueError('Four original samples required')
 for entry in rule['mt_sources']:
  if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]{0,99}',entry['sample']):raise ValueError('Unsafe sample filename')
  for item in [entry['receipt'],entry['vcf'],*entry['qc'],*entry['indices']]:old.verify(item)
 return min(deadline,current+dt.timedelta(hours=rule['max_hours']))
def verify_rows(entry,table):
 counts={};rows=0
 with gzip.open(entry['vcf']['path'],'rt') as original,gzip.open(table,'rt',newline='') as annotated:
  source=(line.rstrip('\n').split('\t') for line in original if not line.startswith('#'));reader=csv.DictReader(annotated,delimiter='\t')
  for call,row in itertools.zip_longest(source,reader):
   if call is None or row is None:raise ValueError('MT export lost/added an original call')
   if [row[f] for f in ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','raw_INFO','raw_FORMAT','raw_sample']]!=call:raise ValueError('Original MT call/quality fields changed')
   if row['sample']!=entry['sample']:raise ValueError('MT sample identity changed')
   data=dict(zip(call[8].split(':'),call[9].split(':')))
   if any(row[f]!=data.get(f,'.') for f in ['GT','AD','AF','DP','FT']):raise ValueError('Original MT sample fields changed')
   for field in ['gnomAD_MT_match_status','Helix_match_status']:counts[field+':'+row[field]]=counts.get(field+':'+row[field],0)+1
   if ',' in call[4] and (row['gnomAD_MT_reference_json'] or row['Helix_reference_json']):raise ValueError('Multiallelic MT source was flattened/joined')
   if row['Helix_reference_json'] and json.loads(row['Helix_reference_json'])['AN'] is not None:raise ValueError('Helix source AN must remain explicitly missing')
   rows+=1
 return dict(rows=rows,status_counts=counts,all_original_calls_fields_order_preserved=True)
def export(rule,deadline=None):
 deadline=min(deadline or END,END);output=old.safe(rule['output']);state=old.safe(rule['state_directory']);state.mkdir(parents=True,exist_ok=True)
 if output.exists():raise ValueError('Existing output cannot be replaced')
 if shutil.disk_usage(state).free<rule['min_free_bytes']:raise ValueError('Required disk reserve reached')
 for item in [rule['source_database'],rule['mt_audit'],rule['source_resource_manifest']]:old.verify(item)
 stage=output.parent/('.'+output.name+'.'+uuid.uuid4().hex);stage.mkdir(mode=0o700);connection=old.connection(rule['source_database']['path']);source=SimpleNamespace(mt=connection)
 # This isolated module instance adds only the independent reserve check to
 # the preserved export guard; the scientific join/call handling is unchanged.
 original_guard=old.output_guard
 def bounded(stage_,cap,end):
  original_guard(stage_,cap,end)
  if shutil.disk_usage(state).free<rule['min_free_bytes']:raise ValueError('Required disk reserve reached')
 old.output_guard=bounded
 try:
  summary=old.mt_exports(rule['mt_sources'],stage,source,deadline,rule['max_output_bytes']);verified=[]
  for entry in rule['mt_sources']:
   verified.append(dict(sample=entry['sample'],**verify_rows(entry,stage/'samples'/entry['sample']/'MT-population-reference.tsv.gz')))
   for item in [entry['vcf'],*entry['indices']]:
    destination=stage/'samples'/entry['sample']/Path(item['path']).name;shutil.copyfile(old.verify(item),destination)
    if old.sha(destination)!=item['sha256']:raise ValueError('Copied original MT VCF/index changed')
    bounded(stage,rule['max_output_bytes'],deadline)
  page='<!doctype html><meta charset="utf-8"><title>Separate mitochondrial reference context</title><h1>Mitochondrial calls and population context</h1><p>Original filtered calls, raw heteroplasmy/quality and NUMT/haplocheck QC retained. Missing population frequency is not zero; multiallelic sources are not split. Research reference context; no automatic clinical class or limit-of-detection claim.</p><ul>'
  for entry in rule['mt_sources']:page+='<li><a href="samples/'+entry['sample']+'/MT-population-reference.tsv.gz">'+entry['sample']+' reference context</a></li>'
  (stage/'index.html').write_text(page+'</ul>');applied=[]
  by_sample={entry['sample']:entry for entry in rule['mt_sources']}
  for path in stage.rglob('*'):
   if not path.is_file():continue
   relative=path.relative_to(stage).as_posix();role='MT_SOURCE_QC' if '/MT-QC/' in relative else 'MT_SOURCE_VCF' if '.vcf.gz' in path.name else 'MT_POPULATION_CONTEXT';sample=path.relative_to(stage).parts[1] if relative.startswith('samples/') else None;entry=by_sample.get(sample);item=dict(path=str(output/relative),relative_path='applied-evidence/MT/'+relative,bytes=path.stat().st_size,sha256=old.sha(path),role=role,source_receipts=[entry['receipt']] if entry else [e['receipt'] for e in rule['mt_sources']])
   if role=='MT_SOURCE_QC':item['source_original_pin']=entry['qc'][int(path.name.split('-',1)[0])]
   elif role=='MT_SOURCE_VCF':item['source_original_pin']=next(p for p in [entry['vcf'],*entry['indices']] if Path(p['path']).name==path.name)
   applied.append(item)
  manifest=dict(schema_version=1,status='PASS_INDEPENDENT_MITOCHONDRIAL_REFERENCE_CONTEXT',source_population_rule=rule['source_population_rule'],resource_ready_manifest=rule['source_resource_manifest'],source_versions=old.load(rule['source_resource_manifest']['path'])['source_versions'],source_database=rule['source_database'],mt_audit=rule['mt_audit'],source_mt_inputs=rule['mt_sources'],source_mt_receipts=[p['receipt'] for p in rule['mt_sources']],source_mt_vcfs=[p['vcf'] for p in rule['mt_sources']],samples=rule['expected_samples'],summary=summary,verified=verified,applied_artifacts=applied,producer_code_sha256=rule['code_sha256'],network_requests=0,no_nuclear_dependency=True,primary_and_original_MT_changed=False,no_automatic_ACMG=True,clinical_validation=False,created_utc=old.now().isoformat());old.dump(stage/'manifest.json',manifest);artifacts={p.relative_to(stage).as_posix():dict(bytes=p.stat().st_size,sha256=old.sha(p)) for p in stage.rglob('*') if p.is_file()};identity=hashlib.sha256(json.dumps(dict(resources=rule['source_resource_manifest'],sources=rule['mt_sources'],code=rule['code_sha256']),sort_keys=True).encode()).hexdigest();old.dump(stage/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='mitochondrial-reference-context-v1',identity=identity,artifacts=artifacts,metadata=dict(no_nuclear_dependency=True,primary_and_original_MT_changed=False)));bounded(stage,rule['max_output_bytes'],deadline);old.promote(stage,output);return manifest
 finally:connection.close();old.output_guard=original_guard
def run(path,digest,execute=False):
 rule=old.load(path);state=Path(rule['state_directory'])
 if not execute:return old.load(state/'status.json') if (state/'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
 end=validate(rule,path,digest);os.umask(0o077);state.mkdir(parents=True,exist_ok=True)
 import fcntl
 with old.safe(state/'controller.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  def status(label,**extra):
   value=dict(status=label,pid=os.getpid(),rule_sha256=digest,deadline_utc=end.isoformat(),updated_utc=old.now().isoformat(),**extra);atomic(state/'status.json',value);return value
  def alarm(signum,frame):raise TimeoutError('Original VM/MT job deadline reached')
  signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.001,(end-old.now()).total_seconds()))
  try:
   import resource
   resource.setrlimit(resource.RLIMIT_AS,(rule['memory_limit_bytes'],rule['memory_limit_bytes']))
   if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
   status('EXPORTING_INDEPENDENT_MT_REFERENCE_CONTEXT');result=export(rule,end);return status('PASS_INDEPENDENT_MITOCHONDRIAL_REFERENCE_CONTEXT',result=result)
  except BaseException as error:status('FAILED_INDEPENDENT_MT_REFERENCE_CONTEXT',error=str(error)[:400]);raise
  finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True);a=sub.add_parser('prepare')
 for name in ['population-rule','rule','output','state-directory']:a.add_argument('--'+name,type=Path,required=True)
 a=sub.add_parser('run');a.add_argument('--rule',type=Path,required=True);a.add_argument('--reviewed-rule-sha256');a.add_argument('--execute',action='store_true');a=p.parse_args();print(json.dumps(prepare(a.population_rule,a.rule,a.output,a.state_directory) if a.command=='prepare' else run(a.rule,a.reviewed_rule_sha256,a.execute),indent=2))
