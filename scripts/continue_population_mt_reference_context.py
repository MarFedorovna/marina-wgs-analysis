#!/usr/bin/env python3
"""Prepare a private review candidate; run only an exactly ROOT-reviewed local-only rule."""
import argparse
import datetime as dt
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import time

ABSOLUTE_DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
EXPECTED_READY_SHA='465ee667cb78baddc3665d06adb15fb27ad54c47415ff14238fda5f0f42c3d85'
def now():return dt.datetime.now(dt.timezone.utc)
def safe(path):
 path=Path(path)
 if any(p.is_symlink() for p in [path,*path.parents]):raise ValueError('Symlink source/state rejected')
 return path
def sha(path):
 digest=hashlib.sha256()
 with safe(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024**2),b''):digest.update(block)
 return digest.hexdigest()
def load(path):return json.loads(safe(path).read_text())
def pin(path):
 path=safe(Path(path).absolute());return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
def verify(item):
 path=safe(item['path'])
 if not path.is_file() or path.stat().st_size!=item['bytes'] or sha(path)!=item['sha256']:raise ValueError('Reviewed pin changed: '+str(path))
def save_new(path,value):
 with safe(path).open('x') as stream:json.dump(value,stream,indent=2);stream.write('\n')
def prepare(binding_path,ready_path,coverage_root,panel_path,mt_root,mt_audit,rule_path,output,state):
 if Path(rule_path).exists() or Path(output).exists():raise ValueError('Choose a new private rule and a new sidecar output')
 ready_pin=pin(ready_path)
 if ready_pin['sha256']!=EXPECTED_READY_SHA:raise ValueError('Exact ready manifest 465ee667 required')
 binding=load(binding_path);package=Path(binding['code_package']['path']);package_pin=pin(package)
 if package_pin['sha256']!=binding['code_package']['sha256']:raise ValueError('Primary sealed package changed')
 files=load(package)['files']
 if len(files)!=15:raise ValueError('Expected unchanged sealed 15-file primary runtime')
 primary_pins=[pin(binding_path),package_pin,pin(binding['annotation_config']['path'])]
 if primary_pins[-1]['sha256']!=binding['annotation_config']['sha256']:raise ValueError('Primary annotation config pin differs')
 for item in files:
  actual=pin(package.parent/item['path'])
  if actual['bytes']!=item['bytes'] or actual['sha256']!=item['sha256']:raise ValueError('Sealed primary code differs')
  primary_pins.append(actual)
 accepted=pin(Path(binding['accepted_filter']['path'])/'workflow-receipt.json')
 if accepted['sha256']!=binding['accepted_filter']['receipt_sha256']:raise ValueError('Accepted original stage5 receipt changed')
 primary_pins.append(accepted)
 panel_pin=pin(panel_path);panel=load(panel_path)
 if panel.get('no_patient_selection') is not True or len(panel['genes'])!=39:raise ValueError('Fixed public39 panel required')
 coverage=[]
 for gene in panel['genes']:
  if gene['symbol'].startswith('MT-'):continue
  source=Path(coverage_root)/(gene['symbol']+'.json');receipt=load(Path(coverage_root)/(gene['symbol']+'-receipt.json'));item=pin(source)
  if item['sha256']!=receipt['sha256'] or item['bytes']!=receipt['bytes']:raise ValueError('Public coverage receipt changed')
  value=load(source)['data']['gene']
  if str(value['hgnc_id'])!=gene['hgnc_id'] or value['reference_genome']!='GRCh38':raise ValueError('Coverage HGNC/build differs')
  coverage.append(dict(item,hgnc_id=gene['hgnc_id'],gene=gene['symbol']))
 audit_pin=pin(mt_audit)
 if load(mt_audit).get('status')!='PASS_RESEARCH_WITH_LIMITATIONS':raise ValueError('MT-v2 technical audit required')
 mt=[]
 for receipt_path in sorted(Path(mt_root).glob('checkpoints/*/calls/workflow-receipt.json')):
  receipt=load(receipt_path)
  if receipt.get('status')!='PASS':raise ValueError('MT source not complete')
  directory=receipt_path.parent;vcf=directory/'mitochondria.filtered.vcf.gz';vcf_pin=pin(vcf)
  expected=receipt['artifacts'][vcf.name]
  if vcf_pin['sha256']!=expected['sha256'] or vcf_pin['bytes']!=expected['bytes']:raise ValueError('Final MT filtered VCF receipt mismatch')
  with gzip.open(vcf,'rt') as stream:
   for line in stream:
    if line.startswith('#CHROM'):samples=line.rstrip('\n').split('\t')[9:];break
   else:raise ValueError('MT sample header missing')
  if len(samples)!=1 or not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]{0,99}',samples[0]):raise ValueError('MT sample header/filename contract')
  qc=[]
  for name in ['research-report.json','numt-assessment.json','haplo-output/output_json','contamination.json']:
   item=pin(directory/name);expected=receipt['artifacts'][name]
   if item['sha256']!=expected['sha256'] or item['bytes']!=expected['bytes']:raise ValueError('MT source QC receipt mismatch')
   qc.append(item)
  mt.append(dict(sample=samples[0],receipt=pin(receipt_path),vcf=vcf_pin,qc=qc))
 if len(mt)!=4 or len({v['sample'] for v in mt})!=4:raise ValueError('Expected four successful MT-v2 sample sources')
 scripts=Path(__file__).resolve().parent;code={str(p):sha(p) for p in [Path(__file__).resolve(),scripts/'export_population_mt_reference_context.py']}
 rule=dict(schema_version=1,review_status='PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING',reviewed_utc=None,code_sha256=code,primary_source_pins=primary_pins,resource_ready_manifest=ready_pin,public_panel=panel_pin,coverage_sources=coverage,mt_sources=mt,mt_audit=audit_pin,expected_samples=sorted(v['sample'] for v in mt),workflow_root=binding['output_root'],output=str(output),state_directory=str(state),deadline_utc=ABSOLUTE_DEADLINE.isoformat(),max_hours=8,interval_seconds=30,max_output_bytes=2*1024**3,no_network=True,no_primary_workflow_mutation=True)
 Path(rule_path).parent.mkdir(parents=True,exist_ok=True);save_new(rule_path,rule)
 return dict(status=rule['review_status'],rule=str(rule_path),rule_sha256=sha(rule_path),MT_sources=4,coverage_genes=len(coverage),no_controller_started=True,no_approval_created=True)
def validate(rule,path,reviewed_sha,current=None):
 if not reviewed_sha or sha(path)!=reviewed_sha or rule.get('review_status')!='ROOT_REVIEWED_EXACT_LOCAL_ONLY_POPULATION_RULE':raise ValueError('Exact ROOT-reviewed rule SHA/status required')
 if rule.get('schema_version')!=1 or rule.get('no_network') is not True or rule.get('no_primary_workflow_mutation') is not True:raise ValueError('Local-only no-primary-mutation contract required')
 if not 0<rule['max_hours']<=8 or not 10<=rule['interval_seconds']<=60 or not 0<rule['max_output_bytes']<=2*1024**3:raise ValueError('Bounded controller required')
 current=current or now();deadline=dt.datetime.fromisoformat(rule['deadline_utc'])
 if deadline>ABSOLUTE_DEADLINE or deadline<=current:raise ValueError('Original absolute deadline exceeded')
 scripts=Path(__file__).resolve().parent;expected={str(Path(__file__).resolve()),str(scripts/'export_population_mt_reference_context.py')}
 if set(rule['code_sha256'])!=expected:raise ValueError('Only this exact controller/exporter pair accepted')
 for name,digest in rule['code_sha256'].items():
  if sha(name)!=digest:raise ValueError('Reviewed code changed')
 for item in [*rule['primary_source_pins'],rule['resource_ready_manifest'],rule['public_panel'],rule['mt_audit'],*rule['coverage_sources']]:verify(item)
 if rule['resource_ready_manifest']['sha256']!=EXPECTED_READY_SHA:raise ValueError('Unreviewed resource manifest')
 for source in rule['mt_sources']:
  for item in [source['receipt'],source['vcf'],*source['qc']]:verify(item)
 if len(rule['mt_sources'])!=4 or {v['sample'] for v in rule['mt_sources']}!=set(rule['expected_samples']):raise ValueError('MT and primary reviewed sample set must match')
 workflow=Path(rule['workflow_root']).resolve();output=safe(rule['output']).resolve();state=safe(rule['state_directory']).resolve()
 if output.exists() or output==workflow or workflow in output.parents or state==workflow or workflow in state.parents or state==output or output in state.parents or state in output.parents:raise ValueError('Separate immutable sidecar output/private state required')
 return min(deadline,current+dt.timedelta(hours=rule['max_hours']))
def readiness(workflow):
 receipts=[]
 for name in ['06_annotate','07_report']:
  path=Path(workflow)/'stages'/name/'workflow-receipt.json'
  if not path.exists():return False
  value=load(path)
  if str(value.get('status','')).startswith('FAIL'):raise ValueError('Primary failed; this controller never restarts it')
  receipts.append(value)
 return all(v.get('status')=='PASS' and v.get('identity') for v in receipts)
def run(path,reviewed_sha,execute=False):
 rule=load(path);state=Path(rule['state_directory'])
 if not execute:return load(state/'status.json') if (state/'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
 deadline=validate(rule,path,reviewed_sha);os.umask(0o077);state.mkdir(mode=0o700,parents=True,exist_ok=True)
 import fcntl
 with safe(state/'controller.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  def status(label,**extra):
   value=dict(status=label,pid=os.getpid(),rule_sha256=reviewed_sha,updated_utc=now().isoformat(),deadline_utc=deadline.isoformat(),**extra);target=safe(state/'status.json');temp=safe(state/'status.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');os.replace(temp,target);return value
  try:
   while not readiness(rule['workflow_root']):
    if now()>=deadline:return status('STOPPED_ORIGINAL_DEADLINE_NO_PRIMARY_ACTUATION')
    status('WAITING_06_07_PASS_NO_PRIMARY_ACTUATION');time.sleep(min(rule['interval_seconds'],(deadline-now()).total_seconds()))
   validate(rule,path,reviewed_sha);status('EXPORTING_LOCAL_REFERENCE_SIDECARS')
   exporter=Path(__file__).resolve().parent/'export_population_mt_reference_context.py';spec=importlib.util.spec_from_file_location('reviewed_population_exporter',exporter);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
   result=module.export(rule,deadline);return status('PASS_LOCAL_POPULATION_MT_CONTEXT_PRIMARY_UNCHANGED',result=result)
  except BaseException as error:status('FAILED_NO_PRIMARY_RETRY_PARTIAL_RETAINED',error=str(error)[:400]);raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
 a=sub.add_parser('prepare');a.add_argument('--binding',type=Path,required=True);a.add_argument('--ready-manifest',type=Path,required=True);a.add_argument('--coverage-root',type=Path,required=True);a.add_argument('--panel',type=Path,required=True);a.add_argument('--mt-root',type=Path,required=True);a.add_argument('--mt-audit',type=Path,required=True);a.add_argument('--rule',type=Path,required=True);a.add_argument('--output',type=Path,required=True);a.add_argument('--state-directory',type=Path,required=True)
 a=sub.add_parser('run');a.add_argument('--rule',type=Path,required=True);a.add_argument('--reviewed-rule-sha256');a.add_argument('--execute',action='store_true');args=p.parse_args()
 result=prepare(args.binding,args.ready_manifest,args.coverage_root,args.panel,args.mt_root,args.mt_audit,args.rule,args.output,args.state_directory) if args.command=='prepare' else run(args.rule,args.reviewed_rule_sha256,args.execute)
 print(json.dumps(result,indent=2))
