#!/usr/bin/env python3
"""Prepare-only by default; ROOT-reviewed bounded continuation to one combined master."""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
READY_SHA='465ee667cb78baddc3665d06adb15fb27ad54c47415ff14238fda5f0f42c3d85'
def now():return dt.datetime.now(dt.timezone.utc)
def safe(path):
 path=Path(path)
 if any(p.is_symlink() for p in [path,*path.parents]):raise ValueError('Symlink source/state refused')
 return path
def load(path):return json.loads(safe(path).read_text())
def sha(path):
 h=hashlib.sha256()
 with safe(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024**2),b''):h.update(block)
 return h.hexdigest()
def pin(path):
 path=safe(Path(path).absolute());return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
def verify(item):
 path=safe(item['path'])
 if not path.is_file() or path.stat().st_size!=item['bytes'] or sha(path)!=item['sha256']:raise ValueError('Pinned source changed: '+str(path))
def prepare(binding_path,tissue_rule_path,population_rule_path,ready_path,panel_path,rule_path,output,state):
 if Path(rule_path).exists() or Path(output).exists():raise ValueError('New private candidate and output required')
 binding=load(binding_path);package=Path(binding['code_package']['path']);package_pin=pin(package)
 if package_pin['sha256']!=binding['code_package']['sha256'] or len(load(package)['files'])!=15:raise ValueError('Unchanged sealed primary 15-file package required')
 primary=[pin(binding_path),package_pin,pin(binding['annotation_config']['path'])]
 if primary[-1]['sha256']!=binding['annotation_config']['sha256']:raise ValueError('Primary config differs')
 for item in load(package)['files']:
  value=pin(package.parent/item['path'])
  if value['bytes']!=item['bytes'] or value['sha256']!=item['sha256']:raise ValueError('Sealed primary code changed')
  primary.append(value)
 accepted=pin(Path(binding['accepted_filter']['path'])/'workflow-receipt.json')
 if accepted['sha256']!=binding['accepted_filter']['receipt_sha256']:raise ValueError('Accepted original stage5 differs')
 primary.append(accepted);tissue=load(tissue_rule_path);population=load(population_rule_path)
 if tissue['workflow_root']!=binding['output_root'] or population['workflow_root']!=binding['output_root']:raise ValueError('Producer workflow roots differ')
 producer=[]
 for rule in [tissue,population]:
  for path,digest in rule['code_sha256'].items():
   item=pin(path)
   if item['sha256']!=digest:raise ValueError('Producer code not exact reviewed candidate')
   producer.append(item)
 ready=pin(ready_path)
 if ready['sha256']!=READY_SHA or population['resource_ready_manifest']!=ready:raise ValueError('Exact population resource release required')
 context=pin(tissue['context_index'])
 if context['sha256']!=tissue['context_index_sha256']:raise ValueError('Tissue context pin differs')
 scripts=Path(__file__).resolve().parent;code={str(path):sha(path) for path in [Path(__file__).resolve(),scripts/'join_combined_annotation_master.py']}
 value=dict(schema_version=1,review_status='PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING',reviewed_utc=None,workflow_root=binding['output_root'],output=str(output),state_directory=str(state),primary_source_pins=primary,producer_code_pins=producer,code_sha256=code,resource_ready_manifest=ready,public_panel=pin(panel_path),tissue_context_index=context,tissue_root=tissue['output'],population_root=population['output'],producer_gates=[dict(kind='tissue',path=str(Path(tissue['state_directory'])/'status.json'),success_status='PASS_SEPARATE_CONTEXT_EXPORT_PRIMARY_UNCHANGED'),dict(kind='population',path=str(Path(population['state_directory'])/'status.json'),success_status='PASS_LOCAL_POPULATION_MT_CONTEXT_PRIMARY_UNCHANGED')],deadline_utc=DEADLINE.isoformat(),max_hours=8,interval_seconds=30,cpu_limit=2,memory_limit_bytes=4*1024**3,max_sqlite_bytes=8*1024**3,max_output_bytes=6*1024**3,min_free_bytes=150*1024**3,no_network=True,no_primary_or_views_mutation=True)
 Path(rule_path).parent.mkdir(parents=True,exist_ok=True)
 with safe(rule_path).open('x') as stream:json.dump(value,stream,indent=2);stream.write('\n')
 return dict(status=value['review_status'],rule=str(rule_path),rule_sha256=sha(rule_path),no_controller_started=True,no_primary_or_views_change=True)
def validate(rule,path,reviewed_sha,current=None):
 if not reviewed_sha or sha(path)!=reviewed_sha or rule.get('review_status')!='ROOT_REVIEWED_EXACT_COMBINED_MASTER_RULE':raise ValueError('Exact ROOT-reviewed combined rule required')
 if rule.get('schema_version')!=1 or rule.get('no_network') is not True or rule.get('no_primary_or_views_mutation') is not True:raise ValueError('Local-only immutable-primary/views contract required')
 if rule['cpu_limit']!=2 or rule['memory_limit_bytes']!=4*1024**3 or not 0<rule['max_hours']<=8 or not 10<=rule['interval_seconds']<=60 or not 0<rule['max_sqlite_bytes']<=8*1024**3 or not 0<rule['max_output_bytes']<=6*1024**3 or rule['min_free_bytes']<150*1024**3:raise ValueError('Reviewed CPU/memory/disk/wait bounds required')
 current=current or now();deadline=dt.datetime.fromisoformat(rule['deadline_utc'])
 if deadline>DEADLINE or deadline<=current:raise ValueError('Original budget deadline cannot extend')
 scripts=Path(__file__).resolve().parent
 if set(rule['code_sha256'])!={str(Path(__file__).resolve()),str(scripts/'join_combined_annotation_master.py')}:raise ValueError('Only exact controller/helper pair allowed')
 for path,digest in rule['code_sha256'].items():
  if sha(path)!=digest:raise ValueError('Reviewed combined code changed')
 for item in [*rule['primary_source_pins'],*rule['producer_code_pins'],rule['resource_ready_manifest'],rule['public_panel'],rule['tissue_context_index']]:verify(item)
 if rule['resource_ready_manifest']['sha256']!=READY_SHA:raise ValueError('Unreviewed resource ready manifest')
 workflow=Path(rule['workflow_root']).resolve();output=safe(rule['output']).resolve();state=safe(rule['state_directory']).resolve()
 if output.exists() or output==workflow or workflow in output.parents or state==workflow or workflow in state.parents or output==state or output in state.parents or state in output.parents:raise ValueError('Separate new output and private state required')
 if output.name!='combined-annotation-master-v1' or output.parent!=workflow.parent:raise ValueError('Exact new combined output root required')
 return min(deadline,current+dt.timedelta(hours=rule['max_hours']))
def readiness(rule):
 for name in ['06_annotate','07_report']:
  path=Path(rule['workflow_root'])/'stages'/name/'workflow-receipt.json'
  if not path.exists():return False
  value=load(path)
  if str(value.get('status','')).startswith('FAIL'):raise ValueError('Primary failed; combined controller cannot restart it')
  if value.get('status')!='PASS' or not value.get('identity'):return False
 for gate in rule['producer_gates']:
  path=Path(gate['path'])
  if not path.exists():return False
  value=load(path);status=str(value.get('status',''))
  if status.startswith(('FAIL','STOPPED')):raise ValueError('Required actual reference enrichment failed/stopped: '+gate['kind'])
  if status!=gate['success_status'] or not value.get('rule_sha256'):return False
 return True
def run(path,reviewed_sha,execute=False):
 rule=load(path);state=Path(rule['state_directory'])
 if not execute:return load(state/'status.json') if (state/'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
 deadline=validate(rule,path,reviewed_sha);os.umask(0o077);state.mkdir(mode=0o700,parents=True,exist_ok=True)
 import fcntl
 with safe(state/'controller.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  def status(label,**extra):
   value=dict(status=label,pid=os.getpid(),rule_sha256=reviewed_sha,deadline_utc=deadline.isoformat(),updated_utc=now().isoformat(),**extra);target=safe(state/'status.json');temporary=safe(state/'status.tmp');temporary.write_text(json.dumps(value,indent=2)+'\n');os.replace(temporary,target);return value
  try:
   while not readiness(rule):
    if now()>=deadline:return status('STOPPED_ORIGINAL_DEADLINE_NO_PRIMARY_ACTUATION')
    status('WAITING_PRIMARY_TISSUE_POPULATION_PASS');time.sleep(min(rule['interval_seconds'],(deadline-now()).total_seconds()))
   validate(rule,path,reviewed_sha);status('JOINING_ONE_ROW_PRESERVING_MASTER')
   helper=Path(__file__).resolve().parent/'join_combined_annotation_master.py';spec=importlib.util.spec_from_file_location('reviewed_combined_master_helper',helper);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);result=module.join(rule,deadline)
   return status('PASS_COMBINED_MASTER_READY_FOR_VIEWS_V2',result=result)
  except BaseException as error:status('FAILED_COMBINED_NO_PRIMARY_OR_VIEWS_RETRY',error=str(error)[:400]);raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
 a=sub.add_parser('prepare')
 for name in ['primary-binding','tissue-rule','population-rule','ready-manifest','panel','rule','output','state-directory']:a.add_argument('--'+name,type=Path,required=True)
 a=sub.add_parser('run');a.add_argument('--rule',type=Path,required=True);a.add_argument('--reviewed-rule-sha256');a.add_argument('--execute',action='store_true');args=p.parse_args()
 result=prepare(args.primary_binding,args.tissue_rule,args.population_rule,args.ready_manifest,args.panel,args.rule,args.output,args.state_directory) if args.command=='prepare' else run(args.rule,args.reviewed_rule_sha256,args.execute);print(json.dumps(result,indent=2))
