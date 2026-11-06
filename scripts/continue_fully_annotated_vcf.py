#!/usr/bin/env python3
"""Prepare-only controller for a ROOT-reviewed, local evidence embedding job."""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
def now():return dt.datetime.now(dt.timezone.utc)
def safe(path):
 path=Path(path)
 if any(p.is_symlink() for p in [path,*path.parents]):raise ValueError('Symlink path refused')
 return path
def load(path):return json.loads(safe(path).read_text())
def sha(path):
 h=hashlib.sha256()
 with safe(path).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()
def pin(path):
 path=safe(Path(path).absolute());return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
def verify(p):
 q=pin(p['path'])
 if q!=p:raise ValueError('Pinned upstream/controller source changed')
def prepare(combined_rule,python,output,state,rule_path):
 if Path(rule_path).exists() or Path(output).exists():raise ValueError('Fresh candidate/output required')
 source=load(combined_rule)
 if source.get('review_status')!='ROOT_REVIEWED_EXACT_COMBINED_MASTER_RULE':raise ValueError('Actual ROOT-reviewed combined producer rule required')
 pins=[pin(combined_rule),*source['primary_source_pins'],*source['producer_code_pins'],source['resource_ready_manifest'],source['public_panel'],source['tissue_context_index']]
 for path,digest in source['code_sha256'].items():
  item=pin(path)
  if item['sha256']!=digest:raise ValueError('Combined producer changed')
  pins.append(item)
 import subprocess
 result=subprocess.run([str(python),'-c','import pysam,json;from pathlib import Path;p=Path(pysam.__file__).parent;print(json.dumps(dict(version=pysam.__version__,files=[str(f) for f in sorted(p.iterdir()) if f.suffix in [".py",".so"]])))'],capture_output=True,text=True,check=True)
 runtime=json.loads(result.stdout)
 if runtime['version']!='0.24.1':raise ValueError('Existing exact pysam 0.24.1 environment required')
 runtime_pins=[pin(p) for p in runtime['files']]
 code={str(p):sha(p) for p in [Path(__file__).resolve(),Path(__file__).resolve().parent/'embed_vcf_evidence.py']}
 value=dict(schema_version=1,review_status='PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING',workflow_root=source['workflow_root'],combined_root=source['output'],combined_state=str(Path(source['state_directory'])/'status.json'),source_pins=pins,tissue_context_index=source['tissue_context_index'],python=pin(Path(python).resolve()),runtime_python=str(python),runtime_pins=runtime_pins,pysam_version='0.24.1',code_sha256=code,output=str(output),state_directory=str(state),deadline_utc=DEADLINE.isoformat(),max_hours=8,interval_seconds=30,cpu_limit=2,memory_limit_bytes=4*1024**3,max_sqlite_bytes=32*1024**3,max_output_bytes=16*1024**3,max_inflation_ratio=24,min_free_bytes=150*1024**3,max_packet_bytes=64*1024**2,max_header_inflated_bytes=64*1024**2,no_network=True,no_primary_mutation=True)
 Path(rule_path).parent.mkdir(parents=True,exist_ok=True)
 with safe(rule_path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')
 return dict(status=value['review_status'],rule=str(rule_path),rule_sha256=sha(rule_path),controller_started=False)
def validate(rule,path,reviewed_sha,current=None):
 if not reviewed_sha or sha(path)!=reviewed_sha or rule.get('review_status')!='ROOT_REVIEWED_EXACT_FULLY_ANNOTATED_VCF_RULE':raise ValueError('Exact ROOT-reviewed VCF rule required')
 if rule.get('schema_version')!=1 or rule.get('no_network') is not True or rule.get('no_primary_mutation') is not True or rule.get('pysam_version')!='0.24.1':raise ValueError('Local-only exact runtime contract required')
 if rule['cpu_limit']!=2 or rule['memory_limit_bytes']!=4*1024**3 or not 0<rule['max_sqlite_bytes']<=32*1024**3 or not 0<rule['max_output_bytes']<=16*1024**3 or not 1<=rule['max_inflation_ratio']<=24 or rule['max_packet_bytes']!=64*1024**2 or rule['max_header_inflated_bytes']!=64*1024**2 or rule['min_free_bytes']<150*1024**3 or not 0<rule['max_hours']<=8 or not 10<=rule['interval_seconds']<=60:raise ValueError('Reviewed bounded resource limits required')
 deadline=dt.datetime.fromisoformat(rule['deadline_utc']);current=current or now()
 if deadline>DEADLINE or deadline<=current:raise ValueError('Original VM budget deadline cannot extend')
 scripts=Path(__file__).resolve().parent
 if set(rule['code_sha256'])!={str(Path(__file__).resolve()),str(scripts/'embed_vcf_evidence.py')}:raise ValueError('Only exact reviewed controller/helper pair allowed')
 for p,digest in rule['code_sha256'].items():
  if sha(p)!=digest:raise ValueError('Reviewed VCF code changed')
 for item in [*rule['source_pins'],*rule['runtime_pins'],rule['python']]:verify(item)
 workflow=safe(rule['workflow_root']).resolve();combined=safe(rule['combined_root']).resolve();output=safe(rule['output']).resolve();state=safe(rule['state_directory']).resolve()
 if output.name!='fully-annotated-vcf-v1' or output.parent!=workflow.parent or combined!=workflow.parent/'combined-annotation-master-v1' or state!=workflow.parent.parent/'private/fully-annotated-vcf-state-v1' or output.exists() or any(a==b or a in b.parents or b in a.parents for a,b in [(output,state),(output,workflow),(state,workflow),(output,combined)]):raise ValueError('Fresh fixed separate VCF output/private state required')
 return min(deadline,current+dt.timedelta(hours=rule['max_hours']))
def ready(rule):
 state=Path(rule['combined_state'])
 if state.exists():
  label=str(load(state).get('status',''))
  if label.startswith(('FAIL','STOPPED')):raise ValueError('Combined producer failed/stopped; cannot restart upstream')
 root=Path(rule['combined_root']);path=root/'workflow-receipt.json'
 if not path.exists():return False
 value=load(path)
 if str(value.get('status','')).startswith('FAIL'):raise ValueError('Combined upstream failed')
 return value.get('status')=='PASS' and value.get('stage')=='combined-annotation-master-v1' and bool(value.get('identity'))
def run(path,reviewed_sha,execute=False):
 rule=load(path);state=Path(rule['state_directory'])
 if not execute:return load(state/'status.json') if (state/'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
 deadline=validate(rule,path,reviewed_sha);os.umask(0o077);state.mkdir(parents=True,exist_ok=True)
 import fcntl
 with safe(state/'controller.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  def status(label,**extra):
   value=dict(status=label,pid=os.getpid(),rule_sha256=reviewed_sha,updated_utc=now().isoformat(),deadline_utc=deadline.isoformat(),progress_path=str(state/'progress.json'),**extra);temporary=safe(state/'status.tmp');temporary.write_text(json.dumps(value,indent=2)+'\n');os.replace(temporary,safe(state/'status.json'));return value
  try:
   import signal
   def budget_alarm(signum,frame):raise TimeoutError('Original VM hard deadline reached')
   signal.signal(signal.SIGALRM,budget_alarm);signal.setitimer(signal.ITIMER_REAL,max(0.001,(deadline-now()).total_seconds()))
   while not ready(rule):
    if now()>=deadline:return status('STOPPED_ORIGINAL_DEADLINE')
    status('WAITING_COMBINED_MASTER_PASS');time.sleep(min(rule['interval_seconds'],max(0,(deadline-now()).total_seconds())))
   validate(rule,path,reviewed_sha);status('EMBEDDING_PORTABLE_VCF_EVIDENCE')
   helper=Path(__file__).resolve().parent/'embed_vcf_evidence.py';spec=importlib.util.spec_from_file_location('reviewed_vcf_embed',helper);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);result=module.embed(rule,deadline)
   return status('PASS_FULLY_ANNOTATED_VCF',result=result)
  except BaseException as error:status('STOPPED_ORIGINAL_DEADLINE' if now()>=deadline else 'FAILED_FULLY_ANNOTATED_VCF',error=str(error)[:400]);raise
  finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True);a=sub.add_parser('prepare')
 for name in ['combined-rule','python','output','state-directory','rule']:a.add_argument('--'+name,type=Path,required=True)
 a=sub.add_parser('run');a.add_argument('--rule',type=Path,required=True);a.add_argument('--reviewed-rule-sha256');a.add_argument('--execute',action='store_true');args=p.parse_args();result=prepare(args.combined_rule,args.python,args.output,args.state_directory,args.rule) if args.command=='prepare' else run(args.rule,args.reviewed_rule_sha256,args.execute);print(json.dumps(result,indent=2))
