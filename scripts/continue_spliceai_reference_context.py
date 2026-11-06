#!/usr/bin/env python3
"""Prepare a bounded ROOT-reviewed local SpliceAI supplement; no default execution."""
import argparse
import datetime as dt
import importlib.util
import json
import os
from pathlib import Path
import signal
import sys

spec=importlib.util.spec_from_file_location('spliceai_local_helper',Path(__file__).resolve().parent/'run_spliceai_reference_context.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
def prepare(binding_path,audit_path,control_path,python,rule_path,output,state):
 if Path(rule_path).exists() or Path(output).exists():raise ValueError('Fresh candidate/output required')
 binding=h.load(binding_path);audit=h.load(audit_path);control=h.load(control_path);root=Path(binding['data_root']);primary=Path(binding['accepted_filter']['path'])/'workflow-receipt.json';primary_pin=h.pin(primary)
 if primary_pin['sha256']!=binding['accepted_filter']['receipt_sha256'] or h.load(primary)['status']!='PASS':raise ValueError('Exact originally accepted filter receipt required')
 if audit.get('patient_inference') is not False or audit['approved_panel_genes']!=39 or audit['unmatched_symbols']!=['MT-TI'] or control.get('patient_inputs') is not False or control['status']!='PASS_PUBLIC_SPLICEAI_CPU_DISPATCH_EQUIVALENCE' or control['distance']!=500 or control['mask']!=0:raise ValueError('Actual fixed39 audit and matching public control required')
 runtime=root/'resources/spliceai-research-runtime-v1';package=runtime/'venv/lib/python3.10/site-packages/spliceai';panel=h.pin(root/'resources/public-hcm-cds-v2/panel-manifest.json')
 if panel['sha256']!='09c518f2659bd2443ee024d65dd83070922864304ed92f705271e250c227fedb':raise ValueError('Exact approved39 panel required')
 models=audit['models'];annotation=audit['annotation'];annotation={k:annotation[k] for k in ['path','sha256']};annotation['bytes']=Path(annotation['path']).stat().st_size
 for item in [panel,annotation,*models]:h.verify(item)
 wheel=h.pin(runtime/'wheel-manifest.json');purpose=h.pin(runtime/'research-purpose-approval.json')
 if h.load(purpose['path']).get('purpose')!='noncommercial_research':raise ValueError('Existing user noncommercial approval required')
 reference_path=Path(binding['reference']['fasta']);reference=dict(path=str(reference_path),bytes=reference_path.stat().st_size,mtime_ns=reference_path.stat().st_mtime_ns,sha256=binding['reference']['sha256'],validation='accepted original reference provenance plus exact per-allele REF checks; no new whole-FASTA hash')
 source=h.pin(Path(binding['accepted_filter']['path'])/'cohort.filtered.vcf.gz')
 if any(source[k]!=audit['source_vcf'][k] for k in ['path','bytes','sha256']):raise ValueError('Audit and actual accepted source differ')
 source_pins=[h.pin(binding_path),primary_pin,h.pin(audit_path),h.pin(control_path),wheel,purpose,panel,annotation,*models,h.pin(str(reference_path)+'.fai'),h.pin(package/'utils.py'),h.pin(package/'__init__.py')]
 code={str(p):h.sha(p) for p in [Path(__file__).resolve(),Path(__file__).resolve().parent/'run_spliceai_reference_context.py']}
 value=dict(schema_version=1,review_status='PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING',source_pins=source_pins,source_vcf=source,primary_receipt=primary_pin,reference=reference,panel=panel,annotation=annotation,models=models,runtime=dict(python=str(python),wheel_manifest=wheel,version='SpliceAI1.3.1 TF2.15.1 numpy1.26.4'),public_control=h.pin(control_path),audit=h.pin(audit_path),expected_variants=audit['unique_variants_in_bundled_public_gene_spans'],distance=500,mask=0,cpu_limit=2,memory_limit_bytes=8*1024**3,max_sqlite_bytes=1024**3,max_output_bytes=1024**3,min_free_bytes=150*1024**3,max_hours=3.5,deadline_utc=h.DEADLINE.isoformat(),code_sha256=code,output=str(output),state_directory=str(state),no_network=True,no_primary_mutation=True)
 Path(rule_path).parent.mkdir(parents=True,exist_ok=True);h.dump(rule_path,value);return dict(status=value['review_status'],rule=str(rule_path),rule_sha256=h.sha(rule_path),controller_started=False,estimated_seconds_public_control=control['median_direct_seconds']*value['expected_variants'])
def validate(rule,path,digest,current=None):
 if not digest or h.sha(path)!=digest or rule.get('review_status')!='ROOT_REVIEWED_EXACT_SPLICEAI_REFERENCE_RULE':raise ValueError('Exact ROOT-reviewed SpliceAI rule required')
 if rule.get('schema_version')!=1 or rule.get('no_network') is not True or rule.get('no_primary_mutation') is not True or rule['distance']!=500 or rule['mask']!=0 or rule['cpu_limit']!=2 or rule['memory_limit_bytes']!=8*1024**3 or rule['max_sqlite_bytes']>1024**3 or rule['max_output_bytes']>1024**3 or rule['min_free_bytes']<150*1024**3 or not 0<rule['max_hours']<=3.5:raise ValueError('Reviewed scientific/resource bounds required')
 deadline=dt.datetime.fromisoformat(rule['deadline_utc']);current=current or h.now()
 if deadline>h.DEADLINE or deadline<=current:raise ValueError('Original budget cannot extend')
 if Path(rule['output'])!=Path('/srv/marina/outputs/spliceai-reference-context-v1') or Path(rule['state_directory'])!=Path('/srv/marina/private/spliceai-reference-context-state-v1') or Path(rule['output']).exists():raise ValueError('Fresh fixed separate output/state required')
 expected={str(Path(__file__).resolve()),str(Path(__file__).resolve().parent/'run_spliceai_reference_context.py')}
 if set(rule['code_sha256'])!=expected:raise ValueError('Exact helper/controller pair required')
 for p,d in rule['code_sha256'].items():
  if h.sha(p)!=d:raise ValueError('Reviewed SpliceAI code changed')
 for item in rule['source_pins']:h.verify(item)
 if str(Path(sys.executable).absolute())!=rule['runtime']['python']:raise ValueError('Launch exact reviewed existing isolated runtime Python')
 return min(deadline,current+dt.timedelta(hours=rule['max_hours']))
def run(path,digest,execute=False):
 rule=h.load(path);state=Path(rule['state_directory'])
 if not execute:return h.load(state/'status.json') if (state/'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
 deadline=validate(rule,path,digest);os.umask(0o077);state.mkdir(parents=True,exist_ok=True)
 import fcntl
 with h.safe(state/'controller.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  def status(label,**extra):
   value=dict(status=label,pid=os.getpid(),rule_sha256=digest,deadline_utc=deadline.isoformat(),updated_utc=h.now().isoformat(),progress_path=str(state/'progress.json'),**extra);h.atomic(state/'status.json',value);return value
  def alarm(signum,frame):raise TimeoutError('Original budget/job deadline reached')
  signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.001,(deadline-h.now()).total_seconds()))
  try:
   import resource
   resource.setrlimit(resource.RLIMIT_AS,(rule['memory_limit_bytes'],rule['memory_limit_bytes']));status('RUNNING_LOCAL_SPLICEAI_PUBLIC39_REFERENCE_PREDICTIONS');result=h.run(rule,deadline);return status('PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT',result=result)
  except BaseException as error:status('STOPPED_ORIGINAL_DEADLINE' if h.now()>=deadline else 'FAILED_LOCAL_SPLICEAI_REFERENCE_CONTEXT',error=str(error)[:400]);raise
  finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True);a=sub.add_parser('prepare')
 for name in ['binding','audit','public-control','python','rule','output','state-directory']:a.add_argument('--'+name,type=Path,required=True)
 a=sub.add_parser('run');a.add_argument('--rule',type=Path,required=True);a.add_argument('--reviewed-rule-sha256');a.add_argument('--execute',action='store_true');a=p.parse_args();result=prepare(a.binding,a.audit,a.public_control,a.python,a.rule,a.output,a.state_directory) if a.command=='prepare' else run(a.rule,a.reviewed_rule_sha256,a.execute);print(json.dumps(result,indent=2))
