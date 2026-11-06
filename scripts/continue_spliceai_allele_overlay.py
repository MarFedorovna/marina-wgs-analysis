#!/usr/bin/env python3
"""ROOT-reviewed waiting controller for the final core33 plus allele-info VCF."""
import argparse
import datetime as dt
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import time

spec=importlib.util.spec_from_file_location('reviewable_allele_overlay',Path(__file__).resolve().parent/'apply_spliceai_allele_overlay.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
FILES=['continue_spliceai_allele_overlay.py','apply_spliceai_allele_overlay.py','spliceai_allele_info_codec.py','spliceai_packet_context.py','embed_vcf_evidence.py']
def prepare(base_path,spliceai_path,python,rule_path,output,state):
 if Path(rule_path).exists() or Path(output).exists():raise ValueError('Fresh private candidate/output required')
 base=h.load(base_path);model=h.load(spliceai_path)
 if base.get('status') not in ['ROOT_REVIEWED_DIRECT_CORE33_RULE','ROOT_REVIEWED_DIRECT_CORE33_V2_RULE'] or base.get('root_reviewed') is not True or model.get('review_status')!='ROOT_REVIEWED_EXACT_SPLICEAI_REFERENCE_RULE':raise ValueError('Actual ROOT-reviewed core33 and SpliceAI source rules required')
 if base['accepted_filter_vcf']!=model['source_vcf']:raise ValueError('Core and model accepted original input differ')
 dependency=[h.pin(base_path),h.pin(spliceai_path),base['driver'],base['core'],base['context_module'],base['codec_module'],*base['context_code']]
 for path,digest in model['code_sha256'].items():
  item=h.pin(path)
  if item['sha256']!=digest:raise ValueError('Reviewed scorer code differs')
  dependency.append(item)
 for item in dependency:h.verify(item)
 metadata=json.loads(subprocess.run([str(python),'-c','import pysam,json;from pathlib import Path;p=Path(pysam.__file__).parent;print(json.dumps(dict(version=pysam.__version__,files=[str(f) for f in sorted(p.iterdir()) if f.suffix in [".py",".so"]])))'],capture_output=True,text=True,check=True).stdout)
 if metadata['version']!='0.24.1':raise ValueError('Existing pinned pysam0.24.1 runtime required')
 code={str(Path(__file__).resolve().parent/name):h.sha(Path(__file__).resolve().parent/name) for name in FILES}
 value=dict(schema_version=1,review_status='PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING',base_root=base['output'],spliceai_root=model['output'],base_state=str(Path(base['state_directory'])/'status.json'),spliceai_state=str(Path(model['state_directory'])/'status.json'),source_vcf=model['source_vcf'],expected_core_driver_sha256=base['driver']['sha256'],expected_model_provenance={k:model[k] for k in ['models','annotation','public_control','code_sha256']},dependency_pins=dependency,runtime_python=str(python),python=h.pin(Path(python).resolve()),runtime_pins=[h.pin(p) for p in metadata['files']],code_sha256=code,output=str(output),state_directory=str(state),cpu_limit=2,memory_limit_bytes=4*1024**3,max_output_bytes=48*1024**3,min_free_bytes=150*1024**3,max_hours=8,interval_seconds=30,deadline_utc=h.DEADLINE.isoformat(),no_network=True,no_primary_or_core_mutation=True,core_packets_reencoded=False)
 Path(rule_path).parent.mkdir(parents=True,exist_ok=True);h.dump(rule_path,value);return dict(status=value['review_status'],rule=str(rule_path),rule_sha256=h.sha(rule_path),controller_started=False)
def validate(rule,path,digest,current=None):
 if not digest or h.sha(path)!=digest or rule.get('review_status')!='ROOT_REVIEWED_EXACT_SPLICEAI_ALLELE_OVERLAY_RULE':raise ValueError('Exact ROOT-reviewed allele overlay rule required')
 if rule.get('schema_version')!=1 or rule.get('no_network') is not True or rule.get('no_primary_or_core_mutation') is not True or rule.get('core_packets_reencoded') is not False or rule['cpu_limit']!=2 or rule['memory_limit_bytes']!=4*1024**3 or not 0<rule['max_output_bytes']<=48*1024**3 or rule['min_free_bytes']<150*1024**3 or not 0<rule['max_hours']<=8 or not 10<=rule['interval_seconds']<=60:raise ValueError('Reviewed scientific/resource limits required')
 deadline=dt.datetime.fromisoformat(rule['deadline_utc']);current=current or h.now()
 if deadline>h.DEADLINE or deadline<=current:raise ValueError('Original VM deadline cannot extend')
 if set(rule['code_sha256'])!={str(Path(__file__).resolve().parent/name) for name in FILES}:raise ValueError('Exact codec/helper/controller source set required')
 for path,expected in rule['code_sha256'].items():
  if h.sha(path)!=expected:raise ValueError('Reviewed overlay source changed')
 for item in [*rule['dependency_pins'],*rule['runtime_pins'],rule['python']]:h.verify(item)
 if Path(rule['base_root'])!=Path('/srv/marina/outputs/direct-core-context-v1') or Path(rule['spliceai_root'])!=Path('/srv/marina/outputs/spliceai-reference-context-v1') or Path(rule['output'])!=Path('/srv/marina/outputs/fully-annotated-vcf-v1') or Path(rule['state_directory'])!=Path('/srv/marina/private/spliceai-allele-overlay-state-v1') or Path(rule['output']).exists():raise ValueError('Fresh fixed separate final root/state required')
 return min(deadline,current+dt.timedelta(hours=rule['max_hours']))
def ready(rule):
 for field in ['base_state','spliceai_state']:
  path=Path(rule[field])
  if path.exists() and str(h.load(path).get('status','')).startswith(('FAIL','STOPPED')):raise ValueError('Required source failed/stopped; final overlay cannot restart it')
 for root,stage in [(rule['base_root'],'direct-core-context-v1'),(rule['spliceai_root'],'spliceai-reference-context-v1')]:
  path=Path(root)/'workflow-receipt.json'
  if not path.exists():return False
  value=h.load(path)
  if value.get('status')!='PASS' or value.get('stage')!=stage:return False
 return True
def run(path,digest,execute=False):
 rule=h.load(path);state=Path(rule['state_directory'])
 if not execute:return h.load(state/'status.json') if (state/'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
 deadline=validate(rule,path,digest);os.umask(0o077);state.mkdir(parents=True,exist_ok=True)
 import fcntl
 with h.safe(state/'controller.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  def status(label,**extra):
   value=dict(status=label,pid=os.getpid(),rule_sha256=digest,updated_utc=h.now().isoformat(),deadline_utc=deadline.isoformat(),progress_path=str(state/'progress.json'),**extra);tmp=state/'status.tmp';tmp.write_text(json.dumps(value,indent=2)+'\n');os.replace(tmp,state/'status.json');return value
  def alarm(signum,frame):raise TimeoutError('Original VM hard deadline reached')
  signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.001,(deadline-h.now()).total_seconds()))
  try:
   while not ready(rule):
    if h.now()>=deadline:return status('STOPPED_ORIGINAL_DEADLINE')
    status('WAITING_CORE33_AND_SPLICEAI_PASS');time.sleep(min(rule['interval_seconds'],max(0,(deadline-h.now()).total_seconds())))
   validate(rule,path,digest);status('APPENDING_ACTUAL_SPLICEAI_ALLELE_INFO_CORE_PACKETS_UNCHANGED');result=h.overlay(rule,deadline);return status('PASS_FULLY_ANNOTATED_VCF_DIRECT',result=result)
  except BaseException as error:status('STOPPED_ORIGINAL_DEADLINE' if h.now()>=deadline else 'FAILED_SPLICEAI_ALLELE_OVERLAY',error=str(error)[:400]);raise
  finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True);a=sub.add_parser('prepare')
 for name in ['base-rule','spliceai-rule','python','rule','output','state-directory']:a.add_argument('--'+name,type=Path,required=True)
 a=sub.add_parser('run');a.add_argument('--rule',type=Path,required=True);a.add_argument('--reviewed-rule-sha256');a.add_argument('--execute',action='store_true');a=p.parse_args();result=prepare(a.base_rule,a.spliceai_rule,a.python,a.rule,a.output,a.state_directory) if a.command=='prepare' else run(a.rule,a.reviewed_rule_sha256,a.execute);print(json.dumps(result,indent=2))
