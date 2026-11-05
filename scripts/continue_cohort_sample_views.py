#!/usr/bin/env python3
"""ROOT-reviewed bounded wait for local cohort/sample report views, no reannotation."""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

ABSOLUTE_DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')


def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def now():return dt.datetime.now(dt.timezone.utc)


def validate(rule,path,reviewed_sha,current_time=None):
    if sha(path)!=reviewed_sha:raise ValueError('Exact ROOT-reviewed rule SHA required')
    expected={'schema_version','review_status','reviewed_utc','code_sha256','ranking_module','ranking_sha256','workflow_root','output',
        'state_directory','python_executable','pysam_version','deadline_utc','max_hours','interval_seconds','no_network','no_primary_workflow_mutation'}
    if set(rule)!=expected or rule['schema_version']!=1:raise ValueError('Unknown/missing reviewed rule fields')
    if rule['review_status']!='ROOT_REVIEWED_EXACT_COHORT_SAMPLE_VIEWS_RULE':raise ValueError('ROOT review required; candidate cannot execute')
    if rule['no_network'] is not True or rule['no_primary_workflow_mutation'] is not True:raise ValueError('Local-only separate report contract required')
    if not 10<=rule['interval_seconds']<=60 or not 0<rule['max_hours']<=8:raise ValueError('Bounded wait required')
    current_time=current_time or now();deadline=dt.datetime.fromisoformat(rule['deadline_utc'])
    if deadline>ABSOLUTE_DEADLINE or deadline<=current_time:raise ValueError('Original absolute budget deadline reached/exceeded')
    scripts=Path(__file__).resolve().parent
    pair={str(Path(__file__).resolve()),str(scripts/'export_cohort_sample_views.py')}
    if set(rule['code_sha256'])!=pair:raise ValueError('Exact controller/exporter pair required')
    for code,digest in rule['code_sha256'].items():
        if Path(code).is_symlink() or sha(code)!=digest:raise ValueError('Reviewed code digest differs')
    if Path(rule['ranking_module']).is_symlink() or sha(rule['ranking_module'])!=rule['ranking_sha256']:raise ValueError('Sealed ranking module digest differs')
    if str(Path(sys.executable).absolute())!=rule['python_executable'] or rule['pysam_version']!='0.24.1':raise ValueError('Reviewed isolated Python/pysam runtime differs')
    workflow=Path(rule['workflow_root']).resolve();output=Path(rule['output']).resolve();state=Path(rule['state_directory']).resolve()
    if output.exists():raise ValueError('New output required; refusing overwrite')
    if output==workflow or workflow in output.parents or state==workflow or workflow in state.parents:raise ValueError('Sidecar/state must be outside workflow')
    if output==state or output in state.parents or state in output.parents:raise ValueError('Immutable output and state must be separate')
    return min(deadline,current_time+dt.timedelta(hours=rule['max_hours']))


def readiness(workflow):
    try:receipts=[load(Path(workflow)/'stages'/stage/'workflow-receipt.json') for stage in ('06_annotate','07_report')]
    except FileNotFoundError:return False
    if any(str(r.get('status','')).startswith('FAIL') for r in receipts):raise ValueError('Upstream failed; no automatic restart or export')
    return all(r.get('status')=='PASS' for r in receipts)


def run(path,reviewed_sha,execute=False):
    rule=load(path);state=Path(rule['state_directory'])
    if not execute:return load(state/'status.json') if (state/'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
    deadline=validate(rule,path,reviewed_sha)
    os.umask(0o077);state.mkdir(mode=0o700,parents=True,exist_ok=True)
    import fcntl
    with (state/'continuation.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Views continuation already owns its lock')
        def status(value,**extra):
            doc=dict(status=value,pid=os.getpid(),rule_sha256=reviewed_sha,updated_utc=now().isoformat(),deadline_utc=deadline.isoformat(),**extra)
            p=state/'status.json';temp=state/'status.tmp'
            if p.is_symlink() or temp.is_symlink():raise ValueError('Symlink state refused')
            temp.write_text(json.dumps(doc,indent=2)+'\n');temp.replace(p);return doc
        try:
            while not readiness(rule['workflow_root']):
                remaining=(deadline-now()).total_seconds()
                if remaining<=0:return status('STOPPED_DEADLINE_NO_EXPORT_NO_PRIMARY_ACTUATION')
                status('WAITING_STAGE6_AND_STAGE7_PASS_NO_PRIMARY_ACTUATION');time.sleep(min(remaining,rule['interval_seconds']))
            validate(rule,path,reviewed_sha)
            import pysam
            if pysam.__version__!=rule['pysam_version']:raise ValueError('Actual pysam runtime release differs')
            status('EXPORTING_LOCAL_COHORT_AND_FOUR_SAMPLE_VIEWS_PRIMARY_UNCHANGED')
            helper=Path(__file__).resolve().parent/'export_cohort_sample_views.py'
            spec=importlib.util.spec_from_file_location('reviewed_cohort_views',helper);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
            result=m.export(rule['workflow_root'],rule['output'],rule['ranking_module'],rule['ranking_sha256'],deadline=deadline)
            return status('PASS_SEPARATE_COHORT_SAMPLE_VIEWS_PRIMARY_UNCHANGED',manifest_sha256=m.sha(Path(rule['output'])/'manifest.json'),counts=result['counts'])
        except BaseException as error:
            status('FAILED_NO_AUTOMATIC_RETRY_PRIMARY_UNCHANGED',error=str(error)[:500]);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rule',type=Path,required=True);p.add_argument('--reviewed-rule-sha256');p.add_argument('--execute',action='store_true')
    a=p.parse_args();print(json.dumps(run(a.rule,a.reviewed_rule_sha256,a.execute),indent=2))
