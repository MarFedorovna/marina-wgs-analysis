#!/usr/bin/env python3
"""Bounded ROOT-reviewed wait: primary PASS plus completed combined enriched master."""
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
def now():return dt.datetime.now(dt.timezone.utc)
def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(rule,path,reviewed_sha,current_time=None):
    if sha(path)!=reviewed_sha:raise ValueError('Exact ROOT-reviewed joined-views rule SHA required')
    fields={'schema_version','review_status','reviewed_utc','code_sha256','ranking_module','ranking_sha256','base_exporter','base_exporter_sha256',
            'joined_manifest','workflow_root','output','state_directory','python_executable','pysam_version','deadline_utc','max_hours','interval_seconds',
            'no_network','no_primary_workflow_mutation'}
    if set(rule)!=fields or rule['schema_version']!=2:raise ValueError('Unknown/missing joined-views rule fields')
    if rule['review_status']!='ROOT_REVIEWED_EXACT_JOINED_COHORT_SAMPLE_VIEWS_V2_RULE':raise ValueError('ROOT review required; candidate cannot execute')
    if rule['no_network'] is not True or rule['no_primary_workflow_mutation'] is not True:raise ValueError('Local-only report contract required')
    current_time=current_time or now();deadline=dt.datetime.fromisoformat(rule['deadline_utc'])
    if deadline>ABSOLUTE_DEADLINE or deadline<=current_time:raise ValueError('Original budget deadline reached/exceeded')
    if not 10<=rule['interval_seconds']<=60 or not 0<rule['max_hours']<=8:raise ValueError('Bounded wait required')
    scripts=Path(__file__).resolve().parent;pair={str(Path(__file__).resolve()),str(scripts/'export_cohort_sample_views_v2.py')}
    if set(rule['code_sha256'])!=pair:raise ValueError('Exact new controller/exporter pair required')
    for code,digest in rule['code_sha256'].items():
        if Path(code).is_symlink() or sha(code)!=digest:raise ValueError('Joined-views code digest differs')
    for field in ('ranking_module','base_exporter'):
        expected=rule['ranking_sha256'] if field=='ranking_module' else rule['base_exporter_sha256']
        if Path(rule[field]).is_symlink() or sha(rule[field])!=expected:raise ValueError('Frozen source digest differs: '+field)
    if str(Path(sys.executable).absolute())!=rule['python_executable'] or rule['pysam_version']!='0.24.1':raise ValueError('Reviewed isolated Python/pysam runtime differs')
    workflow=Path(rule['workflow_root']).resolve();output=Path(rule['output']).resolve();state=Path(rule['state_directory']).resolve()
    if Path(rule['joined_manifest']).resolve()!=workflow.parent/'combined-annotation-master-v1/manifest.json':raise ValueError('Exact combined-master root required')
    if output.exists():raise ValueError('New canonical output required; refusing overwrite')
    if output==workflow or workflow in output.parents or state==workflow or workflow in state.parents or state==output or output in state.parents or state in output.parents:
        raise ValueError('Separate output/state outside primary workflow required')
    return min(deadline,current_time+dt.timedelta(hours=rule['max_hours']))


def readiness(workflow,joined_manifest):
    paths=[Path(workflow)/'stages'/s/'workflow-receipt.json' for s in ('06_annotate','07_report')]
    paths.extend([Path(joined_manifest),Path(joined_manifest).parent/'workflow-receipt.json'])
    try:values=[load(p) for p in paths]
    except FileNotFoundError:return False
    if any(str(v.get('status','')).startswith('FAIL') for v in values):raise ValueError('Upstream failed; no export or automatic restart')
    return all(v.get('status')=='PASS' for v in (values[0],values[1],values[3])) and values[2].get('status')=='PASS_COMBINED_ANNOTATION_MASTER'


def run(path,reviewed_sha,execute=False):
    rule=load(path);state=Path(rule['state_directory'])
    if not execute:return load(state/'status.json') if (state/'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
    deadline=validate(rule,path,reviewed_sha);os.umask(0o077);state.mkdir(mode=0o700,parents=True,exist_ok=True)
    import fcntl
    with (state/'continuation.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Joined views continuation already owns its lock')
        def status(value,**extra):
            doc=dict(status=value,pid=os.getpid(),rule_sha256=reviewed_sha,updated_utc=now().isoformat(),deadline_utc=deadline.isoformat(),**extra)
            p=state/'status.json';temp=state/'status.tmp'
            if p.is_symlink() or temp.is_symlink():raise ValueError('Symlink state refused')
            temp.write_text(json.dumps(doc,indent=2)+'\n');temp.replace(p);return doc
        try:
            while not readiness(rule['workflow_root'],rule['joined_manifest']):
                remaining=(deadline-now()).total_seconds()
                if remaining<=0:return status('STOPPED_DEADLINE_NO_EXPORT_PRIMARY_UNCHANGED')
                status('WAITING_PRIMARY_AND_COMBINED_ENRICHED_MASTER_PASS_NO_PRIMARY_ACTUATION');time.sleep(min(rule['interval_seconds'],remaining))
            validate(rule,path,reviewed_sha)
            import pysam
            if pysam.__version__!=rule['pysam_version']:raise ValueError('Actual pysam release differs')
            helper=Path(__file__).resolve().parent/'export_cohort_sample_views_v2.py'
            spec=importlib.util.spec_from_file_location('joined_views_v2',helper);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
            status('VERIFYING_ORIGINAL_ROW_INVARIANTS_AND_EXPORTING_UNIFIED_RESEARCH_VIEWS')
            result=m.export(rule['workflow_root'],rule['output'],rule['ranking_module'],rule['ranking_sha256'],
                            rule['base_exporter'],rule['base_exporter_sha256'],rule['joined_manifest'],deadline=deadline)
            return status('PASS_UNIFIED_COHORT_AND_SAMPLE_VIEWS_PRIMARY_UNCHANGED',manifest_sha256=sha(Path(rule['output'])/'manifest.json'),counts=result['counts'])
        except BaseException as error:status('FAILED_NO_AUTOMATIC_RETRY_PRIMARY_UNCHANGED',error=str(error)[:500]);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rule',type=Path,required=True);p.add_argument('--reviewed-rule-sha256');p.add_argument('--execute',action='store_true');a=p.parse_args()
    print(json.dumps(run(a.rule,a.reviewed_rule_sha256,a.execute),indent=2))
