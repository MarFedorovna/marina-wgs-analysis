#!/usr/bin/env python3
"""Bounded ROOT-reviewed VCF-first views continuation; no pipeline/process actuation."""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def now():return dt.datetime.now(dt.timezone.utc)


def validate(rule,path,reviewed_sha,current_time=None):
    if sha(path)!=reviewed_sha:raise ValueError('Exact ROOT-reviewed VCF-first rule SHA required')
    keys={'schema_version','review_status','reviewed_utc','code_sha256','ranking_module','ranking_sha256','v1_engine','v1_engine_sha256',
        'v2_helper','v2_helper_sha256','codec_module','codec_sha256','fully_vcf_manifest','fully_manifest_status','control_joined_manifest',
        'packet_append_fields','spliceai_manifest','spliceai_manifest_status','workflow_root','output','state_directory','python_executable','pysam_version','deadline_utc','max_hours','interval_seconds','no_network','no_primary_workflow_mutation'}
    if set(rule)!=keys or rule['schema_version']!=3:raise ValueError('Unknown/missing VCF-first rule fields')
    if rule['review_status']!='ROOT_REVIEWED_EXACT_VCF_FIRST_COHORT_SAMPLE_VIEWS_V3_RULE':raise ValueError('ROOT review required before VCF-first export')
    if rule['no_network'] is not True or rule['no_primary_workflow_mutation'] is not True:raise ValueError('Local-only report contract required')
    current_time=current_time or now();deadline=dt.datetime.fromisoformat(rule['deadline_utc'])
    if deadline>DEADLINE or deadline<=current_time:raise ValueError('Original absolute budget deadline reached/exceeded')
    if not 10<=rule['interval_seconds']<=60 or not 0<rule['max_hours']<=8:raise ValueError('Bounded wait required')
    scripts=Path(__file__).resolve().parent;pair={str(Path(__file__).resolve()),str(scripts/'export_cohort_sample_views_v3.py')}
    if set(rule['code_sha256'])!=pair:raise ValueError('Exact VCF-first controller/exporter pair required')
    for code,digest in rule['code_sha256'].items():
        if Path(code).is_symlink() or sha(code)!=digest:raise ValueError('VCF-first code digest differs')
    for source,digest in [('ranking_module','ranking_sha256'),('v1_engine','v1_engine_sha256'),('v2_helper','v2_helper_sha256'),('codec_module','codec_sha256')]:
        if Path(rule[source]).is_symlink() or sha(rule[source])!=rule[digest]:raise ValueError('Frozen source/codec digest differs: '+source)
    if str(Path(sys.executable).absolute())!=rule['python_executable'] or rule['pysam_version']!='0.24.1':raise ValueError('Reviewed Python/pysam runtime differs')
    workflow=Path(rule['workflow_root']).resolve();output=Path(rule['output']).resolve();state=Path(rule['state_directory']).resolve()
    if Path(rule['fully_vcf_manifest']).resolve()!=workflow.parent/'fully-annotated-vcf-v1/manifest.json':raise ValueError('Exact fully annotated VCF root required')
    if Path(rule['control_joined_manifest']).resolve()!=workflow.parent/'combined-annotation-master-v1/manifest.json':raise ValueError('Exact joined control root required')
    extra=['spliceai_status','spliceai_gene','spliceai_DS_AG','spliceai_DS_AL','spliceai_DS_DG','spliceai_DS_DL','spliceai_DP_AG','spliceai_DP_AL','spliceai_DP_DG','spliceai_DP_DL','spliceai_prediction_json','spliceai_source_json','spliceai_transcript_scope']
    if rule['packet_append_fields']!=extra or Path(rule['spliceai_manifest']).resolve()!=workflow.parent/'spliceai-reference-context-v1/manifest.json' or rule['spliceai_manifest_status']!='PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT':raise ValueError('Exact reviewed SpliceAI supplement contract required')
    if output!=workflow.parent/'cohort-sample-views-v1' or output.exists():raise ValueError('New canonical VCF-first views output required')
    if state==output or state==workflow or output in state.parents or state in output.parents or workflow in state.parents:raise ValueError('Separate private state required')
    return min(deadline,current_time+dt.timedelta(hours=rule['max_hours']))


def readiness(rule):
    paths=[Path(rule['workflow_root'])/'stages'/s/'workflow-receipt.json' for s in ('06_annotate','07_report')]
    paths.extend([Path(rule['fully_vcf_manifest']),Path(rule['fully_vcf_manifest']).parent/'workflow-receipt.json'])
    try:values=[load(p) for p in paths]
    except FileNotFoundError:return False
    if any(str(v.get('status','')).startswith('FAIL') for v in values):raise ValueError('Upstream failed; no automatic restart')
    return all(v.get('status')=='PASS' for v in (values[0],values[1],values[3])) and values[2].get('status')==rule['fully_manifest_status'] and values[3].get('stage')=='fully-annotated-vcf-v1'


def run(path,reviewed_sha,execute=False):
    rule=load(path);state=Path(rule['state_directory'])
    if not execute:return load(state/'status.json') if (state/'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
    deadline=validate(rule,path,reviewed_sha);os.umask(0o077);state.mkdir(mode=0o700,parents=True,exist_ok=True)
    import fcntl
    with (state/'continuation.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('VCF-first continuation already owns its lock')
        def status(value,**extra):
            doc=dict(status=value,pid=os.getpid(),rule_sha256=reviewed_sha,updated_utc=now().isoformat(),deadline_utc=deadline.isoformat(),**extra)
            p=state/'status.json';temp=state/'status.tmp'
            if p.is_symlink() or temp.is_symlink():raise ValueError('Symlink state refused')
            temp.write_text(json.dumps(doc,indent=2)+'\n');temp.replace(p);return doc
        try:
            while not readiness(rule):
                remaining=(deadline-now()).total_seconds()
                if remaining<=0:return status('STOPPED_DEADLINE_NO_EXPORT_PRIMARY_UNCHANGED')
                status('WAITING_PRIMARY_AND_FULLY_ANNOTATED_VCF_PASS_NO_PRIMARY_ACTUATION');time.sleep(min(rule['interval_seconds'],remaining))
            validate(rule,path,reviewed_sha)
            import pysam
            if pysam.__version__!=rule['pysam_version']:raise ValueError('Actual pysam release differs')
            helper=Path(__file__).resolve().parent/'export_cohort_sample_views_v3.py'
            spec=importlib.util.spec_from_file_location('vcf_first_views',helper);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
            status('VERIFYING_VCF_PACKETS_AND_EXPORTING_VIEWS')
            result=module.export(rule,deadline)
            return status('PASS_VCF_FIRST_COHORT_SAMPLE_VIEWS_PRIMARY_UNCHANGED',manifest_sha256=sha(Path(rule['output'])/'manifest.json'),counts=result['counts'])
        except BaseException as error:status('FAILED_NO_AUTOMATIC_RETRY_PRIMARY_UNCHANGED',error=str(error)[:500]);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rule',type=Path,required=True);p.add_argument('--reviewed-rule-sha256');p.add_argument('--execute',action='store_true');a=p.parse_args()
    print(json.dumps(run(a.rule,a.reviewed_rule_sha256,a.execute),indent=2))
