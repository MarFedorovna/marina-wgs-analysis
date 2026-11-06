#!/usr/bin/env python3
"""Bounded own-lock continuation for reviewed direct-packet report views.

Default is status only. No primary workflow process is started, stopped or
modified. Only an exact reviewed rule authorizes the independent report child.
"""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import time

END=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')


def need(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def now():return dt.datetime.now(dt.timezone.utc)


def save(path,value):
    temporary=path.with_suffix('.tmp-'+str(os.getpid()));temporary.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');os.replace(temporary,path)


def checked_rule(path,expected):
    need(sha(path)==expected,'Exact ROOT-reviewed rule SHA256 differs');rule=load(path)
    need(rule.get('root_reviewed') is True and rule.get('schema_version')==1,'Reviewed rule metadata required')
    need(sha(__file__)==rule['controller_sha256'],'Controller code pin differs')
    need(rule['interval_seconds']>=10 and rule['interval_seconds']<=60 and 0<rule['max_hours']<=8,'Bounded polling parameters required')
    need(1<=rule['workers']<=8 and rule.get('synthetic') is False,'Only bounded real report rule accepted')
    need(rule['fully_manifest_path']=='/srv/marina/outputs/fully-annotated-vcf-v1/manifest.json','Exact approved final producer root required')
    need(rule['output']=='/srv/marina/outputs/cohort-sample-views-v1','Exact canonical views root required')
    need(bool(rule.get('authority_expected')),'Source/control/code expectations required')
    for key in rule['authority_expected']:need(key.startswith(('core.','final.')),'Unknown authority expectation namespace')
    for key in ('exporter','proof_module','codec_module','view_policy','ranking_module'):
        item=rule[key];p=Path(item['path']);need(p.is_file() and not p.is_symlink() and p.stat().st_size==item['bytes'] and sha(p)==item['sha256'],'Pinned report helper changed')
    return rule


def continue_views(path,expected,execute=False,watch=False):
    rule=checked_rule(path,expected)
    final=Path(rule['fully_manifest_path'])
    if not execute:return dict(status='STATUS_ONLY_NO_EXECUTION',final_manifest_present=final.is_file(),output_present=Path(rule['output']).exists())
    deadline=min(END,dt.datetime.fromisoformat(rule['deadline_utc']),now()+dt.timedelta(hours=rule['max_hours']))
    need(now()<deadline,'Original deadline exhausted')
    state=Path(rule['state_directory']);state.mkdir(parents=True,exist_ok=True,mode=0o700)
    lock=state/'owns.lock';descriptor=os.open(lock,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    os.write(descriptor,str(os.getpid()).encode());os.close(descriptor)
    def update(status,**extra):
        value=dict(schema_version=1,status=status,pid=os.getpid(),rule_sha256=expected,controller_sha256=rule['controller_sha256'],deadline_utc=deadline.isoformat(),updated_utc=now().isoformat(),primary_process_actuation=False,**extra)
        save(state/'status.json',value);return value
    try:
        while now()<deadline:
            if final.is_file():
                manifest=load(final)
                if manifest.get('status')=='PASS_FULLY_ANNOTATED_VCF_DIRECT':
                    need(not Path(rule['output']).exists(),'Canonical report already exists; refusing overwrite/restart')
                    update('VERIFYING_DIRECT_PACKET_PROOFS_AND_EXPORTING_VIEWS')
                    spec=importlib.util.spec_from_file_location('reviewed_fastviews_exporter',rule['exporter']['path']);exporter=importlib.util.module_from_spec(spec);spec.loader.exec_module(exporter)
                    run=dict(rule,fully_manifest=exporter.pin(final),exporter_sha256=rule['exporter']['sha256'],deadline_utc=deadline.isoformat())
                    result=exporter.export(run)
                    output=Path(rule['output']);manifest_pin=exporter.pin(output/'manifest.json')
                    receipt=dict(schema_version=1,status='PASS',stage='cohort-sample-views-v1',synthetic=False,
                        identity=hashlib.sha256((expected+manifest_pin['sha256']).encode()).hexdigest(),rule_sha256=expected,
                        producer_code_sha256=rule['exporter']['sha256'],report_input_kind='FULLY_ANNOTATED_VCF_PACKETS_V1',
                        artifacts=dict(result['artifacts'],**{'manifest.json':dict(bytes=manifest_pin['bytes'],sha256=manifest_pin['sha256'])}),completed_utc=now().isoformat())
                    with (output/'workflow-receipt.json').open('x',encoding='utf-8') as stream:json.dump(receipt,stream,ensure_ascii=False,indent=2);stream.write('\n')
                    return update('PASS_VCF_FIRST_COHORT_SAMPLE_VIEWS_PRIMARY_UNCHANGED',output_manifest_sha256=manifest_pin['sha256'])
                need(manifest.get('status') not in ('FAIL','FAILED'),'Producer explicitly failed; no automatic restart')
            value=update('WAITING_DIRECT_FULLY_ANNOTATED_VCF_PASS_NO_PRIMARY_ACTUATION')
            if not watch:return value
            time.sleep(min(rule['interval_seconds'],max(0,(deadline-now()).total_seconds())))
        return update('DEADLINE_STOPPED_NO_AUTOMATIC_RESTART')
    except BaseException as error:
        update('FAILED_PARTIALS_RETAINED_NO_AUTOMATIC_RESTART',error_type=type(error).__name__)
        raise
    finally:
        if lock.exists() and lock.read_text()==str(os.getpid()):lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--rule',type=Path,required=True);parser.add_argument('--rule-sha256',required=True)
    parser.add_argument('--execute',action='store_true');parser.add_argument('--watch',action='store_true');args=parser.parse_args()
    print(json.dumps(continue_views(args.rule,args.rule_sha256,args.execute,args.watch),sort_keys=True))
