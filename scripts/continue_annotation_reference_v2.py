#!/usr/bin/env python3
"""ROOT-reviewable conditional continuation. Default prepares a rule; never approves/runs by default."""
import argparse
import copy
import datetime as dt
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def now():return dt.datetime.now(dt.timezone.utc)
def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def save(path,value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    assert not path.is_symlink() and not temporary.is_symlink()
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def pin(path):
    path=Path(path).resolve(strict=True);assert path.is_file() and not path.is_symlink()
    return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))


def verify_pin(item):
    path=Path(item['path']);assert path.is_file() and not path.is_symlink() and path.stat().st_size==item['bytes'] and sha(path)==item['sha256'],'A ROOT-reviewed file changed: '+str(path)


def check_cfg_delta(base,candidate,adapters,sidecar):
    expected=copy.deepcopy(base)
    expected['annotation']['local_adapters']=[*base['annotation'].get('local_adapters',[]),*adapters]
    expected['annotation']['resource_manifest']=str(sidecar)
    assert candidate==expected,'Config changed outside the exact reviewed annotation adapter/resource fields'


def check_binding_delta(draft,final,sidecar,digest):
    expected=copy.deepcopy(draft);expected['external_upstream_code'][str(sidecar)]=digest
    assert final==expected,'Binding changed outside the one reviewed sidecar digest entry'


def expected_sidecar(core,additions):
    result=copy.deepcopy(core)
    for key,value in additions.items():
        assert key not in core['resources'] or core['resources'][key].get('status') not in ('available','bundled_in_verified_cache'),'Cannot replace a ready locked core resource'
        result['resources'][key]=value
    return result


def exact_approval(path,binding_sha,rule_sha):
    path=Path(path)
    identity=dict(binding_sha256=binding_sha,ROOT_reviewed_controller_rule_sha256=rule_sha,
        status='APPROVED_EXACT_BINDING_AFTER_ROOT_REVIEWED_DETERMINISTIC_CHECKS')
    if path.exists():
        actual=load(path);assert all(actual.get(k)==v for k,v in identity.items()),'Unknown/different existing approval'
        return actual
    assert not path.is_symlink()
    value=dict(identity,approved_utc=now().isoformat())
    with path.open('x') as stream:stream.write(json.dumps(value,indent=2)+'\n')
    return value


def annotation_processes(proc=Path('/proc'),exclude=()):
    result=[]
    for folder in proc.glob('[0-9]*'):
        pid=int(folder.name)
        if pid in exclude:continue
        try:
            argv=[s.decode(errors='replace') for s in (folder/'cmdline').read_bytes().split(b'\0') if s]
            names={Path(s).name for s in argv}
            if ('annotation_queue.py' in names and 'run' in argv and '--execute' in argv) or ('workflow.py' in names and any(v in argv for v in ('annotate','report')) and '--execute' in argv):result.append(pid)
        except OSError:pass
    return sorted(result)


def package_pins(binding):
    package=Path(binding['code_package']['path']);assert sha(package)==binding['code_package']['sha256']
    document=load(package);assert len(document['files'])==15,'Expected unchanged 15-file sealed runtime'
    result=[pin(package)]
    for item in document['files']:
        path=package.parent/item['path'];actual=pin(path)
        assert actual['bytes']==item['bytes'] and actual['sha256']==item['sha256'];result.append(actual)
    return result


def prepare(review_dir,materializer,output,max_hours):
    assert 0<max_hours<=8 and not output.exists(),'Choose a new controller rule file'
    review=load(review_dir/'candidate-review.json');draft=load(review_dir/'bindings-public-hcm-v2.json')
    cfg_path=Path(draft['annotation_config']['path']);cfg=load(cfg_path)
    original_binding=load('/srv/marina/private/annotation-queue-bindings.json');base_cfg_path=Path(original_binding['annotation_config']['path'])
    plan_path=review_dir/'resource-overlay-plan.json';plan=load(plan_path);overlay_path=Path(plan['overlay_path']);overlay=load(overlay_path)
    additions=overlay['prepared_resource_manifest_additions'];adapters=overlay['annotation_config_append_local_adapters']
    assert review['public_controls']=='PASS_ACTUAL_SEALED_PARSERS_PUBLIC_SYNTHETIC_CONTROLS'
    assert {i['name'] for i in adapters}=={'gene_evidence','cadd_public_hcm_union','myvariant_dbnsfp48_public_reference','revel_hcm_reference'} and len(adapters)==4 and not overlay.get('pending_reference_adapters')
    assert sha(base_cfg_path)==original_binding['annotation_config']['sha256'] and sha(cfg_path)==draft['annotation_config']['sha256']
    check_cfg_delta(load(base_cfg_path),cfg,adapters,plan['new_manifest'])
    assert not Path(plan['new_manifest']).exists() and not Path(draft['approval_file']).exists(),'Unexpected existing sidecar/approval requires ROOT review'
    assert draft['code_package']==original_binding['code_package'] and draft['accepted_filter']==original_binding['accepted_filter'] and draft['budget']==original_binding['budget']
    expected_draft=copy.deepcopy(original_binding)
    expected_draft.update(annotation_config=dict(path=str(cfg_path),sha256=sha(cfg_path)),
        approval_file=str(review_dir/'ROOT-approval-required.json'),status_file=str(review_dir/'queue-status.json'),
        log_dir=str(review_dir/'logs'),status='SEALED_FOR_ROOT_REVIEW')
    expected_draft['resources']['path']=plan['new_manifest']
    assert draft==expected_draft,'Draft binding contains an unexpected change outside the reviewed builder outputs'
    small=[pin(Path(__file__)),pin(materializer),pin(materializer.parent/'prefetch_myvariant_public_hcm.py'),pin(base_cfg_path),
        pin(cfg_path),pin(plan_path),pin(overlay_path),pin(review_dir/'public-control-receipt.json'),
        pin(Path(draft['accepted_filter']['path'])/'workflow-receipt.json')]
    assert small[-1]['sha256']==draft['accepted_filter']['receipt_sha256'],'Accepted stage5 receipt changed'
    for filename,expected in draft['external_upstream_code'].items():
        item=pin(filename);assert item['sha256']==expected,'A frozen upstream source changed';small.append(item)
    package=package_pins(draft)
    assets=[]
    for item in adapters:
        value=pin(item['path']);assert value['bytes']==item['bytes'] and value['sha256']==item['sha256'];assets.append(value)
    budget=draft['budget'];assert budget['planned_vm_usd']==40 and budget['hourly_usd']==2.39
    origin=dt.datetime.fromisoformat(budget['vm_created_utc'].replace('Z','+00:00'));deadline=origin+dt.timedelta(hours=40/2.39)
    assert deadline.isoformat()=='2026-11-06T03:46:48.046025+00:00','Original VM budget anchor changed'
    value=dict(schema_version=1,status='PREPARED_FOR_ROOT_CODE_REVIEW_NOT_RUNNING',
        fixed_controller_path=str(Path(__file__).resolve()),materializer=str(materializer),review_dir=str(review_dir),
        pins=small,package_pins=package,adapter_pins=assets,approved_adapters=adapters,resource_additions=additions,
        original_annotation_config=str(base_cfg_path),candidate_annotation_config=str(cfg_path),
        draft_binding_path=str(review_dir/'bindings-public-hcm-v2.json'),draft_binding=copy.deepcopy(draft),
        draft_binding_sha256=sha(review_dir/'bindings-public-hcm-v2.json'),base_core_manifest=plan['base_manifest'],
        sidecar_manifest=plan['new_manifest'],expected_core_recipe_sha256=plan['expected_core_recipe_sha256'],
        expected_vep_image=plan['expected_vep_image'],deadline_utc=deadline.isoformat(),max_hours=max_hours,
        unchanged_15_file_package_sha256=draft['code_package']['sha256'],accepted_stage5_receipt_sha256=draft['accepted_filter']['receipt_sha256'],
        post_report_step=None,no_stage3_5_rerun=True,no_patient_API_queries=True,
        ROOT_execution_instruction='Execute only after reviewing this controller and supplying the exact --reviewed-rule-sha256. This is a deterministic ROOT rule approval, not an approval of an unknown future binding.')
    output.parent.mkdir(mode=0o700,parents=True,exist_ok=True);save(output,value)
    return dict(status=value['status'],rule=str(output),rule_sha256=sha(output),controller_sha256=sha(Path(__file__)),no_approval_created=True,no_queue_started=True)


def verify_inputs(rule,full_assets):
    for item in [*rule['pins'],*rule['package_pins']]:verify_pin(item)
    if full_assets:
        for item in rule['adapter_pins']:verify_pin(item)
    check_cfg_delta(load(rule['original_annotation_config']),load(rule['candidate_annotation_config']),rule['approved_adapters'],rule['sidecar_manifest'])


def materialize_and_approve(rule,rule_sha,materialize_call=None):
    verify_inputs(rule,True)
    draft=rule['draft_binding'];binding_path=Path(rule['draft_binding_path']);sidecar=Path(rule['sidecar_manifest']);base=Path(rule['base_core_manifest'])
    assert sha(binding_path)==rule['draft_binding_sha256'] and not sidecar.exists(),'Unexpected sidecar/binding change before materialization'
    assert not Path(draft['approval_file']).exists(),'Unknown existing approval'
    core=load(base);before=sha(base)
    assert core.get('status')=='PASS' and core.get('config_sha256')==rule['expected_core_recipe_sha256'] and core['vep']['image']==rule['expected_vep_image']
    expected=expected_sidecar(core,rule['resource_additions'])
    if materialize_call:materialize_call()
    else:
        parent=str(Path(rule['materializer']).parent);sys.path.insert(0,parent)
        spec=importlib.util.spec_from_file_location('ROOT_pinned_reference_materializer',rule['materializer']);builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
        builder.materialize(Path(rule['review_dir']))
    assert sha(base)==before,'Original core manifest was overwritten'
    assert load(sidecar)==expected,'Sidecar differs from deep-copy core + exact ROOT-reviewed additions'
    digest=sha(sidecar);final=load(binding_path);check_binding_delta(draft,final,sidecar,digest)
    verify_inputs(rule,True)
    assert not annotation_processes(),'Another annotation queue/workflow child is active; refusing approval/launch'
    approval=exact_approval(final['approval_file'],sha(binding_path),rule_sha)
    return dict(binding_sha256=sha(binding_path),sidecar_sha256=digest,approval=approval)


def execute(rule_path,reviewed_sha):
    assert reviewed_sha and sha(rule_path)==reviewed_sha,'ROOT must provide the exact reviewed rule SHA256'
    rule=load(rule_path);assert rule['schema_version']==1 and sha(rule['fixed_controller_path'])==sha(Path(__file__))
    root=rule_path.parent;state_path=root/'controller-status.json';deadline=dt.datetime.fromisoformat(rule['deadline_utc'])
    end=min(deadline,now()+dt.timedelta(hours=rule['max_hours']))
    assert now()<end,'Original VM allowance/controller bound exhausted'
    with (root/'.controller.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return dict(status='ANOTHER_CONDITIONAL_CONTROLLER_ACTIVE')
        prior=load(state_path) if state_path.exists() else {}
        assert not prior or prior.get('rule_sha256')==reviewed_sha,'Unknown controller state'
        def status(name,**extra):
            value=dict(status=name,rule_sha256=reviewed_sha,pid=os.getpid(),updated_utc=now().isoformat(),deadline_utc=rule['deadline_utc'],**extra);save(state_path,value);print(json.dumps(value),flush=True);return value
        try:
            verify_inputs(rule,True)
            assert not annotation_processes(),'Other annotation queue/workflow child is active; ROOT must stop only its waiting controller first'
            assert sha(rule['draft_binding_path'])==rule['draft_binding_sha256'],'Draft changed; prepare/review a new rule'
            assert not Path(rule['sidecar_manifest']).exists() and not Path(rule['draft_binding']['approval_file']).exists(),'Existing sidecar/approval not owned by this fresh conditional run'
            while now()<end:
                verify_inputs(rule,False)
                assert not annotation_processes(),'Another annotation job appeared while waiting'
                base=Path(rule['base_core_manifest'])
                if base.is_file():
                    core=load(base)
                    assert core.get('status')=='PASS','Unexpected non-PASS original core manifest'
                    break
                status('WAITING_ORIGINAL_VEP_CORE_PASS_NO_APPROVAL_NO_QUEUE')
                time.sleep(min(30,max(0,(end-now()).total_seconds())))
            else:return status('STOPPED_ORIGINAL_BUDGET_OR_CONTROLLER_BOUND',queue_started=False)
            status('VERIFYING_PINNED_MATERIALIZATION_BEFORE_EXACT_APPROVAL')
            result=materialize_and_approve(rule,reviewed_sha)
            assert now()<end and not annotation_processes(),'Budget/other annotation job changed before launch'
            binding=load(rule['draft_binding_path']);queue=Path(binding['code_package']['path']).parent/'scripts/annotation_queue.py'
            argv=[sys.executable,'-u',str(queue),'run','--bindings',rule['draft_binding_path'],'--execute']
            log=(root/'sealed-queue-launch.log').open('ab')
            child=subprocess.Popen(argv,stdin=subprocess.DEVNULL,stdout=log,stderr=log,start_new_session=True)
            status('UNCHANGED_SEALED_QUEUE_RUNNING_PUBLIC20_CONTROLS_THEN_6_7',queue_pid=child.pid,**{k:v for k,v in result.items() if k!='approval'})
            while now()<end and child.poll() is None:time.sleep(min(5,max(0,(end-now()).total_seconds())))
            if child.poll() is None:return status('CONTROLLER_OBSERVATION_BOUND_QUEUE_RETAINS_ORIGINAL_VM_DEADLINE_GUARD',queue_pid=child.pid)
            qstatus=load(binding['status_file']) if Path(binding['status_file']).exists() else {}
            assert child.returncode==0 and qstatus.get('status')=='RESEARCH_REPORT_COMPLETE','Sealed queue failed/stopped; inspect its bounded logs, no automatic restart'
            return status('RESEARCH_REPORT_COMPLETE_POSTSTEP_REMAINS_STANDALONE',queue_pid=child.pid,queue_status=binding['status_file'])
        except BaseException as error:
            status('REFUSED_OR_FAILED_NO_AUTOMATIC_RESTART',error=str(error)[:300]);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--review-dir',type=Path,default=Path('/srv/marina/private/annotation-public-hcm-v2-review'));p.add_argument('--materializer',type=Path,default=Path('/srv/marina/private/reviewed-annotation-reference-candidate-v1/scripts/prepare_annotation_reference_candidate.py'));p.add_argument('--rule',type=Path,default=Path('/srv/marina/private/annotation-reference-controller-v1/rule-candidate.json'));p.add_argument('--max-hours',type=float,default=8);g=p.add_mutually_exclusive_group();g.add_argument('--execute',action='store_true');g.add_argument('--status',action='store_true');p.add_argument('--reviewed-rule-sha256');a=p.parse_args()
    if a.status:
        state=a.rule.parent/'controller-status.json';result=load(state) if state.exists() else dict(status='NOT_STARTED',rule_exists=a.rule.exists())
    else:result=execute(a.rule,a.reviewed_rule_sha256) if a.execute else prepare(a.review_dir,a.materializer,a.rule,a.max_hours)
    print(json.dumps(result,indent=2))
