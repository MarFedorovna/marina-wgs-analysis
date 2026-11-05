#!/usr/bin/env python3
"""Prepare recognized annotation config/binding candidates; never approve or run a queue."""
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sqlite3
from prefetch_myvariant_public_hcm import load,save


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value


def public_controls(runtime,adapters,additions,target):
    """Sample exact public MYBPC3 source rows, then exercise actual sealed parsers."""
    a=module('sealed_annotation_public_reference_control',runtime/'scripts/annotation_stages.py')
    control=target/'public-controls';control.mkdir();samples=[]
    for item in adapters:
        original=Path(item['path']);assert original.stat().st_size==item['bytes']
        destination=control/(item['name']+'.tsv');selected=[]
        with original.open() as stream:
            rows=csv.DictReader(stream,delimiter='\t');fields=rows.fieldnames
            for row in rows:
                if item['format']=='gene_evidence_tsv':
                    if row.get('gene') in ('MYBPC3','TNNT2'):selected.append(row)
                elif (row.get('chrom'),row.get('pos'),row.get('ref'),row.get('alt'))==('chr11','47332146','T','G'):
                    selected.append(row)
                    break
        assert selected,'Public control not found in ready adapter: '+item['name']
        with destination.open('x',newline='') as out:
            writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(selected)
        sample=dict(item,path=str(destination),bytes=destination.stat().st_size,sha256=sha(destination))
        samples.append(sample)
    manifest_path=control/'synthetic-core-manifest.json'
    synthetic=dict(schema_version=1,status='PASS',assembly='GRCh38',resources=copy.deepcopy(additions),artifacts=[])
    save(manifest_path,synthetic)
    manifest=a.load_resources(manifest_path,verify_tree=False,local_adapters=samples)
    genes=a.gene_evidence(manifest,manifest_path);assert len(genes['MYBPC3'])==2
    tnnt2=genes['TNNT2'];assert any('GN099:1.0' in e['evidence'].get('cspec_ids_versions','') for e in tnnt2)
    assert all('MISSING_EXACT_HCM_MONDO' in e['evidence'].get('clingen_hcm_validity_source_status','') for e in tnnt2)
    database=a.make_adapter_database(manifest,manifest_path,control/'public-controls.sqlite')
    row=a.annotation_row(dict(chrom='chr11',pos='47332146',ref='T',alt='G',filter='PASS',qual='100'),
        dict(SYMBOL='MYBPC3',HGNC_ID='7551',Feature='ENST00000545968.6',Consequence='missense_variant',Amino_acids='D/A'),
        'PUBLIC_SYNTHETIC_CONTROL',dict(GT='0/1',FT='PASS',DP='30',AD='15,15',GQ='60'),({},manifest,genes,{},database))
    local=json.loads(row['local_evidence_json']);gene_context=json.loads(row['gene_context_json'])
    assert float(row['CADD_PHRED'])==25.0,'Official CADD exact allele did not reach canonical field'
    assert row['FATHMM']=='','Raw, nonaligned API vectors must not become a canonical transcript scalar'
    assert any('FATHMM_dbNSFP48_json' in e and e['FATHMM_dbNSFP48_json'] for e in local)
    assert any(e.get('REVEL_reference_score')=='0.856' for e in local)
    assert {e['evidence']['source_kind'] for e in gene_context}=={'OrphadataScience','ClinGenCSpec'}
    assert row['gene_validity_rank']==0 and 'no_automatic' in row['mechanism_status']
    coverage=a.requirements_coverage({},manifest)
    availability={r['resource_id']:r['status'] for r in coverage['sources'] if r['resource_id'] in ('orphanet','clingen_vcep_cspec','fathmm','revel')}
    assert all(availability.get(key)=='available_local' for key in ('orphanet','clingen_vcep_cspec','fathmm','revel'))
    result=dict(status='PASS_ACTUAL_SEALED_PARSERS_PUBLIC_SYNTHETIC_CONTROLS',runtime_sha256=sha(runtime/'scripts/annotation_stages.py'),
        adapters=len(samples),canonical_CADD_PHRED=row['CADD_PHRED'],canonical_FATHMM=row['FATHMM'],raw_api_predictors_retained=True,
        raw_REVEL_reference_score_retained=True,gene_context_sources=[e['evidence']['source_kind'] for e in gene_context],
        TNNT2_current_HCM_CSpec=True,TNNT2_HCM_validity_missing_preserved=True,availability=availability,
        no_patient_data_read=True,no_ACMG_PVS1_LoF_GoF_inference=True,no_queue_started=True)
    save(target/'public-control-receipt.json',result);return result


def prepare(base_config,base_binding,overlay_path,target):
    assert not target.exists(),'Choose a new private review candidate; do not overwrite'
    original=load(base_config);binding=load(base_binding);overlay=load(overlay_path)
    assert overlay['status']=='ROOT_REVIEW_CANDIDATE_NOT_APPLIED' and not overlay.get('pending_reference_adapters')
    runtime=Path(binding['code_package']['path']).parent
    queue=module('sealed_queue_public_reference_verification',runtime/'scripts/annotation_queue.py')
    queue.verify_bindings(binding)
    assert binding['annotation_config']['sha256']==sha(base_config)
    target.mkdir(mode=0o700,parents=True)
    cfg=copy.deepcopy(original);adapters=overlay['annotation_config_append_local_adapters']
    old_adapters=cfg['annotation'].get('local_adapters',[])
    assert not ({i['name'] for i in old_adapters}&{i['name'] for i in adapters}),'Existing adapter name collision needs review'
    cfg['annotation']['local_adapters']=[*old_adapters,*adapters]
    old_resource=Path(original['annotation']['resource_manifest']);new_resource=old_resource.with_name('resource-manifest-public-hcm-v2.json')
    assert old_resource!=new_resource
    cfg['annotation']['resource_manifest']=str(new_resource)
    assert {k:v for k,v in cfg.items() if k!='annotation'}=={k:v for k,v in original.items() if k!='annotation'}
    allowed={'local_adapters','resource_manifest'}
    assert {k for k in set(cfg['annotation'])|set(original['annotation']) if cfg['annotation'].get(k)!=original['annotation'].get(k)}<=allowed
    config_path=target/'annotation-public-hcm-v2.json';save(config_path,cfg)
    new_binding=copy.deepcopy(binding);new_binding['annotation_config']=dict(path=str(config_path),sha256=sha(config_path))
    new_binding['resources']['path']=str(new_resource)
    new_binding['approval_file']=str(target/'ROOT-approval-required.json')
    new_binding['status_file']=str(target/'queue-status.json');new_binding['log_dir']=str(target/'logs')
    new_binding['status']='SEALED_FOR_ROOT_REVIEW'
    binding_path=target/'bindings-public-hcm-v2.json';save(binding_path,new_binding)
    queue.verify_bindings(new_binding)
    controls=public_controls(runtime,adapters,overlay['prepared_resource_manifest_additions'],target)
    plan=dict(status='WAITING_CORE_VEP_PASS_THEN_ROOT_REVIEW_NOT_APPLIED',base_manifest=str(old_resource),
        new_manifest=str(new_resource),expected_core_recipe_sha256=binding['resources']['expected_recipe_sha256'],
        expected_vep_image=binding['resources']['vep_image'],resource_additions=overlay['prepared_resource_manifest_additions'],
        overlay_path=str(overlay_path),overlay_sha256=sha(overlay_path),new_binding=str(binding_path),
        exact_binding_finalization='After exclusive new sidecar creation, its SHA is inserted into the existing external_upstream_code digest guard; ROOT must review the resulting final binding hash before approval',
        original_core_manifest_not_modified=True,original_15_code_files_unchanged=True)
    save(target/'resource-overlay-plan.json',plan)
    result=dict(status='PRIVATE_ANNOTATION_CONFIG_BINDING_DRAFT_READY_FOR_ROOT_REVIEW',
        annotation_config=dict(path=str(config_path),sha256=sha(config_path)),binding_draft=dict(path=str(binding_path),sha256=sha(binding_path)),
        source_stage5_receipt_sha256=binding['accepted_filter']['receipt_sha256'],
        unchanged_code_package_sha256=binding['code_package']['sha256'],
        changed_config_fields=['annotation.local_adapters','annotation.resource_manifest'],
        recognized_local_adapter_names=[i['name'] for i in adapters],public_controls=controls['status'],
        final_resource_manifest_not_created=True,ROOT_approval_not_created=True,queue_not_started=True,
        note='Run --materialize-resources only after original VEP core PASS; re-review final exact binding hash then approve separately')
    save(target/'candidate-review.json',result);return result


def materialize(target):
    plan=load(target/'resource-overlay-plan.json');base=Path(plan['base_manifest']);new=Path(plan['new_manifest'])
    binding_path=Path(plan['new_binding']);binding=load(binding_path)
    assert not Path(binding['approval_file']).exists(),'Never mutate an approved binding'
    assert new.parent==base.parent and new.name=='resource-manifest-public-hcm-v2.json' and new!=base
    assert sha(Path(plan['overlay_path']))==plan['overlay_sha256']
    core=load(base)
    assert core.get('status')=='PASS' and core.get('config_sha256')==plan['expected_core_recipe_sha256']
    assert core['vep']['image']==plan['expected_vep_image']
    merged=copy.deepcopy(core)
    for key,item in plan['resource_additions'].items():
        assert key not in merged['resources'] or merged['resources'][key].get('status') not in ('available','bundled_in_verified_cache')
        merged['resources'][key]=item
    data=(json.dumps(merged,indent=2,ensure_ascii=False)+'\n').encode()
    assert not new.is_symlink()
    if new.exists():assert new.read_bytes()==data,'Existing sidecar differs'
    else:
        with new.open('xb') as out:out.write(data)
    # This existing map is an actual generic path->digest guard in verify_bindings;
    # including the sidecar binds it without changing any of the 15 sealed files.
    binding['external_upstream_code'][str(new)]=sha(new);save(binding_path,binding)
    runtime=Path(binding['code_package']['path']).parent;queue=module('sealed_queue_final_reference_binding',runtime/'scripts/annotation_queue.py');queue.verify_bindings(binding)
    result=dict(status='FINAL_BINDING_REQUIRES_ROOT_REVIEW_NO_APPROVAL_CREATED',binding=str(binding_path),binding_sha256=sha(binding_path),new_resource_manifest=str(new),new_resource_manifest_sha256=sha(new),original_manifest_sha256=sha(base),unchanged_code_package_sha256=binding['code_package']['sha256'],recognized_existing_digest_guard=True,queue_not_started=True)
    save(target/'final-binding-review.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--base-config',type=Path,default=Path('/srv/marina/private/annotation-v1.json'));p.add_argument('--base-bindings',type=Path,default=Path('/srv/marina/private/annotation-queue-bindings.json'));p.add_argument('--overlay',type=Path,default=Path('/srv/marina/private/public-hcm-union-integration-candidate-v2/config-overlay.json'));p.add_argument('--output',type=Path,default=Path('/srv/marina/private/annotation-public-hcm-v2-review'));p.add_argument('--materialize-resources',action='store_true');a=p.parse_args();print(json.dumps(materialize(a.output) if a.materialize_resources else prepare(a.base_config,a.base_bindings,a.overlay,a.output),indent=2))
