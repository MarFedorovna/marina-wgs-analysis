#!/usr/bin/env python3
"""Expose every audited requirement/resource without claiming unfinished work is done."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import re

CODE=Path(__file__).resolve().parents[1]


def load(path):
    return json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {}


def requirements(path):
    rows=[]
    for line in path.read_text(encoding='utf-8').splitlines():
        if not line.startswith('| '):continue
        columns=[v.strip() for v in line.strip().strip('|').split('|')]
        refs=re.findall(r'\b(?:USER\d{2}|U\d{2}|M\d{2}|G\d{2})\b',columns[0])
        if not refs or len(columns)<3:continue
        policy=columns[-1]
        category=('withdrawn' if 'WITHDRAWN' in policy else 'separate_branch' if 'SEPARATE' in policy or 'SOMATIC' in policy
                  else 'expert_or_missing_input' if 'MANUAL' in policy else 'research_output_required')
        rows.append(dict(id='REQ%03d'%(len(rows)+1),source_refs=refs,source_locator=columns[0],
            requirement=columns[1],required_output=columns[2],policy=policy,category=category,
            actual_status='WITHDRAWN_NOT_RUN' if category=='withdrawn' else 'NOT_PROVEN_COMPLETE_BY_RESOURCE_AVAILABILITY'))
    return rows


def coverage(server_root=None):
    registry=load(CODE/'config/evidence-sources.example.json')['resources']
    policy=load(CODE/'config/resources-inventory.example.json')
    aliases={'clingen_gene_validity':'clingen_validity','ncbi_refseq_gene_mane':'mane',
             'hpo':'hpo_gene','clingen_vcep_cspec':'clingen_vcep_cspec','cadd_uncertain':'cadd'}
    auto={'clinvar','clingen_validity','clingen_dosage','hpo_gene','hgnc','mane','alphamissense'}
    cache={'gnomad','sift','polyphen2','ensembl','dbsnp'}
    profiles=load(CODE/'config/predictor-installation-profiles-v1.json')
    profile_by_id={i['id']:i for key in ['candidates_install_now','not_now'] for i in profiles.get(key,[])}
    observed={}
    supplement={}
    if server_root:
        for p in (server_root/'resources/annotation/.public-source-cache').glob('*/resource.json'):
            item=load(p); actual=p.parent/item.get('path','MISSING')
            if item.get('status')=='available' and actual.is_file() and actual.stat().st_size==item.get('bytes'):
                observed[item.get('name')]=dict(status='RESOURCE_RECEIPT_AVAILABLE_NOT_FRESH_HASH',version=item.get('version'),bytes=item.get('bytes'),sha256=item.get('sha256'))
        alpha=load(server_root/'resources/annotation-supplements/alphamissense-2023-v2-grch38-v3/resource-manifest.json')
        if alpha.get('status')=='PASS':observed['alphamissense']=dict(status='RESOURCE_MANIFEST_PASS_NOT_PATIENT_ANNOTATION',version=alpha.get('dataset_version'))
        supplement=load(server_root/'resources/public-supplements-v1/source-manifest.json')
    sources=[]
    for item in registry:
        name=item['id'];key=aliases.get(name,name)
        if name=='kegg_enrichment_withdrawn':mode='withdrawn'
        elif name=='alphafold_prior_work':mode='historical_superseded_no_legacy_reuse'
        elif 'somatic' in item.get('scope','').lower():mode='separate_branch'
        elif key in auto:mode='automatic_public_resource'
        elif key in cache:mode='cache_115_waiting_verified_metadata_and_controls'
        elif name in ('orphanet','clingen_vcep_cspec'):mode='new_public_supplement_candidate_not_bound'
        elif name in ('omim','mastermind') or key in policy.get('access_required',{}):mode='access_required'
        elif key in profile_by_id:mode='predictor_profile_not_enabled'
        elif key in policy.get('optional',{}):mode='separate_installation_profile_required'
        else:mode='manual_or_offline_export_required'
        status='MISSING_NOT_BENIGN_NOT_ZERO'
        if mode=='withdrawn':status='WITHDRAWN_NOT_RUN'
        if mode=='historical_superseded_no_legacy_reuse':status='NOT_REUSED_NOT_REQUIRED_MODELING'
        row=dict(item,coverage_mode=mode,actual_status=status)
        if key in observed:row.update(actual_status=observed[key]['status'],observed_resource=observed[key])
        if name in ('orphanet','clingen_vcep_cspec') and supplement.get('status')=='PASS_PUBLIC_SOURCES_ONLY':
            wanted=['orphadata_genes','orphadata_disease_ids','orphadata_phenotypes'] if name=='orphanet' else ['clingen_cspec_public','clingen_cspec_page2','clingen_cspec_page3','clingen_cspec_registry_metadata']
            receipts={r['source']['id']:r for r in supplement['resources']}
            if all(k in receipts and receipts[k].get('status')=='PASS' for k in wanted):
                row.update(actual_status='PUBLIC_SOURCE_DATA_PASS_ADAPTER_CANDIDATE_NOT_BOUND',source_receipts=[{k:r[k] for k in ['bytes','sha256','retrieved_utc']} for r in [receipts[w] for w in wanted]])
        if key in profile_by_id:row['installation_profile']=profile_by_id[key]
        sources.append(row)
    extra=[dict(id=key,coverage_mode='predictor_profile_not_enabled',actual_status='MISSING_NOT_BENIGN_NOT_ZERO',installation_profile=value)
           for key,value in profile_by_id.items() if key not in {aliases.get(s['id'],s['id']) for s in sources}]
    rows=requirements(CODE/'docs/requirements-matrix.md')
    source_docs=[]
    for prefix,start,end in [('U',1,19),('M',1,41),('G',42,44),('USER',1,3)]:
        for n in range(start,end+1):
            source_id=prefix+str(n).zfill(2)
            matched=[r['id'] for r in rows if source_id in r['source_refs']]
            source_docs.append(dict(id=source_id,requirements=matched,
                status='EMPTY_FAILED_NOTES_NO_REQUIREMENT_INFERRED' if source_id=='M26' else 'AUDITED_REQUIREMENTS_REGISTERED' if matched else 'AUDIT_REFERENCE_REQUIRES_REVIEW'))
    return dict(schema_version=1,generated_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        scope='All audited Ulugbek/MedInter/GATK requirements; source availability is not completion',
        registry_source_count=len(registry),resources=sources,additional_predictor_profiles=extra,
        audited_sources=source_docs,requirements=rows,
        mandatory_research_outputs=['separate stop_gained','separate stop_lost','start_lost','frameshift','splice',
            'per-transcript consequences/positions/HGVS/version','top10 and top20 per sample and cohort',
            'AlphaMissense exact transcript/protein mapping; N/A vs missing','carrier/intersection matrices','raw and quality evidence'],
        interpretation_gates=['No automatic ACMG/PVS1/LoF/GoF assignment','No KEGG enrichment reinstatement',
            'No legacy candidate reuse','Missing is not benign or AF zero','Public gene evidence is not a patient diagnosis',
            'No patient variant APIs','Sealed annotation runtime approval/binding unchanged; candidate requires root review'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--server-root',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=coverage(a.server_root)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(resources=len(result['resources']),requirements=len(result['requirements']),audited_sources=len(result['audited_sources']),output=str(a.output))))
