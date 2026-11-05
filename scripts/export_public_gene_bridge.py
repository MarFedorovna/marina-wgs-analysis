#!/usr/bin/env python3
"""Compact TSV bridge for the existing gene_evidence_tsv adapter; no patient inputs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
from prefetch_myvariant_public_hcm import public_genes,load

CODE=Path(__file__).resolve().parents[1]
FIELDS=['gene','hgnc_id','mondo_id','source_kind','source_release','assertion_status','orpha_ids',
    'clingen_hcm_validity_source_status',
    'exact_hcm_orpha_ids','association_types','assessment_statuses','hpo_term_count','cspec_ids_versions',
    'cspec_approved_dates','source_manifest_sha256','source_document_paths','source_document_sha256',
    'interpretation_status','source_context_index','limitation']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def rows_for_gene(connection,gene,candidate,source_manifest_sha):
    hgnc,symbol=gene['hgnc_id'],gene['symbol'];connection.row_factory=sqlite3.Row
    associations=[dict(r) for r in connection.execute('SELECT * FROM gene_disease WHERE hgnc_id=?',(hgnc,))]
    orpha_ids=sorted({r['orpha_id'] for r in associations})
    exact=[r[0] for r in connection.execute('''SELECT DISTINCT g.orpha_id FROM gene_disease g JOIN disease_xref x USING(orpha_id)
      WHERE g.hgnc_id=? AND x.target_id='MONDO:0005045' AND x.relation LIKE 'E (%' AND x.validation='Validated' ''',(hgnc,))]
    terms={r[0] for orpha in orpha_ids for r in connection.execute('SELECT hpo_id FROM disease_hpo WHERE orpha_id=?',(orpha,))}
    base=dict(gene=symbol,hgnc_id=hgnc,mondo_id='MONDO:0005045',source_manifest_sha256=source_manifest_sha,
        interpretation_status='PUBLIC_GENE_CONTEXT_ONLY_NO_AUTOMATIC_ACMG_PVS1_LOF_GOF_OR_PHENOTYPE_FIT',
        source_context_index='public-gene-context-index.json',limitation='Existing runtime joins canonical SYMBOL; HGNC is source-validated provenance; per-variant HGNC matching must be checked at integration review')
    orpha=dict(base,source_kind='OrphadataScience',source_release='July2026',
        assertion_status='GENE_ASSOCIATIONS_AVAILABLE_NOT_PATIENT_DISEASE_MATCH' if associations else 'MISSING_GENE_ASSOCIATION_NOT_BENIGN',
        orpha_ids=';'.join(orpha_ids),exact_hcm_orpha_ids=';'.join(sorted(exact)),
        association_types=';'.join(sorted({r['association_type'] for r in associations})),
        assessment_statuses=';'.join(sorted({r['assessment'] for r in associations})),hpo_term_count=str(len(terms)),
        source_document_paths='orphadata-sources.json')
    specs=[json.loads(r[0]) for r in connection.execute('SELECT payload FROM cspec WHERE hgnc_id=? AND mondo_id=?',(hgnc,'MONDO:0005045'))]
    documents=[candidate/(s['spec_id']+'-'+s['version']+'.json') for s in specs]
    cspec=dict(base,source_kind='ClinGenCSpec',source_release='registry20261105;247uniqueentities',
        assertion_status='CURRENT_RELEASED_HCM_SPEC_AVAILABLE_REQUIRES_EXPERT_APPLICATION' if specs else 'MISSING_RELEASED_HCM_SPEC_IN_COMPLETE_REGISTRY',
        cspec_ids_versions=';'.join(s['spec_id']+':'+s['version'] for s in specs),
        cspec_approved_dates=';'.join(s.get('approved_on') or '' for s in specs),
        source_document_paths=';'.join(str(p) for p in documents) or 'clingen-cspec-sources.json',
        source_document_sha256=';'.join(sha(p) for p in documents))
    return [orpha,cspec],dict(gene=gene,orphadata_assertions=associations,exact_hcm_orpha_ids=sorted(exact),
        source_hpo_terms=sorted(terms),cspec=[{k:v for k,v in s.items() if k not in ('rules','rule_sets')} for s in specs],
        full_rules_documents=[dict(path=str(p),sha256=sha(p)) for p in documents])


def export(cfg,source,candidate,target):
    assert not target.exists(),'Choose a new integration candidate version; do not overwrite'
    manifest=load(source/'source-manifest.json');cm=load(candidate/'candidate-manifest.json')
    assert manifest.get('status')=='PASS_PUBLIC_SOURCES_ONLY' and cm.get('status')=='CANDIDATE_PUBLIC_ADAPTER_NOT_BOUND'
    assert cm['source_manifest_sha256']==sha(source/'source-manifest.json')
    database=candidate/'public-gene-evidence.sqlite'
    artifact=next(i for i in cm['files'] if i['path']==database.name)
    assert database.stat().st_size==artifact['bytes'] and sha(database)==artifact['sha256']
    genes=public_genes(cfg)
    validity_cfg=dict(cfg);validity_cfg.pop('additional_released_hcm_cspec',None)
    validity_hgnc={g['hgnc_id'] for g in public_genes(validity_cfg)}
    connection=sqlite3.connect(database.resolve().as_uri()+'?mode=ro',uri=True)
    target.mkdir(mode=0o700,parents=True);tsv=target/'public-hcm-gene-evidence.tsv';index=[];rows=[]
    with tsv.open('x',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=FIELDS,delimiter='\t');writer.writeheader()
        for gene in genes:
            emitted,details=rows_for_gene(connection,gene,candidate,sha(source/'source-manifest.json'))
            for row in emitted:
                row['clingen_hcm_validity_source_status']=('EXACT_HCM_MONDO_CURATIONS_PRESENT_KEEP_ORIGINAL_CLASSIFICATIONS' if gene['hgnc_id'] in validity_hgnc else 'MISSING_EXACT_HCM_MONDO_IN_VALIDITY_SNAPSHOT;HCM_CSPEC_IS_SEPARATE;DO_NOT_PROMOTE_DCM_ARVC')
            writer.writerows(emitted);rows.extend(emitted);index.append(details)
    connection.close()
    (target/'public-gene-context-index.json').write_text(json.dumps(dict(status='PUBLIC_GENE_CONTEXT_NOT_PATIENT_MATCH',genes=index),indent=2)+'\n')
    resources={}
    for resource,ids,filename,version in [
        ('orphanet',['orphadata_genes','orphadata_disease_ids','orphadata_phenotypes'],'orphadata-sources.json','OrphadataScienceJuly2026'),
        ('clingen_vcep_cspec',['clingen_cspec_public','clingen_cspec_page2','clingen_cspec_page3','clingen_cspec_registry_metadata'],'clingen-cspec-sources.json','CSpecRegistry20261105;individual_released_versions')]:
        chosen=[r for r in manifest['resources'] if r['source']['id'] in ids];assert len(chosen)==len(ids)
        p=target/filename;p.write_text(json.dumps(dict(status='PUBLIC_SOURCE_DATA_AVAILABLE_NOT_AUTOMATIC_INTERPRETATION',
            sources=chosen,context_index=str(target/'public-gene-context-index.json'),source_files_not_patient_data=True),indent=2)+'\n')
        resources[resource]=dict(status='available_local_supplement',version=version,path=str(p),bytes=p.stat().st_size,sha256=sha(p),
            local_adapter=True,parsed_into_gene_evidence_tsv=str(tsv),scope='PUBLIC_GENE_CONTEXT_ONLY; source-specific assertions retained')
    adapter=dict(name='gene_evidence',status='available',local_adapter=True,format='gene_evidence_tsv',
        version='OrphadataJuly2026+CSpec20261105_hcm_bridge_v1',path=str(tsv),bytes=tsv.stat().st_size,sha256=sha(tsv),
        license_status='OrphadataCC-BY-4.0 attribution retained;ClinGenCC0 attribution/date retained',
        source_url='https://sciences.orphadata.com/genes/',source_urls=['https://sciences.orphadata.com/genes/','https://sciences.orphadata.com/alignments/','https://sciences.orphadata.com/phenotypes/','https://cspec.genome.network/cspec/ui/svi/help'])
    overlay=dict(status='ROOT_REVIEW_CANDIDATE_NOT_APPLIED',annotation_binding_updated=False,
        annotation_config_append_local_adapters=[adapter],prepared_resource_manifest_additions=resources,
        precautions=['Preserve all existing local_adapters; append only after root review.',
            'Do not replace locked VEP/ClinVar/ClinGen/HGNC/MANE resources.',
            'New manifest copy and queue binding/approval required; old sealed runtime unchanged.',
            'TSV gene context reaches gene_context_json, not pathogenicity or gene_validity_rank.',
            'Current loader joins SYMBOL; exact HGNC is preserved and must be compared against VEP HGNC in controls.'])
    (target/'config-overlay.json').write_text(json.dumps(overlay,indent=2)+'\n')
    result=dict(status='ROOT_REVIEW_INTEGRATION_CANDIDATE_NOT_APPLIED',genes=len(genes),tsv_rows=len(rows),
        released_hcm_spec_rows=sum(bool(r.get('cspec_ids_versions')) for r in rows),
        files=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(target.iterdir()) if p.is_file()],
        no_patient_variants_processed=True,automatic_acmg=False)
    (target/'integration-candidate-manifest.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=CODE/'config/myvariant-public-hcm-v2.json');p.add_argument('--source',type=Path,default=Path('/srv/marina/resources/public-supplements-v1'));p.add_argument('--candidate',type=Path,default=Path('/srv/marina/private/publicsupplements-candidate-v1'));p.add_argument('--output',type=Path,default=Path('/srv/marina/private/publicsupplements-integration-candidate-v1'));a=p.parse_args();print(json.dumps(export(load(a.config),a.source,a.candidate,a.output),indent=2))
