#!/usr/bin/env python3
"""Single read-only bounded status command: main pipeline and every new public resource."""
import argparse
import json
import math
from pathlib import Path
import time
from monitor_pipeline import snapshot as pipeline_snapshot,load,report

FILES={
 'public_supplement_sources':'resources/public-supplements-v1/source-manifest.json',
 'public_hcm_union_panel':'resources/public-hcm-cds-v2/panel-manifest.json',
 'orpha_cspec_bridge39':'private/publicsupplements-integration-candidate-v2/integration-candidate-manifest.json',
 'cadd_public_cds37':'resources/cadd-public-hcm-v1/status.json',
 'cadd_tnnt2_incremental':'resources/cadd-public-hcm-tnnt2-v1/status.json',
 'myvariant_public_genes38':'resources/myvariant-public-hcm-v1/status.json',
 'myvariant_tnnt2_incremental':'resources/myvariant-public-hcm-tnnt2-v1/status.json',
 'myvariant_raw_adapter':'private/myvariant-public-hcm-adapter-candidate-v1/status.json',
 'revel_archive':'resources/revel-source-only-v1/status.json',
 'revel_public_cds39':'resources/revel-public-hcm-v1/status.json',
 'genereviews':'resources/genereviews-source-only-v1/research-reference-manifest.json',
 'spliceai_wheel':'resources/spliceai-source-only-v1/resource-receipt.json',
 'spliceai_runtime':'resources/spliceai-research-runtime-v1/status.json',
 'spliceai_synthetic_control':'resources/spliceai-research-runtime-v1/synthetic-control-receipt.json',
 'refseq_source_queue':'resources/refseq-vep115-source-only-v1/status.json',
 'lovd_public_myh7_gene':'resources/reference-api-public-catalogs-v1/lovd_myh7_gene-receipt.json',
 'omim_public_identifier_map':'resources/reference-api-public-catalogs-v1/omim_mim2gene-receipt.json',
}
FIELDS=['status','pid','bytes','sha256','received_bytes','expected_bytes','updated_utc','retrieved_utc','validated_utc',
 'deadline_utc','error','models_loaded','reference_alignment','reference_mismatch_rejected','default_annotation',
 'coordinate_review','clinical_pipeline_enabled','rows','output_bytes','actual_http_range_bytes','completed_intervals',
 'total_intervals','gene','gene_records','gene_total','completed_genes','total_genes','raw_bytes','records','raw_records',
 'validated_rows','source_records_scanned','kept_public_records','genes','genes_with_cds','cds_union_bases','tsv_rows','released_hcm_spec_rows',
 'patient_variant_queries','patient_scores_enabled','annotation_binding_updated','read_error']


def snapshot(root,local=False):
    result=pipeline_snapshot(root,local)
    if local:return result
    resources={}
    for name,relative in FILES.items():
        value=load(root/relative)
        resources[name]={k:value[k] for k in FIELDS if k in value} if value else dict(status='NOT_YET_PRESENT')
        resources[name]['receipt']=str(root/relative)
    resources['omim_full_content']=dict(status='ACCESS_REQUIRED',enabled=False,public_identifier_map_is_not_full_content=True)
    resources['varsome_api']=dict(status='ACCESS_REQUIRED',enabled=False)
    result['public_reference_resources']=resources
    result['public_reference_scope']='Public fixed gene/CDS reference preparation; available source does not imply patient analysis or queue binding'
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path('/srv/marina'));g=p.add_mutually_exclusive_group();g.add_argument('--once',action='store_true');g.add_argument('--watch',action='store_true');p.add_argument('--local',action='store_true');p.add_argument('--interval',type=float,default=60);p.add_argument('--max-hours',type=float,default=8);p.add_argument('--report-dir',type=Path);a=p.parse_args()
    if not math.isfinite(a.interval) or not 1<=a.interval<=3600 or not math.isfinite(a.max_hours) or not 0<a.max_hours<=8:raise ValueError('Bounded interval1..3600 and max-hours>0..8 required')
    end=time.monotonic()+a.max_hours*3600
    try:
        while True:
            value=snapshot(a.root,a.local);print(json.dumps(value,ensure_ascii=False),flush=True)
            if a.report_dir:report(value,a.report_dir)
            budget=value.get('databases',{}).get('seconds_until_original_deadline',a.max_hours*3600)
            remaining=min(end-time.monotonic(),budget if isinstance(budget,(int,float)) else a.max_hours*3600)
            if not a.watch or remaining<=0:return
            until=time.monotonic()+min(a.interval,remaining)
            while time.monotonic()<until:time.sleep(min(1,max(0,until-time.monotonic())))
    except KeyboardInterrupt:pass


if __name__=='__main__':main()
