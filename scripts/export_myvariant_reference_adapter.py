#!/usr/bin/env python3
"""Convert complete fixed-public-gene API receipts to exact-allele raw-evidence TSV."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from prefetch_myvariant_public_hcm import Reference,load,validate_hit,save

PREDICTORS={'fathmm':'FATHMM_dbNSFP48_json','fathmm-xf':'FATHMM_XF_dbNSFP48_json',
    'fathmm-mkl':'FATHMM_MKL_dbNSFP48_json','cadd':'CADD_dbNSFP48_json','revel':'REVEL_dbNSFP48_json',
    'metalr':'MetaLR_dbNSFP48_json','varity':'VARITY_dbNSFP48_json','sift':'SIFT_dbNSFP48_json',
    'polyphen2':'PolyPhen_dbNSFP48_json','spliceai':'SpliceAI_dbNSFP48_json','maxentscan':'MaxEntScan_dbNSFP48_json'}
FIELDS=['chrom','pos','ref','alt','transcript','public_query_gene','MYVARIANT_API_build','dbNSFP_release',
    'genomic_mapping_status','dbNSFP_source_genes_json','dbNSFP_Ensembl_json','dbNSFP_AA_json',*PREDICTORS.values(),
    'source_page_path','source_page_sha256','source_license_url','transcript_score_status']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def row_for_hit(hit,gene,reference,page,page_sha):
    validation=validate_hit(hit,gene,reference)
    if not validation['status'].startswith('GRCH38_GENOMIC_REFERENCE_MATCH'):return None,validation['status']
    value=hit['dbnsfp'];pos=value['hg38']['start']
    row=dict(chrom='chr'+str(value['chrom']),pos=str(pos),ref=value['ref'],alt=value['alt'],transcript='*',
        public_query_gene=gene,MYVARIANT_API_build='20250624',dbNSFP_release='4.8a',
        genomic_mapping_status=validation['status'],dbNSFP_source_genes_json=json.dumps(value.get('genename')),
        dbNSFP_Ensembl_json=json.dumps(value.get('ensembl',{}),separators=(',',':')),
        dbNSFP_AA_json=json.dumps(value.get('aa',{}),separators=(',',':')),source_page_path=str(page),source_page_sha256=page_sha,
        source_license_url=value.get('_license','https://sites.google.com/site/jpopgen/dbNSFP'),
        transcript_score_status='RAW_SOURCE_VECTORS_ONLY;UNVERSIONED_OR_NONALIGNED_TRANSCRIPT_ARRAYS;DO_NOT_COPY_TO_EVERY_VEP_TRANSCRIPT')
    for key,column in PREDICTORS.items():row[column]=json.dumps(value[key],separators=(',',':')) if key in value else ''
    # No FATHMM or CADD_PHRED canonical columns: transcript arrays and API CADD version
    # cannot silently override the independently pinned official CADD1.7 exact allele score.
    return row,None


def export(sources,target,reference_path):
    assert not target.exists(),'Choose a new candidate version'
    receipts=[]
    for source in sources:
        panel=load(source/'panel-receipt.json')
        assert panel.get('status')=='COMPLETE_PUBLIC_PANEL_SOURCE_NOT_ENABLED','Public gene prefetch is not complete'
        assert panel['api_build']=='20250624' and panel['dbnsfp_version']=='4.8a'
        receipts.append((source,panel))
    target.mkdir(mode=0o700,parents=True);reference=Reference(reference_path);tsv=target/'myvariant-public-hcm-reference.tsv'
    rows=raw_records=0;excluded={};present={k:0 for k in PREDICTORS};missing={k:0 for k in PREDICTORS};genes=[];gene_pages=[]
    with tsv.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=FIELDS,delimiter='\t');writer.writeheader()
        for source,panel in receipts:
            for gene in panel['completed_genes']:
                assert gene not in genes,'Duplicate public gene source in combined panel'
                genes.append(gene);receipt=load(source/gene/'gene-receipt.json');assert receipt['status']=='COMPLETE_PUBLIC_GENE_SOURCE_NOT_ENABLED'
                counted=0
                for item in receipt['pages']:
                    page=source/gene/item['path'];assert page.stat().st_size==item['bytes'] and sha(page)==item['sha256']
                    doc=load(page);gene_pages.append(dict(path=str(page),sha256=item['sha256'],bytes=item['bytes']))
                    for hit in doc['hits']:
                        raw_records+=1;counted+=1;row,reason=row_for_hit(hit,gene,reference,page,item['sha256'])
                        if row is None:excluded[reason]=excluded.get(reason,0)+1;continue
                        writer.writerow(row);rows+=1
                        for key,column in PREDICTORS.items():
                            if row[column]:present[key]+=1
                            else:missing[key]+=1
                    if raw_records%10000<len(doc['hits']):save(target/'status.json',dict(status='NORMALIZING_PUBLIC_SOURCE_ONLY',raw_records=raw_records,validated_rows=rows))
                assert counted==receipt['records'],'Complete gene receipt/page count differs'
    adapter=dict(name='myvariant_dbnsfp48_public_reference',status='available',local_adapter=True,format='variant_evidence_tsv',
        version='MyVariant20250624_dbNSFP4.8a_publicHCM_union_v2_rawvectors',path=str(tsv),bytes=tsv.stat().st_size,sha256=sha(tsv),
        license_status='Non-commercial research confirmed; retain MyVariant and each underlying source terms; no public redistribution of corpus',source_url='https://docs.myvariant.info/en/latest/doc/data.html')
    result=dict(status='PUBLIC_EXACT_ALLELE_RAW_EVIDENCE_ADAPTER_CANDIDATE',adapter=adapter,public_genes=genes,public_gene_count=len(genes),
        raw_records=raw_records,validated_allele_rows=rows,excluded=excluded,predictor_rows_present=present,predictor_rows_missing=missing,
        source_panel_receipts=[dict(path=str(source/'panel-receipt.json'),sha256=sha(source/'panel-receipt.json')) for source,_ in receipts],
        source_pages_manifest_sha256=hashlib.sha256(json.dumps(gene_pages,sort_keys=True).encode()).hexdigest(),
        no_patient_variants_read=True,source_transcript_versions_not_invented=True,canonical_per_transcript_scores_enabled=False,
        CADD_API_subscore_version='NOT_SPECIFIED;DO_NOT_REPLACE_OFFICIAL_CADD1.7',annotation_binding_updated=False,automatic_acmg=False)
    save(target/'candidate-manifest.json',result);save(target/'status.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,nargs='+',default=[Path('/srv/marina/resources/myvariant-public-hcm-v1'),Path('/srv/marina/resources/myvariant-public-hcm-tnnt2-v1')]);p.add_argument('--output',type=Path,default=Path('/srv/marina/private/myvariant-public-hcm-adapter-candidate-v1'));p.add_argument('--reference',type=Path,default=Path('/srv/marina/references/broad-hg38-v0-altaware-g1575676516681666-prebuilt/Homo_sapiens_assembly38.fasta'));a=p.parse_args();print(json.dumps(export(a.sources,a.output,a.reference),indent=2))
