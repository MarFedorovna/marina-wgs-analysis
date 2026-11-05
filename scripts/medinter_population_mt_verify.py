#!/usr/bin/env python3
"""Offline actual public-sidecar controls; never inspect patient outputs or call APIs."""
import hashlib
import json
from pathlib import Path
import re
import sqlite3
from medinter_population_mt_sources import load,save,now,safe
from medinter_population_mt_common import panel,IntervalLookup,population_fields,decode_population_info

BASE=Path('/srv/marina/resources')
SIDE=BASE/'medinter-population-mt-sidecars-v1'
def connection(path):return sqlite3.connect(f'file:{safe(path)}?mode=ro',uri=True)
def artifact(path):
 digest=hashlib.sha256()
 with safe(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024**2),b''):digest.update(block)
 return dict(path=str(path),bytes=path.stat().st_size,sha256=digest.hexdigest(),data_scope='public reference only; no patient inputs')
def ready_manifest():
 controls=SIDE/'actual-public-source-controls.json';data=load(controls)
 if data.get('status')!='PASS_PUBLIC_SOURCE_SIDECAR_CONTROLS_NOT_CLINICAL_VALIDATION' or data.get('patient_inputs') is not False or data.get('annotation_binding_updated') is not False:raise ValueError('Actual public source controls required')
 output=SIDE/'resource-ready-manifest.json'
 if output.exists():print(json.dumps(load(output),indent=2));return
 value=dict(status='AVAILABLE_VERIFIED_PUBLIC_REFERENCE_SIDECARS_NOT_BOUND',validation_receipt=str(controls),validation_receipt_sha256=hashlib.sha256(controls.read_bytes()).hexdigest(),artifacts=data['artifacts'],checks=data['checks'],source_versions=dict(nuclear_counts='gnomAD4.1.1 exomes/genomes; fixed public GENCODE49 CDS',all_sites_AN='gnomAD4.1 exomes only',coverage='gnomad_r4 API, documented v4 exome/v3 genome coverage; exact deployment release not returned',MT='gnomAD3.1 and HelixMTdb20200327; NC_012920.1 validated',ClinVar='2026-10-31/2026-11-01 tabular references; selected submission SCVs all versioned'),population_INFO_storage='SQLite INFO values are zlib-compressed UTF-8 BLOBs; scalar AC/AN/AF/nhomalt/FILTER columns remain queryable',explicit_limitations=['MITOMAP official endpoint HTTP403; not obtained','Genomes all-sites AN not downloaded; absent records never assumed AC=0','Coverage is public gene-bin context with source version/aggregation limitations','No patient variants selected resources or were sent to APIs','ClinVar source preparation is not clinical review; literature and independent cases not established','No CNV/SV calling or analysis; no automatic pathogenicity classification'],annotation_binding_updated=False,patient_inputs=False,created_utc=now().isoformat())
 save(output,value);print(json.dumps(dict(status=value['status'],path=str(output),validation_receipt_sha256=value['validation_receipt_sha256']),indent=2))
def main():
 manifest,intervals=panel();lookup=IntervalLookup(intervals);results={}
 mt=load(SIDE/'mt-resource-receipt.json');db=connection(Path(mt['database']))
 gnomad=db.execute('SELECT COUNT(*) FROM gnomad').fetchone()[0];helix=db.execute('SELECT COUNT(*) FROM helix').fetchone()[0];helix_with_AN=db.execute('SELECT COUNT(AN) FROM helix').fetchone()[0];quarantine=db.execute("SELECT COUNT(*) FROM quarantine WHERE reason='MULTIALLELIC_COUNTS_NOT_ASSIGNED_TO_SINGLE_ALLELE'").fetchone()[0]
 if (gnomad,helix,helix_with_AN,quarantine)!=(18164,14104,0,220) or mt['counts']['ref_quarantine']!=0 or not mt['official_rcrs_matches_local_GRCh38_chrM']:raise ValueError('Actual MT source controls failed')
 results['mt']=dict(status='PASS_ACTUAL_PUBLIC_MT_SOURCE_CONTROLS',gnomad_rows=gnomad,helix_rows=helix,helix_missing_AN_preserved=True,multiallelic_not_flattened=quarantine,rcrs_REF_mismatches=0);db.close()
 receipt=load(SIDE/'clinvar-dossier-resource-receipt.json');db=connection(Path(receipt['database']));submissions=0;unversioned=0
 for (scv,) in db.execute('SELECT SCV FROM submissions'):submissions+=1;unversioned+=not bool(re.fullmatch(r'SCV\d+\.\d+',scv))
 variants=db.execute('SELECT COUNT(DISTINCT variation_id) FROM variants').fetchone()[0];citations=db.execute('SELECT COUNT(*) FROM citations').fetchone()[0]
 if (submissions,unversioned,variants,citations)!=(248946,0,109528,207404):raise ValueError('Actual ClinVar source controls failed')
 if db.execute('SELECT COUNT(*) FROM submissions WHERE variation_id NOT IN (SELECT variation_id FROM variants)').fetchone()[0]:raise ValueError('SCV not in public gene-selected reference variants')
 results['clinvar']=dict(status='PASS_ACTUAL_PUBLIC_GENE_CLINVAR_SOURCE_CONTROLS',variation_ids=variants,SCV_rows=submissions,SCV_version_missing=unversioned,citation_links=citations,literature_read=False,independent_case_counts_established=False,PP5_assigned=False);db.close()
 an=load(SIDE/'an-resource-receipt.json')
 if not an or not an['status'].startswith('READY_'):raise ValueError('AN source subset not ready')
 db=connection(SIDE/'public-hcm39-exome-AN.sqlite');rows=0
 for chrom,pos,value in db.execute('SELECT chrom,pos,AN FROM all_sites_AN'):
  if value<0 or not lookup.contains(chrom,pos):raise ValueError('AN coordinate/denominator outside approved public CDS')
  rows+=1
 if rows!=an['rows'] or rows==0:raise ValueError('AN row count mismatch')
 results['all_sites_AN']=dict(status='PASS_ACTUAL_PUBLIC_CDS_AN_CONTROLS',rows=rows,genomes_AN_not_downloaded=True,missing_not_zero=True);db.close()
 nuclear=BASE/'medinter-population-hcm39-v3';nuclear_receipt=load(nuclear/'resource-receipt.json')
 if not nuclear_receipt or not nuclear_receipt['status'].startswith('READY_'):raise ValueError('Nuclear count subset not ready')
 db=connection(nuclear/'public-population.sqlite');rows=0;null_AN=0
 for source,chrom,pos,ref,alt,filters,AC,AN,AF,nhomalt,info in db.execute('SELECT * FROM variants'):
  info=decode_population_info(info)
  if not lookup.contains(chrom,pos):raise ValueError('Nuclear count coordinate outside public CDS')
  parsed=population_fields([chrom,str(pos),'.',ref,alt,'.',filters,info])
  if (parsed['AC'],parsed['AN'],parsed['AF'],parsed['nhomalt'])!=(AC,AN,AF,nhomalt):raise ValueError('Nuclear allele arrays altered in storage')
  if any(part.split('=',1)[0].lower() in ['vep','csq'] for part in info.split(';')):raise ValueError('Redundant large consequence payload not omitted')
  rows+=1;null_AN+=AN is None
 done=db.execute('SELECT COUNT(*) FROM completed').fetchone()[0];counts=dict(db.execute('SELECT source,COUNT(*) FROM variants GROUP BY source'))
 if done!=2*len(intervals) or rows!=nuclear_receipt['rows'] or not any(k.startswith('exomes-') for k in counts) or not any(k.startswith('genomes-') for k in counts):raise ValueError('Public nuclear source extraction incomplete')
 results['population_counts']=dict(status='PASS_ACTUAL_PUBLIC_CDS_POPULATION_SOURCE_CONTROLS',rows=rows,completed_intervals=done,missing_AN_rows=null_AN,only_omitted_INFO=['vep','CSQ']);db.close()
 coverage=load(BASE/'medinter-population-gene-coverage-v1/resource-receipt.json')
 if coverage.get('downloaded_genes')!=38 or coverage.get('annotation_binding_updated') is not False:raise ValueError('Public gene coverage job not complete')
 results['coverage']=dict(status='PASS_PUBLIC_GENE_CONTEXT_COMPLETENESS',nuclear_genes=38,MT_gene_not_applicable=True,exact_base_depth_claim=False,version_mapping='v4 exome coverage / v3 genome coverage')
 output=SIDE/'actual-public-source-controls.json'
 if output.exists():raise ValueError('Existing actual-control receipt retained')
 artifacts=[artifact(path) for path in [nuclear/'public-population.sqlite',SIDE/'public-hcm39-exome-AN.sqlite',Path(mt['database']),Path(receipt['database']),SIDE/'public-hcm39-literature-review-queue.tsv']]
 save(output,dict(status='PASS_PUBLIC_SOURCE_SIDECAR_CONTROLS_NOT_CLINICAL_VALIDATION',checks=results,artifacts=artifacts,public_panel_genes=len(manifest['genes']),patient_inputs=False,annotation_binding_updated=False,automatic_pathogenicity=False,completed_utc=now().isoformat()))
 print(json.dumps(dict(status='PASS_PUBLIC_SOURCE_SIDECAR_CONTROLS_NOT_CLINICAL_VALIDATION',checks=results),indent=2))
if __name__=='__main__':
 import sys
 if '--write-ready-manifest' in sys.argv:ready_manifest()
 else:main()
