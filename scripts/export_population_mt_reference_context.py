#!/usr/bin/env python3
"""Local-only postannotation population/SCV/MT reference tables; never modifies primary outputs."""
import argparse
from collections import defaultdict
from contextlib import ExitStack
import csv
import ctypes
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import uuid
import zlib

EXPECTED_READY_SHA='465ee667cb78baddc3665d06adb15fb27ad54c47415ff14238fda5f0f42c3d85'
ABSOLUTE_DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
EXTRA=['pop_context_key','pop_gnomad_exact_match_status','pop_exome_AN_v41','pop_exome_AN_status','pop_SCV_map_status','pop_validated_AlleleIDs','pop_validated_VariationIDs','pop_gene_coverage_hgnc','pop_gene_coverage_scope','pop_master_source_row']
def now():return dt.datetime.now(dt.timezone.utc)
def load(path):return json.loads(Path(path).read_text())
def safe(path):
 path=Path(path)
 if any(p.is_symlink() for p in [path,*path.parents]):raise ValueError('Symlink path rejected')
 return path
def sha(path):
 digest=hashlib.sha256()
 with safe(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024**2),b''):digest.update(block)
 return digest.hexdigest()
def pin(path):
 path=safe(Path(path).absolute())
 return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
def verify(item):
 path=safe(Path(item['path']))
 if not path.is_file() or path.stat().st_size!=item['bytes'] or sha(path)!=item['sha256']:raise ValueError('Pinned source changed: '+str(path))
 return path
def dump(path,value):
 with safe(path).open('x') as stream:json.dump(value,stream,indent=2,ensure_ascii=False);stream.write('\n')
def text(path,mode):return gzip.open(path,mode,encoding='utf-8',newline='') if str(path).endswith('.gz') else Path(path).open(mode.replace('t',''),encoding='utf-8',newline='')
def writer(stack,path,fields):
 path.parent.mkdir(parents=True,exist_ok=True);stream=stack.enter_context(text(path,'xt'));result=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n');result.writeheader();return result
def canonical_chrom(value):
 value=str(value).removeprefix('chr')
 if value in ['M','MT']:return 'chrM'
 return 'chr'+value if value in [*map(str,range(1,23)),'X','Y'] else None
def allele_key(chrom,pos,ref,alt):
 chrom=canonical_chrom(chrom)
 if chrom is None or not str(pos).isdigit() or int(pos)<1 or not re.fullmatch('[ACGTN]+',ref or '') or not re.fullmatch('[ACGTN]+',alt or ''):return None
 return chrom,int(pos),ref,alt
def key_text(key):return ':'.join(map(str,key))
def output_guard(stage,cap,deadline):
 if now()>=deadline:raise TimeoutError('Original budget deadline')
 if cap and sum(p.stat().st_size for p in stage.rglob('*') if p.is_file())>cap:raise ValueError('Sidecar output cap reached; partial retained')
def hgnc(value):
 value=str(value or '');value='HGNC:'+value if value.isdigit() else value
 return value if re.fullmatch(r'HGNC:\d+',value) else None
def connection(path):
 db=sqlite3.connect(f'file:{safe(path)}?mode=ro',uri=True);db.row_factory=sqlite3.Row;db.execute('PRAGMA query_only=ON');return db
def promote(source,target):
 if os.name=='nt':os.rename(source,target);return
 libc=ctypes.CDLL(None,use_errno=True);rename=libc.renameat2;rename.argtypes=[ctypes.c_int,ctypes.c_char_p,ctypes.c_int,ctypes.c_char_p,ctypes.c_uint]
 if rename(-100,os.fsencode(source),-100,os.fsencode(target),1):raise OSError(ctypes.get_errno(),'Exclusive output promotion failed')
def ready_artifacts(rule):
 path=verify(rule['resource_ready_manifest'])
 if rule['resource_ready_manifest']['sha256']!=EXPECTED_READY_SHA:raise ValueError('Only exact ROOT-reviewed public resource manifest accepted')
 ready=load(path)
 if ready['status']!='AVAILABLE_VERIFIED_PUBLIC_REFERENCE_SIDECARS_NOT_BOUND' or ready['patient_inputs'] is not False or ready['annotation_binding_updated'] is not False:raise ValueError('Unexpected reference manifest scope')
 artifacts={Path(item['path']).name:item for item in ready['artifacts']}
 for item in artifacts.values():verify(item)
 return ready,artifacts
class Sources:
 def __init__(self,artifacts):
  self.pop=connection(artifacts['public-population.sqlite']['path']);self.an=connection(artifacts['public-hcm39-exome-AN.sqlite']['path']);self.cv=connection(artifacts['public-hcm39-clinvar-dossier-v2.sqlite']['path']);self.mt=connection(artifacts['public-mt-population-v2.sqlite']['path'])
  self.pop_keys=defaultdict(list);self.cv_keys=defaultdict(set)
  for row in self.pop.execute('SELECT rowid,chrom,pos,ref,alt FROM variants'):
   key=allele_key(row['chrom'],row['pos'],row['ref'],row['alt'])
   if key:self.pop_keys[key].append(row['rowid'])
  for row in self.cv.execute('SELECT variation_id,allele_id,chrom,pos,ref,alt FROM variants'):
   key=allele_key(row['chrom'],row['pos'],row['ref'],row['alt'])
   if key and str(row['allele_id']).isdigit() and str(row['variation_id']).isdigit():self.cv_keys[key].add((row['allele_id'],row['variation_id']))
 def counts(self,key):
  result=[]
  for rowid in self.pop_keys.get(key,[]):
   row=dict(self.pop.execute('SELECT * FROM variants WHERE rowid=?',(rowid,)).fetchone());value=row['info'];row['info']=zlib.decompress(value).decode() if isinstance(value,bytes) else value;result.append(row)
  return result
 def scv_map(self,key,reported):
  pairs=sorted(self.cv_keys.get(key,set()));reported=str(reported or '').strip()
  if reported in ['', '.', '-']:return ('EXACT_ALLELES_VALIDATED_SOURCE_ID_MAP' if pairs else 'NO_EXACT_SOURCE_ALLELES_NOT_ID_EQUIVALENCE'),pairs
  values=re.split(r'[&,;|]',reported)
  if not all(re.fullmatch(r'\d+',v) for v in values):return 'UNPARSEABLE_ALLELEID_NO_JOIN',[]
  allowed=set(values);matched=[p for p in pairs if p[0] in allowed]
  if {p[0] for p in matched}!=allowed:return 'ALLELEID_COORDINATE_CONFLICT_OR_SNAPSHOT_DIFFERENCE_NO_JOIN',[]
  return 'EXACT_ALLELES_AND_EXPLICIT_ALLELEID_TO_VARIATIONID',matched
 def close(self):
  for db in [self.pop,self.an,self.cv,self.mt]:db.close()
def coverage_sources(rule):
 result={}
 for item in rule['coverage_sources']:
  data=load(verify(item));gene=data['data']['gene']
  if hgnc(gene['hgnc_id'])!=item['hgnc_id'] or gene['reference_genome']!='GRCh38' or gene['symbol']!=item['gene']:raise ValueError('Coverage gene identity changed')
  result[item['hgnc_id']]=gene
 return result
def nuclear(master,stage,sources,coverage,public_genes,samples,deadline,output_cap=None):
 counters=dict(master_rows=0,selected_rows=0,samples={s:0 for s in samples},variants=0,SCV_map_conflict_rows=0);seen={};links=defaultdict(set);variant_ids=defaultdict(set);matched_genes=set()
 csv.field_size_limit(64*1024**2)
 with ExitStack() as stack,text(master,'rt') as inp:
  reader=csv.DictReader(inp,delimiter='\t');fields=reader.fieldnames or []
  if not {'sample','chrom','pos','ref','alt','GT','FT','DP','AD','GQ','sample_state','site_filter','hgnc_id','clinvar_allele_id'}<=set(fields) or set(EXTRA)&set(fields):raise ValueError('Primary master schema mismatch/already supplemented')
  cohort=writer(stack,stage/'cohort/nuclear-sample-rows.tsv.gz',fields+EXTRA)
  per={s:writer(stack,stage/f'samples/{s}/nuclear-reference-context.tsv.gz',fields+EXTRA) for s in samples}
  population=writer(stack,stage/'cohort/nuclear-population-long.tsv.gz',['variant_key','chrom','pos','ref','alt','dataset','source_version','match_status','AC','AN','AF','nhomalt','population_FILTER','ancestry_quality_INFO'])
  an_out=writer(stack,stage/'cohort/nuclear-exome-AN-v41.tsv.gz',['variant_key','chrom','pos','AN','match_status','capture_flags_json','scope'])
  mappings=writer(stack,stage/'cohort/ClinVar-validated-allele-map.tsv.gz',['variant_key','AlleleID','VariationID','scope'])
  observed=set()
  for number,row in enumerate(reader,1):
   if number%10000==1 and now()>=deadline:raise TimeoutError('Original budget deadline')
   counters['master_rows']+=1;sample=row['sample'];observed.add(sample)
   if sample not in per:raise ValueError('Unexpected primary master sample')
   key=allele_key(row['chrom'],row['pos'],row['ref'],row['alt']);gene=hgnc(row.get('hgnc_id'));scope_gene=gene in public_genes
   if key is None or key[0]=='chrM':continue
   if not (scope_gene or key in sources.pop_keys or key in sources.cv_keys):continue
   if key not in seen:
    pop=sources.counts(key);an=sources.an.execute('SELECT AN,raw_json FROM all_sites_AN WHERE chrom=? AND pos=?',key[:2]).fetchone();context=key_text(key)
    seen[key]=dict(key=context,pop_status='EXACT_GRCh38_ALLELES' if pop else 'NO_RECORD_IN_PUBLIC_CDS_SUBSET_NOT_AF_ZERO',AN=an['AN'] if an else '',AN_status='EXACT_POSITION_EXOME_AN_ONLY' if an else 'MISSING_NOT_ZERO_NOT_GENOME_OR_DEPTH')
    for dataset in ['exomes','genomes']:
     rows=[v for v in pop if v['source'].startswith(dataset+'-')]
     if len(rows)>1:raise ValueError('Duplicate population allele for dataset')
     value=rows[0] if rows else {}
     population.writerow(dict(variant_key=context,chrom=key[0],pos=key[1],ref=key[2],alt=key[3],dataset=dataset,source_version='4.1.1 public HCM39 CDS',match_status='EXACT_ALLELE' if value else 'MISSING_NOT_AF_ZERO',AC=value.get('AC',''),AN=value.get('AN',''),AF=value.get('AF',''),nhomalt=value.get('nhomalt',''),population_FILTER=value.get('filter',''),ancestry_quality_INFO=value.get('info','')))
    an_out.writerow(dict(variant_key=context,chrom=key[0],pos=key[1],AN=seen[key]['AN'],match_status=seen[key]['AN_status'],capture_flags_json=an['raw_json'] if an else '',scope='POSITION_CONTEXT_v4.1_EXOMES_ONLY_NOT_GENOMES_NOT_DEPTH'))
   context=seen[key];status,pairs=sources.scv_map(key,row['clinvar_allele_id'])
   if 'CONFLICT' in status:counters['SCV_map_conflict_rows']+=1
   for allele,variation in pairs:
    if (allele,variation) not in variant_ids[key]:mappings.writerow(dict(variant_key=context['key'],AlleleID=allele,VariationID=variation,scope='EXACT_GRCh38_ALLELES_SOURCE_ID_MAP_NOT_CLINICAL_CLASSIFICATION'));variant_ids[key].add((allele,variation))
    links[(variation,context['key'])].add(sample)
   cov=coverage.get(gene);valid_gene=cov is not None and (not row.get('gene') or row['gene']==cov['symbol'])
   if valid_gene:matched_genes.add(gene)
   extra=dict(pop_context_key=context['key'],pop_gnomad_exact_match_status=context['pop_status'],pop_exome_AN_v41=context['AN'],pop_exome_AN_status=context['AN_status'],pop_SCV_map_status=status,pop_validated_AlleleIDs='|'.join(sorted({p[0] for p in pairs})),pop_validated_VariationIDs='|'.join(sorted({p[1] for p in pairs})),pop_gene_coverage_hgnc=gene if valid_gene else '',pop_gene_coverage_scope='GENE_BIN_CONTEXT_ONLY_NOT_LOCUS_CALLABILITY' if valid_gene else 'MISSING_OR_IDENTITY_MISMATCH_NOT_ZERO',pop_master_source_row=number)
   combined=dict(row,**extra);cohort.writerow(combined);per[sample].writerow(combined);counters['selected_rows']+=1;counters['samples'][sample]+=1
   output_guard(stage,output_cap,deadline)
  if observed!=set(samples):raise ValueError('Master sample set differs from reviewed rule')
  scv=writer(stack,stage/'cohort/ClinVar-SCV-long.tsv.gz',['variant_key','VariationID','SCV_accession_version','submitter','date_evaluated','source_classification','raw_source_json','review_status'])
  sample_scv={s:writer(stack,stage/f'samples/{s}/SCV-reference-links.tsv.gz',['variant_key','VariationID','SCV_accession_version','scope']) for s in samples}
  used={v for v,key in links}
  by_variation=defaultdict(list)
  for variation,key in links:by_variation[variation].append(key)
  for row in sources.cv.execute('SELECT * FROM submissions'):
   if row['variation_id'] not in used:continue
   if not re.fullmatch(r'SCV\d+\.\d+',row['SCV']):raise ValueError('Unversioned SCV refused')
   for context in by_variation[row['variation_id']]:
    scv.writerow(dict(variant_key=context,VariationID=row['variation_id'],SCV_accession_version=row['SCV'],submitter=row['submitter'],date_evaluated=row['date_evaluated'],source_classification=row['classification'],raw_source_json=row['raw_json'],review_status='SOURCE_ASSERTION_ONLY_NO_PP5_NO_INDEPENDENT_CASE_COUNT'))
    for sample in links[(row['variation_id'],context)]:sample_scv[sample].writerow(dict(variant_key=context,VariationID=row['variation_id'],SCV_accession_version=row['SCV'],scope='VALIDATED_LOCAL_SOURCE_LINK_FOR_MANUAL_REVIEW_ONLY'))
    output_guard(stage,output_cap,deadline)
  citations=writer(stack,stage/'cohort/ClinVar-citation-links.tsv.gz',['variant_key','VariationID','citation_source','citation_id','reference_url','primary_article_read','independence_established'])
  for row in sources.cv.execute('SELECT variation_id,citation_source,citation_id FROM citations'):
   if row['variation_id'] not in used:continue
   for context in by_variation[row['variation_id']]:citations.writerow(dict(variant_key=context,VariationID=row['variation_id'],citation_source=row['citation_source'],citation_id=row['citation_id'],reference_url=f"https://pubmed.ncbi.nlm.nih.gov/{row['citation_id']}/" if row['citation_source']=='PubMed' and row['citation_id'].isdigit() else '',primary_article_read=False,independence_established=False))
  covout=writer(stack,stage/'cohort/gene-coverage-context-long.tsv.gz',['hgnc_id','gene','source_cohort','source_context','pos_as_returned','mean','median','over_10','over_20','over_30','scope'])
  for identity in sorted(matched_genes):
   value=coverage[identity]
   for cohort_name in ['exome','genome']:
    for record in value['coverage'][cohort_name]:covout.writerow(dict(hgnc_id=identity,gene=value['symbol'],source_cohort=cohort_name,source_context='gnomad_r4 API documented v4 exome / v3 genome',pos_as_returned=record['pos'],mean=record.get('mean'),median=record.get('median'),over_10=record.get('over_10'),over_20=record.get('over_20'),over_30=record.get('over_30'),scope='GENE_BIN_CONTEXT_ONLY_NO_EXACT_LOCUS_DEPTH_CLAIM'))
 counters.update(variants=len(seen),SCV_linked_variation_ids=len(used),coverage_genes=len(matched_genes),primary_ranking_changed=False)
 return counters
def mt_exports(entries,stage,sources,deadline,output_cap=None):
 summaries=[]
 for entry in entries:
  if now()>=deadline:raise TimeoutError('Original budget deadline')
  sample=entry['sample'];source=verify(entry['vcf']);receipt=load(verify(entry['receipt']))
  if receipt.get('status')!='PASS':raise ValueError('MT technical source receipt not PASS')
  if receipt['artifacts']['mitochondria.filtered.vcf.gz']['sha256']!=entry['vcf']['sha256']:raise ValueError('MT final VCF does not match successful source receipt')
  qcdir=stage/'samples'/sample/'MT-QC';qcdir.mkdir(parents=True,exist_ok=True)
  for index,item in enumerate(entry['qc']):
   target=qcdir/f'{index:02d}-{Path(item["path"]).name}.txt';shutil.copyfile(verify(item),target)
   if sha(target)!=item['sha256']:raise ValueError('Copied original MT QC differs')
  rows=0;multiallelic=0;contigs={};sample_header=None
  with ExitStack() as stack,gzip.open(source,'rt') as inp:
   out=writer(stack,stage/f'samples/{sample}/MT-population-reference.tsv.gz',['sample','CHROM','POS','ID','REF','ALT','QUAL','FILTER','raw_INFO','raw_FORMAT','raw_sample','GT','AD','AF','DP','FT','gnomAD_MT_match_status','gnomAD_MT_reference_json','Helix_match_status','Helix_reference_json','qc_scope'])
   for line in inp:
    if line.startswith('##contig='):
     match=re.search(r'ID=([^,>]+).*length=(\d+)',line)
     if match:contigs[match[1]]=int(match[2])
    if line.startswith('#CHROM'):
     header=line.rstrip('\n').split('\t');sample_header=header[9:]
     if sample_header!=[sample]:raise ValueError('MT source sample differs from reviewed sample')
     continue
    if line.startswith('#'):continue
    if now()>=deadline:raise TimeoutError('Original budget deadline')
    v=line.rstrip('\n').split('\t')
    if len(v)!=10 or sample_header is None or canonical_chrom(v[0])!='chrM' or contigs.get(v[0])!=16569:raise ValueError('MT rCRS contig/sample/format contract mismatch')
    values=dict(zip(v[8].split(':'),v[9].split(':')));key=allele_key(v[0],v[1],v[3],v[4]);pop=None;helix=None
    if key:
     record=sources.mt.execute('SELECT * FROM gnomad WHERE pos=? AND ref=? AND alt=?',key[1:]).fetchone();pop=dict(record) if record else None
     record=sources.mt.execute('SELECT * FROM helix WHERE pos=? AND ref=? AND alt=?',key[1:]).fetchone();helix=dict(record) if record else None
    else:multiallelic+=1
    missing='NO_EXACT_SOURCE_RECORD_NOT_AF_ZERO' if key else 'MULTIALLELIC_OR_SYMBOLIC_NO_FLATTEN_NO_JOIN'
    out.writerow(dict(sample=sample,CHROM=v[0],POS=v[1],ID=v[2],REF=v[3],ALT=v[4],QUAL=v[5],FILTER=v[6],raw_INFO=v[7],raw_FORMAT=v[8],raw_sample=v[9],GT=values.get('GT','.'),AD=values.get('AD','.'),AF=values.get('AF','.'),DP=values.get('DP','.'),FT=values.get('FT','.'),gnomAD_MT_match_status='EXACT_rCRS_ALLELES' if pop else missing,gnomAD_MT_reference_json=json.dumps(pop) if pop else '',Helix_match_status='EXACT_BIALLELIC_rCRS_SOURCE_AN_MISSING' if helix else missing,Helix_reference_json=json.dumps(helix) if helix else '',qc_scope='RAW_FINAL_FILTER_INFO_FORMAT_PRESERVED; linked NUMT/haplocheck source QC; research with limitations, no clinical LoD claim'));rows+=1
    output_guard(stage,output_cap,deadline)
  verify(entry['vcf']);summaries.append(dict(sample=sample,rows=rows,multiallelic_or_symbolic_not_split=multiallelic,raw_call_and_quality_fields_preserved=True))
 with ExitStack() as stack:
  combined=None
  for entry in entries:
   with text(stage/f'samples/{entry["sample"]}/MT-population-reference.tsv.gz','rt') as inp:
    reader=csv.DictReader(inp,delimiter='\t')
    if combined is None:combined=writer(stack,stage/'cohort/MT-population-reference-long.tsv.gz',reader.fieldnames)
    for row in reader:combined.writerow(row);output_guard(stage,output_cap,deadline)
 return summaries
def export(rule,deadline=None):
 deadline=min(ABSOLUTE_DEADLINE,dt.datetime.fromisoformat(rule['deadline_utc']),deadline or now()+dt.timedelta(hours=8))
 if now()>=deadline:raise TimeoutError('Original budget deadline')
 ready,artifacts=ready_artifacts(rule);workflow=safe(Path(rule['workflow_root']));receipts=[]
 for name in ['06_annotate','07_report']:
  path=workflow/'stages'/name/'workflow-receipt.json';value=load(path)
  if value.get('status')!='PASS' or not value.get('identity'):raise ValueError('Both primary stage 06/07 must PASS')
  receipts.append((pin(path),value))
 master=workflow/'stages/07_report/master_evidence.tsv.gz';expected=receipts[1][1]['artifacts'][master.name];master_pin=dict(path=str(master),bytes=expected['bytes'],sha256=expected['sha256']);verify(master_pin)
 public=load(verify(rule['public_panel']));genes={g['hgnc_id'] for g in public['genes']};coverage=coverage_sources(rule);samples=rule['expected_samples']
 if len(set(samples))!=len(samples) or any(not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]{0,99}',s) for s in samples):raise ValueError('Reviewed sample filenames invalid')
 output=safe(Path(rule['output']));state=safe(Path(rule['state_directory']))
 if output.exists() or workflow==output or workflow in output.parents or output==state or output in state.parents or state in output.parents:raise ValueError('Separate new immutable output and private state required')
 output.parent.mkdir(parents=True,exist_ok=True);stage=output.parent/('.'+output.name+'.partial-'+uuid.uuid4().hex);stage.mkdir(mode=0o700)
 sources=Sources(artifacts)
 try:
  summary=nuclear(master,stage,sources,coverage,genes,samples,deadline,rule['max_output_bytes']);mt=mt_exports(rule['mt_sources'],stage,sources,deadline,rule['max_output_bytes'])
  verify(master_pin)
  for item,value in receipts:verify(item)
  for item in rule['coverage_sources']:verify(item)
  for item in artifacts.values():verify(item)
  (stage/'README.md').write_text('Separate local reference context for manual research review. Primary tables, GT/FT/sample quality and ranking are unchanged. Population missingness is not AF=0. Coverage is gene-bin context. SCV assertions and citations are not independent PP5 or independent case evidence. MT raw FILTER/INFO/FORMAT/GT/AD/AF/DP are retained; no multiallelic splitting or clinical LoD claim.\n')
  outputs={str(path.relative_to(stage)):dict(bytes=path.stat().st_size,sha256=sha(path)) for path in stage.rglob('*') if path.is_file()}
  if sum(v['bytes'] for v in outputs.values())>rule['max_output_bytes']:raise ValueError('Sidecar output byte cap reached; partial retained')
  provenance=dict(status='PASS_LOCAL_POPULATION_MT_CONTEXT_EXPORT',source_master=master_pin,primary_receipts=[p for p,v in receipts],resource_ready_manifest=rule['resource_ready_manifest'],nuclear=summary,MT=mt,network_requests=0,primary_outputs_changed=False,automatic_ACMG=False,clinical_validation=False,created_utc=now().isoformat())
  dump(stage/'provenance.json',provenance);outputs['provenance.json']=dict(bytes=(stage/'provenance.json').stat().st_size,sha256=sha(stage/'provenance.json'))
  identity=hashlib.sha256(json.dumps(dict(master=master_pin,resources=rule['resource_ready_manifest'],code=rule['code_sha256'],MT=rule['mt_sources']),sort_keys=True).encode()).hexdigest()
  dump(stage/'workflow-receipt.json',dict(schema_version=1,status='PASS',identity=identity,stage='population-reference-context-v1',completed_utc=now().isoformat(),artifacts=outputs,metadata=provenance,note='Separate local reference sidecar; original primary rankings and files unchanged.'))
  promote(stage,output);return provenance
 finally:sources.close()
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rule',type=Path,required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
 print(json.dumps(export(load(a.rule)) if a.execute else dict(status='DRY_RUN_NO_PATIENT_OR_RESOURCE_READ'),indent=2))
