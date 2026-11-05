#!/usr/bin/env python3
"""One row-preserving local master overlay from completed, actually applied reference sidecars."""
import argparse
from collections import defaultdict
import csv
import ctypes
import datetime as dt
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import uuid

DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
TISSUE_FIELDS=['reference_context_hgnc','reference_context_scope','gtex_v11_gene_heart_medians_json','gtex_v11_transcript_match_status','gtex_v11_transcript_heart_medians_json','pext_gtex_v10_status','pext_source_transcript_match_status','pext_transcript_heart_expression_json','pext_allele_assignment_status','uniprot_gene_protein_context_json','protein_domain_allele_assignment_status','reference_context_index_sha256']
APPENDED=['allref_base_master_row_number','allref_tissue_status',*['allref_tissue_'+f for f in TISSUE_FIELDS], 'allref_population_status','allref_variant_key','allref_gnomad_exomes_summary_json','allref_gnomad_genomes_summary_json','allref_exome_AN_v41','allref_exome_AN_status','allref_SCV_map_status','allref_ClinVar_AlleleIDs','allref_ClinVar_VariationIDs','allref_SCV_reference_count','allref_SCV_classifications_json','allref_SCV_version_examples_json','allref_gene_coverage_summary_json','allref_source_versions_json','allref_population_evidence_ref','allref_SCV_evidence_ref','allref_MT_appendix_ref','allref_review_flags','allref_review_flag_details_json']
def now():return dt.datetime.now(dt.timezone.utc)
def safe(path):
 path=Path(path)
 if any(p.is_symlink() for p in [path,*path.parents]):raise ValueError('Symlink source/output rejected')
 return path
def load(path):return json.loads(safe(path).read_text())
def sha(path):
 value=hashlib.sha256()
 with safe(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024**2),b''):value.update(block)
 return value.hexdigest()
def pin(path):
 path=safe(Path(path).absolute());return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
def verify(item):
 path=safe(item['path'])
 if not path.is_file() or path.stat().st_size!=item['bytes'] or sha(path)!=item['sha256']:raise ValueError('Input artifact changed: '+str(path))
 return path
def dump(path,value):
 with safe(path).open('x') as stream:json.dump(value,stream,indent=2,ensure_ascii=False);stream.write('\n')
def text(path,mode='rt'):return gzip.open(path,mode,encoding='utf-8',newline='') if str(path).endswith('.gz') else Path(path).open(mode.replace('t',''),encoding='utf-8',newline='')
def canonical_tsv(values):
 stream=io.StringIO(newline='');csv.writer(stream,delimiter='\t',lineterminator='\n').writerow(values);return stream.getvalue().encode('utf-8')
def row_key(row,fields):
 if any(row.get(f) is None for f in fields):raise ValueError('Truncated original field vector')
 return hashlib.sha256(canonical_tsv([row[f] for f in fields])).hexdigest()
def hgnc(value):
 value=str(value or '');return 'HGNC:'+value if value.isdigit() else value
def variant_key(row):
 chrom=str(row.get('chrom','')).removeprefix('chr')
 if chrom not in [*map(str,range(1,23)),'X','Y'] or not str(row.get('pos','')).isdigit() or not re.fullmatch('[ACGTN]+',row.get('ref','')) or not re.fullmatch('[ACGTN]+',row.get('alt','')):return None
 return f"chr{chrom}:{int(row['pos'])}:{row['ref']}:{row['alt']}"
def promote(source,target):
 if os.name=='nt':os.rename(source,target);return
 libc=ctypes.CDLL(None,use_errno=True);fn=libc.renameat2;fn.argtypes=[ctypes.c_int,ctypes.c_char_p,ctypes.c_int,ctypes.c_char_p,ctypes.c_uint]
 if fn(-100,os.fsencode(source),-100,os.fsencode(target),1):raise OSError(ctypes.get_errno(),'Exclusive output promotion refused')
def limits(rule):
 if os.name!='nt':
  import resource
  resource.setrlimit(resource.RLIMIT_AS,(rule['memory_limit_bytes'],rule['memory_limit_bytes']))
  if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:rule['cpu_limit']]))
def artifact(root,receipt,name):
 expected=receipt['artifacts'][name];return verify(dict(path=str(Path(root)/name),bytes=expected['bytes'],sha256=expected['sha256']))
def applied_inputs(rule):
 workflow=Path(rule['workflow_root']);stage_receipts=[];documents=[]
 for name in ['06_annotate','07_report']:
  path=workflow/'stages'/name/'workflow-receipt.json';value=load(path)
  if value.get('status')!='PASS' or not value.get('identity'):raise ValueError('Both primary stages must PASS')
  stage_receipts.append(pin(path));documents.append(value)
 base=artifact(workflow/'stages/07_report',documents[1],'master_evidence.tsv.gz');transcripts=artifact(workflow/'stages/06_annotate',documents[0],'transcript_samples.tsv.gz')
 tissue_root=Path(rule['tissue_root']);tissue_path=tissue_root/'manifest.json';tissue=load(tissue_path)
 if tissue.get('status')!='PASS_SEPARATE_RESEARCH_REFERENCE_CONTEXT_EXPORT' or tissue.get('source_sha256')!=sha(transcripts) or tissue.get('network_requests')!=0 or tissue.get('original_pipeline_changed') is not False:raise ValueError('Tissue source must be actual matching local stage6 enrichment')
 if any(tissue['stage_receipts'].get(p['path'],{}).get('sha256')!=p['sha256'] for p in stage_receipts):raise ValueError('Tissue stage6/7 receipt pins differ')
 if tissue['context_index_sha256']!=rule['tissue_context_index']['sha256']:raise ValueError('Tissue context index changed')
 tissue_data=verify(tissue['artifact'])
 pop_root=Path(rule['population_root']);pop_path=pop_root/'workflow-receipt.json';pop=load(pop_path);provenance=pop['metadata']
 if pop.get('status')!='PASS' or pop.get('stage')!='population-reference-context-v1' or not pop.get('identity') or provenance.get('status')!='PASS_LOCAL_POPULATION_MT_CONTEXT_EXPORT' or provenance.get('network_requests')!=0 or provenance.get('primary_outputs_changed') is not False:raise ValueError('Population sidecar is not completed actual local enrichment')
 base_pin=pin(base)
 if provenance['source_master']!=base_pin or {p['path']:p['sha256'] for p in provenance['primary_receipts']}!={p['path']:p['sha256'] for p in stage_receipts}:raise ValueError('Population sidecar is for another master/stage context')
 if provenance['resource_ready_manifest']['sha256']!=rule['resource_ready_manifest']['sha256']:raise ValueError('Population source version/pin differs')
 for name,item in pop['artifacts'].items():verify(dict(path=str(pop_root/name),**item))
 pop_data=artifact(pop_root,pop,'cohort/nuclear-sample-rows.tsv.gz')
 return dict(base=base,base_pin=base_pin,tissue=tissue,tissue_data=tissue_data,pop=pop,pop_root=pop_root,pop_data=pop_data,stage_receipts=stage_receipts,upstream_receipts=[pin(tissue_path),pin(pop_path)])
def database(path):
 db=sqlite3.connect(path);db.row_factory=sqlite3.Row;db.executescript('PRAGMA cache_size=-65536;PRAGMA mmap_size=0;PRAGMA temp_store=FILE;CREATE TABLE row_context(kind TEXT,digest TEXT,n INTEGER,used INTEGER,payload TEXT,PRIMARY KEY(kind,digest));CREATE TABLE lookup(kind TEXT,key TEXT,payload TEXT,PRIMARY KEY(kind,key));')
 return db
def index_rows(db,kind,path,fields,required_extras,deadline):
 count=0
 with text(path) as stream,db:
  reader=csv.DictReader(stream,delimiter='\t')
  if not set(fields+required_extras)<=set(reader.fieldnames or []):raise ValueError('Upstream row schema cannot prove original-row identity: '+kind)
  for row in reader:
   if count%10000==0 and now()>=deadline:raise TimeoutError('Original deadline')
   digest=row_key(row,fields);payload=json.dumps({f:row[f] for f in required_extras},sort_keys=True,separators=(',',':'))
   old=db.execute('SELECT payload FROM row_context WHERE kind=? AND digest=?',(kind,digest)).fetchone()
   if old and old['payload']!=payload:raise ValueError('Identical source rows have contradictory supplemental evidence')
   db.execute('INSERT INTO row_context VALUES (?,?,1,0,?) ON CONFLICT(kind,digest) DO UPDATE SET n=n+1',(kind,digest,payload));count+=1
 return count
def index_long(db,kind,path,key_field,deadline=DEADLINE):
 with text(path) as stream,db:
  reader=csv.DictReader(stream,delimiter='\t')
  for number,row in enumerate(reader,1):
   if number%1000==1 and now()>=deadline:raise TimeoutError('Original deadline during applied-source indexing')
   key=row[key_field]
   if kind=='population':
    expected=variant_key(dict(chrom=row['chrom'],pos=row['pos'],ref=row['ref'],alt=row['alt']))
    if key!=expected or row['dataset'] not in ['exomes','genomes'] or not re.search(r'\b4\.1\.1\b',row['source_version']):raise ValueError('Population allele/dataset/version mismatch')
    key+='|'+row['dataset'];value={f:row[f] for f in ['source_version','match_status','AC','AN','AF','nhomalt','population_FILTER']}
   elif kind=='SCV':
    if not re.fullmatch(r'SCV\d+\.\d+',row['SCV_accession_version']):raise ValueError('Unversioned SCV')
    if not row['VariationID'].isdigit():raise ValueError('Invalid source VariationID')
    key=row['variant_key']+'|'+row['VariationID'];old=db.execute('SELECT payload FROM lookup WHERE kind=? AND key=?',(kind,key)).fetchone();value=json.loads(old['payload']) if old else dict(ids=[],classes=[],count=0)
    value['count']+=1
    if len(value['ids'])<5 and row['SCV_accession_version'] not in value['ids']:value['ids'].append(row['SCV_accession_version'])
    if row['source_classification'] not in value['classes']:value['classes'].append(row['source_classification'])
   elif kind=='map':
    if not row['AlleleID'].isdigit() or not row['VariationID'].isdigit() or not row['scope'].startswith('EXACT_GRCh38_ALLELES_SOURCE_ID_MAP'):raise ValueError('Unvalidated AlleleID/VariationID source map')
    old=db.execute('SELECT payload FROM lookup WHERE kind=? AND key=?',(kind,key)).fetchone();value=json.loads(old['payload']) if old else dict(pairs=[]);pair=[row['AlleleID'],row['VariationID']]
    if pair not in value['pairs']:value['pairs'].append(pair)
   elif kind=='AN':
    if row['scope']!='POSITION_CONTEXT_v4.1_EXOMES_ONLY_NOT_GENOMES_NOT_DEPTH' or not key.startswith(row['chrom']+':'+str(int(row['pos']))+':'):raise ValueError('Exome AN position/version/scope mismatch')
    value={f:row[f] for f in ['AN','match_status','scope']}
   elif kind=='coverage':
    key=row['hgnc_id'];old=db.execute('SELECT payload FROM lookup WHERE kind=? AND key=?',(kind,key)).fetchone();value=json.loads(old['payload']) if old else dict(hgnc_id=key,exome_bins=0,genome_bins=0,scope='GENE_BIN_ONLY_NOT_LOCUS_DEPTH',versions='v4 exome coverage / v3 genome coverage')
    cohort=row['source_cohort']
    if cohort not in ['exome','genome'] or row['scope']!='GENE_BIN_CONTEXT_ONLY_NO_EXACT_LOCUS_DEPTH_CLAIM':raise ValueError('Coverage scope invalid')
    value[cohort+'_bins']+=1
   else:raise ValueError('Unknown index kind')
   encoded=json.dumps(value,sort_keys=True,separators=(',',':'))
   old=db.execute('SELECT payload FROM lookup WHERE kind=? AND key=?',(kind,key)).fetchone()
   if kind=='population' and old and old['payload']!=encoded:raise ValueError('Duplicate contradictory population allele/dataset')
   db.execute('INSERT OR REPLACE INTO lookup VALUES (?,?,?)',(kind,key,encoded))
def take(db,kind,digest):
 row=db.execute('SELECT n,used,payload FROM row_context WHERE kind=? AND digest=?',(kind,digest)).fetchone()
 if row is None:return None
 if row['used']>=row['n']:raise ValueError('Supplemental row multiplicity exhausted')
 db.execute('UPDATE row_context SET used=used+1 WHERE kind=? AND digest=?',(kind,digest));return json.loads(row['payload'])
def lookup(db,kind,key):
 row=db.execute('SELECT payload FROM lookup WHERE kind=? AND key=?',(kind,key)).fetchone();return json.loads(row['payload']) if row else None
def extras(db,row,digest,number,genes,versions):
 tissue=take(db,'tissue',digest);population=take(db,'population',digest);gene=hgnc(row.get('hgnc_id'));key=variant_key(row);in_scope=gene in genes
 if in_scope and tissue is None:raise ValueError('Public gene row missing exact stage6-to-stage7 tissue correspondence')
 if in_scope and key is not None and population is None:raise ValueError('Public gene row missing actually applied population context')
 result={f:'' for f in APPENDED};result.update(allref_base_master_row_number=number,allref_tissue_status='VERIFIED_FULL_ORIGINAL_ROW_MATCH' if tissue else 'OUTSIDE_APPROVED_PUBLIC39_CONTEXT',allref_population_status='VERIFIED_FULL_ORIGINAL_ROW_MATCH' if population else 'OUTSIDE_PUBLIC_SUBSET_OR_UNSUPPORTED_ALLELE_NOT_AF_ZERO',allref_source_versions_json=versions,allref_MT_appendix_ref='applied-evidence/MT-appendix-index.json')
 flags=[];details={}
 if tissue:
  for field in TISSUE_FIELDS:result['allref_tissue_'+field]=tissue[field]
  if any('VERSION_MISMATCH' in tissue[f] for f in ['gtex_v11_transcript_match_status','pext_source_transcript_match_status']):flags.append('TRANSCRIPT_VERSION_CONTEXT_REVIEW')
 if population:
  if key!=population['pop_context_key']:raise ValueError('Population exact allele key differs from original row')
  an_context=lookup(db,'AN',key)
  if an_context is None or an_context['AN']!=population['pop_exome_AN_v41'] or an_context['match_status']!=population['pop_exome_AN_status']:raise ValueError('Applied row AN differs from actual separate exome position context')
  result.update(allref_variant_key=key,allref_exome_AN_v41=population['pop_exome_AN_v41'],allref_exome_AN_status=population['pop_exome_AN_status'],allref_SCV_map_status=population['pop_SCV_map_status'],allref_ClinVar_AlleleIDs=population['pop_validated_AlleleIDs'],allref_ClinVar_VariationIDs=population['pop_validated_VariationIDs'],allref_population_evidence_ref=json.dumps(dict(path='applied-evidence/population/cohort/nuclear-population-long.tsv.gz',variant_key=key)),allref_SCV_evidence_ref=json.dumps(dict(path='applied-evidence/population/cohort/ClinVar-SCV-long.tsv.gz',variant_key=key)))
  for dataset,field in [('exomes','allref_gnomad_exomes_summary_json'),('genomes','allref_gnomad_genomes_summary_json')]:
   value=lookup(db,'population',key+'|'+dataset)
   if value is None:raise ValueError('Applied population row has no matching actual long-table record')
   result[field]=json.dumps(value,sort_keys=True,separators=(',',':'))
   if value['match_status']=='EXACT_ALLELE' and value['population_FILTER'] not in ['PASS','.','']:flags.append('GNOMAD_SOURCE_FILTER_REVIEW')
  if 'CONFLICT' in population['pop_SCV_map_status']:flags.append('ALLELEID_SOURCE_MAP_REVIEW')
  allele_ids=set(filter(None,population['pop_validated_AlleleIDs'].split('|')));variation_ids=set(filter(None,population['pop_validated_VariationIDs'].split('|')));mapping=lookup(db,'map',key) or {'pairs':[]};pairs=[p for p in mapping['pairs'] if p[0] in allele_ids and p[1] in variation_ids]
  if {p[0] for p in pairs}!=allele_ids or {p[1] for p in pairs}!=variation_ids:raise ValueError('Row SCV IDs do not agree with exact source AlleleID/VariationID mapping')
  upstream_id=str(row.get('clinvar_allele_id','')).strip()
  if allele_ids and upstream_id not in ['','.','-']:
   parsed=re.split(r'[&,;|]',upstream_id)
   if not all(v.isdigit() for v in parsed) or set(parsed)!=allele_ids:raise ValueError('Upstream AlleleID was conflated with another identifier')
  scvs=[value for value in (lookup(db,'SCV',key+'|'+identity) for identity in sorted(variation_ids)) if value]
  scv=dict(count=sum(v['count'] for v in scvs),classes=sorted({c for v in scvs for c in v['classes']}),ids=sorted({i for v in scvs for i in v['ids']})[:5]) if scvs else None
  if scv:
   result.update(allref_SCV_reference_count=scv['count'],allref_SCV_classifications_json=json.dumps(scv['classes']),allref_SCV_version_examples_json=json.dumps(scv['ids']))
   if len(scv['classes'])>1:flags.append('SCV_CLASSIFICATIONS_DIFFER_CONDITION_AND_DATE_REVIEW');details['SCV_scope']='Distinct source assertion texts; not an automatic conflict resolution or independent PP5/case evidence'
  cov=lookup(db,'coverage',population['pop_gene_coverage_hgnc']) if population['pop_gene_coverage_hgnc'] else None
  if cov:result['allref_gene_coverage_summary_json']=json.dumps(cov,sort_keys=True,separators=(',',':'))
 result['allref_review_flags']=';'.join(sorted(set(flags)));result['allref_review_flag_details_json']=json.dumps(details,sort_keys=True,separators=(',',':'));return result
def join(rule,deadline=None):
 deadline=min(DEADLINE,dt.datetime.fromisoformat(rule['deadline_utc']),deadline or now()+dt.timedelta(hours=8));limits(rule)
 if now()>=deadline:raise TimeoutError('Original budget deadline')
 source=applied_inputs(rule);base=source['base'];public=load(verify(rule['public_panel']));genes={g['hgnc_id'] for g in public['genes']};ready=load(verify(rule['resource_ready_manifest']));versions=json.dumps(ready['source_versions'],sort_keys=True,separators=(',',':'))
 output=safe(rule['output']);state=safe(rule['state_directory'])
 if output.exists() or output==Path(rule['workflow_root']) or Path(rule['workflow_root']) in output.parents:raise ValueError('A new separate combined output is required')
 output.parent.mkdir(parents=True,exist_ok=True);stage=output.parent/('.'+output.name+'.partial-'+uuid.uuid4().hex);stage.mkdir(mode=0o700);state.mkdir(parents=True,exist_ok=True)
 index=state/('row-index-'+uuid.uuid4().hex+'.sqlite');db=database(index);csv.field_size_limit(64*1024**2)
 try:
  with text(base) as stream:original_fields=next(csv.reader(stream,delimiter='\t'))
  if any(f.startswith('allref_') for f in original_fields) or len(set(original_fields))!=len(original_fields):raise ValueError('Baseline original fields collide with overlay')
  pop_fields=['pop_context_key','pop_gnomad_exact_match_status','pop_exome_AN_v41','pop_exome_AN_status','pop_SCV_map_status','pop_validated_AlleleIDs','pop_validated_VariationIDs','pop_gene_coverage_hgnc','pop_gene_coverage_scope']
  loaded=dict(tissue=index_rows(db,'tissue',source['tissue_data'],original_fields,TISSUE_FIELDS,deadline),population=index_rows(db,'population',source['pop_data'],original_fields,pop_fields,deadline))
  pop_root=source['pop_root'];index_long(db,'population',pop_root/'cohort/nuclear-population-long.tsv.gz','variant_key',deadline);index_long(db,'AN',pop_root/'cohort/nuclear-exome-AN-v41.tsv.gz','variant_key',deadline);index_long(db,'map',pop_root/'cohort/ClinVar-validated-allele-map.tsv.gz','variant_key',deadline);index_long(db,'SCV',pop_root/'cohort/ClinVar-SCV-long.tsv.gz','variant_key',deadline);index_long(db,'coverage',pop_root/'cohort/gene-coverage-context-long.tsv.gz','hgnc_id',deadline)
  if index.stat().st_size>rule['max_sqlite_bytes']:raise ValueError('Private join index size cap')
  target=stage/'master_evidence.enriched.tsv.gz';base_digest=hashlib.sha256();output_digest=hashlib.sha256();header_bytes=canonical_tsv(original_fields);base_digest.update(header_bytes);output_digest.update(header_bytes);count=0;matched={'tissue':0,'population':0};flags=0
  with text(base) as inp,text(target,'xt') as out,db:
   reader=csv.DictReader(inp,delimiter='\t');writer=csv.DictWriter(out,original_fields+APPENDED,delimiter='\t',lineterminator='\n');writer.writeheader()
   for count,row in enumerate(reader,1):
    if count%1000==1:
     if now()>=deadline:raise TimeoutError('Original deadline')
     if index.stat().st_size>rule['max_sqlite_bytes'] or target.stat().st_size>rule['max_output_bytes']:raise ValueError('Join index/output cap')
     if shutil.disk_usage(stage).free<rule['min_free_bytes']:raise ValueError('Keep analysis disk reserve')
    raw=canonical_tsv([row[f] for f in original_fields]);base_digest.update(raw);new=extras(db,row,row_key(row,original_fields),count,genes,versions);combined=dict(row,**new);writer.writerow(combined);output_digest.update(canonical_tsv([combined[f] for f in original_fields]));flags+=bool(new['allref_review_flags'])
    matched['tissue']+=new['allref_tissue_status']=='VERIFIED_FULL_ORIGINAL_ROW_MATCH';matched['population']+=new['allref_population_status']=='VERIFIED_FULL_ORIGINAL_ROW_MATCH'
  orphaned=dict(db.execute('SELECT kind,SUM(n-used) FROM row_context GROUP BY kind'))
  if any(orphaned.values()) or base_digest.hexdigest()!=output_digest.hexdigest():raise ValueError('Original row/order/field or supplemental multiplicity invariant failed')
  auxiliary=[];mt=[]
  for name,item in source['pop']['artifacts'].items():
   if now()>=deadline:raise TimeoutError('Original deadline during applied-result auxiliary copy')
   nuclear_aux={'cohort/nuclear-population-long.tsv.gz','cohort/nuclear-exome-AN-v41.tsv.gz','cohort/ClinVar-validated-allele-map.tsv.gz','cohort/ClinVar-SCV-long.tsv.gz','cohort/ClinVar-citation-links.tsv.gz','cohort/gene-coverage-context-long.tsv.gz'}
   mt_aux=name=='cohort/MT-population-reference-long.tsv.gz' or re.fullmatch(r'samples/[A-Za-z0-9][A-Za-z0-9_.-]{0,99}/MT-population-reference.tsv.gz',name) or re.fullmatch(r'samples/[A-Za-z0-9][A-Za-z0-9_.-]{0,99}/MT-QC/\d\d-[A-Za-z0-9_.-]+\.txt',name)
   selected=name in nuclear_aux or mt_aux
   if not selected:continue
   original=artifact(pop_root,source['pop'],name);relative=Path('applied-evidence/population')/name;destination=stage/relative;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(original,destination)
   if sha(destination)!=item['sha256']:raise ValueError('Copied auxiliary artifact changed')
   entry=dict(role='MT_APPENDIX' if mt_aux else 'APPLIED_NUCLEAR_REFERENCE_LONG_TABLE',path=str(output/relative),relative_path=relative.as_posix(),source_relative_path=name,source_upstream_receipt=source['upstream_receipts'][1],bytes=item['bytes'],sha256=item['sha256']);auxiliary.append(entry)
   if entry['role']=='MT_APPENDIX':mt.append(entry)
  (stage/'applied-evidence').mkdir(exist_ok=True)
  dump(stage/'applied-evidence/MT-appendix-index.json',dict(status='SEPARATE_MT_APPENDIX_NOT_NUCLEAR_GENOTYPE_ROWS',artifacts=mt,no_automatic_classification=True));dump(stage/'auxiliary-index.json',dict(artifacts=auxiliary,source_evidence_keys='variant_key; validated AlleleID/VariationID/SCV; HGNC gene-bin context',no_external_requests=True))
  (stage/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Combined annotation master</title><h1>Combined annotation master</h1><p>Original nuclear rows and ranking retained. Supplemental evidence is for review.</p><ul><li><a href="master_evidence.enriched.tsv.gz">One enriched master</a></li><li><a href="auxiliary-index.json">Applied source long tables and keys</a></li><li><a href="applied-evidence/MT-appendix-index.json">Separate MT appendix</a></li></ul>')
  if sum(p.stat().st_size for p in stage.rglob('*') if p.is_file())>rule['max_output_bytes']:raise ValueError('Combined output cap')
  verify(source['base_pin'])
  for item in [*source['stage_receipts'],*source['upstream_receipts']]:verify(item)
  enriched=dict(path=str(output/target.name),bytes=target.stat().st_size,sha256=sha(target))
  manifest=dict(schema_version=1,status='PASS_COMBINED_ANNOTATION_MASTER',base_master=source['base_pin'],enriched_master=enriched,original_fields=original_fields,appended_fields=APPENDED,row_count=count,invariant_digest=base_digest.hexdigest(),invariant_digest_contract='SHA256 of csv.writer(delimiter=TAB,lineterminator=LF) original ordered header and field vectors, UTF-8',original_rows_fields_order_preserved=True,ranking_eligibility_unchanged=True,stage_receipts=source['stage_receipts'],upstream_receipts=source['upstream_receipts'],producer_code_sha256=rule['code_sha256'],context_rows_loaded=loaded,context_rows_matched=matched,orphan_context_rows=orphaned,review_flag_rows=flags,auxiliary_artifacts=auxiliary,MT_rows_added_to_nuclear_master=0,network_requests=0,created_utc=now().isoformat())
  dump(stage/'manifest.json',manifest);artifacts={p.relative_to(stage).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in stage.rglob('*') if p.is_file()};identity=hashlib.sha256(json.dumps(dict(base=source['base_pin'],upstream=source['upstream_receipts'],code=rule['code_sha256']),sort_keys=True).encode()).hexdigest()
  dump(stage/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='combined-annotation-master-v1',identity=identity,completed_utc=now().isoformat(),artifacts=artifacts,metadata=manifest));promote(stage,output);return manifest
 finally:db.close()
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rule',type=Path,required=True);p.add_argument('--execute',action='store_true');a=p.parse_args();print(json.dumps(join(load(a.rule)) if a.execute else dict(status='DRY_RUN_NO_PATIENT_DATA_READ'),indent=2))
