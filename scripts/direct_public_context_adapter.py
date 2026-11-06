#!/usr/bin/env python3
"""Public-source index and pure row/allele context for a separately reviewed producer."""
import argparse
from collections import defaultdict
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sqlite3
import zlib

READY_SHA='465ee667cb78baddc3665d06adb15fb27ad54c47415ff14238fda5f0f42c3d85'
def module(name):
 p=Path(__file__).resolve().parent/(name+'.py');spec=importlib.util.spec_from_file_location('direct_'+name,p);value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value
base=module('embed_vcf_evidence');tissue=module('export_tissue_protein_context');joined=module('join_combined_annotation_master');splice=module('spliceai_packet_context')
APPENDED=joined.APPENDED
SPLICE_FIELDS=splice.SPLICE_FIELDS
def pack(value):return zlib.compress(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode(),3)
def unpack(value):return json.loads(zlib.decompress(value))
def connect(path):
 db=sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro&immutable=1',uri=True);db.row_factory=sqlite3.Row;db.execute('PRAGMA cache_size=-16384');return db
def hgnc(value):return 'HGNC:'+str(value) if str(value).isdigit() else str(value or '')
def key(row):return joined.variant_key(row)
def build(ready_path,pop_rule_path,tissue_path,output):
 if base.sha(ready_path)!=READY_SHA:raise ValueError('Exact approved source release required')
 ready=base.load(ready_path);rule=base.load(pop_rule_path);context_pin=base.pin(tissue_path)
 if context_pin['sha256']!='231dba7b57acbb5c07bf16779e6af34bb48ab945b016017a0554b9c418e942c6':raise ValueError('Exact approved tissue context required')
 sources={Path(p['path']).name:p for p in ready['artifacts']};used=[sources[n] for n in ['public-population.sqlite','public-hcm39-exome-AN.sqlite','public-hcm39-clinvar-dossier-v2.sqlite']]
 for item in used:base.verify(item)
 output=base.safe(output)
 if output.exists():raise ValueError('Fresh public derived index required')
 output.mkdir(mode=0o700,parents=True);path=output/'public-context.sqlite';db=sqlite3.connect(path);db.executescript('PRAGMA cache_size=-65536;CREATE TABLE reference(kind TEXT,key TEXT,subkey TEXT,payload BLOB,PRIMARY KEY(kind,key,subkey));CREATE TABLE metadata(key TEXT PRIMARY KEY,payload BLOB);')
 counts=defaultdict(int)
 def put(kind,identity,subkey,value):
  raw=pack(value);old=db.execute('SELECT payload FROM reference WHERE kind=? AND key=? AND subkey=?',(kind,identity,subkey)).fetchone()
  if old and old[0]!=raw:raise ValueError('Conflicting public source identity')
  db.execute('INSERT OR IGNORE INTO reference VALUES(?,?,?,?)',(kind,identity,subkey,raw));counts[kind]+=1
 pop=connect(used[0]['path']);an=connect(used[1]['path']);cv=connect(used[2]['path'])
 try:
  with db:
   for row in pop.execute('SELECT * FROM variants'):
    row=dict(row);identity=key(row)
    if identity is None:continue
    dataset=row['source'].split('-',1)[0]
    if dataset not in ['exomes','genomes']:raise ValueError('Unknown public population dataset')
    info=row['info'];info=zlib.decompress(info).decode() if isinstance(info,bytes) else info
    put('population',identity,dataset,dict(variant_key=identity,chrom=row['chrom'],pos=str(row['pos']),ref=row['ref'],alt=row['alt'],dataset=dataset,source_version='4.1.1 public HCM39 CDS',match_status='EXACT_ALLELE',AC=str(row['AC']) if row['AC'] is not None else '',AN=str(row['AN']) if row['AN'] is not None else '',AF=str(row['AF']) if row['AF'] is not None else '',nhomalt=str(row['nhomalt']) if row['nhomalt'] is not None else '',population_FILTER=row['filter'],ancestry_quality_INFO=info))
   for row in an.execute('SELECT * FROM all_sites_AN'):
    put('AN',row['chrom']+':'+str(row['pos']),'',dict(chrom=row['chrom'],pos=str(row['pos']),AN=str(row['AN']),match_status='EXACT_POSITION_EXOME_AN_ONLY',capture_flags_json=row['raw_json'],scope='POSITION_CONTEXT_v4.1_EXOMES_ONLY_NOT_GENOMES_NOT_DEPTH'))
   for row in cv.execute('SELECT * FROM variants'):
    identity=key(dict(row))
    if identity is None:continue
    if not row['variation_id'].isdigit() or not row['allele_id'].isdigit():raise ValueError('Invalid source ID')
    put('map',identity,row['allele_id']+'|'+row['variation_id'],dict(variant_key=identity,AlleleID=row['allele_id'],VariationID=row['variation_id'],scope='EXACT_GRCh38_ALLELES_SOURCE_ID_MAP_NOT_CLINICAL_CLASSIFICATION'))
   for number,row in enumerate(cv.execute('SELECT * FROM submissions'),1):
    if not re.fullmatch(r'SCV[0-9]+\.[0-9]+',row['SCV']):raise ValueError('Versioned SCV required')
    put('SCV',row['variation_id'],row['SCV']+'|'+str(number).zfill(12),dict(VariationID=row['variation_id'],SCV_accession_version=row['SCV'],submitter=row['submitter'],date_evaluated=row['date_evaluated'],source_classification=row['classification'],raw_source_json=row['raw_json'],review_status='SOURCE_ASSERTION_ONLY_NO_PP5_NO_INDEPENDENT_CASE_COUNT'))
   for row in cv.execute('SELECT * FROM citations'):
    put('citation',row['variation_id'],row['citation_source']+'|'+row['citation_id'],dict(VariationID=row['variation_id'],citation_source=row['citation_source'],citation_id=row['citation_id'],raw_source_json=row['raw_json'],reference_url='https://pubmed.ncbi.nlm.nih.gov/'+row['citation_id']+'/' if row['citation_source']=='PubMed' and row['citation_id'].isdigit() else '',primary_article_read=False,independence_established=False))
   coverage={}
   for item in rule['coverage_sources']:
    raw=base.load(base.verify(item));gene=raw['data']['gene']
    if hgnc(gene['hgnc_id'])!=item['hgnc_id'] or gene['symbol']!=item['gene'] or gene['reference_genome']!='GRCh38':raise ValueError('Public coverage gene identity changed')
    coverage[item['hgnc_id']]=gene;put('coverage',item['hgnc_id'],'',gene)
   public=base.public_bundle(dict(tissue_context_index=context_pin,max_header_inflated_bytes=64*1024**2));gene_bundle=dict(v=1,scope='PUBLIC_HGNC_GENE_BIN_CONTEXT_NOT_LOCUS_DEPTH',genes=coverage)
   if sum(len(json.dumps(v,ensure_ascii=False).encode()) for v in [public,gene_bundle])>64*1024**2:raise ValueError('Total inflated public header cap')
   summaries={identity:dict(hgnc_id=identity,scope='GENE_BIN_ONLY_NOT_LOCUS_DEPTH',versions='v4 exome coverage / v3 genome coverage',exome_bins=len(gene['coverage']['exome']),genome_bins=len(gene['coverage']['genome'])) for identity,gene in coverage.items()}
   for name,value in [('tissue',base.load(tissue_path)),('source_versions',ready['source_versions']),('public_bundle',public),('gene_bundle',gene_bundle),('coverage_summaries',summaries)]:db.execute('INSERT INTO metadata VALUES(?,?)',(name,pack(value)))
  db.execute('PRAGMA optimize');db.close();metadata=dict(schema_version=1,status='READY_PUBLIC_DIRECT_CONTEXT_INDEX_NOT_APPLIED',source_ready_manifest=base.pin(ready_path),source_databases=used,tissue_context_index=context_pin,population_rule=base.pin(pop_rule_path),counts=dict(counts),index=base.pin(path),public_header_sha256=hashlib.sha256(base.encode_packet(public).encode()).hexdigest(),gene_header_sha256=hashlib.sha256(base.encode_packet(gene_bundle).encode()).hexdigest(),patient_inputs=False,network_requests=0)
  base.dump(output/'manifest.json',metadata);return metadata
 finally:
  for connection in [pop,an,cv]:connection.close()
class Context:
 def __init__(self,index_manifest,spliceai_manifest=None,verify_index=True,expected_source_vcf=None):
  self.manifest=base.load(index_manifest)
  if self.manifest.get('status')!='READY_PUBLIC_DIRECT_CONTEXT_INDEX_NOT_APPLIED' or self.manifest['source_ready_manifest']['sha256']!=READY_SHA:raise ValueError('Exact prepared public context index required')
  self.db=connect(base.verify(self.manifest['index']) if verify_index else self.manifest['index']['path']);self.tissue=unpack(self.db.execute('SELECT payload FROM metadata WHERE key=?',('tissue',)).fetchone()[0]);self.versions=json.dumps(unpack(self.db.execute('SELECT payload FROM metadata WHERE key=?',('source_versions',)).fetchone()[0]),sort_keys=True,separators=(',',':'));self.context_sha=self.manifest['tissue_context_index']['sha256'];self.coverage=unpack(self.db.execute('SELECT payload FROM metadata WHERE key=?',('coverage_summaries',)).fetchone()[0]);self.public_genes=set(self.tissue['genes']);self.last_key=None;self.last=None;self.splices={};self.splice_source='';self.splice_manifest=None
  if spliceai_manifest:
   self.splice_manifest=base.load(spliceai_manifest)
   receipt=base.load(Path(spliceai_manifest).parent/'workflow-receipt.json');manifest_pin=base.pin(spliceai_manifest)
   if self.splice_manifest.get('status')!='PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT' or receipt.get('status')!='PASS' or receipt.get('stage')!='spliceai-reference-context-v1' or not re.fullmatch('[a-f0-9]{64}',receipt.get('identity','')):raise ValueError('Actual SpliceAI completion receipt required')
   for name,item in [('manifest.json',manifest_pin),('predictions.tsv.gz',self.splice_manifest['predictions'])]:
    if any(receipt['artifacts'].get(name,{}).get(k)!=item[k] for k in ['bytes','sha256']):raise ValueError('SpliceAI completion artifact differs')
   if self.splice_manifest['distance']!=500 or self.splice_manifest['mask']!=0 or self.splice_manifest['gene_scope']!=splice.SCOPE or len(self.splice_manifest['models'])!=5 or self.splice_manifest.get('network_requests')!=0 or self.splice_manifest.get('primary_outputs_changed') is not False:raise ValueError('Reviewed SpliceAI model/scope parameters required')
   if expected_source_vcf is not None and self.splice_manifest['source_vcf']!=expected_source_vcf:raise ValueError('SpliceAI was calculated for a different accepted input')
   with gzip.open(base.verify(self.splice_manifest['predictions']),'rt') as f:
    import csv
    rows=0
    for row in csv.DictReader(f,delimiter='\t'):
     source=json.loads(row['source_json'])
     if source['version']!='1.3.1' or source['distance']!=500 or source['mask']!=0 or source['model_sha256']!=[p['sha256'] for p in self.splice_manifest['models']] or source['annotation_sha256']!=self.splice_manifest['annotation']['sha256'] or source['reference_sha256']!=self.splice_manifest['reference']['sha256']:raise ValueError('Prediction row model provenance differs')
     if row['variant_key']!=splice.key(row['chrom'],row['pos'],row['ref'],row['alt']):raise ValueError('Prediction allele key differs')
     self.splices.setdefault(row['variant_key'],[]).append(row);self.splice_source=row['source_json'];rows+=1
    if rows!=self.splice_manifest['prediction_rows'] or not self.splices:raise ValueError('SpliceAI prediction count differs/empty')
 def get(self,kind,identity):return [unpack(r[0]) for r in self.db.execute('SELECT payload FROM reference WHERE kind=? AND key=? ORDER BY subkey',(kind,identity))]
 def allele_evidence(self,chrom,pos,ref,alt):
  identity=key(dict(chrom=chrom,pos=str(pos),ref=ref,alt=alt))
  if identity==self.last_key and self.last is not None:return self.last
  population=self.get('population',identity) if identity else [];mapping=self.get('map',identity) if identity else [];ans=self.get('AN','chr'+str(chrom).removeprefix('chr')+':'+str(int(pos))) if identity else [];scvs=[];citations=[]
  for variation in sorted({r['VariationID'] for r in mapping}):
   scvs.extend(dict(r,variant_key=identity) for r in self.get('SCV',variation));citations.extend(dict(r,variant_key=identity) for r in self.get('citation',variation))
  result=dict(population=population,exome_AN=[dict(r,variant_key=identity) for r in ans],allele_map=mapping,SCV=scvs,citations=citations);self.last_key=identity;self.last=result;return result
 def augment(self,row,ordinal):
  result={f:'' for f in APPENDED};gene=hgnc(row.get('hgnc_id'));t=tissue.supplement(row,self.tissue,ordinal);identity=key(row);evidence=self.allele_evidence(row['chrom'],row['pos'],row['ref'],row['alt']);flags=[];details={};result.update(allref_base_master_row_number=str(ordinal),allref_tissue_status='APPLIED_DIRECT_EXACT_HGNC_AND_TRANSCRIPT_VERSION_STATUS' if t else 'OUTSIDE_APPROVED_PUBLIC39_CONTEXT',allref_population_status='APPLIED_DIRECT_EXACT_ALLELE_SOURCE_CONTEXT' if identity and (gene in self.tissue['genes'] or evidence['population'] or evidence['allele_map'] or evidence['exome_AN']) else 'OUTSIDE_PUBLIC_SUBSET_OR_UNSUPPORTED_ALLELE_NOT_AF_ZERO',allref_source_versions_json=self.versions,allref_MT_appendix_ref='SEPARATE_MT_APPENDIX_NOT_NUCLEAR_GENOTYPE_ROWS')
  if t:
   for field in joined.TISSUE_FIELDS:result['allref_tissue_'+field]=self.context_sha if field=='reference_context_index_sha256' else str(t[field])
   if any('VERSION_MISMATCH' in t[f] for f in ['gtex_v11_transcript_match_status','pext_source_transcript_match_status']):flags.append('TRANSCRIPT_VERSION_CONTEXT_REVIEW')
  if identity and result['allref_population_status']=='APPLIED_DIRECT_EXACT_ALLELE_SOURCE_CONTEXT':
   result.update(allref_variant_key=identity,allref_population_evidence_ref='INFO/WGS_EVIDENCE/evidence/population',allref_SCV_evidence_ref='INFO/WGS_EVIDENCE/evidence/SCV')
   for dataset,field in [('exomes','allref_gnomad_exomes_summary_json'),('genomes','allref_gnomad_genomes_summary_json')]:
    values=[r for r in evidence['population'] if r['dataset']==dataset]
    if len(values)>1:raise ValueError('Ambiguous exact population dataset')
    value={k:values[0][k] for k in ['source_version','match_status','AC','AN','AF','nhomalt','population_FILTER']} if values else dict(source_version='4.1.1 public HCM39 CDS',match_status='MISSING_NOT_AF_ZERO',AC='',AN='',AF='',nhomalt='',population_FILTER='');result[field]=json.dumps(value,sort_keys=True,separators=(',',':'))
    if value['match_status']=='EXACT_ALLELE' and value['population_FILTER'] not in ['PASS','.','']:flags.append('GNOMAD_SOURCE_FILTER_REVIEW')
   an=evidence['exome_AN'];result.update(allref_exome_AN_v41=an[0]['AN'] if an else '',allref_exome_AN_status=an[0]['match_status'] if an else 'MISSING_NOT_ZERO_NOT_GENOME_OR_DEPTH')
   pairs=evidence['allele_map'];reported=str(row.get('clinvar_allele_id','')).strip()
   if reported in ['', '.', '-']:status='EXACT_ALLELES_VALIDATED_SOURCE_ID_MAP' if pairs else 'NO_EXACT_SOURCE_ALLELES_NOT_ID_EQUIVALENCE'
   else:
    parsed=re.split('[&,;|]',reported)
    if not all(v.isdigit() for v in parsed):status='UNPARSEABLE_ALLELEID_NO_JOIN';pairs=[]
    else:
     pairs=[p for p in pairs if p['AlleleID'] in parsed]
     if {p['AlleleID'] for p in pairs}!=set(parsed):status='ALLELEID_COORDINATE_CONFLICT_OR_SNAPSHOT_DIFFERENCE_NO_JOIN';pairs=[]
     else:status='EXACT_ALLELES_AND_EXPLICIT_ALLELEID_TO_VARIATIONID'
   result.update(allref_SCV_map_status=status,allref_ClinVar_AlleleIDs='|'.join(sorted({p['AlleleID'] for p in pairs})),allref_ClinVar_VariationIDs='|'.join(sorted({p['VariationID'] for p in pairs})))
   allowed={p['VariationID'] for p in pairs};scvs=[r for r in evidence['SCV'] if r['VariationID'] in allowed]
   if scvs:
    classes=sorted({r['source_classification'] for r in scvs});result.update(allref_SCV_reference_count=str(len(scvs)),allref_SCV_classifications_json=json.dumps(classes),allref_SCV_version_examples_json=json.dumps(sorted({r['SCV_accession_version'] for r in scvs})[:5]))
    if len(classes)>1:flags.append('SCV_CLASSIFICATIONS_DIFFER_CONDITION_AND_DATE_REVIEW');details['SCV_scope']='Distinct source assertion texts; not automatic classification or independent evidence'
   if 'CONFLICT' in status:flags.append('ALLELEID_SOURCE_MAP_REVIEW')
   coverage=self.coverage.get(gene)
   if coverage and (not row.get('gene') or row['gene']==self.tissue['genes'][gene]['gene']):result['allref_gene_coverage_summary_json']=json.dumps(coverage,sort_keys=True,separators=(',',':'))
  result.update(allref_review_flags=';'.join(sorted(set(flags))),allref_review_flag_details_json=json.dumps(details,sort_keys=True,separators=(',',':')))
  if self.splice_manifest is not None:result.update(splice.enrich_row(row,self.splices.get(splice.key(row['chrom'],row['pos'],row['ref'],row['alt']),[]),self.splice_source))
  return {field:str(value) for field,value in result.items()}
 def header_bundles(self):return {name:unpack(self.db.execute('SELECT payload FROM metadata WHERE key=?',(key_,)).fetchone()[0]) for name,key_ in [('public_reference','public_bundle'),('gene_context','gene_bundle')]}
 def close(self):self.db.close()
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ready-manifest',type=Path,required=True);p.add_argument('--population-rule',type=Path,required=True);p.add_argument('--tissue-index',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();print(json.dumps(build(a.ready_manifest,a.population_rule,a.tissue_index,a.output),indent=2))
