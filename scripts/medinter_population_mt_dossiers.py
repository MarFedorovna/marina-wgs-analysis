#!/usr/bin/env python3
"""Build local public-reference MT, AN and ClinVar dossier sidecars; no patient input."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import time
from medinter_population_mt_sources import load,save,safe,now,HARD_DEADLINE
from medinter_population_mt_common import panel,PublicFasta,IntervalLookup,info_fields

SOURCES=Path('/srv/marina/resources/medinter-population-mt-sources-v1')
ROOT=Path('/srv/marina/resources/medinter-population-mt-sidecars-v1')
def source(name):
 receipt=load(SOURCES/(name+'-receipt.json'))
 if receipt.get('status')!='DOWNLOADED_VERIFIED_PUBLIC_SOURCE_ONLY_NOT_BOUND':raise ValueError('Source not ready: '+name)
 path=SOURCES/(name+'__'+Path(receipt['source']['url']).name)
 if not path.is_file() or path.stat().st_size!=receipt['bytes']:raise ValueError('Source receipt size mismatch')
 return path,receipt
def tsv(path,header_prefix):
 csv.field_size_limit(16*1024*1024)  # Public ClinVar contains documented long descriptions/gene lists.
 opener=gzip.open if path.name.endswith(('.gz','.bgz')) else open
 with opener(path,'rt',encoding='utf-8-sig',newline='') as stream:
  for line in stream:
   if line.startswith(header_prefix):
    heading=line.lstrip('#').rstrip('\r\n').split('\t');break
  else:raise ValueError('Expected tabular header absent: '+header_prefix)
  reader=csv.DictReader(stream,fieldnames=heading,delimiter='\t')
  for index,row in enumerate(reader):
   if index%10000==0 and now()>=HARD_DEADLINE:raise TimeoutError('Original budget deadline')
   yield row
def open_db(name):
 path=safe(ROOT/name)
 if path.exists():raise ValueError('Existing sidecar retained; no overwrite')
 return sqlite3.connect(path)
def mt():
 fasta=PublicFasta();rcrs=Path('/srv/marina/resources/medinter-population-mt-metadata-v1/rcrs_NC_012920_1.verbatim')
 if not rcrs.exists():raise ValueError('Official rCRS metadata not downloaded')
 sequence=''.join(line.strip() for line in rcrs.read_text().splitlines() if not line.startswith('>')).upper()
 if len(sequence)!=16569 or fasta.sequence('chrM',1,16569)!=sequence:raise ValueError('Official NC_012920.1 and local public reference differ')
 db=open_db('public-mt-population-v2.sqlite');db.executescript('CREATE TABLE gnomad(pos INTEGER,ref TEXT,alt TEXT,filter TEXT,AN INTEGER,AC_hom INTEGER,AC_het INTEGER,AF_hom TEXT,AF_het TEXT,info TEXT,PRIMARY KEY(pos,ref,alt));CREATE TABLE helix(pos INTEGER,ref TEXT,alt TEXT,counts_hom INTEGER,AF_hom TEXT,counts_het INTEGER,AF_het TEXT,AN INTEGER,raw_json TEXT,PRIMARY KEY(pos,ref,alt));CREATE TABLE quarantine(source TEXT,pos INTEGER,ref TEXT,alt TEXT,reason TEXT,raw TEXT);')
 counts={'gnomad':0,'helix':0,'ref_quarantine':0,'multiallelic_quarantine':0};headers=[]
 gnomad,gnomad_receipt=source('gnomad_mt_v3_1')
 with gzip.open(gnomad,'rt') as stream,db:
  for line in stream:
   if line.startswith('#'):headers.append(line.rstrip());continue
   row=line.rstrip('\n').split('\t');chrom,pos,ref,alt=row[0],int(row[1]),row[3],row[4];info=info_fields(row[7])
   if chrom!='chrM' or ',' in alt:raise ValueError('Unexpected MT coordinate/allele encoding')
   if fasta.sequence('chrM',pos,len(ref))!=ref:
    db.execute('INSERT INTO quarantine VALUES (?,?,?,?,?,?)',('gnomad',pos,ref,alt,'rCRS_REF_MISMATCH',line));counts['ref_quarantine']+=1;continue
   for field in ['AN','AC_hom','AC_het','AF_hom','AF_het']:
    if field not in info:raise ValueError('Missing MT field '+field)
   an=int(info['AN']);hom=int(info['AC_hom']);het=int(info['AC_het'])
   if min(an,hom,het)<0 or hom+het>an:raise ValueError('Invalid mtDNA individual counts')
   db.execute('INSERT INTO gnomad VALUES (?,?,?,?,?,?,?,?,?,?)',(pos,ref,alt,row[6],an,hom,het,info['AF_hom'],info['AF_het'],row[7]));counts['gnomad']+=1
 helix,helix_receipt=source('helix_mtdb_20200327')
 with db:
  for row in tsv(helix,'locus\t'):
   chrom,value=row['locus'].split(':');pos=int(value);alleles=json.loads(row['alleles'])
   if len(alleles)!=2:
    db.execute('INSERT INTO quarantine VALUES (?,?,?,?,?,?)',('helix',pos,alleles[0],json.dumps(alleles[1:]),'MULTIALLELIC_COUNTS_NOT_ASSIGNED_TO_SINGLE_ALLELE',json.dumps(row)));counts['multiallelic_quarantine']+=1;continue
   ref,alt=alleles
   if chrom!='chrM':raise ValueError('Unexpected Helix contig')
   if fasta.sequence('chrM',pos,len(ref))!=ref:
    db.execute('INSERT INTO quarantine VALUES (?,?,?,?,?,?)',('helix',pos,ref,alt,'rCRS_REF_MISMATCH',json.dumps(row)));counts['ref_quarantine']+=1;continue
   hom=int(row['counts_hom']);het=int(row['counts_het'])
   if min(hom,het)<0:raise ValueError('Invalid Helix aggregate counts')
   db.execute('INSERT INTO helix VALUES (?,?,?,?,?,?,?,?,?)',(pos,ref,alt,hom,row['AF_hom'],het,row['AF_het'],None,json.dumps(row)));counts['helix']+=1
 db.close();fasta.close()
 (ROOT/'gnomad-mt-source-headers.txt').write_text('\n'.join(headers)+'\n')
 result=dict(status='READY_PUBLIC_MT_FREQUENCY_SIDECAR_NOT_BOUND',database=str(ROOT/'public-mt-population-v2.sqlite'),counts=counts,reference_accession='NC_012920.1',rcrs_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),official_rcrs_matches_local_GRCh38_chrM=True,helix_site_AN_not_supplied_kept_NULL=True,helix_healthy_control_claim=False,gnomad_heteroplasmy_and_haplogroup_fields_preserved_in_INFO=True,mitomap_status='HTTP_403_NO_BYPASS_NOT_DOWNLOADED',no_automatic_pathogenicity=True,annotation_binding_updated=False,patient_inputs=False,source_sha256={'gnomad':gnomad_receipt['sha256'],'helix':helix_receipt['sha256']},completed_utc=now().isoformat())
 save(ROOT/'mt-resource-receipt.json',result)
def an_subset():
 manifest,intervals=panel();lookup=IntervalLookup(intervals);path,receipt=source('gnomad_exomes_all_sites_an_v4_1');db=open_db('public-hcm39-exome-AN.sqlite')
 db.execute('CREATE TABLE all_sites_AN(chrom TEXT,pos INTEGER,AN INTEGER,raw_json TEXT,PRIMARY KEY(chrom,pos))');count=0;zero=0
 with db:
  for row in tsv(path,'locus\t'):
   chrom,value=row['locus'].split(':');pos=int(value)
   if not lookup.contains(chrom,pos):continue
   an=int(row['AN'])
   if an<0:raise ValueError('Negative AN')
   db.execute('INSERT INTO all_sites_AN VALUES (?,?,?,?)',(chrom,pos,an,json.dumps(row)));count+=1;zero+=an==0
 db.close()
 if not count:raise ValueError('Empty public AN subset')
 result=dict(status='READY_PUBLIC_HCM39_EXOME_ALL_SITES_AN_V4_1_NOT_BOUND',rows=count,zero_AN_rows=zero,CDS_union_bases=sum(end-start for chrom,start,end in intervals),missing_positions_not_assumed_AC_zero=True,genome_all_sites_AN_downloaded=False,depth_coverage_downloaded=False,capture_region_flags_preserved=True,panel_genes=len(manifest['genes']),source_sha256=receipt['sha256'],patient_inputs=False,annotation_binding_updated=False,completed_utc=now().isoformat())
 save(ROOT/'an-resource-receipt.json',result)
def clinvar():
 manifest,intervals=panel();symbols={x['symbol'] for x in manifest['genes']};hgnc={x['hgnc_id'] for x in manifest['genes']};db=open_db('public-hcm39-clinvar-dossier-v2.sqlite')
 db.executescript('CREATE TABLE variants(variation_id TEXT,allele_id TEXT,chrom TEXT,pos TEXT,ref TEXT,alt TEXT,raw_json TEXT,PRIMARY KEY(variation_id,allele_id,chrom,pos,ref,alt));CREATE TABLE submissions(variation_id TEXT,SCV TEXT,SCV_version_present INTEGER,submitter TEXT,date_evaluated TEXT,classification TEXT,raw_json TEXT);CREATE TABLE citations(variation_id TEXT,citation_source TEXT,citation_id TEXT,raw_json TEXT,PRIMARY KEY(variation_id,citation_source,citation_id));')
 variants,variant_receipt=source('clinvar_variant_summary_20261101');ids=set();rows=0
 with db:
  for row in tsv(variants,'#AlleleID\t'):
   genes=set(re.split('[;,|]',row['GeneSymbol']));hgncs=set(re.split('[;,|]',row['HGNC_ID']))
   if row['Assembly']!='GRCh38' or not (genes&symbols or hgncs&hgnc):continue
   variation=row['VariationID'];ids.add(variation)
   db.execute('INSERT OR IGNORE INTO variants VALUES (?,?,?,?,?,?,?)',(variation,row['AlleleID'],row['Chromosome'],row['PositionVCF'],row['ReferenceAlleleVCF'],row['AlternateAlleleVCF'],json.dumps(row)));rows+=1
 submissions,submission_receipt=source('clinvar_submission_summary_20261101');scv_count=0;unversioned=0
 with db:
  for row in tsv(submissions,'#VariationID\t'):
   if row['VariationID'] not in ids:continue
   versioned=bool(re.fullmatch(r'SCV\d+\.\d+',row['SCV']));unversioned+=not versioned;scv_count+=1
   db.execute('INSERT INTO submissions VALUES (?,?,?,?,?,?,?)',(row['VariationID'],row['SCV'],int(versioned),row['Submitter'],row['DateLastEvaluated'],row['ClinicalSignificance'],json.dumps(row)))
 citations,citation_receipt=source('clinvar_var_citations_20261031');citation_count=0
 with db:
  for row in tsv(citations,'#AlleleID\t'):
   if row['VariationID'] not in ids:continue
   before=db.total_changes;db.execute('INSERT OR IGNORE INTO citations VALUES (?,?,?,?)',(row['VariationID'],row['citation_source'],row['citation_id'],json.dumps(row)));citation_count+=db.total_changes-before
 with (ROOT/'public-hcm39-literature-review-queue.tsv').open('x',newline='') as out:
  writer=csv.writer(out,delimiter='\t');writer.writerow(['VariationID','citation_source','citation_id','reference_url','full_text_read','independent_observations_established'])
  for variation,kind,identifier in db.execute('SELECT variation_id,citation_source,citation_id FROM citations ORDER BY variation_id,citation_source,citation_id'):
   url=f'https://pubmed.ncbi.nlm.nih.gov/{identifier}/' if kind=='PubMed' and identifier.isdigit() else ''
   writer.writerow([variation,kind,identifier,url,False,False])
 db.close()
 result=dict(status='READY_PUBLIC_HCM39_CLINVAR_REFERENCE_DOSSIER_PREPARATION_NOT_CLINICAL_REVIEW',database=str(ROOT/'public-hcm39-clinvar-dossier-v2.sqlite'),public_genes=len(symbols),variation_ids=len(ids),variant_rows=rows,SCV_rows=scv_count,SCV_rows_without_version=unversioned,citation_records=citation_count,SCV_version_gap='Unversioned public submission_summary SCVs remain explicitly unversioned; exact-version XML/manual verification still required',literature_full_text_downloaded=False,literature_read=False,independent_case_counts_established=False,PP5_not_assigned=True,current_aggregate_reference_not_replacing_core_ClinVar20261025=True,source_sha256={'variants':variant_receipt['sha256'],'submissions':submission_receipt['sha256'],'citations':citation_receipt['sha256']},patient_inputs=False,annotation_binding_updated=False,completed_utc=now().isoformat())
 save(ROOT/'clinvar-dossier-resource-receipt.json',result)
def main():
 import fcntl
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--task',choices=['mt','an','clinvar'],required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
 if not a.execute:print('DRY_RUN_NO_SOURCE_READ');return
 os.umask(0o077);safe(ROOT);ROOT.mkdir(parents=True,exist_ok=True)
 with safe(ROOT/('.'+a.task+'.lock')).open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  save(ROOT/(a.task+'-status.json'),dict(status='BUILDING_PUBLIC_REFERENCE_SIDECAR',pid=os.getpid(),updated_utc=now().isoformat(),patient_inputs=False))
  try:
   {'mt':mt,'an':an_subset,'clinvar':clinvar}[a.task]()
   save(ROOT/(a.task+'-status.json'),dict(status='COMPLETE_PUBLIC_REFERENCE_SIDECAR_NOT_BOUND',pid=os.getpid(),updated_utc=now().isoformat(),patient_inputs=False))
  except Exception as exc:
   save(ROOT/(a.task+'-status.json'),dict(status='FAIL_SIDECAR_REVIEW_REQUIRED_PARTIAL_RETAINED',error=str(exc)[:240],pid=os.getpid(),updated_utc=now().isoformat()));raise
if __name__=='__main__':main()
