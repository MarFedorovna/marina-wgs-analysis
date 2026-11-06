#!/usr/bin/env python3
"""Local, lossless applied-evidence VCF overlay; never reruns an annotator."""
import argparse
import base64
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import uuid
import zlib
from urllib.parse import quote

CODEC='json-utf8+zlib+base64url-nopad'
INFO_IDS=['WGS_EVIDENCE','WGS_GENE','WGS_CONSEQUENCE','WGS_CLINVAR','WGS_CLINVAR_REVIEW','WGS_GNOMAD_EXOMES','WGS_GNOMAD_GENOMES','WGS_CACHE_AF','WGS_TX_SCORES','WGS_LOCAL_SCORES','WGS_TISSUE_STATUS','WGS_ANN_SCOPE','WGS_ANN_STATUS']
SCALAR_TYPES={'WGS_GNOMAD_'+dataset+'_'+field:('Float' if field=='AF' else 'Integer') for dataset in ['EXOMES','GENOMES'] for field in ['AC','AN','AF','nhomalt']}
INFO_IDS.extend(SCALAR_TYPES)
NAMESPACED_INFO_IDS=INFO_IDS
MAX_PACKET_BYTES=64*1024**2
csv.field_size_limit(MAX_PACKET_BYTES)
DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
ACTIVE_DEADLINE=None
def now():return dt.datetime.now(dt.timezone.utc)
def safe(path):
 path=Path(path)
 if any(p.is_symlink() for p in [path,*path.parents]):raise ValueError('Symlink path refused')
 return path
def sha(path):
 h=hashlib.sha256()
 with safe(path).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):
   if ACTIVE_DEADLINE is not None and now()>=ACTIVE_DEADLINE:raise TimeoutError('Original deadline during checksum')
   h.update(b)
 return h.hexdigest()
def pin(path):
 path=safe(Path(path).absolute());return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
def verify(p):
 path=safe(p['path'])
 if path.stat().st_size!=p['bytes'] or sha(path)!=p['sha256']:raise ValueError('Pinned input changed: '+str(path))
 return path
def load(path):return json.loads(safe(path).read_text(encoding='utf-8'))
def dump(path,value):
 with safe(path).open('x',encoding='utf-8') as f:json.dump(value,f,ensure_ascii=False,indent=2);f.write('\n')
def canonical_tsv(values):
 f=io.StringIO(newline='');csv.writer(f,delimiter='\t',lineterminator='\n').writerow(values);return f.getvalue().encode('utf-8')
def encode_packet(value):
 raw=json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode('utf-8')
 if len(raw)>MAX_PACKET_BYTES:raise ValueError('One allele evidence exceeds bounded 64MiB packet')
 return base64.urlsafe_b64encode(zlib.compress(raw,6)).decode('ascii').rstrip('=')
def decode_packet(token):
 if not isinstance(token,str) or len(token)>2*MAX_PACKET_BYTES or not re.fullmatch('[A-Za-z0-9_-]+',token):raise ValueError('Invalid/oversize evidence encoding')
 encoded=base64.b64decode(token+'='*(-len(token)%4),altchars=b'-_',validate=True);decoder=zlib.decompressobj();raw=decoder.decompress(encoded,MAX_PACKET_BYTES+1)
 if len(raw)>MAX_PACKET_BYTES or decoder.unconsumed_tail or not decoder.eof or decoder.unused_data:raise ValueError('Oversize/truncated/concatenated evidence packet')
 value=json.loads(raw)
 if not isinstance(value,dict) or value.get('v')!=1:raise ValueError('Unsupported evidence schema')
 return value
encode_schema=encode_packet
decode_schema=decode_packet
def select_sample(packet,fields,sample):
 if 'sample' not in fields:raise ValueError('Sample field required')
 result=dict(packet);result['rows']=[r for r in packet['rows'] if r['values'][fields.index('sample')]==sample];result['scope']='sample';result['sample']=sample
 # Applied allele evidence and HGNC context are public and contain no cohort genotypes.
 return result
def esc(value):return quote(str(value),safe='._:-')
def rebuild_readable(packet,fields):
 rows=[dict(zip(fields,r['values'])) for r in packet['rows']]
 def unique(values):return ','.join(sorted({esc(v) for v in values if v not in ['',None,'.']})) or '.'
 result={'WGS_GENE':unique(r.get('gene','') for r in rows),'WGS_CONSEQUENCE':unique(r.get('consequence','') for r in rows),'WGS_CLINVAR':unique('|'.join([r.get('clinvar_significance',''),r.get('clinvar_ids',''),r.get('clinvar_snapshot','')]) for r in rows),'WGS_CLINVAR_REVIEW':unique('|'.join([r.get('clinvar_review',''),r.get('allref_SCV_map_status',''),r.get('allref_review_flags','')]) for r in rows)}
 for cohort in ['exomes','genomes']:
  values=[]
  for r in packet.get('evidence',{}).get('population',[]):
   if r['dataset']==cohort:values.append('|'.join(k+':'+r.get(k,'') for k in ['source_version','match_status','AC','AN','AF','nhomalt','population_FILTER']))
  result['WGS_GNOMAD_'+cohort.upper()]=unique(values)
  exact=[r for r in packet.get('evidence',{}).get('population',[]) if r.get('dataset')==cohort and r.get('match_status')=='EXACT_ALLELE' and '4.1.1' in r.get('source_version','')]
  for field in ['AC','AN','AF','nhomalt']:
   candidates={str(r.get(field,'')) for r in exact};value=next(iter(candidates)) if len(candidates)==1 else '.'
   if field=='AF':
    try:valid=math.isfinite(float(value)) and 0<=float(value)<=1
    except ValueError:valid=False
   else:valid=bool(re.fullmatch(r'\d+',value)) and int(value)<=2147483647
   result['WGS_GNOMAD_'+cohort.upper()+'_'+field]=value if valid else '.'
 # Explicit exact transcript/sample labels retain score families and missing values.
 score_fields=['SIFT','PolyPhen','FATHMM','CADD_PHRED','AlphaMissense_score','AlphaMissense_match_status','fourth_predictor_identity']
 result['WGS_TX_SCORES']=unique('|'.join([r.get('sample',''),r.get('transcript',''),r.get('protein_id',''),*[k+':'+r.get(k,'') for k in score_fields]]) for r in rows)
 result['WGS_LOCAL_SCORES']=unique('|'.join([r.get('sample',''),r.get('transcript',''),'source_original_local_evidence_json:'+r.get('local_evidence_json','')]) for r in rows)
 result['WGS_CACHE_AF']=unique('|'.join(['af_max_observed:'+r.get('af_max_observed',''),'af_status:'+r.get('af_status',''),'gnomad_version:'+r.get('gnomad_version',''),'source_af_fields_json:'+r.get('af_fields_json',''),'scope:original_primary_cache_AF_not_whole_genome_AC_AN']) for r in rows)
 result['WGS_TISSUE_STATUS']=unique('|'.join([r.get('gene',''),r.get('transcript',''),r.get('allref_tissue_gtex_v11_transcript_match_status',''),r.get('allref_tissue_pext_source_transcript_match_status',''),'GTEx_v11_pext_GTEx_v10_UniProt_gene_context']) for r in rows)
 result['WGS_ANN_SCOPE']=esc('GRCh38_exact_allele;public39_CDS_subset;exome_AN_position_only;HGNC_coverage_gene_bins_not_depth;MT_separate_appendix;no_automatic_ACMG')
 result['WGS_ANN_STATUS']=unique('|'.join([r.get('allref_tissue_status',''),r.get('allref_population_status','')]) for r in rows) if rows else 'NO_MASTER_ANNOTATION_ROWS'
 return result
def key(chrom,pos,ref,alt):
 chrom=str(chrom).removeprefix('chr');return 'chr'+chrom+':'+str(int(pos))+':'+ref+':'+alt
def pack(value):return zlib.compress(json.dumps(value,ensure_ascii=False,separators=(',',':')).encode('utf-8'),3)
def unpack(value):return json.loads(zlib.decompress(value))
def text(path):return gzip.open(path,'rt',encoding='utf-8',newline='') if str(path).endswith('.gz') else Path(path).open(encoding='utf-8',newline='')
def guard(rule,deadline,paths):
 if now()>=deadline:raise TimeoutError('Original VM budget deadline reached')
 if shutil.disk_usage(Path(rule['state_directory'])).free<rule['min_free_bytes']:raise ValueError('Required disk reserve reached')
 for path,limit in paths:
  if Path(path).exists() and Path(path).stat().st_size>limit:raise ValueError('Bounded intermediate/output limit reached')
def progress(rule,phase,**counts):
 target=safe(Path(rule['state_directory'])/'progress.json');tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(dict(phase=phase,pid=os.getpid(),updated_utc=now().isoformat(),**counts),indent=2)+'\n');os.replace(tmp,target)
def schema_header(schema):return f'##WGS_EVIDENCE_SCHEMA=<Version=1,Codec="{CODEC}",Data="{encode_schema(schema)}">\n'.encode()
def chunk_headers(name,token):
 chunks=[token[i:i+60000] for i in range(0,len(token),60000)];digest=hashlib.sha256(token.encode()).hexdigest()
 return [f'##{name}=<Version=1,Part={i+1},Parts={len(chunks)},SHA256="{digest}",Data="{chunk}">\n'.encode() for i,chunk in enumerate(chunks)]
def decode_header_bundle(lines,name):
 chunks={};parts=None;digest=None
 for raw in lines:
  line=raw.decode() if isinstance(raw,bytes) else raw
  if not line.startswith('##'+name+'='):continue
  match=re.fullmatch(r'##'+re.escape(name)+r'=<Version=1,Part=(\d+),Parts=(\d+),SHA256="([a-f0-9]{64})",Data="([A-Za-z0-9_-]+)">\s*',line)
  if not match:raise ValueError('Invalid chunk header')
  i,n,d,chunk=match.groups();i=int(i);n=int(n)
  if n>4096 or i<1 or i>n or i in chunks or (parts is not None and (parts!=n or digest!=d)):raise ValueError('Conflicting/oversize bundle chunks')
  parts=n;digest=d;chunks[i]=chunk
 if not parts or len(chunks)!=parts:raise ValueError('Missing bundle chunks')
 token=''.join(chunks[i] for i in range(1,parts+1))
 if hashlib.sha256(token.encode()).hexdigest()!=digest:raise ValueError('Header bundle checksum differs')
 return decode_packet(token)
def public_bundle(rule):
 index_path=verify(rule['tissue_context_index']);index=load(index_path);sources=[];total=index_path.stat().st_size
 for path,item in index['verified_raw_sources'].items():
  value=verify(dict(path=path,**item));total+=value.stat().st_size
  if total>rule['max_header_inflated_bytes']:raise ValueError('Public applied source header cap reached')
  sources.append(dict(source=value.parent.name+'/'+value.name,bytes=item['bytes'],sha256=item['sha256'],content_utf8=value.read_bytes().decode('utf-8')))
 for source,item in index['source_manifests'].items():
  path=Path(item['path'])
  if sha(path)!=item['sha256']:raise ValueError('Public tissue source manifest changed')
  total+=path.stat().st_size
  if total>rule['max_header_inflated_bytes']:raise ValueError('Public header cap reached')
  sources.append(dict(source=source+'/resource-manifest.json',**{k:v for k,v in pin(path).items() if k!='path'},content_utf8=path.read_bytes().decode('utf-8')))
 return dict(v=1,scope='PUBLIC39_REFERENCE_OBJECTS_NOT_PATIENT_GENOTYPES_OR_VERIFIED_ALLELE_DOMAIN_PEXT_ASSIGNMENTS',context_index_pin=rule['tissue_context_index'],context_index=index,sources=sources)
def headers(schema,coverage_token,public_token):
 out=[]
 for field in INFO_IDS:
  description={'WGS_EVIDENCE':'Lossless ordered full master row vectors and applied exact-allele public evidence; decode WGS_EVIDENCE_SCHEMA','WGS_TX_SCORES':'Percent-encoded sample|exact transcript|protein|source-field:score entries; no isoform aggregation','WGS_ANN_SCOPE':'Research context scope; gene coverage is not locus depth'}.get(field,'Percent-encoded applied reference context; missing remains missing; see lossless WGS_EVIDENCE')
  if field in SCALAR_TYPES:description='Exact matched gnomAD v4.1.1 public39 CDS global '+field.rsplit('_',1)[-1]+'; missing/outside subset is missing; source FILTER and release in dataset String summary'
  out.append(f'##INFO=<ID={field},Number={1 if field=="WGS_EVIDENCE" or field in SCALAR_TYPES else "."},Type={SCALAR_TYPES.get(field,"String")},Description="{description}">\n'.encode())
 out.append(schema_header(schema))
 out.extend(chunk_headers('WGS_GENE_CONTEXT',coverage_token));out.extend(chunk_headers('WGS_PUBLIC_REFERENCE',public_token))
 return out
def strip_overlay(line):
 ending=b'\r\n' if line.endswith(b'\r\n') else b'\n' if line.endswith(b'\n') else b'';parts=line[:-len(ending)].split(b'\t') if ending else line.split(b'\t')
 items=[v for v in parts[7].split(b';') if v.split(b'=',1)[0].decode() not in INFO_IDS];parts[7]=b';'.join(items) or b'.';return b'\t'.join(parts)+ending
def index_inputs(db,manifest,deadline,rule):
 fields=manifest['original_fields']+manifest['appended_fields'];original=hashlib.sha256();full=hashlib.sha256();original.update(canonical_tsv(manifest['original_fields']));full.update(canonical_tsv(fields));count=0
 with text(verify(manifest['enriched_master'])) as f,db:
  reader=csv.reader(f,delimiter='\t')
  if next(reader)!=fields:raise ValueError('Ordered combined master schema changed')
  for count,values in enumerate(reader,1):
   if count%10000==1:guard(rule,deadline,[(rule['_db_path'],rule['max_sqlite_bytes'])]);progress(rule,'INDEXING_MASTER',rows_indexed=count-1,rows_total=manifest['row_count'])
   if len(values)!=len(fields):raise ValueError('Truncated master row')
   row=dict(zip(fields,values));allele=key(row['chrom'],row['pos'],row['ref'],row['alt']);db.execute('INSERT INTO rows VALUES(?,?,?,0)',(count,allele,pack(values)));original.update(canonical_tsv(values[:len(manifest['original_fields'])]));full.update(canonical_tsv(values))
 if count!=manifest['row_count'] or original.hexdigest()!=manifest['invariant_digest']:raise ValueError('Combined master original-row invariant failed')
 coverage={'fields':[],'genes':{}};evidence_counts={};wanted={'nuclear-population-long.tsv.gz':'population','nuclear-exome-AN-v41.tsv.gz':'exome_AN','ClinVar-validated-allele-map.tsv.gz':'allele_map','ClinVar-SCV-long.tsv.gz':'SCV','ClinVar-citation-links.tsv.gz':'citations','gene-coverage-context-long.tsv.gz':'coverage'}
 ordered=sorted(manifest['auxiliary_artifacts'],key=lambda a:list(wanted).index(Path(a['path']).name) if Path(a['path']).name in wanted else 99)
 for artifact in ordered:
  kind=wanted.get(Path(artifact['path']).name)
  if not kind:continue
  path=verify(artifact);n=0
  with text(path) as f,db:
   reader=csv.DictReader(f,delimiter='\t')
   if kind=='coverage':coverage['fields']=reader.fieldnames
   for n,row in enumerate(reader,1):
    if n%10000==1:guard(rule,deadline,[(rule['_db_path'],rule['max_sqlite_bytes'])]);progress(rule,'INDEXING_APPLIED_EVIDENCE',table=kind,rows_indexed=n-1)
    if kind=='coverage':
     if row['scope']!='GENE_BIN_CONTEXT_ONLY_NO_EXACT_LOCUS_DEPTH_CLAIM':raise ValueError('Coverage scope invalid')
     coverage['genes'].setdefault(row['hgnc_id'],[]).append([row[field] for field in coverage['fields']])
    else:
     allele=row['variant_key']
     if not db.execute('SELECT 1 FROM rows WHERE allele=? LIMIT 1',(allele,)).fetchone():raise ValueError('Applied allele evidence has no combined master row')
     if kind=='population' and (allele!=key(row['chrom'],row['pos'],row['ref'],row['alt']) or row['dataset'] not in ['exomes','genomes'] or '4.1.1' not in row['source_version']):raise ValueError('Population allele/release contract')
     if kind=='allele_map' and (not row['AlleleID'].isdigit() or not row['VariationID'].isdigit() or not row['scope'].startswith('EXACT_GRCh38_ALLELES_SOURCE_ID_MAP')):raise ValueError('Unvalidated source allele mapping')
     if kind in ['SCV','citations'] and not db.execute('SELECT 1 FROM evidence WHERE kind=? AND allele=? AND variation=?',('allele_map',allele,row['VariationID'])).fetchone():raise ValueError('SCV/citation lacks exact validated allele map')
     db.execute('INSERT INTO evidence VALUES(?,?,?,?)',(kind,allele,row.get('VariationID',''),pack(row)))
  evidence_counts[kind]=n
 if set(wanted.values())!=set(evidence_counts):raise ValueError('All six actually applied nuclear evidence tables required')
 return fields,full.hexdigest(),coverage,evidence_counts
def embed(rule,deadline=None):
 global ACTIVE_DEADLINE
 import pysam
 if pysam.__version__!='0.24.1':raise ValueError('Exact reviewed pysam 0.24.1 required')
 deadline=min(deadline or DEADLINE,DEADLINE);ACTIVE_DEADLINE=deadline;output=safe(rule['output']);state=safe(rule['state_directory']);state.mkdir(parents=True,exist_ok=True)
 if output.exists():raise ValueError('Existing output cannot be replaced')
 if os.name!='nt':
  import resource
  resource.setrlimit(resource.RLIMIT_AS,(rule['memory_limit_bytes'],rule['memory_limit_bytes']))
  if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:rule['cpu_limit']]))
 combined=Path(rule['combined_root']);receipt_path=combined/'workflow-receipt.json';receipt=load(receipt_path)
 if receipt.get('status')!='PASS' or receipt.get('stage')!='combined-annotation-master-v1':raise ValueError('Combined master PASS required')
 manifest_path=combined/'manifest.json';manifest=load(manifest_path)
 if manifest.get('status')!='PASS_COMBINED_ANNOTATION_MASTER' or manifest.get('ranking_eligibility_unchanged') is not True:raise ValueError('Exact combined master invariant contract required')
 for name in ['manifest.json','master_evidence.enriched.tsv.gz']:
  expected=receipt['artifacts'][name];verify(dict(path=str(combined/name),bytes=expected['bytes'],sha256=expected['sha256']))
 source_root=Path(rule['workflow_root'])/'stages/06_annotate';primary_path=source_root/'workflow-receipt.json';primary=load(primary_path)
 if primary.get('status')!='PASS' or pin(primary_path) not in manifest['stage_receipts']:raise ValueError('Matching primary stage6 PASS required')
 base=verify(dict(path=str(source_root/'annotated.vcf.gz'),**primary['artifacts']['annotated.vcf.gz']));rule=dict(rule,max_output_bytes=min(rule['max_output_bytes'],max(1024**2,base.stat().st_size*rule['max_inflation_ratio'])))
 db_path=state/('evidence-index-'+uuid.uuid4().hex+'.sqlite');db=sqlite3.connect(db_path);db.executescript('PRAGMA cache_size=-65536;PRAGMA mmap_size=0;PRAGMA temp_store=FILE;CREATE TABLE rows(i INTEGER PRIMARY KEY,allele TEXT,payload BLOB,used INTEGER);CREATE INDEX row_allele ON rows(allele);CREATE TABLE evidence(kind TEXT,allele TEXT,variation TEXT,payload BLOB);CREATE INDEX evidence_key ON evidence(allele,kind,variation);CREATE TABLE seen(allele TEXT PRIMARY KEY);')
 stage=output.parent/('.'+output.name+'.'+uuid.uuid4().hex);stage.mkdir(mode=0o700)
 try:
  rule=dict(rule,_db_path=str(db_path));fields,full_digest,coverage,evidence_counts=index_inputs(db,manifest,deadline,rule);schema=dict(v=1,codec=CODEC,original_fields=manifest['original_fields'],appended_fields=manifest['appended_fields'],full_fields=fields,row_count=manifest['row_count'],invariant_digest=manifest['invariant_digest'],full_fields_digest=full_digest,scope='cohort',sample=None,ordinal='1-based global joined-master order',gene_context_scope='HGNC gene bins NOT exact locus depth',readable_scores='sample|exact transcript|protein|original source-field:value; URL-percent encoded entries; no transcript reassignment')
  gene_object=dict(v=1,scope='PUBLIC_HGNC_GENE_BIN_CONTEXT_NOT_LOCUS_DEPTH',fields=coverage['fields'],genes=coverage['genes']);public_object=public_bundle(rule)
  if sum(len(json.dumps(o,ensure_ascii=False).encode()) for o in [schema,gene_object,public_object])>rule['max_header_inflated_bytes']:raise ValueError('Total inflated applied header evidence exceeds64MiB')
  gene_token=encode_packet(gene_object);public_token=encode_packet(public_object);out=stage/'cohort.fully-annotated.vcf.gz';source_digest=hashlib.sha256();source_headers=hashlib.sha256();count=0;embedded=0;header_seen=False;sample_names=[]
  with gzip.open(base,'rb') as src,pysam.BGZFile(str(out),'wb') as dst,db:
   for line in src:
    if line.startswith(b'#'):
     if line.startswith(b'##INFO=<ID=WGS_') or line.startswith((b'##WGS_EVIDENCE_SCHEMA=',b'##WGS_GENE_CONTEXT=',b'##WGS_PUBLIC_REFERENCE=')):raise ValueError('Existing overlay namespace collision')
     if line.startswith(b'#CHROM'):
      sample_names=line.rstrip(b'\r\n').decode().split('\t')[9:]
      if len(sample_names)!=4 or len(set(sample_names))!=4:raise ValueError('Exactly four original cohort samples required')
      for h in headers(schema,gene_token,public_token):dst.write(h)
      header_seen=True
     source_headers.update(line);dst.write(line);continue
    if not header_seen:raise ValueError('VCF header required')
    count+=1
    if count%1000==1:guard(rule,deadline,[(db_path,rule['max_sqlite_bytes']),(out,rule['max_output_bytes'])]);progress(rule,'EMBEDDING_VCF',records_embedded=count-1,master_rows_embedded=embedded,master_rows_total=manifest['row_count'])
    ending=b'\r\n' if line.endswith(b'\r\n') else b'\n' if line.endswith(b'\n') else b'';parts=(line[:-len(ending)] if ending else line).split(b'\t')
    if len(parts)!=13 or b',' in parts[4]:raise ValueError('Four-sample biallelic source VCF contract required')
    chrom,pos,ref,alt=parts[0].decode(),int(parts[1]),parts[3].decode(),parts[4].decode();allele=key(chrom,pos,ref,alt);db.execute('INSERT INTO seen VALUES(?)',(allele,))
    rows=[dict(i=i,values=unpack(payload)) for i,payload in db.execute('SELECT i,payload FROM rows WHERE allele=? ORDER BY i',(allele,))]
    if not rows or {r['values'][fields.index('sample')] for r in rows}!=set(sample_names):raise ValueError('Every primary allele must have explicit master rows for all four samples, including empty annotations')
    if any(r['values'][fields.index('sample')] not in sample_names for r in rows):raise ValueError('Master sample absent from source VCF')
    format_fields=parts[8].decode().split(':')
    for r in rows:
     row=dict(zip(fields,r['values']));sample_values=parts[9+sample_names.index(row['sample'])].decode().split(':');sample_format=dict(zip(format_fields,sample_values))
     for field in ['GT','FT']:
      if field in fields and field in sample_format and row[field]!=sample_format[field]:raise ValueError('Master GT/FT differs from source VCF')
    evidence={kind:[unpack(r[0]) for r in db.execute('SELECT payload FROM evidence WHERE allele=? AND kind=? ORDER BY rowid',(allele,kind))] for kind in ['population','exome_AN','allele_map','SCV','citations']}
    genes=sorted({str(r['values'][fields.index('hgnc_id')]) for r in rows if 'hgnc_id' in fields});genes=['HGNC:'+g if g.isdigit() else g for g in genes];packet=dict(v=1,allele=[chrom,pos,ref,alt],rows=rows,evidence=evidence,gene_context_keys=[g for g in genes if g in coverage['genes']],scope='cohort',sample=None)
    tags={'WGS_EVIDENCE':encode_packet(packet),**rebuild_readable(packet,fields)}
    if any(p.split(b'=',1)[0].decode() in INFO_IDS for p in parts[7].split(b';')):raise ValueError('Existing INFO namespace collision')
    parts[7]=(parts[7]+b';' if parts[7]!=b'.' else b'')+';'.join(k+'='+v for k,v in tags.items()).encode();new=b'\t'.join(parts)+ending
    if strip_overlay(new)!=line:raise ValueError('Original VCF field/INFO byte invariant violated')
    source_digest.update(line);dst.write(new);db.execute('UPDATE rows SET used=used+1 WHERE allele=?',(allele,));embedded+=len(rows)
  if not header_seen or db.execute('SELECT count(*) FROM rows WHERE used!=1').fetchone()[0] or embedded!=manifest['row_count']:raise ValueError('Every master row must be embedded exactly once')
  guard(rule,deadline,[(out,rule['max_output_bytes'])]);pysam.tabix_index(str(out),preset='vcf',force=False);guard(rule,deadline,[(out,rule['max_output_bytes'])])
  # Independent decode of final BGZF proves every packet and ordered full-vector digest.
  db.execute('CREATE TABLE decoded(i INTEGER PRIMARY KEY,payload BLOB)');decoded_count=0;restored=hashlib.sha256();restored_headers=hashlib.sha256();seek_loci={};new_headers=[];verified_records=0
  with gzip.open(out,'rb') as f,db:
   for line in f:
    if line.startswith(b'#'):
     if line.startswith((b'##WGS_PUBLIC_REFERENCE=',b'##WGS_GENE_CONTEXT=')):new_headers.append(line)
     elif not line.startswith((b'##INFO=<ID=WGS_',b'##WGS_EVIDENCE_SCHEMA=')):restored_headers.update(line)
     continue
    parts=line.rstrip(b'\r\n').split(b'\t');tags=dict(p.split(b'=',1) for p in parts[7].split(b';') if b'=' in p);packet=decode_packet(tags[b'WGS_EVIDENCE'].decode())
    if packet['allele']!=[parts[0].decode(),int(parts[1]),parts[3].decode(),parts[4].decode()]:raise ValueError('Embedded allele differs from VCF')
    verified_records+=1
    if verified_records%1000==1:guard(rule,deadline,[(db_path,rule['max_sqlite_bytes']),(out,rule['max_output_bytes'])]);progress(rule,'VERIFYING_DECODE',rows_verified=decoded_count,rows_total=manifest['row_count'])
    for r in packet['rows']:
     if len(r['values'])!=len(fields):raise ValueError('Decoded row schema mismatch')
     db.execute('INSERT INTO decoded VALUES(?,?)',(r['i'],pack(r['values'])));decoded_count+=1
    restored.update(strip_overlay(line))
    seek_loci.setdefault(parts[0].decode(),int(parts[1]))
  check=hashlib.sha256();check.update(canonical_tsv(fields))
  for expected,(i,payload) in enumerate(db.execute('SELECT i,payload FROM decoded ORDER BY i'),1):
   if i!=expected:raise ValueError('Decoded row ordinal gap')
   if expected%10000==1:guard(rule,deadline,[(db_path,rule['max_sqlite_bytes'])]);progress(rule,'VERIFYING_MASTER_ORDER',rows_verified=expected-1,rows_total=manifest['row_count'])
   values=unpack(payload);source=db.execute('SELECT payload FROM rows WHERE i=?',(i,)).fetchone()
   if source is None or unpack(source[0])!=values:raise ValueError('Decoded field vector differs from joined master')
   check.update(canonical_tsv(values))
  if decoded_count!=manifest['row_count'] or check.hexdigest()!=full_digest or restored.hexdigest()!=source_digest.hexdigest() or restored_headers.hexdigest()!=source_headers.hexdigest():raise ValueError('Lossless VCF decode/original field/header proof failed')
  if decode_header_bundle(new_headers,'WGS_GENE_CONTEXT')!=gene_object or decode_header_bundle(new_headers,'WGS_PUBLIC_REFERENCE')!=public_object:raise ValueError('Self-contained public header roundtrip failed')
  for source in public_object['sources']:
   raw=source['content_utf8'].encode('utf-8')
   if len(raw)!=source['bytes'] or hashlib.sha256(raw).hexdigest()!=source['sha256']:raise ValueError('Embedded public source verbatim checksum failed')
  seeks=0
  with pysam.TabixFile(str(out)) as tabix:
   for chrom,pos in seek_loci.items():
    guard(rule,deadline,[(out,rule['max_output_bytes'])])
    if not any(int(s.split('\t')[1])==pos for s in tabix.fetch(chrom,pos-1,pos)):raise ValueError('Indexed random seek failed')
    seeks+=1
  def final_pin(path):return dict(path=str(output/path.name),bytes=path.stat().st_size,sha256=sha(path))
  result=dict(schema_version=1,status='PASS_FULLY_ANNOTATED_VCF',base_vcf=pin(base),master_manifest=pin(manifest_path),master=pin(verify(manifest['enriched_master'])),combined_receipt=pin(receipt_path),enriched_vcf=final_pin(out),index=final_pin(Path(str(out)+'.tbi')),schema=schema,original_fields=manifest['original_fields'],appended_fields=manifest['appended_fields'],full_fields=fields,row_count=embedded,record_count=count,invariant_digest=manifest['invariant_digest'],full_fields_digest=full_digest,original_records_fields_preserved=True,original_headers_preserved=True,original_headers_digest=source_headers.hexdigest(),original_records_digest=source_digest.hexdigest(),decoded_original_records_digest=restored.hexdigest(),packet_decode_verified=True,decoded_master_vectors_verified=True,public_header_bundle_verified=True,indexed_seek_checks=seeks,applied_evidence_counts=evidence_counts,gene_context_sha256=hashlib.sha256(gene_token.encode()).hexdigest(),public_reference_sha256=hashlib.sha256(public_token.encode()).hexdigest(),MT_scope='SEPARATE_MT_APPENDIX_NO_NUCLEAR_ROWS',MT_appendix_artifacts=[a for a in manifest['auxiliary_artifacts'] if a['role']=='MT_APPENDIX'],stage_receipts=manifest['stage_receipts'],upstream_receipts=[pin(receipt_path)],producer_code_sha256=rule['code_sha256'],network_requests=0,primary_outputs_changed=False,created_utc=now().isoformat())
  dump(stage/'manifest.json',result);artifacts={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in stage.iterdir() if p.is_file()};identity=hashlib.sha256(json.dumps(dict(base=result['base_vcf'],master=result['master'],code=rule['code_sha256']),sort_keys=True).encode()).hexdigest();dump(stage/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='fully-annotated-vcf-v1',identity=identity,artifacts=artifacts,metadata=result))
  if output.exists():raise ValueError('Concurrent output exists')
  # Same-filesystem exclusive promotion, shared helper has Linux RENAME_NOREPLACE.
  import ctypes
  if os.name=='nt':os.rename(stage,output)
  else:
   libc=ctypes.CDLL(None,use_errno=True)
   if libc.renameat2(-100,os.fsencode(stage),-100,os.fsencode(output),1):raise OSError(ctypes.get_errno(),'Exclusive promotion refused')
  return result
 finally:db.close();ACTIVE_DEADLINE=None
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['decode']);p.add_argument('--vcf',type=Path,required=True);p.add_argument('--master',type=Path,required=True);a=p.parse_args()
 # Standalone portable decoder: no external master/database necessary.
 import tempfile
 temporary=tempfile.TemporaryDirectory(prefix='portable-vcf-decode-');db=sqlite3.connect(str(Path(temporary.name)/'rows.sqlite'));db.executescript('PRAGMA cache_size=-65536;PRAGMA mmap_size=0;PRAGMA temp_store=FILE;CREATE TABLE rows(i INTEGER PRIMARY KEY,payload BLOB);');schema=None
 with gzip.open(a.vcf,'rt') as f:
  for line in f:
   if line.startswith('##WGS_EVIDENCE_SCHEMA='):schema=decode_schema(re.search('Data="([A-Za-z0-9_-]+)"',line).group(1))
   elif not line.startswith('#'):
    token=dict(p.split('=',1) for p in line.split('\t')[7].split(';') if '=' in p)['WGS_EVIDENCE']
    for r in decode_packet(token)['rows']:db.execute('INSERT INTO rows VALUES(?,?)',(r['i'],pack(r['values'])))
 if schema is None:raise ValueError('Embedded schema required')
 if a.master.exists():raise ValueError('Decoder never replaces an existing master')
 destination=safe(a.master);staging=safe(destination.with_name('.'+destination.name+'.'+uuid.uuid4().hex));full=hashlib.sha256();original=hashlib.sha256();full.update(canonical_tsv(schema['full_fields']));original.update(canonical_tsv(schema['original_fields']));count=0
 with staging.open('x',encoding='utf-8',newline='') as f:
  writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(schema['full_fields'])
  for count,(i,payload) in enumerate(db.execute('SELECT i,payload FROM rows ORDER BY i'),1):
   values=unpack(payload)
   if len(values)!=len(schema['full_fields']) or any(not isinstance(v,str) for v in values) or (schema['scope']=='cohort' and i!=count):raise ValueError('Portable decoded schema/order differs')
   writer.writerow(values);full.update(canonical_tsv(values));original.update(canonical_tsv(values[:len(schema['original_fields'])]))
 if count!=schema['row_count'] or full.hexdigest()!=schema['full_fields_digest'] or original.hexdigest()!=schema['invariant_digest']:raise ValueError('Portable decoded master count/digest failed')
 os.link(staging,destination);staging.unlink()
 db.close();temporary.cleanup()
