#!/usr/bin/env python3
"""Local official SpliceAI predictions on the fixed, approved public39 gene scope."""
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import time
from types import SimpleNamespace
import uuid

DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
SPLICE_FIELDS=['spliceai_status','spliceai_gene','spliceai_DS_AG','spliceai_DS_AL','spliceai_DS_DG','spliceai_DS_DL','spliceai_DP_AG','spliceai_DP_AL','spliceai_DP_DG','spliceai_DP_DL','spliceai_prediction_json','spliceai_source_json','spliceai_transcript_scope']
SCORE_FIELDS=['DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG','DP_DL']
SCOPE='GENCODE24_GENE_UNION_NOT_ENST49_MAPPING'
PRED_FIELDS=['variant_key','chrom','pos','ref','alt','gene','status',*SCORE_FIELDS,'raw_annotation','source_json','transcript_scope']
def now():return dt.datetime.now(dt.timezone.utc)
def safe(path):
 path=Path(path)
 if any(p.is_symlink() for p in [path,*path.parents]):raise ValueError('Symlink source/state refused')
 return path
def sha(path):
 h=hashlib.sha256()
 with safe(path).open('rb') as f:
  for block in iter(lambda:f.read(1024**2),b''):
   if now()>=DEADLINE:raise TimeoutError('Original deadline during checksum')
   h.update(block)
 return h.hexdigest()
def pin(path):
 path=safe(Path(path).absolute());return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
def verify(item):
 if pin(item['path'])!=item:raise ValueError('Pinned source changed: '+item['path'])
 return Path(item['path'])
def load(path):return json.loads(safe(path).read_text())
def dump(path,value):
 with safe(path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')
def atomic(path,value):
 path=safe(path);temporary=safe(path.with_suffix('.tmp'));temporary.write_text(json.dumps(value,indent=2)+'\n');os.replace(temporary,path)
def variant_key(chrom,pos,ref,alt):return 'chr'+chrom.removeprefix('chr')+':'+str(int(pos))+':'+ref+':'+alt
def parse_result(raw,record,distance):
 p=raw.split('|')
 if len(p)!=10 or p[0]!=record.alts[0] or not p[1]:raise ValueError('Official SpliceAI result identity/schema mismatch')
 values=dict(zip(SCORE_FIELDS,p[2:]));available=all(v!='.' for v in values.values())
 if available:
  if any(not math.isfinite(float(values[f])) or not 0<=float(values[f])<=1 for f in SCORE_FIELDS[:4]) or any(not re.fullmatch('-?[0-9]+',values[f]) or abs(int(values[f]))>distance for f in SCORE_FIELDS[4:]):raise ValueError('Official scores/offsets invalid')
 elif any(v!='.' for v in values.values()):raise ValueError('Partly missing score vector refused')
 return p[1],values,available
def classify(record,distance):
 alt=record.alts[0] if record.alts and len(record.alts)==1 else ''
 if not alt or set(record.ref+alt)-set('ACGT'):return 'UNSUPPORTED_NON_ACGT_OR_SYMBOLIC_ALLELE'
 if len(record.ref)>1 and len(alt)>1:return 'UNSUPPORTED_COMPLEX_ALLELE'
 if len(record.ref)>2*distance or len(alt)>2*distance:return 'UNSUPPORTED_ALLELE_LENGTH'
 return None
class DirectModel:
 def __init__(self,model):self.model=model
 def predict(self,x):return self.model(x,training=False).numpy()
def subset_annotation(ann,symbols):
 indices=[i for i,g in enumerate(ann.genes) if g in symbols]
 for field in ['genes','chroms','strands','tx_starts','tx_ends']:setattr(ann,field,getattr(ann,field)[indices])
 for field in ['exon_starts','exon_ends']:setattr(ann,field,[getattr(ann,field)[i] for i in indices])
 return indices
def guard(rule,deadline,path=None):
 if now()>=deadline:raise TimeoutError('Original VM budget deadline reached')
 if shutil.disk_usage(rule['state_directory']).free<rule['min_free_bytes']:raise ValueError('Required disk reserve reached')
 if path and Path(path).exists() and Path(path).stat().st_size>rule['max_sqlite_bytes']:raise ValueError('SpliceAI checkpoint cap reached')
def run(rule,deadline):
 os.environ.update(CUDA_VISIBLE_DEVICES='-1',TF_NUM_INTEROP_THREADS='1',TF_NUM_INTRAOP_THREADS=str(rule['cpu_limit']),OMP_NUM_THREADS=str(rule['cpu_limit']),TF_CPP_MIN_LOG_LEVEL='2')
 import numpy as np
 import pysam
 import tensorflow as tf
 import importlib.metadata
 from spliceai.utils import Annotator,get_delta_scores
 if tf.__version__!='2.15.1' or np.__version__!='1.26.4' or importlib.metadata.version('spliceai')!='1.3.1':raise ValueError('Exact reviewed SpliceAI/TF/numpy runtime required')
 for item in rule['source_pins']:verify(item)
 state=safe(rule['state_directory']);output=safe(rule['output']);state.mkdir(parents=True,exist_ok=True)
 if output.exists():raise ValueError('Output already exists')
 if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:rule['cpu_limit']]))
 panel=load(verify(rule['panel']));symbols={g['symbol'] for g in panel['genes']};annotation=verify(rule['annotation']);reference=rule['reference'];refpath=safe(reference['path'])
 if refpath.stat().st_size!=reference['bytes'] or refpath.stat().st_mtime_ns!=reference['mtime_ns']:raise ValueError('Accepted FASTA metadata changed')
 control=load(verify(rule['public_control']))
 if control.get('status')!='PASS_PUBLIC_SPLICEAI_CPU_DISPATCH_EQUIVALENCE' or control.get('patient_inputs') is not False or control['distance']!=rule['distance'] or control['mask']!=rule['mask'] or control['max_probability_absolute_error']>1e-7:raise ValueError('Matching actual public equivalence control required')
 ann=Annotator(str(refpath),str(annotation));subset_annotation(ann,symbols);ann.models=[DirectModel(model) for model in ann.models]
 source_json=json.dumps(dict(version='1.3.1',model_ensemble='five official SpliceAI h5 models',model_sha256=[p['sha256'] for p in rule['models']],annotation_sha256=rule['annotation']['sha256'],annotation_scope=SCOPE,reference_sha256=reference['sha256'],assembly='GRCh38',distance=rule['distance'],mask=rule['mask'],tensorflow='2.15.1',dispatch='official get_delta_scores, public-control-equivalent direct Keras inference',no_transcript_reassignment=True),sort_keys=True,separators=(',',':'))
 path=state/'checkpoint.sqlite';db=sqlite3.connect(path);db.executescript('PRAGMA cache_size=-65536;PRAGMA mmap_size=0;CREATE TABLE IF NOT EXISTS metadata(identity TEXT);CREATE TABLE IF NOT EXISTS variants(key TEXT PRIMARY KEY);CREATE TABLE IF NOT EXISTS predictions(key TEXT,gene TEXT,payload TEXT,PRIMARY KEY(key,gene));')
 identity=hashlib.sha256(json.dumps(dict(source=rule['source_vcf'],models=rule['models'],annotation=rule['annotation'],panel=rule['panel'],distance=rule['distance'],mask=rule['mask'],code=rule['code_sha256']),sort_keys=True).encode()).hexdigest();old=db.execute('SELECT identity FROM metadata').fetchone()
 if old and old[0]!=identity:raise ValueError('Resume identity changed')
 if not old:db.execute('INSERT INTO metadata VALUES(?)',(identity,));db.commit()
 completed=db.execute('SELECT count(*) FROM variants').fetchone()[0];start=time.monotonic();last_progress=0;seen=set();inference=0;statuses={}
 try:
  with pysam.VariantFile(str(verify(rule['source_vcf']))) as vcf:
   for chrom,start_pos,end_pos in sorted({('chr'+str(c).removeprefix('chr'),int(s)-1,int(e)) for c,s,e in zip(ann.chroms,ann.tx_starts,ann.tx_ends)}):
    if chrom not in vcf.header.contigs:continue
    for raw in vcf.fetch(chrom,start_pos,end_pos):
     guard(rule,deadline,path);record=SimpleNamespace(chrom=raw.chrom,pos=raw.pos,ref=raw.ref,alts=raw.alts);alt=record.alts[0] if record.alts and len(record.alts)==1 else '';key=variant_key(record.chrom,record.pos,record.ref,alt)
     if key in seen:continue
     seen.add(key)
     if db.execute('SELECT 1 FROM variants WHERE key=?',(key,)).fetchone():continue
     genes,strands,idxs=ann.get_name_and_strand(record.chrom,record.pos);error=classify(record,rule['distance']);records=[]
     if not error:
      if ann.ref_fasta[chrom][record.pos-1:record.pos-1+len(record.ref)].seq.upper()!=record.ref:error='REFERENCE_MISMATCH_NO_SCORE'
      elif record.pos<=5000+rule['distance'] or record.pos+5000+rule['distance']>len(ann.ref_fasta[chrom]):error='REFERENCE_WINDOW_UNAVAILABLE_NO_SCORE'
     if not error:
      values=get_delta_scores(record,ann,rule['distance'],rule['mask']);inference+=1
      if not values:error='NO_MODEL_RESULT_NO_SCORE'
      for result in values:
       gene,scores,available=parse_result(result,record,rule['distance'])
       if gene not in genes:raise ValueError('Unexpected model gene')
       records.append(dict(variant_key=key,chrom=record.chrom,pos=record.pos,ref=record.ref,alt=alt,gene=gene,status='SCORED_EXACT_GRCH38_ALLELE_GENE' if available else 'UNSUPPORTED_MODEL_RESULT_MISSING',**{f:scores[f] if available else '' for f in SCORE_FIELDS},raw_annotation=result,source_json=source_json,transcript_scope=SCOPE))
     if error:
      for gene in genes:records.append(dict(variant_key=key,chrom=record.chrom,pos=record.pos,ref=record.ref,alt=alt,gene=str(gene),status=error,**{f:'' for f in SCORE_FIELDS},raw_annotation='',source_json=source_json,transcript_scope=SCOPE))
     if {r['gene'] for r in records}!=set(genes):raise ValueError('Each overlapping approved model gene requires score or explicit failure status')
     with db:
      for r in records:db.execute('INSERT INTO predictions VALUES(?,?,?)',(key,r['gene'],json.dumps(r,separators=(',',':'))));statuses[r['status']]=statuses.get(r['status'],0)+1
      db.execute('INSERT INTO variants VALUES(?)',(key,))
     completed+=1
     if time.monotonic()-last_progress>=15:
      atomic(state/'progress.json',dict(phase='SCORING_LOCAL_SPLICEAI',pid=os.getpid(),updated_utc=now().isoformat(),variants_completed=completed,variants_estimated=rule['expected_variants'],inferences_this_run=inference,elapsed_seconds=time.monotonic()-start,seconds_per_variant=(time.monotonic()-start)/max(1,inference),statuses_this_run=statuses));last_progress=time.monotonic()
  guard(rule,deadline,path)
  if len(seen)!=rule['expected_variants'] or completed!=len(seen):raise ValueError('Actual approved-scope variant count differs from reviewed audit')
  stage=output.parent/('.'+output.name+'.'+uuid.uuid4().hex);stage.mkdir(mode=0o700);target=stage/'predictions.tsv.gz';counts={};rows=0
  with gzip.open(target,'wt',encoding='utf-8',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=PRED_FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
   for (payload,) in db.execute('SELECT payload FROM predictions ORDER BY key,gene'):
    if rows%1000==0:
     guard(rule,deadline,path)
     if target.stat().st_size>rule['max_output_bytes']:raise ValueError('SpliceAI output cap reached')
    row=json.loads(payload);writer.writerow(row);counts[row['status']]=counts.get(row['status'],0)+1;rows+=1
  if target.stat().st_size>rule['max_output_bytes']:raise ValueError('SpliceAI output cap reached')
  scored=counts.get('SCORED_EXACT_GRCH38_ALLELE_GENE',0)
  if not scored:raise ValueError('No actual predictions; cannot publish placeholder success')
  manifest=dict(schema_version=1,status='PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT',source_vcf=rule['source_vcf'],primary_receipt=rule['primary_receipt'],reference=reference,panel=rule['panel'],annotation=rule['annotation'],models=rule['models'],runtime=rule['runtime'],public_control=rule['public_control'],distance=rule['distance'],mask=rule['mask'],gene_scope=SCOPE,unmatched_panel_genes=['MT-TI'],outside_scope_status='OUTSIDE_APPROVED_PUBLIC39_MODEL_GENE_SPANS_NOT_ZERO',variant_count=completed,prediction_rows=rows,status_counts=counts,predictions=dict(path=str(output/target.name),bytes=target.stat().st_size,sha256=sha(target)),producer_code_sha256=rule['code_sha256'],network_requests=0,no_automatic_classification=True,primary_outputs_changed=False,created_utc=now().isoformat())
  dump(stage/'manifest.json',manifest);artifacts={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in stage.iterdir() if p.is_file()};dump(stage/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='spliceai-reference-context-v1',identity=identity,artifacts=artifacts,metadata=manifest))
  import ctypes
  if os.name=='nt':os.rename(stage,output)
  else:
   libc=ctypes.CDLL(None,use_errno=True)
   if libc.renameat2(-100,os.fsencode(stage),-100,os.fsencode(output),1):raise OSError(ctypes.get_errno(),'Exclusive output promotion failed')
  return manifest
 finally:db.close()
