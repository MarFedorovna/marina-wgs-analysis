#!/usr/bin/env python3
"""Append actual allele-level SpliceAI without decoding or rewriting core vectors."""
import datetime as dt
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import uuid

DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
spec=importlib.util.spec_from_file_location('allele_info_codec',Path(__file__).resolve().parent/'spliceai_allele_info_codec.py');codec=importlib.util.module_from_spec(spec);spec.loader.exec_module(codec)
ACTIVE_DEADLINE=None
def now():return dt.datetime.now(dt.timezone.utc)
def safe(path):
 path=Path(path)
 if any(p.is_symlink() for p in [path,*path.parents]):raise ValueError('Symlink source/output refused')
 return path
def sha(path):
 h=hashlib.sha256()
 with safe(path).open('rb') as f:
  for block in iter(lambda:f.read(1024**2),b''):
   if ACTIVE_DEADLINE and now()>=ACTIVE_DEADLINE:raise TimeoutError('Original deadline during checksum')
   h.update(block)
 return h.hexdigest()
def pin(path):
 path=safe(Path(path).absolute());return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
def verify(item):
 path=safe(item['path'])
 if not path.is_file() or path.stat().st_size!=item['bytes'] or sha(path)!=item['sha256']:raise ValueError('Pinned overlay input changed: '+str(path))
 return path
def load(path):return json.loads(safe(path).read_text(encoding='utf-8'))
def dump(path,value):
 with safe(path).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2);f.write('\n')
def guard(rule,deadline,path=None):
 if now()>=deadline:raise TimeoutError('Original VM allowance reached')
 if shutil.disk_usage(rule['state_directory']).free<rule['min_free_bytes']:raise ValueError('Required remaining disk reserve reached')
 if path and Path(path).exists() and Path(path).stat().st_size>rule['max_output_bytes']:raise ValueError('Reviewed output ceiling reached')
def progress(rule,phase,**counts):
 path=Path(rule['state_directory'])/'progress.json';tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(dict(phase=phase,pid=os.getpid(),updated_utc=now().isoformat(),**counts),indent=2)+'\n');os.replace(tmp,path)
def check_receipt(root,stage,names):
 root=Path(root);path=root/'workflow-receipt.json';receipt=load(path)
 if receipt.get('status')!='PASS' or receipt.get('stage')!=stage or not re.fullmatch('[a-f0-9]{64}',receipt.get('identity','')):raise ValueError('Actual source PASS receipt required: '+stage)
 for name in names:verify(dict(path=str(root/name),**{k:receipt['artifacts'][name][k] for k in ['bytes','sha256']}))
 return receipt,pin(path)
def sources(rule,public_control=False):
 root=Path(rule['base_root']);base_receipt,base_receipt_pin=check_receipt(root,'direct-core-context-v1',['manifest.json','core-manifest.json','cohort.core-context.vcf.gz','cohort.core-context.vcf.gz.tbi']);base=load(root/'manifest.json');core=load(root/'core-manifest.json')
 if base.get('status')!='INTERMEDIATE_CORE33_PASS' or base.get('final_delivery_complete') is not False or core.get('status')!='DIRECT_PACKET_BASELINE_PASS' or any(base.get(f) is not True for f in ['original_records_fields_preserved','packet_decode_verified','direct_core_prefix_verified']):raise ValueError('Actual intermediate core33 authority required')
 if bool(base.get('synthetic_validation'))!=bool(public_control) or bool(core.get('synthetic_validation'))!=bool(public_control):raise ValueError('Public-control and production sources cannot be mixed')
 for key_,name in [('enriched_vcf','cohort.core-context.vcf.gz'),('index','cohort.core-context.vcf.gz.tbi'),('core_manifest','core-manifest.json')]:
  item=base[key_]
  if Path(item['path'])!=root/name:raise ValueError('Unexpected intermediate artifact path')
  expected=base_receipt['artifacts'][name]
  if any(item[k]!=expected[k] for k in ['bytes','sha256']):raise ValueError('Intermediate manifest and PASS receipt artifact pins differ')
 if core['schema']!=base['core_schema'] or core['schema']['full_fields']!=base['schema']['original_fields'] or base['row_count']!=base['schema']['row_count'] or len(base['schema']['appended_fields'])!=33:raise ValueError('Preserved core prefix/schema authority mismatch')
 if base['producer_code_sha256']['final_driver']!=rule['expected_core_driver_sha256'] or core['producer_code_sha256']!=base['producer_code_sha256']:raise ValueError('Completed core does not match the reviewed direct producer')
 splice_root=Path(rule['spliceai_root']);model_receipt,model_receipt_pin=check_receipt(splice_root,'spliceai-reference-context-v1',['manifest.json','predictions.tsv.gz']);model=load(splice_root/'manifest.json')
 if model.get('status')!='PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT' or model.get('distance')!=500 or model.get('mask')!=0 or model.get('gene_scope')!=codec.SCOPE or len(model['models'])!=5 or model.get('network_requests')!=0 or model.get('primary_outputs_changed') is not False:raise ValueError('Actual complete local SpliceAI model source required')
 for field,expected in rule['expected_model_provenance'].items():
  actual=model['producer_code_sha256'] if field=='code_sha256' else model[field]
  if actual!=expected:raise ValueError('Completed model does not match reviewed weights/annotation/control/producer')
 if model['source_vcf']!=rule['source_vcf'] or core['source_pins']['filter_vcf']!=rule['source_vcf']:raise ValueError('Core and model must derive from the exact accepted source05 VCF')
 verify(rule['source_vcf']);predictions={};source_json='';rows=0
 import csv
 with gzip.open(verify(model['predictions']),'rt',encoding='utf-8',newline='') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   source=codec.source(row['source_json'])
   if source['model_sha256']!=[p['sha256'] for p in model['models']] or source['annotation_sha256']!=model['annotation']['sha256'] or source['reference_sha256']!=model['reference']['sha256']:raise ValueError('Model row provenance differs from actual completion')
   predictions.setdefault(row['variant_key'],[]).append(row);source_json=row['source_json'];rows+=1
 if rows!=model['prediction_rows'] or len(predictions)!=model['variant_count'] or not any(p['status']=='SCORED_EXACT_GRCH38_ALLELE_GENE' for values in predictions.values() for p in values):raise ValueError('No complete actual score set')
 for values in predictions.values():codec.validate_predictions(values,[values[0]['chrom'],int(values[0]['pos']),values[0]['ref'],values[0]['alt']],source_json)
 projection=dict(v=1,projection_contract=codec.CONTRACT,projection_fields=codec.SPLICE_FIELDS,scope=codec.SCOPE,model=codec.source(source_json),physical_core_packets_unchanged=True,missing_status=codec.OUTSIDE)
 return dict(base=base,core=core,base_receipt=base_receipt_pin,model=model,model_receipt=model_receipt_pin,projection=projection,predictions=predictions,source_json=source_json)
def data_parts(line):
 ending=b'\r\n' if line.endswith(b'\r\n') else b'\n' if line.endswith(b'\n') else b'';parts=(line[:-len(ending)] if ending else line).split(b'\t',8)
 if len(parts)!=9:raise ValueError('Original VCF data columns required')
 return parts,ending
def info(parts):
 values={}
 for item in parts[7].split(b';'):
  if item==b'.':continue
  key_,separator,value=item.partition(b'=')
  if key_ in values:raise ValueError('Duplicate original INFO key')
  values[key_]=value if separator else None
 return values
def overlay(rule,deadline=None,public_control=False):
 global ACTIVE_DEADLINE
 import pysam
 if pysam.__version__!='0.24.1':raise ValueError('Pinned pysam0.24.1 required')
 deadline=min(deadline or DEADLINE,DEADLINE);ACTIVE_DEADLINE=deadline;state=safe(rule['state_directory']);state.mkdir(parents=True,exist_ok=True);output=safe(rule['output'])
 if output.exists():raise ValueError('Existing final output cannot be replaced')
 if os.name!='nt':
  import resource
  resource.setrlimit(resource.RLIMIT_AS,(rule['memory_limit_bytes'],rule['memory_limit_bytes']))
  if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:rule['cpu_limit']]))
 source=sources(rule,public_control);base=source['base'];stage=output.parent/('.'+output.name+'.overlay-'+uuid.uuid4().hex);stage.mkdir(mode=0o700);target=stage/'cohort.fully-annotated.vcf.gz';count=0;selected=0;seen=set();before=hashlib.sha256();before_headers=hashlib.sha256();before_tokens=hashlib.sha256();sample_names=None
 try:
  with gzip.open(base['enriched_vcf']['path'],'rb') as src,pysam.BGZFile(str(target),'wb') as dst:
   for line in src:
    if line.startswith(b'#'):
     if line.startswith(b'##WGS_SPLICEAI_SCHEMA=') or any(line.startswith(('##INFO=<ID='+field+',').encode()) for field in codec.OVERLAY_INFO_IDS):raise ValueError('Input already has a SpliceAI overlay namespace')
     if line.startswith(b'#CHROM'):
      sample_names=line.rstrip(b'\r\n').decode().split('\t')[9:]
      if sample_names!=base['samples'] or len(sample_names)!=4:raise ValueError('Original four-sample identity differs')
      for header in codec.overlay_headers(source['projection']):dst.write(header)
     before_headers.update(line);dst.write(line);continue
    count+=1
    if count%1000==1:guard(rule,deadline,target);progress(rule,'APPENDING_ALLELE_SPLICEAI_INFO',records=count-1,records_total=base['record_count'],matched_model_alleles=selected)
    parts,ending=data_parts(line);old=info(parts)
    if b'WGS_EVIDENCE' not in old or any(field.encode() in old for field in codec.OVERLAY_INFO_IDS):raise ValueError('Required unchanged core packet or unused overlay namespace missing')
    allele=[parts[0].decode(),int(parts[1]),parts[3].decode(),parts[4].decode()];identity=codec.projection.key(*allele);predictions=source['predictions'].get(identity,[])
    if identity in seen:raise ValueError('Duplicate source allele in overlay')
    if predictions:seen.add(identity);selected+=1
    tags=codec.overlay_info(allele,predictions,source['source_json'] if predictions else '');parts[7]=(parts[7]+b';' if parts[7]!=b'.' else b'')+';'.join(k+'='+v for k,v in tags.items()).encode();new=b'\t'.join(parts)+ending
    if codec.strip_overlay(new,only_spliceai=True)!=line:raise ValueError('An original field/token changed while appending scores')
    before.update(line);before_tokens.update(old[b'WGS_EVIDENCE']+b'\n');dst.write(new)
  if sample_names is None or count!=base['record_count'] or seen!=set(source['predictions']):raise ValueError('All source records and selected model alleles must be retained exactly once')
  guard(rule,deadline,target);progress(rule,'INDEXING_FINAL_VCF',records=count,matched_model_alleles=selected);pysam.tabix_index(str(target),preset='vcf',force=False);guard(rule,deadline,target)
  restored=hashlib.sha256();restored_headers=hashlib.sha256();tokens=hashlib.sha256();verified=0;model_verified=0;seek={};projection_seen=None
  with gzip.open(target,'rb') as f:
   for line in f:
    if line.startswith(b'#'):
     if line.startswith(b'##WGS_SPLICEAI_SCHEMA='):projection_seen=codec.decode_projection_schema(line)
     elif not any(line.startswith(('##INFO=<ID='+field+',').encode()) for field in codec.OVERLAY_INFO_IDS):restored_headers.update(line)
     continue
    verified+=1
    if verified%1000==1:guard(rule,deadline,target);progress(rule,'VERIFYING_ORIGINAL_RECORDS_AND_MODEL_INFO',records_verified=verified-1,records_total=count)
    parts,ending=data_parts(line);fields=info(parts);allele=[parts[0].decode(),int(parts[1]),parts[3].decode(),parts[4].decode()];predictions,metadata=codec.projection_from_info({k.decode():v.decode() if v is not None else None for k,v in fields.items() if k.decode() in codec.OVERLAY_INFO_IDS},allele)
    if predictions:
     identity=codec.projection.key(*allele)
     if predictions!=source['predictions'].get(identity):raise ValueError('Final applied model data differs from completed local inference')
     model_verified+=1
    restored.update(codec.strip_overlay(line,only_spliceai=True));tokens.update(fields[b'WGS_EVIDENCE']+b'\n');seek.setdefault(allele[0],allele[1])
  if verified!=count or model_verified!=selected or restored.hexdigest()!=before.hexdigest() or tokens.hexdigest()!=before_tokens.hexdigest() or restored_headers.hexdigest()!=before_headers.hexdigest() or projection_seen!=source['projection']:raise ValueError('Independent original-record/token/header/model proof failed')
  with pysam.TabixFile(str(target)) as tabix:
   for chrom,pos in seek.items():
    guard(rule,deadline,target)
    if not any(int(row.split('\t')[1])==pos for row in tabix.fetch(chrom,pos-1,pos)):raise ValueError('Indexed final locus check failed')
  def final_pin(path):return dict(path=str(output/path.name),bytes=path.stat().st_size,sha256=sha(path))
  applied=list(base.get('applied_artifacts',[]))+[dict(source['model']['predictions'],relative_path='applied-evidence/spliceai/predictions.tsv.gz',role='APPLIED_SPLICEAI_ALLELE_PREDICTIONS',scope=codec.SCOPE)]
  named={Path(p).name:value for p,value in rule['code_sha256'].items()};producer=dict(base['producer_code_sha256'],spliceai_overlay_helper=named['apply_spliceai_allele_overlay.py'],spliceai_overlay_controller=named['continue_spliceai_allele_overlay.py'],spliceai_allele_info_codec=named['spliceai_allele_info_codec.py'],spliceai_projection=named['spliceai_packet_context.py'])
  result=dict(schema_version=1,status='PASS_FULLY_ANNOTATED_VCF_DIRECT',final_delivery_complete=True,enriched_vcf=final_pin(target),index=final_pin(Path(str(target)+'.tbi')),core_manifest=base['core_manifest'],intermediate_manifest=pin(Path(rule['base_root'])/'manifest.json'),intermediate_receipt=source['base_receipt'],spliceai_manifest=pin(Path(rule['spliceai_root'])/'manifest.json'),spliceai_receipt=source['model_receipt'],spliceai_predictions=source['model']['predictions'],schema=base['schema'],core_schema=base['core_schema'],projection=source['projection'],spliceai_storage=codec.CONTRACT,projection_fields=codec.SPLICE_FIELDS,samples=base['samples'],row_count=base['row_count'],record_count=count,applied_artifacts=applied,direct_core_prefix_verified=True,packet_decode_verified=True,original_records_fields_preserved=True,core33_packets_byte_identical=True,original_headers_preserved=True,packet_decode_verification='Inherited verified core33 packets with exact byte-token proof; no core vector re-encoding',original_record_digest=before.hexdigest(),restored_original_record_digest=restored.hexdigest(),core_packet_tokens_sha256=tokens.hexdigest(),original_headers_digest=before_headers.hexdigest(),selected_model_alleles=selected,indexed_contigs_verified=sorted(seek),producer_code_sha256=producer,core_producer_code_sha256=base['producer_code_sha256'],overlay_code_sha256=rule['code_sha256'],network_requests=0,no_legacy_receipt_created=True,synthetic_validation=public_control,created_utc=now().isoformat())
  dump(stage/'manifest.json',result);artifacts={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in stage.iterdir() if p.is_file()};identity=hashlib.sha256(json.dumps(dict(base=result['intermediate_manifest'],model=result['spliceai_manifest'],code=rule['code_sha256']),sort_keys=True).encode()).hexdigest();dump(stage/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='fully-annotated-vcf-v1',identity=identity,artifacts=artifacts,metadata=dict(projection_contract=codec.CONTRACT,core33_packets_byte_identical=True,no_legacy_receipt_created=True)))
  import ctypes
  if os.name=='nt':os.rename(stage,output)
  else:
   libc=ctypes.CDLL(None,use_errno=True)
   if libc.renameat2(-100,os.fsencode(stage),-100,os.fsencode(output),1):raise OSError(ctypes.get_errno(),'Exclusive final output promotion refused')
  return result
 finally:ACTIVE_DEADLINE=None
