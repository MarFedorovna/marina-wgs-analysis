#!/usr/bin/env python3
"""Prepare private native-download links for receipt-completed outputs over pinned SSH.

Creates only same-filesystem hardlinks inside an existing authenticated Jupyter
root. Does not change authentication, open ports, copy raw/public databases or
read browser cookies/tokens. Keep the generated plan private and outside Git.
"""
import argparse
import base64
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import shlex
import subprocess
import time
from urllib.parse import quote


HELPER = r'''
import base64,datetime,hashlib,html,json,os,pathlib,re,sys
from urllib.parse import quote
os.umask(0o077)
cfg=json.loads(sys.argv[1]);root=pathlib.Path(cfg['root']).resolve(strict=True)
compare_only=cfg.get('compare_only',False)
workspace=pathlib.Path(cfg['workspace']).resolve(strict=True);jroot=pathlib.Path(cfg['jupyter']).resolve(strict=True)
workspace.relative_to(jroot)
links=workspace/'marina-verified-downloads'
if not compare_only:links.mkdir(mode=0o700,exist_ok=True)
if links.is_symlink():raise RuntimeError('Native link directory must not be a symlink')
view=workspace/'RESULTS'
if not compare_only:view.mkdir(mode=0o700,exist_ok=True)
if view.is_symlink():raise RuntimeError('RESULTS must not be a symlink')
if not compare_only:
 for group in ('BAM','gVCF','VCF','QC','Reports'):(view/group).mkdir(mode=0o700,exist_ok=True)
files={};receipts=[];names={};view_rows={}
def human_path(p):
 parts=p.relative_to(root).parts;name=p.name
 if name.endswith(('.bam','.bam.bai','.bai')) and len(parts)>=4 and parts[-2]=='align':
  return pathlib.Path('BAM',parts[-3],name),'full-reference alignment'
 if '.g.vcf.gz' in name:
  if 'gvcf' in parts:
   i=parts.index('gvcf');scope=parts[i+1];sample=parts[i+2]
   return (pathlib.Path('gVCF',sample,name) if scope=='autosomes' else pathlib.Path('gVCF',scope,sample,name)),scope
  return pathlib.Path('Pilot',*parts[1:]),'technical pilot; separate scope'
 if '.vcf' in name or name.endswith(('.bcf','.tbi','.csi')):
  if 'stages' in parts:return pathlib.Path('VCF',name),'completed cohort stage'
  return pathlib.Path('VCF','checkpoints',*parts[1:]),'partition/checkpoint'
 if name.endswith(('.html','.pdf','.xlsx','.docx','.csv','.tsv','.csv.gz','.tsv.gz','.png','.jpg','.jpeg','.svg')):
  return pathlib.Path('Reports',*parts[1:]),'completed report/table'
 return pathlib.Path('QC',*parts[1:]),'completed QC/log/provenance'
def expose(p,entry):
 rel,scope=human_path(p);target=view/rel
 if target.is_symlink():raise RuntimeError('Human view cannot contain symlinks')
 if target.exists() and not os.path.samefile(p,target):
  rel=pathlib.Path('Versions',entry['Sha256'][:16],rel);target=view/rel
 if not compare_only:target.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
 if target.exists():
  if not os.path.samefile(p,target):raise RuntimeError('Refusing different human-view file')
 elif not compare_only:
  try:os.link(p,target)
  except FileExistsError:
   if not os.path.samefile(p,target):raise
 view_rows[entry['RelativePath']]={'ViewPath':rel.as_posix(),'Scope':scope,'Source':entry['RelativePath'],'SizeBytes':entry['SizeBytes'],'ExpectedServerSha256':entry['Sha256'],'Receipt':entry['SourceReceipt'],'DownloadFileName':entry['DownloadFileName']}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def safe_artifact(base,name):
 rel=pathlib.PurePosixPath(name)
 if not rel.parts or rel.is_absolute() or '..' in rel.parts or chr(92) in name or rel.as_posix()!=name:raise RuntimeError('Unsafe artifact path')
 p=base
 for part in rel.parts:
  p=p/part
  if p.is_symlink():raise RuntimeError('Symlink artifact or parent refused')
 if not p.is_file() or base.resolve() not in p.resolve(strict=True).parents:raise RuntimeError('Unsafe completed artifact')
 return p

def primary_receipt(stage):
 p=root/'outputs/workflow-v1/stages'/stage/'workflow-receipt.json'
 if p.is_symlink():raise RuntimeError('Symlink primary receipt')
 d=json.loads(p.read_bytes())
 if d.get('status')!='PASS' or not d.get('identity'):raise RuntimeError('Sidecar requires genuine completed primary workflow receipts')
 return p,d

def source_pin(pin,stage,name):
 rp,d=primary_receipt(stage);p=rp.parent/name;expected=d['artifacts'][name]
 if pin.get('path')!=str(p) or pin.get('sha256')!=expected['sha256'] or pin.get('bytes')!=expected['bytes']:raise RuntimeError('Sidecar source pin differs from completed primary artifact')
 if p.is_symlink() or p.stat().st_size!=expected['bytes']:raise RuntimeError('Primary artifact size changed')

def receipt_pins(pins):
 expected={}
 for stage in ('06_annotate','07_report'):
  p,d=primary_receipt(stage);expected[str(p)]=digest(p)
 if set(pins)!=set(expected) or any(pins[p].get('sha256')!=expected[p] for p in expected):raise RuntimeError('Sidecar primary receipt digest differs')

def combined_manifest(folder):
 p=folder/'manifest.json';rp=folder/'workflow-receipt.json'
 if p.is_symlink() or rp.is_symlink() or p.stat().st_size>32*1024**2:raise RuntimeError('Unsafe combined-master provenance')
 d=json.loads(p.read_bytes());receipt=json.loads(rp.read_bytes())
 required={'schema_version','status','base_master','enriched_master','original_fields','appended_fields','row_count','invariant_digest','original_rows_fields_order_preserved','ranking_eligibility_unchanged','stage_receipts','upstream_receipts','created_utc'}
 if not required<=set(d) or d.get('schema_version')!=1 or d.get('status')!='PASS_COMBINED_ANNOTATION_MASTER' or d['original_rows_fields_order_preserved'] is not True or d['ranking_eligibility_unchanged'] is not True:raise RuntimeError('Unknown combined-master schema/invariant contract')
 if receipt.get('schema_version')!=1 or receipt.get('status')!='PASS' or receipt.get('stage')!='combined-annotation-master-v1' or not re.fullmatch('[0-9a-f]{64}',receipt.get('identity','')):raise RuntimeError('Combined-master completion receipt required')
 source_pin(d['base_master'],'07_report','master_evidence.tsv.gz')
 receipt_pins({x['path']:x for x in d['stage_receipts']})
 if len(d['stage_receipts'])!=2 or type(d['row_count']) is not int or d['row_count']<0 or not re.fullmatch('[0-9a-f]{64}',d['invariant_digest']):raise RuntimeError('Invalid combined-master count/digest')
 original,appended=d['original_fields'],d['appended_fields']
 if not original or not appended or len(set(original+appended))!=len(original+appended) or any(not re.fullmatch('allref_[A-Za-z0-9_]{1,120}',f) for f in appended):raise RuntimeError('Combined master is not append-only namespaced context')
 expected_upstream={str(root/'outputs/tissue-protein-context-v1/manifest.json'),str(root/'outputs/population-reference-context-v1/workflow-receipt.json')}
 if len(d['upstream_receipts'])!=2 or {x.get('path') for x in d['upstream_receipts']}!=expected_upstream:raise RuntimeError('Unknown combined upstream receipt roots')
 for item in d['upstream_receipts']:
  source=pathlib.Path(item['path'])
  if source.is_symlink() or source.stat().st_size!=item.get('bytes') or digest(source)!=item.get('sha256'):raise RuntimeError('Combined upstream receipt digest differs')
  if sidecar_artifacts(source,json.loads(source.read_bytes())) is None:raise RuntimeError('Combined upstream is not complete')
 target=folder/'master_evidence.enriched.tsv.gz';expected=d['enriched_master']
 if expected.get('path')!=str(target):raise RuntimeError('Unexpected combined master path')
 if safe_artifact(folder,target.name).stat().st_size!=expected.get('bytes') or not re.fullmatch('[0-9a-f]{64}',expected.get('sha256','')):raise RuntimeError('Combined master artifact size/digest format differs')
 artifacts=receipt.get('artifacts',{})
 allowed={'manifest.json',target.name,'index.html','auxiliary-index.json','applied-evidence/MT-appendix-index.json'}
 population_path=root/'outputs/population-reference-context-v1/workflow-receipt.json';population=json.loads(population_path.read_bytes())
 population_pin={'path':str(population_path),'bytes':population_path.stat().st_size,'sha256':digest(population_path)}
 for entry in d.get('auxiliary_artifacts',[]):
  relative=entry.get('relative_path','');source=entry.get('source_relative_path','')
  if relative!='applied-evidence/population/'+source or source=='cohort/nuclear-sample-rows.tsv.gz' or entry.get('path')!=str(folder/relative):raise RuntimeError('Unknown combined applied evidence path')
  nuclear={'cohort/'+name+'.tsv.gz' for name in ('nuclear-population-long','nuclear-exome-AN-v41','ClinVar-validated-allele-map','ClinVar-SCV-long','ClinVar-citation-links','gene-coverage-context-long')}
  mt=source=='cohort/MT-population-reference-long.tsv.gz' or re.fullmatch(r'samples/[A-Za-z0-9][A-Za-z0-9_.-]{0,99}/(MT-population-reference\.tsv\.gz|MT-QC/[0-9]{2}-[A-Za-z0-9_.-]+\.txt)',source)
  if source not in nuclear and not mt:raise RuntimeError('Unrecognized applied result-table type')
  role='MT_APPENDIX' if mt else 'APPLIED_NUCLEAR_REFERENCE_LONG_TABLE'
  if entry.get('source_upstream_receipt')!=population_pin or entry.get('role')!=role:raise RuntimeError('Combined applied evidence producer/role differs')
  population_item=population.get('artifacts',{}).get(source,{})
  if not population_item or any(entry.get(k)!=population_item.get(k) or artifacts.get(relative,{}).get(k)!=entry.get(k) for k in ('bytes','sha256')):raise RuntimeError('Combined applied evidence differs from population completion artifact')
  allowed.add(relative)
 if not {'manifest.json',target.name}<=set(artifacts) or not set(artifacts)<=allowed:raise RuntimeError('Unknown combined-master artifacts')
 for name,pin_value in [(target.name,expected),('manifest.json',dict(bytes=p.stat().st_size,sha256=digest(p)))]:
  if any(artifacts[name].get(k)!=pin_value[k] for k in ('bytes','sha256')):raise RuntimeError('Combined receipt artifact identity differs')
 return d,receipt

def direct_vcf_manifest(folder):
 path=folder/'manifest.json';rp=folder/'workflow-receipt.json'
 d=json.loads(safe_artifact(folder,'manifest.json').read_bytes());receipt=json.loads(safe_artifact(folder,'workflow-receipt.json').read_bytes())
 if path.stat().st_size>32*1024**2 or d.get('schema_version')!=1 or d.get('status')!='PASS_FULLY_ANNOTATED_VCF_DIRECT':raise RuntimeError('Unknown direct VCF authority')
 if any(d.get(k) is not True for k in ('direct_core_prefix_verified','packet_decode_verified','original_records_fields_preserved')):raise RuntimeError('Direct VCF scientific checks incomplete')
 if receipt.get('schema_version')!=1 or receipt.get('status')!='PASS' or receipt.get('stage')!='fully-annotated-vcf-v1' or not re.fullmatch('[0-9a-f]{64}',receipt.get('identity','')):raise RuntimeError('Direct VCF completion receipt required')
 core_pin=d['core_manifest'];core_path=pathlib.Path(core_pin['path'])
 if core_path.parent not in {folder,root/'outputs/direct-core-context-v1'} or core_path.name not in {'core-manifest.json','core_manifest.json'} or core_path.is_symlink() or core_path.stat().st_size>32*1024**2 or core_path.stat().st_size!=core_pin['bytes'] or digest(core_path)!=core_pin['sha256']:raise RuntimeError('Direct baseline authority identity differs')
 core=json.loads(core_path.read_bytes())
 if core.get('status')!='DIRECT_PACKET_BASELINE_PASS' or core.get('no_legacy_receipt_created') is not True:raise RuntimeError('Direct baseline authority is not genuine new contract')
 fields=core['schema']['full_fields'];schema=d['schema']
 if schema.get('full_fields',[])[:len(fields)]!=fields or schema.get('ordinal_in_vcf_order') is not True or core['schema'].get('producer_kind')!='DIRECT_FROZEN_ANNOTATION_ROW_V1_NO_LEGACY_06_07_RECEIPTS':raise RuntimeError('Direct field-prefix/ordinal contract differs')
 for value in (schema,core['schema']):
  if value.get('vector_digest_contract')!='ordered_batch_canonical_TSV_utf8_LF_sha256_v1' or not re.fullmatch('[0-9a-f]{64}',value.get('ordered_batch_merkle_root','')):raise RuntimeError('Unknown direct batch proof contract')
  if type(value.get('row_count')) is not int or value['row_count']<1 or not value.get('batches'):raise RuntimeError('Missing direct batch counters')
  import csv,io
  header=io.StringIO(newline='');csv.writer(header,delimiter='\t',lineterminator='\n').writerow(value['full_fields'])
  descriptors=[];ordinal=1
  for batch in value['batches']:
   batch={k:batch[k] for k in ('id','first_i','last_i','row_count','record_count','canonical_row_sha256')}
   if batch['first_i']!=ordinal or batch['last_i']!=ordinal+batch['row_count']-1 or batch['row_count']<1 or batch['record_count']<1 or not re.fullmatch('[0-9a-f]{64}',batch['canonical_row_sha256']):raise RuntimeError('Invalid direct batch boundary')
   ordinal=batch['last_i']+1;descriptors.append(batch)
  if ordinal-1!=value['row_count']:raise RuntimeError('Direct batch total differs')
  commitment={'contract':value['vector_digest_contract'],'canonical_header_sha256':hashlib.sha256(header.getvalue().encode()).hexdigest(),'batches':descriptors}
  if hashlib.sha256(json.dumps(commitment,sort_keys=True,ensure_ascii=False,separators=(',',':')).encode()).hexdigest()!=value['ordered_batch_merkle_root']:raise RuntimeError('Direct metadata batch commitment differs')
 if schema['row_count']!=core['schema']['row_count']:raise RuntimeError('Direct prefix row count differs')
 source=core['source_pins'];needed={'binding','config','resources','snapshot_database','snapshot_manifest','completed_vep_vcf','completed_vep_index','vep_component_producer','filter_receipt','filter_vcf','alpha_manifest','alpha_canonical_matches','alpha_isoform_matches','core_public_proof','core_module'}
 if set(source)!=needed:raise RuntimeError('Unknown direct source lineage')
 source_pin(source['filter_vcf'],'05_filter','cohort.filtered.vcf.gz')
 p,_=primary_receipt('05_filter')
 if source['filter_receipt']!={'path':str(p),'bytes':p.stat().st_size,'sha256':digest(p)}:raise RuntimeError('Direct accepted05 receipt differs')
 for name,item in source.items():
  p=pathlib.Path(item['path'])
  if p.is_symlink() or root not in p.resolve(strict=True).parents or p.stat().st_size!=item['bytes'] or not re.fullmatch('[0-9a-f]{64}',item.get('sha256','')):raise RuntimeError('Unsafe direct source pin')
 producer=root/'private/direct-evidence-candidate-v1/scripts/fast_vcf_evidence_candidate.py'
 if not producer.is_file() or core.get('producer_code_sha256',{}).get('core_module')!=digest(producer) or any(d.get('producer_code_sha256',{}).get(k)!=v for k,v in core.get('producer_code_sha256',{}).items()):raise RuntimeError('Direct reviewed producer identity differs')
 if d.get('projection'):
  if d.get('core_producer_code_sha256')!=core.get('producer_code_sha256') or d.get('core33_packets_byte_identical') is not True:raise RuntimeError('SpliceAI overlay changed the stored core33 authority')
  names={'apply_spliceai_allele_overlay.py','continue_spliceai_allele_overlay.py','spliceai_allele_info_codec.py','spliceai_packet_context.py','embed_vcf_evidence.py'}
  matched=False
  for version in ('v1','v2','v3'):
   package=root/('private/spliceai-allele-overlay-candidate-'+version+'/scripts')
   if all((package/n).is_file() for n in names):
    expected={str(package/n):digest(package/n) for n in names}
    if d.get('overlay_code_sha256')==expected:matched=True
  if not matched:raise RuntimeError('Exact reviewed allele overlay producer differs')
 for key,name in [('enriched_vcf','cohort.fully-annotated.vcf.gz'),('index','cohort.fully-annotated.vcf.gz.tbi')]:
  item=d[key];p=safe_artifact(folder,name)
  if item.get('path')!=str(p) or item.get('bytes')!=p.stat().st_size or any(receipt.get('artifacts',{}).get(name,{}).get(k)!=item[k] for k in ('bytes','sha256')):raise RuntimeError('Direct VCF receipt artifact differs')
 if receipt.get('artifacts',{}).get('manifest.json')!={'bytes':path.stat().st_size,'sha256':digest(path)}:raise RuntimeError('Direct manifest is not receipt-bound')
 return d,receipt,core

def fast_view_artifacts(rp,d):
 if d.get('schema_version')!=1 or d.get('status')!='PASS' or d.get('synthetic') is not False:return None
 if d.get('report_input_kind')!='FULLY_ANNOTATED_VCF_PACKETS_V1' or any(d.get(k) is not True for k in ('sample_packet_scope_verified','public_header_bundles_verified','no_legacy_receipt_created','no_external_requests')) or d.get('automatic_ACMG') is not False or d.get('primary_workflow_changed') is not False:raise RuntimeError('Fast views contract differs')
 fully,_,core=direct_vcf_manifest(root/'outputs/fully-annotated-vcf-v1')
 expected={'fully_manifest':{'path':str(root/'outputs/fully-annotated-vcf-v1/manifest.json'),'bytes':(root/'outputs/fully-annotated-vcf-v1/manifest.json').stat().st_size,'sha256':digest(root/'outputs/fully-annotated-vcf-v1/manifest.json')},'core_manifest':fully['core_manifest'],'vcf':fully['enriched_vcf']}
 if d.get('source_pins')!=expected or d.get('batch_merkle_root')!=fully['schema']['ordered_batch_merkle_root'] or d.get('core_prefix_merkle_root')!=core['schema']['ordered_batch_merkle_root']:raise RuntimeError('Fast views source proof differs')
 version=d.get('report_producer_kind')=='FAST_DIRECT_PACKET_VIEWS_V2'
 producer=root/('private/reviewed-packet-fastviews-v2/scripts/export_packet_views_fast_v2.py' if version else 'private/reviewed-packet-fastviews-v1/scripts/export_packet_views_fast.py')
 producers=[producer]
 if version:producers.append(root/'private/reviewed-packet-fastviews-v2-3/scripts/export_packet_views_fast_v2_1.py')
 ranking=root/'private/annotation-runtime-v1/scripts/annotation_stages.py'
 if not any(p.is_file() and d.get('producer_code_sha256')==digest(p) for p in producers) or d['ranking_module'].get('path')!=str(ranking) or d['ranking_module'].get('sha256')!=digest(ranking) or d['ranking_module'].get('function')!='rank_key':raise RuntimeError('Fast views frozen producer/rank differs')
 samples=d.get('samples',[]);labels=d.get('sample_file_labels',{})
 if len(set(samples))!=4 or len(samples)!=4 or set(labels)!=set(samples) or set(labels.values())!={'sample_0001','sample_0002','sample_0003','sample_0004'}:raise RuntimeError('Fast views four-sample scope differs')
 artifacts=dict(d['artifacts']);required={'index.html','masters/cohort.master.tsv.gz','tables/cohort.top20.tsv','tables/cohort.functional.tsv.gz','tables/functional-category-index.tsv','tables/allele-membership-matrix.tsv.gz','tables/ACMG-manual-top20-worksheet.tsv.gz','tables/gene-overlap-SEPARATE-NOT-SHARED-ALLELES.tsv.gz','vcf/cohort.fully-annotated.vcf.gz','vcf/cohort.fully-annotated.vcf.gz.tbi','provenance/fully-annotated-vcf-manifest.json','provenance/direct-core-manifest.json'}
 for label in labels.values():required.update({'masters/'+label+'.master.tsv.gz','tables/'+label+'.top20.tsv','tables/'+label+'.functional.tsv.gz','vcf/'+label+'.annotated.vcf.gz','vcf/'+label+'.annotated.vcf.gz.tbi','tables/spliceai/'+label+'.top20.tsv'})
 required.add('tables/spliceai/cohort.top20.tsv')
 if not required<=set(artifacts):raise RuntimeError('Missing final canonical view artifacts')
 for key,name in [('enriched_vcf','cohort.fully-annotated.vcf.gz'),('index','cohort.fully-annotated.vcf.gz.tbi')]:
  if any(artifacts['vcf/'+name].get(k)!=fully[key][k] for k in ('bytes','sha256')):raise RuntimeError('Canonical final VCF copy differs')
 permitted=lambda n:n in required or bool(re.fullmatch(r'tables/(cohort|shared_ge2|all4|sample_000[1-4])\.(functional\.tsv\.gz|top20\.tsv)|tables/(cohort|shared_ge2|all4|sample_000[1-4])/[A-Za-z0-9_-]+\.tsv\.gz|tables/spliceai/(cohort|sample_000[1-4])\.scored\.tsv\.gz',n))
 applied=fully.get('applied_artifacts',[])
 for item in applied:
  relative='applied-evidence/spliceai/predictions.tsv.gz';source=root/'outputs/spliceai-reference-context-v1/predictions.tsv.gz'
  if item.get('relative_path')!=relative or item.get('role')!='APPLIED_SPLICEAI_ALLELE_PREDICTIONS' or item.get('scope')!='GENCODE24_GENE_UNION_NOT_ENST49_MAPPING' or item.get('path')!=str(source):raise RuntimeError('Unknown applied output/resource source')
  if any(item.get(k)!=fully.get('spliceai_predictions',{}).get(k) or artifacts.get(relative,{}).get(k)!=item.get(k) for k in ('bytes','sha256')):raise RuntimeError('Applied SpliceAI prediction identity differs')
  receipt_source=root/'outputs/spliceai-reference-context-v1/workflow-receipt.json';sp=json.loads(receipt_source.read_bytes())
  if sp.get('status')!='PASS' or sp.get('stage')!='spliceai-reference-context-v1' or any(sp.get('artifacts',{}).get('predictions.tsv.gz',{}).get(k)!=item[k] for k in ('bytes','sha256')):raise RuntimeError('Applied SpliceAI producer receipt differs')
  required.add(relative)
 for name,item in artifacts.items():
  if not permitted(name) or not isinstance(item,dict) or type(item.get('bytes')) is not int or not re.fullmatch('[0-9a-f]{64}',item.get('sha256','')):raise RuntimeError('Unknown fast output/resource artifact')
  if safe_artifact(rp.parent,name).stat().st_size!=item['bytes']:raise RuntimeError('Fast completed output size changed')
 receipt_path=rp.parent/'workflow-receipt.json';receipt=json.loads(safe_artifact(rp.parent,'workflow-receipt.json').read_bytes())
 expected_receipt=dict(artifacts,**{'manifest.json':{'bytes':rp.stat().st_size,'sha256':digest(rp)}})
 if receipt.get('schema_version')!=1 or receipt.get('status')!='PASS' or receipt.get('stage')!='cohort-sample-views-v1' or receipt.get('artifacts')!=expected_receipt or receipt.get('producer_code_sha256')!=d['producer_code_sha256']:raise RuntimeError('Fast completed view receipt differs')
 # The active immutable importer binds each artifact to SourceReceipt. The
 # manifest cannot list its subsequently-created workflow receipt as an output.
 # Preserve that authentic receipt in captured metadata, without inventing a
 # self-referential manifest artifact or changing the running importer.
 raw=receipt_path.read_bytes();receipts.append({'RelativePath':receipt_path.relative_to(root).as_posix(),'SizeBytes':len(raw),'Sha256':hashlib.sha256(raw).hexdigest(),'ContentBase64':base64.b64encode(raw).decode()})
 artifacts.update({'manifest.json':expected_receipt['manifest.json']})
 return artifacts

def fully_vcf_manifest(folder):
 path=folder/'manifest.json';rp=folder/'workflow-receipt.json'
 if path.is_symlink() or rp.is_symlink() or path.stat().st_size>32*1024**2:raise RuntimeError('Unsafe fully annotated VCF provenance')
 d=json.loads(path.read_bytes());receipt=json.loads(rp.read_bytes());combined,joined_receipt=combined_manifest(root/'outputs/combined-annotation-master-v1')
 if d.get('schema_version')!=1 or d.get('status')!='PASS_FULLY_ANNOTATED_VCF' or d.get('network_requests')!=0 or d.get('primary_outputs_changed') is not False:raise RuntimeError('Unknown fully annotated VCF schema/contract')
 if receipt.get('schema_version')!=1 or receipt.get('status')!='PASS' or receipt.get('stage')!='fully-annotated-vcf-v1' or not re.fullmatch('[0-9a-f]{64}',receipt.get('identity','')):raise RuntimeError('Fully annotated VCF completion receipt required')
 if any(d.get(flag) is not True for flag in ('original_records_fields_preserved','packet_decode_verified','decoded_master_vectors_verified','public_header_bundle_verified')):raise RuntimeError('Fully annotated VCF lossless/public-bundle checks incomplete')
 source_pin(d['base_vcf'],'06_annotate','annotated.vcf.gz');joined=root/'outputs/combined-annotation-master-v1'
 for key,file in [('master_manifest',joined/'manifest.json'),('combined_receipt',joined/'workflow-receipt.json')]:
  if d.get(key)!={'path':str(file),'bytes':file.stat().st_size,'sha256':digest(file)}:raise RuntimeError('Fully annotated VCF joined lineage differs')
 if d.get('master')!=combined['enriched_master']:raise RuntimeError('Fully annotated control master differs')
 schema=d['schema']
 if schema.get('v')!=1 or schema.get('codec')!='json-utf8+zlib+base64url-nopad' or schema.get('scope')!='cohort' or schema.get('sample') is not None:raise RuntimeError('Unknown fully annotated cohort packet schema')
 for key in ('original_fields','appended_fields','row_count','invariant_digest'):
  if d.get(key)!=combined[key] or schema.get(key)!=combined[key]:raise RuntimeError('Fully annotated original-row invariant differs')
 if d.get('full_fields')!=combined['original_fields']+combined['appended_fields'] or schema.get('full_fields')!=d['full_fields'] or not re.fullmatch('[0-9a-f]{64}',d.get('full_fields_digest','')) or schema.get('full_fields_digest')!=d['full_fields_digest']:raise RuntimeError('Fully annotated complete-row schema/digest differs')
 for key,name in [('enriched_vcf','cohort.fully-annotated.vcf.gz'),('index','cohort.fully-annotated.vcf.gz.tbi')]:
  item=d[key];file=safe_artifact(folder,name)
  if item.get('path')!=str(file) or item.get('bytes')!=file.stat().st_size or not re.fullmatch('[0-9a-f]{64}',item.get('sha256','')):raise RuntimeError('Fully annotated VCF artifact identity differs')
  if any(receipt.get('artifacts',{}).get(name,{}).get(k)!=item[k] for k in ('bytes','sha256')):raise RuntimeError('Fully annotated VCF receipt artifact differs')
 expected={'cohort.fully-annotated.vcf.gz','cohort.fully-annotated.vcf.gz.tbi','manifest.json'}
 if set(receipt.get('artifacts',{}))!=expected or receipt['artifacts']['manifest.json']!={'bytes':path.stat().st_size,'sha256':digest(path)}:raise RuntimeError('Unknown fully annotated VCF output set')
 producers=[root/'private'/package/'scripts/embed_vcf_evidence.py' for package in ('fully-annotated-vcf-candidate-v1','reviewed-fully-annotated-vcf-v1')]
 if not any(p.is_file() and d.get('producer_code_sha256',{}).get(str(p))==digest(p) for p in producers):raise RuntimeError('Fully annotated VCF producer identity differs')
 return d,receipt,combined,joined_receipt

def sidecar_artifacts(rp,d):
 name=rp.parent.name
 if (rp.parent/'INCOMPLETE.txt').exists():raise RuntimeError('Completed sidecar has INCOMPLETE marker')
 if name=='cohort-sample-views-v1' and d.get('report_producer_kind') in ('FAST_DIRECT_PACKET_VIEWS_V1','FAST_DIRECT_PACKET_VIEWS_V2'):return fast_view_artifacts(rp,d)
 if name=='consequence-tables-v1':
  if d.get('status')!='PASS' or d.get('schema_version')!=1 or d.get('synthetic') is not False:return None
  producer=root/'private/reviewed-consequence-export-v1/scripts/export_consequence_tables.py'
  if not producer.is_file() or d.get('exporter_sha256')!=digest(producer):raise RuntimeError('Consequence producer identity differs')
  source_pin(d['source'],'06_annotate','transcript_samples.tsv.gz')
  if d.get('workflow_receipt')!=primary_receipt('06_annotate')[1]:raise RuntimeError('Consequence embedded primary receipt differs')
  artifacts=dict(d['artifacts'])
  permitted=lambda n:n in ('index.tsv','index.html','continuation-receipts.json','evidence/transcript_effects.tsv.gz') or bool(re.fullmatch(r'tables/(cohort|sample_[0-9]{4})/[A-Za-z0-9_-]+\.tsv\.gz',n))
 elif name=='tissue-protein-context-v1':
  if d.get('status')!='PASS_SEPARATE_RESEARCH_REFERENCE_CONTEXT_EXPORT':return None
  keys={'status','source_artifact','source_sha256','source_bytes','context_index_sha256','stage_receipts','counters','artifact','scope','original_pipeline_changed','network_requests','created_utc'}
  if set(d)!=keys or d['original_pipeline_changed'] is not False or d['network_requests']!=0:raise RuntimeError('Unknown tissue sidecar schema/contract')
  context=root/'resources/tissue-protein-context-public39-v1/context-index.json'
  if not context.is_file() or context.is_symlink() or d['context_index_sha256']!=digest(context):raise RuntimeError('Public tissue context identity differs; source is not copied')
  source_pin(dict(path=d['source_artifact'],sha256=d['source_sha256'],bytes=d['source_bytes']),'06_annotate','transcript_samples.tsv.gz');receipt_pins(d['stage_receipts'])
  target=rp.parent/'public39-reference-context.tsv.gz'
  if d['artifact'].get('path')!=str(target):raise RuntimeError('Unexpected tissue artifact path')
  artifacts={'public39-reference-context.tsv.gz':{k:d['artifact'][k] for k in ('bytes','sha256')}}
  permitted=lambda n:n=='public39-reference-context.tsv.gz'
 elif name=='cohort-sample-views-v1':
  if d.get('status')!='PASS' or d.get('schema_version')!=1 or d.get('synthetic') is not False:return None
  if d.get('primary_workflow_changed') is not False or d.get('no_external_requests') is not True or d.get('automatic_ACMG') is not False:raise RuntimeError('Unknown cohort sidecar contract')
  ranking=root/'private/annotation-runtime-v1/scripts/annotation_stages.py'
  if d['ranking_module'].get('path')!=str(ranking) or d['ranking_module'].get('sha256')!=digest(ranking) or d['ranking_module'].get('function')!='rank_key':raise RuntimeError('Sealed ranking producer pin differs')
  expected_paths=set()
  for stage,artifact in [('06_annotate','annotated.vcf.gz'),('06_annotate','annotated.vcf.gz.tbi'),('07_report','master_evidence.tsv.gz')]:
   p=root/'outputs/workflow-v1/stages'/stage/artifact;expected_paths.add(str(p));source_pin(dict(d['source_pins'][str(p)],path=str(p)),stage,artifact)
  pins={str(p):d['source_pins'][str(p)] for p in [root/'outputs/workflow-v1/stages'/s/'workflow-receipt.json' for s in ('06_annotate','07_report')]}
  receipt_pins(pins);expected_paths.update(pins)
  if d.get('report_input_kind')=='COMBINED_ANNOTATION_MASTER_V1':
   folder=root/'outputs/combined-annotation-master-v1';combined,joined_receipt=combined_manifest(folder)
   meta=d.get('joined_master',{});manifest=folder/'manifest.json'
   if meta.get('manifest')!={'path':str(manifest),'bytes':manifest.stat().st_size,'sha256':digest(manifest)}:raise RuntimeError('Cohort joined manifest identity differs')
   for field in ('enriched_master','base_master','original_fields','appended_fields','row_count','invariant_digest'):
    if meta.get(field)!=combined[field]:raise RuntimeError('Cohort joined provenance differs: '+field)
   if meta.get('original_rows_fields_order_verified') is not True or meta.get('rank_eligibility_columns_unchanged_verified') is not True:raise RuntimeError('Cohort joined invariants were not verified')
   engine=root/'private/reviewed-cohort-sample-views-v1/scripts/export_cohort_sample_views.py'
   producers=[root/'private'/version/'scripts/export_cohort_sample_views_v2.py' for version in ('reviewed-cohort-sample-views-v2','reviewed-cohort-sample-views-v2-1')]
   if d.get('views_v1_engine')!={'path':str(engine),'sha256':digest(engine)} or not any(p.is_file() and d.get('views_v2_exporter_sha256')==digest(p) for p in producers):raise RuntimeError('Joined views producer/engine identity differs')
   applied={entry['relative_path'] for entry in combined.get('auxiliary_artifacts',[])}
   if applied:applied.update({'auxiliary-index.json','applied-evidence/MT-appendix-index.json'})
   if set(d.get('copied_applied_evidence',[]))!=applied:raise RuntimeError('Joined views applied evidence set differs')
   for relative in applied:
    if any(d.get('artifacts',{}).get(relative,{}).get(k)!=joined_receipt['artifacts'][relative].get(k) for k in ('bytes','sha256')):raise RuntimeError('Joined views applied evidence artifact differs')
   extras=[dict(path=str(manifest),bytes=manifest.stat().st_size,sha256=digest(manifest)),combined['enriched_master'],
           dict(path=str(folder/'workflow-receipt.json'),bytes=(folder/'workflow-receipt.json').stat().st_size,sha256=digest(folder/'workflow-receipt.json')),*combined['upstream_receipts']]
   for item in extras:
    path=item['path'];expected_paths.add(path)
    if any(d['source_pins'].get(path,{}).get(k)!=item[k] for k in ('bytes','sha256')):raise RuntimeError('Cohort joined source pin differs')
  elif d.get('report_input_kind')=='FULLY_ANNOTATED_VCF_PACKETS_V1':
   folder=root/'outputs/fully-annotated-vcf-v1';fully,fully_receipt,combined,joined_receipt=fully_vcf_manifest(folder)
   meta=d.get('fully_annotated_vcf',{});manifest=folder/'manifest.json'
   if meta.get('manifest')!={'path':str(manifest),'bytes':manifest.stat().st_size,'sha256':digest(manifest)} or meta.get('vcf')!=fully['enriched_vcf'] or meta.get('index')!=fully['index'] or meta.get('schema')!=fully['schema']:raise RuntimeError('VCF-first report source differs')
   control=d.get('joined_control',{})
   if control.get('manifest')!=fully['master_manifest'] or control.get('master')!=fully['master'] or control.get('all_full_field_vectors_order_verified') is not True or control.get('used_as_report_input') is not False:raise RuntimeError('VCF-first joined control contract differs')
   if d.get('sample_packet_scope_verified') is not True or d.get('public_header_bundles_verified') is not True:raise RuntimeError('VCF-first individual/public scope not verified')
   checks=d.get('sample_packet_checks',[])
   if len(checks)!=4 or {x.get('sample') for x in checks}!=set(d['samples']) or any(x.get('scope')!='sample' or x.get('other_sample_rows')!=0 or x.get('original_source_ordinals_retained') is not True for x in checks):raise RuntimeError('Per-sample VCF packet pruning not proven')
   engine=root/'private/reviewed-cohort-sample-views-v1/scripts/export_cohort_sample_views.py';helper=root/'private/reviewed-cohort-sample-views-v2-1/scripts/export_cohort_sample_views_v2.py';producer=root/'private/reviewed-cohort-sample-views-v3/scripts/export_cohort_sample_views_v3.py'
   if d.get('views_v1_engine')!={'path':str(engine),'sha256':digest(engine)} or d.get('views_v2_control_helper')!={'path':str(helper),'sha256':digest(helper)} or d.get('views_v3_exporter_sha256')!=digest(producer):raise RuntimeError('VCF-first reviewed report producer differs')
   codec=d.get('packet_codec',{})
   if fully['producer_code_sha256'].get(codec.get('path'))!=codec.get('sha256'):raise RuntimeError('VCF-first packet codec differs from source producer')
   extracted=d.get('extracted_master',{});extracted_path=pathlib.Path(extracted.get('path',''))
   try:relative_extracted=extracted_path.relative_to(root/'private/cohort-sample-views-continuation-v3')
   except ValueError:raise RuntimeError('VCF-derived master root differs')
   if len(relative_extracted.parts)!=2 or not re.fullmatch('extraction-[0-9a-f]{32}',relative_extracted.parts[0]) or relative_extracted.parts[1]!='master.from-vcf.tsv.gz' or extracted_path.is_symlink() or extracted_path.stat().st_size!=extracted.get('bytes') or not re.fullmatch('[0-9a-f]{64}',extracted.get('sha256','')):raise RuntimeError('VCF-derived master pin/path differs')
   if extracted.get('source')!='VCF_INFO_WGS_EVIDENCE_ONLY' or extracted.get('row_count')!=fully['row_count'] or extracted.get('full_fields_digest')!=fully['full_fields_digest'] or extracted.get('original_fields_digest')!=fully['invariant_digest']:raise RuntimeError('VCF-derived complete master invariant differs')
   if any(d.get('artifacts',{}).get('masters/cohort.master.tsv.gz',{}).get(k)!=extracted[k] for k in ('bytes','sha256')):raise RuntimeError('Canonical master differs from VCF packet extraction')
   applied={entry['relative_path'] for entry in combined.get('auxiliary_artifacts',[])}
   if applied:applied.update({'auxiliary-index.json','applied-evidence/MT-appendix-index.json'})
   if set(d.get('copied_applied_evidence',[]))!=applied:raise RuntimeError('VCF-first applied evidence set differs')
   for relative in applied:
    if any(d.get('artifacts',{}).get(relative,{}).get(k)!=joined_receipt['artifacts'][relative].get(k) for k in ('bytes','sha256')):raise RuntimeError('VCF-first applied evidence artifact differs')
   joined=root/'outputs/combined-annotation-master-v1'
   extras=[fully['master_manifest'],fully['combined_receipt'],fully['master'],*combined['upstream_receipts'],
           dict(path=str(manifest),bytes=manifest.stat().st_size,sha256=digest(manifest)),
           dict(path=str(folder/'workflow-receipt.json'),bytes=(folder/'workflow-receipt.json').stat().st_size,sha256=digest(folder/'workflow-receipt.json')),
           fully['enriched_vcf'],fully['index'],extracted]
   for item in extras:
    source=item['path'];expected_paths.add(source)
    if any(d['source_pins'].get(source,{}).get(k)!=item[k] for k in ('bytes','sha256')):raise RuntimeError('VCF-first source pin differs')
   for key,name in [('enriched_vcf','cohort.fully-annotated.vcf.gz'),('index','cohort.fully-annotated.vcf.gz.tbi')]:
    if any(d.get('artifacts',{}).get('vcf/'+name,{}).get(k)!=fully[key][k] for k in ('bytes','sha256')):raise RuntimeError('Canonical cohort VCF copy differs')
   provenance={'fully-annotated-vcf-manifest':dict(path=str(manifest),bytes=manifest.stat().st_size,sha256=digest(manifest)),
       'fully-annotated-vcf-receipt':dict(path=str(folder/'workflow-receipt.json'),bytes=(folder/'workflow-receipt.json').stat().st_size,sha256=digest(folder/'workflow-receipt.json')),
       'joined-control-manifest':fully['master_manifest'],'joined-control-receipt':fully['combined_receipt']}
   for name,item in provenance.items():
    if any(d.get('artifacts',{}).get('provenance/'+name+'.json',{}).get(k)!=item[k] for k in ('bytes','sha256')):raise RuntimeError('Canonical VCF-first provenance copy differs')
  elif d.get('report_input_kind') is not None:raise RuntimeError('Unknown cohort report input kind')
  if set(d['source_pins'])!=expected_paths:raise RuntimeError('Unknown cohort source pin')
  labels=d.get('sample_file_labels',{})
  if len(d.get('samples',[]))!=4 or len(set(d['samples']))!=4 or set(labels)!=set(d['samples']) or set(labels.values())!={'sample_0001','sample_0002','sample_0003','sample_0004'}:raise RuntimeError('Cohort sidecar must contain exact four sample views')
  artifacts=dict(d['artifacts']);scope=r'(cohort|shared_ge2|all4|sample_000[1-4])'
  applied=set(d.get('copied_applied_evidence',[])) if d.get('report_input_kind') in ('COMBINED_ANNOTATION_MASTER_V1','FULLY_ANNOTATED_VCF_PACKETS_V1') else set()
  extra_v3={'vcf/cohort.fully-annotated.vcf.gz','vcf/cohort.fully-annotated.vcf.gz.tbi',*['provenance/'+x+'.json' for x in ('fully-annotated-vcf-manifest','fully-annotated-vcf-receipt','joined-control-manifest','joined-control-receipt')]} if d.get('report_input_kind')=='FULLY_ANNOTATED_VCF_PACKETS_V1' else set()
  permitted=lambda n:n in applied or n in extra_v3 or n in ('index.html','progress.json','tables/allele-membership-matrix.tsv.gz','tables/ACMG-manual-top20-worksheet.tsv.gz','tables/gene-overlap-SEPARATE-NOT-SHARED-ALLELES.tsv.gz','tables/functional-category-index.tsv') or bool(re.fullmatch(r'masters/(cohort|sample_000[1-4])\.master\.tsv\.gz|vcf/sample_000[1-4]\.annotated\.vcf\.gz(\.tbi)?|tables/'+scope+r'\.(functional\.tsv\.gz|top20\.tsv)|tables/'+scope+r'/[A-Za-z0-9_-]+\.tsv\.gz',n))
  required={'vcf/sample_%04d.annotated.vcf.gz%s'%(i,suffix) for i in range(1,5) for suffix in ('','.tbi')}
  if not required<=set(artifacts):raise RuntimeError('Missing sample annotated VCF/TBI pair')
  if extra_v3 and not extra_v3<=set(artifacts):raise RuntimeError('Missing self-contained cohort VCF/provenance artifacts')
 elif name=='fully-annotated-vcf-v1':
  if d.get('status')!='PASS':return None
  direct=json.loads((rp.parent/'manifest.json').read_bytes())
  if direct.get('status')=='PASS_FULLY_ANNOTATED_VCF_DIRECT':
   fully,receipt,core=direct_vcf_manifest(rp.parent)
   return {'manifest.json':{'bytes':(rp.parent/'manifest.json').stat().st_size,'sha256':digest(rp.parent/'manifest.json')},'workflow-receipt.json':{'bytes':rp.stat().st_size,'sha256':digest(rp)}}
  fully,receipt,combined,joined_receipt=fully_vcf_manifest(rp.parent);artifacts=dict(receipt['artifacts']);permitted=lambda n:n in receipt['artifacts']
 elif name=='combined-annotation-master-v1':
  if d.get('status')!='PASS':return None
  combined,receipt=combined_manifest(rp.parent);artifacts=dict(receipt['artifacts'])
  permitted=lambda n:n in receipt['artifacts']
 elif name=='population-reference-context-v1':
  if d.get('status')!='PASS' or d.get('schema_version')!=1 or d.get('stage')!=name or not re.fullmatch('[0-9a-f]{64}',d.get('identity','')):return None
  meta=d['metadata']
  if meta.get('status')!='PASS_LOCAL_POPULATION_MT_CONTEXT_EXPORT' or meta.get('network_requests')!=0 or meta.get('primary_outputs_changed') is not False or meta.get('automatic_ACMG') is not False:raise RuntimeError('Unknown population sidecar schema/contract')
  source_pin(meta['source_master'],'07_report','master_evidence.tsv.gz')
  receipt_pins({x['path']:x for x in meta['primary_receipts']})
  artifacts=dict(d['artifacts'])
  cohort={'nuclear-sample-rows','nuclear-population-long','nuclear-exome-AN-v41','ClinVar-validated-allele-map','ClinVar-SCV-long','ClinVar-citation-links','gene-coverage-context-long','MT-population-reference-long'}
  permitted=lambda n:n in ('README.md','provenance.json') or n in {'cohort/'+s+'.tsv.gz' for s in cohort} or bool(re.fullmatch(r'samples/[A-Za-z0-9][A-Za-z0-9_.-]{0,99}/(nuclear-reference-context\.tsv\.gz|SCV-reference-links\.tsv\.gz|MT-population-reference\.tsv\.gz|MT-QC/[0-9]{2}-[A-Za-z0-9_.-]+\.txt)',n))
 elif name=='mitochondrial-reference-context-v1':
  if d.get('schema_version')!=1 or d.get('status')!='PASS' or d.get('stage')!=name or not re.fullmatch('[0-9a-f]{64}',d.get('identity','')):return None
  mp=safe_artifact(rp.parent,'manifest.json');meta=json.loads(mp.read_bytes())
  if meta.get('schema_version')!=1 or meta.get('status')!='PASS_INDEPENDENT_MITOCHONDRIAL_REFERENCE_CONTEXT' or meta.get('network_requests')!=0 or any(meta.get(k) is not True for k in ('no_nuclear_dependency','no_automatic_ACMG')) or meta.get('primary_and_original_MT_changed') is not False or meta.get('clinical_validation') is not False:raise RuntimeError('Unknown independent MT contract')
  package=root/'private/mitochondrial-reference-context-candidate-v1/scripts';names=('export_mitochondrial_reference_context.py','export_population_mt_reference_context.py')
  expected={str(package/n):digest(package/n) for n in names}
  if meta.get('producer_code_sha256')!=expected or expected[str(package/names[1])]!='4fc954c6cacce92a199fa477fcb75a3c662ecb1fcd0b0dce0ef65b6e13301832':raise RuntimeError('Independent MT reviewed producer differs')
  for key in ('source_population_rule','resource_ready_manifest','mt_audit'):
   item=meta[key];source=pathlib.Path(item['path'])
   if source.is_symlink() or root not in source.resolve(strict=True).parents or source.stat().st_size>32*1024**2 or source.stat().st_size!=item['bytes'] or digest(source)!=item['sha256']:raise RuntimeError('Independent MT source metadata identity differs')
  if meta['resource_ready_manifest']['sha256']!='465ee667cb78baddc3665d06adb15fb27ad54c47415ff14238fda5f0f42c3d85' or meta['source_database'].get('sha256')!='f069164654316043d8a96f47ad5b29413d7a28b0b9ba320648206110e42739b0':raise RuntimeError('Independent MT reference release differs')
  audit=json.loads(pathlib.Path(meta['mt_audit']['path']).read_bytes())
  if audit.get('status')!='PASS_RESEARCH_WITH_LIMITATIONS':raise RuntimeError('Current MT-v2 technical audit required')
  samples=meta.get('samples',[]);inputs=meta.get('source_mt_inputs',[])
  if len(samples)!=4 or len(set(samples))!=4 or len(inputs)!=4 or {x.get('sample') for x in inputs}!=set(samples) or any(not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]{0,99}',s) for s in samples):raise RuntimeError('Exact four independent MT sample sources required')
  if meta.get('source_mt_receipts')!=[e['receipt'] for e in inputs] or meta.get('source_mt_vcfs')!=[e['vcf'] for e in inputs]:raise RuntimeError('Independent MT source aliases differ')
  receipts_by_sample={};originals={}
  for entry in inputs:
   item=entry['receipt'];source=pathlib.Path(item['path'])
   if source.is_symlink() or root/'outputs/mitochondria-recovery-v2' not in source.resolve(strict=True).parents or source.stat().st_size!=item['bytes'] or digest(source)!=item['sha256']:raise RuntimeError('Current MT-v2 source receipt differs')
   receipt=json.loads(source.read_bytes())
   if receipt.get('status')!='PASS' or not receipt.get('identity'):raise RuntimeError('Current MT-v2 source is not PASS')
   receipts_by_sample[entry['sample']]=item;originals[entry['sample']]=[entry['vcf'],*entry['qc'],*entry['indices']]
   for original in originals[entry['sample']]:
    path=pathlib.Path(original['path'])
    if path.is_symlink() or root/'outputs/mitochondria-recovery-v2' not in path.resolve(strict=True).parents or path.stat().st_size!=original['bytes']:raise RuntimeError('Unexpected independent MT original source')
    if path.parent==source.parent and any(receipt.get('artifacts',{}).get(path.name,{}).get(k)!=original[k] for k in ('bytes','sha256')):raise RuntimeError('MT original source differs from call receipt')
  artifacts=dict(d['artifacts']);allowed={'manifest.json','index.html','cohort/MT-population-reference-long.tsv.gz'}
  for sample in samples:allowed.add('samples/'+sample+'/MT-population-reference.tsv.gz')
  required_mt=set(allowed);copied_originals={s:[] for s in samples}
  for entry in meta.get('applied_artifacts',[]):
   relative=entry.get('relative_path','').removeprefix('applied-evidence/MT/');parts=pathlib.PurePosixPath(relative).parts
   if entry.get('relative_path')!='applied-evidence/MT/'+relative or entry.get('path')!=str(rp.parent/relative):raise RuntimeError('Independent MT applied artifact path differs')
   own=parts[1] if len(parts)>2 and parts[0]=='samples' else None
   if own is not None and own not in samples:raise RuntimeError('Unknown MT sample output')
   expected_receipts=[receipts_by_sample[own]] if own else [e['receipt'] for e in inputs]
   if entry.get('source_receipts')!=expected_receipts or any(artifacts.get(relative,{}).get(k)!=entry.get(k) for k in ('bytes','sha256')):raise RuntimeError('Independent MT applied source/output identity differs')
   role=entry.get('role')
   if role in ('MT_SOURCE_QC','MT_SOURCE_VCF'):
    original=entry.get('source_original_pin')
    if own is None or original not in originals[own] or any(entry[k]!=original[k] for k in ('bytes','sha256')):raise RuntimeError('Copied independent MT original changed')
    copied_originals[own].append(original)
    pattern=r'samples/'+re.escape(own)+r'/(MT-QC/[0-9]{2}-[A-Za-z0-9_.-]+\.txt|mitochondria\.filtered\.vcf\.gz(\.tbi|\.csi)?)'
    if not re.fullmatch(pattern,relative):raise RuntimeError('Unknown independent MT VCF/QC type')
   elif role!='MT_POPULATION_CONTEXT' or relative not in allowed:raise RuntimeError('Unknown independent MT applied role')
   allowed.add(relative)
  if not required_mt<=set(artifacts) or any(sorted(json.dumps(p,sort_keys=True) for p in copied_originals[s])!=sorted(json.dumps(p,sort_keys=True) for p in originals[s]) for s in samples):raise RuntimeError('Independent MT omitted a required table/original QC/VCF/index')
  if artifacts.get('manifest.json')!={'bytes':mp.stat().st_size,'sha256':digest(mp)}:raise RuntimeError('Independent MT manifest is not receipt-bound')
  permitted=lambda n:n in allowed
 else:raise RuntimeError('Unrecognized sidecar root')
 for n,v in artifacts.items():
  if not permitted(n):raise RuntimeError('Unrecognized sidecar output/resource path: '+n)
  if not isinstance(v,dict) or type(v.get('bytes')) is not int or v['bytes']<0 or not re.fullmatch('[0-9a-f]{64}',v.get('sha256','')):raise RuntimeError('Invalid sidecar artifact identity')
  p=safe_artifact(rp.parent,n)
  if p.stat().st_size!=v['bytes']:raise RuntimeError('Completed sidecar artifact changed size')
 artifacts[rp.name]={'bytes':rp.stat().st_size,'sha256':digest(rp)}
 return artifacts

def completed_receipts():
 primary_roots={'chr1-technical-pilot','remaining-align','workflow-v1','gene-coverage-v1','mitochondria-recovery-v2'}
 for p in sorted((root/'outputs').rglob('*receipt.json')):
  rel=p.relative_to(root/'outputs')
  if any(x.lower() in {'resources','databases','references','cache','raw','inputs'} for x in rel.parts):continue
  if len(rel.parts)<2 or not (rel.parts[0] in primary_roots or re.fullmatch(r'technical-xy-qc-mq20-[0-9]{8}T[0-9]{6}Z',rel.parts[0])):continue
  yield p,None
 for name,file in [('consequence-tables-v1','manifest.json'),('tissue-protein-context-v1','manifest.json'),('cohort-sample-views-v1','manifest.json'),('population-reference-context-v1','workflow-receipt.json'),('combined-annotation-master-v1','workflow-receipt.json'),('fully-annotated-vcf-v1','workflow-receipt.json'),('mitochondrial-reference-context-v1','workflow-receipt.json')]:
  p=root/'outputs'/name/file
  if p.exists():yield p,'completed_sidecar_v1'

for rp,sidecar_kind in completed_receipts():
 if sidecar_kind is None and rp.name not in ('receipt.json','workflow-receipt.json'):continue
 if any(p.startswith('.') or 'partial' in p.lower() for p in rp.relative_to(root).parts):continue
 if rp.is_symlink() or root not in rp.resolve().parents:raise RuntimeError('Unsafe stage receipt')
 parent=rp.parent
 while parent!=root:
  if parent.is_symlink():raise RuntimeError('Symlink receipt parent refused')
  parent=parent.parent
 if rp.stat().st_size>32*1024**2:raise RuntimeError('Receipt/manifest exceeds bounded metadata size')
 raw=rp.read_bytes();d=json.loads(raw)
 if sidecar_kind:
  artifacts=sidecar_artifacts(rp,d)
 elif rp.name=='workflow-receipt.json':
  if d.get('status')!='PASS' or not d.get('identity'):continue
  artifacts=d.get('artifacts',{})
 else:
  if not d.get('completed_utc') or not d.get('fingerprint'):continue
  artifacts=d.get('outputs',{})
 if not artifacts:continue
 receipts.append({'RelativePath':rp.relative_to(root).as_posix(),'SizeBytes':len(raw),'Sha256':hashlib.sha256(raw).hexdigest(),'ContentBase64':base64.b64encode(raw).decode()})
 for name,value in artifacts.items():
  p=safe_artifact(rp.parent,name)
  if any(x.lower() in {'resources','databases','references','cache','raw','inputs'} for x in p.relative_to(root).parts):continue
  # Needed biological outputs, indices, QC, tables and reports; large public
  # resource caches and FASTQ/FASTA are deliberately absent.
  if not p.name.endswith(('.bam','.bai','.vcf','.vcf.gz','.bcf','.tbi','.csi','.tsv','.csv','.tsv.gz','.csv.gz','.txt','.html','.pdf','.xlsx','.docx','.md','.json','.log','.png','.jpg','.jpeg','.svg')):continue
  sha=value['sha256'] if isinstance(value,dict) else value
  if not re.fullmatch('[0-9a-f]{64}',sha):raise RuntimeError('Invalid receipt digest')
  size=p.stat().st_size
  if isinstance(value,dict) and size!=value.get('bytes',value.get('size_bytes',size)):raise RuntimeError('Completed artifact changed size')
  source=p.relative_to(root).as_posix()
  if source in files and files[source]['Sha256']!=sha:raise RuntimeError('Conflicting receipts')
  filename='__'.join(p.relative_to(root).parts[1:])
  if len(filename.encode())>220:filename=hashlib.sha256(source.encode()).hexdigest()[:16]+'__'+p.name
  if filename in names and names[filename]!=source:raise RuntimeError('Native download filename collision')
  names[filename]=source;target=links/filename
  if target.is_symlink():raise RuntimeError('Native download target must be a hardlink, not symlink')
  if target.exists():
   if not os.path.samefile(p,target):raise RuntimeError('Refusing differing native download link')
  elif not compare_only:
   try:os.link(p,target)
   except FileExistsError:
    if not os.path.samefile(p,target):raise
  jrel=target.relative_to(jroot).as_posix()
  files[source]={'RelativePath':source,'ServerPath':str(p),'SizeBytes':size,'Sha256':sha,'SourceReceipt':rp.relative_to(root).as_posix(),'DownloadFileName':filename,'JupyterRelativePath':jrel}
  expose(p,files[source])
excluded_gvcfs=sum('.g.vcf' in pathlib.PurePosixPath(x['RelativePath']).name.lower() for x in files.values())
def local_selected(entry):
 name=pathlib.PurePosixPath(entry['RelativePath']).name.lower()
 parts=pathlib.PurePosixPath(entry['RelativePath']).parts
 if len(parts)>1 and parts[1].startswith('mitochondria') and parts[1]!='mitochondrial-reference-context-v1':
  if parts[1]!='mitochondria-recovery-v2' or len(parts)!=6 or parts[-2]!='calls':return False
  return name in {'mitochondria.filtered.vcf.gz','mitochondria.filtered.vcf.gz.tbi','mitochondria.pass.vcf.gz','mitochondria.pass.vcf.gz.tbi','mitochondria.raw.vcf.gz','mitochondria.raw.vcf.gz.tbi','research-report.json','numt-assessment.json','contamination.json','heteroplasmy.tsv','circular-coverage.tsv'}
 if '.g.vcf' in name:return False
 if '/checkpoints/' in entry['RelativePath'] and ('.vcf' in name or name.endswith(('.bcf','.tbi','.csi'))):return False
 # One final annotation view is the local deliverable. Intermediate matrices and
 # duplicate context/FVCF copies remain on the server; no existing local file is deleted.
 if len(parts)>1 and parts[1] in {'consequence-tables-v1','tissue-protein-context-v1','population-reference-context-v1','combined-annotation-master-v1','fully-annotated-vcf-v1'}:
  return name in {'manifest.json','workflow-receipt.json','receipt.json','provenance.json'} and len(parts)==3
 if len(parts)>3 and parts[:3]==('outputs','workflow-v1','stages') and parts[3] in {'06_annotate','07_report'}:
  return name.endswith(('.json','.log','.md','.txt','.png','.jpg','.jpeg','.svg'))
 return True
selected=[x for x in files.values() if local_selected(x)]
excluded_redundant=sum(not local_selected(x) and '.g.vcf' not in pathlib.PurePosixPath(x['RelativePath']).name.lower() for x in files.values())
ordered=sorted(selected,key=lambda x:(not x['RelativePath'].endswith('.bam'),-x['SizeBytes']))
rows=sorted(view_rows.values(),key=lambda x:x['ViewPath'])
bams=sum(x['ViewPath'].startswith('BAM/') and x['ViewPath'].endswith('.bam') for x in rows)
gvcfs=sum(x['ViewPath'].startswith('gVCF/') and len(pathlib.PurePosixPath(x['ViewPath']).parts)==3 and x['ViewPath'].endswith('.g.vcf.gz') for x in rows)
cohort=any('/stages/04_joint/' in x['Source'] and x['ViewPath']=='VCF/cohort.raw.vcf.gz' for x in rows)
filtered=any('/stages/05_filter/' in x['Source'] and x['ViewPath']=='VCF/cohort.quality.vcf.gz' for x in rows)
status='РћР±С‰РёР№ VCF РіРѕС‚РѕРІ' if cohort else 'РћР±С‰РµРіРѕ VCF РїРѕРєР° РЅРµС‚: СЃС‚Р°РґРёСЏ РѕР±СЉРµРґРёРЅРµРЅРёСЏ РµС‰С‘ РЅРµ Р·Р°РІРµСЂС€РµРЅР°'
body=['<!doctype html><html lang="ru"><meta charset="utf-8"><meta http-equiv="refresh" content="60"><title>Р“РѕС‚РѕРІС‹Рµ СЂРµР·СѓР»СЊС‚Р°С‚С‹</title><style>body{font:16px system-ui;max-width:1300px;margin:32px auto;padding:0 20px}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}table{border-collapse:collapse;width:100%}small{color:#555}</style>','<h1>Р“РѕС‚РѕРІС‹Рµ СЂРµР·СѓР»СЊС‚Р°С‚С‹ РЅР° СЃРµСЂРІРµСЂРµ</h1>',f'<p>BAM: {bams}; gVCF РѕР±СЂР°Р·С†РѕРІ: {gvcfs} (С‚РµРєСѓС‰РёР№ scope вЂ” chr1вЂ“22/autosomes). '+status+'.</p>','<p>РџРѕРєР°Р·Р°РЅС‹ С„Р°Р№Р»С‹ С‚РѕР»СЊРєРѕ РёР· Р·Р°РІРµСЂС€С‘РЅРЅС‹С… СЃС‚Р°РґРёР№ СЃ РєРІРёС‚Р°РЅС†РёСЏРјРё. BAM вЂ” РІС‹СЂР°РІРЅРёРІР°РЅРёСЏ; gVCF вЂ” РїРѕРѕР±СЂР°Р·С†РѕРІС‹Рµ РІС‹Р·РѕРІС‹; VCF вЂ” РѕР±С‰РёР№ РЅР°Р±РѕСЂ. РўРµС…РЅРёС‡РµСЃРєРёРµ СЂРµР·СѓР»СЊС‚Р°С‚С‹ РµС‰С‘ С‚СЂРµР±СѓСЋС‚ Р°РЅРЅРѕС‚Р°С†РёРё Рё РёРЅС‚РµСЂРїСЂРµС‚Р°С†РёРё.</p>','<p><a href="BAM/">BAM</a> В· <a href="gVCF/">gVCF</a> В· <a href="VCF/">VCF</a> В· <a href="QC/">QC</a> В· <a href="Reports/">Reports</a></p>','<table><tr><th>Р¤Р°Р№Р»</th><th>GiB</th><th>РћР±Р»Р°СЃС‚СЊ</th></tr>']
for row in rows:body.append('<tr><td><a href="../marina-verified-downloads/'+quote(row['DownloadFileName'],safe='')+'?download=1">'+html.escape(row['ViewPath'])+'</a></td><td>'+f'{row["SizeBytes"]/2**30:.4f}'+'</td><td>'+html.escape(row['Scope'])+'</td></tr>')
body.append('</table></html>')
if not compare_only:
 temporary=view/('INDEX.html.'+str(os.getpid())+'.tmp');temporary.write_text('\n'.join(body),encoding='utf-8');os.replace(temporary,view/'INDEX.html')
view_manifest={'GeneratedUtc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'CompletedBAMs':bams,'AutosomalGVCFs':gvcfs,'CohortVCFReady':cohort,'FilteredCohortVCFReady':filtered,'Files':rows}
if not compare_only:
 temporary=view/('manifest.json.'+str(os.getpid())+'.tmp');temporary.write_text(json.dumps(view_manifest,ensure_ascii=False,indent=2),encoding='utf-8');os.replace(temporary,view/'manifest.json')
print(json.dumps({'SchemaVersion':1,'PlanGeneratorVersion':3,'Mode':'COMPARE_ONLY_NO_LINKS_WRITTEN' if compare_only else 'PUBLISH_HARDLINK_VIEW','GeneratedUtc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'ReceiptPolicy':'Completed known-root receipts and strict source-bound sidecars; expected producer SHA256 and actual source size checked; no new BAM/table byte hashes; local backup verification mode recorded separately','LocalSelection':'Canonical VCF-first cohort/sample VCFs, masters/top/types/applied evidence, BAM/index, current MT/QC and provenance; redundant intermediate matrices and duplicate enriched VCF copies stay on server','ExcludedIndividualGVCFArtifacts':excluded_gvcfs,'ExcludedRedundantArtifacts':excluded_redundant,'Files':ordered,'Receipts':receipts,'ServerView':{'Path':str(view),'CompletedBAMs':bams,'AutosomalGVCFs':gvcfs,'CohortVCFReady':cohort,'FilteredCohortVCFReady':filtered}}))
'''


def local_capacity(plan, destination, staging, reserve=15*1024**3):
    """Conservative remaining bytes; only source identity/size metadata is read."""
    destination,staging=Path(destination).resolve(strict=True),Path(staging).resolve(strict=True)
    records_path=destination/'size-only-files.json'
    records=json.loads(records_path.read_text(encoding='utf-8-sig')) if records_path.is_file() else {}
    pending=0;selected=0;canonical=0
    for entry in plan['Files']:
        relative=entry['RelativePath'];name=PurePosixPath(relative).name.lower()
        if name.endswith(('.bam','.bai')):continue
        selected+=1;canonical+=relative.startswith('outputs/cohort-sample-views-v1/')
        identity=json.dumps([relative,entry['SizeBytes'],entry['Sha256']],separators=(',',':'))
        record=records.get('source-'+hashlib.sha256(identity.encode()).hexdigest(),{})
        complete=False
        if record.get('ExpectedServerSha256')==entry['Sha256']:
            stored=(destination/record.get('StoragePath','')).resolve();friendly=(destination/record.get('LocalPath','')).resolve()
            if destination in stored.parents and destination in friendly.parents and stored.is_file() and friendly.is_file() and not stored.is_symlink() and not friendly.is_symlink():
                complete=stored.stat().st_size==entry['SizeBytes'] and stored.stat().st_mtime_ns==record.get('MtimeNs') and os.path.samefile(stored,friendly)
        if complete:continue
        flat=staging/entry['DownloadFileName'];partial=flat.with_name(flat.name+'.part');offset=0
        if not flat.is_symlink() and flat.is_file() and flat.stat().st_size==entry['SizeBytes']:offset=entry['SizeBytes']
        elif not partial.is_symlink() and partial.is_file():offset=min(entry['SizeBytes'],partial.stat().st_size)
        pending+=entry['SizeBytes']-offset
    free=shutil.disk_usage(destination).free
    return dict(status='PASS' if free>=pending+reserve else 'BLOCKED_INSUFFICIENT_SPACE_PLAN_NOT_PUBLISHED',
        required_remaining_bytes=pending,free_bytes=free,reserve_bytes=reserve,selected_small_result_files=selected,
        completed_canonical_files=canonical,all_required_outputs_retained=True,patient_file_hashes_computed=False,
        partial_offsets_retained=True,existing_files_deleted=False)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', required=True)
    parser.add_argument('--ssh-config', required=True)
    parser.add_argument('--known-hosts', required=True)
    parser.add_argument('--remote-root', required=True)
    parser.add_argument('--workspace-root', required=True)
    parser.add_argument('--jupyter-root', required=True)
    parser.add_argument('--output', required=True, help='Private ignored expected-identity plan')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--compare', action='store_true', help='One read-only SSH inventory; no server hardlinks/view writes; new private output only')
    parser.add_argument('--compare-with-plan', type=Path, help='Existing private plan; compare metadata only, never replace it')
    parser.add_argument('--watch', action='store_true')
    parser.add_argument('--interval', type=int, default=60)
    parser.add_argument('--max-hours', type=float, default=8)
    parser.add_argument('--deadline-utc', help='Absolute stop time for local plan refresh')
    parser.add_argument('--pause-flag', help='Persistent local user-pause marker to respect')
    parser.add_argument('--local-results-root',type=Path,help='Private result destination for metadata-only required-byte capacity gate')
    parser.add_argument('--local-results-staging',type=Path,help='Existing independent small-result staging folder; no credentials are opened')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.host):
        raise ValueError('Confirmed SSH alias required')
    if args.interval < 60 or not 0 < args.max_hours <= 8:
        raise ValueError('Interval >=60 seconds and max-hours <=8 required')
    for value in (args.remote_root, args.workspace_root, args.jupyter_root):
        p = PurePosixPath(value)
        if not p.is_absolute() or '..' in p.parts or len(p.parts) < 3 or any(ord(c) < 32 for c in value):
            raise ValueError('Dedicated absolute remote directories required')
    if args.compare and (args.watch or args.execute):
        raise ValueError('--compare is one-shot and cannot be combined with --watch/--execute')
    if args.compare and Path(args.output).exists():
        raise ValueError('Compare requires a new private output; active plan is never replaced')
    if not args.execute and not args.compare:
        print('DRY RUN: no SSH, hardlinks or plan files created')
        return
    config, known = Path(args.ssh_config).expanduser().resolve(strict=True), Path(args.known_hosts).expanduser().resolve(strict=True)
    ssh = ['ssh', '-T', '-F', str(config), '-o', 'BatchMode=yes', '-o', 'StrictHostKeyChecking=yes',
           '-o', 'UserKnownHostsFile=' + str(known), '-o', 'ConnectTimeout=30',
           '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=4']
    command = shlex.join(['python3', '-', json.dumps({'root': args.remote_root, 'workspace': args.workspace_root,
                                                    'jupyter': args.jupyter_root, 'compare_only': args.compare})])
    deadline=time.monotonic()+args.max_hours*3600
    wall_deadline=time.time()+args.max_hours*3600
    import datetime as dt
    budget_deadline=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00').timestamp()
    wall_deadline=min(wall_deadline,budget_deadline)
    deadline=min(deadline,time.monotonic()+max(0,budget_deadline-time.time()))
    if args.deadline_utc:
        import datetime as dt
        absolute=dt.datetime.fromisoformat(args.deadline_utc.replace('Z','+00:00'))
        if absolute.tzinfo is None:raise ValueError('Deadline must contain a UTC offset')
        deadline=min(deadline,time.monotonic()+max(0,absolute.timestamp()-time.time()))
        wall_deadline=min(wall_deadline,absolute.timestamp())
    while True:
        if time.monotonic()>=deadline or time.time()>=wall_deadline or (args.pause_flag and Path(args.pause_flag).exists()):
            print('PAUSED local plan refresh; server outputs unchanged',flush=True);break
        try:
            result = subprocess.run(ssh + [args.host, command], input=HELPER, text=True, capture_output=True,
                                    timeout=min(120,max(1,deadline-time.monotonic())))
        except subprocess.TimeoutExpired:
            if time.monotonic()>=deadline:break
            raise
        if result.returncode:
            raise RuntimeError(result.stderr)
        plan = json.loads(result.stdout)
        if args.compare_with_plan:
            previous=json.loads(args.compare_with_plan.read_text(encoding='utf-8-sig'))
            old={x['RelativePath']:(x['SizeBytes'],x['Sha256']) for x in previous['Files']}
            new={x['RelativePath']:(x['SizeBytes'],x['Sha256']) for x in plan['Files']}
            plan['Comparison']={'Added':len(set(new)-set(old)),'Removed':len(set(old)-set(new)),
                                'ExpectedIdentityChanged':sum(old[p]!=new[p] for p in set(old)&set(new)),
                                'PreviousPlanChanged':False}
        for receipt in plan['Receipts']:
            data = base64.b64decode(receipt['ContentBase64'], validate=True)
            if len(data) != receipt['SizeBytes'] or hashlib.sha256(data).hexdigest() != receipt['Sha256']:
                raise ValueError('Receipt metadata integrity check failed')
        for entry in plan['Files']:
            rel = quote(entry['JupyterRelativePath'], safe='/')
            entry['JupyterDownloadPath'] = '/files/' + rel + '?download=1'
            entry['JupyterUIPath'] = '/lab/tree/' + rel
        canonical=any(e['RelativePath'].startswith('outputs/cohort-sample-views-v1/') for e in plan['Files'])
        if canonical and args.execute and not (args.local_results_root and args.local_results_staging):
            raise ValueError('Final canonical results require explicit local destination/staging capacity inventory before publication')
        if args.local_results_root or args.local_results_staging:
            if not (args.local_results_root and args.local_results_staging):raise ValueError('Both local capacity directories are required')
            plan['LocalCapacity']=local_capacity(plan,args.local_results_root,args.local_results_staging)
            if plan['LocalCapacity']['status']!='PASS' and args.execute:
                blocked=Path(args.output).with_suffix('.capacity-blocked.json')
                blocked.write_text(json.dumps(plan['LocalCapacity'],indent=2)+'\n',encoding='utf-8')
                print('BLOCKED local result capacity: required remaining bytes '+str(plan['LocalCapacity']['required_remaining_bytes'])+'; original plan and all required files retained',flush=True)
                if not args.watch:raise ValueError('Insufficient local space above the unchanged 15GiB reserve; plan not published')
                time.sleep(min(args.interval,max(0,deadline-time.monotonic())));continue
        output = Path(args.output).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        temp = output.with_name(output.name + '.tmp')
        temp.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        temp.replace(output)
        print(f'PASS native plan: {len(plan["Files"])} completed artifacts; {len(plan["Receipts"])} receipts; view {plan["ServerView"]}',flush=True)
        if not args.watch or time.monotonic()+args.interval>=deadline:break
        time.sleep(min(args.interval,max(0,deadline-time.monotonic()),max(0,wall_deadline-time.time())))


if __name__ == '__main__':
    main()
