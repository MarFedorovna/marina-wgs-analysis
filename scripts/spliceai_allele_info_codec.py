#!/usr/bin/env python3
"""Portable allele-level SpliceAI evidence; original core row packets unchanged."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re

def pinned(name,digest):
 path=Path(__file__).resolve().parent/(name+'.py')
 if path.is_symlink() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:raise ValueError('Preserved codec/projection source changed')
 spec=importlib.util.spec_from_file_location('allele_'+name,path);value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value
base=pinned('embed_vcf_evidence','18e7a2102aceb5b962892242763d7d331523d7669c93f16b12f9fe71a1eb238e')
projection=pinned('spliceai_packet_context','b6e0b87feaf763cdfc9388f6fd49e09b3beccc983adbaf9355356342a2d90d48')
CODEC=base.CODEC
CONTRACT='ALLELE_INFO_V1_ROW_PROJECTION_ON_READ'
SCOPE=projection.SCOPE
SPLICE_FIELDS=projection.SPLICE_FIELDS
SCORE_FIELDS=projection.SCORE_FIELDS
OUTSIDE='OUTSIDE_APPROVED_PUBLIC39_MODEL_GENE_SPANS_NOT_ZERO'
OVERLAY_INFO_IDS=['SpliceAI','WGS_SPLICEAI_STATUS','WGS_SPLICEAI_EVIDENCE']
INFO_IDS=[*base.INFO_IDS,*OVERLAY_INFO_IDS]
NAMESPACED_INFO_IDS=INFO_IDS
SAMPLE_PACKET_KEYS={'v','allele','rows','evidence','gene_context_keys','scope','sample'}
for name in ['encode_packet','decode_packet','encode_schema','decode_schema','schema_header','decode_header_bundle','chunk_headers','public_bundle']:globals()[name]=getattr(base,name)
enrich_row=projection.enrich_row
def source(value):
 data=json.loads(value)
 if data.get('version')!='1.3.1' or data.get('distance')!=500 or data.get('mask')!=0 or data.get('annotation_scope')!=SCOPE or data.get('assembly')!='GRCh38' or len(data.get('model_sha256',[]))!=5:raise ValueError('Unsupported SpliceAI source model/parameters')
 for digest in [*data['model_sha256'],data.get('annotation_sha256',''),data.get('reference_sha256','')]:
  if not re.fullmatch('[a-f0-9]{64}',digest):raise ValueError('Full model/annotation/reference hashes required')
 return data
def validate_predictions(predictions,allele,source_json):
 if not predictions:return
 metadata=source(source_json);seen=set();expected=projection.key(*allele)
 for row in predictions:
  if row['variant_key']!=expected or [str(row['chrom']),int(row['pos']),row['ref'],row['alt']]!=allele or not row['gene'] or row['gene'] in seen or row['transcript_scope']!=SCOPE or source(row['source_json'])!=metadata:raise ValueError('SpliceAI allele/gene/model identity mismatch')
  seen.add(row['gene'])
  if row['status']=='SCORED_EXACT_GRCH38_ALLELE_GENE':
   if any(not math.isfinite(float(row[f])) or not 0<=float(row[f])<=1 for f in SCORE_FIELDS[:4]) or any(not re.fullmatch('-?[0-9]+',row[f]) or abs(int(row[f]))>500 for f in SCORE_FIELDS[4:]):raise ValueError('SpliceAI scores/offsets invalid')
   if row['raw_annotation']!='|'.join([row['alt'],row['gene'],*[row[f] for f in SCORE_FIELDS]]):raise ValueError('Raw standard score vector differs')
  elif any(row.get(f) not in ['',None,'.'] for f in SCORE_FIELDS):raise ValueError('Unscored status must have missing scores')
def overlay_info(allele,predictions,source_json):
 if not predictions:return dict(SpliceAI='.',WGS_SPLICEAI_STATUS=OUTSIDE,WGS_SPLICEAI_EVIDENCE='.')
 validate_predictions(predictions,allele,source_json);scored=[p for p in predictions if p['status']=='SCORED_EXACT_GRCH38_ALLELE_GENE'];status='SCORED' if scored else 'NOT_SCORED_EXPLICIT_FAILURE';packet=dict(v=1,allele=allele,status=status,scope=SCOPE,predictions=predictions,source_json=source_json)
 return dict(SpliceAI=','.join(sorted(p['raw_annotation'] for p in scored)) or '.',WGS_SPLICEAI_STATUS=status,WGS_SPLICEAI_EVIDENCE=encode_packet(packet))
def projection_from_info(info,allele=None):
 status=info.get('WGS_SPLICEAI_STATUS');token=info.get('WGS_SPLICEAI_EVIDENCE');standard=info.get('SpliceAI')
 if status==OUTSIDE:
  if token not in [None,'','.'] or standard not in [None,'','.']:raise ValueError('Outside-scope status contains invented scores')
  return [],''
 if status not in ['SCORED','NOT_SCORED_EXPLICIT_FAILURE'] or not token or token=='.':raise ValueError('Actual encoded SpliceAI evidence or explicit outside status required')
 packet=decode_packet(token)
 if packet.get('status')!=status or packet.get('scope')!=SCOPE or (allele is not None and packet.get('allele')!=allele):raise ValueError('Allele SpliceAI packet scope/status mismatch')
 expected=overlay_info(packet['allele'],packet['predictions'],packet['source_json'])
 if any(info.get(k)!=v for k,v in expected.items()):raise ValueError('Standard SpliceAI INFO and full evidence disagree')
 return packet['predictions'],packet['source_json']
def select_sample(packet,fields,sample):return {key:value for key,value in base.select_sample(packet,fields,sample).items() if key in SAMPLE_PACKET_KEYS}
def rebuild_readable(packet,fields):
 predictions={};sources=set()
 for item in packet['rows']:
  row=dict(zip(fields,item['values']))
  for p in json.loads(row.get('spliceai_prediction_json') or '[]'):predictions[(p['variant_key'],p['gene'])]=p
  if row.get('spliceai_source_json'):sources.add(row['spliceai_source_json'])
 if len(sources)>1:raise ValueError('Contradictory row model source metadata')
 return dict(base.rebuild_readable(packet,fields),**overlay_info(packet['allele'],list(predictions.values()),next(iter(sources)) if sources else ''))
def overlay_headers(schema):
 return [b'##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI1.3.1 D500 mask0; GENCODE24 gene union NOT VEP GENCODE49 transcript mapping; ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">\n',b'##INFO=<ID=WGS_SPLICEAI_STATUS,Number=1,Type=String,Description="SCORED, NOT_SCORED_EXPLICIT_FAILURE, or OUTSIDE_APPROVED_PUBLIC39_MODEL_GENE_SPANS_NOT_ZERO; never infer missing scores as zero">\n',b'##INFO=<ID=WGS_SPLICEAI_EVIDENCE,Number=1,Type=String,Description="Full exact-allele gene-model predictions and source hashes, JSON UTF8 zlib base64url without padding; WGS_SPLICEAI_SCHEMA declares row projection">\n',f'##WGS_SPLICEAI_SCHEMA=<Version=1,Codec="{CODEC}",Projection="{CONTRACT}",Data="{encode_schema(schema)}">\n'.encode()]
def headers(schema,coverage_token,public_token):return base.headers(schema,coverage_token,public_token)
def decode_projection_schema(line):
 if isinstance(line,bytes):line=line.decode()
 match=re.fullmatch(r'##WGS_SPLICEAI_SCHEMA=<Version=1,Codec="'+re.escape(CODEC)+r'",Projection="'+CONTRACT+r'",Data="([A-Za-z0-9_-]+)">\s*',line)
 if not match:raise ValueError('Unsupported SpliceAI projection header')
 value=decode_schema(match[1])
 if value.get('projection_contract')!=CONTRACT or value.get('projection_fields')!=SPLICE_FIELDS:raise ValueError('Projection field schema differs')
 return value
def strip_overlay(line,only_spliceai=False):
 ids=OVERLAY_INFO_IDS if only_spliceai else INFO_IDS;ending=b'\r\n' if line.endswith(b'\r\n') else b'\n' if line.endswith(b'\n') else b'';parts=(line[:-len(ending)] if ending else line).split(b'\t');parts[7]=b';'.join(p for p in parts[7].split(b';') if p.split(b'=',1)[0].decode() not in ids) or b'.';return b'\t'.join(parts)+ending
