#!/usr/bin/env python3
"""Pure exact-allele/gene SpliceAI field projection; no transcript reassignment."""
import json
from urllib.parse import quote

SPLICE_FIELDS=['spliceai_status','spliceai_gene','spliceai_DS_AG','spliceai_DS_AL','spliceai_DS_DG','spliceai_DS_DL','spliceai_DP_AG','spliceai_DP_AL','spliceai_DP_DG','spliceai_DP_DL','spliceai_prediction_json','spliceai_source_json','spliceai_transcript_scope']
SCORE_FIELDS=['DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG','DP_DL']
SCOPE='GENCODE24_GENE_UNION_NOT_ENST49_MAPPING'
def key(chrom,pos,ref,alt):return 'chr'+str(chrom).removeprefix('chr')+':'+str(int(pos))+':'+ref+':'+alt
def enrich_row(row,predictions,source_json):
 result={field:'' for field in SPLICE_FIELDS};result.update(spliceai_source_json=source_json if predictions else '',spliceai_transcript_scope=SCOPE)
 expected=key(row['chrom'],row['pos'],row['ref'],row['alt'])
 if any(p['variant_key']!=expected or p['chrom'].removeprefix('chr')!=str(row['chrom']).removeprefix('chr') or int(p['pos'])!=int(row['pos']) or p['ref']!=row['ref'] or p['alt']!=row['alt'] for p in predictions):raise ValueError('SpliceAI exact allele projection mismatch')
 result['spliceai_prediction_json']=json.dumps(predictions,sort_keys=True,separators=(',',':'))
 matching=[p for p in predictions if p['gene']==row.get('gene')]
 if not predictions:result['spliceai_status']='OUTSIDE_APPROVED_PUBLIC39_MODEL_GENE_SPANS_NOT_ZERO'
 elif not matching:result['spliceai_status']='GENCODE24_MODEL_GENE_NOT_MATCHED_TO_VEP_GENE_NO_TRANSCRIPT_ASSIGNMENT'
 elif len(matching)!=1:result['spliceai_status']='AMBIGUOUS_GENE_MODEL_RESULTS_NO_SCALAR_ASSIGNMENT'
 else:
  p=matching[0];result.update(spliceai_status=p['status'],spliceai_gene=p['gene'],spliceai_source_json=p['source_json'])
  if p['status']=='SCORED_EXACT_GRCH38_ALLELE_GENE':
   for field in SCORE_FIELDS:
    if p[field] in ['',None,'.']:raise ValueError('Scored status with missing value')
    result['spliceai_'+field]=str(p[field])
 return result
def readable(packet,fields):
 raw=set();statuses=set();source=set()
 for item in packet['rows']:
  row=dict(zip(fields,item['values']));statuses.add('|'.join([row.get('sample',''),row.get('gene',''),row.get('spliceai_status',''),SCOPE]))
  for p in json.loads(row.get('spliceai_prediction_json') or '[]'):
   if p['status']=='SCORED_EXACT_GRCH38_ALLELE_GENE':raw.add(p['raw_annotation'])
  if row.get('spliceai_source_json'):source.add(row['spliceai_source_json'])
 return {'SpliceAI':','.join(sorted(raw)) or '.', 'WGS_SPLICEAI_STATUS':','.join(quote(v,safe='._:-') for v in sorted(statuses)) or '.', 'WGS_SPLICEAI_SOURCE':','.join(quote(v,safe='._:-') for v in sorted(source)) or '.'}
