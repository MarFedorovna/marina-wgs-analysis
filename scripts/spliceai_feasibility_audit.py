#!/usr/bin/env python3
"""Read-only local source/runtime/scope audit; no patient inference or network."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

def main(root,output):
 import pysam
 runtime=root/'resources/spliceai-research-runtime-v1';package=runtime/'venv/lib/python3.10/site-packages/spliceai';panel=json.loads((root/'resources/public-hcm-cds-v2/panel-manifest.json').read_text());genes=panel['genes'];symbols={g.get('symbol',g.get('gene_symbol')) for g in genes};symbols.discard(None)
 annotation=package/'annotations/grch38.txt'
 with annotation.open() as f:annotations=[r for r in csv.DictReader(f,delimiter='\t') if r['#NAME'] in symbols]
 primary=root/'outputs/workflow-v1/stages/05_filter';receipt=json.loads((primary/'workflow-receipt.json').read_text())
 if receipt.get('status')!='PASS':raise ValueError('Accepted source must PASS')
 source=primary/'cohort.filtered.vcf.gz';entry=receipt['artifacts'][source.name]
 if source.stat().st_size!=entry['bytes']:raise ValueError('Source size changed')
 counts={};seen=set();types={'SNV':0,'indel':0,'complex_or_unsupported':0};start=time.monotonic()
 with pysam.VariantFile(str(source)) as vcf:
  for r in annotations:
   chrom='chr'+r['CHROM'].removeprefix('chr');number=0
   if chrom not in vcf.header.contigs:continue
   for record in vcf.fetch(chrom,int(r['TX_START']),int(r['TX_END'])):
    k=(record.chrom,record.pos,record.ref,record.alts)
    if k in seen:continue
    seen.add(k);number+=1
    alt=record.alts[0] if record.alts and len(record.alts)==1 else ''
    kind='SNV' if len(record.ref)==len(alt)==1 and alt in 'ACGT' else 'indel' if alt and set(record.ref+alt)<=set('ACGT') and min(len(record.ref),len(alt))==1 else 'complex_or_unsupported';types[kind]+=1
   counts[r['#NAME']]=counts.get(r['#NAME'],0)+number
 scan_seconds=time.monotonic()-start;pages=0;hits=0;splice=0;fields=set()
 for folder in [root/'resources/myvariant-public-hcm-v1',root/'resources/myvariant-public-hcm-tnnt2-v1']:
  for path in folder.glob('*/page-*.json'):
   pages+=1
   for hit in json.loads(path.read_text()).get('hits',[]):
    hits+=1;value=hit.get('dbnsfp',{}).get('spliceai')
    if value is not None:
     splice+=1
     if isinstance(value,dict):fields.update(value)
 result=dict(status='READONLY_SPLICEAI_FEASIBILITY_AUDIT',patient_inference=False,network_requests=0,source_vcf=dict(path=str(source),bytes=entry['bytes'],sha256=entry['sha256'],hash_from_accepted_receipt_not_rehashed=True),approved_panel_genes=len(genes),panel_symbols=len(symbols),bundled_annotation_genes=len(annotations),unmatched_symbols=sorted(symbols-{r['#NAME'] for r in annotations}),unique_variants_in_bundled_public_gene_spans=len(seen),variant_types=types,gene_counts=counts,variant_count_seconds=scan_seconds,prefetch=dict(pages=pages,hits=hits,records_with_spliceai=splice,spliceai_fields=sorted(fields)),annotation=dict(path=str(annotation),sha256=hashlib.sha256(annotation.read_bytes()).hexdigest(),scope='official bundled GRCh38 GENCODE24 gene-union annotations, NOT VEP GENCODE49 transcript mapping'),models=[dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted((package/'models').glob('*.h5'))])
 with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path('/srv/marina'));p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.root,a.output)
