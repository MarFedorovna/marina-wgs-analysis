#!/usr/bin/env python3
"""Public ClinGen HCM genes + public GENCODE49 CDS only; accepts no patient data."""
import argparse
import collections
import gzip
import hashlib
import json
from pathlib import Path
import re
from prefetch_myvariant_public_hcm import public_genes,load

CODE=Path(__file__).resolve().parents[1]


def union(ranges):
    result=[]
    for start,end in sorted(ranges):
        if result and start<=result[-1][1]:result[-1][1]=max(end,result[-1][1])
        else:result.append([start,end])
    return result


def build(cfg,gtf,target):
    assert not target.exists(),'Public panel version already exists'
    sha=hashlib.sha256()
    with gtf.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):sha.update(block)
    assert sha.hexdigest()=='d6e6fe0515c95b2a8cd36a853c1989cee9115c736c60237c56ae92b9daaaf7c4','GENCODE49 public source differs'
    genes=public_genes(cfg);selected={g['hgnc_id']:g for g in genes};gene_ids={};ranges=collections.defaultdict(list)
    with gzip.open(gtf,'rt') as stream:
        for line in stream:
            if line.startswith('#'):continue
            row=line.rstrip('\n').split('\t')
            if row[2] not in ('gene','CDS'):continue
            attrs=dict(re.findall(r'(\w+) "([^"]*)";',row[8]));gid=attrs.get('gene_id')
            if row[2]=='gene' and attrs.get('hgnc_id') in selected:
                gene_ids[gid]=attrs['hgnc_id']
            if row[2]=='CDS' and gid in gene_ids and re.fullmatch(r'chr(?:[1-9]|1\d|2[0-2]|X|Y)',row[0]):
                ranges[(gene_ids[gid],row[0])].append((int(row[3])-1,int(row[4])))
    target.mkdir(mode=0o700,parents=True)
    by_chrom=collections.defaultdict(list);rows=[]
    for (hgnc,chrom),values in sorted(ranges.items()):
        for start,end in union(values):rows.append((chrom,start,end,hgnc,selected[hgnc]['symbol']));by_chrom[chrom].append((start,end))
    (target/'hcm-gene-cds.bed').write_text(''.join('\t'.join(map(str,r))+'\n' for r in rows))
    merged=[]
    for chrom,values in sorted(by_chrom.items()):
        for start,end in union(values):merged.append((chrom.removeprefix('chr'),start,end))
    bed=target/'hcm-cds-union-cadd.bed';bed.write_text(''.join('\t'.join(map(str,r))+'\n' for r in merged))
    delta_rows=[r for r in rows if r[3]=='HGNC:11949'] if cfg.get('additional_released_hcm_cspec') else []
    delta_bed=target/'hcm-cds-delta-cadd.bed'
    delta_bed.write_text(''.join('%s\t%d\t%d\n'%(chrom.removeprefix('chr'),start,end) for chrom,start,end,_,_ in delta_rows))
    present={key[0] for key in ranges}
    result=dict(status='PUBLIC_REFERENCE_PANEL_ONLY',source_clingen_sha256=cfg['clingen_csv_sha256'],
        source_gencode_sha256=sha.hexdigest(),gencode_release=49,assembly='GRCh38',
        disease_id=cfg['disease_id'],curation_rows=40,genes=genes,total_genes=len(genes),genes_with_cds=len(present),
        genes_without_supported_cds=[g for g in genes if g['hgnc_id'] not in present],
        cds_union_intervals=len(merged),cds_union_bases=sum(r[2]-r[1] for r in merged),
        bed_sha256=hashlib.sha256(bed.read_bytes()).hexdigest(),coordinate_system='0-based half-open BED',
        source_cspec_addition=cfg.get('additional_released_hcm_cspec'),delta_gene='TNNT2' if delta_rows else None,
        delta_bed_sha256=hashlib.sha256(delta_bed.read_bytes()).hexdigest(),
        source_priority_note='CSV exact-HCM validity and released CSpec HCM applicability remain separate; no DCM/ARVC classification is promoted to HCM',
        no_patient_selection=True,no_patient_variants_read=True)
    (target/'panel-manifest.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=CODE/'config/myvariant-public-hcm-v1.json');p.add_argument('--gtf',type=Path,default=Path('/srv/marina/resources/gene-coverage/gencode49/gencode.v49.annotation.gtf.gz'));p.add_argument('--output',type=Path,default=Path('/srv/marina/resources/public-hcm-cds-v1'));a=p.parse_args();print(json.dumps(build(load(a.config),a.gtf,a.output),indent=2))
