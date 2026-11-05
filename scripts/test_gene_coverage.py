#!/usr/bin/env python3
"""Synthetic interval, denominator and resource-subset regression tests (no patients)."""
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import io
from gene_coverage import union, model, write_bed, aggregate, gene_rows, hcm_curations, subset_hcm, download_gtf


class CoverageTests(unittest.TestCase):
    def test_bounded_resource_ranges(self):
        data=b'A'*(2**21)+b'TEST'
        class Response(io.BytesIO):
            status=206
        def response(request,timeout):
            start,end=map(int,request.headers['Range'].split('=')[1].split('-'))
            result=Response(data[start:end+1]);result.headers={'Content-Range':f'bytes {start}-{end}/{len(data)}'};return result
        with tempfile.TemporaryDirectory() as d,patch('gene_coverage.urllib.request.urlopen',side_effect=response):
            p=Path(d)/'resource.gz'
            download_gtf(p,{'url':'https://public.invalid/fixture','bytes':len(data),'publisher_md5':hashlib.md5(data).hexdigest()})
            self.assertEqual(p.read_bytes(),data)

    def test_union_global_thresholds_and_na(self):
        self.assertEqual(union([(10,20),(0,10),(5,7),(30,40)]),[(0,20),(30,40)])
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);gtf=d/'fixture.gtf.gz'
            rows=[]
            for chrom,gid,name,hgnc,typ,start,end in [('chr1','g1.1','A','HGNC:1','protein_coding',1,20),('chr1','g2.1','B','HGNC:2','protein_coding',11,30),('chrX','g3.1','C','HGNC:3','protein_coding',1,20),('chr1','g4.1','D','HGNC:4','lncRNA',41,50)]:
                attr=f'gene_id "{gid}"; gene_name "{name}"; hgnc_id "{hgnc}"; gene_type "{typ}";'
                rows.append(f'{chrom}\ttest\tgene\t{start}\t{end}\t.\t+\t.\t{attr}')
                if typ=='protein_coding':
                    rows.append(f'{chrom}\ttest\tCDS\t{start}\t{end}\t.\t+\t0\t{attr}')
            with gzip.open(gtf,'wt') as f:f.write('\n'.join(rows)+'\n')
            lengths={'chr1':100,'chrX':100};genes,intervals=model(gtf,lengths)
            bedrows=write_bed(d/'model.bed',genes,intervals,lengths)
            self.assertEqual(sum(b-a for c,a,b,n in bedrows if n=='__GLOBAL__'),30)
            self.assertEqual(sum(b-a for c,a,b,n in bedrows if n!='__GLOBAL__'),40)
            prefix=d/'coverage'
            with gzip.open(str(prefix)+'.regions.bed.gz','wt') as means,gzip.open(str(prefix)+'.thresholds.bed.gz','wt') as thresholds:
                thresholds.write('#chrom\tstart\tend\tregion\t1X\t10X\t20X\n')
                for c,a,b,n in bedrows:
                    means.write(f'{c}\t{a}\t{b}\t{n}\t5.25\n')
                    thresholds.write(f'{c}\t{a}\t{b}\t{n}\t{b-a}\t0\t0\n')
            stats=aggregate(prefix,bedrows);self.assertEqual(stats['__GLOBAL__']['CDSBases'],30)
            self.assertEqual(stats['__GLOBAL__']['MeanDepth'],5.25)
            rows=gene_rows(genes,intervals,stats)
            x=next(r for r in rows if r['Chromosome']=='chrX');self.assertEqual(x['CDSBases'],'NA')
            self.assertEqual(x['Status'],'NOT_ASSESSED_CHROMOSOME')
            curation=[{'GENE ID (HGNC)':'HGNC:1','GENE SYMBOL':'OLD_ALIAS','CLASSIFICATION':'Disputed'},{'GENE ID (HGNC)':'HGNC:4','GENE SYMBOL':'D','CLASSIFICATION':'Limited'}]
            subset=subset_hcm(curation,rows)
            self.assertEqual(subset[0]['MatchMethod'],'exact_HGNC_id');self.assertEqual(subset[0]['CLASSIFICATION'],'Disputed')
            self.assertEqual(subset[1]['Status'],'NA_NO_PROTEIN_CODING_MODEL_MATCH')

    def test_csv_exact_mondo_not_substring(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'curations.csv';p.write_text('metadata\n++++\n"GENE SYMBOL","GENE ID (HGNC)","DISEASE ID (MONDO)"\nA,HGNC:1,MONDO:0005045\nB,HGNC:2,MONDO:00050450\n')
            lock={'csv_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'disease_id':'MONDO:0005045','expected_curation_rows':1}
            self.assertEqual(hcm_curations(p,lock)[0]['GENE SYMBOL'],'A')


if __name__=='__main__':unittest.main()
