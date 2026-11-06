"""Synthetic MT vectors remain raw; no nuclear quality or ploidy classification."""
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('mt_annotation',ROOT/'scripts/annotate_mitochondrial_consequences.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
q=m.helper(os.environ['SUPPLEMENTARY_HELPER'])
codec=q.module(os.environ['SUPPLEMENTARY_CODEC'],'_mt_public_test_codec')


class MTTests(unittest.TestCase):
    def fixture(self,root):
        header=['##fileformat=VCFv4.2','##contig=<ID=chrM,length=16569>',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Source Mutect2 GT">',
            '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Source AF">',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPUBLIC_MT']
        raw=['chrM','100','original','A','C,G','25','low_allele_frac','TLOD=3,4','GT:AF:AD:DP:HL','0/1:0.1,0.2:70,10,20:100:0.3']
        source=root/'source.vcf.gz';annotated=root/'annotated.vcf.gz';frequency=root/'frequency.tsv.gz'
        with gzip.open(source,'wt') as stream:stream.write('\n'.join(header)+'\n'+'\t'.join(raw)+'\n')
        with gzip.open(annotated,'wt') as stream:
            stream.write('\n'.join(header[:-1]+['##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|SYMBOL|Feature|ALLELE_NUM|HGVSc|HGVSp">',header[-1]])+'\n')
            added=list(raw);added[7]+=';CSQ=C|missense_variant|MT-PUBLIC|TX_PUBLIC|1|c.1A>C|p.Public,G|synonymous_variant|MT-PUBLIC|TX_PUBLIC|2|c.1A>G|p.Public'
            stream.write('\t'.join(added)+'\n')
        row=dict(sample='PUBLIC_MT',CHROM=raw[0],POS=raw[1],REF=raw[3],ALT=raw[4],raw_INFO=raw[7],raw_FORMAT=raw[8],raw_sample=raw[9],
                 gnomAD_MT_match_status='MULTIALLELIC_OR_SYMBOLIC_NO_FLATTEN_NO_JOIN',Helix_match_status='MULTIALLELIC_OR_SYMBOLIC_NO_FLATTEN_NO_JOIN')
        with gzip.open(frequency,'wt',newline='') as stream:
            writer=csv.DictWriter(stream,list(row),delimiter='\t');writer.writeheader();writer.writerow(row)
        return source,annotated,frequency,raw

    def test_raw_multiallelic_heteroplasmy_and_filters_remain_unchanged(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);source,annotated,frequency,raw=self.fixture(root);table=root/'functional.tsv.gz'
            self.assertEqual(m.make_mt_table(q,source,annotated,'PUBLIC_MT',frequency,{'QC':'original'},table),1)
            with gzip.open(table,'rt') as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
            self.assertEqual(len(rows),2)
            for row in rows:
                self.assertEqual(row['source_sample'],raw[9]);self.assertEqual(row['AF'],'0.1,0.2');self.assertEqual(row['HL'],'0.3')
                self.assertEqual(row['GT'],'0/1');self.assertEqual(row['FILTER'],'low_allele_frac')
                self.assertEqual(row['alt'],'C,G');self.assertEqual(row['clinical_classification'],'NOT_ASSIGNED')
                self.assertNotIn('GQ',row)
            self.assertEqual([r['matched_alt'] for r in rows],['C','G'])

    def test_actual_mt_packet_bgzf_persisted_roundtrip(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);source,annotated,frequency,_=self.fixture(root);table=root/'functional.tsv.gz'
            m.make_mt_table(q,source,annotated,'PUBLIC_MT',frequency,{},table)
            proof=q.embed(source,annotated,table,root/'mt.vcf.gz',codec,lambda:None)
            self.assertEqual((proof['records'],proof['rows']),(1,2))
            with gzip.open(root/'mt.vcf.gz','rt') as stream:contents=stream.read()
            self.assertIn('WGS_MT_CONSEQUENCE_EVIDENCE=',contents)
            self.assertNotIn('WGS_SUPPLEMENTARY_EVIDENCE=',contents)

    def test_frequency_for_other_sample_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);source,annotated,frequency,_=self.fixture(root)
            with gzip.open(frequency,'rt') as stream:value=stream.read()
            with gzip.open(frequency,'wt') as stream:stream.write(value.replace('PUBLIC_MT','OTHER_SAMPLE'))
            with self.assertRaises(ValueError):m.make_mt_table(q,source,annotated,'PUBLIC_MT',frequency,{},root/'bad.tsv.gz')

    def test_second_sample_not_silently_merged(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);source,annotated,frequency,_=self.fixture(root)
            with self.assertRaises(ValueError):m.make_mt_table(q,source,annotated,'OTHER_SAMPLE',frequency,{},root/'bad.tsv.gz')


if __name__=='__main__':unittest.main()
