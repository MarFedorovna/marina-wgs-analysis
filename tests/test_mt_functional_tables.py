"""Public MT sharing controls: record presence, multiallelic identity and missing semantics."""
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('mt_views',ROOT/'scripts/export_mt_functional_tables.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
q=m.helper(os.environ['SUPPLEMENTARY_HELPER'])
codec=q.module(os.environ['SUPPLEMENTARY_CODEC'],'_mt_views_test_codec')
SAMPLES=['PUBLIC_A','PUBLIC_B','PUBLIC_C','PUBLIC_D']


class MTFunctionalTests(unittest.TestCase):
    def test_genotype_zero_does_not_erase_observed_heteroplasmy_record(self):
        values={};key=('chrM',100,'A','G')
        m.observed_allele(values,key,SAMPLES[0],dict(GT='0/0',AF_vector='.12',HL_raw='.12',FILTER='PASS'))
        m.observed_allele(values,key,SAMPLES[1],dict(GT='0/1',AF_vector='.03',HL_raw='.',FILTER='low_allele_frac'))
        row=m.summarize(key,values[key],SAMPLES)
        self.assertEqual(row['observed_call_count'],2);self.assertEqual(row['site_PASS_count'],1)
        self.assertTrue(m.belongs('shared_ge2','observed',row,SAMPLES))
        self.assertFalse(m.belongs('shared_ge2','site_PASS',row,SAMPLES))
        context=json.loads(row['sample_source_context_json'])
        self.assertEqual(context[SAMPLES[0]]['original_records'][0]['GT'],'0/0')
        self.assertEqual(context[SAMPLES[2]]['state'],'NOT_OBSERVED_NOT_ASSESSED_NO_REFERENCE_INFERENCE')

    def test_transcript_effects_are_unioned_per_exact_allele_not_two_alleles(self):
        values={};key=('chrM',100,'A','G')
        value=m.observed_allele(values,key,SAMPLES[0],dict(FILTER='PASS'))
        for term in ['missense_variant','synonymous_variant']:
            m.attach_effect(value,dict(sample=SAMPLES[0],gene='MT-PUBLIC',consequence=term))
        row=m.summarize(key,value,SAMPLES)
        self.assertEqual(row['observed_call_count'],1);self.assertEqual(row['transcript_effect_rows'],2)
        self.assertEqual(json.loads(row['consequences_json']),['missense_variant','synonymous_variant'])

    def fixture(self,root,sample='PUBLIC_A'):
        source=root/'source.vcf.gz';table=root/'table.tsv.gz';annotated=root/'annotated.vcf.gz'
        headers=['##fileformat=VCFv4.2','##contig=<ID=chrM,length=16569>',
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="SourceGT">',
                 '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample]
        raw=['chrM','100','PUBLIC_VARIANT','A','C,G','20','PASS','TLOD=3,4','GT:AF:HL','0/0:0.1,0.2:0.3']
        for path in [source,annotated]:
            with gzip.open(path,'wt') as stream:stream.write('\n'.join(headers)+'\n'+'\t'.join(raw)+'\n')
        fields=['sample','chrom','pos','ref','alt','source_FORMAT','source_sample','allele_number','matched_alt','consequence','gene','frequency_context_json']
        with gzip.open(table,'wt',newline='') as stream:
            writer=csv.DictWriter(stream,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
            for num,alt,term in [('1','C','missense_variant'),('2','G','synonymous_variant')]:
                writer.writerow(dict(sample=sample,chrom='chrM',pos='100',ref='A',alt='C,G',source_FORMAT=raw[8],source_sample=raw[9],
                                     allele_number=num,matched_alt=alt,consequence=term,gene='MT-PUBLIC',frequency_context_json='{"no_flatten":true}'))
        q.TAG='WGS_MT_CONSEQUENCE_EVIDENCE';q.SCHEMA='WGS_MT_CONSEQUENCE_SCHEMA'
        output=root/'final.vcf.gz';proof=q.embed(source,annotated,table,output,codec,lambda:None)
        return output,table,proof

    def test_actual_packet_multialt_identity_and_all_vectors_preserved(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);vcf,table,proof=self.fixture(root)
            values={};master=root/'master.tsv'
            with master.open('w',newline='') as stream:
                writer=csv.writer(stream,delimiter='\t',lineterminator='\n');writer.writerow(proof['fields'])
                result=m.collect(q,codec,vcf,table,'PUBLIC_A',proof['ordered_vector_sha256'],values,[writer],lambda:None)
            self.assertEqual(set(values),{('chrM',100,'A','C'),('chrM',100,'A','G')})
            self.assertEqual(values[('chrM',100,'A','C')]['terms'],{'missense_variant'})
            self.assertEqual(values[('chrM',100,'A','G')]['terms'],{'synonymous_variant'})
            self.assertEqual(result['transcript_rows'],2)
            self.assertEqual(values[('chrM',100,'A','G')]['samples']['PUBLIC_A'][0]['AF_vector'],'0.1,0.2')
            with gzip.open(table,'rt') as a,master.open() as b:
                self.assertEqual(list(csv.reader(a,delimiter='\t')),list(csv.reader(b,delimiter='\t')))

    def test_missing_table_vector_is_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);vcf,table,proof=self.fixture(root)
            with gzip.open(table,'wt') as stream:stream.write('\t'.join(proof['fields'])+'\n')
            with self.assertRaises(ValueError):m.collect(q,codec,vcf,table,'PUBLIC_A',proof['ordered_vector_sha256'],{},[],lambda:None)

    def test_different_alt_and_different_position_are_not_shared(self):
        values={}
        for sample,key in zip(SAMPLES,[('chrM',100,'A','G'),('chrM',100,'A','C'),('chrM',101,'A','G'),('chrM',102,'A','G')]):
            m.observed_allele(values,key,sample,dict(FILTER='PASS'))
        self.assertTrue(all(m.summarize(k,v,SAMPLES)['observed_call_count']==1 for k,v in values.items()))

    def test_actual_four_sample_export_all_types_and_shared_scopes(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);source=root/'source';source.mkdir();samples=[];artifacts={}
            for sample in SAMPLES:
                folder=source/sample;folder.mkdir();vcf,table,proof=self.fixture(folder,sample)
                final=folder/'mitochondria.fully-consequence-annotated.vcf.gz';vcf.rename(final)
                Path(str(vcf)+'.tbi').rename(Path(str(final)+'.tbi'))
                long=folder/'MT-functional-all-transcripts.tsv.gz';table.rename(long)
                samples.append(dict(sample=sample,packet_proof=proof))
                for path in [final,Path(str(final)+'.tbi'),long]:artifacts[str(path.relative_to(source))]=q.pin(path)
            manifest=dict(status='PUBLIC_SYNTHETIC_INPUT_FIXTURE',samples=samples,patient_data_used=False)
            receipt=dict(status='PUBLIC_SYNTHETIC_INPUT_FIXTURE',artifacts=artifacts)
            q.write(source/'manifest.json',manifest);q.write(source/'workflow-receipt.json',receipt)
            p=dict(input_root=str(source),codec=q.pin(Path(os.environ['SUPPLEMENTARY_CODEC'])),
                   producer=q.pin(ROOT/'scripts/export_mt_functional_tables.py'),reviewed_profile_sha256='PUBLIC_CONTROL_ONLY')
            output=root/'export'
            result=m.produce(q,p,manifest,receipt,output,lambda:None)
            self.assertEqual(result['unique_exact_alleles'],2)
            with (output/'tables/functional-category-index.tsv').open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
            all4=next(r for r in rows if r['scope']=='all4' and r['level']=='site_PASS' and r['consequence']=='all.functional')
            self.assertEqual(all4['observed_unique_alleles'],'2')
            empty=next(r for r in rows if r['scope']=='cohort' and r['level']=='observed' and r['consequence']=='stop_gained')
            self.assertEqual(empty['observed_unique_alleles'],'0')
            for item in result['artifacts'].values():q.verify(item)
            final_receipt=q.load(output/'workflow-receipt.json')
            self.assertEqual(final_receipt['artifacts']['manifest.json'],q.pin(output/'manifest.json'))

    def test_unresolved_multialt_effect_is_not_guessed(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);vcf,table,proof=self.fixture(root)
            with gzip.open(table,'rt') as stream:
                reader=csv.DictReader(stream,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
            for row in rows:row['allele_number']='';row['matched_alt']=''
            with gzip.open(table,'wt',newline='') as stream:
                writer=csv.DictWriter(stream,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
            ambiguous=root/'ambiguous.vcf.gz'
            proof=q.embed(root/'source.vcf.gz',root/'annotated.vcf.gz',table,ambiguous,codec,lambda:None)
            alleles={};result=m.collect(q,codec,ambiguous,table,'PUBLIC_A',proof['ordered_vector_sha256'],alleles,[],lambda:None)
            self.assertEqual(result['unresolved_multiallelic_transcript_rows'],2)
            self.assertTrue(all(not value['terms'] for value in alleles.values()))


if __name__=='__main__':unittest.main()
