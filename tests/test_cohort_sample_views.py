import csv
import datetime as dt
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('views',ROOT/'scripts/export_cohort_sample_views.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
ranking=ROOT/'scripts/annotation_stages.py'
rs=importlib.util.spec_from_file_location('rank',ranking);rank=importlib.util.module_from_spec(rs);rs.loader.exec_module(rank)
try:import pysam
except ImportError:pysam=None


class ViewPureTests(unittest.TestCase):
    def test_expired_deadline_cannot_read_or_actuate(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            with self.assertRaisesRegex(ValueError,'deadline'):
                m.export(root/'missing-workflow',root/'output',root/'missing-ranking','0'*64,deadline=m.now()-dt.timedelta(seconds=1))
            self.assertFalse((root/'output').exists())
    def test_partial_unknown_quality_counts_and_source_info(self):
        flags=m.call_flags('1/.','PASS','PASS')
        self.assertTrue(flags['raw']);self.assertFalse(flags['qc']);self.assertTrue(flags['unknown'])
        self.assertEqual((flags['AC'],flags['AN']),(1,1))
        info=m.parse_info(m.derived_info(dict(AC='5',AN='8',AF='.625',DP='100',CSQ='raw'),flags))
        self.assertEqual((info['COHORT_AC'],info['COHORT_AN'],info['COHORT_AF']),('5','8','.625'))
        self.assertEqual((info['AC'],info['AN'],info['AF']),('1','1','1'))
        self.assertEqual(info['COHORT_DP'],'100');self.assertEqual(info['CSQ'],'raw')
        self.assertTrue(m.call_flags('1|0','PASS','PASS')['qc'])
        with self.assertRaises(ValueError):m.call_flags('0/2','PASS','PASS')

    def test_exact_allele_and_shared_membership_not_gene_overlap(self):
        self.assertNotEqual(m.exact_key('chr1',1,'A','T'),m.exact_key('chr1',2,'A','T'))
        calls={s:dict(flags=m.call_flags(gt,ft,'PASS')) for s,gt,ft in [('A','0/1','PASS'),('B','1/.','PASS'),('C','0/1','LowDP'),('D','./.','PASS')]}
        values=m.membership(calls,list(calls),['A'])
        self.assertEqual((values['carrier_count'],values['pct_of4'],values['raw_observed_alt_count']),(1,25,3))
        self.assertEqual(m.scope_members('shared_ge2',calls,{'A'}),(True,False,False))
        self.assertEqual(values['assessed_count'],1)


@unittest.skipUnless(pysam,'pysam required for actual BGZF/TBI integration controls')
class ViewIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.workflow=self.root/'workflow';self.stage6=self.workflow/'stages/06_annotate';self.stage7=self.workflow/'stages/07_report'
        self.stage6.mkdir(parents=True);self.stage7.mkdir()
        self.samples=['SYNTHETIC_A','SYNTHETIC_B','SYNTHETIC_C','SYNTHETIC_D']
        self.records=[(10,['0/1','0/1','0/0','./.'],['PASS']*4,'stop_gained'),
            (20,['1/.','0/0','0/1','0/0'],['PASS','PASS','LowDP','PASS'],'stop_lost'),
            (30,['0/1']*4,['PASS']*4,'frameshift_variant'),
            (40,['0/1','0/0','0/0','0/0'],['PASS']*4,'missense_variant'),
            (50,['0/0','0/1','0/0','0/0'],['PASS']*4,'missense_variant'),
            (60,['./.']*4,['PASS']*4,'intron_variant')]
        vcf=self.stage6/'annotated.vcf.gz'
        header='##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        for tag,number,kind in [('AC','A','Integer'),('AN','1','Integer'),('AF','A','Float'),('DP','1','Integer')]:header+='##INFO=<ID=%s,Number=%s,Type=%s,Description="source cohort">\n'%(tag,number,kind)
        for tag,number,kind in [('GT','1','String'),('FT','1','String'),('DP','1','Integer'),('AD','R','Integer'),('GQ','1','Integer')]:header+='##FORMAT=<ID=%s,Number=%s,Type=%s,Description="source">\n'%(tag,number,kind)
        header+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(self.samples)+'\n'
        self.rows=[]
        with pysam.BGZFile(str(vcf),'wb') as stream:
            stream.write(header.encode())
            for pos,gts,fts,effect in self.records:
                ac=sum(gt.count('1') for gt in gts);an=sum(c in '01' for gt in gts for c in gt)
                stream.write(('chr1\t%d\t.\tA\tT\t90\tPASS\tAC=%d;AN=%d;AF=%s;DP=80\tGT:FT:DP:AD:GQ\t'%(pos,ac,an,str(ac/an) if an else '.')+'\t'.join(gt+':'+ft+':20:10,10:50' for gt,ft in zip(gts,fts))+'\n').encode())
                for tx in ['ENST_SYNTHETIC.1','ENST_SYNTHETIC.2']:
                    for sample,gt,ft in zip(self.samples,gts,fts):
                        row={f:'' for f in rank.FIELDS};state=rank.sample_state(gt,ft,'PASS')
                        row.update(sample=sample,variant_id=rank.variant_key('chr1',pos,'A','T'),chrom='chr1',pos=str(pos),ref='A',alt='T',
                            GT=gt,FT=ft,DP='20',AD='10,10',GQ='50',site_filter='PASS',sample_state=state,
                            gene='SYNTHETIC_GENE',transcript=tx,consequence=effect,
                            candidate_eligible='1' if state=='alternate' and effect!='intron_variant' else '0',
                            consequence_rank='5' if tx.endswith('.1') else '10',local_evidence_json='{"FATHMM_source":"RAW_UNALIGNED"}',
                            classification_status='not_classified',nmd_assessment='unassessed',mechanism_status='not_inferred')
                        self.rows.append(row)
        pysam.tabix_index(str(vcf),preset='vcf')
        self.write_master()
        self.receipts()

    def write_master(self):
        with gzip.open(self.stage7/'master_evidence.tsv.gz','wt',newline='') as stream:
            writer=csv.DictWriter(stream,rank.FIELDS,delimiter='\t');writer.writeheader();writer.writerows(self.rows)

    def receipts(self):
        for stage in [self.stage6,self.stage7]:
            artifacts={p.name:dict(bytes=p.stat().st_size,sha256=m.sha(p)) for p in stage.iterdir() if p.is_file() and p.name!='workflow-receipt.json'}
            (stage/'workflow-receipt.json').write_text(json.dumps(dict(status='PASS',artifacts=artifacts)))

    def run_export(self):
        return m.export(self.workflow,self.root/'output',ranking,m.sha(ranking),synthetic=True)

    def read_table(self,path):
        opener=gzip.open if str(path).endswith('.gz') else open
        with opener(path,'rt',newline='') as stream:return list(csv.DictReader(stream,delimiter='\t'))

    def test_bgzf_index_counts_partial_and_cohort_provenance(self):
        result=self.run_export();out=self.root/'output'
        self.assertEqual(result['unique_annotated_alleles'],6);self.assertEqual(result['master_transcript_rows'],48)
        with pysam.VariantFile(str(out/'vcf/sample_0001.annotated.vcf.gz')) as vcf:
            records=list(vcf.fetch());self.assertEqual([r.pos for r in records],[10,20,30,40])
            partial=next(r for r in records if r.pos==20)
            self.assertEqual(partial.samples[self.samples[0]]['GT'],(1,None));self.assertEqual(partial.info['AN'],1)
            self.assertEqual(partial.info['AF'],(1.0,));self.assertIn('COHORT_AC',partial.info)
        matrix=self.read_table(out/'tables/allele-membership-matrix.tsv.gz')
        first=next(r for r in matrix if r['pos']=='10');self.assertEqual((first['carrier_count'],first['pct_of4'],first['assessed_count']),('2','50','3'))
        self.assertIn(self.samples[3],json.loads(first['unknown_samples_json']))
        self.assertEqual(m.sha(out/'masters/cohort.master.tsv.gz'),m.sha(self.stage7/'master_evidence.tsv.gz'))
        sample_rows=self.read_table(out/'masters/sample_0001.master.tsv.gz')
        self.assertEqual(len(sample_rows),8);self.assertTrue(all(r['sample']==self.samples[0] for r in sample_rows))
        self.assertIn('1/.',{r['GT'] for r in sample_rows})
        with self.assertRaisesRegex(ValueError,'refusing overwrite'):self.run_export()

    def test_top_unique_exact_sealed_order_shared_and_no_padding(self):
        result=self.run_export();out=self.root/'output'
        top=self.read_table(out/'tables/cohort.top20.tsv')
        self.assertEqual(len(top),4);self.assertEqual(len({r['allele_id'] for r in top}),4)
        self.assertTrue(all(r['transcript']=='ENST_SYNTHETIC.2' for r in top))
        shared=self.read_table(out/'tables/shared_ge2.top20.tsv');self.assertEqual({r['pos'] for r in shared},{'10','30'})
        all4=self.read_table(out/'tables/all4.top20.tsv');self.assertEqual([r['pos'] for r in all4],['30'])
        sample_top=self.read_table(out/'tables/sample_0002.top20.tsv')
        self.assertTrue(all(r['sample']==self.samples[1] and r['GT']=='0/1' for r in sample_top))
        self.assertTrue(all(r['ACMG_class']=='NOT_ASSIGNED' for r in top))
        self.assertEqual(result['counts']['cohort']['consequences']['stop_lost'],dict(raw=1,qc=0,eligible=0))
        self.assertFalse(result['automatic_ACMG'])
        source_best={}
        for row in self.rows:
            if row['candidate_eligible']=='1':
                key=row['variant_id']
                if key not in source_best or rank.rank_key(row)<rank.rank_key(source_best[key]):source_best[key]=row
        expected=sorted(source_best.values(),key=rank.rank_key)[:20]
        self.assertEqual([(r['variant_id'],r['sample'],r['transcript']) for r in top],[(r['variant_id'],r['sample'],r['transcript']) for r in expected])

    def test_source_mismatch_nonpass_and_no_overwrite(self):
        receipt=self.stage7/'workflow-receipt.json';value=json.loads(receipt.read_text());value['status']='WAITING';receipt.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError,'NOT_READY'):self.run_export()
        self.assertFalse((self.root/'output').exists())
        self.receipts();self.rows[0]['GT']='0/0';self.write_master();self.receipts()
        with self.assertRaisesRegex(ValueError,'genotype/quality'):self.run_export()
        self.assertFalse((self.root/'output/manifest.json').exists())

    def test_empty_sample_has_real_readable_bgzf_and_index(self):
        vcf=self.stage6/'annotated.vcf.gz'
        with gzip.open(vcf,'rt') as stream:lines=stream.readlines()
        with pysam.BGZFile(str(vcf),'wb') as stream:
            for line in lines:
                if not line.startswith('#'):
                    fields=line.rstrip('\n').split('\t');fields[-1]='0/0:PASS:20:10,10:50';line='\t'.join(fields)+'\n'
                stream.write(line.encode())
        Path(str(vcf)+'.tbi').unlink();pysam.tabix_index(str(vcf),preset='vcf')
        for row in self.rows:
            if row['sample']==self.samples[3]:row.update(GT='0/0',FT='PASS',sample_state='reference',candidate_eligible='0')
        self.write_master();self.receipts();result=self.run_export()
        self.assertEqual(result['sample_vcf_records'][self.samples[3]],0)
        with pysam.VariantFile(str(self.root/'output/vcf/sample_0004.annotated.vcf.gz')) as check:self.assertEqual(list(check.fetch()),[])
        self.assertEqual(self.read_table(self.root/'output/tables/sample_0004.top20.tsv'),[])
        self.assertEqual(self.read_table(self.root/'output/masters/sample_0004.master.tsv.gz'),[])

    def test_receipt_hash_tamper_rejected_before_output(self):
        path=self.stage6/'annotated.vcf.gz.tbi';path.write_bytes(path.read_bytes()+b'X')
        with self.assertRaisesRegex(ValueError,'artifact differs'):self.run_export()
        self.assertFalse((self.root/'output').exists())


if __name__=='__main__':unittest.main()
