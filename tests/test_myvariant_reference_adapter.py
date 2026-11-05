import importlib.util
from pathlib import Path
import sys
import unittest

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location('myvariant_reference_adapter',SCRIPTS/'export_myvariant_reference_adapter.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class Reference:
    def base(self,c,p):return 'T'


class ApiAdapterTests(unittest.TestCase):
    def hit(self):return dict(_id='chr11:g.47332146T>G',dbnsfp=dict(genename='MYBPC3',chrom='11',ref='T',alt='G',hg38=dict(start=47332146,end=47332146),ensembl=dict(transcriptid=['ENST1','ENST2','ENST3']),fathmm=dict(score=[-1.52,-1.52],pred=['D','D']),cadd=dict(phred=25.0)))
    def test_nonaligned_vectors_remain_raw_and_do_not_override_cadd(self):
        row,reason=m.row_for_hit(self.hit(),'MYBPC3',Reference(),Path('/public/page.json'),'a'*64)
        self.assertIsNone(reason);self.assertEqual(row['transcript'],'*')
        self.assertIn('-1.52',row['FATHMM_dbNSFP48_json']);self.assertNotIn('FATHMM',row);self.assertNotIn('CADD_PHRED',row)
        self.assertEqual(row['SpliceAI_dbNSFP48_json'],'');self.assertIn('NONALIGNED',row['transcript_score_status'])
    def test_hg19_id_mismatch_excluded_even_with_valid_structured_hg38(self):
        hit=self.hit();hit['_id']='chr11:g.47353697T>G'
        row,reason=m.row_for_hit(hit,'MYBPC3',Reference(),Path('/public/page.json'),'a'*64)
        self.assertIsNone(row);self.assertEqual(reason,'RAW_ID_ASSEMBLY_MISMATCH')


if __name__=='__main__':unittest.main()
