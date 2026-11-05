import importlib.util
from pathlib import Path
import sys
import unittest

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location('revel_public_panel',SCRIPTS/'prepare_revel_public_hcm.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class Reference:
    def base(self,c,p):return 'T'


class RevelPanelTests(unittest.TestCase):
    def row(self):return dict(chr='11',hg19_pos='47353697',grch38_pos='47332146',ref='T',alt='G',aaref='D',aaalt='A',REVEL='0.856',Ensembl_transcriptid='ENST00000545968')
    def test_lifted_exact_allele_retains_raw_aa_and_unversioned_transcript(self):
        row,reason=m.candidate_row(self.row(),{'11':([(47332145,47332146)],[47332145])},Reference())
        self.assertIsNone(reason);self.assertEqual(row['pos'],'47332146');self.assertEqual(row['transcript'],'*')
        self.assertEqual(row['REVEL_source_hg19_pos'],'47353697');self.assertIn('VERSION_UNAVAILABLE',row['REVEL_mapping_status']);self.assertNotIn('REVEL',row)
    def test_unmapped_or_wrong_ref_not_promoted(self):
        value=self.row();value['grch38_pos']='.';row,reason=m.candidate_row(value,{'11':([],[])},Reference());self.assertIsNone(row)
        value=self.row();value['ref']='A';row,reason=m.candidate_row(value,{'11':([(47332145,47332146)],[47332145])},Reference())
        self.assertIsNone(row);self.assertIn('REF_MISMATCH',reason)


if __name__=='__main__':unittest.main()
