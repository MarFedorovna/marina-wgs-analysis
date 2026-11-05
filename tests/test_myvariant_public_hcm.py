import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('public_myvariant',Path(__file__).resolve().parents[1]/'scripts/prefetch_myvariant_public_hcm.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class Reference:
    def base(self,chrom,pos):return 'T' if (chrom,pos)==('11',47332146) else 'N'


class PublicMyVariantTests(unittest.TestCase):
    def hit(self):return dict(_id='chr11:g.47332146T>G',dbnsfp=dict(genename=['MYBPC3'],chrom='11',ref='T',alt='G',hg38=dict(start=47332146,end=47332146),fathmm=dict(score=[-1.52,-1.52])))
    def test_exact_grch38_ref_alt_and_public_gene_control(self):
        result=m.validate_hit(self.hit(),'MYBPC3',Reference())
        self.assertEqual(result['variant_key'],'GRCh38:chr11:47332146:T:G')
        self.assertIn('TRANSCRIPT_VERSION_UNRESOLVED',result['status'])
    def test_default_hg19_raw_id_not_treated_as_hg38(self):
        hit=self.hit();hit['_id']='chr11:g.47353697T>G'
        self.assertEqual(m.validate_hit(hit,'MYBPC3',Reference())['status'],'RAW_ID_ASSEMBLY_MISMATCH')
    def test_gene_alias_ref_or_ambiguous_position_not_guessed(self):
        hit=self.hit();self.assertEqual(m.validate_hit(hit,'FHC',Reference())['status'],'GENE_ID_MISMATCH_NO_SCORE_JOIN')
        hit['dbnsfp']['ref']='A';hit['_id']='chr11:g.47332146A>G'
        self.assertEqual(m.validate_hit(hit,'MYBPC3',Reference())['status'],'REFERENCE_BASE_MISMATCH_NO_SCORE_JOIN')
        hit['dbnsfp']['hg38']['start']=[47332146]
        self.assertEqual(m.validate_hit(hit,'MYBPC3',Reference())['status'],'MISSING_OR_AMBIGUOUS_GRCH38_SNV_MAPPING')


if __name__=='__main__':unittest.main()
