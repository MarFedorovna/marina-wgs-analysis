import importlib.util
from pathlib import Path
import sys
import unittest

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location('public_pext',SCRIPTS/'prefetch_pext_public39.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class PextTests(unittest.TestCase):
    def test_hgnc_build_and_canonical_symbol_required(self):
        gene=dict(symbol='MYBPC3',hgnc_id='HGNC:7551');value=dict(reference_genome='GRCh38',symbol='MYBPC3',hgnc_id='7551')
        self.assertEqual(m.validate_gene(value,gene),value)
        for key,bad in [('reference_genome','GRCh37'),('symbol','FHC'),('hgnc_id','HGNC:7577')]:
            altered=dict(value);altered[key]=bad
            with self.assertRaises(AssertionError):m.validate_gene(altered,gene)
    def test_public_query_contains_no_variant_or_case_fields(self):
        self.assertIn('pext',m.QUERY);self.assertNotIn('variants',m.QUERY);self.assertNotIn('subject',m.QUERY);self.assertNotIn('genotype',m.QUERY)


if __name__=='__main__':unittest.main()
