import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('tissue_protein',Path(__file__).resolve().parents[1]/'scripts/prepare_tissue_protein_references.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class TissueProteinTests(unittest.TestCase):
    def test_uniprot_exact_hgnc_crossref_and_reviewed_status(self):
        entry=dict(organism=dict(taxonId=9606),entryType='UniProtKB reviewed (Swiss-Prot)',uniProtKBCrossReferences=[dict(database='HGNC',id='HGNC:7551'),dict(database='GeneID',id='4607')])
        self.assertEqual(m.uniprot_gene_ids(entry),{'HGNC:7551'})
        entry['entryType']='UniProtKB unreviewed (TrEMBL)'
        with self.assertRaises(AssertionError):m.uniprot_gene_ids(entry)
    def test_gtex_needs_hgnc_and_canonical_symbol_not_alias_guess(self):
        genes=[dict(hgnc_id='HGNC:7551',symbol='MYBPC3')];row=dict(description='protein [Source:HGNC Symbol;Acc:HGNC:7551]',geneSymbol='MYBPC3')
        self.assertEqual(m.gtex_identity(row,genes),genes[0]);row['geneSymbol']='FHC';self.assertIsNone(m.gtex_identity(row,genes))


if __name__=='__main__':unittest.main()
