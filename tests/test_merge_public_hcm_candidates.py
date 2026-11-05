import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location('merge_public_sources',SCRIPTS/'merge_public_hcm_candidates.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class MergeTests(unittest.TestCase):
    def test_identical_source_overlap_deduplicated(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);header='chrom\tpos\tref\talt\tCADD_RAW\tCADD_PHRED\n';body='chr1\t1\tA\tG\t1.1\t20\n'
            a=root/'a.tsv';b=root/'b.tsv';a.write_text(header+body);b.write_text(header+body)
            result=m.merge_cadd([a,b],root/'union.tsv');self.assertEqual(result,dict(unique_alleles=1,exact_duplicate_rows=1))
    def test_conflicting_source_values_fail_instead_of_max_or_average(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);header='chrom\tpos\tref\talt\tCADD_RAW\tCADD_PHRED\n';a=root/'a.tsv';b=root/'b.tsv'
            a.write_text(header+'chr1\t1\tA\tG\t1.1\t20\n');b.write_text(header+'chr1\t1\tA\tG\t2.2\t30\n')
            with self.assertRaises(AssertionError):m.merge_cadd([a,b],root/'union.tsv')
            self.assertFalse((root/'union.tsv').exists())


if __name__=='__main__':unittest.main()
