import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from recover_mitochondria import merge_assessment


class AssessmentTests(unittest.TestCase):
    def test_boundary_and_existing_filters_and_genotypes_retained(self):
        header = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n'
        rows = [f'chrM\t{p}\t.\tA\tC\t.\t{f}\t.\tGT:DP:AD:AF\t0/1:20:10,10:0.5\n'
                for p, f in ((100, 'PASS'), (5000, 'weak_evidence'), (16560, 'PASS'))]
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / 'a').write_text(header + ''.join(rows))
            (root / 'b').write_text(header + rows[1].replace('weak_evidence', 'alignment'))
            counts = merge_assessment(root / 'a', root / 'b', root / 'out')
            actual = [x.split('\t') for x in (root / 'out').read_text().splitlines() if not x.startswith('#')]
            self.assertEqual(counts, {'assessed': 1, 'unavailable': 2, 'alignment_artifact': 1, 'rows': 3})
            self.assertEqual(actual[1][6], 'weak_evidence;alignment')
            self.assertTrue(all(r[9] == '0/1:20:10,10:0.5' for r in actual))
            self.assertTrue(all(actual[i][6] == 'NUMTArtifactCheckUnavailable' for i in (0, 2)))


if __name__ == '__main__':
    unittest.main()
