"""Planning safety only; no patient calling or external services."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/sex_chromosome_plan.py'
spec = importlib.util.spec_from_file_location('xy_plan', SOURCE)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class SexChromosomePlanTests(unittest.TestCase):
    def test_unknown_ploidy_blocks_all_without_diploid_default(self):
        row = dict.fromkeys(m.COLUMNS, 'unknown')
        result = m.sample_plan(row, [])
        self.assertEqual(len(result), 3)
        self.assertTrue(all(item['status'] == 'BLOCKED_BEFORE_CALLING' for item in result))
        self.assertTrue(all(item['declared_ploidy'] is None and item['command_argv'] is None for item in result))

    def test_y_zero_is_skipped_and_never_manufactures_reference(self):
        row = dict(sample_id='PUBLIC_SYNTH', karyotype='XX', x_nonpar_ploidy='2', y_nonpar_ploidy='0', x_par_ploidy='2', evidence_source='synthetic_fixture')
        result = m.sample_plan(row, ['PAR_contract_unverified'])
        y = next(item for item in result if item['partition'] == 'Y_NONPAR')
        self.assertEqual(y['status'], 'SKIP_DECLARED_ZERO_COPY_NOT_REFERENCE')
        self.assertIsNone(y['command_argv'])
        self.assertIn('not_assessed_no_GT', y['genotype'])
        self.assertTrue(all(item['status'] == 'BLOCKED_BEFORE_CALLING' for item in result if item['partition'] != 'Y_NONPAR'))

    def test_official_par_boundaries_and_complement_have_no_gap_or_overlap(self):
        lock = m.load(m.CODE / 'config/sex-chromosomes.resources.lock.json')
        parts = m.partitions(lock)
        self.assertEqual(parts['X_PAR'][0], ('chrX', 10001, 2781479))
        self.assertEqual(parts['X_NONPAR'][0], ('chrX', 1, 10000))
        self.assertEqual(parts['X_NONPAR'][1][1], 2781480)
        self.assertEqual(parts['Y_NONPAR'][-1][1], 57217416)
        for chrom in ['chrX', 'chrY']:
            nonpar = parts['X_NONPAR' if chrom == 'chrX' else 'Y_NONPAR']
            par = lock['par_one_based_closed'][chrom]
            self.assertEqual(sum(b-a+1 for _, a, b in nonpar) + sum(b-a+1 for a,b in par), lock['primary_contig_lengths'][chrom])
            self.assertFalse(any(a <= d and c <= b for _,a,b in nonpar for c,d in par))

    def test_fai_mask_probe_distinguishes_n_from_sequence(self):
        with tempfile.TemporaryDirectory() as temporary:
            p = Path(temporary) / 'public.fa'; p.write_bytes(b'>chrY\nNNNNACGT\n')
            row = ['chrY', '8', '6', '8', '9']
            self.assertEqual(m.count_range(p, row, 1, 4), dict(bases=4, N=4, ACGT=0))
            self.assertEqual(m.count_range(p, row, 5, 8), dict(bases=4, N=0, ACGT=4))


if __name__ == '__main__':
    unittest.main()
