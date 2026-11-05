"""Synthetic identity tests; a score must never broadcast across transcripts."""
import importlib.util
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/alphamissense_adapter.py'
spec = importlib.util.spec_from_file_location('alpha_adapter', SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class AlphaTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        header = '#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n'
        canonical = root / 'canonical.tsv'
        canonical.write_text(header + 'chr1\t10\tC\tA\thg38\tP_SYNTH\tENST00000000001.2\tR10L\t0.9\tlikely_pathogenic\n')
        isoforms = root / 'isoforms.tsv'
        isoforms.write_text(header + 'chr1\t10\tC\tA\thg38\t\tENST00000000002.3\tV2L\t0.2\tlikely_benign\n')
        self.db = sqlite3.connect(':memory:')
        self.lookup = m.Lookup(self.db, {'canonical': {'matches': canonical}, 'isoforms': {'matches': isoforms}},
                               {'dataset_version': 'synthetic_fixture'})
        self.csq = dict(Feature='ENST00000000001.2', Amino_acids='R/L', Protein_position='10/100',
                        SWISSPROT='P_SYNTH', Consequence='missense_variant')

    def tearDown(self):
        self.db.close()
        self.temporary.cleanup()

    def test_exact_allele_transcript_version_protein_and_accession(self):
        row = self.lookup.match('1:10:C:A', self.csq)
        self.assertEqual(row['AlphaMissense_score'], .9)
        self.assertEqual(row['AlphaMissense_model_label'], 'likely_pathogenic')
        self.assertEqual(row['AlphaMissense_research_rank'], 1)
        self.assertIn('model_not_ACMG', row['AlphaMissense_match_status'])

    def test_version_mismatch_retains_source_without_assignment(self):
        row = self.lookup.match('1:10:C:A', dict(self.csq, Feature='ENST00000000001.3'))
        self.assertEqual(row['AlphaMissense_score'], '')
        self.assertIn('version_mismatch', row['AlphaMissense_match_status'])
        self.assertEqual(len(json.loads(row['AlphaMissense_allele_predictions_json'])), 2)

    def test_wrong_protein_or_accession_is_not_a_match(self):
        row = self.lookup.match('1:10:C:A', dict(self.csq, Amino_acids='R/V'))
        self.assertEqual(row['AlphaMissense_score'], '')
        self.assertIn('protein_substitution_mismatch', row['AlphaMissense_match_status'])
        row = self.lookup.match('1:10:C:A', dict(self.csq, SWISSPROT='P_OTHER'))
        self.assertEqual(row['AlphaMissense_score'], '')
        self.assertIn('uniprot', row['AlphaMissense_match_status'])

    def test_isoform_without_uniprot_still_requires_exact_transcript_and_protein(self):
        row = self.lookup.match('1:10:C:A', dict(self.csq, Feature='ENST00000000002.3', Amino_acids='V/L', Protein_position='2/100', SWISSPROT=''))
        self.assertEqual(row['AlphaMissense_score'], .2)
        self.assertEqual(row['AlphaMissense_source_dataset'], 'isoforms')
        self.assertEqual(row['AlphaMissense_source_uniprot'], '')

    def test_nonmissense_and_uncovered_are_neither_scores_nor_benign(self):
        row = self.lookup.match('1:10:C:A', dict(self.csq, Consequence='synonymous_variant'))
        self.assertEqual(row['AlphaMissense_score'], '')
        self.assertIn('not_applicable', row['AlphaMissense_match_status'])
        row = self.lookup.match('1:99:C:A', self.csq)
        self.assertEqual(row['AlphaMissense_score'], '')
        self.assertIn('not_benign', row['AlphaMissense_match_status'])

    def test_conflicting_exact_predictions_remain_unassigned(self):
        self.db.execute('INSERT INTO alpha VALUES (?,?,?,?,?,?,?)', ('1:10:C:A', 'ENST00000000001.2', 'R10L', 'P_SYNTH', .1, 'likely_benign', 'isoforms'))
        row = self.lookup.match('1:10:C:A', self.csq)
        self.assertEqual(row['AlphaMissense_score'], '')
        self.assertIn('conflicting', row['AlphaMissense_match_status'])


if __name__ == '__main__':
    unittest.main()
