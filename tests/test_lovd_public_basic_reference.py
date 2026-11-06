"""Offline safety/schema controls; no network and no patient data."""
import importlib.util
import json
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[1] / 'scripts' / 'collect_lovd_public_basic_reference.py'
SPEC = importlib.util.spec_from_file_location('lovd_basic', PATH)
lovd = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lovd)


class PublicBasicTests(unittest.TestCase):
    def variant(self, **extra):
        return dict({'symbol': 'TEST1', 'id': 1, 'position_genomic': {'hg38': 'chr1:100'},
                     'Variant/DNA': ['NM_000001.1:c.1A>G'], 'Variant/DBID': 'TEST1_1'}, **extra)

    def test_object_and_array(self):
        row = self.variant()
        for value in [row, [row]]:
            self.assertEqual(lovd.rows(json.dumps(value).encode(), 'variants', 'TEST1'), [row])

    def test_wrong_gene_rejected(self):
        with self.assertRaises(ValueError):
            lovd.rows(json.dumps(self.variant()).encode(), 'variants', 'OTHER')

    def test_case_level_payload_rejected(self):
        with self.assertRaises(ValueError):
            lovd.rows(json.dumps(self.variant(individuals=[])).encode(), 'variants', 'TEST1')

    def test_position_is_not_refalt(self):
        self.assertEqual(lovd.mapping_status(self.variant(), {'refseq_build': 'hg38'}),
                         'UNMAPPED_NO_EXPLICIT_GRCH38_REF_ALT')

    def test_build_not_relabelled(self):
        self.assertEqual(lovd.mapping_status(self.variant(), {'refseq_build': 'hg19'}),
                         'UNMAPPED_SOURCE_BUILD_hg19')

    def test_explicit_grch38_still_requires_fasta(self):
        row = self.variant(variant_genomic={'hg19': 'chr1:g.200A>G', 'hg38': 'chr1:g.100A>G'})
        self.assertEqual(lovd.mapping_status(row, None),
                         'GRCH38_SNV_CANDIDATE_PENDING_LOCAL_FASTA_VALIDATION')

    def test_complex_description_not_guessed(self):
        row = self.variant(variant_genomic={'hg38': 'chr1:g.100_101del'})
        self.assertEqual(lovd.mapping_status(row, None), 'UNMAPPED_EXPLICIT_GRCH38_COMPLEX_DESCRIPTION')

    def test_arbitrary_endpoint_and_variant_input_rejected(self):
        for kind, symbol in [('individuals', 'TEST1'), ('variants', '../TEST1'), ('variants', 'chr1:100A>G')]:
            with self.assertRaises(ValueError):
                lovd.url(kind, symbol)
        with self.assertRaises(ValueError):
            lovd.fetch('https://example.com/private', 10, lovd.END)


if __name__ == '__main__':
    unittest.main()
