import hashlib
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('nmd', Path(__file__).resolve().parents[1] / 'scripts/research_nmd_geometry.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def model():
    sequence = 'ATG' + 'AAA' * 39 + 'TAA' + 'CCC' * 30
    return dict(transcript='ENST000001.2', hgnc_id='HGNC:1', source_url_or_local_reference='synthetic-control',
        reference_sha256=hashlib.sha256(sequence.encode()).hexdigest(), spliced_cdna=sequence, cds_start=1, cds_end=123,
        exon_ends=[90, len(sequence)], edit=dict(cdna_start=10, ref='A', alt='T'), consequence_kind='stop_gained',
        complete_cds=True, reference_genome='GRCh38', orientation='transcript_5prime_to_3prime', biological_relevance_review='UNREVIEWED_SYNTHETIC',
        chromosome='1', genetic_code='standard_nuclear')


class NmdTests(unittest.TestCase):
    def test_new_stop_from_sequence_not_variant_position(self):
        value = model(); result = m.assess(value)
        self.assertEqual(result['mutant_stop_end_cdna'], 12)
        self.assertEqual(result['distance_stop_end_to_final_junction_nt'], 78)
        self.assertTrue(result['more_than_55nt_upstream'])
        self.assertEqual(result['nmd_prediction'], 'NOT_ASSIGNED'); self.assertFalse(result['automatic_PVS1'])
        # Deleting one A changes the reading frame but does not introduce a stop
        # before the supplied sequence ends. Position 10 is never called the PTC.
        value['edit']['alt'] = ''; value['consequence_kind'] = 'frameshift'
        shifted = m.assess(value)
        self.assertIsNone(shifted['mutant_stop_end_cdna'])
        self.assertEqual(shifted['mutant_final_junction_cdna'], 89)

    def test_final_exon_and_threshold_boundary(self):
        value = model(); value['exon_ends'] = [62, len(value['spliced_cdna'])]
        result = m.assess(value); self.assertEqual(result['distance_stop_end_to_final_junction_nt'], 50)
        self.assertFalse(result['more_than_50nt_upstream'])
        value['exon_ends'] = [9, len(value['spliced_cdna'])]
        self.assertEqual(m.assess(value)['geometry'], 'STOP_IN_FINAL_EXON_OR_AT_FINAL_JUNCTION')

    def test_reference_and_exon_contract_fail_closed(self):
        value = model(); value['reference_sha256'] = '0' * 64
        with self.assertRaises(ValueError): m.assess(value)
        value = model(); value['edit']['ref'] = 'C'
        with self.assertRaises(ValueError): m.assess(value)
        value = model(); value['exon_ends'] = [10, len(value['spliced_cdna'])]; value['edit']['ref'] = 'AAA'
        with self.assertRaises(ValueError): m.assess(value)
        self.assertEqual(m.assess({})['status'], 'UNASSESSED_REQUIRES_VALIDATED_TRANSCRIPT_INPUTS')
        value = model(); value['chromosome'] = 'MT'
        self.assertEqual(m.assess(value)['status'], 'UNASSESSED_NON_NUCLEAR_GENETIC_CODE')


if __name__ == '__main__': unittest.main()
