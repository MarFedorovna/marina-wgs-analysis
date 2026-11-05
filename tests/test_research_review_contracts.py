import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('contracts', ROOT / 'scripts/research_review_contracts.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class ReviewContractTests(unittest.TestCase):
    def test_revel_allele_aa_transcript_no_version_invention(self):
        source = dict(assembly='GRCh38', chrom='chr11', pos=1, ref='T', alt='G', REVEL_version='1.3',
            REVEL_source_transcript_ids='ENST000001', REVEL_source_aa_ref='D', REVEL_source_aa_alt='A', REVEL_reference_score='.856')
        current = dict(assembly='GRCh38', chrom='11', pos=1, ref='T', alt='G', transcript='ENST000001.2', amino_acids='D/A')
        result = m.revel_alignment(source, current)
        self.assertEqual(result['raw_reference_score'], .856)
        self.assertIsNone(result['canonical_transcript_score']); self.assertFalse(result['automatic_PP3_BP4'])
        current['amino_acids'] = 'D/V'; self.assertEqual(m.revel_alignment(source, current)['status'], 'AMINO_ACID_MISMATCH')
        current['alt'] = 'A'; self.assertEqual(m.revel_alignment(source, current)['status'], 'ALLELE_OR_ASSEMBLY_MISMATCH')

    def test_actual_source_cardinality_three_transcripts_two_fathmm_not_zip(self):
        raw = dict(reference_genome='GRCh38', chrom='11', pos=1, ref='T', alt='G', dbnsfp=dict(
            ensembl=dict(transcriptid=['ENST1', 'ENST2', 'ENST3']), aa=dict(pos=[1247, 1247, 1247]),
            fathmm=dict(score=[-1.52, -1.52]), revel=dict(score=[.856, .856, .856])))
        current = dict(assembly='GRCh38', chrom='11', pos=1, ref='T', alt='G')
        self.assertEqual(m.dbnsfp_alignment(raw, 'fathmm', current)['status'], 'NONALIGNED_CARDINALITY_DO_NOT_ZIP')
        result = m.dbnsfp_alignment(raw, 'revel', current)
        self.assertIn('CONTRACT_AND_VERSIONS_UNVERIFIED', result['status']); self.assertIsNone(result['canonical_transcript_score'])
        self.assertEqual(m.dbnsfp_alignment(raw, 'spliceai', current)['status'], 'MISSING_SOURCE_NOT_ZERO_OR_BENIGN')

    def test_blank_template_and_no_inferred_family_or_phenotype(self):
        template = json.loads((ROOT / 'templates/phenotype-pedigree-input-v1.json').read_text())
        result = m.validate_phenotype_pedigree(template)
        self.assertEqual(result['status'], 'PASS_INPUT_STRUCTURE_ONLY_NO_CLINICAL_INFERENCE')
        self.assertEqual(result['inferred_relationships'], 0)
        changed = copy.deepcopy(template)
        changed['subjects'][0]['phenotypes'] = [dict(hpo_id='HP:0001639', status='present', source=None)]
        self.assertTrue(m.validate_phenotype_pedigree(changed)['errors'])

    def test_denovo_phase_need_evidence(self):
        template = json.loads((ROOT / 'templates/phenotype-pedigree-input-v1.json').read_text())
        template['segregation_records'] = [dict(subject_id='SUBJECT_A', de_novo_status='confirmed', phase='trans')]
        result = m.validate_phenotype_pedigree(template)
        self.assertTrue(any('de novo' in e for e in result['errors'])); self.assertTrue(any('cis/trans' in e for e in result['errors']))


if __name__ == '__main__': unittest.main()
