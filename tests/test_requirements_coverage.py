import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('coverage_builder',Path(__file__).resolve().parents[1]/'scripts/build_requirements_coverage.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class CoverageTests(unittest.TestCase):
    def test_all_audited_resources_and_sources_visible(self):
        result=m.coverage()
        self.assertEqual(result['registry_source_count'],69)
        self.assertEqual(len({r['id'] for r in result['resources']}),69)
        docs={d['id']:d for d in result['audited_sources']}
        self.assertTrue(all('M%02d'%i in docs for i in range(1,42)))
        self.assertEqual(docs['M26']['status'],'EMPTY_FAILED_NOTES_NO_REQUIREMENT_INFERRED')
        self.assertTrue(all('U%02d'%i in docs for i in range(1,20)))
        self.assertTrue(all('G%d'%i in docs for i in range(42,45)))

    def test_withdrawn_kegg_not_reenabled_by_all(self):
        result=m.coverage()
        kegg=next(r for r in result['resources'] if r['id']=='kegg_enrichment_withdrawn')
        self.assertEqual(kegg['actual_status'],'WITHDRAWN_NOT_RUN')
        self.assertTrue(any(r['source_refs']==['U07'] and r['category']=='withdrawn' for r in result['requirements']))

    def test_missing_resources_not_marked_installed(self):
        result=m.coverage()
        self.assertFalse(any(r['actual_status'].startswith(('PASS','READY','INSTALLED')) for r in result['resources']))
        self.assertIn('top10 and top20 per sample and cohort',result['mandatory_research_outputs'])
        self.assertIn('No automatic ACMG/PVS1/LoF/GoF assignment',result['interpretation_gates'])


if __name__=='__main__':unittest.main()
