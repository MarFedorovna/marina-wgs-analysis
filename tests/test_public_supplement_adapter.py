import importlib.util
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('supplement_adapter',Path(__file__).resolve().parents[1]/'scripts/public_supplement_adapter.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class SupplementTests(unittest.TestCase):
    def test_exact_hgnc_not_symbol_alias_or_substring(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'genes.xml'
            path.write_text('<JDBOR><Disorder><OrphaCode>123</OrphaCode><Name>Disease</Name><DisorderGeneAssociationList><DisorderGeneAssociation><Gene><Symbol>MYH7</Symbol><ExternalReferenceList><ExternalReference><Source>HGNC</Source><Reference>7577</Reference></ExternalReference></ExternalReferenceList></Gene><DisorderGeneAssociationType><Name>Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType></DisorderGeneAssociation></DisorderGeneAssociationList></Disorder></JDBOR>')
            self.assertEqual(list(m.gene_rows(path))[0][0],'HGNC:7577')
        for invalid in ['MYH7','HGNC:75770x','7577']:
            with self.assertRaises(ValueError):m.require_id(invalid,'HGNC')

    def spec(self, current=True, disease='MONDO:0005045'):
        return dict(entId='GN095',entContent=dict(version='1.0',states=[dict(current=current,name='Released')]),
            ld=dict(Disease=[dict(entId=disease)],Gene=[dict(entContent=dict(HGNC=dict(hgnc_id='HGNC:7551',symbol='MYBPC3')))],CriteriaCode=[dict(entId='PVS1',entContent=dict(applicability='requires review'))]))

    def test_only_current_released_exact_mondo_specs_with_rules(self):
        good=self.spec();wrong=self.spec(disease='MONDO:00050450');draft=self.spec(current=False)
        rows=list(m.released_cspec({'data':[good,wrong,draft]}))
        self.assertEqual(len(rows),1);self.assertEqual(rows[0]['hgnc_id'],'HGNC:7551')
        self.assertEqual(rows[0]['version'],'1.0');self.assertFalse(rows[0]['automatic_acmg_application'])
        self.assertEqual(len(rows[0]['rules']),1)

    def test_duplicate_spec_refused(self):
        with self.assertRaises(ValueError):list(m.released_cspec({'data':[self.spec(),self.spec()]}))

    def test_no_false_installed_without_pass_manifest(self):
        with tempfile.TemporaryDirectory() as t:
            source=Path(t)/'source';source.mkdir();target=Path(t)/'candidate'
            (source/'source-manifest.json').write_text('{"status":"DOWNLOADING","resources":[]}')
            with self.assertRaises(ValueError):m.build(source,target)
            self.assertFalse(target.exists())

    def test_broad_mapping_does_not_become_exact_and_provenance_preserved(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'evidence.sqlite';c=sqlite3.connect(path)
            c.executescript('CREATE TABLE gene_disease(hgnc_id,symbol,orpha_id,disease); CREATE TABLE disease_xref(orpha_id,target_id,relation,validation); CREATE TABLE cspec(hgnc_id,mondo_id,payload); CREATE TABLE provenance(resource_id,payload);')
            c.execute('INSERT INTO gene_disease VALUES(?,?,?,?)',('HGNC:7551','MYBPC3','ORPHA:1','HCM'))
            c.execute('INSERT INTO disease_xref VALUES(?,?,?,?)',('ORPHA:1','MONDO:0005045','NTBT (narrower)','Validated'))
            c.execute('INSERT INTO provenance VALUES(?,?)',('orphadata',json.dumps(dict(sha256='a'*64,source=dict(release='July 2026')))))
            c.commit();c.close()
            result=m.gene_evidence(path,'HGNC:7551','MONDO:0005045')
            self.assertEqual(result['status'],'MISSING_EXACT_GENE_DISEASE_CONTEXT')
            self.assertEqual(result['provenance']['orphadata']['source']['release'],'July 2026')
            self.assertIsNone(result['pathogenicity_classification'])


if __name__=='__main__':unittest.main()
