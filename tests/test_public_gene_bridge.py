import csv
import importlib.util
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location('public_bridge',SCRIPTS/'export_public_gene_bridge.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class BridgeTests(unittest.TestCase):
    def test_exact_hgnc_compact_rows_preserve_specs_without_embedding_rules(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);c=sqlite3.connect(':memory:')
            c.executescript('CREATE TABLE gene_disease(hgnc_id,symbol,orpha_id,disease,association_type,assessment,validation);CREATE TABLE disease_xref(orpha_id,target_id,relation,validation);CREATE TABLE disease_hpo(orpha_id,hpo_id);CREATE TABLE cspec(hgnc_id,mondo_id,payload);')
            c.execute('INSERT INTO gene_disease VALUES(?,?,?,?,?,?,?)',('HGNC:1','PUBLIC','ORPHA:1','Public disease','Disease-causing','Assessed','PMID:1'))
            c.execute('INSERT INTO disease_xref VALUES(?,?,?,?)',('ORPHA:1','MONDO:0005045','E (Exact mapping)','Validated'))
            c.execute('INSERT INTO disease_hpo VALUES(?,?)',('ORPHA:1','HP:0001639'))
            payload=dict(spec_id='GN001',version='1.0',approved_on='2026-01-01',rules=[{'text':'FULL_RULE_NOT_IN_VARIANT_ROW'}],rule_sets=[])
            c.execute('INSERT INTO cspec VALUES(?,?,?)',('HGNC:1','MONDO:0005045',json.dumps(payload)))
            (root/'GN001-1.0.json').write_text(json.dumps(payload))
            rows,index=m.rows_for_gene(c,dict(hgnc_id='HGNC:1',symbol='PUBLIC'),root,'a'*64)
            self.assertEqual(rows[0]['exact_hcm_orpha_ids'],'ORPHA:1');self.assertEqual(rows[1]['cspec_ids_versions'],'GN001:1.0')
            self.assertNotIn('FULL_RULE_NOT_IN_VARIANT_ROW',json.dumps(rows));self.assertEqual(rows[0]['hpo_term_count'],'1')
            missing,_=m.rows_for_gene(c,dict(hgnc_id='HGNC:2',symbol='PUBLIC'),root,'a'*64)
            self.assertTrue(missing[0]['assertion_status'].startswith('MISSING'))

    def test_existing_runtime_gene_adapter_reads_new_context(self):
        annotation=SCRIPTS/'annotation_stages.py'
        if not annotation.exists():self.skipTest('Full existing runtime supplied separately for server integration control')
        spec=importlib.util.spec_from_file_location('existing_annotation',annotation);a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);tsv=root/'bridge.tsv'
            with tsv.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=m.FIELDS,delimiter='\t');writer.writeheader()
                for source in ['OrphadataScience','ClinGenCSpec']:
                    writer.writerow(dict(gene='PUBLIC',hgnc_id='HGNC:1',mondo_id='MONDO:0005045',source_kind=source,source_release='synthetic',assertion_status='SOURCE_CONTEXT_ONLY',cspec_ids_versions='GN001:1.0' if source=='ClinGenCSpec' else ''))
            item=dict(path=str(tsv),local_adapter=True,status='available',version='synthetic_bridge_v1')
            genes=a.gene_evidence(dict(resources={'gene_evidence':item}),root/'manifest.json')
            self.assertEqual(len(genes['PUBLIC']),2)
            db=sqlite3.connect(':memory:');db.execute('CREATE TABLE evidence(variant,transcript,evidence)')
            row=a.annotation_row(dict(chrom='chr1',pos='1',ref='A',alt='G',filter='PASS',qual='100'),
                dict(SYMBOL='PUBLIC',HGNC_ID='1',Feature='ENST_PUBLIC',Consequence='missense_variant'),
                'PUBLIC_SYNTHETIC_CONTROL',dict(GT='0/1',FT='PASS',DP='30',AD='15,15',GQ='60'),
                ({},dict(resources={'gene_evidence':item}),genes,{},db))
            context=json.loads(row['gene_context_json']);self.assertEqual({r['evidence']['source_kind'] for r in context},{'OrphadataScience','ClinGenCSpec'})
            self.assertEqual(row['gene_validity_rank'],0);self.assertIn('no_automatic',row['mechanism_status'])


if __name__=='__main__':unittest.main()
