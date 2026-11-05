import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('joined_views',ROOT/'scripts/export_cohort_sample_views_v2.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
sys.path.insert(0,str(Path(__file__).resolve().parent))
import test_cohort_sample_views as baseline


class Guard:
    def check(self):pass


def dump(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value));return value
def p(path):return m.pin(path,Guard())


def joined(workflow,*,alter=None,extra='allref_public_gene_context'):
    base=workflow/'stages/07_report/master_evidence.tsv.gz'
    for stage in ('06_annotate','07_report'):
        path=workflow/'stages'/stage/'workflow-receipt.json';value=json.loads(path.read_text());value['identity']='1'*64;dump(path,value)
    with gzip.open(base,'rt',newline='') as stream:
        values=list(csv.reader(stream,delimiter='\t'));fields,rows=values[0],values[1:]
    tissue=workflow.parent/'tissue-protein-context-v1/manifest.json'
    population=workflow.parent/'population-reference-context-v1/workflow-receipt.json'
    dump(tissue,dict(status='PASS_SEPARATE_RESEARCH_REFERENCE_CONTEXT_EXPORT',original_pipeline_changed=False,network_requests=0))
    dump(population,dict(status='PASS',stage='population-reference-context-v1',metadata=dict(status='PASS_LOCAL_POPULATION_MT_CONTEXT_EXPORT')))
    folder=workflow.parent/'combined-annotation-master-v1';folder.mkdir(exist_ok=True)
    enriched=folder/'master_evidence.enriched.tsv.gz'
    with gzip.open(enriched,'wt',newline='') as stream:
        writer=csv.writer(stream,delimiter='\t',lineterminator='\n');writer.writerow(fields+[extra])
        for number,row in enumerate(rows):
            source=list(row)
            if alter and number==0:source=alter(source)
            writer.writerow(source+['reference-context-only'])
    digest=io.StringIO(newline='');writer=csv.writer(digest,delimiter='\t',lineterminator='\n');writer.writerow(fields);writer.writerows(rows)
    manifest=dict(schema_version=1,status='PASS_COMBINED_ANNOTATION_MASTER',base_master=p(base),enriched_master=p(enriched),
        original_fields=fields,appended_fields=[extra],row_count=len(rows),invariant_digest=hashlib.sha256(digest.getvalue().encode()).hexdigest(),
        original_rows_fields_order_preserved=True,ranking_eligibility_unchanged=True,
        stage_receipts=[p(workflow/'stages'/s/'workflow-receipt.json') for s in ('06_annotate','07_report')],upstream_receipts=[p(tissue),p(population)],created_utc='2026-11-05T20:00:00Z')
    path=folder/'manifest.json';dump(path,manifest)
    dump(folder/'workflow-receipt.json',dict(status='PASS',stage='combined-annotation-master-v1',identity='2'*64,
        artifacts={file.name:{k:v for k,v in p(file).items() if k!='path'} for file in (path,enriched)}))
    return path


def applied(path):
    folder=path.parent;population_path=folder.parent/'population-reference-context-v1/workflow-receipt.json'
    population=json.loads(population_path.read_text());source_relative='cohort/nuclear-population-long.tsv.gz'
    source=population_path.parent/source_relative;source.parent.mkdir(parents=True,exist_ok=True);source.write_bytes(b'synthetic applied result table')
    expected=p(source);population['artifacts']={source_relative:{k:expected[k] for k in ('bytes','sha256')}};dump(population_path,population)
    relative='applied-evidence/population/'+source_relative;target=folder/relative;target.parent.mkdir(parents=True);target.write_bytes(source.read_bytes())
    dump(folder/'auxiliary-index.json',dict(scope='applied result table links'));dump(folder/'applied-evidence/MT-appendix-index.json',dict(scope='separate MT appendix'))
    value=json.loads(path.read_text());value['upstream_receipts']=[p(Path(item['path'])) for item in value['upstream_receipts']]
    value['auxiliary_artifacts']=[dict(role='APPLIED_NUCLEAR_REFERENCE_LONG_TABLE',path=str(target),relative_path=relative,bytes=expected['bytes'],sha256=expected['sha256'],source_relative_path=source_relative,source_upstream_receipt=p(population_path))]
    dump(path,value);receipt=json.loads((folder/'workflow-receipt.json').read_text())
    for file in [path,target,folder/'auxiliary-index.json',folder/'applied-evidence/MT-appendix-index.json']:
        receipt['artifacts'][file.relative_to(folder).as_posix()]={k:v for k,v in p(file).items() if k!='path'}
    dump(folder/'workflow-receipt.json',receipt);return relative


class JoinedInvariantTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.workflow=Path(self.temp.name).resolve()/'workflow';stage7=self.workflow/'stages/07_report';stage7.mkdir(parents=True)
        stage6=self.workflow/'stages/06_annotate';stage6.mkdir()
        base=stage7/'master_evidence.tsv.gz'
        with gzip.open(base,'wt',newline='') as stream:
            writer=csv.writer(stream,delimiter='\t',lineterminator='\n');writer.writerow(['sample','GT','FT','candidate_eligible','raw_context'])
            writer.writerow(['SYNTHETIC_A','1/.','LowDP','0','literal\ttab\nand newline'])
        dump(stage7/'workflow-receipt.json',dict(status='PASS',identity='1'*64,artifacts={base.name:{k:v for k,v in p(base).items() if k!='path'}}))
        dump(stage6/'workflow-receipt.json',dict(status='PASS',identity='1'*64,artifacts={}))

    def test_quote_aware_original_invariants_and_namespace(self):
        path=joined(self.workflow);value,pins=m.verify_joined(path,self.workflow,Guard())
        self.assertEqual(value['row_count'],1);self.assertEqual(value['appended_fields'],['allref_public_gene_context'])
        self.assertTrue(any(item['path'].endswith('master_evidence.enriched.tsv.gz') for item in pins))
        path=joined(self.workflow,extra='GT_overwrite')
        with self.assertRaisesRegex(ValueError,'namespaced'):m.verify_joined(path,self.workflow,Guard())

    def test_true_bools_and_fresh_hashes_cannot_hide_mutated_original_gt(self):
        def change(row):row[1]='0/0';return row
        path=joined(self.workflow,alter=change)
        with self.assertRaisesRegex(ValueError,'changed an original'):m.verify_joined(path,self.workflow,Guard())
        path=joined(self.workflow);value=json.loads(path.read_text());value['upstream_receipts'][0]['sha256']='0'*64;dump(path,value)
        with self.assertRaisesRegex(ValueError,'pin differs'):m.verify_joined(path,self.workflow,Guard())


@unittest.skipUnless(baseline.pysam and sys.platform=='linux','Linux pysam required for full atomic joined views proof')
class JoinedViewIntegrationTests(unittest.TestCase):
    def test_joined_columns_all_outputs_preserve_sealed_rank_and_original_rows(self):
        fixture=baseline.ViewIntegrationTests();fixture.setUp();self.addCleanup(fixture.doCleanups)
        path=joined(fixture.workflow);applied_relative=applied(path)
        base_exporter=ROOT/'scripts/export_cohort_sample_views.py';ranking=ROOT/'scripts/annotation_stages.py'
        base_result=baseline.m.export(fixture.workflow,fixture.root/'baseline-output',ranking,baseline.m.sha(ranking),synthetic=True)
        output=fixture.root/'joined-output'
        result=m.export(fixture.workflow,output,ranking,baseline.m.sha(ranking),base_exporter,baseline.m.sha(base_exporter),path,synthetic=True)
        self.assertEqual(result['counts'],base_result['counts']);self.assertEqual(result['report_input_kind'],'COMBINED_ANNOTATION_MASTER_V1')
        for filename in ['cohort.top20.tsv','sample_0001.top20.tsv']:
            before=fixture.read_table(fixture.root/'baseline-output/tables'/filename);after=fixture.read_table(output/'tables'/filename)
            self.assertEqual([(r['rank'],r['variant_id'],r['sample'],r['transcript']) for r in before],[(r['rank'],r['variant_id'],r['sample'],r['transcript']) for r in after])
            self.assertTrue(all(r['allref_public_gene_context']=='reference-context-only' for r in after))
        for filename in ['cohort.master.tsv.gz','sample_0001.master.tsv.gz']:
            self.assertTrue(all(r['allref_public_gene_context']=='reference-context-only' for r in fixture.read_table(output/'masters'/filename)))
        self.assertTrue(all(r['allref_public_gene_context']=='reference-context-only' for r in fixture.read_table(output/'tables/sample_0001.functional.tsv.gz')))
        self.assertTrue((output/applied_relative).is_file())
        self.assertEqual(p(output/applied_relative)['sha256'],p(path.parent/applied_relative)['sha256'])
        self.assertIn(applied_relative,(output/'index.html').read_text())
        self.assertIn(applied_relative,result['copied_applied_evidence'])
        bad=json.loads(path.read_text());bad['auxiliary_artifacts'][0]['relative_path']='resources/database.tsv.gz'
        with self.assertRaisesRegex(ValueError,'Only applied'):
            m.copy_applied_evidence(bad,fixture.root/'must-not-copy-resources',Guard())
        self.assertFalse(any(fixture.root.glob('.joined-output.joined-partial-*')))


if __name__=='__main__':unittest.main()
