import base64
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('native_plan',ROOT/'transfer/prepare_native_plan.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class NativePlanSidecarTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.base=Path(self.temp.name).resolve()
        self.root=self.base/'srv/marina';self.root.mkdir(parents=True)
        self.jroot=self.base/'jupyter';self.workspace=self.jroot/'workspace';self.workspace.mkdir(parents=True)
        self.stage6=self.root/'outputs/workflow-v1/stages/06_annotate';self.stage7=self.root/'outputs/workflow-v1/stages/07_report'
        for stage,names in [(self.stage6,['annotated.vcf.gz','annotated.vcf.gz.tbi','transcript_samples.tsv.gz']),
                            (self.stage7,['master_evidence.tsv.gz','table.tsv.gz','table.csv.gz','plots/overview.png','plots/overview.jpg','plots/overview.jpeg','plots/overview.svg'])]:
            stage.mkdir(parents=True)
            artifacts={name:self.file(stage/name,b'synthetic completed table') for name in names}
            self.json(stage/'workflow-receipt.json',dict(status='PASS',identity='1'*64,artifacts=artifacts))
        for i in range(4):
            stage=self.root/'outputs/remaining-align'/('SYNTHETIC_%d'%i)/'align';stage.mkdir(parents=True)
            # Expected historical BAM SHA deliberately differs from file bytes:
            # the planner must retain expected identity and never rehash BAMs.
            artifacts={name:self.file(stage/name,b'synthetic alignment') for name in ['marked.bam','marked.bam.bai']}
            artifacts['marked.bam']['sha256']='a'*64
            self.json(stage/'receipt.json',dict(completed_utc='2026-11-05T20:00:00Z',fingerprint='legacy-complete',outputs=artifacts))
        stage=self.root/'outputs/workflow-v1/checkpoints/03_call';stage.mkdir(parents=True)
        artifacts={name:self.file(stage/name,b'gvcf stays server-side') for name in ['SYNTHETIC.g.vcf.gz','SYNTHETIC.g.vcf.gz.tbi']}
        self.json(stage/'workflow-receipt.json',dict(status='PASS',identity='2'*64,artifacts=artifacts))
        self.sidecars()

    def file(self,path,raw):
        path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
        return dict(bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())

    def json(self,path,value):
        path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value));return value

    def pin(self,path):
        raw=path.read_bytes();return dict(path=str(path),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())

    def primary_pins(self):
        return {str(stage/'workflow-receipt.json'):dict(sha256=self.pin(stage/'workflow-receipt.json')['sha256']) for stage in [self.stage6,self.stage7]}

    def sidecars(self):
        source=self.pin(self.stage6/'transcript_samples.tsv.gz');master=self.pin(self.stage7/'master_evidence.tsv.gz')
        consequent=self.root/'outputs/consequence-tables-v1'
        producer=self.root/'private/reviewed-consequence-export-v1/scripts/export_consequence_tables.py';self.file(producer,b'reviewed consequence producer')
        artifacts={name:self.file(consequent/name,b'synthetic consequence result') for name in ['index.html','index.tsv','tables/cohort/stop_gained.tsv.gz','evidence/transcript_effects.tsv.gz']}
        self.json(consequent/'manifest.json',dict(schema_version=1,status='PASS',synthetic=False,exporter_sha256=self.pin(producer)['sha256'],source=source,
            workflow_receipt=json.loads((self.stage6/'workflow-receipt.json').read_text()),artifacts=artifacts))
        tissue=self.root/'outputs/tissue-protein-context-v1';target=tissue/'public39-reference-context.tsv.gz';artifact=self.file(target,b'synthetic tissue result')
        context=self.root/'resources/tissue-protein-context-public39-v1/context-index.json';self.file(context,b'synthetic public gene context never selected for PC')
        self.json(tissue/'manifest.json',dict(status='PASS_SEPARATE_RESEARCH_REFERENCE_CONTEXT_EXPORT',source_artifact=source['path'],source_sha256=source['sha256'],
            source_bytes=source['bytes'],context_index_sha256=self.pin(context)['sha256'],stage_receipts=self.primary_pins(),counters={},artifact=dict(artifact,path=str(target)),
            scope='public gene context',original_pipeline_changed=False,network_requests=0,created_utc='2026-11-05T20:00:00Z'))
        cohort=self.root/'outputs/cohort-sample-views-v1';ranking=self.root/'private/annotation-runtime-v1/scripts/annotation_stages.py';self.file(ranking,b'sealed ranking')
        artifacts={name:self.file(cohort/name,b'synthetic view') for name in ['index.html','masters/cohort.master.tsv.gz','tables/cohort.top20.tsv','tables/cohort.functional.tsv.gz']+
            ['vcf/sample_%04d.annotated.vcf.gz%s'%(i,suffix) for i in range(1,5) for suffix in ('','.tbi')]}
        pins={str(p):{k:v for k,v in self.pin(p).items() if k!='path'} for p in [self.stage6/'workflow-receipt.json',self.stage7/'workflow-receipt.json',self.stage6/'annotated.vcf.gz',self.stage6/'annotated.vcf.gz.tbi',self.stage7/'master_evidence.tsv.gz']}
        samples=['SYNTHETIC_%d'%i for i in range(4)]
        self.json(cohort/'manifest.json',dict(schema_version=1,status='PASS',synthetic=False,primary_workflow_changed=False,no_external_requests=True,automatic_ACMG=False,
            ranking_module=dict(path=str(ranking),sha256=self.pin(ranking)['sha256'],function='rank_key'),source_pins=pins,samples=samples,
            sample_file_labels={s:'sample_%04d'%(i+1) for i,s in enumerate(samples)},artifacts=artifacts,completed_utc='2026-11-05T20:00:00Z'))
        population=self.root/'outputs/population-reference-context-v1'
        artifacts={name:self.file(population/name,b'synthetic population result') for name in ['README.md','provenance.json','cohort/nuclear-population-long.tsv.gz','samples/SYNTHETIC_0/MT-population-reference.tsv.gz']}
        meta=dict(status='PASS_LOCAL_POPULATION_MT_CONTEXT_EXPORT',source_master=master,primary_receipts=[self.pin(p/'workflow-receipt.json') for p in [self.stage6,self.stage7]],
            network_requests=0,primary_outputs_changed=False,automatic_ACMG=False)
        self.json(population/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='population-reference-context-v1',identity='4'*64,
            completed_utc='2026-11-05T20:00:00Z',metadata=meta,artifacts=artifacts))

    def inventory(self):
        cfg=dict(root=str(self.root),workspace=str(self.workspace),jupyter=str(self.jroot),compare_only=True)
        old=sys.argv;out=io.StringIO();environment={}
        try:
            sys.argv=['helper',json.dumps(cfg)]
            with contextlib.redirect_stdout(out):exec(compile(m.HELPER,'<native-plan-helper>','exec'),environment)
        finally:sys.argv=old
        return json.loads(out.getvalue()),environment

    def test_all_four_completed_sidecars_compressed_tables_and_four_sample_vcfs(self):
        plan,env=self.inventory();paths={x['RelativePath'] for x in plan['Files']}
        for name in ['consequence-tables-v1','tissue-protein-context-v1','population-reference-context-v1','cohort-sample-views-v1']:
            self.assertTrue(any(p.startswith('outputs/'+name+'/') for p in paths),name)
        self.assertEqual(sum('/cohort-sample-views-v1/vcf/' in p for p in paths),8)
        self.assertEqual(sum(p.endswith('.bam') for p in paths),4)
        self.assertFalse(any('.g.vcf' in p for p in paths))
        for suffix in ['table.tsv.gz','table.csv.gz','plots/overview.png','plots/overview.jpg','plots/overview.jpeg','plots/overview.svg']:
            entry=next(x for x in env['files'].values() if x['RelativePath'].endswith(suffix))
            self.assertTrue(env['view_rows'][entry['RelativePath']]['ViewPath'].startswith('Reports/'))
            if suffix in ('table.tsv.gz','table.csv.gz'):self.assertNotIn(entry['RelativePath'],paths)
        self.assertTrue(all(x['Sha256']=='a'*64 for x in plan['Files'] if x['RelativePath'].endswith('.bam')))
        self.assertEqual(plan['Mode'],'COMPARE_ONLY_NO_LINKS_WRITTEN');self.assertEqual(list(self.workspace.iterdir()),[])
        for receipt in plan['Receipts']:
            raw=base64.b64decode(receipt['ContentBase64']);self.assertEqual(hashlib.sha256(raw).hexdigest(),receipt['Sha256'])

    def test_failed_partial_unknown_dataset_and_resource_trees_not_selected(self):
        for name in ['unrecognized-dataset','resources','databases','.population-reference-context-v1.partial-test']:
            folder=self.root/'outputs'/name;artifacts={'db.tsv.gz':self.file(folder/'db.tsv.gz',b'never copy dataset')}
            self.json(folder/'workflow-receipt.json',dict(status='PASS',identity='9'*64,artifacts=artifacts))
        path=self.root/'outputs/tissue-protein-context-v1/manifest.json';value=json.loads(path.read_text());value['status']='FAILED';self.json(path,value)
        plan,_=self.inventory();paths=[x['RelativePath'] for x in plan['Files']]
        self.assertFalse(any('tissue-protein-context-v1' in p for p in paths))
        self.assertFalse(any('db.tsv.gz' in p or '.partial-' in p for p in paths))

    def test_sidecar_bad_path_and_source_identity_rejected(self):
        path=self.root/'outputs/consequence-tables-v1/manifest.json';value=json.loads(path.read_text())
        value['artifacts']['../outside.tsv.gz']=dict(bytes=1,sha256='a'*64);self.json(path,value)
        with self.assertRaises(RuntimeError):self.inventory()
        del value['artifacts']['../outside.tsv.gz'];value['source']['sha256']='0'*64;self.json(path,value)
        with self.assertRaisesRegex(RuntimeError,'source pin'):self.inventory()
        self.assertEqual(list(self.workspace.iterdir()),[])

    def test_recognized_sidecar_cannot_smuggle_resources_or_changed_size(self):
        path=self.root/'outputs/cohort-sample-views-v1/manifest.json';value=json.loads(path.read_text())
        value['artifacts']['resources/reference.tsv.gz']=self.file(path.parent/'resources/reference.tsv.gz',b'public database')
        self.json(path,value)
        with self.assertRaisesRegex(RuntimeError,'resource path'):self.inventory()
        del value['artifacts']['resources/reference.tsv.gz'];self.json(path,value)
        (path.parent/'tables/cohort.functional.tsv.gz').write_bytes(b'changed')
        with self.assertRaisesRegex(RuntimeError,'changed size'):self.inventory()

    def test_producer_and_context_identity_not_just_pass_anywhere(self):
        path=self.root/'outputs/consequence-tables-v1/manifest.json';value=json.loads(path.read_text());value['exporter_sha256']='0'*64;self.json(path,value)
        with self.assertRaisesRegex(RuntimeError,'producer identity'):self.inventory()
        value['exporter_sha256']=self.pin(self.root/'private/reviewed-consequence-export-v1/scripts/export_consequence_tables.py')['sha256'];self.json(path,value)
        context=self.root/'resources/tissue-protein-context-public39-v1/context-index.json';context.write_bytes(b'changed public context')
        with self.assertRaisesRegex(RuntimeError,'context identity'):self.inventory()

    def test_sidecar_symlink_cannot_escape_fixed_output_root(self):
        path=self.root/'outputs/consequence-tables-v1/tables/cohort/stop_gained.tsv.gz'
        outside=self.base/'outside.tsv.gz';outside.write_bytes(path.read_bytes());path.unlink()
        try:path.symlink_to(outside)
        except OSError:self.skipTest('Symlinks unavailable in this local test environment')
        with self.assertRaisesRegex(RuntimeError,'Symlink'):self.inventory()

    def combined(self):
        folder=self.root/'outputs/combined-annotation-master-v1';enriched=folder/'master_evidence.enriched.tsv.gz';self.file(enriched,b'synthetic enriched patient-result fixture')
        upstream=[self.pin(self.root/'outputs/tissue-protein-context-v1/manifest.json'),self.pin(self.root/'outputs/population-reference-context-v1/workflow-receipt.json')]
        manifest=dict(schema_version=1,status='PASS_COMBINED_ANNOTATION_MASTER',base_master=self.pin(self.stage7/'master_evidence.tsv.gz'),enriched_master=self.pin(enriched),
            original_fields=['sample','GT','FT'],appended_fields=['allref_context'],row_count=1,invariant_digest='5'*64,
            original_rows_fields_order_preserved=True,ranking_eligibility_unchanged=True,stage_receipts=[self.pin(p/'workflow-receipt.json') for p in [self.stage6,self.stage7]],
            upstream_receipts=upstream,created_utc='2026-11-05T20:00:00Z')
        self.json(folder/'manifest.json',manifest)
        self.json(folder/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='combined-annotation-master-v1',identity='6'*64,
            artifacts={p.name:{k:v for k,v in self.pin(p).items() if k!='path'} for p in [folder/'manifest.json',enriched]}))
        path=self.root/'outputs/cohort-sample-views-v1/manifest.json';views=json.loads(path.read_text())
        engine=self.root/'private/reviewed-cohort-sample-views-v1/scripts/export_cohort_sample_views.py';producer=self.root/'private/reviewed-cohort-sample-views-v2/scripts/export_cohort_sample_views_v2.py'
        self.file(engine,b'frozen views engine');self.file(producer,b'new reviewed joined views producer')
        views.update(report_input_kind='COMBINED_ANNOTATION_MASTER_V1',views_v1_engine=dict(path=str(engine),sha256=self.pin(engine)['sha256']),views_v2_exporter_sha256=self.pin(producer)['sha256'],
            joined_master=dict(manifest=self.pin(folder/'manifest.json'),**{k:manifest[k] for k in ['enriched_master','base_master','original_fields','appended_fields','row_count','invariant_digest']},
                original_rows_fields_order_verified=True,rank_eligibility_columns_unchanged_verified=True))
        for pin in [self.pin(folder/'manifest.json'),self.pin(folder/'workflow-receipt.json'),self.pin(enriched),*upstream]:views['source_pins'][pin['path']]={k:pin[k] for k in ('bytes','sha256')}
        self.json(path,views);return folder,path

    def test_combined_master_and_joined_views_provenance_recognized_without_copying_resources(self):
        folder,path=self.combined();plan,env=self.inventory();paths={x['RelativePath'] for x in plan['Files']}
        self.assertNotIn('outputs/combined-annotation-master-v1/master_evidence.enriched.tsv.gz',paths)
        self.assertIn('outputs/combined-annotation-master-v1/master_evidence.enriched.tsv.gz',env['files'])
        self.assertIn('outputs/combined-annotation-master-v1/manifest.json',paths)
        self.assertEqual(sum('/cohort-sample-views-v1/vcf/' in p for p in paths),8)
        self.assertFalse(any(p.startswith('resources/') or '/resources/' in p for p in paths))
        views=json.loads(path.read_text());views['joined_master']['rank_eligibility_columns_unchanged_verified']=False;self.json(path,views)
        with self.assertRaisesRegex(RuntimeError,'invariants'):self.inventory()

    def test_joined_views_unknown_source_pin_or_bad_receipt_fails_closed(self):
        folder,path=self.combined();views=json.loads(path.read_text());views['source_pins'][str(self.root/'resources/database.tsv.gz')]=dict(bytes=1,sha256='0'*64);self.json(path,views)
        with self.assertRaisesRegex(RuntimeError,'Unknown cohort source pin'):self.inventory()
        del views['source_pins'][str(self.root/'resources/database.tsv.gz')];self.json(path,views)
        receipt=json.loads((folder/'workflow-receipt.json').read_text());receipt['artifacts']['manifest.json']['sha256']='0'*64;self.json(folder/'workflow-receipt.json',receipt)
        with self.assertRaisesRegex(RuntimeError,'receipt artifact identity'):self.inventory()

    def test_applied_result_evidence_has_dual_receipts_and_remains_reports(self):
        folder,path=self.combined();manifest=json.loads((folder/'manifest.json').read_text())
        population=self.root/'outputs/population-reference-context-v1';source_relative='cohort/nuclear-population-long.tsv.gz'
        source=population/source_relative;relative='applied-evidence/population/'+source_relative
        expected=self.pin(source);self.file(folder/relative,source.read_bytes())
        indices=['auxiliary-index.json','applied-evidence/MT-appendix-index.json']
        for name in indices:self.json(folder/name,dict(scope='completed applied result index'))
        entry=dict(role='APPLIED_NUCLEAR_REFERENCE_LONG_TABLE',path=str(folder/relative),relative_path=relative,
            source_relative_path=source_relative,source_upstream_receipt=self.pin(population/'workflow-receipt.json'),bytes=expected['bytes'],sha256=expected['sha256'])
        manifest['auxiliary_artifacts']=[entry];self.json(folder/'manifest.json',manifest)
        receipt=json.loads((folder/'workflow-receipt.json').read_text())
        for name in ['manifest.json',relative,*indices]:receipt['artifacts'][name]={k:self.pin(folder/name)[k] for k in ('bytes','sha256')}
        self.json(folder/'workflow-receipt.json',receipt)
        views=json.loads(path.read_text());views['joined_master']['manifest']=self.pin(folder/'manifest.json')
        for item in [self.pin(folder/'manifest.json'),self.pin(folder/'workflow-receipt.json')]:views['source_pins'][item['path']]={k:item[k] for k in ('bytes','sha256')}
        views['copied_applied_evidence']=[relative,*indices]
        for name in [relative,*indices]:views['artifacts'][name]=self.file(path.parent/name,(folder/name).read_bytes())
        self.json(path,views);plan,env=self.inventory()
        for branch in ['combined-annotation-master-v1','cohort-sample-views-v1']:
            source_path='outputs/'+branch+'/'+relative
            if branch=='cohort-sample-views-v1':self.assertIn(source_path,{x['RelativePath'] for x in plan['Files']})
            else:self.assertNotIn(source_path,{x['RelativePath'] for x in plan['Files']})
            self.assertTrue(env['view_rows'][source_path]['ViewPath'].startswith('Reports/'))
        manifest['auxiliary_artifacts'][0]['source_upstream_receipt']['sha256']='0'*64;self.json(folder/'manifest.json',manifest)
        with self.assertRaises(RuntimeError):self.inventory()

    def vcf_first(self):
        joined,views_path=self.combined();combined=json.loads((joined/'manifest.json').read_text());folder=self.root/'outputs/fully-annotated-vcf-v1'
        vcf=folder/'cohort.fully-annotated.vcf.gz';index=folder/'cohort.fully-annotated.vcf.gz.tbi';self.file(vcf,b'synthetic self contained evidence VCF');self.file(index,b'synthetic index')
        codec=self.root/'private/fully-annotated-vcf-candidate-v1/scripts/embed_vcf_evidence.py';self.file(codec,b'reviewed packet codec')
        schema=dict(v=1,codec='json-utf8+zlib+base64url-nopad',scope='cohort',sample=None,**{k:combined[k] for k in ['original_fields','appended_fields','row_count','invariant_digest']},
            full_fields=combined['original_fields']+combined['appended_fields'],full_fields_digest='7'*64)
        fully=dict(schema_version=1,status='PASS_FULLY_ANNOTATED_VCF',network_requests=0,primary_outputs_changed=False,
            original_records_fields_preserved=True,packet_decode_verified=True,decoded_master_vectors_verified=True,public_header_bundle_verified=True,
            base_vcf=self.pin(self.stage6/'annotated.vcf.gz'),master_manifest=self.pin(joined/'manifest.json'),combined_receipt=self.pin(joined/'workflow-receipt.json'),
            master=combined['enriched_master'],enriched_vcf=self.pin(vcf),index=self.pin(index),schema=schema,producer_code_sha256={str(codec):self.pin(codec)['sha256']},
            **{k:schema[k] for k in ['original_fields','appended_fields','full_fields','row_count','invariant_digest','full_fields_digest']})
        self.json(folder/'manifest.json',fully)
        self.json(folder/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='fully-annotated-vcf-v1',identity='8'*64,
            artifacts={p.name:{k:self.pin(p)[k] for k in ('bytes','sha256')} for p in [vcf,index,folder/'manifest.json']}))
        views=json.loads(views_path.read_text());extracted=self.root/'private/cohort-sample-views-continuation-v3'/('extraction-'+'a'*32)/'master.from-vcf.tsv.gz'
        self.file(extracted,b'synthetic lossless full row master');views['artifacts']['masters/cohort.master.tsv.gz']=self.file(views_path.parent/'masters/cohort.master.tsv.gz',extracted.read_bytes())
        helper=self.root/'private/reviewed-cohort-sample-views-v2-1/scripts/export_cohort_sample_views_v2.py';producer=self.root/'private/reviewed-cohort-sample-views-v3/scripts/export_cohort_sample_views_v3.py'
        self.file(helper,b'frozen control and applied evidence helper');self.file(producer,b'VCF first view producer')
        views.update(report_input_kind='FULLY_ANNOTATED_VCF_PACKETS_V1',fully_annotated_vcf=dict(manifest=self.pin(folder/'manifest.json'),vcf=fully['enriched_vcf'],index=fully['index'],schema=schema),
            joined_control=dict(manifest=fully['master_manifest'],master=fully['master'],all_full_field_vectors_order_verified=True,used_as_report_input=False),
            extracted_master=dict(self.pin(extracted),row_count=fully['row_count'],full_fields_digest=fully['full_fields_digest'],original_fields_digest=fully['invariant_digest'],source='VCF_INFO_WGS_EVIDENCE_ONLY'),
            sample_packet_scope_verified=True,public_header_bundles_verified=True,sample_packet_checks=[dict(sample=s,scope='sample',other_sample_rows=0,original_source_ordinals_retained=True) for s in views['samples']],
            views_v2_control_helper=dict(path=str(helper),sha256=self.pin(helper)['sha256']),packet_codec=dict(path=str(codec),sha256=self.pin(codec)['sha256']),views_v3_exporter_sha256=self.pin(producer)['sha256'],copied_applied_evidence=[])
        for source in [folder/'manifest.json',folder/'workflow-receipt.json',vcf,index,extracted]:views['source_pins'][str(source)]={k:self.pin(source)[k] for k in ('bytes','sha256')}
        for source,name in [(vcf,'vcf/cohort.fully-annotated.vcf.gz'),(index,'vcf/cohort.fully-annotated.vcf.gz.tbi'),
            (folder/'manifest.json','provenance/fully-annotated-vcf-manifest.json'),(folder/'workflow-receipt.json','provenance/fully-annotated-vcf-receipt.json'),
            (joined/'manifest.json','provenance/joined-control-manifest.json'),(joined/'workflow-receipt.json','provenance/joined-control-receipt.json')]:
            views['artifacts'][name]=self.file(views_path.parent/name,source.read_bytes())
        self.json(views_path,views);return folder,views_path

    def test_vcf_first_final_selection_one_cohort_copy_and_sample_scope_proof(self):
        folder,path=self.vcf_first();plan,env=self.inventory();paths={entry['RelativePath'] for entry in plan['Files']}
        self.assertIn('outputs/cohort-sample-views-v1/vcf/cohort.fully-annotated.vcf.gz',paths)
        self.assertNotIn('outputs/fully-annotated-vcf-v1/cohort.fully-annotated.vcf.gz',paths)
        self.assertIn('outputs/fully-annotated-vcf-v1/cohort.fully-annotated.vcf.gz',env['files'])
        self.assertEqual(sum('/cohort-sample-views-v1/vcf/' in p for p in paths),10)
        self.assertIn('outputs/cohort-sample-views-v1/masters/cohort.master.tsv.gz',paths)
        self.assertNotIn('outputs/workflow-v1/stages/07_report/master_evidence.tsv.gz',paths)
        self.assertEqual(sum(p.endswith('.bam') for p in paths),4)
        views=json.loads(path.read_text());views['sample_packet_checks'][0]['other_sample_rows']=1;self.json(path,views)
        with self.assertRaisesRegex(RuntimeError,'packet pruning'):self.inventory()


if __name__=='__main__':unittest.main()
