import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from types import SimpleNamespace

sys.path.insert(0,str(Path(__file__).parent))
from test_native_plan_sidecars import NativePlanSidecarTests as Harness
from test_native_plan_sidecars import m as planner


class DirectPlanTests(unittest.TestCase):
    setUp=Harness.setUp;file=Harness.file;json=Harness.json;pin=Harness.pin
    sidecars=Harness.sidecars;primary_pins=Harness.primary_pins;inventory=Harness.inventory

    def direct(self):
        for path in [self.root/'outputs/consequence-tables-v1/manifest.json',self.root/'outputs/tissue-protein-context-v1/manifest.json',
            self.root/'outputs/population-reference-context-v1/workflow-receipt.json',self.stage6/'workflow-receipt.json',self.stage7/'workflow-receipt.json']:
            value=json.loads(path.read_text());value['status']='FAILED';self.json(path,value)
        stage=self.root/'outputs/workflow-v1/stages/05_filter'
        self.json(stage/'workflow-receipt.json',dict(status='PASS',identity='5'*64,artifacts={'cohort.filtered.vcf.gz':self.file(stage/'cohort.filtered.vcf.gz',b'accepted05')}))
        needed={'binding','config','resources','snapshot_database','snapshot_manifest','completed_vep_vcf','completed_vep_index','vep_component_producer','filter_receipt','filter_vcf','alpha_manifest','alpha_canonical_matches','alpha_isoform_matches','core_public_proof','core_module'}
        sources={}
        for name in needed:
            path=self.root/'private/direct-lineage'/name;self.file(path,b'finite source');sources[name]=self.pin(path)
        sources['filter_receipt']=self.pin(stage/'workflow-receipt.json');sources['filter_vcf']=self.pin(stage/'cohort.filtered.vcf.gz')
        producer=self.root/'private/direct-evidence-candidate-v1/scripts/fast_vcf_evidence_candidate.py';self.file(producer,b'core source')
        sources['core_module']=self.pin(producer);producers={'core_module':self.pin(producer)['sha256'],'final_driver':'f'*64}
        def schema(fields):
            batch=dict(id=0,first_i=1,last_i=4,row_count=4,record_count=1,canonical_row_sha256='c'*64)
            header=io.StringIO(newline='');csv.writer(header,delimiter='\t',lineterminator='\n').writerow(fields)
            c=dict(contract='ordered_batch_canonical_TSV_utf8_LF_sha256_v1',canonical_header_sha256=hashlib.sha256(header.getvalue().encode()).hexdigest(),batches=[batch])
            return dict(full_fields=fields,row_count=4,batches=[batch],vector_digest_contract=c['contract'],ordinal_in_vcf_order=True,
                producer_kind='DIRECT_FROZEN_ANNOTATION_ROW_V1_NO_LEGACY_06_07_RECEIPTS',ordered_batch_merkle_root=hashlib.sha256(json.dumps(c,sort_keys=True,separators=(',',':')).encode()).hexdigest())
        folder=self.root/'outputs/fully-annotated-vcf-v1';core=folder/'core-manifest.json'
        self.json(core,dict(status='DIRECT_PACKET_BASELINE_PASS',no_legacy_receipt_created=True,schema=schema(['sample','GT']),source_pins=sources,producer_code_sha256=producers))
        for name in ['cohort.fully-annotated.vcf.gz','cohort.fully-annotated.vcf.gz.tbi']:self.file(folder/name,b'one final '+name.encode())
        full=dict(schema_version=1,status='PASS_FULLY_ANNOTATED_VCF_DIRECT',direct_core_prefix_verified=True,packet_decode_verified=True,original_records_fields_preserved=True,
            core_manifest=self.pin(core),schema=schema(['sample','GT','allref_context']),producer_code_sha256=producers,
            enriched_vcf=self.pin(folder/'cohort.fully-annotated.vcf.gz'),index=self.pin(folder/'cohort.fully-annotated.vcf.gz.tbi'))
        self.json(folder/'manifest.json',full)
        self.json(folder/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='fully-annotated-vcf-v1',identity='6'*64,
            artifacts={n:{k:v for k,v in self.pin(folder/n).items() if k!='path'} for n in ['manifest.json','cohort.fully-annotated.vcf.gz','cohort.fully-annotated.vcf.gz.tbi']}))
        views=self.root/'outputs/cohort-sample-views-v1';meta=json.loads((views/'manifest.json').read_text());artifacts=meta['artifacts']
        required={'tables/functional-category-index.tsv','tables/allele-membership-matrix.tsv.gz','tables/ACMG-manual-top20-worksheet.tsv.gz','tables/gene-overlap-SEPARATE-NOT-SHARED-ALLELES.tsv.gz',
            'provenance/fully-annotated-vcf-manifest.json','provenance/direct-core-manifest.json','tables/spliceai/cohort.top20.tsv'}
        for label in meta['sample_file_labels'].values():required.update({'masters/'+label+'.master.tsv.gz','tables/'+label+'.top20.tsv','tables/'+label+'.functional.tsv.gz','tables/spliceai/'+label+'.top20.tsv'})
        for name in required:artifacts[name]=self.file(views/name,b'final result')
        for name in ['cohort.fully-annotated.vcf.gz','cohort.fully-annotated.vcf.gz.tbi']:artifacts['vcf/'+name]=self.file(views/'vcf'/name,(folder/name).read_bytes())
        exporter=self.root/'private/reviewed-packet-fastviews-v1/scripts/export_packet_views_fast.py';self.file(exporter,b'reviewed fast source')
        meta.update(report_producer_kind='FAST_DIRECT_PACKET_VIEWS_V1',report_input_kind='FULLY_ANNOTATED_VCF_PACKETS_V1',sample_packet_scope_verified=True,
            public_header_bundles_verified=True,no_legacy_receipt_created=True,producer_code_sha256=self.pin(exporter)['sha256'],
            source_pins=dict(fully_manifest=self.pin(folder/'manifest.json'),core_manifest=self.pin(core),vcf=full['enriched_vcf']),
            batch_merkle_root=full['schema']['ordered_batch_merkle_root'],core_prefix_merkle_root=json.loads(core.read_text())['schema']['ordered_batch_merkle_root'])
        self.json(views/'manifest.json',meta);self.finish(views,meta);return views

    def finish(self,views,meta):
        self.json(views/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='cohort-sample-views-v1',producer_code_sha256=meta['producer_code_sha256'],
            artifacts=dict(meta['artifacts'],**{'manifest.json':{k:v for k,v in self.pin(views/'manifest.json').items() if k!='path'}})))

    def test_direct_one_canonical_vcf_four_samples_without_06_07(self):
        self.direct();plan,_=self.inventory();paths={x['RelativePath'] for x in plan['Files']}
        self.assertEqual(sum('/vcf/cohort.fully-annotated.vcf.gz' in p for p in paths),2)
        self.assertEqual(sum('/vcf/sample_' in p for p in paths),8)
        self.assertFalse(any(p.startswith('outputs/fully-annotated-vcf-v1/cohort.') or '.g.vcf' in p for p in paths))
        self.assertTrue(any('/tables/spliceai/' in p for p in paths))

    def test_resource_smuggling_rejected(self):
        views=self.direct();meta=json.loads((views/'manifest.json').read_text());meta['artifacts']['resources/db.tsv.gz']=self.file(views/'resources/db.tsv.gz',b'NOT a result')
        self.json(views/'manifest.json',meta);self.finish(views,meta)
        with self.assertRaisesRegex(RuntimeError,'output/resource'):self.inventory()

    def test_exact_queue_recovery_v2_1_exporter_path_is_accepted(self):
        views=self.direct();meta=json.loads((views/'manifest.json').read_text())
        producer=self.root/'private/reviewed-packet-fastviews-v2-3/scripts/export_packet_views_fast_v2_1.py';self.file(producer,b'exact reviewed queue-only release')
        meta.update(report_producer_kind='FAST_DIRECT_PACKET_VIEWS_V2',producer_code_sha256=self.pin(producer)['sha256'])
        self.json(views/'manifest.json',meta);self.finish(views,meta)
        plan,_=self.inventory();self.assertEqual(sum('/vcf/sample_' in e['RelativePath'] for e in plan['Files']),8)

    def test_capacity_counts_resume_bytes_preserves_all_required_files(self):
        destination=self.base/'results';staging=self.base/'small-staging';destination.mkdir();staging.mkdir()
        entry=dict(RelativePath='outputs/cohort-sample-views-v1/masters/cohort.master.tsv.gz',SizeBytes=100,DownloadFileName='cohort.master.tsv.gz',Sha256='1'*64)
        (staging/'cohort.master.tsv.gz.part').write_bytes(b'x'*40)
        with patch.object(planner.shutil,'disk_usage',return_value=SimpleNamespace(free=159)):
            value=planner.local_capacity(dict(Files=[entry]),destination,staging,reserve=100)
        self.assertEqual(value['required_remaining_bytes'],60);self.assertTrue(value['status'].startswith('BLOCKED'))
        self.assertTrue(value['all_required_outputs_retained']);self.assertFalse(value['patient_file_hashes_computed'])
        self.assertEqual((staging/'cohort.master.tsv.gz.part').stat().st_size,40)

    def overlay_applied_scores_need_real_source_receipt(self,version):
        views=self.direct();folder=self.root/'outputs/fully-annotated-vcf-v1';full=json.loads((folder/'manifest.json').read_text())
        if version=='v3':
            core_path=Path(full['core_manifest']['path']);core=json.loads(core_path.read_text())
            core['producer_code_sha256']['resume_parent_driver']='8'*64
            provenance=folder/'resume-provenance.json';self.json(provenance,dict(parent_reviewed_rule='synthetic old reviewed recipe',source_origin_mapping=[]))
            core['resume_provenance']=self.pin(provenance);self.json(core_path,core)
            full['core_manifest']=self.pin(core_path);full['producer_code_sha256']=core['producer_code_sha256']
        package=self.root/('private/spliceai-allele-overlay-candidate-'+version+'/scripts')
        names=['apply_spliceai_allele_overlay.py','continue_spliceai_allele_overlay.py','spliceai_allele_info_codec.py','spliceai_packet_context.py','embed_vcf_evidence.py']
        for name in names:self.file(package/name,b'reviewed overlay code fixture')
        prediction=self.root/'outputs/spliceai-reference-context-v1/predictions.tsv.gz';self.file(prediction,b'actual applied score rows fixture')
        item=dict(self.pin(prediction),relative_path='applied-evidence/spliceai/predictions.tsv.gz',role='APPLIED_SPLICEAI_ALLELE_PREDICTIONS',scope='GENCODE24_GENE_UNION_NOT_ENST49_MAPPING')
        score_receipt=prediction.parent/'workflow-receipt.json'
        self.json(score_receipt,dict(status='PASS',stage='spliceai-reference-context-v1',artifacts={'predictions.tsv.gz':{k:item[k] for k in ('bytes','sha256')}}))
        full.update(projection=dict(projection_contract='ALLELE_INFO_V1_ROW_PROJECTION_ON_READ'),core33_packets_byte_identical=True,
            core_producer_code_sha256=full['producer_code_sha256'],overlay_code_sha256={str(package/n):self.pin(package/n)['sha256'] for n in names},
            spliceai_predictions=self.pin(prediction),applied_artifacts=[item])
        self.json(folder/'manifest.json',full);wf=json.loads((folder/'workflow-receipt.json').read_text());wf['artifacts']['manifest.json']={k:self.pin(folder/'manifest.json')[k] for k in ('bytes','sha256')};self.json(folder/'workflow-receipt.json',wf)
        meta=json.loads((views/'manifest.json').read_text());meta['source_pins']['fully_manifest']=self.pin(folder/'manifest.json');meta['source_pins']['core_manifest']=full['core_manifest']
        meta['artifacts'][item['relative_path']]=self.file(views/item['relative_path'],prediction.read_bytes())
        self.json(views/'manifest.json',meta);self.finish(views,meta)
        plan,_=self.inventory();self.assertTrue(any(e['RelativePath'].endswith(item['relative_path']) for e in plan['Files']))
        wf=json.loads(score_receipt.read_text());wf['status']='FAILED';self.json(score_receipt,wf)
        with self.assertRaisesRegex(RuntimeError,'producer receipt'):self.inventory()

    def test_v2_overlay_applied_scores_need_real_source_receipt(self):
        self.overlay_applied_scores_need_real_source_receipt('v2')

    def test_v3_resume_parent_and_overlay_paths_preserve_strict_backup_contract(self):
        self.overlay_applied_scores_need_real_source_receipt('v3')


if __name__=='__main__':unittest.main()
