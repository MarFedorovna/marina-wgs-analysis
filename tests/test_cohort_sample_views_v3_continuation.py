import datetime as dt
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('vcf_first_controller',ROOT/'scripts/continue_cohort_sample_views_v3.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class VcfFirstContinuationTests(unittest.TestCase):
    def test_primary_and_joined_alone_do_not_start_views(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);workflow=root/'workflow';fully=root/'fully-annotated-vcf-v1/manifest.json'
            rule=dict(workflow_root=str(workflow),fully_vcf_manifest=str(fully),fully_manifest_status='PASS_FULLY_ANNOTATED_VCF')
            for stage in ('06_annotate','07_report'):
                path=workflow/'stages'/stage;path.mkdir(parents=True);(path/'workflow-receipt.json').write_text(json.dumps(dict(status='PASS')))
            self.assertFalse(m.readiness(rule));fully.parent.mkdir();fully.write_text(json.dumps(dict(status='PASS_FULLY_ANNOTATED_VCF')))
            self.assertFalse(m.readiness(rule));(fully.parent/'workflow-receipt.json').write_text(json.dumps(dict(status='PASS',stage='fully-annotated-vcf-v1')));self.assertTrue(m.readiness(rule))

    def test_rule_codec_and_engine_pins_and_original_deadline(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp).resolve();scripts=ROOT/'scripts';codecs=scripts/'embed_vcf_evidence.py'
            rule=dict(schema_version=3,review_status='ROOT_REVIEWED_EXACT_VCF_FIRST_COHORT_SAMPLE_VIEWS_V3_RULE',reviewed_utc='2026-11-05T22:00:00Z',
                code_sha256={str(p):m.sha(p) for p in (scripts/'export_cohort_sample_views_v3.py',scripts/'continue_cohort_sample_views_v3.py')},
                ranking_module=str(scripts/'annotation_stages.py'),ranking_sha256=m.sha(scripts/'annotation_stages.py'),
                v1_engine=str(scripts/'export_cohort_sample_views.py'),v1_engine_sha256=m.sha(scripts/'export_cohort_sample_views.py'),
                v2_helper=str(scripts/'export_cohort_sample_views_v2.py'),v2_helper_sha256=m.sha(scripts/'export_cohort_sample_views_v2.py'),codec_module=str(codecs),codec_sha256=m.sha(codecs),
                fully_vcf_manifest=str(root/'fully-annotated-vcf-v1/manifest.json'),fully_manifest_status='PASS_FULLY_ANNOTATED_VCF',
                control_joined_manifest=str(root/'combined-annotation-master-v1/manifest.json'),workflow_root=str(root/'workflow'),output=str(root/'cohort-sample-views-v1'),state_directory=str(root/'private-state'),
                packet_append_fields=['spliceai_status','spliceai_gene','spliceai_DS_AG','spliceai_DS_AL','spliceai_DS_DG','spliceai_DS_DL','spliceai_DP_AG','spliceai_DP_AL','spliceai_DP_DG','spliceai_DP_DL','spliceai_prediction_json','spliceai_source_json','spliceai_transcript_scope'],
                spliceai_manifest=str(root/'spliceai-reference-context-v1/manifest.json'),spliceai_manifest_status='PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT',
                python_executable=str(Path(sys.executable).absolute()),pysam_version='0.24.1',deadline_utc=m.DEADLINE.isoformat(),max_hours=8,interval_seconds=60,no_network=True,no_primary_workflow_mutation=True)
            path=root/'rule.json';path.write_text(json.dumps(rule));current=m.DEADLINE-dt.timedelta(hours=1)
            m.validate(rule,path,m.sha(path),current)
            with self.assertRaisesRegex(ValueError,'deadline'):m.validate(rule,path,m.sha(path),m.DEADLINE)
            rule['codec_sha256']='0'*64
            with self.assertRaisesRegex(ValueError,'codec'):m.validate(rule,path,m.sha(path),current)


if __name__=='__main__':unittest.main()
