import datetime as dt
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('joined_controller',ROOT/'scripts/continue_cohort_sample_views_v2.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class JoinedContinuationTests(unittest.TestCase):
    def test_primary_pass_alone_never_starts_baseline_only_export(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);workflow=root/'workflow';joined=root/'combined-annotation-master-v1/manifest.json'
            for stage in ('06_annotate','07_report'):
                p=workflow/'stages'/stage;p.mkdir(parents=True);(p/'workflow-receipt.json').write_text(json.dumps(dict(status='PASS')))
            self.assertFalse(m.readiness(workflow,joined))
            joined.parent.mkdir();joined.write_text(json.dumps(dict(status='PASS_COMBINED_ANNOTATION_MASTER')))
            self.assertFalse(m.readiness(workflow,joined))
            (joined.parent/'workflow-receipt.json').write_text(json.dumps(dict(status='PASS')))
            self.assertTrue(m.readiness(workflow,joined))

    def test_rule_code_base_pins_and_original_deadline(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp).resolve();scripts=ROOT/'scripts'
            rule=dict(schema_version=2,review_status='ROOT_REVIEWED_EXACT_JOINED_COHORT_SAMPLE_VIEWS_V2_RULE',reviewed_utc='2026-11-05T20:00:00Z',
                code_sha256={str(p.resolve()):m.sha(p) for p in (scripts/'continue_cohort_sample_views_v2.py',scripts/'export_cohort_sample_views_v2.py')},
                ranking_module=str(scripts/'annotation_stages.py'),ranking_sha256=m.sha(scripts/'annotation_stages.py'),
                base_exporter=str(scripts/'export_cohort_sample_views.py'),base_exporter_sha256=m.sha(scripts/'export_cohort_sample_views.py'),
                joined_manifest=str(root/'combined-annotation-master-v1/manifest.json'),workflow_root=str(root/'workflow'),output=str(root/'output'),state_directory=str(root/'state'),
                python_executable=str(Path(sys.executable).absolute()),pysam_version='0.24.1',deadline_utc=m.ABSOLUTE_DEADLINE.isoformat(),max_hours=8,interval_seconds=60,no_network=True,no_primary_workflow_mutation=True)
            path=root/'rule.json';path.write_text(json.dumps(rule));current=m.ABSOLUTE_DEADLINE-dt.timedelta(hours=1)
            m.validate(rule,path,m.sha(path),current)
            self.assertEqual(m.run(path,None)['status'],'NOT_STARTED_ROOT_REVIEW_REQUIRED')
            with self.assertRaisesRegex(ValueError,'deadline'):m.validate(rule,path,m.sha(path),m.ABSOLUTE_DEADLINE)
            rule['base_exporter_sha256']='0'*64
            with self.assertRaisesRegex(ValueError,'Frozen source'):m.validate(rule,path,m.sha(path),current)


if __name__=='__main__':unittest.main()
