import datetime as dt
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('controller',ROOT/'scripts/continue_cohort_sample_views.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class ViewsContinuationTests(unittest.TestCase):
    def rule(self,root):
        scripts=ROOT/'scripts';ranking=scripts/'annotation_stages.py'
        rule=dict(schema_version=1,review_status='ROOT_REVIEWED_EXACT_COHORT_SAMPLE_VIEWS_RULE',reviewed_utc='2026-11-05T19:00:00Z',
            code_sha256={str(p.resolve()):m.sha(p) for p in (scripts/'continue_cohort_sample_views.py',scripts/'export_cohort_sample_views.py')},
            ranking_module=str(ranking.resolve()),ranking_sha256=m.sha(ranking),workflow_root=str(root/'workflow'),output=str(root/'output'),state_directory=str(root/'state'),
            python_executable=str(Path(sys.executable).absolute()),pysam_version='0.24.1',deadline_utc=m.ABSOLUTE_DEADLINE.isoformat(),max_hours=8,
            interval_seconds=60,no_network=True,no_primary_workflow_mutation=True)
        path=root/'rule.json';path.write_text(json.dumps(rule));return rule,path

    def test_root_review_code_ranking_and_unknown_fields(self):
        with tempfile.TemporaryDirectory() as temp:
            rule,path=self.rule(Path(temp));time=m.ABSOLUTE_DEADLINE-dt.timedelta(hours=1)
            m.validate(rule,path,m.sha(path),time)
            with self.assertRaises(ValueError):m.validate(rule,path,'0'*64,time)
            rule['ranking_sha256']='0'*64
            with self.assertRaisesRegex(ValueError,'ranking'):m.validate(rule,path,m.sha(path),time)
            rule,path=self.rule(Path(temp));rule['unknown']=True
            with self.assertRaisesRegex(ValueError,'Unknown'):m.validate(rule,path,m.sha(path),time)

    def test_no_default_execution_and_expired_rule(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);rule,path=self.rule(root)
            self.assertEqual(m.run(path,None)['status'],'NOT_STARTED_ROOT_REVIEW_REQUIRED');self.assertFalse((root/'output').exists())
            with self.assertRaisesRegex(ValueError,'deadline'):m.validate(rule,path,m.sha(path),m.ABSOLUTE_DEADLINE)
            self.assertFalse(m.readiness(root/'workflow'))
            for stage in ('06_annotate','07_report'):
                p=root/'workflow/stages'/stage;p.mkdir(parents=True);(p/'workflow-receipt.json').write_text(json.dumps(dict(status='PASS')))
            self.assertTrue(m.readiness(root/'workflow'))
            (root/'workflow/stages/07_report/workflow-receipt.json').write_text(json.dumps(dict(status='FAILED')))
            with self.assertRaises(ValueError):m.readiness(root/'workflow')


if __name__=='__main__':unittest.main()
