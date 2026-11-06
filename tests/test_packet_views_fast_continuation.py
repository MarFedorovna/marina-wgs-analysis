import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import continue_packet_views_fast as m


class ContinuationGuards(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.path=self.root/'rule.json'
        def item(name):
            path=ROOT/'scripts'/name
            return dict(path=str(path),bytes=path.stat().st_size,sha256=m.sha(path))
        self.rule=dict(schema_version=1,root_reviewed=True,synthetic=False,controller_sha256=m.sha(m.__file__),
            interval_seconds=10,max_hours=8,workers=2,deadline_utc=m.END.isoformat(),state_directory=str(self.root/'state'),
            fully_manifest_path='/srv/marina/outputs/fully-annotated-vcf-v1/manifest.json',output='/srv/marina/outputs/cohort-sample-views-v1',
            authority_expected={'core.source_pins':{'synthetic_test_only':'not_execution_authority'}},exporter=item('export_packet_views_fast.py'),
            proof_module=item('fast_vcf_evidence_candidate.py'),codec_module=item('embed_vcf_evidence.py'),view_policy=item('export_cohort_sample_views.py'),ranking_module=item('annotation_stages.py'))

    def write(self):
        self.path.write_text(json.dumps(self.rule));return m.sha(self.path)

    def test_changed_rule_digest_rejected(self):
        digest=self.write();self.path.write_text(self.path.read_text()+' ')
        with self.assertRaisesRegex(ValueError,'rule SHA256'):m.continue_views(self.path,digest)

    def test_status_only_does_not_create_lock_or_work(self):
        digest=self.write();result=m.continue_views(self.path,digest)
        self.assertEqual(result['status'],'STATUS_ONLY_NO_EXECUTION');self.assertFalse((self.root/'state').exists())

    def test_existing_own_lock_refused_no_restarts(self):
        state=self.root/'state';state.mkdir();(state/'owns.lock').write_text('another-controller')
        digest=self.write()
        with self.assertRaises(FileExistsError):m.continue_views(self.path,digest,execute=True)
        self.assertEqual((state/'owns.lock').read_text(),'another-controller')

    def test_unknown_authority_or_extra_scope_rejected(self):
        self.rule['authority_expected']={'legacy07.status':'PASS'}
        with self.assertRaisesRegex(ValueError,'namespace'):m.checked_rule(self.path,self.write())
        self.rule['authority_expected']={'core.source_pins':{}}
        self.rule['fully_manifest_path']='/srv/marina/resources/fake/manifest.json'
        with self.assertRaisesRegex(ValueError,'producer root'):m.checked_rule(self.path,self.write())

    def test_changed_helper_rejected(self):
        self.rule['proof_module']['sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'helper'):m.checked_rule(self.path,self.write())


if __name__=='__main__':unittest.main()
