import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('component_download',Path(__file__).resolve().parents[1]/'transfer/download_completed_component.py')
component=importlib.util.module_from_spec(spec);spec.loader.exec_module(component)

class CompletedComponentDownload(unittest.TestCase):
    def plan(self,root):
        return dict(schema_version=1,status='ROOT_REVIEWED_COMPLETED_COMPONENT_DOWNLOAD',
            scientific_status='COMPLETED_VEP_COMPONENT_NOT_FINAL_ALL_SOURCES_OR_REPORT',
            deadline_utc='2026-11-06T03:46:48+00:00',destination=str(root/'not_created'),
            state_directory=str(root/'also_not_created'),pause_file=str(root/'pause'),source_provenance={},
            files=[dict(DownloadFileName='cohort.vep.vcf.gz'+s,ServerPath='/srv/private/completed.vcf.gz'+s,
                        SizeBytes=n,Sha256='a'*64) for s,n in [('',100),('.tbi',50)]])

    def test_dry_run_never_reads_auth_or_creates_destination(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);result=component.run(self.plan(root),root/'MUST_NOT_READ')
            self.assertEqual('PLAN_ONLY_NO_AUTH_OR_NETWORK',result['status'])
            self.assertEqual(150,result['bytes']);self.assertEqual([],list(root.iterdir()))

    def test_unpaired_index_is_rejected(self):
        plan=self.plan(Path('/unused'));plan['files'][1]['ServerPath']='/srv/private/different.vcf.gz.tbi'
        with self.assertRaisesRegex(ValueError,'same completed component pair'):component.run(plan)

    def test_wrong_scope_or_unsafe_output_refused(self):
        plan=self.plan(Path('/unused'));wrong=copy.deepcopy(plan);wrong['scientific_status']='ALL_SOURCES_COMPLETE'
        with self.assertRaises(ValueError):component.run(wrong)
        wrong=copy.deepcopy(plan);wrong['files'][0]['DownloadFileName']='../other.vcf.gz'
        with self.assertRaises(ValueError):component.run(wrong)
        wrong=copy.deepcopy(plan);wrong['files'][0]['Sha256']='not-a-source-digest'
        with self.assertRaises(ValueError):component.run(wrong)

if __name__=='__main__':unittest.main()
