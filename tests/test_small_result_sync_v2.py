import base64
import datetime as dt
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'transfer'))
import sync_completed_small_results as old
import sync_completed_small_results_v2 as m
sys.path.insert(0,str(Path(__file__).parent))
from test_small_result_sync import FakeSftp


class FinalRootSyncTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.profile=dict(schema_version=1,plan=str(self.root/'plan.json'),auth_config=str(self.root/'AUTH_MUST_NOT_BE_READ'),download_directory=str(self.root/'downloads'),
            destination=str(self.root/'results'),state_directory=str(self.root/'new-state'),remote_root='/srv/marina',wsl_distribution='Ubuntu-22.04',max_hours=8,interval_seconds=60,
            deadline_utc=m.DEADLINE.isoformat(),max_file_bytes=32*1024**3,max_total_bytes=64*1024**3,min_free_bytes=15*1024**3,required_completed_roots=['cohort-sample-views-v1'],allow_onedrive=True,pause_file=str(self.root/'PAUSE'))
        source='outputs/cohort-sample-views-v1/manifest.json'
        names=['vcf/cohort.fully-annotated.vcf.gz','vcf/cohort.fully-annotated.vcf.gz.tbi',*[f'vcf/sample_{i:04d}.annotated.vcf.gz{s}' for i in range(1,5) for s in ('','.tbi')]]
        receipt=dict(status='PASS',schema_version=1,synthetic=False,artifacts={name:dict(bytes=12,sha256='a'*64) for name in names})
        raw=json.dumps(receipt).encode()
        self.plan=dict(SchemaVersion=1,PlanGeneratorVersion=3,Receipts=[dict(RelativePath=source,ContentBase64=base64.b64encode(raw).decode(),SizeBytes=len(raw),Sha256=hashlib.sha256(raw).hexdigest())],
            Files=[dict(RelativePath='outputs/cohort-sample-views-v1/'+name,ServerPath='/srv/marina/outputs/cohort-sample-views-v1/'+name,SizeBytes=12,Sha256='a'*64,SourceReceipt=source,DownloadFileName=name.replace('/','__')) for name in names])

    def test_new_final_root_does_not_wait_for_four_unexecuted_legacy_roots(self):
        m.validate_profile(self.profile);_,summary=m.inventory(self.plan,self.profile)
        self.assertEqual(summary['missing_completed_roots'],[]);self.assertEqual(summary['individual_annotated_vcf_artifacts'],10)
        legacy=dict(self.profile,required_completed_roots=sorted(old.SIDE_ROOTS));_,old_summary=old.inventory(self.plan,legacy)
        self.assertEqual(len(old_summary['missing_completed_roots']),4)
        with self.assertRaisesRegex(ValueError,'final VCF-first'):m.validate_profile(legacy)

    def test_dry_inventory_never_reads_auth_or_creates_transfer_state(self):
        Path(self.profile['plan']).write_text(json.dumps(self.plan));before=set(self.root.iterdir())
        result=m.run(self.profile)
        self.assertEqual(result['status'],'DRY_RUN_NO_AUTH_CONFIG_KEY_NETWORK_OR_IMPORT');self.assertEqual(set(self.root.iterdir()),before)
        self.assertFalse(Path(self.profile['auth_config']).exists());self.assertFalse(result['patient_file_hashes_computed'])

    def test_partial_resumes_without_new_hash_and_optional_mt_does_not_block(self):
        entry=self.plan['Files'][0];data=b'ACGT'*10000;entry=dict(entry,SizeBytes=len(data));directory=self.root/'downloads';directory.mkdir()
        partial=directory/(entry['DownloadFileName']+'.part');partial.write_bytes(data[:1234]);sftp=FakeSftp(data)
        identity=dict(RelativePath=entry['RelativePath'],RemotePath=entry['ServerPath'],ExpectedBytes=len(data),ExpectedServerSha256=entry['Sha256'],RemoteMtime=7,SourceReceipt=entry['SourceReceipt'],Verification='size_only')
        partial.with_name(partial.name+'.source.json').write_text(json.dumps(identity))
        final=m.download_one(sftp,sftp.blocks,entry,directory,m.now()+dt.timedelta(minutes=1),Path(self.profile['pause_file']),reserve=0)
        self.assertEqual(sftp.offsets,[1234]);self.assertEqual(final.read_bytes(),data)
        self.assertIsNone(json.loads(final.with_name(final.name+'.completed.json').read_text())['ActualSha256'])
        self.assertTrue(m.selected('outputs/mitochondrial-reference-context-v1/cohort/MT-population-reference-long.tsv.gz'))
        _,summary=m.inventory(self.plan,self.profile);self.assertEqual(summary['missing_completed_roots'],[])


if __name__=='__main__':unittest.main()
