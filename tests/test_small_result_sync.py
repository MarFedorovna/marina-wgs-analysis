import base64
import datetime as dt
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('small_sync',ROOT/'transfer/sync_completed_small_results.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class FakeSftp:
    def __init__(self,data,mtime=7):self.data=data;self.mtime=mtime;self.stats=0;self.offsets=[]
    def stat(self,path):self.stats+=1;return SimpleNamespace(st_size=len(self.data),st_mtime=self.mtime)
    def open(self,path,mode):return io.BytesIO(self.data)
    def blocks(self,stream,start,size):
        self.offsets.append(start);stream.seek(start)
        for offset in range(0,size,32768):yield stream.read(min(32768,size-offset))


class SmallResultSyncTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name).resolve()
        self.profile=dict(schema_version=1,plan=str(self.root/'plan.json'),auth_config=str(self.root/'MUST_NOT_READ_AUTH_CONFIG'),
            download_directory=str(self.root/'downloads'),destination=str(self.root/'results'),state_directory=str(self.root/'state'),remote_root='/srv/marina',
            wsl_distribution='Ubuntu-22.04',max_hours=8,interval_seconds=60,deadline_utc=m.DEADLINE.isoformat(),max_file_bytes=32*1024**3,max_total_bytes=64*1024**3,
            min_free_bytes=15*1024**3,required_completed_roots=sorted(m.SIDE_ROOTS),allow_onedrive=True,pause_file=str(self.root/'PAUSE'))
        source='outputs/workflow-v1/stages/07_report/workflow-receipt.json'
        receipt=dict(status='PASS',identity='1'*64,artifacts={'master_evidence.tsv.gz':dict(bytes=12,sha256='a'*64)})
        raw=json.dumps(receipt).encode()
        self.entry=dict(RelativePath='outputs/workflow-v1/stages/07_report/master_evidence.tsv.gz',ServerPath='/srv/marina/outputs/workflow-v1/stages/07_report/master_evidence.tsv.gz',
            SizeBytes=12,Sha256='a'*64,SourceReceipt=source,DownloadFileName='completed_master.tsv.gz')
        self.plan=dict(SchemaVersion=1,PlanGeneratorVersion=2,Files=[self.entry],Receipts=[dict(RelativePath=source,ContentBase64=base64.b64encode(raw).decode(),SizeBytes=len(raw),Sha256=hashlib.sha256(raw).hexdigest())])

    def test_dry_inventory_never_reads_auth_or_creates_transfer_paths(self):
        Path(self.profile['plan']).write_text(json.dumps(self.plan));before=set(self.root.iterdir())
        result=m.run(self.profile)
        self.assertEqual(result['status'],'DRY_RUN_NO_AUTH_CONFIG_KEY_NETWORK_OR_IMPORT');self.assertEqual(result['selected_files'],1)
        self.assertEqual(set(self.root.iterdir()),before);self.assertFalse(Path(self.profile['auth_config']).exists())

    def test_strict_selection_no_bam_gvcf_resources_and_sample_vcf_yes(self):
        for path in ['outputs/remaining-align/S/align/marked.bam','outputs/workflow-v1/stages/03_call/S.g.vcf.gz',
            'resources/cadd.tsv.gz','outputs/database-cache-v1/source.tsv.gz','outputs/cohort-sample-views-v1/resources/db.tsv.gz']:
            self.assertFalse(m.selected(path),path)
        for i in range(1,5):
            for suffix in ('','.tbi'):self.assertTrue(m.selected('outputs/cohort-sample-views-v1/vcf/sample_%04d.annotated.vcf.gz%s'%(i,suffix)))
        self.assertTrue(m.selected('outputs/combined-annotation-master-v1/applied-evidence/population/cohort/ClinVar-SCV-long.tsv.gz'))
        with self.assertRaises(ValueError):m.selected('outputs/workflow-v1/../resources/db.tsv.gz')

    def test_receipt_tamper_cap_and_future_completion_not_claimed(self):
        files,summary=m.inventory(self.plan,self.profile)
        self.assertEqual(summary['missing_completed_roots'],sorted(m.SIDE_ROOTS));self.assertEqual(files[0]['Sha256'],'a'*64)
        self.plan['Files'][0]['Sha256']='b'*64
        with self.assertRaisesRegex(ValueError,'producer receipt'):m.inventory(self.plan,self.profile)
        self.plan['Files'][0]['Sha256']='a'*64;self.profile['max_file_bytes']=1
        with self.assertRaisesRegex(ValueError,'cap'):m.inventory(self.plan,self.profile)

    def test_resume_reads_only_remaining_prefix_and_completion_is_size_only(self):
        data=b'ACGT'*40000;entry=dict(self.entry,SizeBytes=len(data));directory=self.root/'downloads';directory.mkdir();sftp=FakeSftp(data)
        part=directory/(entry['DownloadFileName']+'.part');part.write_bytes(data[:12345])
        identity=dict(RelativePath=entry['RelativePath'],RemotePath=entry['ServerPath'],ExpectedBytes=len(data),ExpectedServerSha256=entry['Sha256'],RemoteMtime=7,SourceReceipt=entry['SourceReceipt'],Verification='size_only')
        (directory/(part.name+'.source.json')).write_text(json.dumps(identity))
        final=m.download_one(sftp,sftp.blocks,entry,directory,m.now()+dt.timedelta(minutes=1),Path(self.profile['pause_file']),reserve=0)
        self.assertEqual(final.read_bytes(),data);self.assertEqual(sftp.offsets[0],12345)
        complete=json.loads(final.with_name(final.name+'.completed.json').read_text());self.assertIsNone(complete['ActualSha256'])
        self.assertNotIn('auth',json.dumps(complete).lower())
        self.assertEqual(m.download_one(sftp,sftp.blocks,entry,directory,m.now()+dt.timedelta(minutes=1),Path(self.profile['pause_file']),reserve=0),final)

    def test_changed_remote_identity_and_short_chunk_preserve_partial(self):
        data=b'A'*70000;entry=dict(self.entry,SizeBytes=len(data));directory=self.root/'downloads';directory.mkdir();sftp=FakeSftp(data)
        def short(stream,start,size):yield b'A'*10
        with self.assertRaisesRegex(ValueError,'contiguous'):m.download_one(sftp,short,entry,directory,m.now()+dt.timedelta(minutes=1),Path(self.profile['pause_file']),reserve=0)
        part=directory/(entry['DownloadFileName']+'.part');self.assertEqual(part.stat().st_size,0)
        sftp.mtime=8
        with self.assertRaisesRegex(ValueError,'source identity'):m.download_one(sftp,sftp.blocks,entry,directory,m.now()+dt.timedelta(minutes=1),Path(self.profile['pause_file']),reserve=0)
        self.assertEqual(part.stat().st_size,0)

    def test_deadline_and_pause_before_any_source_access(self):
        directory=self.root/'downloads';directory.mkdir();sftp=FakeSftp(b'A')
        with self.assertRaisesRegex(ValueError,'deadline'):m.download_one(sftp,sftp.blocks,self.entry,directory,m.now()-dt.timedelta(seconds=1),self.root/'pause',reserve=0)
        (self.root/'pause').touch()
        with self.assertRaisesRegex(ValueError,'pause'):m.download_one(sftp,sftp.blocks,self.entry,directory,m.now()+dt.timedelta(minutes=1),self.root/'pause',reserve=0)
        self.assertEqual(sftp.stats,0);self.assertEqual(list(directory.iterdir()),[])


if __name__=='__main__':unittest.main()
