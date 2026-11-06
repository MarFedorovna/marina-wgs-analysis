import base64
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'transfer'))
import download_cohort_views_dataset as m


def plan():
    artifacts={'tables/cohort.top20.tsv':dict(bytes=12,sha256='a'*64),'masters/cohort.master.tsv.gz':dict(bytes=30*1024**2,sha256='b'*64)}
    manifest=dict(schema_version=1,status='PASS',synthetic=False,report_input_kind='FULLY_ANNOTATED_VCF_PACKETS_V1',report_producer_kind='FAST_DIRECT_PACKET_VIEWS_V2',
        sample_packet_scope_verified=True,public_header_bundles_verified=True,no_external_requests=True,no_legacy_receipt_created=True,primary_workflow_changed=False,automatic_ACMG=False,producer_code_sha256='c'*64,artifacts=artifacts)
    def captured(name,value):
        raw=json.dumps(value).encode();return dict(path=m.REMOTE+name,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),content_base64=base64.b64encode(raw).decode())
    mp=captured('manifest.json',manifest);receipt=dict(schema_version=1,status='PASS',synthetic=False,stage='cohort-sample-views-v1',identity='d'*64,producer_code_sha256='c'*64,artifacts=dict(artifacts,**{'manifest.json':{k:mp[k] for k in ('bytes','sha256')}}))
    return dict(schema_version=1,status='ROOT_REVIEWED_COHORT_VIEWS_DATASET_DOWNLOAD',deadline_utc=m.END.isoformat(),destination='NOT_OPENED',state_directory='NOT_CREATED',pause_file='NOT_OPENED',source_provenance=dict(manifest=mp,receipt=captured('workflow-receipt.json',receipt)),source_aliases=[])


class FakeSftp:
    def __init__(self,data):self.data=data
    def open(self,path,mode):return io.BytesIO(self.data)
    def stat(self,path):return SimpleNamespace(st_size=len(self.data),st_mtime=7)
    def blocks(self,stream,start,size):
        stream.seek(start)
        for offset in range(0,size,32768):yield stream.read(min(32768,size-offset))


class DatasetTests(unittest.TestCase):
    def test_dry_and_small_phase_are_explicit_no_auth(self):
        p=plan();full=m.run(p,auth='MUST_NOT_READ');small=m.run(p,auth='MUST_NOT_READ',small_only=True)
        self.assertIn('NO_AUTH',full['status']);self.assertEqual(small['pending_large_files'],1);self.assertFalse(small['whole_dataset_complete'])
    def test_wrong_receipt_fake_pass_and_resource_paths_rejected(self):
        p=plan();p['source_provenance']['receipt']['sha256']='f'*64
        with self.assertRaises(ValueError):m.validate(p)
        for name in ['resources/db.tsv.gz','x/../database.tsv.gz','calls/sample.g.vcf.gz','x.bam','/tmp/file.tsv','C:/evil.tsv','x\x00.tsv','NUL.tsv','a/COM1.tsv','bad. /file.tsv']:
            with self.assertRaises(ValueError):m.safe_relative(name)
    def test_real_physical_preallocation_does_not_create_final(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'own.part';m.preallocate(p,1024**2)
            self.assertEqual(p.stat().st_size,1024**2);self.assertGreaterEqual(m.allocation_bytes(p),1024**2)
            with self.assertRaises(ValueError):m.preallocate(p,1)
    def fixture(self,root):
        name='tables/public.tsv.gz';data=bytes(range(251))*211;expected=dict(bytes=len(data),sha256='a'*64);target=root/'public.tsv.gz';part=root/'public.tsv.gz.part';meta=root/'chunks.json';note=root/'complete.json';m.preallocate(part,len(data))
        source=dict(artifact=name,bytes=len(data),sha256='a'*64,remote_mtime=7,manifest_sha256='c'*64);lock=threading.Lock()
        m.atomic(meta,dict(schema_version=1,artifact=name,expected=expected,source=source,chunk_bytes=m.CHUNK,file_id=part.stat().st_ino,completed=[]))
        tasks=[(name,expected,target,part,meta,note,source,lock,i,i*m.CHUNK,min(m.CHUNK,len(data)-i*m.CHUNK)) for i in range((len(data)+m.CHUNK-1)//m.CHUNK)]
        return data,target,part,meta,note,tasks
    def test_sixteen_workers_nonoverlapping_exact_bytes_and_durable_resume(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(m,'CHUNK',8192):
            data,target,part,meta,note,tasks=self.fixture(Path(temp));stop=threading.Event();deadline=m.now()+dt.timedelta(minutes=1);pause=Path(temp)/'PAUSE'
            s=FakeSftp(data);m.transfer_chunk(s,s.blocks,tasks[0],stop,deadline,pause,lambda n:None)
            self.assertEqual(m.load(meta)['completed'],[0]);self.assertFalse(target.exists())
            # Owned but uncommitted bytes are rewritten from the chunk boundary.
            with part.open('r+b') as stream:stream.seek(8192);stream.write(b'partial-failed-chunk')
            def work(task):
                remote=FakeSftp(data);m.transfer_chunk(remote,remote.blocks,task,stop,deadline,pause,lambda n:None)
            with ThreadPoolExecutor(max_workers=16) as pool:list(pool.map(work,tasks[1:]))
            self.assertEqual(target.read_bytes(),data);self.assertFalse(part.exists());self.assertIsNone(m.load(note)['actual_sha256'])
            self.assertEqual(m.load(meta)['completed'],list(range(len(tasks))))
    def test_failed_chunk_keeps_uncommitted_bitmap_and_prevents_final(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(m,'CHUNK',8192):
            data,target,part,meta,note,tasks=self.fixture(Path(temp));remote=FakeSftp(data)
            def short(stream,start,size):yield b'bad'
            with self.assertRaises(ValueError):m.transfer_chunk(remote,short,tasks[0],threading.Event(),m.now()+dt.timedelta(minutes=1),Path(temp)/'PAUSE',lambda n:None)
            self.assertEqual(m.load(meta)['completed'],[]);self.assertTrue(part.exists());self.assertFalse(target.exists())
    def test_short_local_write_cannot_commit_preallocated_zero_tail(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(m,'CHUNK',8192):
            data,target,part,meta,note,tasks=self.fixture(Path(temp));remote=FakeSftp(data);original=Path.open
            class ShortWriter:
                def __init__(self,stream):self.stream=stream;self.first=True
                def __enter__(self):return self
                def __exit__(self,*args):self.stream.close()
                def seek(self,offset):return self.stream.seek(offset)
                def write(self,data):
                    if self.first:self.first=False;return self.stream.write(data[:3])
                    return 0
                def fileno(self):return self.stream.fileno()
            def opened(path,*args,**kwargs):
                stream=original(path,*args,**kwargs)
                return ShortWriter(stream) if path==part and args and args[0]=='r+b' else stream
            with patch.object(Path,'open',opened):
                with self.assertRaisesRegex(ValueError,'local write'):m.transfer_chunk(remote,remote.blocks,tasks[0],threading.Event(),m.now()+dt.timedelta(minutes=1),Path(temp)/'PAUSE',lambda n:None)
            self.assertEqual(part.stat().st_size,len(data));self.assertEqual(m.load(meta)['completed'],[]);self.assertFalse(target.exists())
    def test_alias_is_same_inode_with_source_completion_identity(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'main.vcf.gz';source.write_bytes(b'completed source');expected=dict(bytes=source.stat().st_size,sha256='a'*64)
            receipt=root/'source.completed.json';m.atomic(receipt,dict(Status='COMPLETED_SIZE_ONLY',SourceIdentity=dict(ExpectedBytes=expected['bytes'],ExpectedServerSha256=expected['sha256']),LocalMtimeNs=source.stat().st_mtime_ns))
            target=root/'mirror/cohort.vcf.gz';note=root/'alias.json';m.alias_file(dict(local_path=str(source),completion_receipt=str(receipt)),expected,target,note)
            self.assertTrue(target.samefile(source));self.assertTrue(m.completed(target,note,expected))


if __name__=='__main__':unittest.main()
