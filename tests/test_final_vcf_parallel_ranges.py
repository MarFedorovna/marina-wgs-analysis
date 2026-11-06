import datetime as dt
import io
import json
from pathlib import Path
import sys
import tempfile
import threading
from types import SimpleNamespace
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'transfer'));sys.path.insert(0,str(Path(__file__).parent))
import download_final_vcf_parallel_ranges as m
import sync_completed_small_results as helper
from test_download_final_annotated_vcf import FinalPairTests


class FakeSftp:
    def __init__(self,data):self.data=data;self.starts=[]
    def open(self,path,mode):return io.BytesIO(self.data)
    def stat(self,path):return SimpleNamespace(st_size=len(self.data),st_mtime=7)
    def blocks(self,stream,start,size):
        self.starts.append(start);stream.seek(start)
        for offset in range(0,size,32768):yield stream.read(min(32768,size-offset))


class RangeTests(unittest.TestCase):
    def fixture(self,root):
        data=bytes(range(251))*211;prefix=97
        identity=dict(ExpectedBytes=len(data),RemoteMtime=7,ExpectedServerSha256='a'*64,RemotePath='/synthetic',RelativePath='outputs/synthetic.vcf.gz',SourceReceipt='synthetic',Verification='size_only')
        return data,dict(source_identity=identity,prefix_bytes=prefix,ranges=m.partition(prefix,len(data)))

    def test_partition_and_capacity_bounded_sixteen(self):
        ranges=m.partition(97,9999);self.assertEqual(len(ranges),16)
        self.assertEqual(ranges[0]['offset'],97);self.assertEqual(sum(d['length'] for d in ranges),9902)
        self.assertTrue(all(a['offset']+a['length']==b['offset'] for a,b in zip(ranges,ranges[1:])))
        with self.assertRaises(ValueError):m.partition(0,9999,17)
        with self.assertRaisesRegex(ValueError,'two copies'):m.check_space(Path('.'),100,free=m.base.RESERVE+199)
        m.check_space(Path('.'),100,free=m.base.RESERVE+200)

    def test_ranges_resume_and_ordered_append_exact_without_hashes(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);data,layout=self.fixture(root);prefix=root/'cohort.part';prefix.write_bytes(data[:97]);paths=[]
            for desc in layout['ranges']:
                path=root/('range-%02d.part'%desc['index']);paths.append(path);identity=dict(layout['source_identity'],**desc);sftp=FakeSftp(data)
                if desc['index']==0:
                    path.write_bytes(data[desc['offset']:desc['offset']+5]);helper.atomic(path.with_name(path.name+'.source.json'),identity)
                m.transfer_range(sftp,sftp.blocks,dict(ServerPath='/synthetic'),identity,path,helper,helper.now()+dt.timedelta(minutes=1),threading.Event(),lambda value:None)
                if desc['index']==0:self.assertEqual(sftp.starts,[desc['offset']+5])
                self.assertIsNone(helper.load(path.with_name(path.name+'.completed.json'))['actual_sha256'])
            m.assemble(prefix,layout,paths,len(data),helper,helper.now()+dt.timedelta(minutes=1))
            self.assertEqual(prefix.read_bytes(),data);self.assertTrue(all(p.is_file() for p in paths))

    def test_wrong_identity_and_missing_range_cannot_mutate_prefix(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);data,layout=self.fixture(root);prefix=root/'cohort.part';prefix.write_bytes(data[:97]);path=root/'range-00.part';path.write_bytes(b'ownpartial')
            identity=dict(layout['source_identity'],**layout['ranges'][0]);helper.atomic(path.with_name(path.name+'.source.json'),dict(identity,offset=999))
            sftp=FakeSftp(data)
            with self.assertRaisesRegex(ValueError,'identity'):m.transfer_range(sftp,sftp.blocks,dict(ServerPath='/synthetic'),identity,path,helper,helper.now()+dt.timedelta(minutes=1),threading.Event(),lambda value:None)
            self.assertEqual(path.read_bytes(),b'ownpartial')
            with self.assertRaises((ValueError,FileNotFoundError)):m.assemble(prefix,layout,[root/('missing-%d'%i) for i in range(16)],len(data),helper,helper.now()+dt.timedelta(minutes=1))
            self.assertEqual(prefix.read_bytes(),data[:97]);self.assertFalse((root/'cohort.vcf.gz').exists())

    def test_dry_plan_opens_no_auth(self):
        plan=FinalPairTests().plan();result=m.run(plan,auth='AUTH_MUST_NOT_OPEN')
        self.assertEqual(result['connections'],16);self.assertIn('NO_AUTH',result['status'])


if __name__=='__main__':unittest.main()
