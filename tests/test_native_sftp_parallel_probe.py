"""Test the bounded probe using in-memory synthetic reads and no real SSH key."""
from pathlib import Path
import sys
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'transfer'))
import native_sftp_parallel_probe as probe

class Remote:
 def __enter__(self):return self
 def __exit__(self,*args):pass
class Client:
 def close(self):pass
class SFTP:
 def __init__(self,size=8*2**20):self.size=size
 def stat(self,path):return SimpleNamespace(st_size=self.size,st_mtime=1234)
 def open(self,path,mode):
  if mode!='rb':raise AssertionError('Probe must be read-only')
  return Remote()
 def close(self):pass

class ParallelProbeTests(unittest.TestCase):
 def test_four_and_eight_connections_read_exactly_four_mib(self):
  for count in [4,8]:
   with self.subTest(count=count):
    requests=[];lock=threading.Lock()
    def blocks(remote,start,size):
     with lock:requests.append((start,size))
     yield b'x'*size
    cfg={'files':[{'remote_path':'/synthetic/result','size_bytes':8*2**20}]}
    with patch.object(probe,'connect',return_value=(Client(),SFTP())),patch.object(probe,'blocks',side_effect=blocks):result=probe.run(cfg,count)
    self.assertEqual(result['total_source_bytes_read'],4*2**20)
    self.assertEqual(sum(size for start,size in requests),4*2**20)
    self.assertEqual(len(requests),count+1)
    intervals=sorted((start,start+size) for start,size in requests)
    self.assertEqual(intervals[0][0],0);self.assertEqual(intervals[-1][1],4*2**20)
    self.assertTrue(all(a[1]==b[0] for a,b in zip(intervals,intervals[1:])))
    self.assertFalse(result['patient_hashes_computed'])
 def test_invalid_connection_count_does_not_connect(self):
  with patch.object(probe,'connect') as connect,self.assertRaises(ValueError):probe.run({},16)
  connect.assert_not_called()
 def test_size_mismatch_prevents_any_source_reads(self):
  cfg={'files':[{'remote_path':'/synthetic/result','size_bytes':8*2**20}]}
  with patch.object(probe,'connect',return_value=(Client(),SFTP(9*2**20))),patch.object(probe,'blocks') as blocks,self.assertRaises(ValueError):probe.run(cfg,4)
  blocks.assert_not_called()
 def test_truncated_baseline_fails_without_parallel_connections(self):
  cfg={'files':[{'remote_path':'/synthetic/result','size_bytes':8*2**20}]}
  with patch.object(probe,'connect',return_value=(Client(),SFTP())) as connect,patch.object(probe,'blocks',return_value=iter([b'too short'])),self.assertRaises(ValueError):probe.run(cfg,4)
  self.assertEqual(connect.call_count,1)
if __name__=='__main__':unittest.main()
