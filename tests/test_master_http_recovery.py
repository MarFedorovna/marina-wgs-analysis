from pathlib import Path
import datetime as dt,hashlib,importlib.util,json,os,sys,tempfile,unittest,zipfile,zlib
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'transfer'));import recover_master_http_segments as r
class Recovery(unittest.TestCase):
 def fixture(self,root):
  base,physical=r.modules();part=root/'cohort.master.tsv.gz.part';base.preallocate(part,32)
  with part.open('r+b') as f:f.write(b'A'*16);f.flush();os.fsync(f.fileno())
  expected=dict(bytes=32,sha256='a'*64);source=dict(artifact=r.ARTIFACT,bytes=32,sha256='a'*64,remote_mtime=7,manifest_sha256='b'*64)
  snap=dict(schema_version=1,artifact=r.ARTIFACT,expected=expected,source=source,chunk_bytes=16,file_id=part.stat().st_ino,completed=[0],missing=[1]);meta=root/'ledger.json';base.atomic(meta,{k:v for k,v in snap.items() if k!='missing'});progress=dict(identity={'public':'control'},segments={});pp=root/'http-progress.json'
  ctx=dict(base=base,snapshot=snap,part=part,meta=meta,physical=physical.physical_allocation,deadline=r.now()+dt.timedelta(minutes=2),progress=progress,progress_path=pp);return ctx
 def archive(self,root,ctx,sub=0,mutate=None,corrupt=False):
  payload=b'B'*8;value=r.spec_for(ctx['snapshot'],1,sub);value['payload_crc32']=f'{zlib.crc32(payload)&0xffffffff:08x}'
  if mutate:mutate(value)
  p=root/f'part-{sub}.zip'
  with zipfile.ZipFile(p,'w',zipfile.ZIP_STORED) as z:z.writestr('block.bin',payload);z.writestr('segment.json',json.dumps(value))
  if corrupt:
   data=bytearray(p.read_bytes());n=data.find(payload);self.assertGreater(n,0);data[n]^=1;p.write_bytes(data)
  return p,dict(name=p.name,bytes=p.stat().st_size,index=1,suboffset=sub,offset=16+sub,payload_bytes=8,payload_crc32=f'{zlib.crc32(payload)&0xffffffff:08x}')
 def test_actual_windows_only_missing_offsets_fsync_checkpoints_and_bitmap(self):
  with tempfile.TemporaryDirectory() as d,patch.object(r,'CHUNK',16),patch.object(r,'SUB',8):
   root=Path(d);ctx=self.fixture(root);p,s=self.archive(root,ctx);self.assertTrue(r.apply_segment(ctx,p,s));self.assertEqual(json.loads(ctx['meta'].read_text())['completed'],[0]);self.assertEqual(len(ctx['progress']['segments']),1)
   p2,s2=self.archive(root,ctx,8);self.assertTrue(r.apply_segment(ctx,p2,s2));self.assertEqual(ctx['part'].read_bytes(),b'A'*16+b'B'*16);self.assertEqual(json.loads(ctx['meta'].read_text())['completed'],[0,1]);self.assertFalse(r.apply_segment(ctx,p,s));self.assertEqual(len(ctx['progress']['segments']),2)
 def test_wrong_offset_source_crc_unowned_refused_before_any_write(self):
  cases=[('offset',lambda v:v.update(offset=v['offset']+1),False,False),('source',lambda v:v['source'].update(remote_mtime=99),False,False),('CRC',None,True,False),('unowned',None,False,True)]
  for label,change,corrupt,unowned in cases:
   with self.subTest(label=label),tempfile.TemporaryDirectory() as d,patch.object(r,'CHUNK',16),patch.object(r,'SUB',8):
    root=Path(d);ctx=self.fixture(root);before=ctx['part'].read_bytes()
    # Keep the immutable snapshot independent from mutated serialized ZIP metadata.
    if label=='source':change=lambda v:v.update(source=dict(v['source'],remote_mtime=99))
    p,s=self.archive(root,ctx,mutate=change,corrupt=corrupt)
    if unowned:ctx['snapshot']['file_id']+=1
    with self.assertRaises((ValueError,zipfile.BadZipFile)):r.apply_segment(ctx,p,s)
    self.assertEqual(ctx['part'].read_bytes(),before);self.assertEqual(ctx['progress']['segments'],{});self.assertEqual(json.loads(ctx['meta'].read_text())['completed'],[0])
 def test_existing_completed_chunk_cannot_be_targeted(self):
  with tempfile.TemporaryDirectory() as d,patch.object(r,'CHUNK',16),patch.object(r,'SUB',8):
   ctx=self.fixture(Path(d));value=r.spec_for(ctx['snapshot'],0,0);value['payload_crc32']='00000000'
   with self.assertRaises(ValueError):r.segment_contract(value,ctx['snapshot'])
if __name__=='__main__':unittest.main()
