import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

MODULE=Path(__file__).resolve().parents[1]/'transfer/assemble_bam_chunks.py'
spec=importlib.util.spec_from_file_location('chunk_assembler',MODULE)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class ChunkBackupTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.downloads=self.root/'downloads';self.destination=self.root/'destination'
        self.downloads.mkdir();self.destination.mkdir()
        self.payload=b'abcdefghijklmnopq'
        self.plan=dict(schema_version=1,export_id='a'*24,verification='size_only',content_hashes_computed=False,
                       download_filename='synthetic.bam',size_bytes=len(self.payload),chunk_bytes=8,chunk_count=3,
                       source_identity={'size_bytes':len(self.payload),'mtime_ns':123},chunks=[])
        for i,offset in enumerate(range(0,len(self.payload),8)):
            content=self.payload[offset:offset+8]
            self.plan['chunks'].append(dict(index=i,offset=offset,size_bytes=len(content),filename=f'public{i}.bin'))
        self.space=patch.object(mod.shutil,'disk_usage',return_value=type('Disk',(),{'free':100*2**30})())
        self.space.start()

    def tearDown(self):self.space.stop();self.tmp.cleanup()

    def write(self,index):
        row=self.plan['chunks'][index];(self.downloads/row['filename']).write_bytes(self.payload[row['offset']:row['offset']+row['size_bytes']])

    def run_assembly(self):return mod.assemble(self.plan,self.downloads,self.destination,True,0)

    def test_resume_and_atomic_exact_reconstruction(self):
        self.write(0);r=self.run_assembly();self.assertEqual(r['completed_chunks'],1)
        self.assertFalse((self.destination/'synthetic.bam').exists())
        self.write(1);self.write(2);r=self.run_assembly()
        self.assertEqual(r['status'],'COMPLETED_SIZE_ONLY')
        self.assertEqual((self.destination/'synthetic.bam').read_bytes(),self.payload)
        self.assertFalse(r['PatientHashesComputed'])

    def test_gap_duplicate_and_traversal_rejected(self):
        for mutation in ['gap','duplicate','traversal']:
            p=copy.deepcopy(self.plan)
            if mutation=='gap':p['chunks'][1]['offset']+=1
            elif mutation=='duplicate':p['chunks'][1]['filename']=p['chunks'][0]['filename']
            else:p['chunks'][0]['filename']='../escape.bin'
            with self.assertRaises(ValueError):mod.validate_manifest(p)

    def test_native_writer_is_not_touched(self):
        native=self.destination/'synthetic.bam.part';native.write_bytes(b'writer')
        with self.assertRaises(ValueError):self.run_assembly()
        self.assertEqual(native.read_bytes(),b'writer')

    def test_browser_partial_not_assembled(self):
        self.write(0);(self.downloads/'public0.bin.crdownload').write_bytes(b'')
        result=self.run_assembly();self.assertEqual(result['completed_chunks'],0)

    def test_wrong_final_never_overwritten(self):
        final=self.destination/'synthetic.bam';final.write_bytes(b'unrelated')
        with self.assertRaises(ValueError):self.run_assembly()
        self.assertEqual(final.read_bytes(),b'unrelated')

    def test_uncommitted_tail_recovers_only_owned_checkpoint(self):
        self.write(0);self.run_assembly()
        partial=self.destination/'synthetic.bam.chunk-assembly.part'
        with partial.open('ab') as stream:stream.write(b'crash')
        self.write(1);self.write(2);self.run_assembly()
        self.assertEqual((self.destination/'synthetic.bam').read_bytes(),self.payload)

    def test_free_space_reserve_enforced(self):
        with patch.object(mod.shutil,'disk_usage',return_value=type('Disk',(),{'free':15*2**30})()):
            with self.assertRaises(ValueError):self.run_assembly()
        self.assertFalse((self.destination/'synthetic.bam').exists())


if __name__=='__main__':unittest.main()
