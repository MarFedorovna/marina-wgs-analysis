import copy
import datetime as dt
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'transfer'));sys.path.insert(0,str(Path(__file__).parent))
import download_final_vcf_parallel_ranges_v2 as m
import download_final_annotated_vcf as original
import sync_completed_small_results as helper
from test_download_final_annotated_vcf import FinalPairTests
from test_final_vcf_parallel_ranges import RangeTests


class AssemblyResume(unittest.TestCase):
    def fixture(self,root):
        data,layout=RangeTests().fixture(root);prefix=root/'cohort.part';prefix.write_bytes(data[:layout['prefix_bytes']]);paths=[]
        for desc in layout['ranges']:
            path=root/('range-%02d.part'%desc['index']);path.write_bytes(data[desc['offset']:desc['offset']+desc['length']]);paths.append(path)
            helper.atomic(path.with_name(path.name+'.completed.json'),dict(status='COMPLETED_RANGE_SIZE_ONLY',identity=dict(layout['source_identity'],**desc),bytes=desc['length'],local_mtime_ns=path.stat().st_mtime_ns,actual_sha256=None))
        return data,layout,prefix,paths,root/'assembly-phase.json'

    def test_injected_midappend_deadline_resumes_exact_bytes_without_truncate(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);data,layout,prefix,paths,marker=self.fixture(root);deadline=helper.now()+dt.timedelta(minutes=1);calls=[0]
            def injected_now():
                calls[0]+=1;return deadline+dt.timedelta(seconds=1) if calls[0]>3 else helper.now()
            failing=SimpleNamespace(load=helper.load,atomic=helper.atomic,now=injected_now)
            with self.assertRaisesRegex(ValueError,'Assembly deadline'):m.assemble(prefix,layout,paths,len(data),failing,deadline,marker)
            partial=prefix.read_bytes();self.assertGreater(len(partial),layout['prefix_bytes']);self.assertLess(len(partial),len(data));self.assertTrue(marker.exists())
            m.assemble(prefix,layout,paths,len(data),helper,deadline,marker)
            self.assertEqual(prefix.read_bytes(),data);self.assertEqual(prefix.read_bytes()[:len(partial)],partial);self.assertTrue(all(p.exists() for p in paths))

    def test_changed_download_prefix_requires_marker_and_wrong_layout_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            data,layout,prefix,paths,marker=self.fixture(Path(temp));prefix.write_bytes(data[:layout['prefix_bytes']+7]);before=prefix.read_bytes()
            with self.assertRaisesRegex(ValueError,'without authorized assembly marker'):m.assemble(prefix,layout,paths,len(data),helper,helper.now()+dt.timedelta(minutes=1),marker)
            self.assertEqual(prefix.read_bytes(),before)
            helper.atomic(marker,m.assembly_marker(layout,len(data)));wrong=copy.deepcopy(layout);wrong['ranges'][1]['offset']+=1
            with self.assertRaises(ValueError):m.assemble(prefix,wrong,paths,len(data),helper,helper.now()+dt.timedelta(minutes=1),marker)
            self.assertEqual(prefix.read_bytes(),before)

    def test_corrupt_appended_tail_and_missing_completed_range_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            data,layout,prefix,paths,marker=self.fixture(Path(temp));helper.atomic(marker,m.assembly_marker(layout,len(data)));prefix.write_bytes(data[:layout['prefix_bytes']]+b'wrong')
            with self.assertRaisesRegex(ValueError,'differs from completed'):m.assemble(prefix,layout,paths,len(data),helper,helper.now()+dt.timedelta(minutes=1),marker)
            self.assertEqual(prefix.read_bytes(),data[:layout['prefix_bytes']]+b'wrong')
        with tempfile.TemporaryDirectory() as temp:
            data,layout,prefix,paths,marker=self.fixture(Path(temp));helper.atomic(marker,m.assembly_marker(layout,len(data)))
            paths[-1].with_name(paths[-1].name+'.completed.json').unlink()
            with self.assertRaises(FileNotFoundError):m.assemble(prefix,layout,paths,len(data),helper,helper.now()+dt.timedelta(minutes=1),marker)
            self.assertEqual(prefix.read_bytes(),data[:layout['prefix_bytes']])

    def test_new_window_is_transfer_only_original_base_unchanged(self):
        self.assertEqual(original.END.isoformat(),'2026-11-06T03:46:48.046025+00:00')
        self.assertEqual(m.base.END.isoformat(),'2026-11-06T04:50:00+00:00')
        plan=FinalPairTests().plan();plan['deadline_utc']=m.base.END.isoformat();m.base.validate(plan)
        with self.assertRaises(ValueError):original.validate(plan)


if __name__=='__main__':unittest.main()
