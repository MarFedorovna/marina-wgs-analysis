from pathlib import Path
import sys,tempfile,unittest,json
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'transfer'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import download_cohort_views_dataset_v2 as v
from test_cohort_dataset_download import plan

class ResumeOnly(unittest.TestCase):
    def fixture(self,root):
        base=v.base_module();p=plan();p.update(destination=str(root/'results'),state_directory=str(root/'state'),pause_file=str(root/'PAUSE'))
        destination=Path(p['destination']);state=Path(p['state_directory']);destination.mkdir();(state/'files').mkdir(parents=True)
        _,files=base.validate(p)
        for name,expected in files.items():
            target,part,meta,note=base.local_paths(destination,state,name);target.parent.mkdir(parents=True,exist_ok=True);base.preallocate(part,expected['bytes'])
            source=dict(artifact=name,bytes=expected['bytes'],sha256=expected['sha256'],remote_mtime=7,manifest_sha256=p['source_provenance']['manifest']['sha256'])
            base.atomic(meta,dict(schema_version=1,artifact=name,expected=expected,source=source,chunk_bytes=base.CHUNK,file_id=part.stat().st_ino,completed=[]))
        return base,p,files
    def test_actual_windows_file_standard_allocation_and_dry_no_auth(self):
        with tempfile.TemporaryDirectory() as d:
            base,p,files=self.fixture(Path(d));proof=v.run(p,auth='MUST_NOT_READ',execute=False)
            self.assertEqual(proof['owned_parts'],len(files));self.assertEqual(proof['additional_allocation_bytes'],0);self.assertGreaterEqual(proof['physical_allocated_bytes'],proof['expected_bytes'])
    def test_missing_part_refused_before_auth(self):
        with tempfile.TemporaryDirectory() as d:
            base,p,files=self.fixture(Path(d));name=next(iter(files));target,part,meta,note=base.local_paths(Path(p['destination']),Path(p['state_directory']),name);part.unlink()
            with self.assertRaisesRegex(ValueError,'existing owned'):v.run(p,auth='MUST_NOT_READ')
    def test_full_logical_but_short_physical_allocation_refused(self):
        with tempfile.TemporaryDirectory() as d:
            base,p,files=self.fixture(Path(d))
            with self.assertRaisesRegex(ValueError,'physical allocation'):v.preflight(base,p,allocation=lambda f:max(0,f.stat().st_size-1))
    def test_changed_source_or_inode_ledger_refused(self):
        with tempfile.TemporaryDirectory() as d:
            base,p,files=self.fixture(Path(d));name=next(iter(files));target,part,meta,note=base.local_paths(Path(p['destination']),Path(p['state_directory']),name);value=base.load(meta);value['source']['manifest_sha256']='f'*64;base.atomic(meta,value)
            with self.assertRaisesRegex(ValueError,'ledger identity'):v.preflight(base,p)
    def test_runtime_floor_and_safe_reasons_do_not_leak(self):
        base=v.base_module();self.assertEqual(v.safe_reason(ValueError('secret user endpoint'),base),'UNCLASSIFIED_ValueError');self.assertEqual(v.safe_reason(ValueError('15GiB runtime reserve reached'),base),'5GiB preallocated-only runtime floor reached')
        with tempfile.TemporaryDirectory() as d:
            base,p,files=self.fixture(Path(d))
            with patch.object(v.shutil,'disk_usage',return_value=type('Space',(),{'free':v.FLOOR-1})()):
                with self.assertRaisesRegex(ValueError,'5GiB'):v.preflight(base,p)

if __name__=='__main__':unittest.main()
