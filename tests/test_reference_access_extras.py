"""Safety gates for public reference preparation, with no external network."""
import importlib.util
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
import zipfile

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
import download_external_reference_source as source
import queue_refseq_source_download as queue

class ReferenceDownloadSafety(unittest.TestCase):
 def test_resume_requires_exact_start_end_total(self):
  source.validate_range('bytes 20-99/100',20,100)
  for bad in ['bytes 0-99/100','bytes 20-98/100','bytes 20-99/101','bytes */100',None]:
   with self.subTest(bad=bad),self.assertRaises(ValueError):source.validate_range(bad,20,100)
 def test_symlink_ancestor_rejected(self):
  with tempfile.TemporaryDirectory() as directory:
   base=Path(directory);(base/'real').mkdir()
   try:(base/'link').symlink_to(base/'real',target_is_directory=True)
   except OSError:self.skipTest('Host does not allow symlink creation')
   with self.assertRaises(ValueError):source.safe(base/'link'/'file')
 def test_wheel_requires_all_five_models_and_crc(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'test.whl'
   with zipfile.ZipFile(path,'w') as archive:
    for i in range(1,6):archive.writestr(f'spliceai/models/spliceai{i}.h5',b'synthetic test only')
   self.assertEqual(source.validate_archive(path,'spliceai')['member_count'],5)
   with zipfile.ZipFile(path,'w') as archive:archive.writestr('spliceai/models/spliceai1.h5',b'only one')
   with self.assertRaises(ValueError):source.validate_archive(path,'spliceai')
 def test_tar_links_rejected_before_extraction(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'test.tar.gz'
   with tarfile.open(path,'w:gz') as archive:
    entry=tarfile.TarInfo('bad');entry.type=tarfile.SYMTYPE;entry.linkname='/outside';archive.addfile(entry)
   with self.assertRaises(ValueError):source.validate_archive(path,'genereviews')
 def test_refseq_needs_primary_pass_matching_assembly_version(self):
  good=dict(status='PASS',assembly='GRCh38',vep=dict(cache_version=115))
  self.assertTrue(queue.valid_primary(good))
  for bad in [{},dict(good,status='RUNNING'),dict(good,assembly='GRCh37'),dict(good,vep=dict(cache_version=114))]:self.assertFalse(queue.valid_primary(bad))
 def test_refseq_disk_guard_keeps_analysis_and_unmeasured_expansion(self):
  self.assertEqual(queue.space_requirement(100), (150+128+10)*1024**3+100)
if __name__=='__main__':unittest.main()
