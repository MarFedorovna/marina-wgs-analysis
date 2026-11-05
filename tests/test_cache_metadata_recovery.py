import importlib.util
import io
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(SCRIPTS))
import prepare_annotation_resources as p
import recover_annotation_extracted_stage as r

INFO='species\thomo_sapiens\nassembly\tGRCh38\nsource_assembly\tGRCh38.p14\nsource_gencode\tGENCODE 49\nsource_sift\t6.2.1\nsource_polyphen\t2.2.3\nsource_gnomADe\tv4.1\nsource_gnomADg\tv4.1\nsource_dbSNP\t156\n'


class CacheRecoveryTests(unittest.TestCase):
    def metadata(self,text):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);info=root/'info.txt';info.write_text(text);return p.cache_metadata_records(info,root)
    def test_real_style_missing_mane_is_unreported_not_fabricated_14(self):
        result=self.metadata(INFO+'source_ClinVar\t202502\n')
        self.assertEqual(result['mane_cache']['status'],'source_version_unreported_in_cache_metadata');self.assertIsNone(result['mane_cache']['version'])
        self.assertEqual(result['sift']['status'],'bundled_in_verified_cache');self.assertIn('v4.1',result['gnomad']['version'])
    def test_present_mane_is_reported_verbatim_not_external_summary_version(self):
        result=self.metadata(INFO+'source_MANE\t1.3\n');self.assertEqual(result['mane_cache']['version'],'1.3')
    def test_wrong_nonmane_versions_and_duplicate_metadata_still_fail(self):
        for old,new in [('6.2.1','6.2.2'),('GENCODE 49','GENCODE 48'),('GRCh38.p14','GRCh37'),('v4.1','v3.1.2'),('\t156','\t155')]:
            with self.subTest(old=old),self.assertRaises(ValueError):self.metadata(INFO.replace(old,new))
        with self.assertRaises(ValueError):self.metadata(INFO+'source_sift\t6.2.1\n')
    def fixture(self,root):
        archive=root/'cache.tar.gz';cache=root/'cache';cache.mkdir();files={'homo_sapiens/115_GRCh38/info.txt':INFO.encode(),'homo_sapiens/115_GRCh38/1/data.gz':b'PUBLIC FIXTURE'}
        with tarfile.open(archive,'w:gz') as tar:
            for name,body in files.items():
                item=tarfile.TarInfo(name);item.size=len(body);tar.addfile(item,io.BytesIO(body));path=cache/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(body)
        return archive,cache
    def test_every_extracted_file_matches_archive_exactly(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);archive,cache=self.fixture(root);before={str(x):x.read_bytes() for x in cache.rglob('*') if x.is_file()};result=r.verify_extracted(archive,cache,100000)
            self.assertEqual(result['files'],2);self.assertEqual(result['status'],'PASS_EXACT_ARCHIVE_TO_EXTRACTED_FILES');self.assertEqual(before,{str(x):x.read_bytes() for x in cache.rglob('*') if x.is_file()})
    def test_corruption_or_extra_extracted_file_is_not_resume_pass(self):
        for kind in ['corrupt','extra']:
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as t:
                root=Path(t);archive,cache=self.fixture(root)
                if kind=='corrupt':(cache/'homo_sapiens/115_GRCh38/1/data.gz').write_bytes(b'PUBLIC BROKEN!')
                else:(cache/'extra').write_bytes(b'not in publisher archive')
                with self.assertRaises(ValueError):r.verify_extracted(archive,cache,100000)
    def test_review_only_does_not_read_or_mutate_stage(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);before=list(root.iterdir());value=r.recover(root/'cfg.json',root,root/'missing','a'*64,root/'status',execute=False)
            self.assertEqual(value['status'],'REVIEW_ONLY_NOT_EXECUTED');self.assertEqual(before,list(root.iterdir()))


if __name__=='__main__':unittest.main()
