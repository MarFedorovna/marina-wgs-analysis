import hashlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

spec=importlib.util.spec_from_file_location('revel_source',Path(__file__).resolve().parents[1]/'scripts/download_revel_source.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class Response(io.BytesIO):
    def __init__(self,body,url,size):super().__init__(body);self.url=url;self.status=200;self.headers={'Content-Length':str(size),'Last-Modified':'fixture'}


class RevelTests(unittest.TestCase):
    def fixture(self):
        stream=io.BytesIO()
        with zipfile.ZipFile(stream,'w') as z:z.writestr('revel-v1.3.csv','chr,hg19_pos,grch38_pos,ref,alt,score\n')
        body=stream.getvalue();cfg=dict(bundle_id='revel-source-only-v1',download_authorized_source_only=True,enable_variant_scores=False,
            url='https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip',file='revel.zip',
            expected_bytes=len(body),max_bytes=700000000,publisher_md5=hashlib.md5(body).hexdigest(),deadline_utc='2026-11-06T03:46:48.046025+00:00')
        return cfg,body

    def test_archive_checksum_retained_but_scores_never_enabled(self):
        cfg,body=self.fixture()
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)/cfg['bundle_id']
            def open_(request,timeout):return Response(b'' if request.get_method()=='HEAD' else body,cfg['url'],len(body))
            with patch.object(m.urllib.request,'urlopen',side_effect=open_):result=m.run(cfg,root,True)
            self.assertEqual(result['status'],'DOWNLOADED_NOT_ENABLED_LICENSE_AND_ADAPTER_REVIEW_REQUIRED')
            self.assertEqual(result['publisher_md5'],cfg['publisher_md5']);self.assertFalse(result['patient_scores_enabled']);self.assertFalse(result['extracted'])

    def test_bad_publisher_checksum_preserves_partial_without_ready(self):
        cfg,body=self.fixture();cfg['publisher_md5']='0'*32
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)/cfg['bundle_id']
            def open_(request,timeout):return Response(b'' if request.get_method()=='HEAD' else body,cfg['url'],len(body))
            with patch.object(m.urllib.request,'urlopen',side_effect=open_):
                with self.assertRaises(AssertionError):m.run(cfg,root,True)
            self.assertTrue((root/'revel.zip.partial').exists());self.assertFalse((root/'revel.zip').exists());self.assertFalse((root/'resource-receipt.json').exists())


if __name__=='__main__':unittest.main()
