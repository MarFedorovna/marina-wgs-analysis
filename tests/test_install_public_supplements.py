"""Linux-only installer fixtures; no network and no real public downloads."""
import datetime as dt
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
import urllib.error
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('public_installer',Path(__file__).resolve().parents[1]/'scripts/install_public_supplements.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class Response(io.BytesIO):
    def __init__(self,body,url,status=200,headers=None):
        super().__init__(body);self.url=url;self.status=status;self.headers=headers or {}


class InstallerTests(unittest.TestCase):
    def item(self):
        return dict(id='orphadata_genes',file='genes.xml',url='https://www.orphadata.com/data/xml/test.xml',max_bytes=1000,format='xml',release='fixture')

    def test_verified_download_then_idempotent_reuse_no_network(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);item=self.item();body=b'<JDBOR date="2026-06-23"><Disorder/></JDBOR>'
            def open_(request,timeout):return Response(b'' if request.get_method()=='HEAD' else body,item['url'],headers={'Content-Length':str(len(body)),'ETag':'fixture'})
            with patch.object(m.urllib.request,'urlopen',side_effect=open_):
                first=m.download(item,folder,m.utc()+dt.timedelta(seconds=30))
            self.assertEqual(first['status'],'PASS');self.assertEqual(first['bytes'],len(body))
            with patch.object(m.urllib.request,'urlopen',side_effect=AssertionError('Ready source must not redownload')):
                self.assertEqual(m.download(item,folder,m.utc()+dt.timedelta(seconds=30))['sha256'],first['sha256'])

    def test_unsupported_head_uses_documented_get(self):
        with tempfile.TemporaryDirectory() as temp:
            item=self.item();item.update(id='clingen_cspec_public',file='cspec.json',format='json');body=b'{"data":[{"entId":"GN001"}]}'
            def open_(request,timeout):
                if request.get_method()=='HEAD':raise urllib.error.HTTPError(item['url'],400,'HEAD unsupported',{},None)
                return Response(body,item['url'],headers={'Content-Length':str(len(body))})
            with patch.object(m.urllib.request,'urlopen',side_effect=open_):
                result=m.download(item,Path(temp),m.utc()+dt.timedelta(seconds=30))
            self.assertEqual(result['validation']['entities'],1)

    def test_ignored_range_preserves_partial(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);item=self.item();partial=folder/'genes.xml.partial';partial.write_bytes(b'PREFIX')
            with patch.object(m.urllib.request,'urlopen',return_value=Response(b'',item['url'])):
                # Separate Response objects are needed because HEAD closes its stream.
                with patch.object(m.urllib.request,'urlopen',side_effect=lambda *a,**kw:Response(b'WRONG',item['url'])):
                    with self.assertRaises(AssertionError):m.download(item,folder,m.utc()+dt.timedelta(seconds=30))
            self.assertEqual(partial.read_bytes(),b'PREFIX');self.assertFalse((folder/'genes.xml').exists())

    def test_expired_deadline_never_starts_network(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(m.urllib.request,'urlopen',side_effect=AssertionError('No outbound')) as network:
            with self.assertRaises(AssertionError):m.download(self.item(),Path(temp),m.utc()-dt.timedelta(seconds=1))
            network.assert_not_called()


if __name__=='__main__':unittest.main()

