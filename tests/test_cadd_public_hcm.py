import datetime as dt
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
import urllib.request
from unittest.mock import patch

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location('public_cadd',SCRIPTS/'prefetch_cadd_public_hcm.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class Response(io.BytesIO):
    def __init__(self,body,url,headers,status=200):super().__init__(body);self.url=url;self.headers=headers;self.status=status


class CaddTests(unittest.TestCase):
    def config(self):return dict(source_url='https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz',source_bytes=32,tile_bytes=8,max_network_bytes=16)
    def test_actual_http_range_cap_and_cache_reuse(self):
        cfg=self.config();calls=[];body=bytes(range(32))
        def open_(request,timeout):
            if request.get_method()=='HEAD':return Response(b'',cfg['source_url'],{'Content-Length':'32','ETag':'fixture'})
            start,end=map(int,request.headers['Range'][6:].split('-'));calls.append((start,end))
            return Response(body[start:end+1],cfg['source_url'],{'Content-Range':'bytes %d-%d/32'%(start,end),'ETag':'fixture'},206)
        with tempfile.TemporaryDirectory() as t,patch.object(m.urllib.request,'urlopen',side_effect=open_):
            source=m.BoundedSource(cfg,Path(t),m.now()+dt.timedelta(seconds=30))
            self.assertEqual(source.tile(2)[2],body[2:8]);self.assertEqual(source.tile(2)[2],body[2:8])
            self.assertEqual(len(calls),1);source.tile(8);self.assertEqual(source.network_bytes,16)
            with self.assertRaises(AssertionError):source.tile(16)
            self.assertEqual(len(calls),2)

    def test_proxy_preserves_whole_file_length_for_virtual_seeks(self):
        class Source:
            cfg={'source_bytes':32}
            def tile(self,start):return start,min(start+7,31),bytes(range(start,min(start+8,32)))
        server,url=m.proxy(Source())
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'Range':'bytes=8-15'})) as response:
                self.assertEqual(response.headers['Content-Range'],'bytes 8-15/32');self.assertEqual(response.read(),bytes(range(8,16)))
            with urllib.request.urlopen(url) as response:
                self.assertEqual(response.status,200);self.assertEqual(response.headers['Content-Length'],'32');self.assertEqual(response.read(),bytes(range(32)))
        finally:server.shutdown();server.server_close()


if __name__=='__main__':unittest.main()
