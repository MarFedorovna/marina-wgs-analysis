#!/usr/bin/env python3
"""Fetch fixed public HCM CDS CADD ranges with an actual network-byte cap; no patient input."""
import argparse
import datetime as dt
import fcntl
import hashlib
import http.server
import json
import math
import os
from pathlib import Path
import re
import threading
import time
import urllib.request
from prefetch_myvariant_public_hcm import Reference,load,save

CODE=Path(__file__).resolve().parents[1]


def now():return dt.datetime.now(dt.timezone.utc)


class BoundedSource:
    """Only one pinned public file. Local tile cache enforces source network bytes."""
    def __init__(self,cfg,root,deadline):
        self.cfg=cfg;self.root=root;self.deadline=deadline
        self.network_bytes=0;self.requests=0;self.cache=root/'range-cache';self.cache.mkdir(exist_ok=True)
        self.cached_bytes=sum(p.stat().st_size for p in self.cache.glob('*.part'))
        self.network_bytes=self.cached_bytes
        self.lock=threading.Lock()
        with urllib.request.urlopen(urllib.request.Request(cfg['source_url'],method='HEAD'),timeout=20) as response:
            assert response.url==cfg['source_url'] and int(response.headers['Content-Length'])==cfg['source_bytes']
            self.identity=dict(url=response.url,bytes=cfg['source_bytes'],etag=response.headers.get('ETag'),last_modified=response.headers.get('Last-Modified'))
        old=load(root/'source-identity.json');assert not old or old==self.identity,'CADD publisher identity changed'
        save(root/'source-identity.json',self.identity)

    def tile(self,start):
        with self.lock:return self._tile(start)

    def _tile(self,start):
        assert now()<self.deadline,'Original budget/runtime deadline reached'
        size=self.cfg['tile_bytes'];left=start//size*size;right=min(left+size,self.cfg['source_bytes'])-1
        path=self.cache/('%012d.part'%left)
        if path.exists():
            data=path.read_bytes();assert len(data)==right-left+1
        else:
            assert self.network_bytes+(right-left+1)<=self.cfg['max_network_bytes'],'CADD actual HTTP byte ceiling reached'
            headers={'Range':'bytes=%d-%d'%(left,right),'Accept-Encoding':'identity'}
            with urllib.request.urlopen(urllib.request.Request(self.cfg['source_url'],headers=headers),timeout=20) as response:
                assert response.status==206 and response.url==self.cfg['source_url']
                assert response.headers.get('Content-Range')=='bytes %d-%d/%d'%(left,right,self.cfg['source_bytes'])
                if self.identity['etag']:assert response.headers.get('ETag')==self.identity['etag']
                data=response.read(size+1)
            assert len(data)==right-left+1,'Unexpected bounded HTTP range length'
            self.network_bytes+=len(data);self.requests+=1;path.write_bytes(data)
        return start,right,data[start-left:]


def proxy(source):
    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self,*args):pass
        def do_HEAD(self):
            self.send_response(200);self.send_header('Content-Length',str(source.cfg['source_bytes']));self.send_header('Accept-Ranges','bytes');self.end_headers()
        def do_GET(self):
            try:
                range_header=self.headers.get('Range')
                match=re.fullmatch(r'bytes=(\d+)-(\d*)',range_header or 'bytes=0-')
                if not match:raise ValueError('Only bounded byte range reads are supported')
                start=int(match[1]);assert 0<=start<source.cfg['source_bytes']
                end=min(int(match[2]) if match[2] else source.cfg['source_bytes']-1,source.cfg['source_bytes']-1)
                self.send_response(206 if range_header else 200);self.send_header('Content-Type','application/octet-stream');self.send_header('Accept-Ranges','bytes')
                if range_header:self.send_header('Content-Range','bytes %d-%d/%d'%(start,end,source.cfg['source_bytes']))
                self.send_header('Content-Length',str(end-start+1));self.end_headers()
                position=start
                # HTSlib needs the actual file/response length for virtual seeks. Stream
                # capped tiles only while its reader keeps this connection open.
                while position<=end:
                    left,right,data=source.tile(position);data=data[:end-position+1]
                    self.wfile.write(data);self.wfile.flush();position+=len(data)
            except (BrokenPipeError,ConnectionResetError):pass
            except BaseException as error:
                source.failure=str(error);self.send_error(509,'Bounded source read stopped')
    server=http.server.ThreadingHTTPServer(('127.0.0.1',0),Handler)
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    return server,'http://127.0.0.1:%d/scores.tsv.gz'%server.server_port


def run(cfg,root,execute=False):
    assert cfg['patient_variant_inputs_allowed'] is False and root.is_absolute() and root.name==cfg['bundle_id'] and not root.is_symlink()
    assert cfg['max_network_bytes']+cfg['max_output_bytes']+cfg['index_bytes']<=cfg['max_resource_bytes']<=300000000
    assert cfg['source_url']=='https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz'
    status=root/'status.json'
    if not execute:return load(status) or dict(status='NOT_STARTED_PUBLIC_CDS_ONLY')
    assert load(Path(cfg['purpose_approval'])).get('status')=='USER_CONFIRMED_NONCOMMERCIAL_RESEARCH'
    os.umask(0o077);root.mkdir(mode=0o700,parents=True,exist_ok=True)
    with (root/'.prefetch.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return dict(status='ALREADY_RUNNING')
        import pysam
        panel=Path(cfg['panel']);panel_manifest=load(panel/'panel-manifest.json');bed=panel/cfg.get('bed_name','hcm-cds-union-cadd.bed')
        assert panel_manifest['no_patient_selection'] is True and panel_manifest['assembly']=='GRCh38'
        assert hashlib.sha256(bed.read_bytes()).hexdigest()==panel_manifest[cfg.get('bed_sha_field','bed_sha256')]
        intervals=[(v[0],int(v[1]),int(v[2])) for v in (line.split('\t') for line in bed.read_text().splitlines())]
        deadline=min(dt.datetime.fromisoformat(cfg['deadline_utc']),now()+dt.timedelta(minutes=cfg['max_minutes']))
        index=root/'whole_genome_SNVs.tsv.gz.tbi'
        if not index.exists():
            if cfg.get('reuse_index'):
                data=Path(cfg['reuse_index']).read_bytes()
            else:
                with urllib.request.urlopen(cfg['index_source_url'],timeout=20) as response:data=response.read(cfg['index_bytes']+1)
            assert len(data)==cfg['index_bytes'] and hashlib.md5(data).hexdigest()==cfg['index_publisher_md5'];index.write_bytes(data)
        assert index.stat().st_size==cfg['index_bytes'] and hashlib.md5(index.read_bytes()).hexdigest()==cfg['index_publisher_md5']
        source=BoundedSource(cfg,root,deadline);server,url=proxy(source);reference=Reference(cfg['reference_fasta'])
        rows=0;completed=0;output_bytes=0;network_error=None;final=root/'cadd-public-hcm.tsv';partial=root/'cadd-public-hcm.tsv.partial'
        assert not final.exists() and not partial.exists(),'Preserve existing subset; inspect receipt before choosing a new output version'
        try:
            tabix=pysam.TabixFile(url,index=str(index));headers=list(tabix.header)
            assert any('CADD GRCh38-v1.7' in h for h in headers),'CADD header/version differs'
            controls=list(tabix.fetch('11',47332145,47332146))
            assert len(controls)==3 and {v.split('\t')[3] for v in controls}=={'A','C','G'} and all(v.split('\t')[2]=='T' for v in controls)
            assert reference.base('11',47332146)=='T'
            save(root/'public-control-receipt.json',dict(status='PASS_PUBLIC_REFERENCE_CONTROL_ONLY',source_headers=headers,rows=controls,index_md5=cfg['index_publisher_md5'],pysam=pysam.__version__,no_patient_alleles_used=True))
            with partial.open('x') as out:
                heading='chrom\tpos\tref\talt\tCADD_RAW\tCADD_PHRED\n';out.write(heading);output_bytes=len(heading)
                for chrom,start,end in intervals:
                    assert now()<deadline,'Original budget/runtime deadline reached'
                    for line in tabix.fetch(chrom,start,end):
                        v=line.split('\t');assert len(v)==6 and v[0]==chrom and start<int(v[1])<=end
                        assert v[2] in 'ACGT' and v[3] in 'ACGT' and v[2]!=v[3]
                        assert reference.base(chrom,int(v[1]))==v[2],'Public GRCh38 reference mismatch'
                        assert math.isfinite(float(v[4])) and math.isfinite(float(v[5])) and float(v[5])>=0
                        record='chr'+line+'\n';output_bytes+=len(record);assert output_bytes<=cfg['max_output_bytes'],'CADD output ceiling reached'
                        out.write(record);rows+=1
                    completed+=1
                    save(status,dict(status='PREFETCHING_PUBLIC_HCM_CDS_NOT_PATIENT_ANNOTATION',pid=os.getpid(),completed_intervals=completed,total_intervals=len(intervals),rows=rows,output_bytes=output_bytes,actual_http_range_bytes=source.network_bytes,index_bytes=cfg['index_bytes'],range_requests=source.requests,updated_utc=now().isoformat()))
            tabix.close();partial.replace(final)
        except BaseException as error:
            network_error=getattr(source,'failure',str(error)) or str(error)
        finally:server.shutdown();server.server_close()
        result=dict(status='COMPLETE_PUBLIC_HCM_CDS_SOURCE_ADAPTER_CANDIDATE' if completed==len(intervals) and not network_error else 'PARTIAL_PUBLIC_CDS_SOURCE_NOT_ENABLED',
            source=cfg,source_identity=source.identity,public_panel_sha256=hashlib.sha256((panel/'panel-manifest.json').read_bytes()).hexdigest(),
            completed_intervals=completed,total_intervals=len(intervals),rows=rows,output_bytes=output_bytes,actual_http_range_bytes=source.network_bytes,
            no_patient_input=True,annotation_binding_updated=False,automatic_acmg=False,error=network_error,updated_utc=now().isoformat())
        if final.exists():result['adapter']=dict(path=str(final),bytes=final.stat().st_size,sha256=hashlib.sha256(final.read_bytes()).hexdigest(),format='variant_evidence_tsv',scope='PUBLIC_HCM_CDS_SNV_ONLY; no genome-wide/indel/noncoding claim')
        save(root/'source-receipt.json',result);save(status,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=CODE/'config/cadd-public-hcm-v1.json');p.add_argument('--root',type=Path,default=Path('/srv/marina/resources/cadd-public-hcm-v1'));p.add_argument('--execute',action='store_true');a=p.parse_args();print(json.dumps(run(load(a.config),a.root,a.execute),indent=2))
