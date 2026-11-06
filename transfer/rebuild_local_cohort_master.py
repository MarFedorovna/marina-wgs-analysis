#!/usr/bin/env python3
"""Offline semantic master recovery from completed VCF packets; default preview."""
import argparse
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import multiprocessing
import os
from pathlib import Path
import re
import shutil
import sys
import time
import uuid
import preview_vcf_master_recovery as preview

MAX_BATCH=64*1024**2
MAX_PENDING=384*1024**2
DEFAULT_WORKERS=max(1,min(6,(os.cpu_count() or 1)-4))

def need(ok,msg):
    if not ok:raise ValueError(msg)
def atomic(path,value):
    path=Path(path);temp=path.with_name(path.name+'.tmp-'+uuid.uuid4().hex)
    with temp.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2);f.flush();os.fsync(f.fileno())
    os.replace(temp,path)
def sha_small(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hash_generated(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
def check_time(end):need(time.monotonic()<end,'Local bounded recovery deadline reached; own parts retained')
def fields_authority(vcf,view,synthetic=False):
    codec=preview.module('spliceai_allele_info_codec');proof=preview.module('fast_vcf_evidence_candidate');stored=projection=None
    with gzip.open(vcf,'rt',encoding='utf-8',newline='') as f:
        while True:
            line=f.readline(MAX_BATCH+1);need(line and len(line)<=MAX_BATCH,'Invalid bounded VCF header')
            if line.startswith('##WGS_EVIDENCE_SCHEMA='):
                m=re.search(r'Data="([A-Za-z0-9_-]+)"',line);need(m is not None,'Stored schema missing');stored=codec.decode_schema(m[1])
            elif line.startswith('##WGS_SPLICEAI_SCHEMA='):projection=codec.decode_projection_schema(line)
            elif line.startswith('#CHROM'):samples=line.rstrip().split('\t')[9:];break
    need(stored and projection and stored.get('codec')==codec.CODEC and stored.get('scope')=='cohort' and stored.get('sample') is None and stored.get('ordinal_in_vcf_order') is True,'Portable ordered row authority required')
    fields=stored['full_fields'];logical=fields+projection['projection_fields']
    need(len(set(logical))==len(logical) and logical==view['logical_full_fields'],'Exact logical field order differs')
    need(proof.merkle_digest(fields,stored['batches'])==stored['ordered_batch_merkle_root']==view['batch_merkle_root'],'Stored header/batch descriptor commitment differs')
    ordinal=1;seen_ids=set()
    for batch in stored['batches']:
        need(type(batch['id']) is int and batch['id']>=0 and batch['id'] not in seen_ids and re.fullmatch('[0-9a-f]{64}',batch['canonical_row_sha256']),'Invalid/duplicate source batch identifier');seen_ids.add(batch['id'])
        need(batch['first_i']==ordinal and batch['row_count']>0 and batch['last_i']==ordinal+batch['row_count']-1 and batch['record_count']>0,'Noncontiguous source batch boundaries');ordinal=batch['last_i']+1
    need(ordinal-1==stored['row_count']==view['master_transcript_rows'] and sum(b['record_count'] for b in stored['batches'])==view['vcf_records'],'Declared full row/record counts differ')
    need(len(samples)==4 and len(set(samples))==4,'Four sample columns required')
    return stored,logical

class HashedOutput:
    def __init__(self,raw):self.raw=raw;self.h=hashlib.sha256()
    def write(self,data):
        n=self.raw.write(data);need(n==len(data),'Short generated gzip write');self.h.update(data);return n
    def flush(self):self.raw.flush()

class CanonicalSink:
    def __init__(self,output):self.output=output;self.h=hashlib.sha256();self.rows=0
    def write(self,text):
        b=text.encode('utf-8');self.output.write(b);self.h.update(b);return len(text)
    def row(self,values):self.writer.writerow(values);self.rows+=1

def process_batch(job):
    check_time(job['end']);codec=preview.module('spliceai_allele_info_codec');proof=preview.module('fast_vcf_evidence_candidate');base=job['batch'];fields=job['fields'];start=base['first_i'];records=0;stored=proof.VectorDigest();target=Path(job['part']);temp=target.with_name(target.name+'.attempt-'+uuid.uuid4().hex)
    need(not target.exists(),'Generated batch already exists without verified checkpoint')
    with temp.open('xb') as raw:
        hashed=HashedOutput(raw)
        with gzip.GzipFile(filename='',fileobj=hashed,mode='wb',compresslevel=1,mtime=0) as output:
            sink=CanonicalSink(output);sink.writer=csv.writer(sink,delimiter='\t',lineterminator='\n')
            for line in job['lines']:
                check_time(job['end']);rows,start=preview.project_record(codec,line,fields,start)
                for ordinal,values,logical in rows:stored.row(values);sink.row(logical)
                records+=1
        raw.flush();os.fsync(raw.fileno());output_sha=hashed.h.hexdigest()
    need(start-1==base['last_i'] and sink.rows==base['row_count'] and records==base['record_count'] and stored.digest.hexdigest()==base['canonical_row_sha256'],'Stored packet batch vectors differ from original authority')
    check_time(job['end']);derived=dict(base,canonical_row_sha256=sink.h.hexdigest());need(not target.exists(),'Generated batch appeared; no overwrite');os.rename(temp,target)
    result=dict(schema_version=1,status='PASS_GENERATED_BATCH',source_identity=job['identity'],source_batch=base,derived_batch=derived,generated_bytes=target.stat().st_size,generated_sha256=output_sha)
    atomic(target.with_suffix(target.suffix+'.json'),result);return result

def source_jobs(vcf,batches):
    with gzip.open(vcf,'rt',encoding='utf-8',newline='') as f:
        for line in f:
            if line.startswith('#CHROM'):break
        for batch in batches:
            lines=[];amount=0
            for _ in range(batch['record_count']):
                line=f.readline(MAX_BATCH+1);need(line and not line.startswith('#') and len(line)<=MAX_BATCH,'Truncated source batch');amount+=len(line.encode('utf-8'));need(amount<=MAX_BATCH,'Source batch exceeds64MiB bounded input');lines.append(line)
            yield batch,lines,amount
        need(not f.read(1),'Extra source VCF record outside closed batch authority')

def recover(vcf,manifest_path,output,workers=DEFAULT_WORKERS,max_hours=8,synthetic=False):
    need(1<=workers<=6 and 0<max_hours<=8,'Bounded local workers/runtime required');end=time.monotonic()+max_hours*3600
    vcf=Path(vcf).resolve(strict=True);manifest_path=Path(manifest_path).resolve(strict=True);view=json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    need(view.get('status')=='PASS' and view.get('synthetic') is synthetic and view.get('report_input_kind')=='FULLY_ANNOTATED_VCF_PACKETS_V1','Completed semantic views control required')
    if not synthetic:
        need(view['code_pins']['codec_module']['sha256']==preview.PINS['spliceai_allele_info_codec'] and view['code_pins']['proof_module']['sha256']==preview.PINS['fast_vcf_evidence_candidate'],'Scientific helper pins differ from completed views')
    expected=view['source_pins']['vcf'];note=json.loads(vcf.with_name(vcf.name+'.completed.json').read_text(encoding='utf-8-sig'))
    need(note.get('Status')=='COMPLETED_SIZE_ONLY' and note['SourceIdentity']['ExpectedBytes']==expected['bytes']==vcf.stat().st_size and note['SourceIdentity']['ExpectedServerSha256']==expected['sha256'] and note['LocalMtimeNs']==vcf.stat().st_mtime_ns,'Completed source identity differs; source VCF is never rehashed')
    schema,logical=fields_authority(vcf,view,synthetic);proof=preview.module('fast_vcf_evidence_candidate')
    output=Path(output).resolve();need(not output.is_relative_to(manifest_path.parent) and not manifest_path.parent.is_relative_to(output) and output!=vcf.parent,'Reconstruction must be separate from source/current dataset');need(not output.is_symlink() and not any((x/'.git').exists() for x in (output,*output.parents)),'New local recovery directory outside Git required');output.mkdir(parents=True,exist_ok=True)
    identity=dict(vcf_bytes=expected['bytes'],expected_source_sha256=expected['sha256'],local_source_mtime_ns=vcf.stat().st_mtime_ns,views_manifest_sha256=sha_small(manifest_path),stored_merkle=schema['ordered_batch_merkle_root'],derived_merkle=view['derived_batch_merkle_root'],producer_sha256=sha_small(__file__),preview_sha256=sha_small(preview.__file__))
    guard=output/'source-identity.json'
    if guard.exists():need(json.loads(guard.read_text())==identity,'Existing reconstruction source/producer identity differs')
    else:atomic(guard,identity)
    final=output/'cohort.master.reconstructed.tsv.gz';receipt=output/'LOCAL-RECONSTRUCTION-RECEIPT.json'
    if final.exists():
        r=json.loads(receipt.read_text());need(r['status']=='PASS_SEMANTIC_MASTER_RECONSTRUCTION' and r['source_identity']==identity and r['output_bytes']==final.stat().st_size and r['output_mtime_ns']==final.stat().st_mtime_ns,'Unrelated final output; no overwrite');return r
    parts=output/'verified-batches';parts.mkdir(exist_ok=True);need(not parts.is_symlink(),'Batch directory symlink refused')
    import msvcrt
    own=(output/'owns.lock').open('a+b');own.write(b'0');own.flush();own.seek(0);msvcrt.locking(own.fileno(),msvcrt.LK_NBLCK,1)
    pending={};pending_bytes=0;results={};pool=None
    def status(label):atomic(output/'status.json',dict(status=label,completed_batches=len(results),batch_count=len(schema['batches']),expected_rows=schema['row_count'],no_source_file_hash=True,workers=workers,updated_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
    def collect(done):
        nonlocal pending_bytes
        for f in done:
            pending_bytes-=pending.pop(f);r=f.result();results[r['source_batch']['id']]=r
        status('REBUILDING_LOCAL_CANONICAL_MASTER')
    try:
        reserve=15*1024**3;need(shutil.disk_usage(output).free>=reserve,'Local15GiB reserve required')
        pool=ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'))
        for batch,lines,amount in source_jobs(vcf,schema['batches']):
            check_time(end);part=parts/(str(batch['id']).zfill(6)+'.tsv.gz');record=part.with_suffix(part.suffix+'.json')
            if part.exists() or record.exists():
                need(part.is_file() and record.is_file() and not part.is_symlink() and not record.is_symlink(),'Incomplete generated checkpoint retained; review before reuse');r=json.loads(record.read_text());need(r['status']=='PASS_GENERATED_BATCH' and r['source_identity']==identity and r['source_batch']==batch and part.stat().st_size==r['generated_bytes'] and hash_generated(part)==r['generated_sha256'],'Generated checkpoint identity/content changed');results[batch['id']]=r;continue
            while pending and (len(pending)>=workers or pending_bytes+amount>MAX_PENDING):collect(wait(pending,return_when=FIRST_COMPLETED,timeout=max(.1,end-time.monotonic()))[0]);check_time(end)
            need(shutil.disk_usage(output).free>=reserve+2*MAX_BATCH,'Local reconstruction reserve reached; checkpoints retained')
            job=dict(batch=batch,lines=lines,fields=schema['full_fields'],identity=identity,part=str(part),end=end);future=pool.submit(process_batch,job);pending[future]=amount;pending_bytes+=amount
        while pending:collect(wait(pending,return_when=FIRST_COMPLETED,timeout=max(.1,end-time.monotonic()))[0]);check_time(end)
        need(vcf.stat().st_size==identity['vcf_bytes'] and vcf.stat().st_mtime_ns==identity['local_source_mtime_ns'],'Source size/mtime changed during local read')
        ordered=[results[b['id']] for b in schema['batches']];derived=[r['derived_batch'] for r in ordered]
        need(proof.merkle_digest(logical,derived)==view['derived_batch_merkle_root'],'Full logical151-field vectors differ from completed views authority')
        need(sum(x['row_count'] for x in derived)==view['master_transcript_rows'] and sum(x['record_count'] for x in derived)==view['vcf_records'],'Full source counts differ')
        assembly=output/'cohort.master.reconstructed.tsv.gz.assembling';need(not assembly.exists(),'Owned interrupted assembly retained; use a new output or review it')
        need(shutil.disk_usage(output).free>=reserve+sum(r['generated_bytes'] for r in ordered),'Assembly requires space for generated parts plus final file')
        status('ASSEMBLING_VERIFIED_GZIP_MEMBERS')
        with assembly.open('xb') as out:
            out.write(gzip.compress(proof.canonical_row_bytes(logical),compresslevel=1,mtime=0))
            for batch in schema['batches']:
                check_time(end)
                with (parts/(str(batch['id']).zfill(6)+'.tsv.gz')).open('rb') as inp:shutil.copyfileobj(inp,out,8*1024**2)
            out.flush();os.fsync(out.fileno())
        check_time(end);need(vcf.stat().st_size==identity['vcf_bytes'] and vcf.stat().st_mtime_ns==identity['local_source_mtime_ns'],'Source changed before local promotion')
        need(not final.exists(),'Final appeared; no overwrite');os.rename(assembly,final)
        result=dict(schema_version=1,status='PASS_SEMANTIC_MASTER_RECONSTRUCTION',source_identity=identity,rows=view['master_transcript_rows'],records=view['vcf_records'],fields=len(logical),derived_merkle_root=view['derived_batch_merkle_root'],stored_merkle_root=schema['ordered_batch_merkle_root'],output_bytes=final.stat().st_size,output_mtime_ns=final.stat().st_mtime_ns,source_file_sha_recomputed=False,compressed_byte_identity_to_server_claimed=False,canonical_TSV_UTF8_LF_vectors_equal=True,automatic_ACMG=False,created_utc=dt.datetime.now(dt.timezone.utc).isoformat());atomic(receipt,result);status(result['status']);return result
    except BaseException:status('FAILED_OWN_CHECKPOINTS_RETAINED');raise
    finally:
        if pool is not None:pool.shutdown(wait=True,cancel_futures=True)
        own.seek(0);msvcrt.locking(own.fileno(),msvcrt.LK_UNLCK,1);own.close()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--vcf',type=Path,required=True);p.add_argument('--views-manifest',type=Path,required=True);p.add_argument('--output',type=Path);p.add_argument('--workers',type=int,default=DEFAULT_WORKERS);p.add_argument('--max-hours',type=float,default=8);p.add_argument('--execute',action='store_true');a=p.parse_args()
    if not a.execute:print(json.dumps(preview.preview(a.vcf,a.views_manifest)))
    else:need(os.name=='nt' and a.output is not None,'Explicit new Windows output required');print(json.dumps(recover(a.vcf,a.views_manifest,a.output,a.workers,a.max_hours)))
