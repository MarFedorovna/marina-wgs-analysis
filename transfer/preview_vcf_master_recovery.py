#!/usr/bin/env python3
"""Bounded offline feasibility preview only; never writes/rebuilds a full master."""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import time

CODE=Path(__file__).resolve().parents[1]/'scripts'
PINS={'spliceai_allele_info_codec':'59adb1e99081cca1291f5807e5094e3c638637ef1ff84f066c47bac03f32c761','fast_vcf_evidence_candidate':'0c4532fdfe316441cf1eb762297a63a725450a7a84ba082db1fda92407b42fce'}
def need(ok,message):
    if not ok:raise ValueError(message)
def module(name):
    path=CODE/(name+'.py');need(not path.is_symlink() and hashlib.sha256(path.read_bytes()).hexdigest()==PINS[name],'Pinned portable scientific helper changed')
    spec=importlib.util.spec_from_file_location(name,path);value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value

def info_dict(text):
    result={}
    for entry in text.split(';'):
        if entry in ('','.'):continue
        key,sep,value=entry.partition('=');need(key not in result,'Duplicate INFO key');result[key]=value if sep else None
    return result

def project_record(codec,line,fields,expected_ordinal):
    columns=line.rstrip('\r\n').split('\t');need(len(columns)>=10 and ',' not in columns[4],'Biallelic sample VCF required')
    allele=[columns[0],int(columns[1]),columns[3],columns[4]];info=info_dict(columns[7]);packet=codec.decode_packet(info['WGS_EVIDENCE'])
    need(packet.get('v')==1 and packet.get('allele')==allele and packet.get('scope','cohort')=='cohort' and packet.get('sample') is None,'Cohort packet identity differs')
    predictions,source_json=codec.projection_from_info(info,allele);ix={k:i for i,k in enumerate(fields)};cache={};out=[]
    for row in packet['rows']:
        values=row['values'];need(row['i']==expected_ordinal and len(values)==len(fields) and all(isinstance(v,str) for v in values),'Lossless row ordinal/vector differs');expected_ordinal+=1
        need([values[ix['chrom']],int(values[ix['pos']]),values[ix['ref']],values[ix['alt']]]==allele,'Source allele/vector mismatch')
        gene=values[ix['gene']]
        if gene not in cache:
            extra=codec.enrich_row({k:values[ix[k]] for k in ('chrom','pos','ref','alt','gene')},predictions,source_json);cache[gene]=[extra[k] for k in codec.SPLICE_FIELDS]
        out.append((row['i'],values,values+cache[gene]))
    need(out,'Explicit sample rows required');return out,expected_ordinal

def preview(vcf,manifest_path,record_limit=10):
    need(type(record_limit) is int and 1<=record_limit<=10,'Preview is bounded to one through ten alleles')
    vcf=Path(vcf);manifest_path=Path(manifest_path);need(vcf.is_file() and not vcf.is_symlink() and manifest_path.is_file() and manifest_path.stat().st_size<8*1024**2,'Local completed files required')
    view=json.loads(manifest_path.read_text(encoding='utf-8-sig'));need(view.get('status')=='PASS' and view.get('synthetic') is False and view.get('report_input_kind')=='FULLY_ANNOTATED_VCF_PACKETS_V1','Completed direct views authority required')
    note=json.loads(vcf.with_name(vcf.name+'.completed.json').read_text(encoding='utf-8-sig'));expected=view['source_pins']['vcf']
    need(note.get('Status')=='COMPLETED_SIZE_ONLY' and note['SourceIdentity']['ExpectedBytes']==expected['bytes']==vcf.stat().st_size and note['SourceIdentity']['ExpectedServerSha256']==expected['sha256'] and note['LocalMtimeNs']==vcf.stat().st_mtime_ns,'Completed local VCF source identity differs')
    codec=module('spliceai_allele_info_codec');proof=module('fast_vcf_evidence_candidate');stored=projection=None;records=rows=0;ordinal=1;started=time.monotonic()
    with vcf.open('rb') as raw:
        with gzip.open(raw,'rt',encoding='utf-8',newline='') as stream:
            while True:
                line=stream.readline(64*1024**2+1);need(line and len(line)<=64*1024**2,'Truncated or oversized prefix line');need(raw.tell()<=32*1024**2,'Compressed preview read ceiling32MiB reached')
                if line.startswith('##WGS_EVIDENCE_SCHEMA='):
                    token=re.search(r'Data="([A-Za-z0-9_-]+)"',line);need(token is not None,'Stored schema missing');stored=codec.decode_schema(token[1])
                elif line.startswith('##WGS_SPLICEAI_SCHEMA='):projection=codec.decode_projection_schema(line)
                elif line.startswith('#CHROM'):
                    need(stored is not None and projection is not None,'Both stored and allele projection schemas required');fields=stored['full_fields'];logical=fields+projection['projection_fields'];need(logical==view['logical_full_fields'],'Logical full master fields differ')
                    need(proof.merkle_digest(fields,stored['batches'])==stored['ordered_batch_merkle_root']==view['batch_merkle_root'],'Stored batch descriptor commitment differs')
                    need(stored['row_count']==sum(x['row_count'] for x in stored['batches'])==view['master_transcript_rows'] and sum(x['record_count'] for x in stored['batches'])==view['vcf_records'],'Full declared row/record counts differ');samples=len(line.rstrip().split('\t')[9:])
                elif not line.startswith('#'):
                    projected,ordinal=project_record(codec,line,fields,ordinal);rows+=len(projected);records+=1
                    if records>=record_limit:break
        consumed=raw.tell()
    return dict(status='PASS_BOUNDED_LOCAL_MASTER_RECOVERY_PREVIEW',preview_records=records,preview_rows=rows,compressed_bytes_read_upper_bound=consumed,elapsed_seconds=time.monotonic()-started,stored_fields=len(fields),logical_fields=len(logical),samples=samples,full_source_records=view['vcf_records'],full_source_rows=view['master_transcript_rows'],batch_count=len(stored['batches']),stored_merkle_authority=view['batch_merkle_root'],derived_merkle_authority=view['derived_batch_merkle_root'],full_dataset_scanned=False,full_derived_Merkle_verified=False,full_file_sha_recomputed=False,output_written=False,compressed_master_byte_identity_claimed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--vcf',type=Path,required=True);p.add_argument('--views-manifest',type=Path,required=True);p.add_argument('--preview-records',type=int,default=10);a=p.parse_args();print(json.dumps(preview(a.vcf,a.views_manifest,a.preview_records)))
