#!/usr/bin/env python3
"""One packet pass, bounded batch workers, immutable research views.

This NEW producer consumes a direct packet authority, never fabricates legacy
06/07 receipts and never reruns VEP. Execution requires an exact reviewed rule.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
import csv
import datetime as dt
import gzip
import hashlib
import html
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import sys
import time
import uuid

END = dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
CONTRACT = 'ordered_batch_canonical_TSV_utf8_LF_sha256_v1'
PRODUCER = 'DIRECT_FROZEN_ANNOTATION_ROW_V1_NO_LEGACY_06_07_RECEIPTS'
RANK_FIELDS = ['phenotype_context_rank','target_gene_context_rank','clinvar_review_rank','gene_validity_rank',
 'consequence_rank','population_rank','predictors_damaging','AlphaMissense_research_rank','mane_select',
 'variant_id','gene','transcript','sample']
SPLICE_FIELDS = ['spliceai_status','spliceai_gene','spliceai_DS_AG','spliceai_DS_AL','spliceai_DS_DG','spliceai_DS_DL',
 'spliceai_DP_AG','spliceai_DP_AL','spliceai_DP_DG','spliceai_DP_DL','spliceai_prediction_json','spliceai_source_json','spliceai_transcript_scope']
BGZF_EOF = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
G = {}


def need(ok, message):
    if not ok: raise ValueError(message)


def load(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def js(value): return json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':'))
def now(): return dt.datetime.now(dt.timezone.utc)
def guard(deadline): need(now() < min(END,dt.datetime.fromisoformat(deadline)), 'Original deadline reached; no automatic restart')


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''): digest.update(block)
    return digest.hexdigest()


def pin(path):
    path=Path(path).resolve(strict=True)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))


def verify(item,content=True):
    path=Path(item['path'])
    need(not path.is_symlink() and path.is_file() and path.stat().st_size==item['bytes'],'Pinned input size/path changed')
    need(re.fullmatch('[0-9a-f]{64}',item['sha256']) is not None,'Invalid expected SHA256')
    if content: need(sha(path)==item['sha256'],'Pinned input content changed')
    return path


def module(item,name):
    path=verify(item);spec=importlib.util.spec_from_file_location(name,path)
    if str(path.parent) not in sys.path:sys.path.insert(0,str(path.parent))
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value


def info_dict(value):
    output={}
    for entry in value.split(';'):
        if not entry or entry=='.':continue
        key,sep,text=entry.partition('=');need(key not in output,'Duplicate INFO key')
        output[key]=text if sep else None
    return output


def info_text(value):return ';'.join(k if v is None else k+'='+str(v) for k,v in value.items()) or '.'


def init_worker(rule,fields,samples,stored_fields=None):
    G.clear();G.update(rule=rule,fields=fields,samples=samples,ix={f:i for i,f in enumerate(fields)})
    G['stored_fields']=stored_fields or fields
    G['engine']=module(rule['view_policy'],'fastviews_frozen_policy')
    G['rank']=module(rule['ranking_module'],'fastviews_frozen_rank')
    G['proof']=module(rule['proof_module'],'fastviews_direct_proof')
    G['codec']=module(rule['codec_module'],'fastviews_final_codec')
    if stored_fields and stored_fields!=fields:G['projection']=G['codec']
    if hasattr(os,'sched_getaffinity'):
        allowed=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,set(allowed[:rule['workers']]))
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    except ImportError:pass
    if hasattr(signal,'SIGALRM'):
        def expired(signum,frame):raise TimeoutError('Original deadline reached in bounded report worker')
        signal.signal(signal.SIGALRM,expired)
        seconds=(min(END,dt.datetime.fromisoformat(rule['deadline_utc']))-now()).total_seconds()
        signal.alarm(max(1,math.ceil(seconds)))


def rank_vector(values):
    ix=G['ix'];return G['rank'].rank_key({f:values[ix[f]] for f in RANK_FIELDS})


def predictions(values):
    """Actual model values only; no missing-to-zero or clinical threshold."""
    ix=G['ix']
    if 'spliceai_prediction_json' not in ix:return []
    raw=values[ix['spliceai_prediction_json']]
    result=json.loads(raw) if raw else []
    need(isinstance(result,list),'SpliceAI evidence must be a list')
    for row in result:
        need(isinstance(row,dict),'Malformed SpliceAI prediction row')
    return result


def scored_max(rows):
    scores=[]
    for row in rows:
        if row.get('status')!='SCORED_EXACT_GRCH38_ALLELE_GENE':continue
        for field in ('DS_AG','DS_AL','DS_DG','DS_DL'):
            value=row.get(field,'')
            if value=='':continue
            number=float(value);need(math.isfinite(number) and 0<=number<=1,'Invalid actual SpliceAI delta score')
            scores.append(number)
    return max(scores) if scores else None


def keep_top(items,value,key):
    items.append(value);items.sort(key=key);del items[20:]


def category_name(term):return re.sub('[^A-Za-z0-9_-]','_',term)[:80]+'-'+hashlib.sha256(term.encode()).hexdigest()[:8]


def process_batch(job):
    """Decode each source packet once; all rows of an allele stay together."""
    import pysam
    rule,fields,samples=G['rule'],G['fields'],G['samples'];ix=G['ix'];engine=G['engine'];codec=G['codec'];proof=G['proof']
    guard(rule['deadline_utc']);folder=Path(job['folder']);folder.mkdir()
    labels={s:'sample_'+str(i+1).zfill(4) for i,s in enumerate(samples)}
    scopes=['cohort','shared_ge2','all4',*samples]
    summary_fields=list(dict.fromkeys([*engine.SUMMARY,*fields,'scope_raw_member','scope_qc_member','scope_candidate_member']))
    splice_fields=[*summary_fields,'spliceai_max_delta','spliceai_all_gene_predictions_json','spliceai_prediction_ranking_scope']
    table_fields={'masters/cohort.master.tsv.gz':fields,'tables/allele-membership-matrix.tsv.gz':engine.SUMMARY[:24]}
    tops={s:[] for s in scopes};splice_tops={s:[] for s in ['cohort',*samples]}
    counts={s:Counter() for s in scopes};terms_counts={s:defaultdict(Counter) for s in scopes}
    genes=Counter();sample_rows=Counter();sample_records=Counter();sample_digests={s:proof.VectorDigest() for s in samples}
    digest=proof.VectorDigest();prefix_digest=proof.VectorDigest();derived_digest=proof.VectorDigest();ordinal=job['proof']['first_i'];records=0
    writers={};raw_streams={};sample_payload_hash={s:hashlib.sha256() for s in samples};schema_prefix=job['prefix_count']
    with ExitStack() as stack:
        def writer(relative,columns):
            if relative not in writers:
                path=folder/relative;path.parent.mkdir(parents=True,exist_ok=True)
                stream=stack.enter_context(gzip.open(path,'xt',encoding='utf-8',newline='',compresslevel=3))
                writers[relative]=csv.writer(stream,delimiter='\t',lineterminator='\n');table_fields[relative]=columns
            return writers[relative]
        def emit(relative,columns,row):
            writer(relative,columns).writerow([row.get(f,'') for f in columns])
        def master(relative,values):writer(relative,fields).writerow(values)
        for line in job['lines']:
            guard(rule['deadline_utc']);columns=line.rstrip('\r\n').split('\t');need(len(columns)==9+len(samples),'VCF sample columns changed')
            if records%1000==0:need(shutil.disk_usage(folder).free>=rule.get('min_free_bytes',32*1024**3),'VM report disk reserve reached; partials retained')
            allele=[columns[0],int(columns[1]),columns[3],columns[4]];need(',' not in columns[4],'Biallelic source required')
            info=info_dict(columns[7]);need('WGS_EVIDENCE' in info,'Missing evidence packet')
            packet=codec.decode_packet(info['WGS_EVIDENCE'])
            need(packet.get('v')==1 and packet.get('allele')==allele and packet.get('scope','cohort')=='cohort' and packet.get('sample') is None,'Packet allele/scope changed')
            fmt=columns[8].split(':');need(len(fmt)==len(set(fmt)) and 'GT' in fmt and 'FT' in fmt,'Explicit unique GT/FT fields required')
            calls={}
            sequence=columns[4] not in ('.','*') and not columns[4].startswith('<') and '[' not in columns[4] and ']' not in columns[4]
            for index,sample in enumerate(samples):
                raw=dict(zip(fmt,columns[9+index].split(':')))
                calls[sample]={f:raw.get(f,'.') for f in ('GT','FT','DP','AD','GQ')}
                calls[sample].update(site_filter=columns[6],flags=engine.call_flags(raw.get('GT','.'),raw.get('FT','.'),columns[6],sequence))
            best_by={};eligible={};representative=None;terms=set();gene_set=set();seen=set();prediction_map={}
            rows=packet.get('rows');need(isinstance(rows,list) and rows,'Every allele requires explicit sample rows, including no-CSQ placeholders')
            projected_rows=[];projection_cache={}
            if 'projection' in G:source_predictions,source_json=G['projection'].projection_from_info(info,allele)
            for row in rows:
                values=row.get('values');need(row.get('i')==ordinal and len(values)==len(G['stored_fields']) and all(isinstance(v,str) for v in values),'Noncontiguous ordinal or malformed complete field vector')
                ordinal+=1;digest.row(values);prefix_digest.row(values[:schema_prefix])
                if 'projection' in G:
                    gene=values[ix['gene']]
                    if gene not in projection_cache:
                        small={f:values[ix[f]] for f in ('chrom','pos','ref','alt','gene')}
                        extra=G['projection'].enrich_row(small,source_predictions,source_json)
                        projection_cache[gene]=[extra[f] for f in SPLICE_FIELDS]
                    values=values+projection_cache[gene]
                projected_rows.append(dict(i=row['i'],values=values));derived_digest.row(values);master('masters/cohort.master.tsv.gz',values)
                sample=values[ix['sample']];need(sample in calls,'Unknown packet sample');seen.add(sample);call=calls[sample]
                need([values[ix['chrom']],int(values[ix['pos']]),values[ix['ref']],values[ix['alt']]]==allele,'Row allele identity differs')
                for field in ('GT','FT','DP','AD','GQ','site_filter'):
                    need(values[ix[field]]==str(call[field]),'Packet genotype/quality differs from VCF')
                need(values[ix['variant_id']]==G['rank'].variant_key(*allele),'Frozen variant identity differs')
                key=rank_vector(values)
                if sample not in best_by or key<best_by[sample][0]:best_by[sample]=(key,values)
                if representative is None or (call['flags']['raw'] and not calls[representative[1][ix['sample']]]['flags']['raw']) or (call['flags']['raw']==calls[representative[1][ix['sample']]]['flags']['raw'] and key<representative[0]):representative=(key,values)
                state=values[ix['candidate_eligible']];need(state in ('0','1'),'Unknown frozen eligibility')
                if state=='1':
                    need(call['flags']['qc'],'Frozen eligibility contradicts complete PASS ALT')
                    if sample not in eligible or key<eligible[sample][0]:eligible[sample]=(key,values)
                if call['flags']['raw']:
                    master('masters/'+labels[sample]+'.master.tsv.gz',values);sample_rows[sample]+=1;sample_digests[sample].row(values)
                terms.update(t for t in re.split('[&,;|]',values[ix['consequence']]) if t)
                if values[ix['gene']]:gene_set.add(values[ix['gene']])
                # JSON repeats across transcript/sample rows: parse a distinct token once per allele.
                if 'spliceai_prediction_json' in ix:
                    token=values[ix['spliceai_prediction_json']]
                    if token and token not in prediction_map:prediction_map[token]=predictions(values)
            need(seen==set(samples),'Missing explicit sample row; missing is not reference')
            terms=terms or {'unannotated'}
            if terms & engine.PLOF:terms.add('pLoF_candidates')
            best=min(eligible.values(),key=lambda v:v[0]) if eligible else representative
            rowdict=dict(zip(fields,best[1]));summary=dict(allele_id=engine.allele_id(*allele),chrom=allele[0],pos=allele[1],ref=allele[2],alt=allele[3],
                **engine.membership(calls,samples,eligible),cohort_AC_source=info.get('AC','.'),cohort_AN_source=info.get('AN','.'),cohort_AF_source=info.get('AF','.'),
                annotation_transcript_rows=len(rows),union_consequences_json=js(sorted(terms)),representative_sample=rowdict['sample'],
                functional_membership_basis='ANY_SOURCE_TRANSCRIPT; distinct alleles; overlapping categories',ACMG_class='NOT_ASSIGNED',
                ACMG_criteria_evidence_status='NOT_ASSESSED_MANUAL_EVIDENCE_REQUIRED',ACMG_provisional_review='PENDING_MANUAL_REVIEW',
                all_transcript_evidence='masters/cohort.master.tsv.gz',evidence_join='exact allele + sample; complete transcript vectors retained',
                interpretation_status='RESEARCH_VIEWS_NOT_CLINICAL_BURDEN_ACMG_OR_MECHANISM')
            summary.update(rowdict);emit('tables/allele-membership-matrix.tsv.gz',engine.SUMMARY[:24],summary)
            for gene in gene_set:
                for sample,call in calls.items():
                    if call['flags']['raw']:genes[(gene,sample,'raw')]+=1
                    if call['flags']['qc']:genes[(gene,sample,'qc')]+=1
            for scope in scopes:
                raw,qc,candidate=engine.scope_members(scope,calls,eligible)
                scoped=dict(summary)
                chosen=eligible.get(scope,best_by.get(scope)) if scope in samples else best
                if scope in samples:scoped.update(dict(zip(fields,chosen[1])));scoped['representative_sample']=scope
                scoped.update(scope_raw_member=int(raw),scope_qc_member=int(qc),scope_candidate_member=int(candidate))
                label=labels.get(scope,scope)
                if raw:
                    emit('tables/'+label+'.functional.tsv.gz',summary_fields,scoped)
                    for term in terms:emit('tables/'+label+'/'+category_name(term)+'.tsv.gz',summary_fields,scoped)
                for level,present in [('raw',raw),('qc',qc),('eligible',candidate)]:
                    if present:
                        counts[scope][level+'_alleles']+=1
                        for term in terms:terms_counts[scope][term][level]+=1
                if candidate:
                    chosen=eligible[scope] if scope in samples else best
                    selected=dict(scoped);selected.update(dict(zip(fields,chosen[1])));selected['representative_sample']=selected['sample']
                    keep_top(tops[scope],(chosen[0],selected),lambda v:v[0])
            # Preserve raw observed ALT records, but prune every private packet/readable tag.
            for index,sample in enumerate(samples):
                if not calls[sample]['flags']['raw']:continue
                logical_packet=dict(packet,rows=projected_rows)
                selected=codec.select_sample(logical_packet,fields,sample)
                selected={k:v for k,v in selected.items() if k in {'v','allele','rows','evidence','gene_context_keys','scope','sample'}}
                need(selected['rows'] and all(v['values'][ix['sample']]==sample for v in selected['rows']),'Sample packet contains another person')
                sample_info={k:v for k,v in info.items() if k not in codec.INFO_IDS}
                sample_info.update(codec.rebuild_readable(selected,fields));sample_info['WGS_EVIDENCE']=codec.encode_packet(selected)
                subset=list(columns[:9])+[columns[9+index]];subset[7]=engine.derived_info(sample_info,calls[sample]['flags'])
                if sample not in raw_streams:
                    path=folder/'vcf'/(labels[sample]+'.body.bgz');path.parent.mkdir(exist_ok=True)
                    raw_streams[sample]=stack.enter_context(pysam.BGZFile(str(path),'wb'))
                payload=('\t'.join(subset)+'\n').encode();sample_payload_hash[sample].update(payload)
                raw_streams[sample].write(payload);sample_records[sample]+=1
            all_predictions={js(p):p for rows_ in prediction_map.values() for p in rows_}
            score=scored_max(all_predictions.values())
            if score is not None:
                for scope in ['cohort',*samples]:
                    qc=any(c['flags']['qc'] for c in calls.values()) if scope=='cohort' else calls[scope]['flags']['qc']
                    if not qc:continue
                    chosen=min((best_by[s] for s in samples if calls[s]['flags']['qc']),key=lambda v:v[0]) if scope=='cohort' else best_by[scope]
                    record=dict(summary);record.update(dict(zip(fields,chosen[1])));record.update(spliceai_max_delta=format(score,'.12g'),
                        spliceai_all_gene_predictions_json=js(list(all_predictions.values())),
                        spliceai_prediction_ranking_scope='PUBLIC39_GENE_MODELS; GENCODE24_GENE_UNION_NOT_ENST49_MAPPING; prediction only, no automatic ACMG')
                    label=labels.get(scope,scope);emit('tables/spliceai/'+label+'.scored.tsv.gz',splice_fields,record)
                    keep_top(splice_tops[scope],((-score,engine.allele_id(*allele)),record),lambda v:v[0])
            records+=1
    actual=dict(id=job['proof']['id'],first_i=job['proof']['first_i'],last_i=ordinal-1,row_count=ordinal-job['proof']['first_i'],record_count=records,canonical_row_sha256=digest.digest.hexdigest())
    need(proof.portable_batch(actual)==proof.portable_batch(job['proof']),'Final packet batch proof differs')
    prefix=dict(actual,canonical_row_sha256=prefix_digest.digest.hexdigest())
    need(proof.portable_batch(prefix)==proof.portable_batch(job['prefix_proof']),'Frozen original field prefix proof differs')
    derived=dict(actual,canonical_row_sha256=derived_digest.digest.hexdigest())
    persisted_parts={}
    for sample in raw_streams:
        path=folder/'vcf'/(labels[sample]+'.body.bgz');observed=hashlib.sha256()
        with gzip.open(path,'rb') as stream:
            for block in iter(lambda:stream.read(4*1024**2),b''):observed.update(block)
        need(observed.hexdigest()==sample_payload_hash[sample].hexdigest(),'Persisted sample payload differs from every scope-checked emitted byte')
        persisted_parts[sample]=dict(bytes=path.stat().st_size,sha256=sha(path),uncompressed_sha256=observed.hexdigest())
    return dict(folder=str(folder),proof=actual,prefix_proof=prefix,derived_proof=derived,table_fields=table_fields,tops=tops,splice_tops=splice_tops,
        counts={s:dict(v) for s,v in counts.items()},terms_counts={s:{k:dict(v) for k,v in value.items()} for s,value in terms_counts.items()},
        genes=[[g,s,k,n] for (g,s,k),n in genes.items()],sample_rows=dict(sample_rows),sample_records=dict(sample_records),
        sample_digests={s:sample_digests[s].digest.hexdigest() for s in samples},persisted_sample_parts=persisted_parts)


def validate_schema(schema,prefix,proof):
    fields=schema.get('full_fields',[]);original=prefix.get('full_fields',[])
    need(schema.get('scope')=='cohort' and schema.get('sample') is None and schema.get('codec')=='json-utf8+zlib+base64url-nopad','Wrong packet codec/scope')
    need(prefix.get('producer_kind')==PRODUCER and schema.get('ordinal_in_vcf_order') is True and prefix.get('ordinal_in_vcf_order') is True,'New direct ordinal authority required; legacy receipts are not accepted')
    need(fields[:len(original)]==original and len(fields)==len(set(fields)),'Frozen original field prefix differs')
    need(schema.get('vector_digest_contract')==prefix.get('vector_digest_contract')==CONTRACT,'Unknown vector proof contract')
    need(proof.merkle_digest(fields,schema['batches'])==schema['ordered_batch_merkle_root'],'Final ordered batch commitment differs')
    need(proof.merkle_digest(original,prefix['batches'])==prefix['ordered_batch_merkle_root'],'Core ordered batch commitment differs')
    need(len(schema['batches'])==len(prefix['batches'])>0,'Core/final batch count differs')
    ordinal=1
    for final,base in zip(schema['batches'],prefix['batches']):
        need(all(final[k]==base[k] for k in ('id','first_i','last_i','row_count','record_count')),'Core/final batch boundaries differ')
        need(final['first_i']==ordinal and final['row_count']>0 and final['last_i']==ordinal+final['row_count']-1 and final['record_count']>0,'Noncontiguous batch proof')
        ordinal=final['last_i']+1
    need(ordinal-1==schema['row_count']==prefix['row_count'],'Proof total row count differs')
    required=set(RANK_FIELDS+['chrom','pos','ref','alt','GT','FT','DP','AD','GQ','site_filter','candidate_eligible','consequence'])
    need(required<=set(fields),'Required frozen source fields missing')
    if schema.get('spliceai_storage')=='ALLELE_INFO_V1_ROW_PROJECTION_ON_READ':
        need(schema.get('projection_fields')==SPLICE_FIELDS and not set(SPLICE_FIELDS)&set(fields),'Exact approved allele INFO projection required')
        return fields+SPLICE_FIELDS
    need(fields[-len(SPLICE_FIELDS):]==SPLICE_FIELDS,'Exact approved SpliceAI extension required')
    return fields


class HashReader:
    def __init__(self,stream):self.stream=stream;self.digest=hashlib.sha256()
    def read(self,size=-1):
        data=self.stream.read(size);self.digest.update(data);return data
    def seekable(self):return False


def source_jobs(path,schema,prefix,headers_out,folder,expected_sha):
    """Read/hash compressed source once; batch boundaries come from authority."""
    with Path(path).open('rb') as raw:
        hashed=HashReader(raw)
        with gzip.GzipFile(fileobj=hashed,mode='rb') as compressed,io.TextIOWrapper(compressed,encoding='utf-8',newline='') as stream:
            first=None;seen=set()
            for line in stream:
                if line.startswith('#'):headers_out.append(line)
                else:first=line;break
            need(first is not None,'No VCF records')
            for number,(final,base) in enumerate(zip(schema['batches'],prefix['batches'])):
                lines=[];batch_bytes=0
                for _ in range(final['record_count']):
                    line=first if first is not None else next(stream,None);first=None
                    need(line is not None and not line.startswith('#'),'VCF record count/header differs')
                    identity=line.split('\t',5);key=(identity[0],int(identity[1]),identity[3],identity[4])
                    need(key not in seen,'Duplicate exact normalized source allele; unique-allele views would be ambiguous');seen.add(key)
                    batch_bytes+=len(line.encode('utf-8'))
                    need(batch_bytes<=64*1024**2,'One authority batch exceeds the reviewed 64 MiB input memory ceiling')
                    lines.append(line)
                yield dict(folder=str(folder/str(number).zfill(5)),proof=final,prefix_proof=base,prefix_count=len(prefix['full_fields']),lines=lines,input_bytes=batch_bytes)
            need(next(stream,None) is None,'Unexpected extra VCF record')
        need(hashed.digest.hexdigest()==expected_sha,'Source VCF SHA differs; incomplete staging is not promoted')


def concatenate_tables(results,output):
    fields={}
    for result in results:
        for name,columns in result['table_fields'].items():
            need(name not in fields or fields[name]==columns,'Worker table field order differs');fields[name]=columns
    for name,columns in sorted(fields.items()):
        target=output/name;target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('xb') as out:
            header=io.StringIO(newline='');csv.writer(header,delimiter='\t',lineterminator='\n').writerow(columns)
            out.write(gzip.compress(header.getvalue().encode(),compresslevel=3,mtime=0))
            for result in results:
                part=Path(result['folder'])/name
                if part.exists():
                    with part.open('rb') as inp:shutil.copyfileobj(inp,out,8*1024**2)
    return fields


def write_table(path,fields,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with (gzip.open(path,'xt',encoding='utf-8',newline='',compresslevel=3) if str(path).endswith('.gz') else path.open('x',encoding='utf-8',newline='')) as stream:
        writer=csv.DictWriter(stream,fields,delimiter='\t',extrasaction='ignore',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def sample_vcfs(rule,schema,headers,samples,results,output,codec,engine):
    import pysam
    labels={s:'sample_'+str(i+1).zfill(4) for i,s in enumerate(samples)};records={};rows={}
    for sample in samples:
        count=sum(r['sample_records'].get(sample,0) for r in results);row_count=sum(r['sample_rows'].get(sample,0) for r in results)
        logical_fields=schema['full_fields']+schema.get('projection_fields',[])
        own=dict(schema,scope='sample',sample=sample,row_count=row_count,ordinal_in_vcf_order=False,full_fields=logical_fields,
            source_ordered_batch_merkle_root=schema['ordered_batch_merkle_root'],sample_packet_contract='ONLY_SELECTED_SAMPLE_ROWS_SOURCE_GLOBAL_ORDINALS')
        own.pop('batches',None);own.pop('ordered_batch_merkle_root',None)
        own['appended_fields']=schema.get('appended_fields',[])+schema.get('projection_fields',[])
        own['source_projection_contract']=schema.get('spliceai_storage');own.pop('projection_fields',None);own.pop('spliceai_storage',None)
        subset=[line for line in headers if not line.startswith('##WGS_EVIDENCE_SCHEMA=')]
        insertion=next(i for i,line in enumerate(subset) if line.startswith('#CHROM'))
        subset.insert(insertion,codec.schema_header(own).decode())
        target=output/'vcf'/(labels[sample]+'.annotated.vcf.gz');target.parent.mkdir(exist_ok=True)
        temporary=target.with_suffix('.header.bgz')
        with pysam.BGZFile(str(temporary),'wb') as stream:stream.write(engine.subset_header(subset,sample))
        parts=[(temporary,sha(temporary)),*[(Path(r['folder'])/'vcf'/(labels[sample]+'.body.bgz'),r['persisted_sample_parts'][sample]['sha256']) for r in results if r['sample_records'].get(sample,0)]]
        with target.open('xb') as out:
            for part,expected_sha in parts:
                copied=hashlib.sha256()
                with part.open('rb') as inp:
                    need(part.stat().st_size>=len(BGZF_EOF),'Truncated BGZF part');remaining=part.stat().st_size-len(BGZF_EOF)
                    while remaining:
                        block=inp.read(min(8*1024**2,remaining));need(block,'Truncated BGZF part');out.write(block);copied.update(block);remaining-=len(block)
                    need(inp.read()==BGZF_EOF,'Unexpected BGZF end marker')
                    copied.update(BGZF_EOF);need(copied.hexdigest()==expected_sha,'A persisted scope-verified sample BGZF block changed during assembly')
            out.write(BGZF_EOF)
        temporary.unlink();pysam.tabix_index(str(target),preset='vcf')
        seen=0
        with pysam.VariantFile(str(target)) as stream:
            need(list(stream.header.samples)==[sample],'Sample VCF header scope differs')
            for record in stream.fetch():
                gt=record.samples[sample]['GT'];need(1 in gt,'Sample VCF contains an unobserved ALT')
                ac=sum(t==1 for t in gt);an=sum(t is not None for t in gt)
                need(tuple(record.info['AC'])==(ac,) and record.info['AN']==an and abs(record.info['AF'][0]-ac/an)<1e-6,'Derived subset AC/AN/AF differs')
                # Decode persisted first/last independently; each emitted packet was
                # scope-checked in the worker, without a second whole-source pass.
                if seen==0 or seen==count-1:
                    packet=codec.decode_packet(record.info['WGS_EVIDENCE'])
                    need(packet.get('sample')==sample and all(r['values'][logical_fields.index('sample')]==sample for r in packet['rows']),'Persisted packet privacy audit failed')
                seen+=1
        need(seen==count,'Indexed sample VCF count differs');records[sample]=count;rows[sample]=row_count
    return records,rows


def export(rule):
    """Reviewed execution entry point. No subprocess/process actuation of primary."""
    guard(rule['deadline_utc']);need(1<=rule['workers']<=8,'At most eight bounded workers')
    need(rule.get('root_reviewed') is True,'Exact ROOT-reviewed rule required')
    need(sha(__file__)==rule['exporter_sha256'],'Exporter code pin differs')
    output=Path(rule['output']);need(not output.exists(),'Refusing existing canonical output')
    manifest_path=verify(rule['fully_manifest']);manifest=load(manifest_path)
    need(manifest.get('status')=='PASS_FULLY_ANNOTATED_VCF_DIRECT' and manifest.get('schema_version')==1,'Final direct VCF authority is not PASS')
    need(all(manifest.get(flag) is True for flag in ('original_records_fields_preserved','packet_decode_verified','direct_core_prefix_verified')),'Final producer proof flags missing')
    core_path=verify(manifest['core_manifest']);core=load(core_path)
    need(core.get('status')=='DIRECT_PACKET_BASELINE_PASS' and core.get('no_legacy_receipt_created') is True,'Authentic new direct baseline authority required')
    # Deterministic conditions are pinned in the ROOT-reviewed rule. They bind
    # the completed producer to the approved source/control/code identities.
    expected=rule.get('authority_expected',{})
    need(rule.get('synthetic') is True or bool(expected),'Exact reviewed producer/source authority expectations required')
    for dotted,value in expected.items():
        need(dotted.startswith(('core.','final.')),'Unknown authority expectation namespace')
        target=core if dotted.startswith('core.') else manifest
        for key in dotted.split('.')[1:]:
            need(isinstance(target,dict) and key in target,'Expected source authority field missing');target=target[key]
        need(target==value,'Completed producer source/control/code identity differs from reviewed rule')
    proof=module(rule['proof_module'],'parent_direct_proof');codec=module(rule['codec_module'],'parent_final_codec');engine=module(rule['view_policy'],'parent_frozen_policy')
    rank=module(rule['ranking_module'],'parent_frozen_rank');need(core['schema']['full_fields']==rank.FIELDS,'Direct original fields differ from sealed annotation_row')
    stored_schema=manifest['schema'];schema=dict(stored_schema);prefix=core['schema']
    if manifest.get('projection'):
        projection=manifest['projection']
        need(projection.get('projection_contract')=='ALLELE_INFO_V1_ROW_PROJECTION_ON_READ' and projection.get('projection_fields')==SPLICE_FIELDS,'Unknown same-VCF projection contract')
        need(manifest.get('core33_packets_byte_identical') is True,'Stored core33 packet preservation proof required')
        schema.update(spliceai_storage=projection['projection_contract'],projection_fields=projection['projection_fields'])
    fields=validate_schema(schema,prefix,proof)
    need(schema.get('appended_fields',[])+schema.get('projection_fields',[])==rule['approved_appended_fields'],'Unapproved appended/projected packet fields')
    source=verify(manifest['enriched_vcf'],content=False);verify(manifest['index'])
    samples=manifest['samples'];need(len(samples)==len(set(samples))==4,'Exactly four private source sample IDs required')
    if hasattr(os,'sched_getaffinity'):
        os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:rule['workers']]))
    staging=output.parent/('.'+output.name+'-'+uuid.uuid4().hex);staging.mkdir(mode=0o700)
    work=staging/'.batch-work';work.mkdir();headers=[];results=[]
    accumulated_tops=defaultdict(list);accumulated_splice=defaultdict(list);accumulated_counts=defaultdict(Counter)
    accumulated_terms=defaultdict(lambda:defaultdict(Counter));gene_counts=defaultdict(Counter)
    def consume(result):
        for scope,values in result.pop('tops').items():
            for value in values:keep_top(accumulated_tops[scope],value,lambda item:item[0])
        for scope,values in result.pop('splice_tops').items():
            for value in values:keep_top(accumulated_splice[scope],value,lambda item:item[0])
        for scope,value in result.pop('counts').items():accumulated_counts[scope].update(value)
        for scope,terms in result.pop('terms_counts').items():
            for term,value in terms.items():accumulated_terms[scope][term].update(value)
        for gene,sample,kind,n in result.pop('genes'):gene_counts[gene][(sample,kind)]+=n
        results.append(result)
    try:
        jobs=source_jobs(source,schema,prefix,headers,work,manifest['enriched_vcf']['sha256'])
        with ProcessPoolExecutor(max_workers=rule['workers'],initializer=init_worker,initargs=(rule,fields,samples,schema['full_fields'])) as pool:
            pending=[];queued_bytes=0
            for job in jobs:
                guard(rule['deadline_utc'])
                while pending and (len(pending)>=rule['workers'] or queued_bytes+job['input_bytes']>128*1024**2):
                    future,n=pending.pop(0);consume(future.result());queued_bytes-=n
                pending.append((pool.submit(process_batch,job),job['input_bytes']));queued_bytes+=job['input_bytes']
            for future,n in pending:consume(future.result())
        need(proof.merkle_digest(schema['full_fields'],[r['proof'] for r in results])==schema['ordered_batch_merkle_root'],'Decoded stored-vector Merkle proof differs')
        need(proof.merkle_digest(prefix['full_fields'],[r['prefix_proof'] for r in results])==prefix['ordered_batch_merkle_root'],'Decoded frozen prefix Merkle proof differs')
        header_schema=[line for line in headers if line.startswith('##WGS_EVIDENCE_SCHEMA=')]
        need(len(header_schema)==1,'Exactly one packet schema required')
        token=re.search('Data="([A-Za-z0-9_-]+)"',header_schema[0]);need(token and codec.decode_schema(token[1])==stored_schema,'Stored header schema differs from authority')
        if manifest.get('projection'):
            declared=[line for line in headers if line.startswith('##WGS_SPLICEAI_SCHEMA=')]
            need(len(declared)==1 and codec.decode_projection_schema(declared[0])==manifest['projection'],'Same-VCF projection header differs from final authority')
        need(next(line for line in headers if line.startswith('#CHROM')).rstrip('\r\n').split('\t')[9:]==samples,'VCF private sample identity differs')
        codec.decode_header_bundle(headers,'WGS_PUBLIC_REFERENCE');codec.decode_header_bundle(headers,'WGS_GENE_CONTEXT')
        columns=concatenate_tables(results,staging)
        for sample in samples:
            target=staging/'masters'/('sample_'+str(samples.index(sample)+1).zfill(4)+'.master.tsv.gz')
            if not target.exists():write_table(target,fields,[])
        counts,rows=sample_vcfs(rule,schema,headers,samples,results,staging,codec,engine)
        labels={s:'sample_'+str(i+1).zfill(4) for i,s in enumerate(samples)};scopes=['cohort','shared_ge2','all4',*samples]
        final_counts={};category_rows=[];acmg=[]
        for scope in scopes:
            top=accumulated_tops[scope]
            label=labels.get(scope,scope);fields_out=['rank',*list(dict.fromkeys([*engine.SUMMARY,*fields,'scope_raw_member','scope_qc_member','scope_candidate_member']))]
            if not (staging/'tables'/(label+'.functional.tsv.gz')).exists():write_table(staging/'tables'/(label+'.functional.tsv.gz'),fields_out[1:],[])
            write_table(staging/'tables'/(label+'.top20.tsv'),fields_out,[dict(v[1],rank=i+1) for i,v in enumerate(top)])
            count=accumulated_counts[scope];tc=accumulated_terms[scope]
            final_counts[scope]=dict(count,top20_unique_alleles=len(top))
            for term in sorted(set(engine.CORE)|set(tc)):
                category_rows.append(dict(scope=label,category=term,raw_unique_alleles=tc[term]['raw'],qc_unique_alleles=tc[term]['qc'],eligible_unique_alleles=tc[term]['eligible'],path=label+'/'+category_name(term)+'.tsv.gz'))
                target=staging/'tables'/label/(category_name(term)+'.tsv.gz')
                if not target.exists():write_table(target,fields_out[1:],[])
            if scope in ['cohort',*samples]:
                for i,(_,row) in enumerate(top,1):
                    for criterion in rank.ACMG:acmg.append(dict(scope=label,rank=i,allele_id=row['allele_id'],sample=row['sample'],gene=row.get('gene',''),transcript=row.get('transcript',''),criterion=criterion,status='NOT_ASSESSED_MANUAL_REVIEW_REQUIRED',strength='',source='',rationale=''))
                splice=accumulated_splice[scope]
                sf=[*fields_out,'spliceai_max_delta','spliceai_all_gene_predictions_json','spliceai_prediction_ranking_scope']
                write_table(staging/'tables/spliceai'/(label+'.top20.tsv'),sf,[dict(v[1],rank=i+1) for i,v in enumerate(splice)])
        write_table(staging/'tables/functional-category-index.tsv',['scope','category','raw_unique_alleles','qc_unique_alleles','eligible_unique_alleles','path'],category_rows)
        write_table(staging/'tables/ACMG-manual-top20-worksheet.tsv.gz',['scope','rank','allele_id','sample','gene','transcript','criterion','status','strength','source','rationale'],acmg)
        write_table(staging/'tables/gene-overlap-SEPARATE-NOT-SHARED-ALLELES.tsv.gz',['gene','raw_sample_allele_counts_json','qc_sample_allele_counts_json','scope'],[
            dict(gene=g,raw_sample_allele_counts_json=js({s:c[(s,'raw')] for s in samples}),qc_sample_allele_counts_json=js({s:c[(s,'qc')] for s in samples}),scope='GENE_OVERLAP_DIFFERENT_ALLELES_ALLOWED_NOT_ALLELE_INTERSECTION') for g,c in sorted(gene_counts.items())])
        # Canonical single inbound VCF copy: hardlink on the VM; one PC backup.
        (staging/'vcf').mkdir(exist_ok=True)
        for item,name in [(manifest['enriched_vcf'],'cohort.fully-annotated.vcf.gz'),(manifest['index'],'cohort.fully-annotated.vcf.gz.tbi')]:os.link(item['path'],staging/'vcf'/name)
        provenance=staging/'provenance';provenance.mkdir()
        shutil.copyfile(manifest_path,provenance/'fully-annotated-vcf-manifest.json');shutil.copyfile(core_path,provenance/'direct-core-manifest.json')
        # Only explicitly receipt-bound applied outputs; never raw reference DBs.
        for item in manifest.get('applied_artifacts',[]):
            relative=Path(item['relative_path'])
            need(not relative.is_absolute() and '..' not in relative.parts and relative.parts[0]=='applied-evidence','Applied output path outside allowlist')
            origin=verify(item);target=staging/relative;target.parent.mkdir(parents=True,exist_ok=True);os.link(origin,target)
        shutil.rmtree(work)
        page=['<!doctype html><meta charset="utf-8"><title>One annotation: cohort and samples</title><h1>One annotation: cohort and four sample views</h1>',
            '<p>Exact-allele research prioritization. No automatic ACMG class, NMD, LoF/GoF or clinical burden assignment. Missing and filtered calls remain explicit.</p>']
        for scope in ['cohort',*samples]:
            label=labels.get(scope,scope);page.append('<h2>'+html.escape(scope)+'</h2><p>'+ ' · '.join('<a href="'+p+'">'+n+'</a>' for p,n in [
                ('masters/'+label+'.master.tsv.gz','Complete transcript master'),('tables/'+label+'.top20.tsv','Top20'),('tables/'+label+'.functional.tsv.gz','Functional consequences'),('tables/spliceai/'+label+'.top20.tsv','SpliceAI prediction top20')])+'</p>')
            if scope in samples:page.append('<p><a href="vcf/'+label+'.annotated.vcf.gz">Individual annotated VCF</a> · <a href="vcf/'+label+'.annotated.vcf.gz.tbi">TBI</a></p>')
        page.append('<p><a href="vcf/cohort.fully-annotated.vcf.gz">Complete self-contained cohort VCF</a> · <a href="tables/functional-category-index.tsv">Functional type files/counts</a> · <a href="tables/allele-membership-matrix.tsv.gz">Exact allele memberships</a> · <a href="tables/ACMG-manual-top20-worksheet.tsv.gz">Manual ACMG worksheet</a></p>')
        (staging/'index.html').write_text('\n'.join(page),encoding='utf-8')
        artifacts={str(p.relative_to(staging)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(staging.rglob('*')) if p.is_file()}
        result=dict(schema_version=1,status='PASS',synthetic=bool(rule.get('synthetic',False)),report_input_kind='FULLY_ANNOTATED_VCF_PACKETS_V1',
            report_producer_kind='FAST_DIRECT_PACKET_VIEWS_V1',samples=samples,sample_file_labels=labels,source_pins=dict(fully_manifest=rule['fully_manifest'],core_manifest=manifest['core_manifest'],vcf=manifest['enriched_vcf']),
            ranking_module=dict(rule['ranking_module'],function='rank_key'),batch_merkle_root=schema['ordered_batch_merkle_root'],core_prefix_merkle_root=prefix['ordered_batch_merkle_root'],
            producer_code_sha256=sha(__file__),code_pins={k:rule[k] for k in ('proof_module','codec_module','view_policy','ranking_module','projection_module') if k in rule},
            logical_full_fields=fields,derived_batch_merkle_root=proof.merkle_digest(fields,[r['derived_proof'] for r in results]),
            master_transcript_rows=schema['row_count'],vcf_records=sum(r['proof']['record_count'] for r in results),sample_vcf_records=counts,sample_master_rows=rows,
            counts=final_counts,sample_packet_scope_verified=True,public_header_bundles_verified=True,top20_padding=False,automatic_ACMG=False,
            sample_packet_checks=[dict(sample=s,scope='sample',other_sample_rows=0,original_source_ordinals_retained=True,row_count=rows[s],
                persisted_payload_verified_by='all emitted bytes SHA256 == BGZF reread bytes; every copied compressed block SHA256 verified; indexed GT/AC/AN/AF audit') for s in samples],
            no_external_requests=True,primary_workflow_changed=False,no_legacy_receipt_created=True,artifacts=artifacts,completed_utc=now().isoformat())
        (staging/'manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        # Linux renameat2 NOREPLACE closes the final existence-check race.
        import ctypes
        library=ctypes.CDLL(None,use_errno=True)
        need(hasattr(library,'renameat2'),'Atomic no-overwrite promotion requires renameat2')
        if library.renameat2(-100,os.fsencode(staging),-100,os.fsencode(output),1):
            error=ctypes.get_errno();raise OSError(error,'Atomic no-overwrite report promotion refused')
        return result
    except BaseException:
        (staging/'INCOMPLETE.txt').write_text('No completed manifest was promoted. Original inputs and active workflows remain unchanged.\n');raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--rule',type=Path,required=True);parser.add_argument('--execute',action='store_true')
    args=parser.parse_args();rule=load(args.rule)
    if args.execute:print(js(export(rule)))
    else:print(js(dict(status='PREPARE_ONLY_NO_PATIENT_EXECUTION',workers=rule.get('workers'),deadline_utc=rule.get('deadline_utc'),output=rule.get('output'))))
