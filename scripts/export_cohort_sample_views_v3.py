#!/usr/bin/env python3
"""VCF packets are the report input; joined TSV is only an exact lossless control."""
import csv
import gzip
import hashlib
import html
import importlib.util
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import uuid
import zlib

INFO_ID='WGS_EVIDENCE'
CODEC='json-utf8+zlib+base64url-nopad'
SPLICE_FIELDS=['spliceai_status','spliceai_gene','spliceai_DS_AG','spliceai_DS_AL','spliceai_DS_DG','spliceai_DS_DL',
    'spliceai_DP_AG','spliceai_DP_AL','spliceai_DP_DG','spliceai_DP_DL','spliceai_prediction_json','spliceai_source_json','spliceai_transcript_scope']


def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def module(path,expected,name):
    path=Path(path)
    if path.is_symlink() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('Pinned views/codec module digest differs')
    spec=importlib.util.spec_from_file_location(name,path);value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value


class DigestText:
    def __init__(self,fields):
        self.digest=hashlib.sha256();self.writer=csv.writer(self,delimiter='\t',lineterminator='\n');self.writer.writerow(fields)
    def write(self,value):self.digest.update(value.encode('utf-8'));return len(value)
    def row(self,values):self.writer.writerow(values)


def header_schema(line,codec):
    match=re.fullmatch(r'##WGS_EVIDENCE_SCHEMA=<Version=1,Codec="'+re.escape(CODEC)+r'",Data="([A-Za-z0-9_-]+)">\r?\n?',line)
    if not match:raise ValueError('Unknown VCF packet schema declaration')
    return codec.decode_schema(match[1])


def schema_line(schema,codec):
    return '##WGS_EVIDENCE_SCHEMA=<Version=1,Codec="'+CODEC+'",Data="'+codec.encode_schema(schema)+'">\n'


def info_dict(text):
    result={}
    for part in text.split(';'):
        if not part or part=='.':continue
        key,sep,value=part.partition('=')
        if key in result:raise ValueError('Duplicate VCF INFO key')
        result[key]=value if sep else None
    return result


def info_text(info):return ';'.join(k if v is None else k+'='+v for k,v in info.items()) or '.'


def validate_fully(rule,engine,helper,guard):
    workflow=Path(rule['workflow_root']).resolve();path=Path(rule['fully_vcf_manifest']).resolve();folder=path.parent
    if path!=workflow.parent/'fully-annotated-vcf-v1/manifest.json':raise ValueError('Exact fully annotated VCF root required')
    value=load(path);receipt_path=folder/'workflow-receipt.json';receipt=load(receipt_path)
    if value.get('schema_version')!=1 or value.get('status')!=rule['fully_manifest_status']:raise ValueError('Completed fully annotated VCF manifest required')
    if any(value.get(flag) is not True for flag in ('original_records_fields_preserved','packet_decode_verified','decoded_master_vectors_verified')):
        raise ValueError('Completed producer lossless/original-record verification required')
    if receipt.get('schema_version')!=1 or receipt.get('status')!='PASS' or receipt.get('stage')!='fully-annotated-vcf-v1' or not re.fullmatch('[0-9a-f]{64}',receipt.get('identity','')):raise ValueError('Fully annotated VCF completion receipt required')
    sources,primary_pins=engine.checked_inputs(workflow,guard)
    joined,control_pins=helper.verify_joined(rule['control_joined_manifest'],workflow,guard)
    for item,expected in [(value['base_vcf'],sources['vcf']),(value['master_manifest'],Path(rule['control_joined_manifest'])),(value['master'],Path(joined['enriched_master']['path'])),
                         (value['combined_receipt'],Path(rule['control_joined_manifest']).parent/'workflow-receipt.json')]:
        if item.get('path')!=str(expected):raise ValueError('Fully annotated VCF lineage path differs')
        helper.check_pin(item,guard)
    if value['master']!=joined['enriched_master']:raise ValueError('Fully annotated VCF control master identity differs')
    control_fields=joined['original_fields']+joined['appended_fields'];extra=rule['packet_append_fields']
    if extra!=SPLICE_FIELDS:raise ValueError('Only the exact reviewed SpliceAI field extension is allowed')
    expected_fields=control_fields+extra
    schema=value['schema']
    if schema.get('v')!=1 or schema.get('codec')!=CODEC or schema.get('scope')!='cohort' or schema.get('sample') is not None:raise ValueError('Cohort packet scope/codec differs')
    if schema.get('original_fields')!=joined['original_fields'] or schema.get('appended_fields')!=joined['appended_fields']+extra or schema.get('full_fields')!=expected_fields:
        raise ValueError('VCF packet complete field schema differs from joined control')
    if value.get('control_prefix_fields')!=control_fields or value.get('supplement_fields')!=extra or value.get('master_role')!='ORIGINAL_PREFIX_CONTROL':raise ValueError('VCF original-prefix/supplement scope differs')
    if schema.get('row_count')!=joined['row_count'] or schema.get('invariant_digest')!=joined['invariant_digest']:raise ValueError('VCF original-row invariant/count differs')
    for field in ('original_fields','appended_fields','full_fields','row_count','invariant_digest','full_fields_digest'):
        if value.get(field)!=schema.get(field):raise ValueError('Producer manifest/schema field differs: '+field)
    pins=list(control_pins)+[helper.pin(path,guard),helper.pin(receipt_path,guard)]
    splice_root=workflow.parent/'spliceai-reference-context-v1'
    if Path(rule['spliceai_manifest']).resolve()!=splice_root/'manifest.json':raise ValueError('Exact local SpliceAI result root required')
    for key,name in [('spliceai_manifest','manifest.json'),('spliceai_receipt','workflow-receipt.json'),('spliceai_predictions','predictions.tsv.gz')]:
        item=value[key]
        if item.get('path')!=str(splice_root/name):raise ValueError('Unknown SpliceAI source path')
        helper.check_pin(item,guard);pins.append(item)
    splice_manifest=load(splice_root/'manifest.json');splice_receipt=load(splice_root/'workflow-receipt.json')
    if splice_manifest.get('status')!=rule['spliceai_manifest_status'] or splice_receipt.get('status')!='PASS' or splice_receipt.get('stage')!='spliceai-reference-context-v1':raise ValueError('Completed reviewed local SpliceAI scores required')
    for name,key in [('manifest.json','spliceai_manifest'),('predictions.tsv.gz','spliceai_predictions')]:
        if any(splice_receipt.get('artifacts',{}).get(name,{}).get(k)!=value[key][k] for k in ('bytes','sha256')):raise ValueError('SpliceAI completion artifact differs')
    for key,name in [('enriched_vcf','cohort.fully-annotated.vcf.gz'),('index','cohort.fully-annotated.vcf.gz.tbi')]:
        item=value[key]
        if item.get('path')!=str(folder/name):raise ValueError('Unknown fully annotated VCF artifact path')
        helper.check_pin(item,guard);pins.append(item)
    for name,item in [('manifest.json',helper.pin(path,guard)),('cohort.fully-annotated.vcf.gz',value['enriched_vcf']),('cohort.fully-annotated.vcf.gz.tbi',value['index'])]:
        expected=receipt.get('artifacts',{}).get(name,{})
        if any(expected.get(k)!=item[k] for k in ('bytes','sha256')):raise ValueError('Fully annotated completion artifact digest differs')
    primary_pins.update({p['path']:{k:p[k] for k in ('bytes','sha256')} for p in pins})
    return value,joined,pins,primary_pins


def extract_master(value,codec,target,control,guard):
    fields=value['schema']['full_fields'];original=value['schema']['original_fields'];sample_col=fields.index('sample')
    positions={name:fields.index(name) for name in ('chrom','pos','ref','alt')};control_fields=value['control_prefix_fields']
    full_digest=DigestText(fields);original_digest=DigestText(original);row_count=0;record_count=0;samples=None
    sample_metrics={};namespace=set();schema_seen=None;public_headers=[]
    csv.field_size_limit(64*1024**2)
    db_path=target.with_suffix('.ordinal-index.sqlite');db=sqlite3.connect(db_path)
    db.executescript('PRAGMA cache_size=-32768;PRAGMA mmap_size=0;PRAGMA temp_store=FILE;CREATE TABLE rows(i INTEGER PRIMARY KEY,payload BLOB,sample TEXT,observed INTEGER);')
    with gzip.open(value['enriched_vcf']['path'],'rt',encoding='utf-8',newline='') as vcf,gzip.open(target,'wt',encoding='utf-8',newline='') as out,gzip.open(control,'rt',encoding='utf-8',newline='') as comparison:
        expected=csv.reader(comparison,delimiter='\t')
        if next(expected)!=control_fields or fields!=control_fields+value['supplement_fields']:raise ValueError('Control prefix/approved full header differs')
        writer=csv.writer(out,delimiter='\t',lineterminator='\n');writer.writerow(fields)
        for line in vcf:
            if line.startswith('##WGS_EVIDENCE_SCHEMA='):
                if schema_seen is not None:raise ValueError('Duplicate VCF packet schema header')
                schema_seen=header_schema(line,codec)
            elif line.startswith('##INFO=<ID='):
                tag=line.split('ID=',1)[1].split(',',1)[0]
                if tag in codec.INFO_IDS:namespace.add(tag)
            elif line.startswith(('##WGS_PUBLIC_REFERENCE=','##WGS_GENE_CONTEXT=')):public_headers.append(line)
            elif line.startswith('#CHROM'):
                samples=line.rstrip('\r\n').split('\t')[9:]
                if len(samples)!=4 or len(set(samples))!=4:raise ValueError('Four exact cohort VCF sample IDs required')
                sample_metrics={s:dict(rows=0,full=DigestText(fields),original=DigestText(original)) for s in samples}
                for name in ('WGS_PUBLIC_REFERENCE','WGS_GENE_CONTEXT'):codec.decode_header_bundle(public_headers,name)
            elif line.startswith('#'):continue
            else:
                if samples is None or schema_seen!=value['schema'] or namespace!=set(codec.INFO_IDS):raise ValueError('VCF header/schema differs from completed manifest')
                record_count+=1
                if record_count%1000==1:guard.check()
                columns=line.rstrip('\r\n').split('\t')
                if len(columns)!=13:raise ValueError('Cohort VCF column count differs')
                allele=[columns[0],int(columns[1]),columns[3],columns[4]]
                info=info_dict(columns[7]);packet=codec.decode_packet(info.get(INFO_ID,''))
                if packet.get('v')!=1 or packet.get('allele')!=allele:raise ValueError('Embedded packet allele/version differs')
                formats=columns[8].split(':');gt_index=formats.index('GT')
                observed={s for s,call in zip(samples,columns[9:]) if '1' in re.split(r'[/|]',call.split(':')[gt_index])}
                seen_samples=set()
                for item in packet['rows']:
                    row_count+=1;values=item['values']
                    ordinal=item.get('i')
                    if type(ordinal) is not int or not 1<=ordinal<=value['schema']['row_count'] or len(values)!=len(fields) or any(not isinstance(v,str) for v in values):raise ValueError('Packet global ordinal/field vector differs')
                    if [values[positions['chrom']],int(values[positions['pos']]),values[positions['ref']],values[positions['alt']]]!=allele:raise ValueError('Packet row genomic allele differs')
                    sample=values[sample_col]
                    if sample not in sample_metrics:raise ValueError('Packet row has unknown sample')
                    seen_samples.add(sample)
                    if not values[fields.index('spliceai_status')]:raise ValueError('SpliceAI missingness/status cannot be blank')
                    scope=values[fields.index('spliceai_transcript_scope')]
                    if scope not in ('','GENCODE24_GENE_UNION_NOT_ENST49_MAPPING'):raise ValueError('Unreviewed SpliceAI transcript-model mapping')
                    db.execute('INSERT INTO rows VALUES(?,?,?,?)',(ordinal,zlib.compress(json.dumps(values,ensure_ascii=False).encode(),3),sample,int(sample in observed)))
                if seen_samples!=set(samples):raise ValueError('Current sealed report-engine compatibility precondition failed: each source allele needs explicit rows for all four samples; no source records were dropped')
                if record_count%1000==0:
                    db.commit()
                    if shutil.disk_usage(target.parent).free<64*1024**3 or db_path.stat().st_size>64*1024**3:raise ValueError('Bounded ordinal-index/disk reserve reached')
        db.commit()
        for ordinal,payload,sample,observed in db.execute('SELECT i,payload,sample,observed FROM rows ORDER BY i'):
            if ordinal%10000==1:guard.check()
            values=json.loads(zlib.decompress(payload));before=next(expected,None)
            if before!=values[:len(control_fields)]:raise ValueError('VCF packet changed/lost/reordered a complete original/control row')
            writer.writerow(values);full_digest.row(values);original_digest.row(values[:len(original)])
            if observed:
                metric=sample_metrics[sample];metric['rows']+=1;metric['full'].row(values);metric['original'].row(values[:len(original)])
        if next(expected,None) is not None:raise ValueError('VCF packet stream omitted control rows')
    db.close();db_path.unlink()
    if row_count!=value['schema']['row_count'] or original_digest.digest.hexdigest()!=value['schema']['invariant_digest'] or full_digest.digest.hexdigest()!=value['schema']['full_fields_digest']:
        raise ValueError('Extracted master full/original digest or row count differs')
    metrics={s:dict(rows=m['rows'],full_fields_digest=m['full'].digest.hexdigest(),invariant_digest=m['original'].digest.hexdigest()) for s,m in sample_metrics.items()}
    return dict(rows=row_count,records=record_count,full_fields_digest=full_digest.digest.hexdigest(),invariant_digest=original_digest.digest.hexdigest(),samples=samples,
                sample_metrics=metrics,namespaced_info_ids=sorted(namespace),public_header_bundles_verified=True,
                public_header_digest=hashlib.sha256(''.join(public_headers).encode()).hexdigest())


def prune_sample(path,sample,source_schema,metrics,namespace,codec,guard):
    import pysam
    fields=source_schema['full_fields'];sample_col=fields.index('sample');tmp=path.with_name(path.name+'.sample-scoped.tmp.gz')
    schema=dict(source_schema,scope='sample',sample=sample,row_count=metrics['rows'],full_fields_digest=metrics['full_fields_digest'],invariant_digest=metrics['invariant_digest'])
    actual_full=DigestText(fields);actual_original=DigestText(source_schema['original_fields']);count=0;records=0
    db_path=path.with_name(path.name+'.scope-check.sqlite');db=sqlite3.connect(db_path)
    db.executescript('PRAGMA cache_size=-32768;PRAGMA mmap_size=0;PRAGMA temp_store=FILE;CREATE TABLE original(i INTEGER PRIMARY KEY,payload BLOB);CREATE TABLE verified(i INTEGER PRIMARY KEY,payload BLOB);')
    with gzip.open(path,'rt',encoding='utf-8',newline='') as inp,pysam.BGZFile(str(tmp),'wb') as out:
        for line in inp:
            if line.startswith('##WGS_EVIDENCE_SCHEMA='):out.write(schema_line(schema,codec).encode());continue
            if line.startswith('#'):
                if line.startswith('#CHROM') and line.rstrip('\r\n').split('\t')[9:]!=[sample]:raise ValueError('Sample VCF header differs')
                out.write(line.encode());continue
            records+=1
            if records%1000==1:guard.check()
            columns=line.rstrip('\r\n').split('\t');info=info_dict(columns[7]);packet=codec.decode_packet(info[INFO_ID])
            selected=codec.select_sample(packet,fields,sample)
            expected=[item for item in packet['rows'] if item['values'][sample_col]==sample]
            if selected.get('rows')!=expected or not expected or selected.get('scope')!='sample' or selected.get('sample')!=sample:raise ValueError('Sample packet did not retain exact own transcript rows/ordinals/scope')
            for item in selected['rows']:
                if item['values'][sample_col]!=sample:raise ValueError('Another individual remains in sample packet')
                count+=1;db.execute('INSERT INTO original VALUES(?,?)',(item['i'],zlib.compress(json.dumps(item['values'],ensure_ascii=False).encode(),3)))
            replacement=codec.rebuild_readable(selected,fields)
            if set(replacement)-set(namespace):raise ValueError('Codec emitted an unknown INFO tag')
            updated={k:v for k,v in info.items() if k not in namespace}
            updated.update(replacement);updated[INFO_ID]=codec.encode_packet(selected)
            columns[7]=info_text(updated);out.write(('\t'.join(columns)+'\n').encode())
    db.commit()
    for i,payload in db.execute('SELECT i,payload FROM original ORDER BY i'):
        values=json.loads(zlib.decompress(payload));actual_full.row(values);actual_original.row(values[:len(source_schema['original_fields'])])
    if count!=metrics['rows'] or actual_full.digest.hexdigest()!=metrics['full_fields_digest'] or actual_original.digest.hexdigest()!=metrics['invariant_digest']:
        raise ValueError('Sample-only packet row/digest verification differs')
    # Replace only the new unpublished staging file; never alter a primary output.
    Path(str(path)+'.tbi').unlink();tmp.replace(path);pysam.tabix_index(str(path),preset='vcf',force=False)
    with pysam.VariantFile(str(path)) as check:
        if list(check.header.samples)!=[sample] or sum(1 for _ in check.fetch())!=records:raise ValueError('Sample packet BGZF/TBI verification failed')
    verified_rows=0;verified_schema=None;verified_full=DigestText(fields);public_headers=[]
    with gzip.open(path,'rt',encoding='utf-8',newline='') as check:
        for line in check:
            if line.startswith('##WGS_EVIDENCE_SCHEMA='):verified_schema=header_schema(line,codec)
            elif line.startswith(('##WGS_PUBLIC_REFERENCE=','##WGS_GENE_CONTEXT=')):public_headers.append(line)
            elif not line.startswith('#'):
                tags=info_dict(line.rstrip('\r\n').split('\t')[7]);packet=codec.decode_packet(tags[INFO_ID])
                if packet.get('scope')!='sample' or packet.get('sample')!=sample:raise ValueError('Persisted sample packet scope differs')
                for row in packet['rows']:
                    if row['values'][sample_col]!=sample:raise ValueError('Another individual persisted in sample packet')
                    verified_rows+=1;db.execute('INSERT INTO verified VALUES(?,?)',(row['i'],zlib.compress(json.dumps(row['values'],ensure_ascii=False).encode(),3)))
                for key,value in codec.rebuild_readable(packet,fields).items():
                    if tags.get(key)!=value:raise ValueError('Persisted readable sample INFO differs from selected packet')
    db.commit()
    for i,payload in db.execute('SELECT i,payload FROM verified ORDER BY i'):
        values=json.loads(zlib.decompress(payload));verified_full.row(values)
        if db.execute('SELECT payload FROM original WHERE i=?',(i,)).fetchone()!=(payload,):raise ValueError('Persisted source ordinal/field vector differs')
    db.close();db_path.unlink()
    if verified_schema!=schema or verified_rows!=count or verified_full.digest.hexdigest()!=metrics['full_fields_digest']:raise ValueError('Persisted sample header/rows digest differs')
    for name in ('WGS_PUBLIC_REFERENCE','WGS_GENE_CONTEXT'):codec.decode_header_bundle(public_headers,name)
    return dict(sample=sample,records=records,rows=count,scope='sample',original_source_ordinals_retained=True,other_sample_rows=0,
                full_fields_digest=metrics['full_fields_digest'],invariant_digest=metrics['invariant_digest'],public_header_bundles_verified=True)


def copy_file(expected,target,helper,guard):
    source=Path(expected['path']);helper.check_pin(expected,guard)
    with source.open('rb') as inp,target.open('xb') as out:
        for block in iter(lambda:inp.read(8*1024**2),b''):
            guard.check()
            if shutil.disk_usage(target.parent).free<64*1024**3:raise ValueError('VCF-first disk reserve reached')
            out.write(block)
    return {k:helper.pin(target,guard)[k] for k in ('bytes','sha256')}


def export(rule,deadline=None):
    engine=module(rule['v1_engine'],rule['v1_engine_sha256'],'vcf_first_frozen_engine')
    helper=module(rule['v2_helper'],rule['v2_helper_sha256'],'vcf_first_frozen_control_helper')
    codec=module(rule['codec_module'],rule['codec_sha256'],'vcf_first_reviewed_packet_codec')
    guard=engine.Guard(deadline or engine.ABSOLUTE_DEADLINE);guard.check()
    if hasattr(os,'sched_getaffinity'):os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
    try:
        import resource
        soft,hard=resource.getrlimit(resource.RLIMIT_AS);resource.setrlimit(resource.RLIMIT_AS,(min(4*1024**3,soft) if soft>0 else 4*1024**3,hard))
    except ImportError:pass
    value,joined,pins,source_pins=validate_fully(rule,engine,helper,guard)
    output=Path(rule['output']).resolve();workflow=Path(rule['workflow_root']).resolve()
    if output.exists() or output!=workflow.parent/'cohort-sample-views-v1':raise ValueError('New canonical VCF-first views output required')
    state=Path(rule['state_directory']).resolve()
    if state==workflow or workflow in state.parents or state==output or state in output.parents or output in state.parents:raise ValueError('Separate private extraction state outside primary workflow required')
    work=Path(rule['state_directory'])/('extraction-'+uuid.uuid4().hex);work.mkdir(mode=0o700,parents=True)
    extracted=work/'master.from-vcf.tsv.gz'
    extraction=extract_master(value,codec,extracted,Path(joined['enriched_master']['path']),guard)
    extracted_pin=helper.pin(extracted,guard);source_pins[str(extracted)]={k:extracted_pin[k] for k in ('bytes','sha256')}
    def inputs(root,current_guard):
        current_guard.check()
        return dict(vcf=Path(value['enriched_vcf']['path']),vcf_index=Path(value['index']['path']),master=extracted),dict(source_pins)
    engine.checked_inputs=inputs;engine.PROJECTION=list(engine.PROJECTION)+value['schema']['appended_fields'];engine.SUMMARY=list(engine.SUMMARY)+value['schema']['appended_fields']
    staging=output.parent/('.'+output.name+'.vcf-first-partial-'+uuid.uuid4().hex)
    result=engine.export(workflow,staging,rule['ranking_module'],rule['ranking_sha256'],deadline=guard.deadline)
    applied=helper.copy_applied_evidence(joined,staging,guard);result['artifacts'].update(applied)
    scope=[]
    for sample,label in result['sample_file_labels'].items():
        path=staging/'vcf'/(label+'.annotated.vcf.gz')
        scope.append(prune_sample(path,sample,value['schema'],extraction['sample_metrics'][sample],extraction['namespaced_info_ids'],codec,guard))
        for file in (path,Path(str(path)+'.tbi')):result['artifacts'][str(file.relative_to(staging))]={k:helper.pin(file,guard)[k] for k in ('bytes','sha256')}
    for key,name in [('enriched_vcf','cohort.fully-annotated.vcf.gz'),('index','cohort.fully-annotated.vcf.gz.tbi')]:
        result['artifacts']['vcf/'+name]=copy_file(value[key],staging/'vcf'/name,helper,guard)
    provenance=staging/'provenance';provenance.mkdir()
    for source,name in [(rule['fully_vcf_manifest'],'fully-annotated-vcf-manifest.json'),(str(Path(rule['fully_vcf_manifest']).parent/'workflow-receipt.json'),'fully-annotated-vcf-receipt.json'),
                        (rule['control_joined_manifest'],'joined-control-manifest.json'),(str(Path(rule['control_joined_manifest']).parent/'workflow-receipt.json'),'joined-control-receipt.json')]:
        result['artifacts']['provenance/'+name]=copy_file(helper.pin(source,guard),provenance/name,helper,guard)
    with (staging/'index.html').open('a',encoding='utf-8') as page:page.write('<p><a href="vcf/cohort.fully-annotated.vcf.gz">One fully annotated cohort VCF: source of all master/top20/type views</a> · <a href="vcf/cohort.fully-annotated.vcf.gz.tbi">TBI</a></p>\n')
    result['artifacts']['index.html']={k:helper.pin(staging/'index.html',guard)[k] for k in ('bytes','sha256')}
    result.update(report_input_kind='FULLY_ANNOTATED_VCF_PACKETS_V1',fully_annotated_vcf=dict(manifest=helper.pin(rule['fully_vcf_manifest'],guard),vcf=value['enriched_vcf'],index=value['index'],schema=value['schema']),
        extracted_master=dict(extracted_pin,row_count=extraction['rows'],full_fields_digest=extraction['full_fields_digest'],original_fields_digest=extraction['invariant_digest'],source='VCF_INFO_WGS_EVIDENCE_ONLY'),
        joined_control=dict(manifest=helper.pin(rule['control_joined_manifest'],guard),master=joined['enriched_master'],all_full_field_vectors_order_verified=True,used_as_report_input=False),
        approved_packet_append_fields=value['supplement_fields'],spliceai_source={key:value[key] for key in ('spliceai_manifest','spliceai_receipt','spliceai_predictions')},
        sample_packet_scope_verified=True,sample_packet_checks=scope,public_header_bundles_verified=extraction['public_header_bundles_verified'],public_header_digest=extraction['public_header_digest'],
        views_v1_engine=dict(path=rule['v1_engine'],sha256=rule['v1_engine_sha256']),
        views_v2_control_helper=dict(path=rule['v2_helper'],sha256=rule['v2_helper_sha256']),packet_codec=dict(path=rule['codec_module'],sha256=rule['codec_sha256']),
        views_v3_exporter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),copied_applied_evidence=sorted(k for k in applied if k!='index.html'))
    for item in pins:helper.check_pin(item,guard)
    (staging/'manifest.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n');helper.promote(staging,output)
    return result
