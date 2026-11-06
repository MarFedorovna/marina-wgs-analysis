#!/usr/bin/env python3
"""Separate MT functional views: source observations and literal site-PASS, not nuclear GT rules."""
import argparse
from collections import Counter
import csv
import datetime as dt
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import resource
import signal
import time

END=dt.datetime.fromisoformat('2026-11-06T06:25:00+00:00')
CORE=['stop_gained','stop_lost','start_lost','frameshift_variant','splice_donor_variant',
      'splice_acceptor_variant','splice_region_variant','missense_variant','synonymous_variant',
      'inframe_insertion','inframe_deletion','protein_altering_variant']
SUMMARY=['allele_id','chrom','pos','ref','alt','observed_call_count','observed_call_samples',
         'site_PASS_count','site_PASS_samples','not_observed_samples','sample_source_context_json',
         'consequences_json','genes_json','transcript_effect_rows','frequency_contexts_json',
         'membership_basis','clinical_interpretation']


def helper(path):
    spec=importlib.util.spec_from_file_location('_mt_views_helpers',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result


def observed_allele(alleles,key,sample,record):
    value=alleles.setdefault(key,dict(samples={},terms=set(),genes=set(),rows=0,frequencies={}))
    value['samples'].setdefault(sample,[]).append(record)
    return value


def attach_effect(value,row):
    value['terms'].update(filter(None,row['consequence'].split('&')))
    if row.get('gene'):value['genes'].add(row['gene'])
    value['rows']+=1
    if row.get('frequency_context_json'):
        value['frequencies'].setdefault(row['sample'],set()).add(row['frequency_context_json'])


def summarize(key,value,samples):
    observed=[sample for sample in samples if sample in value['samples']]
    passed=[sample for sample in observed if any(r['FILTER']=='PASS' for r in value['samples'][sample])]
    context={sample:dict(state='OBSERVED_MUTECT2_CALL_RECORD' if sample in observed else 'NOT_OBSERVED_NOT_ASSESSED_NO_REFERENCE_INFERENCE',
                         original_records=value['samples'].get(sample,[])) for sample in samples}
    chrom,pos,ref,alt=key
    return dict(allele_id=':'.join(map(str,key)),chrom=chrom,pos=pos,ref=ref,alt=alt,
        observed_call_count=len(observed),observed_call_samples='&'.join(observed),
        site_PASS_count=len(passed),site_PASS_samples='&'.join(passed),not_observed_samples='&'.join(s for s in samples if s not in observed),
        sample_source_context_json=json.dumps(context,ensure_ascii=False,sort_keys=True),
        consequences_json=json.dumps(sorted(value['terms'] or {'unannotated_or_allele_unresolved'})),
        genes_json=json.dumps(sorted(value['genes'])),transcript_effect_rows=value['rows'],
        frequency_contexts_json=json.dumps({s:sorted(values) for s,values in value['frequencies'].items()},ensure_ascii=False),
        membership_basis='SOURCE_CALL_RECORD_PRESENCE; literal_FILTER_PASS_separate; GT_AF_HL_not_reclassified',
        clinical_interpretation='NOT_ASSIGNED; transcript_effect_is_not_LoF_GoF_PVS1_or_clinical_burden')


def belongs(scope,level,row,samples):
    members=set(filter(None,row['observed_call_samples' if level=='observed' else 'site_PASS_samples'].split('&')))
    if scope=='cohort':return bool(members)
    if scope=='shared_ge2':return len(members)>=2
    if scope=='all4':return len(members)==len(samples)
    return scope in members


def collect(q,codec,vcf,table,sample,expected_digest,alleles,master_writers,guard):
    rows_total=records=unresolved=0;digest=hashlib.sha256()
    with q.text(vcf) as stream,q.text(table) as original:
        headers=q.header(stream)
        q.need(headers[-1].split('\t')[9:]==[sample],'MT packet source contains another sample')
        schema_lines=[line for line in headers if line.startswith('##WGS_MT_CONSEQUENCE_SCHEMA=')]
        q.need(len(schema_lines)==1,'One MT consequence schema required')
        token=re.search(r'Data="([A-Za-z0-9_-]+)"',schema_lines[0]).group(1)
        schema=codec.decode_packet(token);fields=schema['fields'];sample_index=fields.index('sample')
        reader=csv.DictReader(original,delimiter='\t');q.need(reader.fieldnames==fields,'Source MT table schema differs from VCF')
        digest.update(codec.canonical_tsv(fields))
        for line in stream:
            guard();record=line.rstrip('\r\n').split('\t')
            encoded=[t.split('=',1)[1] for t in record[7].split(';') if t.startswith('WGS_MT_CONSEQUENCE_EVIDENCE=')]
            q.need(len(encoded)==1,'One complete MT evidence packet required')
            packet=codec.decode_packet(encoded[0]);q.need(tuple(packet['allele'])==q.key(record),'MT packet allele differs')
            chrom='chrM' if record[0] in ['M','MT','chrM'] else record[0]
            q.need(chrom=='chrM','MT views cannot include nuclear rows')
            pos=int(record[1]);ref=record[3];alts=record[4].split(',')
            fmt=dict(zip(record[8].split(':'),record[9].split(':')))
            for index,alt in enumerate(alts,1):
                # Identity comes from the source ALT vector, not diploid GT or inferred AF.
                q.need(re.fullmatch('[ACGT]+',ref) and re.fullmatch('[ACGT]+',alt),'Unsupported source allele: keep source, refuse exact-sharing claim')
                observed_allele(alleles,(chrom,pos,ref,alt),sample,dict(ALT_vector=record[4],ALT_index_1=index,
                    GT=fmt.get('GT','.'),AF_vector=fmt.get('AF','.'),HL_raw=fmt.get('HL','.'),AD_vector=fmt.get('AD','.'),DP_raw=fmt.get('DP','.'),
                    FILTER=record[6],QUAL=record[5],FORMAT=record[8],raw_sample=record[9]))
            for packed in packet['rows']:
                guard();values=packed['values'];source_row=next(reader,None)
                q.need(packed['i']==rows_total and len(values)==len(fields) and values[sample_index]==sample,'MT packet row ordinal/shape/sample differs')
                q.need(source_row is not None and [source_row[f] for f in fields]==values,'Actual VCF vectors differ from original MT functional table')
                row=dict(zip(fields,values));q.need(q.row_key(row)==q.key(record),'MT evidence row allele vector differs')
                q.need(row['source_FORMAT']==record[8] and row['source_sample']==record[9],'MT evidence original heteroplasmy fields differ')
                for writer in master_writers:writer.writerow(values)
                digest.update(codec.canonical_tsv(values));rows_total+=1
                number=row.get('allele_number','');matched=row.get('matched_alt','')
                if number.isdigit() and 1<=int(number)<=len(alts) and matched==alts[int(number)-1]:
                    attach_effect(alleles[(chrom,pos,ref,matched)],row)
                elif len(alts)==1:
                    attach_effect(alleles[(chrom,pos,ref,alts[0])],row)
                else:
                    unresolved+=1
            records+=1
        q.need(next(reader,None) is None and digest.hexdigest()==expected_digest,'MT source table/packet complete-vector digest differs')
    return dict(sample=sample,records=records,transcript_rows=rows_total,unresolved_multiallelic_transcript_rows=unresolved,
                ordered_vector_sha256=digest.hexdigest(),fields=fields)


def produce(q,p,manifest,receipt,output,guard):
    codec=q.module(p['codec']['path'],'_mt_functional_codec');root=Path(p['input_root'])
    samples=[r['sample'] for r in manifest['samples']]
    q.need(len(samples)==len(set(samples))==4,'Exactly4sourceMTsamples required')
    labels={sample:'sample_%04d'%(i+1) for i,sample in enumerate(samples)}
    output.mkdir(mode=0o700,parents=True,exist_ok=False);(output/'masters').mkdir();(output/'tables').mkdir()
    alleles={};proofs=[];fields=None
    cohort_stream=gzip.open(output/'masters/cohort.master.tsv.gz','wt',encoding='utf-8',newline='')
    cohort_writer=csv.writer(cohort_stream,delimiter='\t',lineterminator='\n')
    try:
        for sample_result in manifest['samples']:
            sample=sample_result['sample'];vcf=root/sample/'mitochondria.fully-consequence-annotated.vcf.gz';table=root/sample/'MT-functional-all-transcripts.tsv.gz'
            for path in [vcf,Path(str(vcf)+'.tbi'),table]:
                item=receipt['artifacts'][str(path.relative_to(root))];q.verify(item)
            expected_fields=sample_result['packet_proof']['fields']
            if fields is None:fields=expected_fields;cohort_writer.writerow(fields)
            q.need(fields==expected_fields,'Per-sample MT schemas differ')
            with gzip.open(output/'masters'/(labels[sample]+'.master.tsv.gz'),'wt',encoding='utf-8',newline='') as stream:
                writer=csv.writer(stream,delimiter='\t',lineterminator='\n');writer.writerow(fields)
                proofs.append(collect(q,codec,vcf,table,sample,sample_result['packet_proof']['ordered_vector_sha256'],alleles,[cohort_writer,writer],guard))
    finally:cohort_stream.close()
    summaries=[summarize(key,value,samples) for key,value in sorted(alleles.items())]
    terms=sorted(set(CORE)|set().union(*(value['terms'] for value in alleles.values())))
    if any(not value['terms'] for value in alleles.values()):terms.append('unannotated_or_allele_unresolved')
    index_rows=[];scopes=['cohort','shared_ge2','all4',*samples]
    for scope in scopes:
        for level in ['observed','site_PASS']:
            guard();folder=output/'tables'/labels.get(scope,scope)/level;folder.mkdir(parents=True)
            selected=[row for row in summaries if belongs(scope,level,row,samples)]
            targets=[('all.functional',selected)]
            for term in terms:targets.append((term,[row for row in selected if term in json.loads(row['consequences_json'])]))
            for name,rows in targets:
                slug=name if name=='all.functional' else re.sub('[^A-Za-z0-9_-]','_',name)[:80]+'-'+hashlib.sha256(name.encode()).hexdigest()[:8]
                path=folder/(slug+'.tsv.gz')
                with gzip.open(path,'wt',encoding='utf-8',newline='') as stream:
                    writer=csv.DictWriter(stream,SUMMARY,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
                index_rows.append(dict(scope=labels.get(scope,scope),sample_id=scope if scope in samples else '',level=level,
                    consequence=name,observed_unique_alleles=len(rows),relative_path=str(path.relative_to(output)),
                    zero_semantics='0observed_in_this_source_scope; not_missing_data_or_proven_biological_absence'))
    index=output/'tables/functional-category-index.tsv'
    with index.open('x',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,list(index_rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(index_rows)
    page='<!doctype html><meta charset="utf-8"><title>MT functional views</title><h1>Separate mitochondrial functional views</h1><p>Source Mutect2 call-record observations and literal site PASS are separate. Missing is not reference. No nuclear GQ, diploid carrier model, LoF/GoF/PVS1 or clinical burden assignment.</p><a href="masters/cohort.master.tsv.gz">All original transcript evidence decoded from MT VCFs</a> · <a href="tables/functional-category-index.tsv">All functional types, scopes and explicit zero counts</a>'
    (output/'index.html').write_text(page,encoding='utf-8')
    final=dict(schema_version=1,status='PASS_MT_FUNCTIONAL_VIEWS',samples=samples,sample_file_labels=labels,
        producer=p['producer'],reviewed_profile_sha256=p['reviewed_profile_sha256'],source_manifest=q.pin(root/'manifest.json'),source_receipt=q.pin(root/'workflow-receipt.json'),
        proofs=proofs,unique_exact_alleles=len(alleles),functional_category_entries=len(index_rows),
        source_authority='ACTUAL_MT_VCF_PACKETS_DECODED_AND_ALL_VECTOR_ROWS_MATCH_SOURCE_TABLES',
        shared_identity='EXACT_chrM_POS_REF_matchedALT; ALLELE_NUM_for_transcript_projection; raw_ALT_AF_HL_vectors_retained',
        quality_policy='Observed source call record independent of nuclear GT; literal site FILTER PASS separate; no nuclear GQ or candidate_eligible',
        automatic_ACMG=False,LoF_GoF_PVS1_assigned=False,clinical_burden_assigned=False,
        artifacts={str(f.relative_to(output)):q.pin(f) for f in output.rglob('*') if f.is_file()})
    q.write(output/'manifest.json',final)
    q.write(output/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='mitochondrial-functional-views-v1',
        artifacts={str(f.relative_to(output)):q.pin(f) for f in output.rglob('*') if f.is_file()}))
    return final


def prepare(args):
    q=helper(args.helper);end=dt.datetime.fromisoformat(args.deadline_utc.replace('Z','+00:00'))
    q.need(q.now()<end<=END,'MTfunctional deadline must be<=06:25UTC')
    source=q.load(args.input_profile)
    q.need(source['status']=='ROOT_REVIEWED_MT_CONSEQUENCE_ANNOTATION' and source['root_reviewed'] is True,'Reviewed MTconsequence source profile required')
    profile=dict(schema_version=1,status='CANDIDATE_MT_FUNCTIONAL_VIEWS_NOT_RUNNING',root_reviewed=False,
        producer=q.pin(__file__),helper=q.pin(args.helper),codec=q.pin(args.codec),input_profile=q.pin(args.input_profile),
        expected_input_producer=source['producer'],input_root=str(q.safe(source['output'])),output=str(q.safe(args.output)),
        deadline_utc=end.isoformat(),cpu_limit=2,memory_gib=4,min_free_gib=150,max_output_bytes=512*1024**2)
    Path(args.profile).parent.mkdir(parents=True,exist_ok=True);q.write(args.profile,profile)
    print(json.dumps(dict(status=profile['status'],profile=q.pin(args.profile)),indent=2))


def run(args):
    q=helper(args.helper);q.need(q.sha(args.profile)==args.reviewed_profile_sha256,'Exact reviewed profileSHA required');p=q.load(args.profile)
    if not args.execute:print('PLAN_ONLY_NO_MT_EXPORT');return
    q.need(p['status']=='ROOT_REVIEWED_MT_FUNCTIONAL_VIEWS' and p['root_reviewed'] is True,'Reviewed MTfunctional profile required')
    for name in ['producer','helper','codec','input_profile','expected_input_producer']:q.verify(p[name])
    q.need(p['producer']['sha256']==q.sha(__file__),'Exporter producer differs')
    end=dt.datetime.fromisoformat(p['deadline_utc']);q.need(q.now()<end<=END and p['cpu_limit']==2 and p['memory_gib']==4,'MTfunctional limits differ')
    state=Path(args.profile).parent/'mt-functional-status.json';output=q.safe(p['output']);q.need(not output.exists(),'Fresh MTfunctional output required')
    def guard():
        q.need(q.now()<end and __import__('shutil').disk_usage(output.parent).free>=150*1024**3,'MTfunctional deadline/reserve reached')
        if output.exists():q.need(sum(f.stat().st_size for f in output.rglob('*') if f.is_file())<=p['max_output_bytes'],'MTfunctional output cap reached')
    def expired(*unused):raise TimeoutError('MTfunctional deadline')
    signal.signal(signal.SIGALRM,expired);signal.setitimer(signal.ITIMER_REAL,(end-q.now()).total_seconds())
    os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    root=q.safe(p['input_root'])
    while not (root/'workflow-receipt.json').exists():
        guard();q.write(state,dict(status='WAITING_FOR_MT_CONSEQUENCE_PASS',pid=os.getpid()),replace=True);time.sleep(10)
    receipt=q.load(root/'workflow-receipt.json');manifest=q.load(root/'manifest.json')
    q.need(receipt['status']=='PASS' and receipt['stage']=='mitochondrial-consequence-annotation-v1' and manifest['status']=='PASS_MT_CONSEQUENCE_ANNOTATION','Actual MTconsequence completion required')
    q.verify(receipt['artifacts']['manifest.json'])
    q.need(manifest['producer']==p['expected_input_producer'] and manifest['reviewed_profile']==p['input_profile'],'MTconsequence provenance differs')
    p['reviewed_profile_sha256']=q.sha(args.profile)
    try:
        q.write(state,dict(status='EXPORTING_MT_ALL_FUNCTIONAL_TYPES',pid=os.getpid()),replace=True)
        result=produce(q,p,manifest,receipt,output,guard)
        q.write(state,dict(status=result['status'],pid=os.getpid(),unique_alleles=result['unique_exact_alleles']),replace=True)
    except BaseException as error:
        q.write(state,dict(status='FAILED_OR_BOUNDED_STOP',error_type=type(error).__name__),replace=True);raise
    print(json.dumps(dict(status=result['status'],unique_alleles=result['unique_exact_alleles']),indent=2))


def main():
    os.umask(0o077);parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for name in ['helper','codec','input-profile','output','profile']:a.add_argument('--'+name,type=Path,required=True)
    a.add_argument('--deadline-utc',required=True)
    a=sub.add_parser('run');a.add_argument('--helper',type=Path,required=True);a.add_argument('--profile',type=Path,required=True)
    a.add_argument('--reviewed-profile-sha256',required=True);a.add_argument('--execute',action='store_true')
    args=parser.parse_args();{'prepare':prepare,'run':run}[args.command](args)


if __name__=='__main__':main()
