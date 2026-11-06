#!/usr/bin/env python3
"""Offline MT consequences on unchanged Mutect2 calls; no nuclear genotype rules."""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import signal
import subprocess
import time

END = dt.datetime.fromisoformat('2026-11-06T06:25:00+00:00')
FIELDS = ['sample','chrom','pos','ref','alt','GT','AD','AF','HL','DP','QUAL','FILTER',
          'source_INFO','source_FORMAT','source_sample','allele_number','matched_alt',
          'gene','transcript','consequence','hgvs_c','hgvs_p','csq_json','frequency_context_json',
          'source_qc_json','genetic_code_status','heteroplasmy_scope','clinical_classification']


def helper(path):
    import importlib.util
    spec=importlib.util.spec_from_file_location('_mt_supplementary_helpers',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    result.TAG='WGS_MT_CONSEQUENCE_EVIDENCE';result.SCHEMA='WGS_MT_CONSEQUENCE_SCHEMA'
    return result


def sequence(reference):
    import pysam
    with pysam.FastaFile(str(reference)) as fasta:
        if fasta.get_reference_length('chrM') != 16569:
            raise ValueError('Exact16569bp chrM rCRS required')
        return fasta.fetch('chrM').upper()


def selected_resources(q,manifest_path):
    manifest=q.load(manifest_path);q.need(manifest.get('status')=='PASS' and manifest.get('assembly')=='GRCh38','Verified GRCh38 resources required')
    prefix='cache/homo_sapiens/115_GRCh38/'
    selected=[]
    for item in manifest['artifacts']:
        name=item['path']
        if name.startswith(prefix+'MT/') or (name.startswith(prefix) and '/' not in name[len(prefix):]):
            selected.append(dict(item,path=str(Path(manifest_path).parent/name)))
    q.need(any('/MT/' in item['path'] for item in selected),'Actual verified MT cache subtree required')
    for name in ['clinvar','clinvar_index']:
        item=manifest['resources'][name]
        q.need(item['status']=='available','ClinVar source/index unavailable')
        path=Path(item['path']);path=path if path.is_absolute() else Path(manifest_path).parent/path
        selected.append(dict(path=str(path),bytes=item['bytes'],sha256=item['sha256']))
    return manifest,selected


def prepare(args):
    q=helper(args.helper);cfg=q.load(args.base_config);context=q.load(args.mt_context)
    end=dt.datetime.fromisoformat(args.deadline_utc.replace('Z','+00:00'))
    q.need(q.now()<end<=END,'Separate MT consequence deadline must be<=06:25UTC')
    q.need(context['status']=='PASS_INDEPENDENT_MITOCHONDRIAL_REFERENCE_CONTEXT' and len(context['source_mt_inputs'])==4,'Successful4sampleMT reference context required')
    context_root=Path(args.mt_context).parent;context_receipt=q.load(context_root/'workflow-receipt.json')
    q.need(context_receipt['status']=='PASS','MT reference context completion required')
    seq=sequence(cfg['reference']['fasta'])
    source=q.load(Path(cfg['reference']['fasta']).parent/'reference-receipt.json')
    q.need(any(r.get('sha256')==cfg['reference']['sha256'] for r in source['sources']),'Original reference receipt lineage differs')
    entries=[]
    for entry in context['source_mt_inputs']:
        sample=entry['sample'];q.verify(entry['receipt']);q.verify(entry['vcf'])
        calls=q.load(entry['receipt']['path']);q.need(calls['status']=='PASS','Original MT calling stage not successful')
        q.need(calls['artifacts']['mitochondria.filtered.vcf.gz']['sha256']==entry['vcf']['sha256'],'Unshifted filteredMT source differs')
        table=context_root/'samples'/sample/'MT-population-reference.tsv.gz'
        expected=context_receipt['artifacts'][str(table.relative_to(context_root))]
        table_pin=q.pin(table);q.need(all(table_pin[k]==expected[k] for k in ['bytes','sha256']),'AppliedMT frequency context changed')
        entries.append(dict(sample=sample,vcf=entry['vcf'],receipt=entry['receipt'],qc=entry['qc'],frequency_table=table_pin))
    manifest,selected=selected_resources(q,cfg['annotation']['resource_manifest'])
    profile=dict(schema_version=1,status='CANDIDATE_MT_CONSEQUENCE_NOT_RUNNING',root_reviewed=False,
        producer=q.pin(__file__),helper=q.pin(args.helper),codec=q.pin(args.codec),base_config=q.pin(args.base_config),
        reference=cfg['reference'],fai=q.pin(cfg['reference']['fasta']+'.fai'),mt_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
        reference_receipt=q.pin(Path(cfg['reference']['fasta']).parent/'reference-receipt.json'),
        resource_manifest=q.pin(cfg['annotation']['resource_manifest']),verified_resource_scope='MT_CACHE_SUBTREE_GLOBAL_METADATA_AND_CLINVAR_ONLY',selected_resources=selected,
        vep_image=cfg['tools']['vep_image'],source_context_manifest=q.pin(args.mt_context),source_context_receipt=q.pin(context_root/'workflow-receipt.json'),
        source_versions=context['source_versions'],samples=entries,output=str(q.safe(args.output)),deadline_utc=end.isoformat(),
        cpu_limit=2,memory_gib=8,min_free_gib=150,no_recalling=True,no_nuclear_quality_filter=True,no_ploidy_assignment=True,
        genetic_code_requirement='NCBI_VERTEBRATE_MITOCHONDRIAL_TABLE2_ACTUAL_VEP_PUBLIC_CONTROL_REQUIRED',
        public_control_receipt=None)
    Path(args.profile).parent.mkdir(parents=True,exist_ok=True);q.write(args.profile,profile)
    print(json.dumps(dict(status=profile['status'],profile=q.pin(args.profile)),indent=2))


def vep(q,p,input_vcf,work,guard):
    resources=Path(p['resource_manifest']['path']).parent;ref=Path(p['reference']['fasta'])
    original_samples=[]
    with q.text(input_vcf) as stream:
        headers=q.header(stream);original_samples=headers[-1].split('\t')[9:]
        q.need(len(original_samples)==1,'Single-sample originalMT VCF required')
    copied=work/'input.vcf.gz'
    if str(input_vcf).endswith('.gz'):shutil.copyfile(input_vcf,copied)
    else:
        with Path(input_vcf).open('rb') as source,gzip.open(copied,'wb') as output:shutil.copyfileobj(source,output)
    source_manifest=q.load(p['resource_manifest']['path']);clinvar=source_manifest['resources']['clinvar']
    cv=Path(clinvar['path']);cv=cv if cv.is_absolute() else resources/cv
    name='marina-mt-consequence-'+str(os.getpid())+'-'+str(time.time_ns())
    args=['docker','run','--rm','--pull=never','--name',name,'--network=none','--cpus','2','--memory','8g',
          '--user',str(os.getuid())+':'+str(os.getgid()),
          '--mount','type=bind,src='+str(work.resolve())+',dst=/work',
          '--mount','type=bind,src='+str(resources)+',dst=/resources,readonly',
          '--mount','type=bind,src='+str(ref.parent)+',dst=/reference,readonly',p['vep_image'],
          'vep','--offline','--cache','--dir_cache','/resources/cache','--cache_version','115','--assembly','GRCh38',
          '--species','homo_sapiens','--input_file','/work/input.vcf.gz','--output_file','/work/vep.vcf.gz',
          '--format','vcf','--vcf','--compress_output','bgzip','--fasta','/reference/'+ref.name,'--check_ref',
          '--hgvs','--hgvsg','--symbol','--numbers','--canonical','--transcript_version','--protein','--biotype',
          '--allele_number','--total_length','--variant_class','--check_existing','--no_stats','--fork','2',
          '--custom','file=/resources/'+cv.relative_to(resources).as_posix()+',short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNSIGCONF%CLNDN%CLNDISDB%ALLELEID']
    proc=None
    try:
        with (work/'offline-vep.log').open('xb') as log:
            guard();proc=subprocess.Popen(args,stdout=log,stderr=log,start_new_session=True)
            proc.wait(timeout=min(1800,(dt.datetime.fromisoformat(p['deadline_utc'])-q.now()).total_seconds()))
            q.need(proc.returncode==0,'MT offlineVEP failed; retain originalsource and private log')
    finally:
        if proc is not None and proc.poll() is None:os.killpg(proc.pid,signal.SIGKILL)
        subprocess.run(['docker','rm','-f',name],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=10)
    guard();return work/'vep.vcf.gz'


def make_mt_table(q,source,annotated,sample,frequency_table,qc,table):
    frequencies={}
    if frequency_table:
        with q.text(frequency_table) as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                q.need(row['sample']==sample,'Frequency table contains another sample')
                key=(row['CHROM'].removeprefix('chr'),int(row['POS']),row['REF'],row['ALT'])
                q.need(key not in frequencies,'Duplicate originalMT frequency record');frequencies[key]=row
    with q.text(source) as original,q.text(annotated) as ann,gzip.open(table,'wt',encoding='utf-8',newline='') as output:
        original_header,annotated_header=q.header(original),q.header(ann)
        q.need(original_header[-1].split('\t')[9:]==[sample],'OriginalMT sample differs')
        declaration=next(line for line in annotated_header if line.startswith('##INFO=<ID=CSQ,'))
        csq_fields=re.search(r'Format: ([^"\n]+)',declaration).group(1).split('|')
        writer=csv.DictWriter(output,FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
        count=0
        for line in original:
            src=line.rstrip('\r\n').split('\t');new=next(ann,None)
            q.need(new is not None,'MTVEP skipped originalrecord')
            record=new.rstrip('\r\n').split('\t')
            q.need(q.key(src)==q.key(record) and src[8:]==record[8:],'MT annotation changed GT/AF/HL/FORMAT')
            info={v.split('=',1)[0]:(v.split('=',1)[1] if '=' in v else True) for v in q.info_tokens(record[7])}
            source_fmt=dict(zip(src[8].split(':'),src[9].split(':')))
            frequency=frequencies.get(q.key(src))
            if frequency_table:
                q.need(frequency is not None and frequency['raw_INFO']==src[7] and frequency['raw_FORMAT']==src[8] and frequency['raw_sample']==src[9], 'AppliedMT frequency source no longer matches originalrawfields')
            csqs=info.get('CSQ','').split(',') if info.get('CSQ') else ['']
            for item in csqs:
                csq=dict(zip(csq_fields,item.split('|')))
                number=csq.get('ALLELE_NUM','');alts=src[4].split(',')
                matched=alts[int(number)-1] if number.isdigit() and 1<=int(number)<=len(alts) else ''
                writer.writerow(dict(sample=sample,chrom=src[0],pos=src[1],ref=src[3],alt=src[4],
                    **{k:source_fmt.get(k,'.') for k in ['GT','AD','AF','HL','DP']},QUAL=src[5],FILTER=src[6],
                    source_INFO=src[7],source_FORMAT=src[8],source_sample=src[9],allele_number=number,matched_alt=matched,
                    gene=csq.get('SYMBOL',''),transcript=csq.get('Feature',''),consequence=csq.get('Consequence',''),
                    hgvs_c=csq.get('HGVSc',''),hgvs_p=csq.get('HGVSp',''),csq_json=json.dumps(csq,ensure_ascii=False),
                    frequency_context_json=json.dumps(frequency or {'status':'PUBLIC_CONTROL_NO_FREQUENCY_QUERY'},ensure_ascii=False),
                    source_qc_json=json.dumps(qc,ensure_ascii=False),genetic_code_status='VERTEBRATE_MITOCHONDRIAL_CODE2_PUBLIC_CONTROL_VERIFIED',
                    heteroplasmy_scope='SOURCE_RAW_AF_HL_AD_DP_VECTOR_UNCHANGED; no_nuclear_GQ_or_diploid_reclassification',clinical_classification='NOT_ASSIGNED'))
            count+=1
        q.need(next(ann,None) is None,'MTVEP added unexpectedrecord')
    return count


def run(args):
    q=helper(args.helper);q.need(q.sha(args.profile)==args.reviewed_profile_sha256,'ExactROOT-reviewed profileSHA required');p=q.load(args.profile)
    if not args.execute:print('PLAN_ONLY_NO_MT_ANNOTATION');return
    q.need(p['root_reviewed'] is True and p['status']=='ROOT_REVIEWED_MT_CONSEQUENCE_ANNOTATION','ROOT-reviewed profile required')
    end=dt.datetime.fromisoformat(p['deadline_utc']);q.need(q.now()<end<=END and p['cpu_limit']==2 and p['memory_gib']==8,'MTdeadline/resource caps differ')
    for key in ['producer','helper','codec','base_config','fai','reference_receipt','resource_manifest','source_context_manifest','source_context_receipt']:
        q.verify(p[key])
    q.need(p['producer']['sha256']==q.sha(__file__) and p['helper']['path']==str(Path(args.helper).absolute()),'MTproducer/helper differs')
    proof=q.load(q.verify(p['public_control_receipt']));q.need(proof['status']=='PASS_ACTUAL_MT_CODE2_OFFLINE_CONTROL' and proof['producer']['sha256']==p['producer']['sha256'] and proof['vep_image']==p['vep_image'],'Actualmitochondrialcode2 VEP proof required')
    output=q.safe(p['output']);q.need(not output.exists(),'Fresh MTconsequence output required')
    output.mkdir(mode=0o700,parents=True)
    def guard():
        q.need(q.now()<end and shutil.disk_usage(output).free>=150*1024**3,'MTdeadline/disk reserve reached')
    def expired(*unused):raise TimeoutError('SeparateMT consequence deadline')
    signal.signal(signal.SIGALRM,expired);signal.setitimer(signal.ITIMER_REAL,(end-q.now()).total_seconds())
    os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]));resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    for item in p['selected_resources']:guard();q.verify(item)
    seq=sequence(p['reference']['fasta']);q.need(hashlib.sha256(seq.encode()).hexdigest()==p['mt_sequence_sha256'],'MTreference sequence changed')
    codec=q.module(p['codec']['path'],'_mt_evidence_codec');results=[]
    try:
        for entry in p['samples']:
            guard();source=q.verify(entry['vcf']);q.verify(entry['receipt']);frequency=q.verify(entry['frequency_table'])
            work=output/entry['sample'];work.mkdir(mode=0o700)
            q.write(output/'status.json',dict(status='RUNNING_MT_CONSEQUENCE',sample=entry['sample'],pid=os.getpid()),replace=True)
            qc={}
            (work/'source-QC').mkdir()
            for i,item in enumerate(entry['qc']):
                path=q.verify(item);content=path.read_text();qc[str(i)+'-'+path.name]=content
                shutil.copyfile(path,work/'source-QC'/(str(i)+'-'+path.name+'.txt'))
            annotation=vep(q,p,source,work,guard);table=work/'MT-functional-all-transcripts.tsv.gz'
            count=make_mt_table(q,source,annotation,entry['sample'],frequency,qc,table)
            final=work/'mitochondria.fully-consequence-annotated.vcf.gz'
            evidence=q.embed(source,annotation,table,final,codec,guard)
            results.append(dict(sample=entry['sample'],records=count,packet_proof=evidence,source=entry['vcf']))
        manifest=dict(schema_version=1,status='PASS_MT_CONSEQUENCE_ANNOTATION',samples=results,producer=p['producer'],
            reviewed_profile=q.pin(args.profile),public_control_receipt=p['public_control_receipt'],source_context=p['source_context_manifest'],
            original_mutect2_INFO_FORMAT_GT_AF_HL_FILTER_preserved=True,no_nuclear_GQ_or_ploidy_policy_applied=True,
            no_recalling=True,no_automatic_ACMG=True,AlphaMissense='NOT_APPLIED_NO_MT_MODEL',SpliceAI='NOT_APPLIED_NO_MT_MODEL',
            source_versions=p['source_versions'],resource_scope=p['verified_resource_scope'])
        q.write(output/'manifest.json',manifest)
        q.write(output/'status.json',dict(status=manifest['status'],pid=os.getpid()),replace=True)
        q.write(output/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='mitochondrial-consequence-annotation-v1',
            artifacts={str(f.relative_to(output)):q.pin(f) for f in output.rglob('*') if f.is_file()}))
    except BaseException as error:
        q.write(output/'status.json',dict(status='FAILED_OR_BOUNDED_STOP',error_type=type(error).__name__),replace=True);raise
    print(json.dumps(manifest,indent=2))


def main():
    os.umask(0o077);parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for name in ['helper','base-config','mt-context','codec','output','profile']:a.add_argument('--'+name,type=Path,required=True)
    a.add_argument('--deadline-utc',required=True)
    a=sub.add_parser('run');a.add_argument('--helper',type=Path,required=True);a.add_argument('--profile',type=Path,required=True)
    a.add_argument('--reviewed-profile-sha256',required=True);a.add_argument('--execute',action='store_true')
    args=parser.parse_args();{'prepare':prepare,'run':run}[args.command](args)


if __name__=='__main__':main()
