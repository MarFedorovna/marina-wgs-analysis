#!/usr/bin/env python3
"""Separate, reviewed X/Y germline calling. Default: validate/plan, never call."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import random
import shutil
import signal
import subprocess
import threading
import time
import uuid

END = dt.datetime(2026, 11, 6, 6, 20, tzinfo=dt.timezone.utc)
PARTS = ('X_NONPAR', 'X_PAR', 'Y_NONPAR')
MODEL = 'USER_AUTHORIZED_INFERRED_FROM_QC_NOT_CLINICAL_KARYOTYPE'

def need(ok, message):
    if not ok: raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024**2),b''): h.update(b)
    return h.hexdigest()

def load(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def now(): return dt.datetime.now(dt.timezone.utc)
def atomic(path, value):
    path=Path(path); tmp=path.with_name(path.name+'.tmp-'+uuid.uuid4().hex)
    with tmp.open('x') as f: json.dump(value,f,indent=2);f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)

def pin(path):
    p=Path(path);return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))

def checked_pin(value):
    p=Path(value['path']);need(p.is_absolute() and p.is_file() and not p.is_symlink(),'Pinned file absent or unsafe')
    need(p.stat().st_size==value['bytes'] and sha(p)==value['sha256'],'Pinned small file changed')
    return p

def identity(value):
    p=Path(value['path']);s=p.stat()
    need(p.is_absolute() and p.is_file() and not p.is_symlink() and s.st_size==value['bytes'] and s.st_mtime_ns==value['mtime_ns'],'Large input identity changed; no large input hash performed')
    need(re.fullmatch('[0-9a-f]{64}',value['sha256']),'Historical expected SHA missing')
    return p

def checked_profile(p):
    need(p['schema_version']==1 and p['status']=='ROOT_REVIEWED_XY_GERMLINE_V1' and p['processing_model']==MODEL,'Unreviewed processing model')
    end=dt.datetime.fromisoformat(p['deadline_utc']);need(now()<end<=END and (end-now()).total_seconds()<=7200,'Expired or excessive separate calling window')
    need(p['resources']==dict(call_workers=6,call_cpus=4,call_memory_gib=16,max_cpu=24,max_memory_gib=120),'Unreviewed resource limits')
    need(re.fullmatch(r'[^\s@]+@sha256:[0-9a-f]{64}',p['gatk_image']),'Immutable CPU GATK image required')
    samples=p['samples'];ids=[s['sample'] for s in samples]
    need(len(ids)==4 and len(set(ids))==4 and all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}',s) for s in ids),'Exactly four explicit samples required')
    need(sorted((s['X_NONPAR'],s['X_PAR'],s['Y_NONPAR']) for s in samples)==[(1,2,1),(1,2,1),(2,2,0),(2,2,0)],'Expected reviewed 2+2 processing models')
    qc=load(checked_pin(p['qc_receipt']));need(qc.get('status')=='PASS','Technical QC must have PASS receipt')
    contract=load(checked_pin(p['alignment_contract']))
    need(contract['status']=='PASS' and contract['processing_model']==MODEL and contract['reference_sha256']==p['reference']['sha256'] and contract['par_policy']=='Y_PAR_HARD_MASKED_AT_ALIGNMENT' and contract['mask_present_at_original_alignment'] is True and contract['par_report_coordinates']=='X_only','Original alignment/PAR contract unverified')
    mask=load(checked_pin(p['reference_mask_proof']))
    need(mask['status']=='PASS_REFERENCE_IDENTITY' and mask['fasta_sha256']==p['reference']['sha256'] and mask['Y_PAR_all_hard_masked_N'] is True,'Actual reference/mask proof unverified')
    identity(p['reference']);checked_pin(p['reference_fai']);checked_pin(p['reference_dict'])
    need(p['reference_fai']['path']==p['reference']['path']+'.fai','FAI path mismatch')
    for s in samples:
        identity(s['bam']);checked_pin(s['bai']);r=load(checked_pin(s['alignment_receipt']))
        need(r['sample']==s['sample'] and r['config']['reference']['sha256']==p['reference']['sha256'],'Original BAM alignment reference mismatch')
        need(r['outputs'][Path(s['bam']['path']).name]==s['bam']['sha256'],'Historical BAM SHA mismatch')
        need(contract['sample_alignment_receipt_sha256'][s['sample']]==s['alignment_receipt']['sha256'],'Alignment contract does not bind this sample')
        need(s['bai']['path']==s['bam']['path']+'.bai','BAM index path mismatch')
    module=checked_pin(p['filter_module']);need(module.name=='calling_stages.py','Expected reviewed filter module')
    need(15<=p['min_free_gib']<=150,'Explicit disk reserve required')
    checked_pin(p['par_lock']);checked_pin(p['executor'])
    need(p['executor']['sha256']==sha(__file__),'Executor identity mismatch')
    need(Path(p['reference_dict']['path']).parent==Path(p['reference']['path']).parent,'Reference dictionary must be alongside FASTA')
    lock=load(p['par_lock']['path'])
    need(lock['assembly']=='GRCh38' and lock['primary_contig_lengths']=={'chrX':156040895,'chrY':57227415},'GRCh38 contig contract mismatch')
    need(lock['par_one_based_closed']=={'chrX':[[10001,2781479],[155701383,156030895]],'chrY':[[10001,2781479],[56887903,57217415]]},'PAR coordinates changed')
    need(contract['official_par_gff_sha256']==lock['retrieved_gff_sha256'],'PAR official-source mismatch')
    for key in ('bcftools','samtools'):checked_pin(p[key])
    out=Path(p['output']);need(out.is_absolute() and not out.is_symlink() and '/outputs/' in out.as_posix() and out.name.startswith('xy-germline-'),'New separate XY output required')
    return end

def intervals(lock):
    out={'X_PAR':[('chrX',a,b) for a,b in lock['par_one_based_closed']['chrX']]}
    for chrom,label in [('chrX','X_NONPAR'),('chrY','Y_NONPAR')]:
        start=1;rows=[]
        for a,b in lock['par_one_based_closed'][chrom]:
            if start<a:rows.append((chrom,start,a-1))
            start=b+1
        if start<=lock['primary_contig_lengths'][chrom]:rows.append((chrom,start,lock['primary_contig_lengths'][chrom]))
        out[label]=rows
    return out

class Runner:
    def __init__(self,p,out,end):self.p=p;self.out=out;self.end=end;self.mutex=threading.Lock();self.active=set()
    def command(self,args,label,capture=False):
        need(now()<self.end,'Absolute calling deadline');need(shutil.disk_usage(self.out).free>=self.p['min_free_gib']*1024**3,'Disk reserve exhausted')
        with self.mutex:
            with (self.out/'commands.jsonl').open('a') as f:f.write(json.dumps(dict(label=label,argv=args))+'\n')
        log=self.out/(label+'.log')
        with log.open('ab') as f:
            r=subprocess.run(args,stdout=subprocess.PIPE if capture else f,stderr=f,timeout=max(.01,(self.end-now()).total_seconds()),check=True)
        return r.stdout.decode() if capture else None
    def gatk(self,args,label,cpus=4,memory=16,heap=12):
        roots={str(Path(self.p['reference']['path']).parent),str(self.out)}
        roots.update(str(Path(s['bam']['path']).parent) for s in self.p['samples'])
        name='marina-xy-'+uuid.uuid4().hex;cmd=['docker','run','--rm','--network','none','--name',name,'--user',str(os.getuid())+':'+str(os.getgid()),'--cpus',str(cpus),'--memory',str(memory)+'g']
        for root in sorted(roots):
            need(',' not in root,'Unsupported Docker mount path')
            cmd+=['--mount','type=bind,src='+root+',dst='+root+('' if root==str(self.out) else ',readonly')]
        cmd+=[self.p['gatk_image'],'gatk','--java-options',f'-Xmx{heap}g -XX:ParallelGCThreads=2',*args]
        try:self.command(cmd,label)
        finally:
            # Only a uniquely named container created by this invocation may be stopped.
            subprocess.run(['docker','stop','--time','2',name],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=10)
    def status(self,status,**kw):
        with self.mutex:atomic(self.out/'status.json',dict(status=status,pid=os.getpid(),deadline_utc=self.end.isoformat(),updated_utc=now().isoformat(),**kw))

def vcf_check(path,samples,ploidies,regions,allow_missing=()):
    opener=gzip.open if str(path).endswith('.gz') else open;count=0
    with opener(path,'rt') as f:
        header=None
        for line in f:
            if line.startswith('#CHROM'):header=line.rstrip().split('\t')[9:];need(header==samples,'VCF sample identity/order mismatch')
            elif not line.startswith('#'):
                cols=line.rstrip('\n').split('\t');need(header is not None and len(cols)==9+len(samples),'Malformed VCF')
                pos=int(cols[1]);ends=[int(x[4:]) for x in cols[7].split(';') if x.startswith('END=')];right=ends[0] if ends else pos+len(cols[3])-1
                need(right>=pos and any(cols[0]==c and a<=pos<=right<=b for c,a,b in regions),'Record or reference block outside declared partition')
                keys=cols[8].split(':');need('GT' in keys,'GT absent')
                for sid,value in zip(samples,cols[9:]):
                    gt=value.split(':')[keys.index('GT')];als=re.split('[/|]',gt)
                    if sid in allow_missing:need(all(a=='.' for a in als),'Excluded Y sample has invented genotype')
                    elif any(a!='.' for a in als):need(len(als)==ploidies[sid],'Called GT ploidy mismatch')
                count+=1
    need(header==samples,'Missing VCF header');return count

def check_final(path,ids,samples,regions):
    count=0;last=None;seen_at_position=set();last_coord=None
    with gzip.open(path,'rt') as f:
        header=None
        for line in f:
            if line.startswith('#CHROM'):header=line.rstrip().split('\t')[9:];need(header==ids,'Final sample order mismatch')
            elif not line.startswith('#'):
                c=line.rstrip('\n').split('\t');need(header==ids and len(c)==9+len(ids),'Final malformed record')
                coordinate=(('chrX','chrY').index(c[0]),int(c[1]));need(last_coord is None or coordinate>=last_coord,'Final chromosome/position order mismatch')
                if coordinate!=last_coord:seen_at_position=set()
                key=tuple(c[:2]+c[3:5]);need(key not in seen_at_position,'Duplicate final exact allele');seen_at_position.add(key);last_coord=coordinate
                part=[name for name,rows in regions.items() if any(c[0]==chrom and a<=int(c[1])<=b for chrom,a,b in rows)];need(len(part)==1,'Final record lacks unique partition')
                gtindex=c[8].split(':').index('GT')
                for sample,value in zip(samples,c[9:]):
                    alleles=re.split('[/|]',value.split(':')[gtindex]);ploidy=sample[part[0]]
                    if not ploidy:need(all(a=='.' for a in alleles),'Final excluded Y sample has GT')
                    elif any(a!='.' for a in alleles):need(len(alleles)==ploidy,'Final mixed GT ploidy differs')
                count+=1
    need(header==ids,'Final VCF header absent');return count

def module(path):
    spec=importlib.util.spec_from_file_location('xy_reviewed_filters',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def checkpoint(r,name,identity_value,produce):
    target=r.out/name;receipt=target/'receipt.json'
    if target.exists():
        need(receipt.is_file(),'Existing task incomplete; only a new attempt directory may be used')
        saved=load(receipt);need(saved['status']=='PASS' and saved['identity']==identity_value,'Task receipt identity mismatch')
        for x in saved['artifacts']:checked_pin(x)
        return target
    stage=r.out/(name+'.attempt-'+uuid.uuid4().hex);stage.mkdir(parents=True)
    artifacts=produce(stage)
    need(now()<r.end,'Deadline reached before task promotion')
    atomic(stage/'receipt.json',dict(status='PASS',identity=identity_value,artifacts=[dict(pin(x),path=str(target/x.relative_to(stage))) for x in artifacts]))
    os.rename(stage,target);return target

def execute(p,end):
    import fcntl
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=True);guard=(out/'owns.lock').open('a');fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (out/'manifest.json').exists():
        saved=load(out/'manifest.json');need(saved['status']=='PASS_XY_GERMLINE_CALLING' and saved['profile_sha256']==hashlib.sha256(json.dumps(p,sort_keys=True).encode()).hexdigest(),'Existing completion identity differs')
        for artifact in saved['artifacts']:checked_pin(artifact)
        need(load(out/'workflow-receipt.json')['status']=='PASS','Completion receipt absent');guard.close();return saved
    r=Runner(p,out,end);rows=intervals(load(p['par_lock']['path']));ids=[s['sample'] for s in p['samples']];ref=p['reference']['path'];bc=p['bcftools']['path'];sam=p['samtools']['path'];m=module(p['filter_module']['path'])
    profile_identity=hashlib.sha256(json.dumps(p,sort_keys=True).encode()).hexdigest()
    def alarm(*_):raise TimeoutError('Whole XY executor deadline')
    old=signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.01,(end-now()).total_seconds()))
    try:
        profile_copy=out/'reviewed-profile.json'
        if profile_copy.exists():need(load(profile_copy)==p,'Stored profile changed')
        else:atomic(profile_copy,p)
        r.status('PREFLIGHT')
        version=r.command(['docker','image','inspect',p['gatk_image'],'--format','{{json .RepoDigests}}'],'image-digest',True);need(p['gatk_image'] in json.loads(version),'Installed GATK image differs')
        # An independently executed synthetic mixed-ploidy control must bind this exact image/code.
        control=load(checked_pin(p['public_control']))
        need(control.get('multiallelic_removed_alt_missing_verified') is True and control['status']=='PASS_XY_MIXED_PLOIDY_CONTROL' and control['gatk_image']==p['gatk_image'] and control['executor_sha256']==sha(__file__) and control['bcftools']==p['bcftools'] and control['samtools']==p['samtools'] and control['haploid_gt_verified'] and control['diploid_gt_verified'],'Public mixed-ploidy control not proven')
        for s in p['samples']:
            r.command([sam,'quickcheck','-v',s['bam']['path']],'bam-check-'+s['sample'])
            actual=r.command([sam,'view','-H',s['bam']['path']],'bam-header-'+s['sample'],True)
            sm={piece[3:] for line in actual.splitlines() if line.startswith('@RG\t') for piece in line.split('\t') if piece.startswith('SM:')}
            need(sm=={s['sample']},'BAM read-group sample mismatch')
        for part,regions in rows.items():
            path=out/(part+'.intervals');text=''.join(f'{c}:{a}-{b}\n' for c,a,b in regions)
            if path.exists():need(path.read_text()==text,'Existing intervals differ')
            else:path.write_text(text)
        r.status('HAPLOTYPECALLER_RUNNING',task_count=10)
        def call(s,part):
            name='call-'+part+'-'+s['sample'];ident=dict(profile=profile_identity,sample=s['sample'],partition=part,ploidy=s[part])
            def produce(d):
                f=d/'variants.g.vcf.gz'
                r.gatk(['HaplotypeCaller','-R',ref,'-I',s['bam']['path'],'-L',str(out/(part+'.intervals')),'--interval-padding','0','--sample-ploidy',str(s[part]),'--native-pair-hmm-threads','4','-ERC','GVCF','-O',str(f)],name)
                r.command([bc,'index','--force','--tbi',str(f)],name+'-index');vcf_check(f,[s['sample']],{s['sample']:s[part]},rows[part]);return [f,Path(str(f)+'.tbi')]
            return (part,s['sample'],checkpoint(r,name,ident,produce)/'variants.g.vcf.gz')
        jobs=[(s,part) for part in PARTS for s in p['samples'] if s[part]>0];called={}
        with ThreadPoolExecutor(max_workers=6) as pool:
            fs=[pool.submit(call,s,part) for s,part in jobs]
            for f in as_completed(fs):part,sample,path=f.result();called[part,sample]=path;r.status('HAPLOTYPECALLER_RUNNING',completed_calls=len(called),task_count=10)
        final_parts=[]
        for part in PARTS:
            selected=[s for s in p['samples'] if s[part]>0];ploidy={s['sample']:s[part] for s in selected}
            def joint(d):
                combined=d/'combined.g.vcf.gz';raw=d/'joint.raw.vcf.gz';norm=d/'joint.normalized.vcf.gz';expanded=d/'expanded.vcf';tagged=d/'all.soft-filtered.vcf';chosen=d/'selected.vcf'
                args=['CombineGVCFs','-R',ref,'-L',str(out/(part+'.intervals'))]
                for s in selected:args+=['-V',str(called[part,s['sample']])]
                r.gatk(args+['-O',str(combined)],part+'-combine',4,24,18)
                r.gatk(['GenotypeGVCFs','-R',ref,'-V',str(combined),'-L',str(out/(part+'.intervals')),'--standard-min-confidence-threshold-for-calling','0','-O',str(raw)],part+'-genotype',4,24,18)
                vcf_check(raw,sorted(s['sample'] for s in selected),ploidy,rows[part])
                r.command([bc,'norm','-m','-any','--multi-overlaps','.','-f',ref,'-c','e','-Oz','-o',str(norm),str(raw)],part+'-normalize')
                m.expand_sample_columns(norm,expanded,ids)
                m.soft_filter_vcf(expanded,tagged,chosen)
                result=d/'all.soft-filtered.vcf.gz';r.command([bc,'view','-Oz','-o',str(result),str(tagged)],part+'-compress');r.command([bc,'index','--tbi',str(result)],part+'-index')
                vcf_check(result,ids,ploidy,rows[part],[s for s in ids if s not in ploidy]);return [result,Path(str(result)+'.tbi'),raw,Path(str(raw)+'.tbi')]
            r.status('JOINT_NORMALIZE_SOFTFILTER',partition=part)
            d=checkpoint(r,'joint-'+part,dict(profile=profile_identity,partition=part),joint);final_parts.append(d/'all.soft-filtered.vcf.gz')
        merged=out/'cohort.xy.soft-filtered.vcf.gz';unsorted=out/'cohort.unsorted.vcf.gz'
        need(not merged.exists(),'Final output already exists; verify its completed receipt instead')
        r.command([bc,'concat','-a','-Oz','-o',str(unsorted),*map(str,final_parts)],'concat')
        r.command([bc,'sort','-Oz','-o',str(merged),str(unsorted)],'sort');r.command([bc,'index','--tbi',str(merged)],'final-index')
        count=0
        for part in PARTS:count+=vcf_check(out/('joint-'+part)/'all.soft-filtered.vcf.gz',ids,{s['sample']:s[part] for s in p['samples'] if s[part]},rows[part],[s['sample'] for s in p['samples'] if not s[part]])
        need(check_final(merged,ids,p['samples'],rows)==count,'Final merged record count differs from partition total')
        artifacts=[pin(merged),pin(Path(str(merged)+'.tbi'))]+[pin(out/('joint-'+part)/name) for part in PARTS for name in ('joint.raw.vcf.gz','joint.raw.vcf.gz.tbi')]
        manifest=dict(schema_version=1,status='PASS_XY_GERMLINE_CALLING',processing_model=MODEL,clinical_sex_assigned=False,autosomes_or_mt_modified=False,annotation_performed=False,profile_sha256=profile_identity,samples=ids,partitions=rows,record_count=count,record_count_kind='DECOMPOSED_BIALLELIC_ROWS',decomposition='bcftools norm -m -any --multi-overlaps .; removed ALT maps to missing; original joint retained',artifacts=artifacts,producer_sha256=sha(__file__),filter_module=p['filter_module'],source_profile=pin(profile_copy),source_qc=p['qc_receipt'],alignment_contract=p['alignment_contract'],missing_Y_samples='NOT_ASSESSED_NO_GT_NEVER_HOMREF')
        atomic(out/'manifest.json',manifest);atomic(out/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='xy-germline-v1',artifacts=artifacts+[pin(out/'manifest.json')],identity=profile_identity));r.status('PASS_XY_GERMLINE_CALLING',record_count=count)
        return manifest
    except BaseException as e:r.status('FAILED_PARTIALS_RETAINED',error_type=type(e).__name__);raise
    finally:signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old);guard.close()

def public_control(image,out):
    """Real CPU tools, synthetic independent reads; never opens cohort inputs."""
    need(re.fullmatch(r'[^\s@]+@sha256:[0-9a-f]{64}',image),'Immutable image required')
    out=Path(out).resolve();need(not out.exists(),'Control output must be new');out.mkdir(parents=True)
    rng=random.Random(739);seq=''.join(rng.choice('ACGT') for _ in range(2000));pos=500;refbase=seq[pos-1];alt=next(x for x in 'ACGT' if x!=refbase)
    ref=out/'control.fa';ref.write_text('>chrX\n'+seq+'\n')
    sam=shutil.which('samtools');bc=shutil.which('bcftools');need(sam and bc,'Control requires existing samtools and bcftools')
    p=dict(reference={'path':str(ref)},samples=[],gatk_image=image,min_free_gib=1)
    r=Runner(p,out,min(END,now()+dt.timedelta(minutes=5)));r.command([sam,'faidx',str(ref)],'control-faidx');r.command([sam,'dict','-o',str(out/'control.dict'),str(ref)],'control-dict')
    gvcfs=[]
    for sid,ploidy in [('PUBLIC_H1',1),('PUBLIC_D2',2)]:
        path=out/(sid+'.sam')
        with path.open('x') as f:
            f.write('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrX\tLN:2000\n@RG\tID:'+sid+'\tSM:'+sid+'\tPL:ILLUMINA\n')
            for i in range(100):
                left=400+i%40;bases=seq[left-1:left+149]
                if ploidy==1 or i%2==0:bases=bases[:pos-left]+alt+bases[pos-left+1:]
                f.write(f'{sid}_{i}\t0\tchrX\t{left}\t60\t150M\t*\t0\t0\t{bases}\t'+('I'*150)+'\tRG:Z:'+sid+'\n')
        bam=out/(sid+'.bam');r.command([sam,'sort','-o',str(bam),str(path)],sid+'-sort');r.command([sam,'index',str(bam)],sid+'-index')
        g=out/(sid+'.g.vcf.gz');r.gatk(['HaplotypeCaller','-R',str(ref),'-I',str(bam),'-L','chrX:300-700','--sample-ploidy',str(ploidy),'--native-pair-hmm-threads','2','-ERC','GVCF','-O',str(g)],sid+'-hc',2,8,6);gvcfs.append(g)
    combined=out/'combined.g.vcf.gz';joint=out/'joint.vcf.gz'
    r.gatk(['CombineGVCFs','-R',str(ref),'-V',str(gvcfs[0]),'-V',str(gvcfs[1]),'-O',str(combined)],'control-combine',2,8,6)
    r.gatk(['GenotypeGVCFs','-R',str(ref),'-V',str(combined),'--standard-min-confidence-threshold-for-calling','0','-O',str(joint)],'control-joint',2,8,6)
    # GATK canonical sample order is lexical, independent of Combine input order.
    actual=r.command([bc,'query','-l',str(joint)],'control-samples',True).splitlines()
    vcf_check(joint,actual,{'PUBLIC_H1':1,'PUBLIC_D2':2},[('chrX',1,2000)])
    found=False
    with gzip.open(joint,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            c=line.rstrip().split('\t')
            if int(c[1])!=pos or c[3]!=refbase or alt not in c[4].split(','):continue
            gtindex=c[8].split(':').index('GT');genotypes={sid:value.split(':')[gtindex] for sid,value in zip(actual,c[9:])}
            hap=re.split('[/|]',genotypes['PUBLIC_H1']);dip=re.split('[/|]',genotypes['PUBLIC_D2'])
            need(len(hap)==1 and hap[0].isdigit() and int(hap[0])>0 and len(dip)==2 and set(dip)=={'0','1'},'Control genotypes do not match planted haploid ALT/diploid HET')
            found=True
    need(found,'Planted control allele missing')
    # Exercise exact installed bcftools projection: removed alternate is missing, never REF.
    multi=out/'multiallelic.vcf';split=out/'multiallelic.split.vcf'
    alts=[x for x in 'ACGT' if x!=refbase][:2]
    multi.write_text('##fileformat=VCFv4.2\n##contig=<ID=chrX,length=2000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPUBLIC_H1\tPUBLIC_D2\n'+f'chrX\t{pos}\t.\t{refbase}\t'+','.join(alts)+'\t50\tPASS\t.\tGT:AD:PL\t2:0,0,30:80,80,0\t1/2:0,15,15:90,80,70,60,0,50\n')
    r.command([bc,'norm','-m','-any','--multi-overlaps','.','-f',str(ref),'-c','e','-Ov','-o',str(split),str(multi)],'control-decompose')
    splitrows=[line.split('\t') for line in split.read_text().splitlines() if not line.startswith('#')]
    need(len(splitrows)==2 and [row[9].split(':')[0] for row in splitrows]==['.','1'] and [row[10].split(':')[0] for row in splitrows]==['1/.','./1'],'Installed decomposer invented reference or changed projection')
    proof=dict(schema_version=1,status='PASS_XY_MIXED_PLOIDY_CONTROL',synthetic=True,gatk_image=image,executor_sha256=sha(__file__),haploid_gt_verified=True,diploid_gt_verified=True,multiallelic_removed_alt_missing_verified=True,bcftools=pin(bc),samtools=pin(sam),public_control_source='DETERMINISTIC_SYNTHETIC_REFERENCE_AND_READS_NO_PATIENT_INPUT',joint=pin(joint),created_utc=now().isoformat())
    atomic(out/'public-control.json',proof);return proof

def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--profile',type=Path);a.add_argument('--profile-sha256');a.add_argument('--execute',action='store_true');a.add_argument('--public-control',action='store_true');a.add_argument('--gatk-image');a.add_argument('--control-output',type=Path);v=a.parse_args()
    if v.public_control:
        need(v.gatk_image and v.control_output and v.execute and not v.profile,'Explicit isolated public control inputs required');print(json.dumps(public_control(v.gatk_image,v.control_output)));return
    need(v.profile and v.profile_sha256 and not v.gatk_image and not v.control_output,'Exact reviewed profile required');need(sha(v.profile)==v.profile_sha256,'Exact reviewed profile SHA required');p=load(v.profile);end=checked_profile(p)
    if not v.execute:print(json.dumps(dict(status='VALIDATED_NOT_EXECUTED',calls=10,partitions=list(PARTS),max_cpu=24,max_memory_gib=120,deadline_utc=end.isoformat())));return
    execute(p,end)

if __name__=='__main__':main()
