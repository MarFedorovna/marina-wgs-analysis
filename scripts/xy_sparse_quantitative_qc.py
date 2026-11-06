#!/usr/bin/env python3
"""Indexed, sparse X/Y technical measurements; never assigns sex or ploidy."""
import argparse
import bisect
import collections
import concurrent.futures
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import statistics
import subprocess
import time

PAR={'chrX':[(10001,2781479),(155701383,156030895)],'chrY':[(10001,2781479),(56887903,57217415)]}
FLAGS='0xF04'
PAR_SOURCE_SHA256='cc1b9b6dac63774e8e22ec87cfb71edfea59483534520094a1758a4a509f7311'


def need(ok,message):
    if not ok:raise ValueError(message)


def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return h.hexdigest()


def pin(p):
    p=Path(p).resolve(strict=True);return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))


def verify(p):
    path=Path(p['path']);need(not path.is_symlink() and path.stat().st_size==p['bytes'] and sha(path)==p['sha256'],'Small source/config pin changed.')
    return path


def stat_identity(p):
    p=Path(p).resolve(strict=True);s=p.stat()
    return dict(path=str(p),bytes=s.st_size,dev=s.st_dev,ino=s.st_ino,mtime_ns=s.st_mtime_ns)


def save(p,value):
    p=Path(p);temp=p.with_suffix(p.suffix+'.tmp')
    temp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');os.replace(temp,p)


def in_par(chrom,pos):return any(a<=pos<=b for a,b in PAR.get(chrom,[]))


class Fasta:
    def __init__(self,path):
        self.path=Path(path);self.file=self.path.open('rb')
        self.rows={p[0]:tuple(map(int,p[1:5])) for p in (line.rstrip().split('\t') for line in Path(str(path)+'.fai').read_text().splitlines())}
    def get(self,chrom,start,end):
        length,offset,bases,width=self.rows[chrom];need(0<=start<end<=length,'Out-of-bounds reference request.')
        address=lambda n:offset+(n//bases)*width+n%bases
        self.file.seek(address(start));value=self.file.read(address(end-1)-address(start)+1).replace(b'\n',b'').replace(b'\r',b'').decode().upper()
        need(len(value)==end-start,'FASTA/FAI mismatch.');return value
    def close(self):self.file.close()


def panel_rows(path,reference):
    rows=[];seen=set()
    with Path(path).open(newline='') as stream:
        for r in csv.DictReader(stream,delimiter='\t'):
            chrom,pos,ref,alt,af=r['chrom'],int(r['pos']),r['ref'],r['alt'],float(r['af'])
            need(chrom=='chrX' and not in_par(chrom,pos) and (chrom,pos) not in seen,'Panel must contain unique X-nonPAR sites.')
            need(ref in 'ACGT' and alt in 'ACGT' and len(ref)==len(alt)==1 and ref!=alt and .1<=af<=.9,'Common biallelic SNV AF contract changed.')
            need(reference.get(chrom,pos-1,pos)==ref,'Public SNP REF differs from exact analysis FASTA.')
            seen.add((chrom,pos));rows.append(dict(chrom=chrom,pos=pos,ref=ref,alt=alt,af=af,source_id=r.get('source_id','')))
    need(rows,'No public SNP sites.');return sorted(rows,key=lambda r:r['pos'])


def survey_windows(reference,sex_windows=96,autosome_windows=4,width=1000):
    result=[]
    for chrom in ['chr'+str(n) for n in range(1,23)]+['chrX','chrY']:
        length=reference.rows[chrom][0];count=sex_windows if chrom in PAR else autosome_windows
        for slot in range(count):
            anchor=int((slot+.5)*length/count);chosen=None
            for delta in range(0,100001,width):
                for start in (anchor+delta,anchor-delta):
                    end=start+width
                    if start<0 or end>length or any(start<b and end>a-1 for a,b in PAR.get(chrom,[])):continue
                    if any(r['chrom']==chrom and start<r['end'] and end>r['start'] for r in result):continue
                    sequence=reference.get(chrom,start,end)
                    if set(sequence)<=set('ACGT'):
                        chosen=dict(chrom=chrom,start=start,end=end,bases=width,gc_fraction=(sequence.count('G')+sequence.count('C'))/width,slot=slot);break
                if chosen:break
            if chosen:result.append(chosen)
    need(any(r['chrom']=='chrY' for r in result),'No ACGT-only Y nonPAR survey windows.')
    return result


def merge_regions(windows,sites,order):
    regions=collections.defaultdict(list)
    for r in windows:regions[r['chrom']].append((r['start'],r['end']))
    for r in sites:regions[r['chrom']].append((r['pos']-1,r['pos']))
    merged=[]
    for chrom in order:
        for a,b in sorted(regions[chrom]):
            if merged and merged[-1][0]==chrom and a<=merged[-1][2]:merged[-1]=(chrom,merged[-1][1],max(b,merged[-1][2]))
            else:merged.append((chrom,a,b))
    return merged


def parse_bases(ref,text,qualities,depth):
    if depth==0:
        need(text=='*' and qualities=='*','Unexpected zero-depth pileup encoding.');return collections.Counter()
    counts=collections.Counter();i=0;observations=0
    while i<len(text):
        c=text[i]
        if c=='^':need(i+1<len(text),'Truncated read-start marker.');i+=2;continue
        if c=='$':i+=1;continue
        if c in '+-':
            j=i+1
            while j<len(text) and text[j].isdigit():j+=1
            need(j>i+1,'Indel length absent.');n=int(text[i+1:j]);need(n>0 and j+n<=len(text),'Truncated indel.');i=j+n;continue
        need(c in '.,ACGTNacgtn*#<>','Unknown pileup base marker.')
        counts[ref if c in '.,' else c.upper()]+=1;observations+=1;i+=1
    need(observations==depth==len(qualities),'Pileup depth/quality/base counts differ.')
    need(all(ord(q)-33>=20 for q in qualities),'Unexpected base below reviewed BQ threshold.')
    return counts


def classify(counts,ref,alt):
    dp=sum(counts[b] for b in 'ACGT');paired=counts[ref]+counts[alt];other=sum(counts.values())-paired
    result=dict(depth_acgt=dp,ref_count=counts[ref],alt_count=counts[alt],other_count=other,allele_balance=counts[alt]/paired if paired else None)
    if dp<10:state='LOW_DP'
    elif dp>100:state='HIGH_DP'
    elif other/max(sum(counts.values()),1)>.05:state='OFF_ALLELE_OR_NONBASE'
    elif result['allele_balance'] is None:state='NO_REF_ALT'
    elif .25<=result['allele_balance']<=.75 and min(counts[ref],counts[alt])>=3:state='HET_BALANCED'
    elif result['allele_balance']<=.05:state='HOM_REF_LIKE'
    elif result['allele_balance']>=.95:state='HOM_ALT_LIKE'
    else:state='AMBIGUOUS_AB'
    return dict(result,status=state)


def wilson(k,n):
    if not n:return None
    z=1.959963984540054;p=k/n;den=1+z*z/n;center=(p+z*z/(2*n))/den;delta=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0,center-delta),min(1,center+delta)]


def heterozygosity(rows,expected_he_available=True):
    called=[r for r in rows if r['status'] in ('HET_BALANCED','HOM_REF_LIKE','HOM_ALT_LIKE')]
    het=sum(r['status']=='HET_BALANCED' for r in called);bins=len({r['pos']//10000000 for r in called})
    sufficient=len(called)>=100 and bins>=10
    expected=sum(2*r['af']*(1-r['af']) for r in called) if expected_he_available else None
    return dict(panel_sites=len(rows),status_counts=dict(collections.Counter(r['status'] for r in rows)),n_called=len(called),n_het=het,
        het_fraction=het/len(called) if called else None,het_fraction_site_wilson95=wilson(het,len(called)),covered_10Mb_bins=bins,
        sufficiency_status='DESCRIPTIVE_PANEL_METRICS_ONLY' if sufficient else 'INSUFFICIENT_CALLABLE_SITES_OR_SPREAD',
        sum_expected_heterozygosity=expected,F_descriptive=(1-het/expected) if sufficient and expected else None,
        F_status='F_PROXY_FROM_AB_FILTERED_CALLS_GLOBAL_AF_NOT_ANCESTRY_CALIBRATED' if sufficient and expected else 'NOT_ESTIMATED',
        sex_assigned=False,ploidy_assigned=False)


def deadline_guard(deadline,seconds):
    remaining=min(seconds,(dt.datetime.fromisoformat(deadline)-dt.datetime.now(dt.timezone.utc)).total_seconds());need(remaining>0,'QC deadline expired.')
    def expired(signum,frame):raise TimeoutError('Separate bounded XY QC deadline reached.')
    signal.signal(signal.SIGALRM,expired);signal.setitimer(signal.ITIMER_REAL,remaining)


def child_limit():
    import resource
    resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))


def measure(job):
    import resource
    profile,index=job;p=profile;sample=p['samples'][index];out=Path(p['output']);deadline_guard(p['deadline_utc'],p['max_wall_seconds'])
    os.sched_setaffinity(0,set(p['cpu_ids'][index*4:(index+1)*4]))
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    need(stat_identity(sample['bam']['path'])==sample['bam'],'BAM stat identity changed; no BAM hash is computed.')
    verify(sample['index']);verify(sample['receipt'])
    windows=p['windows'];sites=p['sites'];lookup={(r['chrom'],r['pos']):r for r in sites};observed={}
    by_chrom=collections.defaultdict(list)
    for i,r in enumerate(windows):by_chrom[r['chrom']].append((r['start'],r['end'],i))
    for r in by_chrom.values():r.sort()
    starts={chrom:[r[0] for r in ranges] for chrom,ranges in by_chrom.items()}
    stats=[dict(depth_sum=0,breadth1=0,breadth5=0,breadth10=0,max_depth=0) for _ in windows];capped=0
    prefix='sample_'+str(index+1).zfill(4);view=None;pileup=None;began=time.monotonic()
    with (out/(prefix+'.samtools-view.log')).open('xb') as v_err,(out/(prefix+'.samtools-pileup.log')).open('xb') as p_err:
        try:
            view=subprocess.Popen(['samtools','view','-u','-@','2','-M','-L',p['bed']['path'],'-q','30','-F',FLAGS,sample['bam']['path']],stdout=subprocess.PIPE,stderr=v_err,preexec_fn=child_limit)
            pileup=subprocess.Popen(['samtools','mpileup','-A','-q','30','-Q','20','--ff',FLAGS,'-d','10000','-f',p['reference']['path'],'-l',p['bed']['path'],'-'],stdin=view.stdout,stdout=subprocess.PIPE,stderr=p_err,text=True,encoding='utf-8',preexec_fn=child_limit)
            view.stdout.close();last=None
            for line in pileup.stdout:
                columns=line.rstrip('\n').split('\t');need(len(columns)==6,'Unexpected single-BAM pileup schema.')
                chrom,pos,ref,depth=columns[0],int(columns[1]),columns[2],int(columns[3]);key=(chrom,pos)
                need(key!=last,'Duplicate pileup position.');last=key
                counts=parse_bases(ref,columns[4],columns[5],depth);acgt=sum(counts[b] for b in 'ACGT')
                capped+=int(depth>=10000)
                if key in lookup:
                    need(ref==lookup[key]['ref'] and key not in observed,'SNP reference/position changed.')
                    observed[key]=dict(lookup[key],**classify(counts,ref,lookup[key]['alt']))
                n=bisect.bisect_right(starts.get(chrom,[]),pos-1)-1
                if n>=0 and pos<=by_chrom[chrom][n][1]:
                    s=stats[by_chrom[chrom][n][2]];s['depth_sum']+=acgt;s['max_depth']=max(s['max_depth'],acgt)
                    for cutoff in (1,5,10):s['breadth'+str(cutoff)]+=int(acgt>=cutoff)
            need(pileup.wait()==0 and view.wait()==0,'Indexed BAM selection or pileup failed.')
        finally:
            if pileup is not None and pileup.stdout is not None:pileup.stdout.close()
            for proc in (pileup,view):
                if proc is not None and proc.poll() is None:proc.terminate();proc.wait(timeout=10)
    rows=[observed.get((r['chrom'],r['pos']),dict(r,**classify(collections.Counter(),r['ref'],r['alt']))) for r in sites]
    with (out/(prefix+'.x-site-counts.tsv')).open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    survey=[]
    for w,s in zip(windows,stats):survey.append(dict(w,**s,mean_depth=s['depth_sum']/w['bases']))
    summary={}
    for chrom in by_chrom:
        a=[r for r in survey if r['chrom']==chrom];bases=sum(r['bases'] for r in a)
        summary[chrom]=dict(windows=len(a),target_ACGT_bases=bases,mean_depth=sum(r['depth_sum'] for r in a)/bases,
            median_window_depth=statistics.median(r['mean_depth'] for r in a),windows_with_signal=sum(r['breadth1']>0 for r in a),
            **{'breadth_ge'+str(c):sum(r['breadth'+str(c)] for r in a)/bases for c in (1,5,10)})
    auto=[summary[c]['median_window_depth'] for c in summary if c not in PAR];base=statistics.median(auto)
    ratios={c:(summary[c]['median_window_depth']/base if base else None) for c in PAR if c in summary}
    result=dict(sample_id=sample['sample_id'],status='PASS_DESCRIPTIVE_SPARSE_XY_QC',X_nonPAR=heterozygosity(rows,p['expected_he_available']),
        survey=summary,window_metrics=survey,autosome_median_chromosome_window_depth=base,
        autosome_nonzero_chromosomes=sum(v>0 for v in auto),autosome_max_min_nonzero_ratio=max(auto)/min(v for v in auto if v>0) if any(auto) else None,
        sex_chromosome_to_autosome_median_window_depth=ratios,seconds=time.monotonic()-began,sex_assigned=False,ploidy_changed=False)
    result['positions_at_pileup_depth_cap']=capped
    save(out/(prefix+'.summary.json'),result);return result


def prepare(args):
    cfg=load(args.config);reference=Fasta(cfg['reference']);manifest=load(args.panel_manifest)
    need(sha(cfg['par_source'])==PAR_SOURCE_SHA256,'Pinned official GRC PAR source changed.')
    need(len(cfg['samples'])==4 and reference.rows['chrX'][0]==156040895 and reference.rows['chrY'][0]==57227415,'Exact four-sample GRCh38 BAM QC required.')
    need(manifest.get('status')=='PASS_PUBLIC_X_NONPAR_COMMON_SNP_PANEL' and manifest.get('assembly')=='GRCh38','Reviewed independent public GRCh38 panel required.')
    sites=panel_rows(args.panel,reference);windows=survey_windows(reference);reference.close()
    args.profile.parent.mkdir(parents=True,exist_ok=True);bed=args.profile.with_suffix('.targets.bed')
    need(not args.profile.exists() and not bed.exists() and not args.output.exists(),'Use new private profile/BED/output paths.')
    order=list(reference.rows)
    with bed.open('x') as stream:
        for chrom,a,b in merge_regions(windows,sites,order):stream.write(f'{chrom}\t{a}\t{b}\n')
    samples=[]
    for row in cfg['samples']:
        bam=Path(row['bam']);receipt=load(row['receipt']);need(receipt.get('sample')==row['sample_id'] and receipt.get('completed_utc') and bam.name in receipt.get('outputs',{})
            and receipt.get('config',{}).get('reference',{}).get('sha256')==cfg['reference_sha256'],'Completed BAM receipt/reference mismatch.')
        samples.append(dict(sample_id=row['sample_id'],bam=stat_identity(bam),index=pin(str(bam)+'.bai'),receipt=pin(row['receipt']),historical_bam_sha256_not_recomputed=receipt['outputs'][bam.name]))
    available=sorted(os.sched_getaffinity(0));need(len(available)>=16,'Sixteen bounded CPUs required.')
    need(pin(args.panel)==manifest.get('sites') and manifest.get('row_count')==len(sites),'Public panel/source manifest mismatch.')
    verify(manifest['source_manifest'])
    he_available=manifest.get('allele_frequency_source')=='1000GENOMES_PHASE3_GLOBAL_AF'
    profile=dict(schema_version=1,status='CANDIDATE_SPARSE_XY_QC_NOT_EXECUTED',root_reviewed=False,producer=pin(__file__),
        config=pin(args.config),reference=stat_identity(cfg['reference']),reference_expected_sha256_not_recomputed=cfg.get('reference_sha256'),
        fai=pin(str(cfg['reference'])+'.fai'),par_source=pin(cfg['par_source']),public_panel=pin(args.panel),panel_manifest=pin(args.panel_manifest),
        expected_he_available=he_available,sites=sites,windows=windows,bed=pin(bed),samples=samples,
        output=str(args.output.resolve()),cpu_ids=available[-16:],workers=4,max_memory_bytes=24*1024**3,max_wall_seconds=2400,
        deadline_utc=dt.datetime.fromisoformat(args.deadline_utc).isoformat(),no_full_BAM_hash=True,no_full_BAM_scan=True,no_variant_calling=True,no_ploidy_change=True)
    save(args.profile,profile);print(json.dumps(dict(status=profile['status'],profile=pin(args.profile),panel_sites=len(sites),survey_windows=len(windows)),indent=2))


def run(args):
    import resource
    need(sha(args.profile)==args.reviewed_profile_sha256,'Exact reviewed profile SHA required.');p=load(args.profile)
    verify(p['producer']);need(p['producer']['sha256']==sha(__file__),'Producer changed.')
    if not args.execute:print(json.dumps(dict(status='PLAN_ONLY_NO_BAM_READ',samples=len(p['samples']),sites=len(p['sites']),windows=len(p['windows']))));return
    need(p.get('root_reviewed') is True and p.get('status')=='ROOT_REVIEWED_SPARSE_XY_QC','ROOT-reviewed exact profile required.')
    need(len(p['samples'])==4 and len(p['cpu_ids'])==16 and p['workers']==4 and p['max_memory_bytes']<=24*1024**3 and p['max_wall_seconds']<=2400,'Resource caps changed.')
    for key in ('config','fai','par_source','public_panel','panel_manifest','bed'):verify(p[key])
    need(stat_identity(p['reference']['path'])==p['reference'],'Reference stat identity changed.')
    os.sched_setaffinity(0,set(p['cpu_ids']));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));deadline_guard(p['deadline_utc'],p['max_wall_seconds'])
    out=Path(p['output']);out.mkdir(mode=0o700,parents=True,exist_ok=False);started=dt.datetime.now(dt.timezone.utc).isoformat()
    save(out/'status.json',dict(status='RUNNING_SPARSE_XY_QC',pid=os.getpid(),started_utc=started,deadline_utc=p['deadline_utc']))
    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(measure,[(p,i) for i in range(4)]))
        report=dict(schema_version=1,status='PASS_DESCRIPTIVE_SPARSE_XY_QC',started_utc=started,completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
            producer=p['producer'],reviewed_profile=pin(args.profile),samples=results,sex_assigned=False,ploidy_changed=False,variant_calling_started=False,
            method=dict(MQ=30,BQ_BAQ=20,excluded_flags=FLAGS,pair_overlap_removed=True,anomalous_pairs_included=True,BAQ='samtools default, reference provided',
                Y_scope='ACGT-only surveyed nonPAR windows; not complete Y breadth',unique_mappability_mask=False,F_statistic='descriptive public global AF proxy from AB-filtered classifications only',
                X_classification=dict(min_ACGT_depth=10,max_ACGT_depth=100,min_het_ref_alt_each=3,het_AB=[.25,.75],hom_ref_AB_max=.05,hom_alt_AB_min=.95,max_other_fraction=.05),
                minimum_F_sites=100,minimum_F_10Mb_bins=10,pileup_depth_cap=10000),
            limitations=['Sparse survey can miss covered targets and cannot establish genome-wide Y absence.','NonPAR and high mapping quality do not establish unique mappability.',
                'Target/GC bias, LD, ancestry mismatch, contamination and mosaicism can distort ratios and heterozygosity.',
                'Wilson interval is a site-count approximation; LD/ascertainment are not modeled. It is not diagnostic confidence.',
                'No sex, karyotype, copy number or calling ploidy is assigned. Insufficient callable sites remain unresolved.'])
        save(out/'quantitative-xy-profile.json',report);save(out/'status.json',dict(status=report['status'],completed_utc=report['completed_utc'],pid=os.getpid()))
        artifacts={f.name:pin(f) for f in out.iterdir() if f.is_file() and f.name!='workflow-receipt.json'}
        save(out/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='sparse_quantitative_xy_qc_only',artifacts=artifacts,no_ploidy_change=True,no_variant_calling=True))
    except BaseException as error:
        save(out/'status.json',dict(status='FAILED_OR_BOUNDED_STOP',error_type=type(error).__name__,pid=os.getpid()));raise
    print('PASS descriptive sparse X/Y QC; no sex or ploidy assignment')


def main():
    os.umask(0o077);parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    prep=sub.add_parser('prepare')
    for name in ('config','panel','panel-manifest','output','profile'):prep.add_argument('--'+name,type=Path,required=True)
    prep.add_argument('--deadline-utc',required=True)
    action=sub.add_parser('run');action.add_argument('--profile',type=Path,required=True);action.add_argument('--reviewed-profile-sha256',required=True);action.add_argument('--execute',action='store_true')
    args=parser.parse_args();{'prepare':prepare,'run':run}[args.command](args)


if __name__=='__main__':main()
