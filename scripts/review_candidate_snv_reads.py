#!/usr/bin/env python3
"""Indexed local SNV read support. Research metrics only; no calling or classification."""
import argparse,csv,datetime as dt,gzip,hashlib,json,os,time,uuid
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path

EXCLUDED_FLAGS={4:'UNMAPPED',256:'SECONDARY',512:'QCFAIL',1024:'DUPLICATE',2048:'SUPPLEMENTARY'}

def need(ok,msg):
 if not ok:raise ValueError(msg)
def stamp(path):
 s=Path(path).stat();return dict(bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def js(x):return json.dumps(x,ensure_ascii=False,sort_keys=True,separators=(',',':'))
def atomic(path,value):
 p=Path(path);tmp=p.with_name(p.name+'.tmp-'+uuid.uuid4().hex)
 with tmp.open('x',encoding='utf-8') as f:json.dump(value,f,ensure_ascii=False,indent=2);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def target_base(read,pos0):
 """Return (base, numeric quality, reason); stored BAM sequence already reference-oriented."""
 cigar=read.cigartuples
 if not cigar:return None,None,'MISSING_CIGAR'
 # Reject unsupported operations anywhere, even after the target.
 if any(op not in range(9) or n<=0 for op,n in cigar):return None,None,'UNSUPPORTED_CIGAR'
 seq=read.query_sequence
 if not seq or seq=='*':return None,None,'MISSING_SEQUENCE'
 if sum(n for op,n in cigar if op in (0,1,4,7,8))!=len(seq):return None,None,'MALFORMED_CIGAR_SEQUENCE'
 ref=int(read.reference_start);query=0
 for op,n in cigar:
  if op in (0,7,8):
   if ref<=pos0<ref+n:
    i=query+pos0-ref;quals=getattr(read,'query_qualities',None)
    if quals is None or len(quals)!=len(seq) or quals[i]==255:return None,None,'MISSING_BASE_QUALITY'
    return seq[i].upper(),int(quals[i]),None
   ref+=n;query+=n
  elif op in (1,4):query+=n
  elif op in (2,3):
   if ref<=pos0<ref+n:return None,None,'DELETION_AT_LOCUS' if op==2 else 'REF_SKIP_AT_LOCUS'
   ref+=n
  # hard clip and padding consume neither axis
 return None,None,'NO_ALIGNED_BASE_AT_LOCUS'

def summarize_reads(reads,target,mapq=20,baseq=20,read_cap=100000,end=None):
 ref,alt=target['ref'],target['alt'];counts=Counter();excluded=Counter();templates={};unknown_names=0
 for ordinal,read in enumerate(reads,1):
  need(ordinal<=read_cap,'Per-locus read cap exceeded; no complete metric may be promoted')
  if end is not None:need(time.monotonic()<end,'Read-review deadline exceeded')
  counts['fetched_reads']+=1;flag=int(read.flag);bad=[name for bit,name in EXCLUDED_FLAGS.items() if flag&bit]
  if bad:
   excluded.update(bad);counts['excluded_flag_reads']+=1;continue
  mq=int(read.mapping_quality)
  if mq==255:excluded['MAPQ_UNKNOWN_255']+=1;continue
  if mq<mapq:excluded['LOW_MAPQ']+=1;continue
  base,q,reason=target_base(read,int(target['pos'])-1)
  if reason:excluded[reason]+=1;continue
  if q<baseq:excluded['LOW_BASEQ']+=1;continue
  if base not in 'ACGT':excluded['AMBIGUOUS_BASE']+=1;continue
  kind='REF' if base==ref else 'ALT' if base==alt else 'OTHER';strand='reverse' if flag&16 else 'forward'
  counts['passing_reads']+=1;counts['reads_'+kind]+=1;counts['reads_'+kind+'_'+strand]+=1
  name=getattr(read,'query_name',None)
  if not name or name=='*':unknown_names+=1;identity=('UNKNOWN_UNCOLLAPSED',ordinal)
  else:
   try:rg=read.get_tag('RG')
   except KeyError:rg=None
   if not rg:
    counts['passing_reads_missing_RG_uncollapsed']+=1;identity=('MISSING_RG_UNCOLLAPSED',ordinal)
   else:identity=(rg,name)
  entry=templates.setdefault(identity,{'bases':set(),'strands':set(),'reads':0});entry['bases'].add(base);entry['strands'].add(strand);entry['reads']+=1
 for item in templates.values():
  counts['templates_total']+=1
  if item['reads']>1:counts['multi_read_templates']+=1
  if len(item['bases'])!=1:counts['ambiguous_disagreeing_templates']+=1;continue
  base=next(iter(item['bases']));kind='REF' if base==ref else 'ALT' if base==alt else 'OTHER'
  counts['templates_'+kind]+=1
  orientation='both' if len(item['strands'])==2 else next(iter(item['strands']))
  counts['templates_'+kind+'_'+orientation]+=1
 counts['passing_reads_unknown_name_uncollapsed']=unknown_names
 for name in ('fetched_reads','passing_reads','excluded_flag_reads','templates_total','multi_read_templates','ambiguous_disagreeing_templates','passing_reads_missing_RG_uncollapsed'):counts.setdefault(name,0)
 for prefix in ('reads','templates'):
  for kind in ('REF','ALT','OTHER'):
   for suffix in ('','_forward','_reverse')+ (('_both',) if prefix=='templates' else ()):counts.setdefault(prefix+'_'+kind+suffix,0)
 flags=[]
 if counts['passing_reads']==0:flags.append('NO_PASSING_BASE_COVERAGE_NOT_REFERENCE_GENOTYPE')
 if counts['templates_ALT']<3:flags.append('ALT_TEMPLATES_LT3_DESCRIPTIVE')
 if counts['reads_ALT'] and not(counts['reads_ALT_forward'] and counts['reads_ALT_reverse']):flags.append('ALT_READS_ONE_STRAND_DESCRIPTIVE')
 return dict(status='ASSESSED_SNV_READ_SUPPORT',counts=dict(counts),exclusion_reasons=dict(excluded),descriptive_flags=flags)

def candidates(path,max_loci):
 op=gzip.open if str(path).endswith('.gz') else open;targets=[];seen=set()
 with op(path,'rt',encoding='utf-8-sig',newline='') as f:
  reader=csv.DictReader(f,delimiter='\t');need({'chrom','pos','ref','alt','hcm_phenotype_priority','sample_calls_json'}<=set(reader.fieldnames or []),'Candidate columns missing')
  for row in reader:
   if row['hcm_phenotype_priority']!='1':continue
   pos=int(row['pos']);need(pos>0,'Invalid position');key=(row['chrom'],pos,row['ref'],row['alt']);need(key not in seen,'Duplicate exact candidate allele');seen.add(key)
   need(len(targets)<max_loci,'Maximum target count exceeded');calls=json.loads(row['sample_calls_json']);need(isinstance(calls,dict),'Source sample calls must be a map')
   targets.append(dict(chrom=key[0],pos=pos,ref=key[2],alt=key[3],variant_id=row.get('variant_id',''),gene=row.get('gene',''),source_sample_calls=calls,source_membership={k:row.get(k,'') for k in ('quality_passing_carriers_json','eligible_samples_json','raw_observed_alt_samples_json')}))
 need(bool(targets),'No HCM-positive candidate targets');return targets

def header_identity(bam,sample):
 raw=bam.header if isinstance(bam.header,dict) else getattr(bam.header,'SAMheader',None)
 if raw is None:return dict(status='HEADER_SM_UNAVAILABLE_EXPLICIT_MAPPING_ONLY',read_group_count=0)
 if isinstance(raw,dict):
  groups=raw.get('RG',raw.get('@RG',[]));groups=[groups] if isinstance(groups,dict) else groups
 elif isinstance(raw,(str,bytes)):
  if isinstance(raw,bytes):raw=raw.decode('utf-8')
  groups=[dict(cell.split(':',1) for cell in line.split('\t')[1:] if ':' in cell) for line in raw.splitlines() if line.startswith('@RG\t')]
 else:raise ValueError('Unrecognized BAM SAM header representation')
 need(isinstance(groups,list) and all(isinstance(g,dict) for g in groups),'Invalid header RG structure')
 sms={g['SM'] for g in groups if g.get('SM')}
 need(not sms or sms=={sample},'BAM header SM differs from explicit source sample')
 return dict(status='EXACT_HEADER_SM_MATCH' if sms else 'HEADER_SM_UNAVAILABLE_EXPLICIT_MAPPING_ONLY',read_group_count=len(groups))

def worker(sample,entry,targets,mapq,baseq,cap,seconds):
 import bamnostic
 need(bamnostic.__version__=='1.3','Reviewed bamnostic 1.3 required')
 before={k:stamp(entry[k]) for k in ('bam','index')};out=[];end=time.monotonic()+seconds
 with bamnostic.AlignmentFile(entry['bam'],'rb',index_filename=entry['index'],require_index=True,max_cache=8) as bam:
  # Exact contig identities only: no chr aliases or patient identity inference.
  refs=dict(zip(bam.references,bam.lengths));identity=header_identity(bam,sample)
  for target in targets:
   need(time.monotonic()<end,'Read-review deadline exceeded');base={k:target[k] for k in ('chrom','pos','ref','alt','variant_id','gene')};base.update(sample=sample,source_call=target['source_sample_calls'][sample],source_membership=target['source_membership'])
   if len(target['ref'])!=1 or len(target['alt'])!=1 or target['ref'] not in 'ACGT' or target['alt'] not in 'ACGT' or target['ref']==target['alt']:result=dict(status='NOT_ASSESSED_UNSUPPORTED_NON_SNV',counts={},exclusion_reasons={},descriptive_flags=['NOT_A_REFERENCE_OR_BENIGN_INFERENCE'])
   else:
    need(target['chrom'] in refs and target['pos']<=refs[target['chrom']],'Target contig/position absent from exact BAM dictionary')
    result=summarize_reads(bam.fetch(target['chrom'],target['pos']-1,target['pos']),target,mapq,baseq,cap,end)
   out.append(dict(base,**result))
 need(before=={k:stamp(entry[k]) for k in before},'BAM/index changed during indexed review')
 return dict(sample=sample,rows=out,bam_index_stat=before,header_identity=identity)

def write_table(path,rows):
 fields=['sample','chrom','pos','ref','alt','variant_id','gene','status','source_call','source_membership','counts','exclusion_reasons','descriptive_flags']
 with Path(path).open('x',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
  for row in rows:w.writerow({k:js(row[k]) if isinstance(row[k],(dict,list)) else row[k] for k in fields})

def run(args):
 need(1<=args.workers<=4 and 1<=args.max_loci<=500 and 1<=args.read_cap<=100000 and 0<=args.mapq<=254 and 0<=args.baseq<=254 and 0<args.max_hours<=4,'Execution bounds invalid')
 before={str(Path(p).resolve()):stamp(p) for p in (args.candidates,args.bam_manifest)}
 targets=candidates(args.candidates,args.max_loci);manifest=json.loads(Path(args.bam_manifest).read_text(encoding='utf-8-sig'));need(manifest.get('schema_version')==1,'Unknown BAM manifest schema');samples=manifest.get('samples');need(isinstance(samples,dict) and len(samples)==4,'Exactly four explicitly mapped samples required')
 for target in targets:need(set(target['source_sample_calls'])==set(samples),'Candidate source samples differ from BAM mapping')
 need(before=={p:stamp(p) for p in before},'Candidate/manifest changed while loading')
 for sample,entry in samples.items():
  need(isinstance(sample,str) and sample and set(entry)>={'bam','index'},'Invalid explicit BAM/index mapping')
  for k in ('bam','index'):need(Path(entry[k]).is_file(),'Explicit BAM/index missing')
 output=Path(args.output).resolve();need(not output.exists(),'Output must be new; no overwrite/resume of partial review')
 plan=dict(status='PLAN_ONLY_NO_BAM_OPEN',targets=len(targets),samples=4,indexed_queries_max=4*len(targets),workers=args.workers)
 if not args.execute:return plan
 output.mkdir(parents=True);atomic(output/'status.json',dict(status='RUNNING_INDEXED_REVIEW',**{k:v for k,v in plan.items() if k!='status'}));results={}
 try:
  with ProcessPoolExecutor(max_workers=args.workers) as pool:
   jobs=[pool.submit(worker,s,e,targets,args.mapq,args.baseq,args.read_cap,args.max_hours*3600) for s,e in samples.items()]
   for future in as_completed(jobs):
    result=future.result();results[result['sample']]=result;atomic(output/'status.json',dict(status='RUNNING_INDEXED_REVIEW',completed_samples=len(results),samples=4))
  need(before=={p:stamp(p) for p in before},'Candidate/manifest changed during review');rows=[]
  for i,sample in enumerate(samples,1):write_table(output/('sample_%04d.read-support.tsv'%i),results[sample]['rows']);rows.extend(results[sample]['rows'])
  write_table(output/'cohort.read-support.tsv',rows)
  with (output/'targets.igv.bed').open('x',encoding='utf-8') as f:
   for i,t in enumerate(targets,1):f.write('%s\t%d\t%d\ttarget_%04d\n'%(t['chrom'],t['pos']-1,t['pos'],i))
  summary=dict(schema_version=1,status='PASS_LOCAL_INDEXED_SNV_READ_SUPPORT_REVIEW',created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),producer_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_stats=before,bam_index_stats={s:results[s]['bam_index_stat'] for s in samples},header_identity={s:results[s]['header_identity'] for s in samples},parameters=vars(args),targets=len(targets),rows=len(rows),samples=list(samples),results=rows,limitations=['RESEARCH_READ_AND_TEMPLATE_COUNTS_NOT_GATK_DP','NO_BAQ_NO_CALLING_NO_CLASSIFICATION','NO_QNAME_OR_SEQUENCE_OUTPUT','NO_GENOTYPE_INFERENCE_FROM_MISSING_COVERAGE','BAM_MAPPING_EXPLICIT_NOT_BIOLOGICAL_IDENTITY_PROOF'])
  atomic(output/'summary.json',summary);atomic(output/'status.json',dict(status=summary['status'],targets=len(targets),rows=len(rows)));return {k:summary[k] for k in ('status','targets','rows')}
 except Exception as exc:atomic(output/'status.json',dict(status='FAILED_NO_COMPLETE_REVIEW',error_type=type(exc).__name__));raise

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('candidates','bam-manifest','output'):p.add_argument('--'+name,required=True)
 p.add_argument('--workers',type=int,default=4);p.add_argument('--mapq',type=int,default=20);p.add_argument('--baseq',type=int,default=20);p.add_argument('--max-loci',type=int,default=500);p.add_argument('--read-cap',type=int,default=100000);p.add_argument('--max-hours',type=float,default=2);p.add_argument('--execute',action='store_true');print(js(run(p.parse_args())))
if __name__=='__main__':main()
