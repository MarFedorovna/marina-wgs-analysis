#!/usr/bin/env python3
"""Separate full-transcript phenotype research reranking; default validates metadata only."""
import argparse
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import csv,datetime as dt,gzip,hashlib,html,itertools,json,multiprocessing,os
from pathlib import Path
import re,shutil,sys,time,uuid

HPO='HP:0001639';MONDO='MONDO:0005045'
STRENGTH={'Definitive':3,'Strong':2,'Moderate':1}
OLD_RANK=['clinvar_review_rank','consequence_rank','population_rank','predictors_damaging','AlphaMissense_research_rank']
NEW=['hcm_hpo_match','hcm_clingen_evidence_json','hcm_clingen_validity_rank','hcm_phenotype_priority','hcm_priority_reason','hcm_manual_review_flags_json','hcm_inheritance_mechanism_status']
REQUIRED=['sample','chrom','pos','ref','alt','gene','hgnc_id','transcript','candidate_eligible','quality_rank','sample_state','consequence','variant_id','mane_select',*OLD_RANK]
csv.field_size_limit(64*1024**2)
def need(ok,msg):
 if not ok:raise ValueError(msg)
def js(value):return json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':'))
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def public_sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()
def stamp(p):
 s=Path(p).stat();return {'bytes':s.st_size,'mtime_ns':s.st_mtime_ns}
def atomic(p,value):
 p=Path(p);tmp=p.with_name(p.name+'.tmp-'+uuid.uuid4().hex)
 with tmp.open('x',encoding='utf-8') as f:json.dump(value,f,ensure_ascii=False,indent=2);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def text(p):return gzip.open(p,'rt',encoding='utf-8-sig',newline='') if str(p).endswith('.gz') else Path(p).open(encoding='utf-8-sig',newline='')
def public_table(path,groups):
 with text(path) as stream:
  for n,line in enumerate(stream):
   need(n<100,'Public resource header not found')
   cleaned=line.lstrip('#')
   found=None
   for sep in ['\t',',']:
    fields=next(csv.reader([cleaned],delimiter=sep));normalized={x.strip().lower() for x in fields}
    if all(normalized&set(group) for group in groups):found=sep;break
   if found:
    for row in csv.DictReader(itertools.chain([cleaned],stream),delimiter=found):yield {str(k).strip().lower():v.strip() for k,v in row.items() if k is not None}
    return
  raise ValueError('Public resource missing required header')
def get(row,*keys):return next((row[k] for k in keys if row.get(k)), '')
def hcm_disease(row):
 code=get(row,'disease id (mondo)','disease_id','mondo_id');label=' '.join(get(row,'disease label','disease_name','disease').lower().split())
 return code==MONDO or re.search(r'(^| )hypertrophic cardiomyopathy( |$)',label) is not None
def resources(hpo,clingen):
 hp=set()
 for row in public_table(hpo,[['gene_symbol','gene'],['hpo_id']]):
  if row.get('hpo_id')==HPO:hp.add(get(row,'gene_symbol','gene').upper())
 cg=defaultdict(list)
 for row in public_table(clingen,[['gene symbol','gene_symbol','gene'],['classification','classification status'],['disease label','disease_name','disease']]):
  if hcm_disease(row):cg[get(row,'gene symbol','gene_symbol','gene').upper()].append(row)
 return hp,dict(cg)
def context(row,hp,cg,confirmed):
 gene=row['gene'].upper();hgnc=row.get('hgnc_id','');matches=[]
 for evidence in cg.get(gene,[]):
  eid=get(evidence,'gene id (hgnc)','hgnc_id')
  if hgnc and eid and eid!=hgnc:continue
  matches.append(evidence)
 rank=max([STRENGTH.get(get(e,'classification','classification status'),0) for e in matches]+[0]);hpo=gene in hp
 positive=bool(confirmed and (hpo or rank>0));reasons=[]
 if not confirmed:reasons.append('PHENOTYPE_NOT_CONFIRMED_NO_BOOST')
 if hpo:reasons.append('EXACT_HPO_HP0001639_GENE_ASSOCIATION_NOT_CAUSAL_PROOF')
 if rank:reasons.append('POSITIVE_CLINGEN_HCM_GENE_DISEASE_CONTEXT_NOT_GENOTYPE_AWARE')
 elif matches:reasons.append('WEAK_OR_NEGATIVE_CLINGEN_HCM_CURATION_NO_POSITIVE_BOOST')
 if not hpo and not matches:reasons.append('NO_EXACT_HCM_CONTEXT_IN_THESE_SOURCES_NOT_BENIGN')
 flags=[];af=row.get('af_max_observed','');classification=row.get('clinvar_significance','')
 if af in ('','.'):flags.append('SOURCE_AF_MISSING_NOT_ZERO')
 if 'conflict' in classification.lower():flags.append('CLINVAR_CONFLICT_MANUAL_REVIEW')
 labels={x.strip().lower() for x in re.split(r'[/&|]',classification)}
 if labels&{'pathogenic','likely_pathogenic'} and af not in ('','.'):
  try:
   if float(af)>=.01:flags.append('EXTERNAL_PATHOGENIC_LABEL_WITH_HIGH_SOURCE_AF_MANUAL_REVIEW')
  except ValueError:flags.append('SOURCE_AF_NOT_NUMERIC_MANUAL_REVIEW')
 return dict(hcm_hpo_match=int(hpo),hcm_clingen_evidence_json=js(matches),hcm_clingen_validity_rank=rank,hcm_phenotype_priority=int(positive),hcm_priority_reason=';'.join(reasons),hcm_manual_review_flags_json=js(flags),hcm_inheritance_mechanism_status='NOT_ASSESSED; GENE_DISEASE_CONTEXT_ONLY_NOT_GENOTYPE_AWARE')
def rank_key(row):
 # Same original lexicographic order after removing old phenotype/target/gene-validity dimensions.
 return (-int(row['hcm_phenotype_priority']),-int(row['hcm_clingen_validity_rank']),*[-int(row.get(k) or 0) for k in OLD_RANK],-int(bool(row.get('mane_select'))),row['variant_id'],row.get('gene',''),row.get('transcript',''),row['sample'])
def allele(row):
 chrom=str(row['chrom']);chrom='chr'+chrom.removeprefix('chr');need(',' not in row['alt'],'Decomposed exact allele required');return (chrom,int(row['pos']),row['ref'],row['alt'])
def is_eligible(row):return row.get('candidate_eligible')=='1' and row.get('quality_rank')=='1' and row.get('sample_state')=='alternate'
def scan_sample(job):
 path=Path(job['path']);before=stamp(path);need(before==job['source_stat'],'Source master changed before reading');best={};terms=defaultdict(set);rows=eligible_tx=0;hpkeys=set();cgkeys=set();ctxkeys=set();output=Path(job['aux'])
 with text(path) as inp,output.open('x',encoding='utf-8',newline='\n') as aux:
  reader=csv.DictReader(inp,delimiter='\t');need(reader.fieldnames==job['fields'] and set(REQUIRED)<=set(reader.fieldnames),'Full master field order differs')
  for row in reader:
   need(time.monotonic()<job['end'],'Local bounded rerank deadline reached');rows+=1
   if rows%250000==0:atomic(job['progress'],dict(status='SCANNING_FULL_SAMPLE_MASTER',label=job['label'],source_rows=rows,eligible_transcript_rows=eligible_tx))
   need(row['sample']==job['sample'],'Master contains another sample')
   if row['candidate_eligible']!='1':continue
   need(is_eligible(row),'Source eligibility contradicts quality/state; no silent change');key=allele(row);extra=context(row,job['hp'],job['cg'],job['confirmed']);value=dict(row,**extra);eligible_tx+=1;aux.write(js(value)+'\n')
   if extra['hcm_hpo_match']:hpkeys.add(key)
   if extra['hcm_clingen_validity_rank']>0:cgkeys.add(key)
   if extra['hcm_phenotype_priority']:ctxkeys.add(key)
   terms[key].update(x for x in re.split(r'[&,;|]',row['consequence']) if x)
   if key not in best or rank_key(value)<rank_key(best[key]):best[key]=value
  aux.flush();os.fsync(aux.fileno())
 need(stamp(path)==before,'Master changed during full scan');need(rows==job['expected_rows'],'Full per-sample source row count differs');need(len(best)==job['expected_eligible'],'Per-sample unique eligible allele count differs')
 atomic(job['progress'],dict(status='PASS_FULL_SAMPLE_SCAN',label=job['label'],source_rows=rows,eligible_transcript_rows=eligible_tx,eligible_unique_alleles=len(best)))
 return dict(sample=job['sample'],label=job['label'],rows=rows,eligible_tx=eligible_tx,best=best,terms={k:sorted(v) for k,v in terms.items()},hpo=hpkeys,clingen=cgkeys,positive=ctxkeys)
def clinical(path,samples):
 p=load(path);need(p.get('schema_version')==1,'Clinical profile schema required')
 if isinstance(p.get('samples'),list):
  rows=p['samples'];need(len(rows)==len(samples) and len({r.get('sample_id') for r in rows})==len(rows) and {r.get('sample_id') for r in rows}==set(samples),'Exact four sample clinical records required')
  need(p.get('created_utc') and p.get('assertion_source'),'User confirmation provenance and created_utc required');t=dt.datetime.fromisoformat(p['created_utc']);need(t.tzinfo is not None,'Clinical provenance timestamp needs timezone')
  for row in rows:
   diagnoses=row.get('diagnoses',[]);features=row.get('phenotypic_features',[])
   matched=any(isinstance(d,dict) and d.get('confirmation')=='confirmed_by_user' and (d.get('mondo_id')==MONDO or re.search(r'(^| )hypertrophic cardiomyopathy( |$)',str(d.get('label','')).lower())) for d in diagnoses)
   present=any(isinstance(f,dict) and f.get('hpo_id')==HPO and f.get('status')=='present' and f.get('source')=='user_report' for f in features)
   need(matched and present,'Each sample needs explicit user-confirmed HCM and present user-reported HP0001639')
 elif isinstance(p.get('samples'),dict):
  need(set(p['samples'])==set(samples),'Exact private per-sample clinical profile required')
  for sample in samples:
   row=p['samples'][sample];need(row.get('confirmed_hcm') is True and HPO in row.get('hpo_terms',[]) and row.get('confirmation_source') and row.get('confirmed_utc'),'Explicit confirmed HCM and provenance required');t=dt.datetime.fromisoformat(row['confirmed_utc']);need(t.tzinfo is not None,'Clinical confirmation timestamp needs timezone')
 else:raise ValueError('Clinical samples list/dictionary required')
 return p
def source_matrix(path,samples):
 before=stamp(path);selected={}
 with text(path) as stream:
  for row in csv.DictReader(stream,delimiter='\t'):
   if int(row['eligible_sample_count'])<=0:continue
   key=allele(row);need(key not in selected,'Duplicate exact matrix allele');eligible=json.loads(row['eligible_samples_json']);qc=json.loads(row['quality_passing_carriers_json']);raw=json.loads(row['raw_observed_alt_samples_json']);need(len(set(eligible))==len(eligible)==int(row['eligible_sample_count']) and set(eligible)<=set(qc)<=set(raw)<=set(samples) and len(qc)==int(row['carrier_count']),'Source matrix membership inconsistent');selected[key]=row
 need(stamp(path)==before,'Source matrix changed during read');return selected
def belongs(scope,eligible,qc):
 if scope=='cohort':return bool(eligible)
 if scope=='shared_ge2':return len(qc)>=2 and bool(eligible)
 if scope=='all4':return len(qc)==4 and bool(eligible)
 return scope in eligible

def table(path,fields,rows):
 with gzip.open(path,'wt',encoding='utf-8',newline='') if str(path).endswith('.gz') else path.open('x',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def run(views,profiles,hpo,clingen,output,workers=4,execute=False,max_hours=8):
 views=Path(views).resolve(strict=True);m=load(views/'manifest.json');need(m.get('status')=='PASS' and m.get('synthetic') is False and m.get('report_input_kind')=='FULLY_ANNOTATED_VCF_PACKETS_V1','Closed primary VCF-first views required');samples=m['samples'];labels=m['sample_file_labels'];need(len(samples)==4 and len(set(samples))==4,'Four exact source samples required');cp=clinical(profiles,samples);need(1<=workers<=4 and 0<max_hours<=8,'Bounded four-worker local run required');fields=m['logical_full_fields'];need(not set(NEW)&set(fields) and 'rank' not in fields,'Source contains new namespace/rank; no overwrite')
 masters={s:views/'masters'/(labels[s]+'.master.tsv.gz') for s in samples}
 for sample,p in masters.items():need(p.is_file() and p.stat().st_size==m['artifacts'][p.relative_to(views).as_posix()]['bytes'],'Completed individual master missing/size differs')
 matrix=views/'tables/allele-membership-matrix.tsv.gz';need(matrix.stat().st_size==m['artifacts']['tables/allele-membership-matrix.tsv.gz']['bytes'],'Closed source matrix size mismatch')
 if not execute:return dict(status='PLAN_ONLY_NO_PAYLOAD_SCAN_OR_OUTPUT',samples=4,expected_sample_rows=m['sample_master_rows'],sample_ALT_master_rows=sum(m['sample_master_rows'].values()),source_cohort_all_sample_rows=m['master_transcript_rows'],expected_candidates={s:m['counts'][s]['eligible_alleles'] for s in samples},workers=workers)
 output=Path(output).resolve();need(not output.is_relative_to(views) and not views.is_relative_to(output) and not output.is_symlink(),'New separate output required');need(not output.exists() or (output.is_dir() and not any(output.iterdir())),'Output must be absent or empty; no overwrite');stage=output.with_name('.'+output.name+'-'+uuid.uuid4().hex);stage.mkdir(parents=True);(stage/'transcript-evidence').mkdir();(stage/'progress').mkdir();end=time.monotonic()+max_hours*3600
 source_stats={str(p):stamp(p) for p in [Path(profiles),Path(hpo),Path(clingen),views/'manifest.json',matrix,*masters.values()]};hp,cg=resources(hpo,clingen);matrix_rows=source_matrix(matrix,samples);need(len(matrix_rows)==m['counts']['cohort']['eligible_alleles'],'Source candidate matrix count differs');jobs=[]
 for s in samples:jobs.append(dict(path=str(masters[s]),source_stat=stamp(masters[s]),sample=s,label=labels[s],fields=fields,hp=hp,cg=cg,confirmed=True,end=end,expected_rows=m['sample_master_rows'][s],expected_eligible=m['counts'][s]['eligible_alleles'],aux=str(stage/'transcript-evidence'/(labels[s]+'.eligible-transcripts.jsonl')),progress=str(stage/'progress'/(labels[s]+'.json'))))
 results={}
 with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
  for future in as_completed([pool.submit(scan_sample,j) for j in jobs]):
   result=future.result();results[result['sample']]=result;atomic(stage/'progress/status.json',dict(status='SCANNING_FULL_TRANSCRIPTS',completed_samples=len(results),samples=4))
 need(sum(r['rows'] for r in results.values())==sum(m['sample_master_rows'].values()),'Scanned ALT source rows differ; cohort/reference rows are separate')
 scopes=['cohort','shared_ge2','all4',*samples];ranked={s:[] for s in scopes};newfields=list(dict.fromkeys(['rank',*fields,*NEW,'hcm_candidate_consequence_union_json','hcm_functional_membership_basis',*next(iter(matrix_rows.values())).keys(),'ACMG_class','ACMG_criteria_evidence_status']))
 for key,member in matrix_rows.items():
  eligible={s:results[s]['best'][key] for s in samples if key in results[s]['best']};need(set(eligible)==set(json.loads(member['eligible_samples_json'])),'Full transcript eligibility differs from source membership');qc=json.loads(member['quality_passing_carriers_json'])
  for scope in scopes:
   if not belongs(scope,eligible,qc):continue
   chosen=eligible[scope] if scope in samples else min(eligible.values(),key=rank_key);terms=set(results[scope]['terms'][key]) if scope in samples else {term for s in eligible for term in results[s]['terms'][key]};value=dict(chosen,**member);value.update(hcm_candidate_consequence_union_json=js(sorted(terms)),hcm_functional_membership_basis='UNION_SOURCE_ELIGIBLE_TRANSCRIPTS_FOR_EXACT_ALLELE',ACMG_class='NOT_ASSIGNED',ACMG_criteria_evidence_status='NOT_ASSESSED_MANUAL_EVIDENCE_REQUIRED');ranked[scope].append(value)
 for scope in scopes:need(len(ranked[scope])==m['counts'][scope]['eligible_alleles'],'Source scope candidate count changed')
 (stage/'tables').mkdir();(stage/'functional').mkdir();tops={};summary={};category_index=[]
 for scope,rows in ranked.items():
  rows.sort(key=rank_key);label=labels.get(scope,scope)
  for i,row in enumerate(rows,1):row['rank']=i
  table(stage/'tables'/(label+'.all-candidates.tsv.gz'),newfields,rows);positive=[dict(row,rank=i) for i,row in enumerate((x for x in rows if x['hcm_phenotype_priority']),1)];top=positive[:20];tops[label]=top;table(stage/'tables'/(label+'.hcm-top20.tsv'),newfields,top);(stage/'tables'/(label+'.hcm-top20.json')).write_text(json.dumps(top,ensure_ascii=False,indent=2),encoding='utf-8')
  summary[label]=dict(source_eligible_alleles=len(rows),hcm_context_positive_alleles=len(positive),hpo_matched_alleles=sum(bool(x['hcm_hpo_match']) for x in rows),clingen_positive_alleles=sum(x['hcm_clingen_validity_rank']>0 for x in rows),top_rows=len(top));categories=defaultdict(list)
  for row in rows:
   for term in json.loads(row['hcm_candidate_consequence_union_json']):categories[term].append(row)
  folder=stage/'functional'/label;folder.mkdir()
  for term,rr in sorted(categories.items()):
   name=re.sub('[^A-Za-z0-9_-]','_',term)[:80]+'-'+hashlib.sha256(term.encode()).hexdigest()[:8]+'.tsv.gz';table(folder/name,newfields,rr);category_index.append(dict(scope=label,consequence=term,candidate_unique_alleles=len(rr),path='functional/'+label+'/'+name))
 table(stage/'functional-category-index.tsv',['scope','consequence','candidate_unique_alleles','path'],category_index)
 for path,before in source_stats.items():need(stamp(path)==before,'An original source changed during reranking')
 total=dict(producer_code_sha256=public_sha(__file__),public_resource_sha256={'hpo_gene_associations':public_sha(hpo),'clingen_validity':public_sha(clingen)},status='PASS_FULL_TRANSCRIPT_HCM_RESEARCH_RERANK',source_cohort_all_sample_rows=m['master_transcript_rows'],scanned_sample_ALT_rows=sum(r['rows'] for r in results.values()),eligible_transcript_rows={labels[s]:results[s]['eligible_tx'] for s in samples},source_counts_preserved=True,scopes=summary,hpo_exact_term=HPO,hpo_resource_genes=len(hp),clingen_HCM_resource_genes=len(cg),clinical_classification='NOT_ASSIGNED',inheritance_mechanism='NOT_ASSESSED',no_PP4_PVS1_or_LoF_GoF_assignment=True,top20_positive_context_only=True,worker_count=workers,payload_hashes_computed=False,completed_utc=dt.datetime.now(dt.timezone.utc).isoformat())
 sources=[dict(path=p,**v,verification='SOURCE_SIZE_MTIME_UNCHANGED_NO_PAYLOAD_HASH') for p,v in source_stats.items()];atomic(stage/'workbook-data.json',dict(summary=total,scopes=tops,sources=sources));atomic(stage/'manifest.json',dict(schema_version=1,**total,sample_file_labels=labels,sources=sources,original_rank_fields_unchanged=True,all_eligible_transcript_vectors_retained=True));atomic(stage/'clinical-profile-provenance.json',cp)
 (stage/'index.html').write_text('<!doctype html><meta charset="utf-8"><h1>HCM research reranking</h1><p>Full eligible transcripts; exact phenotype context. Research ranking, no clinical classification.</p><a href="workbook-data.json">Workbook data</a> · <a href="functional-category-index.tsv">Candidate functional views</a>',encoding='utf-8');need(not output.exists() or not any(output.iterdir()),'Output appeared; no overwrite')
 if output.exists():output.rmdir()
 os.rename(stage,output);return total

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--views',type=Path,required=True);p.add_argument('--clinical-profiles',type=Path,required=True);p.add_argument('--hpo-genes',type=Path,required=True);p.add_argument('--clingen-validity',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4);p.add_argument('--max-hours',type=float,default=8);p.add_argument('--execute',action='store_true');a=p.parse_args();print(json.dumps(run(a.views,a.clinical_profiles,a.hpo_genes,a.clingen_validity,a.output,a.workers,a.execute,a.max_hours)))
