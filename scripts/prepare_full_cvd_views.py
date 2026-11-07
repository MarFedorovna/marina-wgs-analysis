import csv,gzip,json,sqlite3,re,collections
from pathlib import Path

import argparse
parser=argparse.ArgumentParser(description="Build full cardiovascular research views from private prepared inputs.")
parser.add_argument('--vcf-qc',type=Path,required=True,help='Autosomal QC-passing bgzip VCF with four samples')
parser.add_argument('--opencravat-sqlite',type=Path,required=True)
parser.add_argument('--clingen-validity-csv',type=Path,required=True)
parser.add_argument('--vep-tables',type=Path,required=True,help='Directory containing cohort functional-category TSV.GZ files')
parser.add_argument('--supplementary-json',type=Path,required=True,help='Private prepared JSON containing previous cohort, X/Y and MT sheets')
parser.add_argument('--xy-quality-vcf',type=Path,required=True)
parser.add_argument('--samples',nargs=4,required=True,metavar=('SAMPLE1','SAMPLE2','SAMPLE3','SAMPLE4'))
parser.add_argument('--exclusions',type=Path,help='Private JSON array of CHROM:POS:REF:ALT keys')
parser.add_argument('--output-json',type=Path,required=True,help='Private output; may contain sensitive variant and sample data')
args=parser.parse_args();SAMPLES=args.samples
source=json.load(args.supplementary_json.open(encoding='utf8'));old={item['name']:item for item in source['sheets']};oldh={v:i for i,v in enumerate(old['Все_исходные_кандидаты']['headers'])};oldrows={r[oldh['Источник варианта']]:r for r in old['Все_исходные_кандидаты']['values']}
excluded=set(json.load(args.exclusions.open(encoding='utf8'))) if args.exclusions else set()
cvfile=args.clingen_validity_csv
regex=re.compile(r'cardio|cardiac|heart|aort|vascular|arterial|artery|arrhyth|long.qt|short.qt|brugada|ventricul|atrial|hypercholesterolemia|lipoprotein|marfan|sudden cardiac',re.I);cv=collections.defaultdict(list)
with cvfile.open(encoding='utf8',newline='') as f:
 lines=list(csv.reader(f));head=next(i for i,r in enumerate(lines) if any('GENE SYMBOL' in x.upper() for x in r));
 for vals in lines[head+1:]:
  if len(vals)<10:continue
  row=dict(zip(lines[head],vals));gene=row.get('GENE SYMBOL','');disease=row.get('DISEASE LABEL','');level=row.get('CLASSIFICATION','')
  if gene and (regex.search(disease) or regex.search(row.get('GCEP',''))) and level.lower() in ('definitive','strong','moderate'):cv[gene].append((disease,level,row.get('DISEASE ID (MONDO)','')))
def norm(a):
 c,p,r,t=a.split(':',3);p=int(p);r='' if r=='-' else r;t='' if t=='-' else t
 while r and t and r[-1]==t[-1]:r=r[:-1];t=t[:-1]
 while r and t and r[0]==t[0]:r=r[1:];t=t[1:];p+=1
 return c.removeprefix('chr'),p,r or '-',t or '-'
def num(x):
 try:return float(x) if x not in ('','.','NA',None) else None
 except (ValueError,TypeError):return None
def pathogenic(x):
 sig=str(x or '').lower().replace('_',' ')
 if 'conflict' in sig:return False
 return bool({z.strip() for z in re.split(r'[/|;&]',sig)}&{'pathogenic','likely pathogenic'})
def qc(call):
 try:return call.get('FT')=='PASS' and float(call['DP'])>=10 and float(call['GQ'])>=20 and any(x.isdigit() and int(x)>0 for x in re.split(r'[/|]',call.get('GT','')))
 except (ValueError,TypeError,KeyError):return False
con=sqlite3.connect('file:'+str(args.opencravat_sqlite.resolve()).replace('\\','/')+'?mode=ro',uri=True);con.row_factory=sqlite3.Row
assert con.execute('pragma quick_check(1)').fetchone()[0]=='ok'
genes={x['base__hugo']:dict(x) for x in con.execute('select * from gene')}
records={};dups=set()
for rec in con.execute('select * from variant'):
 k=norm(':'.join(str(rec[v]) for v in ('base__chrom','base__pos','base__ref_base','base__alt_base')));r=dict(rec)
 if k in records and records[k]!=r:dups.add(k)
 else:records[k]=r
for k in dups:del records[k]
if not records: raise ValueError("OpenCRAVAT database has no unambiguous variant rows")
allowed={'MIS','SYN','STG','STL','FSD','FSI','IND','INI','SPL','MLO','UT3','UT5','NMD','EXL','STR','CSS'}
impact_order=['STG','FSD','FSI','SPL','STL','MLO','EXL','MIS','IND','INI','UT3','UT5','SYN','STR','NMD','CSS']
ontology=json.loads(con.execute("select subdict from variant_reportsub where module='base'").fetchone()[0])['so']
table=args.vep_tables;catnames=['stop_gained','stop_lost','start_lost','frameshift_variant','splice_acceptor_variant','splice_donor_variant','splice_region_variant','splice_donor_5th_base_variant','splice_donor_region_variant','splice_polypyrimidine_tract_variant','missense_variant','synonymous_variant','inframe_deletion','inframe_insertion','pLoF_candidates']
vep_categories=collections.defaultdict(set);vep_af={};vep_best={};vep_sigs=collections.defaultdict(set)
for cat in catnames:
 path=next(table.glob(cat+'-*.tsv.gz'))
 with gzip.open(path,'rt',encoding='utf8',newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r.get('scope_qc_member')=='1':
    allele=r['allele_id'];vep_categories[allele].add(cat);af=num(r.get('af_max_observed'))
    if af is not None:vep_af[allele]=max(vep_af.get(allele,0),af)
    if r.get('clinvar_significance'):vep_sigs[allele].add(r['clinvar_significance'])
    if af is None or af<=.01 or pathogenic(r.get('clinvar_significance')):
     fields=['gene','transcript','hgvs_c','hgvs_p','consequence','impact','mane_select','SIFT','PolyPhen','CADD_PHRED','AlphaMissense_score','AlphaMissense_source_transcript','FATHMM','clinvar_significance','clinvar_review','clinvar_conditions','clinvar_ids','af_max_observed']
     mini={k:r.get(k,'') for k in fields};prev=vep_best.get(allele)
     def vepkey(x):return (0 if x.get('gene') in cv else 1,{'HIGH':0,'MODERATE':1,'LOW':2}.get(x.get('impact'),3),not bool(x.get('mane_select')),x.get('gene',''),x.get('transcript',''))
     if prev is None or vepkey(mini)<vepkey(prev):vep_best[allele]=mini
 print('VEP подготовлен',cat,flush=True)
H=['Ранг','Образцы-носители','Ген','Хромосома','Позиция','REF','ALT','HGVS-c','HGVS-p','Категории по всем транскриптам','Consequence выбранного транскрипта','Транскрипт','Образец метрик','GT','DP','AD','Доля ALT в AD','GQ','QUAL','Фильтр сайта','Фильтр генотипа','Число носителей','Доля из 4','ClinVar','ClinVar review','ClinVar заболевание','ClinVar ID','gnomAD AF max','gnomAD4 AF','gnomAD4 AC','gnomAD4 AN','Отбор','ССЗ-контекст ClinGen','SIFT','SIFT prediction','SIFT transcript','PolyPhen','PolyPhen transcript','CADD PHRED','REVEL','REVEL transcript','AlphaMissense','AlphaMissense transcript','FATHMM','FATHMM transcript','SpliceAI DS_AG','SpliceAI DS_AL','SpliceAI DS_DG','SpliceAI DS_DL','SpliceAI DP_AG','SpliceAI DP_AL','SpliceAI DP_DG','SpliceAI DP_DL','dbSNP','ClinGen заболевание','ClinGen классификация','HPO','OMIM ID','GTEx gene','GTEx tissue','UniProt','NCBI Gene ID','PubMed count','LoFtool','pHI','Prec','Статус OpenCRAVAT',*[f'GT {sample}' for sample in SAMPLES],'Оценки и транскрипты источников','Источник варианта','Критерии рейтинга']
def clean(v):
 if v in ('','.','NA','N/A','[]','{}'):return None
 if isinstance(v,(dict,list)):return json.dumps(v,ensure_ascii=False,separators=(',',':'))
 return v
all_candidates={};missing=[];n=0;quality_eligible=0
with gzip.open(args.vcf_qc,'rt',encoding='utf8') as f:
 for line in f:
  if line.startswith('#CHROM'):
   if line.rstrip('\n').split('\t')[9:]!=SAMPLES:raise ValueError('VCF sample order differs from --samples')
   continue
  if line.startswith('#'):continue
  v=line.rstrip('\n').split('\t');n+=1;a=':'.join((v[0],v[1],v[3],v[4]));
  if a in excluded:continue
  calls={s:dict(zip(v[8].split(':'),x.split(':'))) for s,x in zip(SAMPLES,v[9:])};car=[s for s in SAMPLES if qc(calls[s])]
  if not car or v[6]!='PASS' or num(v[5]) is None or num(v[5])<30:continue
  quality_eligible+=1;r=records.get(norm(a))
  if not r:missing.append(a);continue
  try:mapping=json.loads(r.get('base__all_mappings') or '{}')
  except json.JSONDecodeError:mapping={}
  options=[]
  for gene,ts in mapping.items():
   for t in ts:
    if len(t)>=5:options.append((gene,t))
  codes={r.get('base__so')}|{t[2] for _,t in options}
  functional=bool(codes&allowed or vep_categories.get(a));prior_sig=' | '.join(sorted(vep_sigs.get(a,set())));sig=r.get('clinvar__sig') or prior_sig or None;p=pathogenic(sig)
  af=[num(r.get('gnomad4__af'))]+[num(r.get('gnomad4__af_'+i)) for i in ('afr','ami','amr','asj','eas','fin','mid','nfe','sas','rem')]
  previous=oldrows.get(a)
  if previous:af.append(num(previous[oldh['gnomAD AF max']]))
  af.append(vep_af.get(a))
  af=max((x for x in af if x is not None),default=None)
  if not p and not (functional and (af is None or af<=.01)):continue
  def optkey(x):
   gene,t=x;code=t[2];return (0 if gene in cv else 1,impact_order.index(code) if code in impact_order else 99,not bool(t[1]),gene,t[3])
  gene,t=min(options,key=optkey) if options else (r.get('base__hugo'),[None,None,r.get('base__so'),r.get('base__transcript'),r.get('base__cchange')]);
  vr=vep_best.get(a,{})
  if vr.get('gene') in cv and gene not in cv:gene=vr['gene'];t=[None,vr.get('hgvs_p'),None,vr.get('transcript'),vr.get('hgvs_c')]
  ctx=cv.get(gene,[]);chosen=SAMPLES[next(i for i,s in enumerate(SAMPLES) if s in car)];g=genes.get(gene,{})
  reason='ClinVar P/LP; высокая AF' if p and af is not None and af>.01 else 'ClinVar P/LP' if p else 'Функциональный; AF отсутствует' if af is None else 'Функциональный; AF ≤ 1%'
  categories={ontology.get(x,x) for x in codes if x and x in allowed}|vep_categories.get(a,set());
  z={'Образцы-носители':'; '.join(car),'Ген':gene,'Хромосома':v[0],'Позиция':int(v[1]),'REF':v[3],'ALT':v[4],'HGVS-c':t[4],'HGVS-p':t[1],'Категории по всем транскриптам':'; '.join(sorted(categories)),'Consequence выбранного транскрипта':ontology.get(t[2],t[2]),'Транскрипт':t[3],'Образец метрик':chosen,'GT':calls[chosen].get('GT'),'DP':num(calls[chosen].get('DP')),'AD':calls[chosen].get('AD'),'GQ':num(calls[chosen].get('GQ')),'QUAL':num(v[5]),'Фильтр сайта':v[6],'Фильтр генотипа':calls[chosen].get('FT'),'Число носителей':len(car),'Доля из 4':len(car)/4,'ClinVar':sig,'ClinVar review':r.get('clinvar__rev_stat'),'ClinVar заболевание':r.get('clinvar__disease_names'),'ClinVar ID':r.get('clinvar__id'),'gnomAD AF max':af,'gnomAD4 AF':r.get('gnomad4__af'),'gnomAD4 AC':r.get('gnomad4__ac'),'gnomAD4 AN':r.get('gnomad4__an'),'Отбор':reason,'ССЗ-контекст ClinGen':'; '.join(f'{d} ({cl}; {mo})' for d,cl,mo in ctx),'SIFT':r.get('sift__score'),'SIFT prediction':r.get('sift__prediction'),'SIFT transcript':r.get('sift__transcript'),'PolyPhen':r.get('polyphen2__hvar_pred'),'CADD PHRED':r.get('cadd__phred'),'REVEL':r.get('revel__score'),'REVEL transcript':r.get('revel__transcript'),'AlphaMissense':r.get('alphamissense__am_pathogenicity'),'AlphaMissense transcript':r.get('alphamissense__transcript_id'),'FATHMM':r.get('fathmm__fathmm_score'),'FATHMM transcript':r.get('fathmm__ens_tid'),'dbSNP':r.get('dbsnp__rsid'),'ClinGen заболевание':g.get('clingen__disease'),'ClinGen классификация':g.get('clingen__classification'),'HPO':g.get('hpo__term'),'OMIM ID':r.get('omim__omim_id'),'GTEx gene':r.get('gtex__gtex_gene'),'GTEx tissue':r.get('gtex__gtex_tissue'),'UniProt':g.get('uniprot__acc'),'NCBI Gene ID':g.get('ncbigene__entrez'),'PubMed count':g.get('pubmed__n'),'LoFtool':g.get('loftool__loftool_score'),'pHI':g.get('phi__phi'),'Prec':g.get('prec__prec'),'Статус OpenCRAVAT':'Есть строка результата','Источник варианта':a}
  for suf in ('DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG','DP_DL'):z['SpliceAI '+suf]=r.get('spliceai__'+suf.lower())
  for s in SAMPLES:z['GT '+s]=calls[s].get('GT')
  z['Оценки и транскрипты источников']={'VEP':vr,'OpenCRAVAT':{k:r[k] for k in r if k.startswith(('sift__','polyphen2__','cadd__','revel__','alphamissense__','fathmm__','spliceai__')) and r[k] is not None}}
  if vr and str(vr.get('transcript','')).split('.')[0]==str(t[3] or '').split('.')[0]:
   if z['SIFT'] is None:
    match=re.fullmatch(r'(.*?)\(([0-9.eE+-]+)\)',vr.get('SIFT',''))
    if match:z['SIFT']=num(match[2]);z['SIFT prediction']=match[1];z['SIFT transcript']=vr['transcript']
   if z['PolyPhen'] is None:z['PolyPhen']=vr.get('PolyPhen');z['PolyPhen transcript']=vr['transcript']
   if z['AlphaMissense'] is None and num(vr.get('AlphaMissense_score')) is not None:z['AlphaMissense']=num(vr['AlphaMissense_score']);z['AlphaMissense transcript']=vr.get('AlphaMissense_source_transcript') or vr['transcript']
  if z['CADD PHRED'] is None:z['CADD PHRED']=num(vr.get('CADD_PHRED'))
  if z.get('PolyPhen transcript') is None and r.get('polyphen2__hvar_pred') is not None:z['PolyPhen transcript']=r.get('base__transcript')
  z['_calls']=calls;z['_codes']=list(codes);z['_categories']=list(categories);z['_cv']=bool(ctx);z['_p']=p;z['_af']=af
  if not z['Consequence выбранного транскрипта'] and vr.get('gene')==gene:z['Consequence выбранного транскрипта']=vr.get('consequence')
  if not z['ClinVar review']:z['ClinVar review']=vr.get('clinvar_review')
  if not z['ClinVar заболевание']:z['ClinVar заболевание']=vr.get('clinvar_conditions')
  if not z['ClinVar ID']:z['ClinVar ID']=vr.get('clinvar_ids')
  all_candidates[a]=z
if not n or missing: raise ValueError(f'VCF/annotation mismatch: records={n}, unmatched={len(missing)}; inspect input audit and supply exclusions')
print('VCF scanned',n,'QC genotypes',quality_eligible,'candidates',len(all_candidates),'CVD gene contexts',len(cv),flush=True)
def key(z,sample=None):
 sig=str(z.get('ClinVar') or '').lower();conflict='conflict' in sig;high=z['_af'] is not None and z['_af']>.01;benign='benign' in sig and not z['_p']
 severe=any(x in z['_categories'] for x in ('stop_gained','frameshift_variant','splice_acceptor_variant','splice_donor_variant'))
 c=z['_calls'].get(sample,{}) if sample else {};gq=num(c.get('GQ')) if sample else z['GQ'];dp=num(c.get('DP')) if sample else z['DP']
 model_available=any(z.get(k) is not None for k in ('SIFT','PolyPhen','CADD PHRED','REVEL','AlphaMissense','FATHMM','SpliceAI DS_AG','SpliceAI DS_AL','SpliceAI DS_DG','SpliceAI DS_DL'))
 return (int(high),-int(z['_cv']),0 if z['_p'] and not conflict else 3 if benign and not conflict else 2 if conflict else 1,0 if severe else 1 if 'missense_variant' in z['_categories'] else 2,int(z['_af'] is None),z['_af'] if z['_af'] is not None else 1,-min(gq or 0,99),-min(dp or 0,100),-int(model_available),z['Хромосома'],z['Позиция'],z['REF'],z['ALT'])
ordered=sorted(all_candidates.values(),key=key)
def row(z,sample=None,rank=None):
 z=dict(z);z['Ранг']=rank;z['Категории по всем транскриптам']='; '.join(sorted(set(z['_categories'])))
 if sample:
  c=z['_calls'][sample];z.update({'Образец метрик':sample,'GT':c.get('GT'),'DP':num(c.get('DP')),'AD':c.get('AD'),'GQ':num(c.get('GQ')),'Фильтр генотипа':c.get('FT')})
 ad=str(z.get('AD') or '').split(',')
 if len(ad)==2 and num(ad[0]) is not None and num(ad[1]) is not None and num(ad[0])+num(ad[1])>0:z['Доля ALT в AD']=num(ad[1])/(num(ad[0])+num(ad[1]))
 model_available=any(z.get(k) is not None for k in ('SIFT','PolyPhen','CADD PHRED','REVEL','AlphaMissense','FATHMM','SpliceAI DS_AG','SpliceAI DS_AL','SpliceAI DS_DG','SpliceAI DS_DL'))
 z['Критерии рейтинга']='; '.join(['AF > 1%' if z['_af'] is not None and z['_af']>.01 else 'AF неизвестна' if z['_af'] is None else 'AF ≤ 1%','есть связь гена с ССЗ в ClinGen' if z['_cv'] else 'связь гена с ССЗ в панели не установлена',z['Отбор'],'есть оценки моделей' if model_available else 'нет оценок моделей'])
 return [clean(z.get(h)) for h in H]
sheets=[];counts={}
def add(name,head,values):
 assert len(name)<=31;sheets.append({'name':name,'headers':head,'values':values});counts[name]=len(values)
add('Все_прошедшие',H,[row(z,rank=i) for i,z in enumerate(ordered,1)])
for sheet_name,sel in [('Общие_2_и_более',[z for z in ordered if z['Число носителей']>=2]),('Общие_все_4',[z for z in ordered if z['Число носителей']==4])]:add(sheet_name,H,[row(z,rank=i) for i,z in enumerate(sel,1)])
for sample in SAMPLES:
 pool=sorted([z for z in ordered if sample in z['Образцы-носители'].split('; ')],key=lambda z:key(z,sample));add(sample,H,[row(z,sample,i) for i,z in enumerate(pool,1)])
cats={'Stop_gained':['stop_gained'],'Stop_lost':['stop_lost'],'Start_lost':['start_lost'],'Frameshift':['frameshift_variant','frameshift_truncation','frameshift_elongation'],'Splice_acceptor':['splice_acceptor_variant'],'Splice_donor':['splice_donor_variant'],'Splice_region':['splice_region_variant','splice_site_variant','splice_donor_5th_base_variant','splice_donor_region_variant','splice_polypyrimidine_tract_variant'],'Missense':['missense_variant'],'Synonymous':['synonymous_variant'],'Inframe':['inframe_insertion','inframe_deletion'],'pLoF_candidates':['pLoF_candidates']}
for name,terms in cats.items():
 pool=[z for z in ordered if set(z['_categories'])&set(terms)];add(name,H,[row(z) for z in pool])
for sample in SAMPLES:
 pool=sorted([z for z in ordered if sample in z['Образцы-носители'].split('; ')],key=lambda z:key(z,sample));functional_terms=set(sum(cats.values(),[]));add('Функц_'+sample,H,[row(z,sample) for z in pool if set(z['_codes'])&allowed or set(z['_categories'])&functional_terms])
# X/Y is preserved as all QC calls, with per-sample call metrics restored from the source VCF.
xy=old['X_Y_все']
if 'HCM-контекст' in xy['headers']:
 drop=xy['headers'].index('HCM-контекст');xy={'headers':[v for j,v in enumerate(xy['headers']) if j!=drop],'values':[[v for j,v in enumerate(r) if j!=drop] for r in xy['values']]}
xyh={v:i for i,v in enumerate(xy['headers'])};calls_by={};xy_header_seen=False
with args.xy_quality_vcf.open(encoding='utf8') as f:
 for line in f:
  if line.startswith('#CHROM'):
   if line.rstrip('\n').split('\t')[9:]!=SAMPLES:raise ValueError('X/Y VCF sample order differs from --samples')
   xy_header_seen=True;continue
  if line.startswith('#'):continue
  if not xy_header_seen:raise ValueError('X/Y VCF is missing the #CHROM header')
  v=line.rstrip().split('\t');a=':'.join((v[0],v[1],v[3],v[4]));calls_by[a]={s:dict(zip(v[8].split(':'),x.split(':'))) for s,x in zip(SAMPLES,v[9:])}
if not xy_header_seen:raise ValueError('X/Y VCF is missing the #CHROM header')
add('X_Y_все',xy['headers'],xy['values'])
for sample in SAMPLES:
 out=[]
 for orig in xy['values']:
  if sample not in (orig[xyh['Образцы-носители']] or '').split('; '):continue
  r=list(orig);c=calls_by[r[xyh['Источник варианта']]][sample]
  for label,val in [('Образец метрик',sample),('GT',c.get('GT')),('DP',num(c.get('DP'))),('AD',c.get('AD')),('GQ',num(c.get('GQ'))),('Фильтр генотипа',c.get('FT'))]:r[xyh[label]]=val
  r[xyh['Ранг']]=None;out.append(r)
 add('X_Y_'+sample,xy['headers'],out)
for name in ['MT_все']+['MT_'+s for s in SAMPLES]:add(name,old[name]['headers'],old[name]['values'])
for sheet in sheets:
 if sheet['name'] in cats or sheet['name'].startswith('Функц_'):
  rank_col=sheet['headers'].index('Ранг')
  for i,values in enumerate(sheet['values'],1):values[rank_col]=i
summary=[['Аутосомные записи после QC',n],['QC генотипов после исключений',quality_eligible],['Прошли исследовательский отбор',len(ordered)],['С носителями минимум в двух образцах',sum(z['Число носителей']>=2 for z in ordered)],['С носителями во всех четырёх',sum(z['Число носителей']==4 for z in ordered)],['Гены с проверенной кураторской связью ССЗ',len(cv)]]+[[s+' — полный список',counts[s]] for s in SAMPLES]+[['Строк на листе '+k,v] for k,v in counts.items()]
sheets.insert(0,{'name':'Сводка','headers':['Показатель','Значение'],'values':summary})
ann=[['Фенотип','Сердечно-сосудистые заболевания в целом; конкретные диагнозы и родство не предоставлены.'],['Отбор','Чёткая классификация ClinVar P/LP или функциональный вариант с AF ≤ 1% либо отсутствующей AF.'],['Частота','Максимум доступных gnomAD4 AF и популяционных AF, а для ранее аннотированных вариантов также исходной максимальной AF; источники могут отличаться.'],['ССЗ','ClinGen gene–disease validity, классификации definitive/strong/moderate; поиск кардиологических, сосудистых и аритмогенных терминов в диагнозе или экспертной группе. Список не исчерпывает все гены ССЗ.'],['Рейтинг','AF > 1% понижает ранг; затем связь ClinGen с ССЗ, ClinVar, тип последствия, AF, GQ и DP. Это исследовательская эвристика, не ACMG.'],['Последствия','Объединение всех транскриптов OpenCRAVAT и доступных функциональных классов VEP; выбранный транскрипт указан отдельно.'],['SpliceAI','Оценки уровня аллеля; ген модели не указан в скачанном результате и не приписывается автоматически выбранному гену.'],['Генотип','Носителями считаются только альтернативные генотипы с FT PASS, DP ≥ 10 и GQ ≥ 20; данные строки относятся к указанному образцу метрик.'],['X/Y','Все QC вызовы и полные листы по образцам; без установления кариотипа.'],['MT','PASS вызовы; MT AF/AD/DP отдельно, без ядерного фильтра GQ.'],['ACMG','Итоговый автоматический ACMG-класс не присваивался'],['Репозиторий','https://github.com/MarFedorovna/marina-wgs-analysis']]
add('Аннотация_источники',['Источник или правило','Сведения'],ann)
out={'sheets':sheets,'stats':{'vcf_records':n,'quality_eligible':quality_eligible,'candidates':len(ordered),'cardio_candidates':sum(z['_cv'] for z in ordered),'sample_counts':{s:counts[s] for s in SAMPLES},'shared':sum(z['Число носителей']>=2 for z in ordered),'all4':sum(z['Число носителей']==4 for z in ordered),'sheet_counts':counts}}
args.output_json.parent.mkdir(parents=True,exist_ok=True)
args.output_json.write_text(json.dumps(out,ensure_ascii=False),encoding='utf8')
print('FINISHED',json.dumps(out['stats'],ensure_ascii=True),flush=True)

