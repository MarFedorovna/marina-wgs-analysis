#!/usr/bin/env python3
"""Local, receipt-bound sample VCFs and exact-allele research report views."""
import argparse
from contextlib import ExitStack
import csv
import datetime as dt
import gzip
import hashlib
import html
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3

ABSOLUTE_DEADLINE = dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
CORE = ['stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant', 'splice_acceptor_variant',
        'splice_donor_variant', 'splice_region_variant', 'missense_variant', 'synonymous_variant',
        'inframe_insertion', 'inframe_deletion', 'pLoF_candidates']
PLOF = {'stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant'}
PROJECTION = ['sample', 'GT', 'FT', 'DP', 'AD', 'GQ', 'QUAL', 'sample_state', 'site_filter',
    'gene', 'gene_id', 'hgnc_id', 'transcript', 'protein_id', 'hgvs_c', 'hgvs_p', 'consequence',
    'amino_acids', 'protein_position', 'protein_length', 'protein_region_heuristic', 'codons', 'exon', 'intron',
    'mane_select', 'canonical', 'stop_gain', 'stop_loss', 'start_lost', 'frameshift', 'splice', 'missense', 'synonymous', 'plof_candidate',
    'nmd_assessment', 'mechanism_status', 'classification_status', 'clinvar_significance', 'clinvar_review',
    'clinvar_conflicts', 'clinvar_conditions', 'clinvar_ids', 'clinvar_snapshot', 'clingen_evidence_json',
    'gene_context_json', 'af_max_observed', 'af_status', 'af_fields_json', 'ancestry_af_json', 'gnomad_version',
    'SIFT', 'PolyPhen', 'FATHMM', 'CADD_PHRED', 'AlphaMissense_score', 'AlphaMissense_match_status',
    'AlphaMissense_allele_predictions_json', 'local_evidence_json', 'resource_versions_json',
    'candidate_reason', 'ranking_rule']
SUMMARY = ['allele_id', 'chrom', 'pos', 'ref', 'alt', 'carrier_count', 'pct_of4', 'raw_observed_alt_count',
    'raw_pct_of4', 'assessed_count', 'eligible_sample_count', 'quality_passing_carriers_json', 'raw_observed_alt_samples_json',
    'eligible_samples_json', 'reference_samples_json', 'unknown_samples_json', 'filtered_samples_json',
    'partial_samples_json', 'unassessed_samples_json', 'sample_calls_json', 'cohort_AC_source', 'cohort_AN_source',
    'cohort_AF_source', 'annotation_transcript_rows', 'union_consequences_json', 'representative_sample',
    *PROJECTION, 'functional_membership_basis', 'ACMG_class', 'ACMG_criteria_evidence_status', 'ACMG_provisional_review',
    'all_transcript_evidence', 'evidence_join', 'interpretation_status']


def js(value): return json.dumps(value, ensure_ascii=False, separators=(',', ':'), sort_keys=True)
def load(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def now(): return dt.datetime.now(dt.timezone.utc)


class Guard:
    def __init__(self, deadline=ABSOLUTE_DEADLINE, hours=8):
        self.deadline = min(deadline, ABSOLUTE_DEADLINE, now() + dt.timedelta(hours=hours))
    def check(self):
        if now() >= self.deadline: raise ValueError('Original absolute budget deadline reached; no automatic retry')


def sha(path, guard=None):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            if guard: guard.check()
            h.update(block)
    return h.hexdigest()


def exact_key(chrom, pos, ref, alt):
    if ',' in alt: raise ValueError('Expected already normalized decomposed biallelic source')
    return js([chrom, int(pos), ref, alt])


def allele_id(chrom, pos, ref, alt):
    return ':'.join((chrom, str(int(pos)), ref, alt))


def call_flags(gt, ft, site_filter, sequence=True):
    tokens = re.split(r'[/|]', gt or '.')
    if not tokens or any(t not in ('0', '1', '.') for t in tokens):
        raise ValueError('Non-biallelic or malformed GT in decomposed source')
    raw = '1' in tokens
    partial = '.' in tokens
    filtered = ft not in ('', '.', 'PASS') or site_filter not in ('', '.', 'PASS')
    unknown = partial or ft in ('', '.') or site_filter in ('', '.')
    assessed = sequence and not partial and ft == site_filter == 'PASS'
    return dict(raw=raw, qc=raw and assessed, assessed=assessed, partial=partial, filtered=filtered,
                unknown=unknown, reference=assessed and not raw, unassessed=not assessed,
                sequence_assessed=sequence, AC=tokens.count('1'), AN=sum(t != '.' for t in tokens))


def parse_info(value):
    return dict((field.split('=', 1) + [None])[:2] for field in value.split(';') if field and field != '.')


def derived_info(info, flags):
    result = dict(info)
    for tag in ('AC', 'AN', 'AF', 'NS', 'DP'):
        if tag in result: result['COHORT_' + tag] = result.pop(tag)
    result.update(AC=str(flags['AC']), AN=str(flags['AN']), AF=format(flags['AC']/flags['AN'], '.12g') if flags['AN'] else '.',
                  NS=str(int(flags['AN'] > 0)), VIEW_OBSERVED_ALT=None,
                  VIEW_QUALITY_STATE='PARTIAL_UNKNOWN' if flags['partial'] else 'FILTERED' if flags['filtered'] else 'QUALITY_UNKNOWN' if flags['unknown'] else 'UNASSESSED_NON_SEQUENCE' if not flags['sequence_assessed'] else 'PASS')
    return ';'.join(k if v is None else k+'='+v for k,v in result.items())


def subset_header(headers, sample):
    result = []
    for line in headers:
        match = re.match(r'##INFO=<ID=(AC|AN|AF|NS|DP),', line)
        if match:
            tag = match[1]
            result.append(line.replace('ID='+tag+',', 'ID=COHORT_'+tag+',', 1))
        elif line.startswith('#CHROM'):
            continue
        else:
            if line.startswith('##INFO=<ID=COHORT_') or line.startswith('##INFO=<ID=VIEW_'):
                raise ValueError('Reserved subset provenance tag already in source')
            result.append(line)
    for tag, number, kind, description in [
        ('AC','A','Integer','ALT allele count in this one-sample subset; called GT tokens, including partial/filtered calls'),
        ('AN','1','Integer','Called allele tokens in this one-sample subset; missing tokens excluded'),
        ('AF','A','Float','AC/AN from subset GT; not population gnomAD AF or quality-eligible carrier fraction'),
        ('NS','1','Integer','Samples with at least one called GT allele in this subset'),
        ('VIEW_OBSERVED_ALT','0','Flag','Source GT contains ALT; partial/filtered calls retained'),
        ('VIEW_QUALITY_STATE','1','String','Explicit subset call quality state; source FILTER and FORMAT remain unchanged')]:
        result.append('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % (tag,number,kind,description))
    result.append('##cohort_sample_views="Original INFO counts/depth preserved as COHORT_*; remaining site annotations retain cohort/source provenance; research only"\n')
    result.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')
    return ''.join(result).encode()


def membership(calls, samples, eligible=()):
    selected = lambda name: [s for s in samples if calls[s]['flags'][name]]
    carriers, raw = selected('qc'), selected('raw')
    return dict(carrier_count=len(carriers), pct_of4=25*len(carriers), raw_observed_alt_count=len(raw), raw_pct_of4=25*len(raw),
        assessed_count=len(selected('assessed')), eligible_sample_count=len(eligible), quality_passing_carriers_json=js(carriers),
        raw_observed_alt_samples_json=js(raw), eligible_samples_json=js(sorted(eligible)), reference_samples_json=js(selected('reference')),
        unknown_samples_json=js(selected('unknown')), filtered_samples_json=js(selected('filtered')), partial_samples_json=js(selected('partial')),
        unassessed_samples_json=js(selected('unassessed')), sample_calls_json=js(calls))


def scope_members(scope, calls, eligible):
    raw = {s for s,c in calls.items() if c['flags']['raw']}
    qc = {s for s,c in calls.items() if c['flags']['qc']}
    if scope == 'cohort': return bool(raw), bool(qc), bool(eligible)
    if scope == 'shared_ge2': return len(raw)>=2, len(qc)>=2, len(qc)>=2 and bool(eligible)
    if scope == 'all4': return len(raw)==4, len(qc)==4, len(qc)==4 and bool(eligible)
    return scope in raw, scope in qc, scope in eligible


def ranked_insert(top, row, key):
    top.append(row); top.sort(key=key)
    del top[20:]


def checked_inputs(workflow, guard):
    root = Path(workflow).resolve()
    receipts = {stage: root/'stages'/stage/'workflow-receipt.json' for stage in ('06_annotate','07_report')}
    values = {name:load(path) for name,path in receipts.items()}
    if any(v.get('status') != 'PASS' for v in values.values()): raise ValueError('NOT_READY: stage6 AND stage7 PASS required')
    sources = {'vcf':root/'stages/06_annotate/annotated.vcf.gz', 'vcf_index':root/'stages/06_annotate/annotated.vcf.gz.tbi',
               'master':root/'stages/07_report/master_evidence.tsv.gz'}
    pins = {str(p):dict(bytes=p.stat().st_size,sha256=sha(p,guard)) for p in receipts.values()}
    for name,path in sources.items():
        stage = '07_report' if name=='master' else '06_annotate'
        expected = values[stage].get('artifacts',{}).get(path.name,{})
        if path.is_symlink() or path.stat().st_size != expected.get('bytes') or sha(path,guard) != expected.get('sha256'):
            raise ValueError('Completed artifact differs from receipt: '+str(path))
        pins[str(path)] = expected
    return sources,pins


def export(workflow, output, ranking_module, ranking_sha256, *, deadline=ABSOLUTE_DEADLINE, synthetic=False):
    guard = Guard(deadline); guard.check()
    if hasattr(os,'sched_getaffinity'):
        os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
    try:
        import resource
        cap=4*1024**3;soft,hard=resource.getrlimit(resource.RLIMIT_AS)
        resource.setrlimit(resource.RLIMIT_AS,(min(cap,soft) if soft>0 else cap,hard))
    except ImportError:
        pass
    workflow, output, ranking_module = Path(workflow).resolve(), Path(output).resolve(), Path(ranking_module).resolve()
    if output.exists() or output==workflow or workflow in output.parents: raise ValueError('New sidecar output outside workflow required; refusing overwrite')
    if sha(ranking_module,guard) != ranking_sha256: raise ValueError('Pinned sealed ranking module SHA differs')
    spec=importlib.util.spec_from_file_location('sealed_views_ranking',ranking_module); sealed=importlib.util.module_from_spec(spec);spec.loader.exec_module(sealed)
    sources,pins=checked_inputs(workflow,guard)
    import pysam
    output.mkdir(mode=0o700,parents=True); (output/'vcf').mkdir(); (output/'masters').mkdir(); (output/'tables').mkdir()
    db=sqlite3.connect(output/'.working.sqlite');db.row_factory=sqlite3.Row
    db.executescript('''PRAGMA journal_mode=DELETE;PRAGMA temp_store=FILE;PRAGMA cache_size=-65536;
      CREATE TABLE alleles(aid INTEGER PRIMARY KEY,key TEXT UNIQUE,chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,calls TEXT,info TEXT,processed INTEGER DEFAULT 0);
      CREATE TABLE gene_members(gene TEXT,sample TEXT,aid INTEGER,raw INTEGER,qc INTEGER,PRIMARY KEY(gene,sample,aid)) WITHOUT ROWID;''')
    def progress(status,**extra):
        p=output/'progress.json';temp=output/'progress.tmp';temp.write_text(json.dumps(dict(status=status,updated_utc=now().isoformat(),pid=os.getpid(),**extra),indent=2)+'\n');temp.replace(p)
    try:
        progress('READING_COHORT_VCF_WRITING_RAW_OBSERVED_ALT_SAMPLE_VCFS')
        with gzip.open(sources['vcf'],'rt') as stream, ExitStack() as stack:
            headers=[]
            for line in stream:
                headers.append(line)
                if line.startswith('#CHROM'):break
            samples=headers[-1].rstrip('\r\n').split('\t')[9:]
            if len(samples)!=4 or len(set(samples))!=4 or any(not s or '\t' in s for s in samples):raise ValueError('Exactly four private VCF-header sample IDs required')
            if set(samples)&{'cohort','shared_ge2','all4'}:raise ValueError('Sample ID collides with a reserved report scope; choose an explicitly reviewed scope mapping')
            labels={sample:'sample_%04d'%(i+1) for i,sample in enumerate(samples)}
            handles={s:stack.enter_context(pysam.BGZFile(str(output/'vcf'/(labels[s]+'.annotated.vcf.gz')),'wb')) for s in samples}
            for sample,handle in handles.items():handle.write(subset_header(headers,sample))
            counts={s:0 for s in samples};vcf_count=0
            for line in stream:
                fields=line.rstrip('\r\n').split('\t')
                if len(fields)!=13:raise ValueError('VCF column count differs from four-sample header')
                chrom,pos,_,ref,alt=fields[:5];key=exact_key(chrom,pos,ref,alt);names=fields[8].split(':')
                if len(set(names))!=len(names) or 'GT' not in names:raise ValueError('Malformed/absent FORMAT GT')
                calls={};real=bool(re.fullmatch('[ACGT]+',ref) and re.fullmatch('[ACGT]+',alt))
                for sample,raw in zip(samples,fields[9:]):
                    values=raw.split(':');data={name:(values[i] if i<len(values) else '.') for i,name in enumerate(names)}
                    flags=call_flags(data['GT'],data.get('FT','.'),fields[6],real)
                    calls[sample]=dict(GT=data['GT'],FT=data.get('FT','.'),DP=data.get('DP','.'),AD=data.get('AD','.'),GQ=data.get('GQ','.'),
                        site_filter=fields[6],flags=flags,source_sample_column=raw)
                    if flags['raw']:
                        derived=list(fields[:9]);derived[7]=derived_info(parse_info(fields[7]),flags)
                        handles[sample].write(('\t'.join(derived+[raw])+'\n').encode());counts[sample]+=1
                db.execute('INSERT INTO alleles(key,chrom,pos,ref,alt,calls,info) VALUES(?,?,?,?,?,?,?)',
                    (key,chrom,int(pos),ref,alt,js(calls),js(parse_info(fields[7]))));vcf_count+=1
                if vcf_count%5000==0:guard.check();db.commit();progress('SPLITTING_VCF',cohort_records=vcf_count,sample_records=counts)
            db.commit()
        for sample in samples:
            guard.check();path=output/'vcf'/(labels[sample]+'.annotated.vcf.gz')
            pysam.tabix_index(str(path),preset='vcf',force=False)
            with pysam.VariantFile(str(path)) as check:
                if list(check.header.samples)!=[sample]:raise ValueError('Subset VCF sample header differs')
                read_count=0
                for record in check.fetch():
                    if read_count%10000==0:guard.check()
                    gt=record.samples[sample]['GT'];ac=sum(x==1 for x in gt);an=sum(x is not None for x in gt)
                    if ac<1 or tuple(record.info['AC'])!=(ac,) or record.info['AN']!=an or abs(record.info['AF'][0]-ac/an)>1e-6:
                        raise ValueError('BGZF/TBI subset AC/AN/AF validation failed')
                    read_count+=1
                if read_count!=counts[sample]:raise ValueError('Indexed subset record count differs')
        progress('GROUPING_MASTER_BY_EXACT_ALLELE_BOUNDED_TRANSCRIPT_AGGREGATION',vcf_records=vcf_count)
        scopes=['cohort','shared_ge2','all4',*samples]; tops={scope:[] for scope in scopes}; stats={scope:{'raw_alleles':0,'qc_alleles':0,'eligible_alleles':0,'consequences':{}} for scope in scopes}
        matrix_fields=['allele_id','chrom','pos','ref','alt',*SUMMARY[5:24]]
        with ExitStack() as stack:
            def writer(path,fields,gz=True):
                stream=stack.enter_context(gzip.open(path,'wt',encoding='utf-8',newline='') if gz else Path(path).open('x',encoding='utf-8',newline=''))
                result=csv.DictWriter(stream,fields,delimiter='\t',extrasaction='ignore',lineterminator='\n');result.writeheader();return result
            matrix=writer(output/'tables/allele-membership-matrix.tsv.gz',matrix_fields)
            writers={scope:writer(output/'tables'/((labels.get(scope,scope))+'.functional.tsv.gz'),[*SUMMARY,'scope_raw_member','scope_qc_member','scope_candidate_member']) for scope in scopes}
            category_writers={}
            def category_writer(scope,term):
                key=scope,term
                if key not in category_writers:
                    folder=output/'tables'/labels.get(scope,scope);folder.mkdir(exist_ok=True)
                    slug=re.sub('[^A-Za-z0-9_-]','_',term)[:80]+'-'+hashlib.sha256(term.encode()).hexdigest()[:8]
                    category_writers[key]=writer(folder/(slug+'.tsv.gz'),[*SUMMARY,'scope_raw_member','scope_qc_member','scope_candidate_member'])
                return category_writers[key]
            group=None;source_rows=0;annotated_alleles=0
            def finish(group):
                nonlocal annotated_alleles
                if group is None:return
                allele=group['allele'];calls=group['calls'];eligible=set(group['best_eligible']);terms=group['terms'] or {'unannotated'}
                if terms&PLOF:terms.add('pLoF_candidates')
                if group['seen_samples']!=set(samples):raise ValueError('Master allele missing an explicit sample row; missing is not reference')
                best=min(group['best_eligible'].values(),key=sealed.rank_key) if eligible else group['representative']
                info=json.loads(allele['info']);summary=dict(allele_id=allele_id(allele['chrom'],allele['pos'],allele['ref'],allele['alt']),
                    chrom=allele['chrom'],pos=allele['pos'],ref=allele['ref'],alt=allele['alt'],**membership(calls,samples,eligible),
                    cohort_AC_source=info.get('AC','.'),cohort_AN_source=info.get('AN','.'),cohort_AF_source=info.get('AF','.'),
                    annotation_transcript_rows=group['rows'],union_consequences_json=js(sorted(terms)),representative_sample=best['sample'],
                    functional_membership_basis='ANY_SOURCE_TRANSCRIPT; allele-level QC/candidate membership; categories overlap; see complete transcript appendix',
                    ACMG_class='NOT_ASSIGNED',ACMG_criteria_evidence_status='NOT_ASSESSED_MANUAL_EVIDENCE_REQUIRED',
                    ACMG_provisional_review='PENDING_GENE_DISEASE_TRANSCRIPT_QUALITY_AND_CRITERION_REVIEW',
                    all_transcript_evidence='masters/cohort.master.tsv.gz (relative to export root)',evidence_join='exact chrom,pos,ref,alt + sample; transcript rows not collapsed in appendix',
                    interpretation_status='RESEARCH_VIEWS_NOT_CLINICAL_BURDEN_ACMG_OR_MECHANISM')
                summary.update({field:best.get(field,'') for field in PROJECTION});matrix.writerow(summary)
                for gene in group['genes']:
                    for sample,call in calls.items():
                        if call['flags']['raw']:db.execute('INSERT INTO gene_members VALUES(?,?,?,?,?)',(gene,sample,allele['aid'],1,int(call['flags']['qc'])))
                for scope in scopes:
                    raw,qc,candidate=scope_members(scope,calls,eligible)
                    if raw:
                        scoped=dict(summary)
                        if scope in samples:
                            representative=group['best_eligible'].get(scope,group['best_by_sample'][scope])
                            scoped.update({field:representative.get(field,'') for field in PROJECTION});scoped['representative_sample']=scope
                        scoped.update(scope_raw_member=int(raw),scope_qc_member=int(qc),scope_candidate_member=int(candidate))
                        writers[scope].writerow(scoped)
                        for term in terms:category_writer(scope,term).writerow(scoped)
                    for level,present in [('raw',raw),('qc',qc),('eligible',candidate)]:
                        if present:
                            stats[scope][level+'_alleles']+=1
                            for term in terms:
                                stat=stats[scope]['consequences'].setdefault(term,dict(raw=0,qc=0,eligible=0));stat[level]+=1
                    selected=group['best_eligible'].get(scope) if scope in samples else best if candidate else None
                    if selected is not None:
                        ranked_insert(tops[scope],dict(selected,_summary=summary),lambda r:sealed.rank_key(r))
                db.execute('UPDATE alleles SET processed=1 WHERE aid=?',(allele['aid'],));annotated_alleles+=1
                if annotated_alleles%1000==0:
                    guard.check()
                    if shutil.disk_usage(output).free < 64*1024**3 + sources['master'].stat().st_size:
                        raise ValueError('Sidecar disk reserve reached; primary workflow is not modified')
                    db.commit();progress('GROUPING_MASTER',master_rows=source_rows,annotated_alleles=annotated_alleles)
            with gzip.open(sources['master'],'rt',encoding='utf-8',newline='') as stream:
                reader=csv.DictReader(stream,delimiter='\t');source_fields=reader.fieldnames or []
                required={'sample','variant_id','chrom','pos','ref','alt','GT','FT','site_filter','transcript','consequence','candidate_eligible'}
                if not required<=set(source_fields) or len(source_fields)!=len(set(source_fields)):raise ValueError('Master columns missing/duplicated')
                sample_masters={sample:writer(output/'masters'/(labels[sample]+'.master.tsv.gz'),source_fields) for sample in samples}
                for row in reader:
                    source_rows+=1
                    if source_rows%10000==0:guard.check()
                    if None in row or any(v is None for v in row.values()):raise ValueError('Malformed master TSV row')
                    key=exact_key(row['chrom'],row['pos'],row['ref'],row['alt'])
                    if group is None or group['key']!=key:
                        finish(group)
                        allele=db.execute('SELECT * FROM alleles WHERE key=?',(key,)).fetchone()
                        if allele is None or allele['processed']:raise ValueError('Master allele absent from VCF or appears in noncontiguous groups')
                        group=dict(key=key,allele=allele,calls=json.loads(allele['calls']),terms=set(),genes=set(),rows=0,seen_samples=set(),best_eligible={},best_by_sample={},representative=row)
                    sample=row['sample']
                    if sample not in samples:raise ValueError('Master sample differs from private VCF header')
                    call=group['calls'][sample]
                    if row['variant_id']!=sealed.variant_key(row['chrom'],row['pos'],row['ref'],row['alt']):raise ValueError('Master normalized variant identity differs from sealed contract')
                    if row['GT']!=call['GT'] or row['FT']!=call['FT'] or row['site_filter']!=call['site_filter']:
                        raise ValueError('Master genotype/quality differs from annotated VCF')
                    for field in ('DP','AD','GQ'):
                        if field in row and row[field]!=call[field]:raise ValueError('Master FORMAT field differs: '+field)
                    group['seen_samples'].add(sample);group['rows']+=1
                    previous=group['best_by_sample'].get(sample)
                    if previous is None or sealed.rank_key(row)<sealed.rank_key(previous):group['best_by_sample'][sample]=row
                    if call['flags']['raw']:sample_masters[sample].writerow(row)
                    group['terms'].update(term for term in re.split(r'[&,;|]',row['consequence']) if term)
                    if row.get('gene'):group['genes'].add(row['gene'])
                    if call['flags']['raw'] and (not group['calls'][group['representative']['sample']]['flags']['raw'] or sealed.rank_key(row)<sealed.rank_key(group['representative'])):
                        group['representative']=row
                    if row['candidate_eligible']=='1':
                        if not call['flags']['qc']:raise ValueError('Source candidate eligibility contradicts quality-passing complete ALT')
                        previous=group['best_eligible'].get(sample)
                        if previous is None or sealed.rank_key(row)<sealed.rank_key(previous):group['best_eligible'][sample]=row
                    elif row['candidate_eligible']!='0':raise ValueError('Unknown sealed candidate eligibility state')
                finish(group)
            db.commit()
            if db.execute('SELECT COUNT(*) FROM alleles WHERE processed=0').fetchone()[0]:raise ValueError('VCF allele missing from full master evidence')
            for scope in scopes:
                path=output/'tables'/((labels.get(scope,scope))+'.top20.tsv')
                fields=['rank',*SUMMARY,*[f for f in source_fields if f not in SUMMARY], 'scope_raw_member','scope_qc_member','scope_candidate_member']
                top=writer(path,fields,False)
                for number,row in enumerate(tops[scope],1):
                    record=dict(row['_summary']);record.update({k:v for k,v in row.items() if k!='_summary'})
                    record.update(rank=number,representative_sample=row['sample'],scope_raw_member=1,scope_qc_member=1,scope_candidate_member=1)
                    top.writerow(record)
                stats[scope]['top20_unique_alleles']=len(tops[scope])
            worksheet=writer(output/'tables/ACMG-manual-top20-worksheet.tsv.gz',
                ['scope','rank','allele_id','sample','gene','transcript','criterion','status','strength','source','rationale'])
            for scope in ['cohort',*samples]:
                for number,row in enumerate(tops[scope],1):
                    for criterion in sealed.ACMG:
                        worksheet.writerow(dict(scope=labels.get(scope,scope),rank=number,allele_id=row['_summary']['allele_id'],sample=row['sample'],gene=row.get('gene',''),
                            transcript=row.get('transcript',''),criterion=criterion,status='NOT_ASSESSED_MANUAL_REVIEW_REQUIRED',strength='',source='',rationale=''))
            category_index=writer(output/'tables/functional-category-index.tsv',
                ['scope','category','raw_unique_alleles','qc_unique_alleles','eligible_unique_alleles','path'],False)
            for scope in scopes:
                for term in sorted(set(CORE)|set(stats[scope]['consequences'])):
                    category_writer(scope,term)
                    counts_for_term=stats[scope]['consequences'].get(term,dict(raw=0,qc=0,eligible=0))
                    slug=re.sub('[^A-Za-z0-9_-]','_',term)[:80]+'-'+hashlib.sha256(term.encode()).hexdigest()[:8]
                    category_index.writerow(dict(scope=labels.get(scope,scope),category=term,raw_unique_alleles=counts_for_term['raw'],qc_unique_alleles=counts_for_term['qc'],
                        eligible_unique_alleles=counts_for_term['eligible'],path=labels.get(scope,scope)+'/'+slug+'.tsv.gz'))
            gene_table=writer(output/'tables/gene-overlap-SEPARATE-NOT-SHARED-ALLELES.tsv.gz',
                ['gene','raw_sample_allele_counts_json','qc_sample_allele_counts_json','qc_sample_count','scope'])
            previous=None;raw_counts={};qc_counts={}
            def gene_finish():
                if previous is not None:gene_table.writerow(dict(gene=previous,raw_sample_allele_counts_json=js(raw_counts),qc_sample_allele_counts_json=js(qc_counts),qc_sample_count=sum(n>0 for n in qc_counts.values()),scope='GENE_OVERLAP_DIFFERENT_ALLELES_ALLOWED_NOT_ALLELE_INTERSECTION'))
            for record in db.execute('SELECT gene,sample,COUNT(*) AS n,SUM(qc) AS q FROM gene_members GROUP BY gene,sample ORDER BY gene,sample'):
                if record['gene']!=previous:gene_finish();previous=record['gene'];raw_counts={};qc_counts={}
                raw_counts[record['sample']]=record['n'];qc_counts[record['sample']]=record['q']
            gene_finish()
        progress('PRESERVING_BYTE_IDENTICAL_ALL_TRANSCRIPT_APPENDIX')
        appendix=output/'masters/cohort.master.tsv.gz'
        with sources['master'].open('rb') as inp,appendix.open('xb') as out:
            for block in iter(lambda:inp.read(8*1024*1024),b''):guard.check();out.write(block)
        if sha(appendix,guard)!=pins[str(sources['master'])]['sha256']:raise ValueError('All-transcript appendix differs')
        for path,pin in pins.items():
            if Path(path).stat().st_size!=pin['bytes'] or sha(path,guard)!=pin['sha256']:raise ValueError('Input/receipt changed during export')
        if sha(ranking_module,guard)!=ranking_sha256:raise ValueError('Sealed ranking module changed during export')
        db.close();(output/'.working.sqlite').unlink()
        manifest=dict(schema_version=1,status='PASS',synthetic=synthetic,samples=samples,sample_file_labels=labels,
            vcf_records=vcf_count,sample_vcf_records=counts,master_transcript_rows=source_rows,unique_annotated_alleles=annotated_alleles,
            source_pins=pins,ranking_module=dict(path=str(ranking_module),sha256=ranking_sha256,function='rank_key',eligibility='unchanged source candidate_eligible == 1'),
            counts=stats,shared_identity='Exact already-normalized CHROM/POS/REF/ALT; gene overlap is separate',
            carrier_contract='Quality-passing complete ALT GT; raw observed ALT, partial, filtered and unknown retained separately; PCT denominator exactly4',
            shared_top20_contract='At least2 (all4:4) quality-passing ALT carriers and at least1 sealed eligible sample/transcript row; best eligible row uses unchanged sealed rank_key',
            top20_padding=False,automatic_ACMG=False,individual_functional_burden='Descriptive distinct-allele counts only; transcript effects do not prove LoF/GoF/NMD or clinical severity',
            no_external_requests=True,primary_workflow_changed=False,completed_utc=now().isoformat(),pysam_version=pysam.__version__)
        page=['<!doctype html><meta charset="utf-8"><title>Cohort and sample views</title><h1>One annotation: cohort and four sample views</h1>',
            '<p>Research prioritization; no automatic ACMG class or functional burden conclusion. Counts are distinct alleles; consequence groups overlap. Missing and filtered are explicit.</p>',
            '<p><a href="masters/cohort.master.tsv.gz">Complete cohort master: all transcript evidence</a> · <a href="tables/cohort.top20.tsv">Cohort top20</a> · <a href="tables/cohort.functional.tsv.gz">Cohort functional table</a> · <a href="tables/functional-category-index.tsv">Functional type files and counts</a></p>']
        for sample in samples:
            label=labels[sample]
            page.append('<h2>'+html.escape(sample)+'</h2><p>'+ ' · '.join('<a href="'+path+'">'+name+'</a>' for path,name in [
                ('masters/'+label+'.master.tsv.gz','Sample master: observed ALT, all transcript evidence'),
                ('vcf/'+label+'.annotated.vcf.gz','Annotated sample VCF'),('vcf/'+label+'.annotated.vcf.gz.tbi','TBI'),
                ('tables/'+label+'.top20.tsv','Sample top20'),('tables/'+label+'.functional.tsv.gz','Functional consequences')])+'</p>')
        page.append('<details><summary>Auxiliary membership/shared views</summary><p>Exact allele sharing differs from gene overlap.</p><a href="tables/allele-membership-matrix.tsv.gz">Membership matrix</a> · <a href="tables/shared_ge2.top20.tsv">QC shared≥2 top20</a> · <a href="tables/all4.top20.tsv">QC all4 top20</a> · <a href="tables/gene-overlap-SEPARATE-NOT-SHARED-ALLELES.tsv.gz">Separate gene overlap</a></details>')
        page.append('<p><a href="tables/ACMG-manual-top20-worksheet.tsv.gz">Manual ACMG worksheet</a> · <a href="manifest.json">Manifest and source hashes</a></p>')
        (output/'index.html').write_text('\n'.join(page)+'\n',encoding='utf-8')
        progress('PASS',unique_annotated_alleles=annotated_alleles)
        manifest['artifacts']={str(p.relative_to(output)):dict(bytes=p.stat().st_size,sha256=sha(p,guard)) for p in sorted(output.rglob('*')) if p.is_file()}
        (output/'manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
        return manifest
    except BaseException:
        db.close();(output/'INCOMPLETE.txt').write_text('No PASS manifest. Partial report only; primary workflow untouched.\n');raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow-root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--ranking-module',type=Path,required=True);parser.add_argument('--ranking-sha256',required=True)
    parser.add_argument('--synthetic',action='store_true');args=parser.parse_args()
    if hasattr(os,'sched_getaffinity'):
        os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
    print(js(export(args.workflow_root,args.output,args.ranking_module,args.ranking_sha256,synthetic=args.synthetic)))
