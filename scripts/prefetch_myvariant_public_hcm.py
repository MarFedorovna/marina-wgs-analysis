#!/usr/bin/env python3
"""Download fixed public ClinGen HCM gene reference records. No patient input accepted."""
import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import time
import urllib.parse
import urllib.request

CODE=Path(__file__).resolve().parents[1]


def load(path):return json.loads(path.read_text(encoding='utf-8-sig')) if path.is_file() else {}
def now():return dt.datetime.now(dt.timezone.utc)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    if path.is_symlink() or temporary.is_symlink():raise ValueError('Refusing symlink')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def public_genes(cfg):
    path=Path(cfg['clingen_csv'])
    if digest(path)!=cfg['clingen_csv_sha256']:raise ValueError('Public ClinGen source differs')
    lines=path.read_text().splitlines();start=next(i for i,s in enumerate(lines) if '"GENE SYMBOL","GENE ID (HGNC)"' in s)
    all_rows=list(csv.DictReader(lines[start:]));rows=[r for r in all_rows if r.get('DISEASE ID (MONDO)')==cfg['disease_id']]
    if len(rows)!=cfg['expected_curation_rows']:raise ValueError('Exact public HCM curation count differs')
    genes={}
    for r in rows:
        hgnc,symbol=r['GENE ID (HGNC)'],r['GENE SYMBOL']
        if not re.fullmatch(r'HGNC:\d+',hgnc) or not re.fullmatch(r'[A-Z0-9][A-Z0-9.-]{0,50}',symbol):raise ValueError('Unexpected canonical public gene ID')
        if hgnc in genes and genes[hgnc]!=symbol:raise ValueError('Conflicting public gene symbols')
        genes[hgnc]=symbol
    additional=cfg.get('additional_released_hcm_cspec')
    if additional:
        database=Path(additional['sqlite'])
        if digest(database)!=additional['sha256']:raise ValueError('Public CSpec candidate changed')
        canonical={r.get('GENE ID (HGNC)'):r.get('GENE SYMBOL') for r in all_rows if r.get('GENE ID (HGNC)')}
        connection=sqlite3.connect(database.resolve().as_uri()+'?mode=ro',uri=True)
        for hgnc,raw in connection.execute('SELECT hgnc_id,payload FROM cspec WHERE mondo_id=?',(cfg['disease_id'],)):
            content=json.loads(raw);symbol=content['symbol']
            if canonical.get(hgnc)!=symbol or content.get('state')!='Released':raise ValueError('CSpec/current HGNC source identity mismatch')
            genes[hgnc]=symbol
        connection.close()
    return [dict(hgnc_id=k,symbol=v) for k,v in sorted(genes.items(),key=lambda p:p[1])]


class Reference:
    def __init__(self,path):
        self.stream=Path(path).open('rb')
        self.index={v[0]:list(map(int,v[1:5])) for v in (line.split('\t') for line in Path(str(path)+'.fai').read_text().splitlines())}
    def base(self,chrom,pos):
        row=self.index.get('chr'+chrom)
        if not row or not 1<=pos<=row[0]:return None
        _,offset,bases,width=row;n=pos-1;self.stream.seek(offset+(n//bases)*width+(n%bases));return self.stream.read(1).decode().upper()


def validate_hit(hit,gene,reference):
    source=hit.get('dbnsfp',{});names=source.get('genename',[])
    if isinstance(names,str):names=[names]
    if gene not in names:return dict(status='GENE_ID_MISMATCH_NO_SCORE_JOIN')
    chrom=str(source.get('chrom',''));mapping=source.get('hg38',{});pos=mapping.get('start')
    ref,alt=source.get('ref'),source.get('alt')
    if not re.fullmatch(r'(?:[1-9]|1\d|2[0-2]|X|Y)',chrom) or not isinstance(pos,int) or mapping.get('end')!=pos or ref not in ('A','C','G','T') or alt not in ('A','C','G','T') or ref==alt:
        return dict(status='MISSING_OR_AMBIGUOUS_GRCH38_SNV_MAPPING')
    key='chr%s:g.%d%s>%s'%(chrom,pos,ref,alt)
    if hit.get('_id')!=key:return dict(status='RAW_ID_ASSEMBLY_MISMATCH',structured_grch38_key=key)
    if reference.base(chrom,pos)!=ref:return dict(status='REFERENCE_BASE_MISMATCH_NO_SCORE_JOIN',structured_grch38_key=key)
    return dict(status='GRCH38_GENOMIC_REFERENCE_MATCH_TRANSCRIPT_VERSION_UNRESOLVED',
        variant_key='GRCh38:chr%s:%d:%s:%s'%(chrom,pos,ref,alt),raw_id=hit['_id'],
        note='Transcript/protein arrays preserved raw; no flattened score or patient interpretation enabled')


def run(cfg,root,execute=False):
    assert cfg['patient_variant_inputs_allowed'] is False and cfg['base_url']=='https://myvariant.info/v1' and cfg['assembly']=='hg38'
    assert root.is_absolute() and root.name==cfg['bundle_id'] and not root.is_symlink()
    assert cfg['max_raw_bytes']<=1000000000 and cfg['max_records']<=1000000 and cfg['max_minutes']<=30
    status=root/'status.json'
    if not execute:return load(status) or dict(status='NOT_STARTED_PUBLIC_GENE_REFERENCE_ONLY')
    os.umask(0o077);root.mkdir(mode=0o700,parents=True,exist_ok=True)
    with (root/'.prefetch.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return dict(status='ALREADY_RUNNING')
        deadline=min(dt.datetime.fromisoformat(cfg['deadline_utc']),now()+dt.timedelta(minutes=cfg['max_minutes']))
        panel_genes=public_genes(cfg);genes=[g for g in panel_genes if not cfg.get('only_gene') or g['symbol']==cfg['only_gene']]
        assert genes,'Requested gene must belong to the fixed public source panel'
        save(root/'public-gene-panel.json',dict(source_sha256=cfg['clingen_csv_sha256'],disease_id=cfg['disease_id'],curation_rows=40,genes=panel_genes,download_scope=genes,no_patient_selection=True))
        raw_bytes=sum(p.stat().st_size for p in root.glob('*/page-*.json'));records=0;completed=[];blocked=[]
        def fetch(path,params=None,limit=12000000):
            nonlocal raw_bytes
            assert now()<deadline,'Original budget/runtime deadline reached'
            url=cfg['base_url']+'/'+path+('?' + urllib.parse.urlencode(params) if params else '')
            with urllib.request.urlopen(url,timeout=30) as r:
                assert r.url.startswith(cfg['base_url']+'/'),'Unreviewed redirect'
                data=r.read(min(limit,cfg['max_raw_bytes']-raw_bytes)+1)
            if len(data)>limit or raw_bytes+len(data)>cfg['max_raw_bytes']:raise ValueError('Bounded public cache byte ceiling reached')
            raw_bytes+=len(data);return data,json.loads(data),url
        metadata_raw,metadata,_=fetch('metadata',dict(assembly='hg38'))
        assert metadata.get('build_version')==cfg['expected_api_build'] and metadata.get('src',{}).get('dbnsfp',{}).get('version')==cfg['expected_dbnsfp_version'],'API build changed; review new snapshot'
        (root/'api-metadata.json').write_bytes(metadata_raw)
        fields_raw,fields,_=fetch('metadata/fields',dict(assembly='hg38'));(root/'api-fields.json').write_bytes(fields_raw)
        reference=Reference(cfg['reference_fasta'])
        # Public gene control before any bulk prefetch. No patient coordinate/ALT is read.
        control_raw,control,control_url=fetch('query',dict(q='dbnsfp.genename:MYBPC3',assembly='hg38',size=1,fields=cfg['fields']))
        assert len(control.get('hits',[]))==1
        validation=validate_hit(control['hits'][0],'MYBPC3',reference)
        assert validation['status'].startswith('GRCH38_GENOMIC_REFERENCE_MATCH'),'Public reference control failed'
        save(root/'public-control-receipt.json',dict(status='PASS_PUBLIC_REFERENCE_CONTROL_ONLY',query_url=control_url,validation=validation,raw_sha256=hashlib.sha256(control_raw).hexdigest(),api_build=metadata['build_version'],dbnsfp_version=cfg['expected_dbnsfp_version']))
        (root/'public-control-raw.json').write_bytes(control_raw)
        for gene in genes:
            folder=root/gene['symbol'];receipt=load(folder/'gene-receipt.json')
            if receipt.get('status')=='COMPLETE_PUBLIC_GENE_SOURCE_NOT_ENABLED':
                assert receipt['api_build']==metadata['build_version']
                completed.append(gene['symbol']);records+=receipt['records'];continue
            if records>=cfg['max_records'] or raw_bytes>=cfg['max_raw_bytes'] or now()>=deadline:
                blocked.append(dict(gene=gene['symbol'],reason='BOUNDED_RECORD_BYTE_OR_TIME_CAP'));continue
            folder.mkdir(mode=0o700,exist_ok=True)
            # An interrupted scroll restarts this public gene only. Keep old partial pages as evidence,
            # but never mix duplicate pages into a COMPLETE receipt.
            run_id=now().strftime('%Y%m%dT%H%M%S%f');page=0;count=0;seen=set();pages=[];validation_counts={};expected=None
            params=dict(q='dbnsfp.genename:'+gene['symbol'],assembly='hg38',fields=cfg['fields'],fetch_all='true',size=1000)
            try:
                while True:
                    data,response,url=fetch('query',params)
                    hits=response.get('hits',[])
                    if expected is None:
                        expected=response.get('total');assert isinstance(expected,int),'Unexpected API total'
                        if records+expected>cfg['max_records']:raise ValueError('Gene exceeds remaining bounded record ceiling')
                    if not hits:break
                    path=folder/('page-%s-%04d.json'%(run_id,page));path.write_bytes(data)
                    pages.append(dict(path=path.name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
                    with (folder/('validation-%s.jsonl'%run_id)).open('a') as out:
                        for hit in hits:
                            if hit['_id'] in seen:raise ValueError('Duplicate public source record in scroll')
                            seen.add(hit['_id']);value=validate_hit(hit,gene['symbol'],reference);validation_counts[value['status']]=validation_counts.get(value['status'],0)+1
                            out.write(json.dumps(value)+'\n')
                    count+=len(hits);page+=1
                    save(status,dict(status='PREFETCHING_PUBLIC_GENES_NOT_PATIENT_ANNOTATION',pid=os.getpid(),gene=gene['symbol'],gene_records=count,gene_total=expected,completed_genes=len(completed),total_genes=len(genes),raw_bytes=raw_bytes,total_records=records+count,updated_utc=now().isoformat()))
                    if count>=expected:break
                    scroll=response.get('_scroll_id');assert scroll,'API did not return complete scroll'
                    params=dict(scroll_id=scroll,assembly='hg38');time.sleep(0.5)
                assert count==expected,'Public gene scroll incomplete'
                save(folder/'gene-receipt.json',dict(status='COMPLETE_PUBLIC_GENE_SOURCE_NOT_ENABLED',gene=gene,records=count,api_build=metadata['build_version'],dbnsfp_version=cfg['expected_dbnsfp_version'],pages=pages,validation_counts=validation_counts,transcript_scores_enabled=False,source_specific_terms_retained=True))
                completed.append(gene['symbol']);records+=count
            except (OSError,ValueError,AssertionError) as error:
                blocked.append(dict(gene=gene['symbol'],reason=str(error)[:200]));save(folder/'incomplete-status.json',dict(status='INCOMPLETE_NOT_ENABLED',reason=str(error)[:200]))
            time.sleep(0.5)
        result=dict(status='COMPLETE_PUBLIC_PANEL_SOURCE_NOT_ENABLED' if not blocked else 'PARTIAL_PUBLIC_PANEL_SOURCE_NOT_ENABLED',
            completed_genes=completed,blocked_genes=blocked,total_genes=len(genes),records=records,raw_bytes=raw_bytes,
            api_build=metadata['build_version'],dbnsfp_version=cfg['expected_dbnsfp_version'],cadd_subscore_version='NOT_SPECIFIED_BY_API_METADATA',
            patient_variant_queries=False,patient_scores_enabled=False,annotation_binding_updated=False,updated_utc=now().isoformat())
        save(root/'panel-receipt.json',result);save(status,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=CODE/'config/myvariant-public-hcm-v2.json');p.add_argument('--root',type=Path,default=Path('/srv/marina/resources/myvariant-public-hcm-v1'));p.add_argument('--execute',action='store_true');a=p.parse_args();print(json.dumps(run(load(a.config),a.root,a.execute),indent=2))
