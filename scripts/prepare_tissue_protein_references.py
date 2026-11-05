#!/usr/bin/env python3
"""Compact official public UniProt/GTEx39 summaries; no patient input or pipeline binding."""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import urllib.parse
import urllib.request

CODE=Path(__file__).resolve().parents[1]


def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return dt.datetime.now(dt.timezone.utc)


def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp');assert not path.is_symlink() and not temporary.is_symlink()
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def panel(cfg):
    path=Path(cfg['panel_manifest']);assert sha(path)==cfg['panel_sha256']
    source=load(path);assert source['no_patient_selection'] is True and len(source['genes'])==cfg['expected_genes']==39
    return source['genes']


def uniprot_gene_ids(entry):
    assert entry.get('organism',{}).get('taxonId')==9606 and entry.get('entryType')=='UniProtKB reviewed (Swiss-Prot)'
    return {x['id'] for x in entry.get('uniProtKBCrossReferences',[]) if x.get('database')=='HGNC'}


def gtex_identity(row,genes):
    hgnc=re.search(r'Acc:(HGNC:\d+)',row.get('description',''))
    return next((g for g in genes if hgnc and g['hgnc_id']==hgnc[1] and row.get('geneSymbol')==g['symbol']),None)


def run(cfg,profile,execute=False):
    root=Path(cfg['root'])/profile
    if not execute:return load(root/'status.json') if (root/'status.json').exists() else dict(status='NOT_STARTED_PUBLIC_REFERENCE_ONLY')
    assert cfg['patient_queries_allowed'] is False and cfg['max_bytes_per_profile']<=100000000 and cfg['max_minutes']<=30
    os.umask(0o077);root.mkdir(mode=0o700,parents=True,exist_ok=True)
    if (root/'resource-manifest.json').exists():return load(root/'resource-manifest.json')
    import fcntl
    with (root/'.prepare.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return dict(status='ALREADY_RUNNING')
        genes=panel(cfg);deadline=min(dt.datetime.fromisoformat(cfg['deadline_utc']),now()+dt.timedelta(minutes=cfg['max_minutes']));receipts=[];total=0
        save(root/'status.json',dict(status='FETCHING_FIXED_PUBLIC_GENE_REFERENCES',pid=os.getpid(),profile=profile,genes=39,updated_utc=now().isoformat()))
        def fetch(url,name):
            nonlocal total
            assert now()<deadline
            assert url.startswith(('https://rest.uniprot.org/uniprotkb/search?','https://gtexportal.org/api/v2/reference/gene?','https://gtexportal.org/api/v2/expression/median','https://gtexportal.org/api/v2/openapi.json'))
            path=root/name;receipt_path=path.with_suffix(path.suffix+'.receipt.json')
            if path.exists():
                receipt=load(receipt_path);assert path.stat().st_size==receipt['bytes'] and sha(path)==receipt['sha256'] and receipt['url']==url
                data=load(path);headers=receipt['selected_headers']
            else:
                with urllib.request.urlopen(url,timeout=25) as response:
                    assert response.url==url,'Unreviewed redirect';raw=response.read(cfg['max_bytes_per_profile']-total+1)
                    headers={k:v for k,v in response.headers.items() if k.lower() in ('x-uniprot-release','x-uniprot-release-date','x-total-results','link','etag','last-modified','date')}
                assert total+len(raw)<=cfg['max_bytes_per_profile'],'Compact source cap reached'
                data=json.loads(raw);path.write_bytes(raw)
                receipt=dict(status='PUBLIC_SOURCE_RECEIVED',url=url,bytes=len(raw),sha256=sha(path),selected_headers=headers,retrieved_utc=now().isoformat());save(receipt_path,receipt)
            total+=path.stat().st_size;receipts.append(receipt);return data,headers
        try:
            if profile=='uniprot':
                query='('+' OR '.join('gene_exact:'+g['symbol'] for g in genes)+') AND (organism_id:9606) AND (reviewed:true)'
                url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode(dict(query=query,format='json',size=500))
                doc,headers=fetch(url,'uniprot-reviewed-human-hcm39.json');h={k.lower():v for k,v in headers.items()}
                assert h.get('x-uniprot-release')==cfg['uniprot_release'] and int(h['x-total-results'])==len(doc['results']) and 'rel="next"' not in h.get('link','')
                mapping={g['hgnc_id']:g for g in genes};entries=[];features=[]
                for entry in doc['results']:
                    matching=uniprot_gene_ids(entry)&mapping.keys()
                    for hgnc in matching:
                        base=dict(hgnc_id=hgnc,gene=mapping[hgnc]['symbol'],accession=entry['primaryAccession'],release=h['x-uniprot-release'],entry_audit=entry.get('entryAudit'),sequence_length=entry.get('sequence',{}).get('length'),sequence=entry.get('sequence',{}).get('value'),source_sha256=receipts[-1]['sha256'])
                        entries.append(base)
                        for feature in entry.get('features',[]):features.append(dict(hgnc_id=hgnc,gene=base['gene'],accession=base['accession'],release=base['release'],feature=feature))
                assert entries
                save(root/'protein-summary.json',dict(status='PUBLIC_REVIEWED_PROTEIN_REFERENCE_ONLY',proteins=entries))
                save(root/'protein-features.json',dict(coordinates='UniProt protein sequence/isoform-specific; not genomic; retain complete location/evidence',features=features))
                available={e['hgnc_id'] for e in entries}
                result=dict(status='AVAILABLE_PUBLIC_UNIPROT39_REFERENCE_NOT_BOUND',release=h['x-uniprot-release'],release_date=h.get('x-uniprot-release-date'),proteins=len(entries),features=len(features),genes_with_reviewed_protein=len(available),missing_genes=[g for g in genes if g['hgnc_id'] not in available],license=cfg['uniprot_license'],license_url=cfg['uniprot_license_url'])
            else:
                fetch('https://gtexportal.org/api/v2/openapi.json','gtex-openapi.json')
                params=[('geneId',g['symbol']) for g in genes]+[('gencodeVersion',cfg['gtex_gencode']),('genomeBuild',cfg['gtex_genome']),('itemsPerPage','100000')]
                doc,_=fetch('https://gtexportal.org/api/v2/reference/gene?'+urllib.parse.urlencode(params),'gtex-reference-genes.json')
                references={r['gencodeId']:(r,gtex_identity(r,genes)) for r in doc['data'] if gtex_identity(r,genes)}
                assert references and all(r['gencodeVersion']==cfg['gtex_gencode'] and r['genomeBuild']==cfg['gtex_genome'] for r,g in references.values())
                tables={}
                for endpoint in ('medianGeneExpression','medianTranscriptExpression','medianExonExpression','medianJunctionExpression'):
                    params=[('gencodeId',g) for g in references]+[('datasetId',cfg['gtex_dataset'])]+[('tissueSiteDetailId',t) for t in cfg['tissues']]+[('itemsPerPage','100000')]
                    data,_=fetch('https://gtexportal.org/api/v2/expression/'+endpoint+'?'+urllib.parse.urlencode(params),endpoint+'.json')
                    assert data.get('paging_info',{}).get('numberOfPages',1)<=1,'Compact endpoint needs pagination review'
                    assert all(r.get('tissueSiteDetailId') in cfg['tissues'] and r.get('datasetId',cfg['gtex_dataset'])==cfg['gtex_dataset'] for r in data['data'])
                    tables[endpoint]=len(data['data'])
                result=dict(status='AVAILABLE_PUBLIC_GTEX_HEART_MEDIANS_NOT_PEXT_NOT_BOUND',dataset=cfg['gtex_dataset'],gencode=cfg['gtex_gencode'],genome=cfg['gtex_genome'],tissues=cfg['tissues'],source_reference_genes=len(references),table_rows=tables,
                    gene_reference_records=[dict(hgnc_id=g['hgnc_id'],gene=g['symbol'],source_gencode_id=r['gencodeId']) for r,g in references.values()],
                    version_limitation='GTEx GENCODE47 is not current workflow GENCODE49; no silent transcript/exon merge or patient phenotype assignment',protected_or_donor_data_requested=False,pext_status='SEPARATE_GNOMAD_SOURCE_REQUIRED')
            result.update(profile=profile,source_receipts=receipts,total_download_bytes=total,public_panel_sha256=cfg['panel_sha256'],no_patient_queries=True,annotation_binding_updated=False,completed_utc=now().isoformat())
            save(root/'resource-manifest.json',result);save(root/'status.json',result);return result
        except BaseException as error:
            save(root/'status.json',dict(status='INCOMPLETE_OR_BLOCKED_NOT_BOUND',profile=profile,error=str(error)[:300],download_bytes=total,source_receipts=receipts,updated_utc=now().isoformat()));raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=CODE/'config/tissue-protein-public39-v1.json');p.add_argument('--profile',choices=['uniprot','gtex'],required=True);p.add_argument('--execute',action='store_true');a=p.parse_args();print(json.dumps(run(load(a.config),a.profile,a.execute),indent=2))
