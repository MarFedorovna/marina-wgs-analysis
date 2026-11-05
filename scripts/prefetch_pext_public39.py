#!/usr/bin/env python3
"""Fixed public HGNC39 pext/gene transcript reference queries; no variant query field."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.error
import urllib.request
from prepare_tissue_protein_references import load,save,panel,now

CODE=Path(__file__).resolve().parents[1]
QUERY='''query PublicGeneReference($symbol:String!) { gene(gene_symbol:$symbol,reference_genome:GRCh38) {
 reference_genome gene_id gene_version symbol hgnc_id chrom start stop strand
 transcripts { transcript_id transcript_version gtex_tissue_expression { tissue value } }
 pext { flags regions { start stop mean tissues { tissue value } } }
} }'''


def validate_gene(value,gene):
    assert value and value['reference_genome']=='GRCh38' and value['symbol']==gene['symbol']
    hgnc=value.get('hgnc_id','');hgnc='HGNC:'+hgnc if str(hgnc).isdigit() else hgnc
    assert hgnc==gene['hgnc_id'],'Exact public HGNC source identity differs'
    return value


def run(cfg,execute=False):
    root=Path(cfg['root'])/'pext'
    if not execute:return load(root/'status.json') if (root/'status.json').exists() else dict(status='NOT_STARTED_PUBLIC_GENE_PEXT_ONLY')
    assert cfg['patient_queries_allowed'] is False
    os.umask(0o077);root.mkdir(mode=0o700,parents=True,exist_ok=True)
    import fcntl
    with (root/'.prepare.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return dict(status='ALREADY_RUNNING')
        genes=panel(cfg);deadline=min(dt.datetime.fromisoformat(cfg['deadline_utc']),now()+dt.timedelta(minutes=30));receipts=[];missing=[];total=0;last=0
        save(root/'status.json',dict(status='FETCHING_PUBLIC_GENE_PEXT_SEQUENTIAL_RATE_LIMITED',pid=os.getpid(),completed_genes=0,total_genes=39,updated_utc=now().isoformat()))
        try:
            for gene in genes:
                path=root/(gene['symbol']+'.json');receipt_path=root/(gene['symbol']+'-receipt.json')
                if path.exists():
                    receipt=load(receipt_path);assert path.stat().st_size==receipt['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest()==receipt['sha256'];doc=load(path)
                else:
                    delay=max(0,7-(time.monotonic()-last))
                    if delay:time.sleep(delay)
                    assert now()<deadline
                    payload=dict(query=QUERY,variables=dict(symbol=gene['symbol']))
                    request=urllib.request.Request('https://gnomad.broadinstitute.org/api/',data=json.dumps(payload).encode(),headers={'Content-Type':'application/json'})
                    with urllib.request.urlopen(request,timeout=25) as response:
                        raw=response.read(cfg['max_bytes_per_profile']-total+1)
                    last=time.monotonic();assert total+len(raw)<=cfg['max_bytes_per_profile'],'Compact pext response ceiling reached'
                    doc=json.loads(raw);assert not doc.get('errors'),'GraphQL returned errors; no source marked complete'
                    receipt=dict(status='PUBLIC_GENE_REFERENCE_RECEIVED',endpoint='https://gnomad.broadinstitute.org/api/',query=payload,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),retrieved_utc=now().isoformat())
                    path.write_bytes(raw);save(receipt_path,receipt)
                value=doc['data'].get('gene')
                if value is not None:validate_gene(value,gene)
                total+=receipt['bytes'];receipts.append(receipt)
                if value is None or not value.get('pext') or not value['pext'].get('regions'):missing.append(dict(gene=gene,reason='PEXT_NOT_AVAILABLE_IN_PUBLIC_GENE_MODEL_NOT_ZERO_OR_BENIGN'))
                save(root/'status.json',dict(status='FETCHING_PUBLIC_GENE_PEXT_SEQUENTIAL_RATE_LIMITED',pid=os.getpid(),completed_genes=len(receipts),total_genes=39,download_bytes=total,last_gene=gene['symbol'],updated_utc=now().isoformat()))
            result=dict(status='AVAILABLE_PUBLIC_PEXT39_REFERENCE_NOT_BOUND',genes=39,download_bytes=total,source_receipts=receipts,missing_pext=missing,
                reference_genome='GRCh38',source_context='gnomAD browser GRCh38 gene models; pext based on GTExv10 per primary release documentation',
                source_release_returned_by_API=False,primary_documentation=['https://gnomad.broadinstitute.org/news/2024-11-gnomad-v4-pext/','https://gnomad.broadinstitute.org/news/2026-03-gnomad-v4-1-1/'],
                coordinate_contract='Preserve API start/stop as returned; gene context usable, allele/exon overlay requires reviewed coordinate/transcript-version contract',
                version_limitation='Browser GTExv10 gene/transcript model is not GTExv11/GENCODE47 or workflow GENCODE49',
                patient_variant_queries=False,annotation_binding_updated=False,no_automatic_PVS1_or_benign_call=True,completed_utc=now().isoformat())
            save(root/'resource-manifest.json',result);save(root/'status.json',result);return result
        except BaseException as error:
            save(root/'status.json',dict(status='INCOMPLETE_OR_BLOCKED_PUBLIC_PEXT_NOT_BOUND',completed_genes=len(receipts),total_genes=39,download_bytes=total,error=str(error)[:300],updated_utc=now().isoformat()));raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=CODE/'config/tissue-protein-public39-v1.json');p.add_argument('--execute',action='store_true');a=p.parse_args();print(json.dumps(run(load(a.config),a.execute),indent=2))
