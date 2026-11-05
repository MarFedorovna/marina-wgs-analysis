#!/usr/bin/env python3
"""Merge source-only public reference candidates. Never updates the sealed queue."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from prefetch_myvariant_public_hcm import load,save


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def merge_cadd(paths,destination):
    seen={};duplicates=0
    for path in paths:
        with path.open() as stream:
            rows=csv.DictReader(stream,delimiter='\t')
            assert rows.fieldnames==['chrom','pos','ref','alt','CADD_RAW','CADD_PHRED']
            for row in rows:
                key=row['chrom'],int(row['pos']),row['ref'],row['alt'];value=row['CADD_RAW'],row['CADD_PHRED']
                if key in seen:
                    assert seen[key]==value,'Conflicting official CADD score for identical allele'
                    duplicates+=1
                seen[key]=value
    with destination.open('x',newline='') as out:
        writer=csv.writer(out,delimiter='\t');writer.writerow(['chrom','pos','ref','alt','CADD_RAW','CADD_PHRED'])
        for (chrom,pos,ref,alt),(raw,phred) in sorted(seen.items()):writer.writerow([chrom,pos,ref,alt,raw,phred])
    return dict(unique_alleles=len(seen),exact_duplicate_rows=duplicates)


def build(bridge,cadd_sources,target,myvariant=None,revel=None):
    assert not target.exists(),'Choose a new candidate version'
    bridge_overlay=load(bridge/'config-overlay.json');assert bridge_overlay['status']=='ROOT_REVIEW_CANDIDATE_NOT_APPLIED'
    source_receipts=[];paths=[]
    for source in cadd_sources:
        receipt=load(source/'source-receipt.json');assert receipt.get('status')=='COMPLETE_PUBLIC_HCM_CDS_SOURCE_ADAPTER_CANDIDATE'
        item=receipt['adapter'];path=Path(item['path'])
        assert path.stat().st_size==item['bytes'] and sha(path)==item['sha256']
        assert receipt['source']['version']=='CADD GRCh38 v1.7';paths.append(path)
        source_receipts.append(dict(path=str(source/'source-receipt.json'),sha256=sha(source/'source-receipt.json')))
    target.mkdir(mode=0o700,parents=True);tsv=target/'cadd-public-hcm-union-v2.tsv';statistics=merge_cadd(paths,tsv)
    adapter=dict(name='cadd_public_hcm_union',status='available',local_adapter=True,format='variant_evidence_tsv',
        version='CADD_GRCh38_v1.7_publicHCM39_CDS_SNV_union_v2',path=str(tsv),bytes=tsv.stat().st_size,sha256=sha(tsv),
        license_status='User-confirmed non-commercial research; official CADD attribution retained',source_url='https://cadd.bihealth.org/download',
        limitation='Public fixed HCM nuclear CDS SNVs only; MT-TI, noncoding, indels and other genes remain unassessed')
    adapters=[*bridge_overlay['annotation_config_append_local_adapters'],adapter]
    additions=dict(bridge_overlay['prepared_resource_manifest_additions'])
    additions['cadd']=dict(status='available_local_supplement',version=adapter['version'],path=str(tsv),bytes=adapter['bytes'],sha256=adapter['sha256'],local_adapter=True,scope=adapter['limitation'])
    pending=[]
    for name,path in [('myvariant',myvariant),('revel',revel)]:
        manifest=load(path/'candidate-manifest.json') if path else {}
        expected='PUBLIC_EXACT_ALLELE_RAW_EVIDENCE_ADAPTER_CANDIDATE' if name=='myvariant' else 'PUBLIC_REVEL_ALLELE_REFERENCE_ADAPTER_CANDIDATE'
        if manifest.get('status')!=expected:pending.append(name);continue
        item=manifest['adapter'];file=Path(item['path']);assert file.stat().st_size==item['bytes'] and sha(file)==item['sha256']
        adapters.append(item)
        keys=['dbnsfp','fathmm','metalr','varity'] if name=='myvariant' else ['revel']
        for key in keys:additions[key]=dict(status='available_local_supplement',version=item['version'],path=str(file),bytes=item['bytes'],sha256=item['sha256'],local_adapter=True,scope='Fixed public HCM panel raw source evidence; canonical per-transcript mapping not enabled; see candidate manifest',source_candidate_manifest=str(path/'candidate-manifest.json'))
        source_receipts.append(dict(path=str(path/'candidate-manifest.json'),sha256=sha(path/'candidate-manifest.json')))
    overlay=dict(status='ROOT_REVIEW_CANDIDATE_NOT_APPLIED',annotation_config_append_local_adapters=adapters,
        prepared_resource_manifest_additions=additions,annotation_binding_updated=False,
        clinical_classification_automatic=False,pending_reference_adapters=pending,
        caution='Append adapters to a new reviewed config/manifest and create a new ROOT binding; old sealed runtime/approval remain untouched. Gene bridge preserves missing TNNT2 HCM validity; raw predictor arrays do not become per-transcript scalar scores.')
    save(target/'config-overlay.json',overlay)
    result=dict(status='PUBLIC_HCM_UNION_V2_ROOT_REVIEW_CANDIDATE',public_gene_count=39,source_receipts=source_receipts,
        cadd=statistics,bridge_manifest=dict(path=str(bridge/'integration-candidate-manifest.json'),sha256=sha(bridge/'integration-candidate-manifest.json')),
        adapter_count=len(adapters),pending_adapters=pending,annotation_binding_updated=False,no_patient_variants_read=True,
        files=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(target.iterdir()) if p.is_file()])
    save(target/'candidate-manifest.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bridge',type=Path,default=Path('/srv/marina/private/publicsupplements-integration-candidate-v2'));p.add_argument('--cadd-sources',type=Path,nargs='+',default=[Path('/srv/marina/resources/cadd-public-hcm-v1'),Path('/srv/marina/resources/cadd-public-hcm-tnnt2-v1')]);p.add_argument('--myvariant',type=Path,default=Path('/srv/marina/private/myvariant-public-hcm-adapter-candidate-v1'));p.add_argument('--revel',type=Path,default=Path('/srv/marina/resources/revel-public-hcm-v1'));p.add_argument('--output',type=Path,default=Path('/srv/marina/private/public-hcm-union-integration-candidate-v2'));a=p.parse_args();print(json.dumps(build(a.bridge,a.cadd_sources,a.output,a.myvariant,a.revel),indent=2))
