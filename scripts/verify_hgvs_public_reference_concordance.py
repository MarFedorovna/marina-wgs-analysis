#!/usr/bin/env python3
"""Check the two public HCM control transcript IDs against MANE1.4/GENCODE49.

No patient records, no network and no genomic-to-transcript variant conversion.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hgvs-root',type=Path,required=True)
    parser.add_argument('--gencode-gtf',type=Path,required=True)
    args=parser.parse_args()
    manifest=json.loads((args.hgvs_root/'public-source-manifest.json').read_text())
    pairs=manifest.get('mane14_public_pairs',[])
    wanted={p['Ensembl_nuc']:p for p in pairs}
    observed={}
    with gzip.open(args.gencode_gtf,'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            fields=line.rstrip().split('\t')
            if len(fields)!=9 or fields[2]!='transcript':continue
            attrs=dict(re.findall(r'(\w+) "([^"]+)"',fields[8]))
            tid=attrs.get('transcript_id')
            if tid in wanted:
                observed[tid]=dict(gene_name=attrs.get('gene_name'),chromosome=fields[0],strand=fields[6],start=int(fields[3]),end=int(fields[4]))
                if len(observed)==len(wanted):break
    controls=json.loads((args.hgvs_root/'public-control-receipt.json').read_text())
    rows=[]
    for pair in pairs:
        match=observed.get(pair['Ensembl_nuc'])
        rows.append(dict(**pair,gencode49_transcript=match,
                         identity_matches=bool(match and match['gene_name']==pair['symbol'] and match['strand']==pair['chr_strand'])))
    passed=len(rows)==2 and all(row['identity_matches'] for row in rows) and controls['status']=='PASS'
    result=dict(schema_version=1,status='SELECTED_PUBLIC_ID_CONCORDANCE_PASS' if passed else 'INCOMPLETE_OR_MISMATCH',
                mane_release='1.4',gencode_release=49,expected_vep_cache_release=115,
                hgvs_controls_pass=controls['status']=='PASS',public_pairs=rows,
                gencode_file_sha256=hashlib.sha256(args.gencode_gtf.read_bytes()).hexdigest(),
                genomic_to_transcript_variant_mapping_validated=False,
                vep115_runtime_control_validated=False,all_hcm_transcripts_validated=False,
                network_queries=False,patient_data=False)
    (args.hgvs_root/'public-id-concordance-receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
    raise SystemExit(0 if passed else 1)


if __name__=='__main__':main()
