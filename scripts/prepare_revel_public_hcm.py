#!/usr/bin/env python3
"""Stream the verified REVEL archive into a fixed public CDS reference TSV; no patient input."""
import argparse
import bisect
import csv
import datetime as dt
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import zipfile
from prefetch_myvariant_public_hcm import Reference,load,save

SOURCE_FIELDS=['chr','hg19_pos','grch38_pos','ref','alt','aaref','aaalt','REVEL','Ensembl_transcriptid']
FIELDS=['chrom','pos','ref','alt','transcript','REVEL_reference_score','REVEL_source_transcript_ids',
    'REVEL_source_aa_ref','REVEL_source_aa_alt','REVEL_source_hg19_pos','REVEL_mapping_status','REVEL_version']


def index_bed(path):
    by_chrom={}
    for line in path.read_text().splitlines():
        chrom,start,end=line.split('\t')[:3];by_chrom.setdefault(chrom,[]).append((int(start),int(end)))
    return {c:(sorted(v),[r[0] for r in sorted(v)]) for c,v in by_chrom.items()}


def candidate_row(row,index,reference):
    chrom=row['chr'];value=row['grch38_pos']
    if chrom not in index or not value.isdigit():return None,'NO_LIFTED_POSITION_OR_OUTSIDE_PUBLIC_CHROMOSOMES'
    pos=int(value);ranges,starts=index[chrom];i=bisect.bisect_right(starts,pos-1)-1
    if i<0 or not ranges[i][0]<=pos-1<ranges[i][1]:return None,'OUTSIDE_FIXED_PUBLIC_CDS'
    ref,alt=row['ref'],row['alt']
    if ref not in ('A','C','G','T') or alt not in ('A','C','G','T') or ref==alt:return None,'INVALID_SOURCE_ALLELE'
    if reference.base(chrom,pos)!=ref:return None,'LIFTED_GRCH38_REF_MISMATCH_NOT_ENABLED'
    score=float(row['REVEL']);assert math.isfinite(score) and 0<=score<=1
    transcripts=row['Ensembl_transcriptid']
    # Do not manufacture transcript versions or copy this value into every VEP transcript.
    mapping='EXACT_GRCH38_ALLELE_REFERENCE;SOURCE_TRANSCRIPT_VERSION_UNAVAILABLE;CURRENT_AA_MATCH_REQUIRED'
    if not re.fullmatch(r'ENST\d+(?:;ENST\d+)*',transcripts):mapping+=';SOURCE_TRANSCRIPTS_REQUIRE_REVIEW'
    return dict(chrom='chr'+chrom,pos=str(pos),ref=ref,alt=alt,transcript='*',
        REVEL_reference_score=row['REVEL'],REVEL_source_transcript_ids=transcripts,
        REVEL_source_aa_ref=row['aaref'],REVEL_source_aa_alt=row['aaalt'],REVEL_source_hg19_pos=row['hg19_pos'],
        REVEL_mapping_status=mapping,REVEL_version='1.3'),None


def run(archive_path,panel,target,reference_path,approval):
    assert load(approval).get('status')=='USER_CONFIRMED_NONCOMMERCIAL_RESEARCH'
    receipt=load(archive_path.parent/'resource-receipt.json')
    assert receipt.get('bytes')==667102707 and receipt.get('publisher_md5')=='3ea2bc33e6b5455fc7e9899da863b5fe'
    assert archive_path.stat().st_size==receipt['bytes'] and receipt['source']['version']=='1.3'
    pm=load(panel/'panel-manifest.json');bed=panel/'hcm-cds-union-cadd.bed'
    assert pm['no_patient_selection'] is True and pm['assembly']=='GRCh38' and hashlib.sha256(bed.read_bytes()).hexdigest()==pm['bed_sha256']
    assert not target.exists(),'Choose a new candidate version; preserve any existing partial'
    os.umask(0o077);target.mkdir(mode=0o700,parents=True)
    index=index_bed(bed);reference=Reference(reference_path);deadline=min(dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00'),dt.datetime.now(dt.timezone.utc)+dt.timedelta(minutes=30))
    scanned=kept=out_bytes=0;reasons={};partial=target/'revel-public-hcm.tsv.partial'
    try:
        with zipfile.ZipFile(archive_path) as archive:
            members=archive.infolist();assert len(members)==1 and members[0].filename=='revel_with_transcript_ids' and members[0].file_size==6496023640
            with archive.open(members[0]) as binary,io.TextIOWrapper(binary,encoding='utf-8',newline='') as stream,partial.open('x',newline='') as out:
                reader=csv.DictReader(stream);assert reader.fieldnames==SOURCE_FIELDS,'REVEL source schema changed'
                writer=csv.DictWriter(out,fieldnames=FIELDS,delimiter='\t');writer.writeheader()
                for row in reader:
                    scanned+=1;value,reason=candidate_row(row,index,reference)
                    if value:
                        writer.writerow(value);kept+=1
                    elif reason=='LIFTED_GRCH38_REF_MISMATCH_NOT_ENABLED':reasons[reason]=reasons.get(reason,0)+1
                    if scanned%100000==0:
                        out.flush();out_bytes=partial.stat().st_size
                        assert out_bytes<=300000000,'Derived public reference output ceiling reached'
                        assert dt.datetime.now(dt.timezone.utc)<deadline,'Original budget/runtime cap reached'
                        save(target/'status.json',dict(status='STREAMING_VERIFIED_SOURCE_PUBLIC_CDS_ONLY',pid=os.getpid(),source_records_scanned=scanned,kept_public_records=kept,output_bytes=out_bytes,full_archive_uncompressed_bytes=members[0].file_size,uncompressed_file_not_extracted=True,updated_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
        final=target/'revel-public-hcm.tsv';partial.replace(final)
        adapter=dict(name='revel_hcm_reference',status='available',format='variant_evidence_tsv',local_adapter=True,version='REVEL1.3_publicHCM_CDS_reference_v1',path=str(final),bytes=final.stat().st_size,sha256=hashlib.sha256(final.read_bytes()).hexdigest(),license_status='User confirmed non-commercial research; REVEL author terms/attribution retained',source_url='https://sites.google.com/site/revelgenomics/downloads')
        result=dict(status='PUBLIC_REVEL_ALLELE_REFERENCE_ADAPTER_CANDIDATE',adapter=adapter,source_archive_sha256=receipt['sha256'],publisher_md5=receipt['publisher_md5'],public_panel_sha256=hashlib.sha256((panel/'panel-manifest.json').read_bytes()).hexdigest(),source_records_scanned=scanned,kept_public_records=kept,excluded_reasons=reasons,
            source_transcript_version_available=False,canonical_per_transcript_REVEL_enabled=False,
            raw_fields_retained_for_exact_current_transcript_base_and_amino_acid_review=True,
            historical_hg19_position_not_used_as_GRCh38=True,archive_crc_validated_by_complete_stream=True,
            no_patient_variants_processed=True,annotation_binding_updated=False,automatic_acmg=False)
        save(target/'candidate-manifest.json',result);save(target/'status.json',result);return result
    except BaseException as error:
        save(target/'status.json',dict(status='INCOMPLETE_PUBLIC_REFERENCE_NOT_ENABLED',source_records_scanned=scanned,kept_public_records=kept,error=str(error)[:200]));raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,default=Path('/srv/marina/resources/revel-source-only-v1/revel-v1.3_all_chromosomes.zip'));p.add_argument('--panel',type=Path,default=Path('/srv/marina/resources/public-hcm-cds-v2'));p.add_argument('--output',type=Path,default=Path('/srv/marina/resources/revel-public-hcm-v1'));p.add_argument('--reference',type=Path,default=Path('/srv/marina/references/broad-hg38-v0-altaware-g1575676516681666-prebuilt/Homo_sapiens_assembly38.fasta'));p.add_argument('--approval',type=Path,default=Path('/srv/marina/private/resource-use-approval-20261105-noncommercial.json'));a=p.parse_args();print(json.dumps(run(a.archive,a.panel,a.output,a.reference,a.approval),indent=2))
