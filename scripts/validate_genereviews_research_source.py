#!/usr/bin/env python3
"""New research-reference manifest; retain original download receipt and archive."""
import datetime as dt
import json
from pathlib import Path
from download_external_reference_source import safe,save,sha256,validate_archive

ROOT=Path('/srv/marina/resources/genereviews-source-only-v1')
def main():
 safe(ROOT)
 purpose=Path('/srv/marina/resources/spliceai-research-runtime-v1/research-purpose-approval.json')
 approval=json.loads(purpose.read_text())
 if approval.get('purpose')!='noncommercial_research' or approval.get('source')!='explicit_user_confirmation':raise ValueError('Purpose approval missing')
 receipt=ROOT/'resource-receipt.json';source=json.loads(receipt.read_text());archive=ROOT/'gene_NBK1116.tar.gz'
 if archive.stat().st_size!=640489521 or sha256(archive)!=source['sha256']:raise ValueError('Source receipt mismatch')
 inventory=validate_archive(archive,'genereviews')
 manifest=dict(status='AVAILABLE_NONCOMMERCIAL_RESEARCH_REFERENCE_VERBATIM',source=str(archive),source_url=source['url'],source_bytes=archive.stat().st_size,source_sha256=source['sha256'],original_download_receipt_sha256=sha256(receipt),validated_archive_inventory=inventory,archive_unmodified=True,extracted=False,purpose='noncommercial_research',purpose_approval_sha256=sha256(purpose),attribution='GeneReviews® [Internet]. University of Washington, Seattle; 1993–2026. © 1993–2026 University of Washington, Seattle.',attribution_url='https://www.ncbi.nlm.nih.gov/books/NBK1116/',terms_url='https://www.ncbi.nlm.nih.gov/sites/books/NBK138602/',reuse_requirements=['Retain source and copyright attribution with each copy','Link original materials when published on the web','No further modifications to GeneReviews materials'],clinical_variant_annotator=False,automatic_pathogenicity_evidence=False,validated_utc=dt.datetime.now(dt.timezone.utc).isoformat())
 path=ROOT/'research-reference-manifest.json'
 if path.exists():raise ValueError('Research manifest already exists; history preserved')
 save(path,manifest)
 print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
