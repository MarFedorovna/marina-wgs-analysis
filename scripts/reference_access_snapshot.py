#!/usr/bin/env python3
"""Read-only snapshot of this isolated public-reference delta; no secrets."""
import datetime as dt
import json
import os
from pathlib import Path

BASE=Path('/srv/marina/resources')
FILES={
 'genereviews':'genereviews-source-only-v1/research-reference-manifest.json',
 'spliceai_wheel':'spliceai-source-only-v1/resource-receipt.json',
 'spliceai_runtime':'spliceai-research-runtime-v1/status.json',
 'spliceai_synthetic_control':'spliceai-research-runtime-v1/synthetic-control-receipt.json',
 'refseq_full_cache_source_queue':'refseq-vep115-source-only-v1/status.json',
 'lovd_public_myh7_gene':'reference-api-public-catalogs-v1/lovd_myh7_gene-receipt.json',
 'omim_public_identifier_map':'reference-api-public-catalogs-v1/omim_mim2gene-receipt.json',
}
FIELDS=['status','pid','bytes','source_bytes','sha256','source_sha256','received_bytes','expected_bytes','updated_utc','retrieved_utc','validated_utc','deadline_utc','error','models_loaded','probabilities_center','reference_alignment','reference_mismatch_rejected','default_annotation','coordinate_review','clinical_pipeline_enabled']
def main():
 result={}
 for name,relative in FILES.items():
  path=BASE/relative
  if not path.is_file():result[name]=dict(status='NOT_YET_AVAILABLE',receipt=str(path));continue
  data=json.loads(path.read_text());item={key:data[key] for key in FIELDS if key in data};item['receipt']=str(path)
  if 'pid' in data:
   try:os.kill(data['pid'],0);item['pid_currently_exists']=True
   except (ProcessLookupError,PermissionError):item['pid_currently_exists']=False
  result[name]=item
 result['omim_full_content']=dict(status='ACCESS_REQUIRED',secret_env='OMIM_API_KEY',enabled=False)
 result['varsome_api']=dict(status='ACCESS_REQUIRED',secret_env='VARSOME_API_KEY',enabled=False)
 print(json.dumps(dict(snapshot_utc=dt.datetime.now(dt.timezone.utc).isoformat(),resources=result,patient_variant_requests=False),indent=2))
if __name__=='__main__':main()
