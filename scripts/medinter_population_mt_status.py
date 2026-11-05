#!/usr/bin/env python3
"""Read-only snapshot of this MedInter public-source delta; no patient outputs or secrets."""
import json
import os
from pathlib import Path
from medinter_population_mt_sources import now,load

BASE=Path('/srv/marina/resources')
FILES={
 'gnomad_public39_counts':'medinter-population-hcm39-v3/status.json',
 'gnomad_public39_counts_receipt':'medinter-population-hcm39-v3/resource-receipt.json',
 'gnomad_exome_AN_source':'medinter-population-mt-sources-v1/gnomad_exomes_all_sites_an_v4_1-status.json',
 'gnomad_exome_AN_subset':'medinter-population-mt-sidecars-v1/an-status.json',
 'gnomad_exome_AN_subset_receipt':'medinter-population-mt-sidecars-v1/an-resource-receipt.json',
 'gnomad_gene_coverage':'medinter-population-gene-coverage-v1/status.json',
 'mt_population':'medinter-population-mt-sidecars-v1/mt-resource-receipt.json',
 'clinvar_dossier':'medinter-population-mt-sidecars-v1/clinvar-status.json',
 'clinvar_dossier_receipt':'medinter-population-mt-sidecars-v1/clinvar-dossier-resource-receipt.json',
 'mitomap_access':'medinter-population-mt-metadata-v1/mitomap_resources_access-receipt.json',
}
FIELDS=['status','pid','updated_utc','completed_utc','bytes','actual_http_bytes','max_http_bytes','completed_intervals','total_intervals','rows','counts','SCV_rows','SCV_rows_without_version','citation_records','variation_ids','completed_genes','downloaded_genes','download_bytes','error','http_status','not_applicable','depth_coverage_downloaded','genome_all_sites_AN_downloaded']
def main():
 result={}
 for name,relative in FILES.items():
  path=BASE/relative;raw=load(path)
  if not raw:result[name]={'status':'NOT_YET_AVAILABLE','path':str(path)};continue
  item={key:raw[key] for key in FIELDS if key in raw};item['path']=str(path)
  if 'pid' in item:
   try:os.kill(item['pid'],0);item['pid_exists_now']=True
   except (ProcessLookupError,PermissionError):item['pid_exists_now']=False
  if name=='gnomad_exome_AN_subset' and item.get('pid_exists_now'):
   directory=Path('/proc')/str(item['pid'])/'fd'
   for descriptor in directory.iterdir():
    try:
     target=descriptor.resolve()
     if target.name!='gnomad_exomes_all_sites_an_v4_1__gnomad.exomes.v4.1.allele_number_all_sites.tsv.bgz':continue
     details=Path('/proc')/str(item['pid'])/'fdinfo'/descriptor.name
     offset=next(int(line.split()[1]) for line in details.read_text().splitlines() if line.startswith('pos:'))
     item.update(compressed_source_bytes_read=offset,compressed_source_bytes_total=target.stat().st_size)
    except (FileNotFoundError,PermissionError,StopIteration):continue
  result[name]=item
 print(json.dumps(dict(snapshot_utc=now().isoformat(),resources=result,annotation_binding_updated=False,patient_variant_inputs=False),indent=2))
if __name__=='__main__':main()
