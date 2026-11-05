"""Synthetic integration controls for one row-preserving applied-evidence master."""
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import join_combined_annotation_master as joiner
import continue_combined_annotation_master as controller

class CombinedMasterTests(unittest.TestCase):
 def table(self,path,fields,rows):
  path.parent.mkdir(parents=True,exist_ok=True)
  with joiner.text(path,'xt') as stream:
   writer=csv.DictWriter(stream,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
 def fixture(self,root):
  workflow=root/'workflow';tissue=root/'tissue';population=root/'population';output=root/'combined-annotation-master-v1';state=root/'private';state.mkdir();tissue.mkdir();population.mkdir()
  fields=['sample','chrom','pos','ref','alt','transcript','hgnc_id','GT','FT','AD','DP','GQ','candidate_eligible','quality_rank','clinvar_allele_id','note']
  a=dict(zip(fields,['CONTROL_A','chr1','11','A','C','ENST1.1','1','./.','.','3,2','5','.','0','9','111','embedded\tTAB and\nnewline']))
  b=dict(a,sample='CONTROL_B',transcript='ENST1.2',GT='0/1',FT='PASS',candidate_eligible='1',quality_rank='1',note='B')
  c=dict(a,chrom='chr2',pos='22',ref='G',alt='T',hgnc_id='999',clinvar_allele_id='',note='outside scope')
  base_rows=[a,c,b,a.copy()];transcript_rows=[b,a,a.copy(),c]
  master=workflow/'stages/07_report/master_evidence.tsv.gz';transcripts=workflow/'stages/06_annotate/transcript_samples.tsv.gz';self.table(master,fields,base_rows);self.table(transcripts,fields,transcript_rows)
  stage_pins=[]
  for name,path in [('06_annotate',transcripts),('07_report',master)]:
   receipt=path.parent/'workflow-receipt.json';value=joiner.pin(path);receipt.write_text(json.dumps(dict(status='PASS',identity='1'*64,artifacts={path.name:{k:value[k] for k in ['bytes','sha256']}})));stage_pins.append(joiner.pin(receipt))
  context=root/'context.json';context.write_text('{}');context_pin=joiner.pin(context)
  tissue_rows=[]
  for row in [b,a,a.copy()]:
   extra={f:'' for f in joiner.TISSUE_FIELDS};extra.update(reference_context_hgnc='HGNC:1',reference_context_scope='GENE_CONTEXT',gtex_v11_gene_heart_medians_json='[{"median":3.2}]',gtex_v11_transcript_match_status='EXACT_TRANSCRIPT_VERSION_MATCH',pext_source_transcript_match_status='TRANSCRIPT_VERSION_MISMATCH',reference_context_index_sha256=context_pin['sha256']);tissue_rows.append(dict(row,**extra))
  tissue_data=tissue/'public39-reference-context.tsv.gz';self.table(tissue_data,fields+joiner.TISSUE_FIELDS,tissue_rows)
  (tissue/'manifest.json').write_text(json.dumps(dict(status='PASS_SEPARATE_RESEARCH_REFERENCE_CONTEXT_EXPORT',source_sha256=joiner.sha(transcripts),network_requests=0,original_pipeline_changed=False,stage_receipts={p['path']:{'sha256':p['sha256']} for p in stage_pins},context_index_sha256=context_pin['sha256'],artifact=joiner.pin(tissue_data))))
  key='chr1:11:A:C';pop_extra=dict(pop_context_key=key,pop_gnomad_exact_match_status='EXACT_GRCh38_ALLELES',pop_exome_AN_v41='100',pop_exome_AN_status='EXACT_POSITION_EXOME_AN_ONLY',pop_SCV_map_status='EXACT_ALLELES_AND_EXPLICIT_ALLELEID_TO_VARIATIONID',pop_validated_AlleleIDs='111',pop_validated_VariationIDs='9001',pop_gene_coverage_hgnc='HGNC:1',pop_gene_coverage_scope='GENE_BIN_CONTEXT_ONLY_NOT_LOCUS_CALLABILITY')
  self.table(population/'cohort/nuclear-sample-rows.tsv.gz',fields+list(pop_extra),[dict(row,**pop_extra) for row in [a,b,a.copy()]])
  pop_fields=['variant_key','chrom','pos','ref','alt','dataset','source_version','match_status','AC','AN','AF','nhomalt','population_FILTER']
  common=dict(variant_key=key,chrom='chr1',pos='11',ref='A',alt='C',source_version='4.1.1 public HCM39 CDS')
  self.table(population/'cohort/nuclear-population-long.tsv.gz',pop_fields,[dict(common,dataset='exomes',match_status='EXACT_ALLELE',AC='1',AN='100',AF='0.01',nhomalt='0',population_FILTER='AS_VQSR'),dict(common,dataset='genomes',match_status='MISSING_NOT_AF_ZERO',AC='',AN='',AF='',nhomalt='',population_FILTER='')])
  self.table(population/'cohort/nuclear-exome-AN-v41.tsv.gz',['variant_key','chrom','pos','AN','match_status','scope'],[dict(variant_key=key,chrom='chr1',pos='11',AN='100',match_status='EXACT_POSITION_EXOME_AN_ONLY',scope='POSITION_CONTEXT_v4.1_EXOMES_ONLY_NOT_GENOMES_NOT_DEPTH')])
  self.table(population/'cohort/ClinVar-validated-allele-map.tsv.gz',['variant_key','AlleleID','VariationID','scope'],[dict(variant_key=key,AlleleID='111',VariationID='9001',scope='EXACT_GRCh38_ALLELES_SOURCE_ID_MAP_NOT_CLINICAL_CLASSIFICATION')])
  self.table(population/'cohort/ClinVar-SCV-long.tsv.gz',['variant_key','VariationID','SCV_accession_version','source_classification'],[dict(variant_key=key,VariationID='9001',SCV_accession_version='SCV000001.1',source_classification='Pathogenic'),dict(variant_key=key,VariationID='9001',SCV_accession_version='SCV000002.2',source_classification='Benign')])
  self.table(population/'cohort/gene-coverage-context-long.tsv.gz',['hgnc_id','source_cohort','scope'],[dict(hgnc_id='HGNC:1',source_cohort='exome',scope='GENE_BIN_CONTEXT_ONLY_NO_EXACT_LOCUS_DEPTH_CLAIM')])
  self.table(population/'cohort/ClinVar-citation-links.tsv.gz',['variant_key','PMID'],[dict(variant_key=key,PMID='1')])
  self.table(population/'cohort/MT-population-reference-long.tsv.gz',['sample','REF','ALT','GT'],[dict(sample='CONTROL_A',REF='A',ALT='C,G',GT='1/2')])
  ready=root/'ready.json';ready.write_text(json.dumps({'source_versions':{'counts':'4.1.1','AN':'4.1 exomes','coverage':'gene bins only'}}));ready_pin=joiner.pin(ready)
  artifacts={p.relative_to(population).as_posix():{k:joiner.pin(p)[k] for k in ['bytes','sha256']} for p in population.rglob('*') if p.is_file()}
  provenance=dict(status='PASS_LOCAL_POPULATION_MT_CONTEXT_EXPORT',source_master=joiner.pin(master),primary_receipts=stage_pins,resource_ready_manifest=ready_pin,network_requests=0,primary_outputs_changed=False)
  (population/'workflow-receipt.json').write_text(json.dumps(dict(status='PASS',stage='population-reference-context-v1',identity='2'*64,metadata=provenance,artifacts=artifacts)))
  panel=root/'panel.json';panel.write_text(json.dumps({'genes':[{'hgnc_id':'HGNC:1'}]}))
  rule=dict(workflow_root=str(workflow),tissue_root=str(tissue),population_root=str(population),tissue_context_index=context_pin,resource_ready_manifest=ready_pin,public_panel=joiner.pin(panel),output=str(output),state_directory=str(state),deadline_utc=joiner.DEADLINE.isoformat(),memory_limit_bytes=4*1024**3,cpu_limit=2,max_sqlite_bytes=8*1024**3,max_output_bytes=6*1024**3,min_free_bytes=0,code_sha256={'synthetic':'3'*64})
  return rule,fields,base_rows
 def test_full_combined_master_preserves_every_row_order_field_and_actual_context(self):
  with tempfile.TemporaryDirectory() as directory:
   rule,fields,expected=self.fixture(Path(directory))
   with patch.object(joiner,'limits',return_value=None):result=joiner.join(rule)
   output=Path(rule['output'])
   with joiner.text(output/'master_evidence.enriched.tsv.gz') as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
   self.assertEqual([{f:row[f] for f in fields} for row in rows],expected);self.assertEqual(result['row_count'],4);self.assertEqual(result['context_rows_matched'],{'tissue':3,'population':3});self.assertEqual(result['orphan_context_rows'],{'tissue':0,'population':0})
   digest=hashlib.sha256(joiner.canonical_tsv(fields)+b''.join(joiner.canonical_tsv([row[f] for f in fields]) for row in expected)).hexdigest();self.assertEqual(digest,result['invariant_digest'])
   self.assertEqual(json.loads(rows[0]['allref_gnomad_genomes_summary_json'])['AF'],'');self.assertEqual(rows[0]['allref_SCV_reference_count'],'2');self.assertIn('SCV_CLASSIFICATIONS_DIFFER',rows[0]['allref_review_flags']);self.assertEqual(rows[0]['allref_tissue_gtex_v11_gene_heart_medians_json'],'[{"median":3.2}]');self.assertEqual(rows[1]['allref_population_status'],'OUTSIDE_PUBLIC_SUBSET_OR_UNSUPPORTED_ALLELE_NOT_AF_ZERO')
   self.assertEqual(result['MT_rows_added_to_nuclear_master'],0);self.assertTrue((output/'applied-evidence/MT-appendix-index.json').is_file());self.assertTrue(all(p['relative_path'].startswith('applied-evidence/') and p['source_upstream_receipt'] for p in result['auxiliary_artifacts']))
 def test_changed_original_tissue_field_cannot_join_by_index(self):
  with tempfile.TemporaryDirectory() as directory:
   root=Path(directory);rule,fields,rows=self.fixture(root);db=joiner.database(root/'index.sqlite');source=root/'different.tsv.gz';bad=dict(rows[0],FT='PASS');extra={f:'' for f in joiner.TISSUE_FIELDS};self.table(source,fields+joiner.TISSUE_FIELDS,[dict(bad,**extra)])
   joiner.index_rows(db,'tissue',source,fields,joiner.TISSUE_FIELDS,joiner.DEADLINE)
   self.assertIsNone(joiner.take(db,'tissue',joiner.row_key(rows[0],fields)));db.close()
 def test_AlleleID_VariationID_conflation_rejected_in_combined_overlay(self):
  with tempfile.TemporaryDirectory() as directory:
   root=Path(directory);rule,fields,rows=self.fixture(root);db=joiner.database(root/'index.sqlite');payload={f:'' for f in joiner.TISSUE_FIELDS};joiner.index_rows(db,'tissue',root/'tissue/public39-reference-context.tsv.gz',fields,joiner.TISSUE_FIELDS,joiner.DEADLINE)
   p=root/'population/cohort';joiner.index_long(db,'population',p/'nuclear-population-long.tsv.gz','variant_key');joiner.index_long(db,'AN',p/'nuclear-exome-AN-v41.tsv.gz','variant_key');joiner.index_long(db,'map',p/'ClinVar-validated-allele-map.tsv.gz','variant_key')
   row=rows[0];digest=joiner.row_key(row,fields);pop=dict(pop_context_key='chr1:11:A:C',pop_exome_AN_v41='100',pop_exome_AN_status='EXACT_POSITION_EXOME_AN_ONLY',pop_SCV_map_status='EXACT',pop_validated_AlleleIDs='9001',pop_validated_VariationIDs='111',pop_gene_coverage_hgnc='')
   db.execute('INSERT INTO row_context VALUES (?,?,1,0,?)',('population',digest,json.dumps(pop)))
   with self.assertRaises(ValueError):joiner.extras(db,row,digest,1,{'HGNC:1'},'{}')
   db.close()
 def test_root_gate_and_all_four_upstream_PASS_requirements(self):
  with tempfile.TemporaryDirectory() as directory:
   root=Path(directory);workflow=root/'workflow'
   for name in ['06_annotate','07_report']:
    p=workflow/'stages'/name/'workflow-receipt.json';p.parent.mkdir(parents=True);p.write_text(json.dumps({'status':'PASS','identity':'a'*64}))
   tissue=root/'tissue-state.json';pop=root/'population-state.json';tissue.write_text(json.dumps({'status':'TISSUE_OK','rule_sha256':'b'*64}));pop.write_text(json.dumps({'status':'WAITING'}));rule=dict(workflow_root=str(workflow),producer_gates=[dict(kind='tissue',path=str(tissue),success_status='TISSUE_OK'),dict(kind='population',path=str(pop),success_status='POP_OK')])
   self.assertFalse(controller.readiness(rule));pop.write_text(json.dumps({'status':'POP_OK','rule_sha256':'c'*64}));self.assertTrue(controller.readiness(rule));pop.write_text(json.dumps({'status':'FAILED'}))
   with self.assertRaises(ValueError):controller.readiness(rule)
   path=root/'candidate.json';path.write_text(json.dumps({'review_status':'PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING'}))
   with self.assertRaises(ValueError):controller.validate(load_json(path),path,controller.sha(path))
 def test_CPU_memory_and_original_deadline_are_not_extendable(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'rule.json';base=dict(schema_version=1,review_status='ROOT_REVIEWED_EXACT_COMBINED_MASTER_RULE',no_network=True,no_primary_or_views_mutation=True,cpu_limit=2,memory_limit_bytes=4*1024**3,max_hours=8,interval_seconds=30,max_sqlite_bytes=8*1024**3,max_output_bytes=6*1024**3,min_free_bytes=150*1024**3,deadline_utc=controller.DEADLINE.isoformat())
   for change in [dict(cpu_limit=3),dict(memory_limit_bytes=8*1024**3),dict(deadline_utc='2026-11-07T00:00:00+00:00')]:
    rule=dict(base,**change);path.write_text(json.dumps(rule))
    with self.subTest(change=change),self.assertRaises(ValueError):controller.validate(rule,path,controller.sha(path),dt.datetime(2026,11,5,20,tzinfo=dt.timezone.utc))
def load_json(path):return json.loads(path.read_text())
if __name__=='__main__':unittest.main()
