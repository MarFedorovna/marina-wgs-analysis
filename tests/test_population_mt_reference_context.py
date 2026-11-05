"""Synthetic-only local join, sample quality and review-gate controls."""
from collections import defaultdict
from contextlib import ExitStack
import csv
import datetime as dt
import gzip
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import export_population_mt_reference_context as exporter
import continue_population_mt_reference_context as controller

class LocalContextTests(unittest.TestCase):
 def source(self):
  source=exporter.Sources.__new__(exporter.Sources);source.cv_keys=defaultdict(set,{('chr1',11,'A','C'):{('111','9001')}});source.pop_keys={('chr1',11,'A','C'):[1]}
  source.an=sqlite3.connect(':memory:');source.an.row_factory=sqlite3.Row;source.an.execute('CREATE TABLE all_sites_AN(chrom,pos,AN,raw_json)');source.an.execute('INSERT INTO all_sites_AN VALUES (?,?,?,?)',('chr1',11,100,'{"capture_flag":false}'))
  source.cv=sqlite3.connect(':memory:');source.cv.row_factory=sqlite3.Row
  source.cv.executescript('CREATE TABLE submissions(variation_id,SCV,SCV_version_present,submitter,date_evaluated,classification,raw_json);CREATE TABLE citations(variation_id,citation_source,citation_id);')
  source.cv.execute('INSERT INTO submissions VALUES (?,?,?,?,?,?,?)',('9001','SCV000001.2',1,'PUBLIC_CONTROL','2020-01-01','Uncertain significance','{}'));source.cv.execute('INSERT INTO citations VALUES (?,?,?)',('9001','PubMed','1'))
  source.counts=lambda key:[dict(source='exomes-chr1',AC='0',AN=100,AF='0',nhomalt='0',filter='AC0',info='AC=0;AN=100;AF=0;nhomalt=0')]
  return source
 def test_AlleleID_and_VariationID_never_conflated(self):
  source=self.source();key=('chr1',11,'A','C')
  self.assertEqual(source.scv_map(key,'111')[1],[('111','9001')])
  self.assertEqual(source.scv_map(key,'9001')[1],[])
  self.assertEqual(source.scv_map(('chr1',11,'A','G'),'111')[1],[])
  source.an.close();source.cv.close()
 def test_original_quality_and_rank_columns_preserved_and_missing_AF_not_zero(self):
  with tempfile.TemporaryDirectory() as directory:
   root=Path(directory);source=self.source();master=root/'master.gz'
   fields=['sample','chrom','pos','ref','alt','GT','FT','DP','AD','GQ','sample_state','site_filter','hgnc_id','clinvar_allele_id','quality_rank','gene']
   rows=[dict(zip(fields,['CONTROL_A','chr1','11','A','C','./.','.','5','3,2','.','unknown_no_call','PASS','1','111','9','CONTROL_GENE'])),dict(zip(fields,['CONTROL_B','chr1','11','A','C','0/0','LowGQ','7','7,0','2','filtered_reference','PASS','1','9001','8','CONTROL_GENE']))]
   with gzip.open(master,'wt',newline='') as stream:
    writer=csv.DictWriter(stream,fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
   output=root/'output';output.mkdir()
   result=exporter.nuclear(master,output,source,{}, {'HGNC:1'},['CONTROL_A','CONTROL_B'],exporter.ABSOLUTE_DEADLINE)
   self.assertEqual(result['selected_rows'],2);self.assertFalse(result['primary_ranking_changed'])
   with gzip.open(output/'cohort/nuclear-sample-rows.tsv.gz','rt') as stream:actual=list(csv.DictReader(stream,delimiter='\t'))
   for expected,value in zip(rows,actual):self.assertEqual({k:value[k] for k in fields},expected)
   self.assertEqual(actual[1]['pop_validated_VariationIDs'],'');self.assertIn('CONFLICT',actual[1]['pop_SCV_map_status'])
   with gzip.open(output/'cohort/nuclear-population-long.tsv.gz','rt') as stream:pop=list(csv.DictReader(stream,delimiter='\t'))
   self.assertEqual(pop[0]['AF'],'0');self.assertEqual(pop[1]['AF'],'');self.assertEqual(pop[1]['match_status'],'MISSING_NOT_AF_ZERO')
   source.an.close();source.cv.close()
 def mt_source(self,root,alt='C,G',raw='1/2:10,3,4:0.1,0.2:17:LowDepth'):
  vcf=root/'final.vcf.gz';sample='CONTROL_A'
  with gzip.open(vcf,'wt') as out:out.write('##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCONTROL_A\nchrM\t5\t.\tA\t'+alt+'\t50\tnumt_boundary\tFLAG\tGT:AD:AF:DP:FT\t'+raw+'\n')
  item=exporter.pin(vcf);receipt=root/'receipt.json';receipt.write_text(json.dumps({'status':'PASS','artifacts':{'mitochondria.filtered.vcf.gz':{'bytes':item['bytes'],'sha256':item['sha256']}}}))
  return dict(sample=sample,vcf=item,receipt=exporter.pin(receipt),qc=[])
 def test_MT_multiallelic_raw_call_filters_not_split_or_changed(self):
  with tempfile.TemporaryDirectory() as directory:
   root=Path(directory);entry=self.mt_source(root);stage=root/'output';stage.mkdir();source=type('MTSource',(),{})()
   result=exporter.mt_exports([entry],stage,source,exporter.ABSOLUTE_DEADLINE)
   with gzip.open(stage/'samples/CONTROL_A/MT-population-reference.tsv.gz','rt') as stream:row=next(csv.DictReader(stream,delimiter='\t'))
   self.assertEqual(row['ALT'],'C,G');self.assertEqual(row['FILTER'],'numt_boundary');self.assertEqual(row['AF'],'0.1,0.2');self.assertEqual(row['AD'],'10,3,4');self.assertEqual(row['GT'],'1/2');self.assertEqual(row['gnomAD_MT_reference_json'],'');self.assertIn('NO_FLATTEN',row['gnomAD_MT_match_status']);self.assertEqual(result[0]['rows'],1)
 def test_MT_sample_heteroplasmy_not_replaced_by_population_AF(self):
  with tempfile.TemporaryDirectory() as directory:
   root=Path(directory);entry=self.mt_source(root,'C','0/1:69,31:0.31:100:PASS');stage=root/'output';stage.mkdir();source=type('MTSource',(),{})();source.mt=sqlite3.connect(':memory:');source.mt.row_factory=sqlite3.Row
   source.mt.executescript('CREATE TABLE gnomad(pos,ref,alt,AF_het);CREATE TABLE helix(pos,ref,alt,AN,AF_het);')
   source.mt.execute('INSERT INTO gnomad VALUES (5,"A","C",0.0001)');source.mt.execute('INSERT INTO helix VALUES (5,"A","C",NULL,0.0002)')
   exporter.mt_exports([entry],stage,source,exporter.ABSOLUTE_DEADLINE)
   with gzip.open(stage/'cohort/MT-population-reference-long.tsv.gz','rt') as stream:row=next(csv.DictReader(stream,delimiter='\t'))
   self.assertEqual(row['AF'],'0.31');self.assertEqual(row['AD'],'69,31');self.assertEqual(json.loads(row['gnomAD_MT_reference_json'])['AF_het'],0.0001);self.assertIsNone(json.loads(row['Helix_reference_json'])['AN']);source.mt.close()
 def test_master_integrity_pin_rejects_tampering(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'synthetic';path.write_bytes(b'original');item=exporter.pin(path);path.write_bytes(b'changed!')
   with self.assertRaises(ValueError):exporter.verify(item)
 def test_controller_prepared_rule_is_not_execution_approval(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'rule.json';rule={'review_status':'PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING'};path.write_text(json.dumps(rule))
   with self.assertRaises(ValueError):controller.validate(rule,path,controller.sha(path))
 def test_controller_original_absolute_budget_not_extendable(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'rule.json';rule=dict(schema_version=1,review_status='ROOT_REVIEWED_EXACT_LOCAL_ONLY_POPULATION_RULE',no_network=True,no_primary_workflow_mutation=True,max_hours=8,interval_seconds=30,max_output_bytes=1024,deadline_utc='2026-11-07T00:00:00+00:00');path.write_text(json.dumps(rule))
   with self.assertRaises(ValueError):controller.validate(rule,path,controller.sha(path),dt.datetime(2026,11,5,20,tzinfo=dt.timezone.utc))
 def test_both_primary_stages_need_PASS_and_identity(self):
  with tempfile.TemporaryDirectory() as directory:
   root=Path(directory)
   for stage in ['06_annotate','07_report']:
    target=root/'stages'/stage;target.mkdir(parents=True);(target/'workflow-receipt.json').write_text(json.dumps({'status':'PASS','identity':'synthetic'}))
   self.assertTrue(controller.readiness(root));(root/'stages/07_report/workflow-receipt.json').write_text(json.dumps({'status':'FAIL'}))
   with self.assertRaises(ValueError):controller.readiness(root)
if __name__=='__main__':unittest.main()
