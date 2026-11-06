import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
def module(name):
 spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/(name+'.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
h=module('apply_spliceai_allele_overlay');c=h.codec
try:
 import pysam
 HAVE_PYSAM=pysam.__version__=='0.24.1'
except ImportError:HAVE_PYSAM=False
def metadata():return json.dumps(dict(version='1.3.1',distance=500,mask=0,annotation_scope=c.SCOPE,assembly='GRCh38',model_sha256=[str(i)*64 for i in range(1,6)],annotation_sha256='a'*64,reference_sha256='b'*64),sort_keys=True)
def prediction():return dict(variant_key='chr11:47332146:T:G',chrom='chr11',pos='47332146',ref='T',alt='G',gene='MYBPC3',status='SCORED_EXACT_GRCH38_ALLELE_GENE',DS_AG='0.02',DS_AL='0.60',DS_DG='0.03',DS_DL='0.01',DP_AG='-164',DP_AL='-22',DP_DG='-50',DP_DL='-301',raw_annotation='G|MYBPC3|0.02|0.60|0.03|0.01|-164|-22|-50|-301',source_json=metadata(),transcript_scope=c.SCOPE)
class Codec(unittest.TestCase):
 def test_scored_gene_positions_and_source_roundtrip(self):
  allele=['chr11',47332146,'T','G'];p=prediction();info=c.overlay_info(allele,[p],metadata());self.assertEqual(([p],metadata()),c.projection_from_info(info,allele));self.assertEqual(p['raw_annotation'],info['SpliceAI'])
  bad=dict(info,SpliceAI='G|MYBPC3|0|0|0|0|0|0|0|0')
  with self.assertRaises(ValueError):c.projection_from_info(bad,allele)
  with self.assertRaises(ValueError):c.projection_from_info(info,['chr11',47332146,'T','A'])
 def test_outside_and_unsupported_are_not_zero(self):
  self.assertEqual(([],''),c.projection_from_info(c.overlay_info(['chr1',1,'A','T'],[],'')))
  p=prediction();p.update(status='UNSUPPORTED_COMPLEX_ALLELE',raw_annotation='')
  for f in c.SCORE_FIELDS:p[f]=''
  info=c.overlay_info(['chr11',47332146,'T','G'],[p],metadata());self.assertEqual('.',info['SpliceAI']);self.assertEqual([p],c.projection_from_info(info)[0])
  with self.assertRaises(ValueError):c.projection_from_info(dict(SpliceAI='0',WGS_SPLICEAI_STATUS=c.OUTSIDE,WGS_SPLICEAI_EVIDENCE='.'))
 def test_sample_logical_packets_do_not_leak_membership(self):
  original=['sample','chrom','pos','ref','alt','gene','GT','FT'];fields=original+c.SPLICE_FIELDS;rows=[]
  for i,s in enumerate(['S1','S2'],1):
   row=dict(sample=s,chrom='chr11',pos='47332146',ref='T',alt='G',gene='MYBPC3',GT='0|1',FT='PASS');row.update(c.enrich_row(row,[prediction()],metadata()));rows.append(dict(i=i,values=[row[f] for f in fields]))
  packet=dict(v=1,allele=['chr11',47332146,'T','G'],rows=rows,membership=dict(sample_states_json='{"S1":"alternate","S2":"reference"}'),unknown_private={'S2':'private'},scope='cohort',sample=None,evidence={},gene_context_keys=[]);selected=c.select_sample(packet,fields,'S1');self.assertNotIn('S2',json.dumps(selected));rebuilt=c.rebuild_readable(selected,fields);self.assertEqual(prediction()['raw_annotation'],rebuilt['SpliceAI']);self.assertEqual([prediction()],c.projection_from_info(rebuilt,packet['allele'])[0])

@unittest.skipUnless(HAVE_PYSAM and os.environ.get('WGS_ALLELE_OVERLAY_PUBLIC_MODEL'),'Explicit actual-public-model integration on pinned pysam0.24.1')
class PublicIntegration(unittest.TestCase):
 def fixture(self,root):
  modelroot=Path(os.environ['WGS_ALLELE_OVERLAY_PUBLIC_MODEL']);model=h.load(modelroot/'manifest.json');base=root/'direct-core-context-v1';base.mkdir();state=root/'state';state.mkdir();samples=['PUBLIC_S1','PUBLIC_S2','PUBLIC_S3','PUBLIC_S4'];original=['sample','chrom','pos','ref','alt','GT','FT','gene','hgnc_id','transcript'];appended=['allref_control_'+str(i) for i in range(33)];fields=original+appended;rows=[];records=[];tokens=[]
  for chrom,pos,ref,alt,gene in [('chr1',100,'A','T',''),('chr11',47332146,'T','G','MYBPC3')]:
   packetrows=[]
   for sample,gt,ft in zip(samples,['0|1','./.','0/0','1/1'],['PASS','.','LowQual','PASS']):
    row=dict(sample=sample,chrom=chrom,pos=str(pos),ref=ref,alt=alt,GT=gt,FT=ft,gene=gene,hgnc_id='HGNC:7551' if gene else '',transcript='ENST00000545968.3' if gene else '',**{f:'' for f in appended});values=[row[f] for f in fields];rows.append(values);packetrows.append(dict(i=len(rows),values=values))
   packet=dict(v=1,allele=[chrom,pos,ref,alt],scope='cohort',sample=None,rows=packetrows,membership=dict(sample_states_json=json.dumps(dict(zip(samples,['alternate','unknown','filtered','alternate'])))),evidence={},gene_context_keys=[]);token=c.encode_packet(packet);tokens.append(token);records.append(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t43\tLowQual\tCSQ={alt}|original||;WGS_EVIDENCE={token}\tGT:FT\t0|1:PASS\t./.:.\t0/0:LowQual\t1/1:PASS\n')
  core_schema=dict(v=1,scope='cohort',sample=None,codec=c.CODEC,original_fields=original,appended_fields=[],full_fields=original,row_count=len(rows),batches=[dict(id=0,first_i=1,last_i=len(rows),row_count=len(rows),record_count=2,canonical_row_sha256='c'*64)],ordered_batch_merkle_root='d'*64,vector_digest_contract='ordered_batch_canonical_TSV_utf8_LF_sha256_v1');schema=dict(core_schema,appended_fields=appended,full_fields=fields)
  header='##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr11,length=135086622>\n##FILTER=<ID=LowQual,Description="Synthetic public control">\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Original annotation">\n##INFO=<ID=WGS_EVIDENCE,Number=1,Type=String,Description="Unchanged core33 vectors">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filter">\n'+c.schema_header(schema).decode()+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n';vcf=base/'cohort.core-context.vcf.gz'
  with pysam.BGZFile(str(vcf),'wb') as f:f.write((header+''.join(records)).encode())
  pysam.tabix_index(str(vcf),preset='vcf');producer=dict(final_driver='e'*64,core_module='f'*64,frozen_package_manifest='1'*64,context_adapter='2'*64,codec='18e7a2102aceb5b962892242763d7d331523d7669c93f16b12f9fe71a1eb238e');h.dump(base/'core-manifest.json',dict(status='DIRECT_PACKET_BASELINE_PASS',schema=core_schema,source_pins=dict(filter_vcf=model['source_vcf']),producer_code_sha256=producer,synthetic_validation=True));manifest=dict(schema_version=1,status='INTERMEDIATE_CORE33_PASS',final_delivery_complete=False,enriched_vcf=h.pin(vcf),index=h.pin(str(vcf)+'.tbi'),core_manifest=h.pin(base/'core-manifest.json'),core_schema=core_schema,schema=schema,samples=samples,original_records_fields_preserved=True,packet_decode_verified=True,direct_core_prefix_verified=True,row_count=len(rows),record_count=2,producer_code_sha256=producer,synthetic_validation=True);h.dump(base/'manifest.json',manifest);artifacts={p.name:{k:v for k,v in h.pin(p).items() if k!='path'} for p in base.iterdir()};h.dump(base/'workflow-receipt.json',dict(status='PASS',stage='direct-core-context-v1',identity='9'*64,artifacts=artifacts))
  names=['apply_spliceai_allele_overlay.py','continue_spliceai_allele_overlay.py','spliceai_allele_info_codec.py','spliceai_packet_context.py','embed_vcf_evidence.py'];rule=dict(base_root=str(base),spliceai_root=str(modelroot),source_vcf=model['source_vcf'],expected_core_driver_sha256=producer['final_driver'],expected_model_provenance=dict(models=model['models'],annotation=model['annotation'],public_control=model['public_control'],code_sha256=model['producer_code_sha256']),output=str(root/'fully-annotated-vcf-v1'),state_directory=str(state),cpu_limit=2,memory_limit_bytes=4*1024**3,max_output_bytes=1024**3,min_free_bytes=0,code_sha256={str(ROOT/'scripts'/n):h.sha(ROOT/'scripts'/n) for n in names});return rule,records,tokens,header
 def test_actual_full_scores_and_core_packet_bytes_and_indexes(self):
  with tempfile.TemporaryDirectory() as t:
   rule,original,tokens,header=self.fixture(Path(t));result=h.overlay(rule,public_control=True);self.assertEqual('PASS_FULLY_ANNOTATED_VCF_DIRECT',result['status']);self.assertTrue(result['core33_packets_byte_identical']);self.assertEqual(1,result['selected_model_alleles']);self.assertEqual(['chr1','chr11'],result['indexed_contigs_verified']);self.assertEqual(h.load(Path(rule['base_root'])/'manifest.json')['schema'],result['schema'])
   for item in result['applied_artifacts']:
    relative=Path(item['relative_path']);self.assertFalse(relative.is_absolute());self.assertNotIn('..',relative.parts);self.assertEqual('applied-evidence',relative.parts[0]);self.assertEqual('applied-evidence/spliceai/predictions.tsv.gz',relative.as_posix());destination=Path(t)/'view-staging'/relative;destination.parent.mkdir(parents=True,exist_ok=True);os.link(h.verify(item),destination);self.assertEqual(item['sha256'],h.sha(destination))
   with gzip.open(Path(rule['output'])/'cohort.fully-annotated.vcf.gz','rb') as f:lines=list(f)
   data=[line for line in lines if not line.startswith(b'#')];self.assertEqual([v.encode() for v in original],[c.strip_overlay(line,only_spliceai=True) for line in data]);infos=[{k.decode():v.decode() if v is not None else None for k,v in h.info(h.data_parts(line)[0]).items()} for line in data];self.assertEqual(tokens,[v['WGS_EVIDENCE'] for v in infos]);self.assertEqual(([],''),c.projection_from_info(infos[0],['chr1',100,'A','T']));predictions,source=c.projection_from_info(infos[1],['chr11',47332146,'T','G']);self.assertEqual('0.60',predictions[0]['DS_AL']);self.assertEqual('-22',predictions[0]['DP_AL']);self.assertEqual(result['projection'],c.decode_projection_schema(next(line for line in lines if line.startswith(b'##WGS_SPLICEAI_SCHEMA='))))
   proof=dict(status='PASS_ACTUAL_PUBLIC_MODEL_ALLELE_OVERLAY_CORE33_BYTE_PRESERVATION',patient_inputs=False,synthetic_genotypes_only=True,actual_model_manifest=h.pin(Path(rule['spliceai_root'])/'manifest.json'),actual_annotation=infos[1]['SpliceAI'],all_original_record_bytes_preserved=True,all_core_packet_tokens_identical=True,physical_schema_and_Merkle_unchanged=True,projection_contract=c.CONTRACT,record_count=2,row_count=8,selected_model_alleles=1,indexed_contigs=result['indexed_contigs_verified'],producer_sha256=h.sha(ROOT/'scripts/apply_spliceai_allele_overlay.py'),codec_sha256=h.sha(ROOT/'scripts/spliceai_allele_info_codec.py'))
   if os.environ.get('WGS_ALLELE_OVERLAY_PUBLIC_RECEIPT'):Path(os.environ['WGS_ALLELE_OVERLAY_PUBLIC_RECEIPT']).write_text(json.dumps(proof,indent=2)+'\n')
 def test_wrong_accepted_input_refused_without_final_output(self):
  with tempfile.TemporaryDirectory() as t:
   rule,_,_,_=self.fixture(Path(t));rule['source_vcf']=dict(rule['source_vcf'],sha256='0'*64)
   with self.assertRaises(ValueError):h.overlay(rule,public_control=True)
   self.assertFalse(Path(rule['output']).exists())
if __name__=='__main__':unittest.main()
