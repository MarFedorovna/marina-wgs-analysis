import csv
import datetime as dt
import gzip
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
import tempfile
import unittest
import zlib
import base64

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('embed',ROOT/'scripts/embed_vcf_evidence.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
spec=importlib.util.spec_from_file_location('controller',ROOT/'scripts/continue_fully_annotated_vcf.py');c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
try:
 import pysam
 HAVE_PYSAM=pysam.__version__=='0.24.1'
except ImportError:HAVE_PYSAM=False

class Codec(unittest.TestCase):
 def test_lossless_unicode_empty_json_and_missing(self):
  p=dict(v=1,allele=['chr1',10,'A','G'],rows=[dict(i=2,values=['S1','0|1','','.',json.dumps({'x':'кириллица\t\n"'})])],evidence={})
  self.assertEqual(p,m.decode_packet(m.encode_packet(p)))
  self.assertEqual(m.decode_schema(m.encode_schema(dict(v=1,full_fields=['sample','GT']))),dict(v=1,full_fields=['sample','GT']))
 def test_sample_projection_readables_keep_exact_transcripts(self):
  fields=['sample','transcript','FATHMM','GT'];p=dict(v=1,rows=[dict(i=1,values=['S1','TX.1','deleterious','0|1']),dict(i=2,values=['S2','TX.2','tolerated','./.'])],evidence={'SCV':[{'public':'only'}]},gene_context_keys=[])
  sub=m.select_sample(p,fields,'S1');self.assertEqual([1],[r['i'] for r in sub['rows']]);self.assertNotIn('S2',str(m.rebuild_readable(sub,fields)));self.assertIn('TX.1',m.rebuild_readable(sub,fields)['WGS_TX_SCORES']);self.assertNotIn('TX.2',m.rebuild_readable(sub,fields)['WGS_TX_SCORES'])
 def test_reject_bad_version_and_expansion(self):
  for token in ['!!!',m.encode_packet(dict(v=2))]:
   with self.assertRaises(ValueError):m.decode_packet(token)
  old=m.MAX_PACKET_BYTES
  try:
   m.MAX_PACKET_BYTES=32;token=base64.urlsafe_b64encode(zlib.compress(b'x'*1000)).decode().rstrip('=')
   with self.assertRaises(ValueError):m.decode_packet(token)
  finally:m.MAX_PACKET_BYTES=old
 def test_strip_preserves_every_original_byte(self):
  for info in [b'.',b'X=1;CSQ=A|B||;FLAG']:
   original=b'chr1\t1\t.\tA\tG\t.\tLowQual\t'+info+b'\tGT:FT\t0|1:PASS\t./.:.\t0/0:PASS\t1/1:PASS\r\n';parts=original.rstrip(b'\r\n').split(b'\t');parts[7]=(info+b';' if info!=b'.' else b'')+b'WGS_EVIDENCE=abc;WGS_GENE=G';self.assertEqual(original,m.strip_overlay(b'\t'.join(parts)+b'\r\n'))
 def test_unreviewed_controller_never_executes(self):
  with tempfile.TemporaryDirectory() as t:
   path=Path(t)/'rule.json';path.write_text('{}')
   with self.assertRaises(ValueError):c.validate({},path,c.sha(path),dt.datetime(2026,11,5,tzinfo=dt.timezone.utc))
 def test_scalar_count_fields_exact_only_missing_never_zero(self):
  p=dict(rows=[],evidence={'population':[dict(dataset='exomes',source_version='4.1.1',match_status='EXACT_ALLELE',AC='1',AN='200',AF='0.005',nhomalt='0')]});tags=m.rebuild_readable(p,[]);self.assertEqual('0.005',tags['WGS_GNOMAD_EXOMES_AF']);self.assertEqual('0',tags['WGS_GNOMAD_EXOMES_nhomalt']);self.assertEqual('.',tags['WGS_GNOMAD_GENOMES_AN']);p['evidence']['population'][0]['match_status']='MISSING_NOT_AF_ZERO';self.assertEqual('.',m.rebuild_readable(p,[])['WGS_GNOMAD_EXOMES_AF'])

@unittest.skipUnless(HAVE_PYSAM,'Exact pysam0.24.1 public control integration runs on server')
class Integration(unittest.TestCase):
 def fixture(self,root,include_empty_annotation_rows=True):
  workflow=root/'workflow-v1';stage6=workflow/'stages/06_annotate';stage6.mkdir(parents=True);combined=root/'combined-annotation-master-v1';combined.mkdir();state=root/'state';state.mkdir()
  header='##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Original">\n##FILTER=<ID=LowQual,Description="Synthetic filtered control">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filter">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n'
  records=['chr1\t10\t.\tA\tG\t43\tLowQual\tCSQ=G|missense|TX.1,G|synonymous|TX.2\tGT:FT\t0|1:PASS\t./.:.\t0/0:PASS\t1/1:PASS\n','chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT:FT\t0/0:PASS\t./.:.\t0/0:PASS\t0/0:PASS\n']
  with pysam.BGZFile(str(stage6/'annotated.vcf.gz'),'wb') as f:f.write((header+''.join(records)).encode())
  p6=dict(status='PASS',identity='a'*64,artifacts={'annotated.vcf.gz':{k:v for k,v in m.pin(stage6/'annotated.vcf.gz').items() if k!='path'}});(stage6/'workflow-receipt.json').write_text(json.dumps(p6))
  fields=['sample','chrom','pos','ref','alt','GT','FT','hgnc_id','gene','transcript','FATHMM','local_evidence_json','consequence'];extra=['allref_tissue_status','allref_population_status'];full=fields+extra;rows=[]
  for sample,gt,ft in [('S1','0|1','PASS'),('S2','./.','.'),('S3','0/0','PASS'),('S4','1/1','PASS')]:
   for tx in ['TX.1','TX.2']:rows.append([sample,'chr1','10','A','G',gt,ft,'HGNC:1','GENE',tx,'deleterious' if tx=='TX.1' else '','{"source":"fixed-public-control"}','missense_variant','VERIFIED','VERIFIED'])
  if include_empty_annotation_rows:
   for sample,gt,ft in [('S1','0/0','PASS'),('S2','./.','.'),('S3','0/0','PASS'),('S4','0/0','PASS')]:rows.append([sample,'chr1','20','C','T',gt,ft,'','','','','{}','','OUTSIDE','OUTSIDE'])
  path=combined/'master_evidence.enriched.tsv.gz'
  with gzip.open(path,'wt',newline='') as f:
   writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(full);writer.writerows(rows)
  aux=[]
  def table(name,header,values):
   p=combined/name
   with gzip.open(p,'wt',newline='') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(values)
   aux.append(dict(**m.pin(p),role='APPLIED_NUCLEAR_REFERENCE_LONG_TABLE'))
  allele='chr1:10:A:G'
  table('nuclear-population-long.tsv.gz',['variant_key','chrom','pos','ref','alt','dataset','source_version','match_status','AC','AN','AF','nhomalt','population_FILTER','ancestry_quality_INFO'],[[allele,'chr1',10,'A','G','exomes','4.1.1 public HCM39 CDS','EXACT_ALLELE',1,200,0.005,0,'PASS','AC_afr=1;AN_afr=10']])
  table('nuclear-exome-AN-v41.tsv.gz',['variant_key','chrom','pos','AN','match_status','scope'],[[allele,'chr1',10,0,'EXACT','POSITION_CONTEXT_v4.1_EXOMES_ONLY_NOT_GENOMES_NOT_DEPTH']])
  table('ClinVar-validated-allele-map.tsv.gz',['variant_key','AlleleID','VariationID','scope'],[[allele,'77','88','EXACT_GRCh38_ALLELES_SOURCE_ID_MAP_NOT_CLINICAL_CLASSIFICATION']])
  table('ClinVar-SCV-long.tsv.gz',['variant_key','VariationID','SCV_accession_version','raw_source_json'],[[allele,'88','SCV00001.2','{"public_citation":"PMID1"}']])
  table('ClinVar-citation-links.tsv.gz',['variant_key','VariationID','citation_source','citation_id'],[[allele,'88','PubMed','1']])
  table('gene-coverage-context-long.tsv.gz',['hgnc_id','gene','scope','mean'],[['HGNC:1','GENE','GENE_BIN_CONTEXT_ONLY_NO_EXACT_LOCUS_DEPTH_CLAIM','23']])
  digest=hashlib.sha256();digest.update(m.canonical_tsv(fields))
  for row in rows:digest.update(m.canonical_tsv(row[:len(fields)]))
  manifest=dict(status='PASS_COMBINED_ANNOTATION_MASTER',ranking_eligibility_unchanged=True,original_fields=fields,appended_fields=extra,enriched_master=m.pin(path),row_count=len(rows),invariant_digest=digest.hexdigest(),auxiliary_artifacts=aux,stage_receipts=[m.pin(stage6/'workflow-receipt.json')]);(combined/'manifest.json').write_text(json.dumps(manifest));artifacts={name:{k:v for k,v in m.pin(combined/name).items() if k!='path'} for name in ['manifest.json','master_evidence.enriched.tsv.gz']};(combined/'workflow-receipt.json').write_text(json.dumps(dict(status='PASS',stage='combined-annotation-master-v1',identity='b'*64,artifacts=artifacts)))
  source=root/'public.json';source.write_text('{"features":[{"begin":1,"end":4}],"transcripts":["TX.1"]}')
  context=root/'context-index.json';context.write_text(json.dumps(dict(verified_raw_sources={str(source):{k:v for k,v in m.pin(source).items() if k!='path'}},source_manifests={})))
  rule=dict(output=str(root/'fully-annotated-vcf-v1'),state_directory=str(state),workflow_root=str(workflow),combined_root=str(combined),tissue_context_index=m.pin(context),max_header_inflated_bytes=64*1024**2,memory_limit_bytes=4*1024**3,cpu_limit=2,min_free_bytes=0,max_sqlite_bytes=8*1024**3,max_output_bytes=1024**3,max_inflation_ratio=24,code_sha256={})
  return rule,records,rows
 def test_complete_portable_roundtrip_filter_phase_missing_multi_tx_empty_annotation(self):
  with tempfile.TemporaryDirectory() as t:
   rule,records,rows=self.fixture(Path(t));result=m.embed(rule);self.assertTrue(result['packet_decode_verified']);self.assertEqual(12,result['row_count']);self.assertEqual(2,result['record_count']);self.assertEqual(1,result['indexed_seek_checks']);self.assertTrue(result['original_headers_preserved']);out=Path(rule['output'])/'cohort.fully-annotated.vcf.gz'
   with gzip.open(out,'rb') as f:lines=list(f)
   with pysam.VariantFile(str(out)) as f:
    first=next(f);self.assertEqual(200,first.info['WGS_GNOMAD_EXOMES_AN']);self.assertAlmostEqual(.005,first.info['WGS_GNOMAD_EXOMES_AF']);self.assertIsNone(first.info['WGS_GNOMAD_GENOMES_AF']);self.assertEqual('Integer',f.header.info['WGS_GNOMAD_EXOMES_AN'].type)
   reconstructed=Path(t)/'reconstructed.tsv';subprocess.run([sys.executable,str(ROOT/'scripts/embed_vcf_evidence.py'),'decode','--vcf',str(out),'--master',str(reconstructed)],check=True,capture_output=True)
   with reconstructed.open(newline='') as f:decoded=list(csv.reader(f,delimiter='\t'))
   self.assertEqual(rows,decoded[1:])
   self.assertEqual([r.encode() for r in records],[m.strip_overlay(line) for line in lines if not line.startswith(b'#')]);self.assertTrue(any(b'WGS_GENE_CONTEXT' in line for line in lines));public=m.decode_header_bundle(lines,'WGS_PUBLIC_REFERENCE');self.assertIn('features',public['sources'][0]['content_utf8']);packets=[m.decode_packet(dict(x.split(b'=',1) for x in l.split(b'\t')[7].split(b';'))[b'WGS_EVIDENCE'].decode()) for l in lines if not l.startswith(b'#')];self.assertEqual(4,len(packets[1]['rows']));self.assertEqual('AC_afr=1;AN_afr=10',packets[0]['evidence']['population'][0]['ancestry_quality_INFO']);self.assertEqual('0',packets[0]['evidence']['exome_AN'][0]['AN']);self.assertEqual(rows,[r['values'] for packet in packets for r in packet['rows']]);self.assertEqual('SCV00001.2',packets[0]['evidence']['SCV'][0]['SCV_accession_version'])
 def test_zero_master_row_allele_refused(self):
  with tempfile.TemporaryDirectory() as t:
   rule,_,_=self.fixture(Path(t),False)
   with self.assertRaises(ValueError):m.embed(rule)
   self.assertFalse(Path(rule['output']).exists())
 @unittest.skipUnless(os.environ.get('WGS_PUBLIC_CONTEXT_INDEX_CONTROL'),'Actual approved public index supplied only for explicit server control')
 def test_actual_public_reference_bundle_roundtrip(self):
  with tempfile.TemporaryDirectory() as t:
   rule,_,_=self.fixture(Path(t));rule['tissue_context_index']=m.pin(Path(os.environ['WGS_PUBLIC_CONTEXT_INDEX_CONTROL']));result=m.embed(rule)
   with gzip.open(Path(rule['output'])/'cohort.fully-annotated.vcf.gz','rb') as f:lines=[line for line in f if line.startswith(b'#')]
   bundle=m.decode_header_bundle(lines,'WGS_PUBLIC_REFERENCE');self.assertGreater(len(bundle['sources']),30);self.assertTrue(result['public_header_bundle_verified'])
   evidence=dict(status='PASS_ACTUAL_PUBLIC_REFERENCE_HEADER_BGZF_TABIX_ROUNDTRIP',synthetic_genotypes_only=True,patient_inputs_read=False,source_context_index=rule['tissue_context_index'],source_objects=len(bundle['sources']),source_raw_bytes=sum(s['bytes'] for s in bundle['sources']),public_header_chunks=sum(line.startswith(b'##WGS_PUBLIC_REFERENCE=') for line in lines),public_reference_sha256=result['public_reference_sha256'],record_count=result['record_count'],row_count=result['row_count'],indexed_seek_checks=result['indexed_seek_checks'],original_headers_preserved=result['original_headers_preserved'],public_header_bundle_verified=result['public_header_bundle_verified'],pysam_version=pysam.__version__,producer_sha256=m.sha(ROOT/'scripts/embed_vcf_evidence.py'))
   if os.environ.get('WGS_PUBLIC_CONTEXT_CONTROL_RECEIPT'):Path(os.environ['WGS_PUBLIC_CONTEXT_CONTROL_RECEIPT']).write_text(json.dumps(evidence,indent=2)+'\n')
 def test_orphan_master_refused(self):
  with tempfile.TemporaryDirectory() as t:
   rule,_,_=self.fixture(Path(t));manifest_path=Path(rule['combined_root'])/'manifest.json';manifest=json.loads(manifest_path.read_text());manifest['row_count']+=1;manifest_path.write_text(json.dumps(manifest));receipt_path=Path(rule['combined_root'])/'workflow-receipt.json';receipt=json.loads(receipt_path.read_text());receipt['artifacts']['manifest.json']={k:v for k,v in m.pin(manifest_path).items() if k!='path'};receipt_path.write_text(json.dumps(receipt))
   with self.assertRaises(ValueError):m.embed(rule)
   self.assertFalse(Path(rule['output']).exists())
if __name__=='__main__':unittest.main()
