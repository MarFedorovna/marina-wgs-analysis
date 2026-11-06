import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sqlite3
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
def module(name):
 spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/(name+'.py'));value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value
a=module('direct_public_context_adapter');codec=module('spliceai_vcf_codec')
def fixture():
 c=a.Context.__new__(a.Context);c.db=sqlite3.connect(':memory:');c.db.row_factory=sqlite3.Row;c.db.execute('CREATE TABLE reference(kind TEXT,key TEXT,subkey TEXT,payload BLOB)');c.context_sha='a'*64;c.last_key=None;c.last=None;c.versions='{"population":"4.1.1"}';c.tissue=dict(genes={'HGNC:1':dict(gene='GENE',gtex_gene=[],gtex_transcripts={'ENST1.1':[]},pext_transcripts={'ENST1.1':[]},pext_status='CONTEXT_ONLY',uniprot=[])});c.public_genes={'HGNC:1'};c.coverage={'HGNC:1':dict(hgnc_id='HGNC:1',scope='GENE_BIN_ONLY_NOT_LOCUS_DEPTH',exome_bins=1,genome_bins=2)};c.splices={};c.splice_source='{"version":"1.3.1","distance":500,"mask":0}';c.splice_manifest={}
 identity='chr1:10:A:G'
 def put(kind,key,sub,value):c.db.execute('INSERT INTO reference VALUES(?,?,?,?)',(kind,key,sub,a.pack(value)))
 put('population',identity,'exomes',dict(variant_key=identity,chrom='chr1',pos='10',ref='A',alt='G',dataset='exomes',source_version='4.1.1 public HCM39 CDS',match_status='EXACT_ALLELE',AC='1',AN='200',AF='0.005',nhomalt='0',population_FILTER='AS_VQSR',ancestry_quality_INFO='AC_afr=1;AN_afr=10'))
 put('AN','chr1:10','',dict(AN='0',match_status='EXACT_POSITION_EXOME_AN_ONLY',scope='POSITION_CONTEXT_v4.1_EXOMES_ONLY_NOT_GENOMES_NOT_DEPTH'))
 put('map',identity,'77|88',dict(variant_key=identity,AlleleID='77',VariationID='88',scope='EXACT_GRCh38_ALLELES_SOURCE_ID_MAP_NOT_CLINICAL_CLASSIFICATION'))
 put('SCV','88','SCV1.2',dict(VariationID='88',SCV_accession_version='SCV1.2',source_classification='Uncertain significance',raw_source_json='{"full_assertion":"retained"}'))
 put('citation','88','PubMed|1',dict(VariationID='88',citation_source='PubMed',citation_id='1',primary_article_read=False))
 prediction=dict(variant_key=identity,chrom='chr1',pos='10',ref='A',alt='G',gene='GENE',status='SCORED_EXACT_GRCH38_ALLELE_GENE',DS_AG='0.00',DS_AL='0.60',DS_DG='0.03',DS_DL='0.01',DP_AG='-164',DP_AL='-22',DP_DG='-50',DP_DL='-301',raw_annotation='G|GENE|0.00|0.60|0.03|0.01|-164|-22|-50|-301',source_json=c.splice_source,transcript_scope=a.splice.SCOPE);c.splices[identity]=[prediction]
 row=dict(sample='S1',chrom='chr1',pos='10',ref='A',alt='G',GT='0|1',FT='PASS',hgnc_id='HGNC:1',gene='GENE',transcript='ENST1.1',clinvar_allele_id='77',candidate_eligible='0')
 return c,row
class Unit(unittest.TestCase):
 def test_direct_context_exact_alleles_scvs_and_original_fields(self):
  c,row=fixture();before=copy.deepcopy(row);extra=c.augment(row,1);self.assertEqual(before,row);self.assertEqual(46,len(extra));self.assertEqual('0',extra['allref_exome_AN_v41']);self.assertIn('0.005',extra['allref_gnomad_exomes_summary_json']);self.assertEqual('1',extra['allref_SCV_reference_count']);self.assertIn('NOT_LOCUS_DEPTH',extra['allref_gene_coverage_summary_json']);self.assertEqual('0.00',extra['spliceai_DS_AG']);self.assertEqual('-22',extra['spliceai_DP_AL']);self.assertIn('GNOMAD_SOURCE_FILTER_REVIEW',extra['allref_review_flags']);self.assertEqual('0',row['candidate_eligible']);self.assertIn('full_assertion',c.allele_evidence('chr1',10,'A','G')['SCV'][0]['raw_source_json']);c.close()
 def test_alleleid_is_not_variationid_or_false_zero(self):
  c,row=fixture();row['clinvar_allele_id']='88';extra=c.augment(row,1);self.assertIn('CONFLICT',extra['allref_SCV_map_status']);self.assertEqual('',extra['allref_SCV_reference_count']);row.update(chrom='chr22',pos='123',gene='',hgnc_id='',ref='T',alt='A');extra=c.augment(row,2);self.assertEqual('',extra['allref_exome_AN_v41']);self.assertEqual('',extra['spliceai_DS_AG']);self.assertIn('NOT_ZERO',extra['spliceai_status']);c.close()
 def test_different_gene_preserves_model_scope_not_transcript_score(self):
  c,row=fixture();row.update(gene='OTHER',hgnc_id='HGNC:999');extra=c.augment(row,1);self.assertEqual('',extra['spliceai_DS_AL']);self.assertIn('NO_TRANSCRIPT_ASSIGNMENT',extra['spliceai_status']);self.assertIn('0.60',extra['spliceai_prediction_json']);c.close()
 def test_sample_projection_and_original_info_preservation(self):
  c,row=fixture();fields=list(row)+a.APPENDED+a.SPLICE_FIELDS;rows=[]
  for n,sample in enumerate(['S1','S2'],1):
   r=dict(row,sample=sample);values=dict(r,**c.augment(r,n));rows.append(dict(i=n,values=[values[field] for field in fields]))
  packet=dict(v=1,allele=['chr1',10,'A','G'],rows=rows,evidence=c.allele_evidence('chr1',10,'A','G'),gene_context_keys=['HGNC:1'],scope='cohort',sample=None);selected=codec.select_sample(packet,fields,'S1');tags=codec.rebuild_readable(selected,fields);self.assertNotIn('S2',str(tags));self.assertEqual('G|GENE|0.00|0.60|0.03|0.01|-164|-22|-50|-301',tags['SpliceAI']);self.assertEqual('0.005',tags['WGS_GNOMAD_EXOMES_AF']);self.assertEqual('.',tags['WGS_GNOMAD_GENOMES_AF']);self.assertEqual(selected,codec.decode_packet(codec.encode_packet(selected)))
  original=b'chr1\t10\t.\tA\tG\t.\tPASS\tCSQ=G|TX.1\tGT\t0|1\n';parts=original.rstrip(b'\n').split(b'\t');parts[7]+=b';'+';'.join(k+'='+v for k,v in tags.items()).encode();self.assertEqual(original,codec.strip_overlay(b'\t'.join(parts)+b'\n'));c.close()

@unittest.skipUnless(os.environ.get('WGS_DIRECT_PUBLIC_INDEX_CONTROL'),'Explicit existing public index control')
class ActualPublic(unittest.TestCase):
 @unittest.skipUnless(os.environ.get('WGS_DIRECT_SPLICEAI_MODEL_CONTROL'),'Explicit official-model public adapter control only')
 def test_actual_model_completion_projects_and_encodes_full_scores(self):
  import pysam
  source_path=Path('/srv/marina/private/spliceai-reference-context-candidate-v1/scripts/run_spliceai_reference_context.py');self.assertEqual('c75b3aa5c6ced8854e30bc227a188ba1d82ba9d2b8e00e4307a5f52be72dd6e0',a.base.sha(source_path));spec=importlib.util.spec_from_file_location('frozen_public_model_scorer',source_path);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h);root=Path('/srv/marina');package=root/'resources/spliceai-research-runtime-v1/venv/lib/python3.10/site-packages/spliceai';reference=root/'references/broad-hg38-v0-altaware-g1575676516681666-prebuilt/Homo_sapiens_assembly38.fasta'
  with tempfile.TemporaryDirectory() as t:
   folder=Path(t);state=folder/'state';state.mkdir();source=folder/'public.vcf.gz'
   if os.environ.get('WGS_DIRECT_SPLICEAI_PUBLIC_ARTIFACTS'):
    retained=Path(os.environ['WGS_DIRECT_SPLICEAI_PUBLIC_ARTIFACTS']);source=retained.with_name(retained.name+'.input.vcf.gz')
    if source.exists():raise ValueError('Fresh retained public input required')
   with pysam.BGZFile(str(source),'wb') as f:f.write(b'##fileformat=VCFv4.2\n##contig=<ID=chr11,length=135086622>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr11\t47332146\t.\tT\tG\t.\tPASS\t.\n')
   pysam.tabix_index(str(source),preset='vcf');rule=dict(source_pins=[],output=os.environ.get('WGS_DIRECT_SPLICEAI_PUBLIC_ARTIFACTS',str(folder/'output')),state_directory=str(state),source_vcf=h.pin(source),primary_receipt={},panel=h.pin(root/'resources/public-hcm-cds-v2/panel-manifest.json'),annotation=h.pin(package/'annotations/grch38.txt'),reference=dict(path=str(reference),bytes=reference.stat().st_size,mtime_ns=reference.stat().st_mtime_ns,sha256='93157a161863464c9435062fd67c173fdaf99cb8b32f1455018361387ffa5564'),public_control=h.pin(root/'private/fully-annotated-vcf-candidate-v1/spliceai-public-cpu-control.json'),models=[h.pin(p) for p in sorted((package/'models').glob('*.h5'))],distance=500,mask=0,cpu_limit=2,min_free_bytes=0,max_sqlite_bytes=1024**3,max_output_bytes=1024**3,expected_variants=1,code_sha256={},runtime={});model=h.run(rule,h.DEADLINE)
   context=a.Context(Path(os.environ['WGS_DIRECT_PUBLIC_INDEX_CONTROL']),Path(rule['output'])/'manifest.json',verify_index=False,expected_source_vcf=rule['source_vcf']);row=dict(sample='PUBLIC_CONTROL',chrom='chr11',pos='47332146',ref='T',alt='G',gene='MYBPC3',hgnc_id='HGNC:7551',transcript='ENST00000545968.3',clinvar_allele_id='',GT='0|1',FT='PASS');extra=context.augment(row,1);self.assertEqual('0.60',extra['spliceai_DS_AL']);self.assertEqual('-22',extra['spliceai_DP_AL']);self.assertIn('NOT_ENST49_MAPPING',extra['spliceai_transcript_scope']);fields=list(row)+a.APPENDED+a.SPLICE_FIELDS;values=dict(row,**extra);packet=dict(v=1,allele=['chr11',47332146,'T','G'],rows=[dict(i=1,values=[values[f] for f in fields])],evidence=context.allele_evidence('chr11',47332146,'T','G'),gene_context_keys=['HGNC:7551'],scope='cohort',sample=None);tags=codec.rebuild_readable(packet,fields);self.assertEqual('G|MYBPC3|0.02|0.60|0.03|0.01|-164|-22|-50|-301',tags['SpliceAI']);self.assertEqual(packet,codec.decode_packet(codec.encode_packet(packet)));context.close()
   proof=dict(status='PASS_ACTUAL_OFFICIAL_MODEL_DIRECT_PROJECTION_AND_VCF_CODEC',patient_inputs=False,source_scorer_sha256=a.base.sha(source_path),adapter_sha256=a.base.sha(ROOT/'scripts/direct_public_context_adapter.py'),codec_sha256=a.base.sha(ROOT/'scripts/spliceai_vcf_codec.py'),actual_prediction=tags['SpliceAI'],transcript_scope=extra['spliceai_transcript_scope'],full_packet_roundtrip=True,exact_source_lineage_checked=True)
   if os.environ.get('WGS_DIRECT_SPLICEAI_CONTROL_RECEIPT'):Path(os.environ['WGS_DIRECT_SPLICEAI_CONTROL_RECEIPT']).write_text(json.dumps(proof,indent=2)+'\n')
 def test_actual_public_index_context_and_header_hashes(self):
  path=Path(os.environ['WGS_DIRECT_PUBLIC_INDEX_CONTROL']);c=a.Context(path);context=c.tissue['genes']['HGNC:7551'];transcript=next(iter(context['gtex_transcripts']), 'ENST00000545968.3');row=dict(sample='SYNTHETIC_PUBLIC_CONTROL',chrom='chr11',pos='47332146',ref='T',alt='G',gene='MYBPC3',hgnc_id='HGNC:7551',transcript=transcript,clinvar_allele_id='');extra=c.augment(row,1);self.assertIn('APPLIED_DIRECT_EXACT_HGNC',extra['allref_tissue_status'])
  # A rare public model control need not occur in population reference data.
  identity=c.db.execute("SELECT p.key FROM reference p JOIN reference m ON m.kind='map' AND m.key=p.key WHERE p.kind='population' LIMIT 1").fetchone()[0];chrom,pos,ref,alt=identity.split(':');row.update(chrom=chrom,pos=pos,ref=ref,alt=alt,gene='',hgnc_id='',transcript='');extra=c.augment(row,2);evidence=c.allele_evidence(chrom,pos,ref,alt);self.assertTrue(evidence['population']);self.assertEqual(identity,extra['allref_variant_key']);self.assertTrue(evidence['allele_map']);self.assertTrue(evidence['SCV']);bundles=c.header_bundles()
  for name,bundle in bundles.items():self.assertEqual(bundle,codec.decode_packet(codec.encode_packet(bundle)))
  self.assertEqual(c.manifest['public_header_sha256'],hashlib.sha256(codec.encode_packet(bundles['public_reference']).encode()).hexdigest());self.assertEqual(c.manifest['gene_header_sha256'],hashlib.sha256(codec.encode_packet(bundles['gene_context']).encode()).hexdigest());self.assertIn('ancestry_quality_INFO',evidence['population'][0]);proof=dict(status='PASS_ACTUAL_PUBLIC_DIRECT_CONTEXT_AND_FULL_HEADER_ROUNDTRIP',patient_inputs=False,index_manifest=a.base.pin(path),index=c.manifest['index'],known_public_allele=identity,population_rows=len(evidence['population']),validated_source_maps=len(evidence['allele_map']),full_SCV_records=len(evidence['SCV']),citation_records=len(evidence['citations']),public_header_sha256=c.manifest['public_header_sha256'],gene_header_sha256=c.manifest['gene_header_sha256'],producer_sha256=a.base.sha(ROOT/'scripts/direct_public_context_adapter.py'));c.close()
  if os.environ.get('WGS_DIRECT_PUBLIC_CONTROL_RECEIPT'):Path(os.environ['WGS_DIRECT_PUBLIC_CONTROL_RECEIPT']).write_text(json.dumps(proof,indent=2)+'\n')
if __name__=='__main__':unittest.main()
