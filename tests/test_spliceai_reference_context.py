import datetime as dt
import gzip
import importlib.util
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('scorer',ROOT/'scripts/run_spliceai_reference_context.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
class Unit(unittest.TestCase):
 def test_true_zero_is_available_not_missing(self):
  record=SimpleNamespace(alts=('G',));gene,values,available=h.parse_result('G|MYBPC3|0.00|0.60|0.03|0.01|-164|-22|-50|-301',record,500);self.assertTrue(available);self.assertEqual('0.00',values['DS_AG']);self.assertEqual('-22',values['DP_AL'])
 def test_wrong_allele_bad_positions_and_partial_missing_refused(self):
  record=SimpleNamespace(alts=('G',))
  for raw in ['T|MYBPC3|0|0|0|0|0|0|0|0','G|MYBPC3|0|0|0|0|501|0|0|0','G|MYBPC3|.|0|0|0|0|0|0|0']:
   with self.assertRaises(ValueError):h.parse_result(raw,record,500)
 def test_unsupported_complex_no_fake_zero(self):
  self.assertEqual('UNSUPPORTED_COMPLEX_ALLELE',h.classify(SimpleNamespace(ref='TA',alts=('GC',)),500));self.assertIsNone(h.classify(SimpleNamespace(ref='T',alts=('TG',)),500));self.assertEqual('UNSUPPORTED_NON_ACGT_OR_SYMBOLIC_ALLELE',h.classify(SimpleNamespace(ref='A',alts=('<DEL>',)),500))
 def test_schema_does_not_claim_transcript_mapping(self):
  self.assertEqual(13,len(h.SPLICE_FIELDS));self.assertIn('NOT_ENST49_MAPPING',h.SCOPE)

@unittest.skipUnless(os.environ.get('WGS_SPLICEAI_PUBLIC_STAGE_CONTROL'),'Explicit existing-runtime public control only')
class PublicIntegration(unittest.TestCase):
 def test_actual_official_model_stage_and_resume(self):
  import pysam
  root=Path('/srv/marina');package=root/'resources/spliceai-research-runtime-v1/venv/lib/python3.10/site-packages/spliceai';control=root/'private/fully-annotated-vcf-candidate-v1/spliceai-public-cpu-control.json';reference=root/'references/broad-hg38-v0-altaware-g1575676516681666-prebuilt/Homo_sapiens_assembly38.fasta'
  with tempfile.TemporaryDirectory() as t:
   folder=Path(t);state=folder/'state';state.mkdir();source=folder/'public.vcf.gz'
   header='##fileformat=VCFv4.2\n##contig=<ID=chr11,length=135086622>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n';records='chr11\t47332146\t.\tT\tG\t.\tPASS\t.\nchr11\t47332147\t.\tAA\tGC\t.\tPASS\t.\n'
   with pysam.BGZFile(str(source),'wb') as f:f.write((header+records).encode())
   pysam.tabix_index(str(source),preset='vcf');models=[h.pin(p) for p in sorted((package/'models').glob('*.h5'))];rule=dict(source_pins=[],output=str(folder/'output'),state_directory=str(state),source_vcf=h.pin(source),primary_receipt={},panel=h.pin(root/'resources/public-hcm-cds-v2/panel-manifest.json'),annotation=h.pin(package/'annotations/grch38.txt'),reference=dict(path=str(reference),bytes=reference.stat().st_size,mtime_ns=reference.stat().st_mtime_ns,sha256='93157a161863464c9435062fd67c173fdaf99cb8b32f1455018361387ffa5564'),public_control=h.pin(control),models=models,distance=500,mask=0,cpu_limit=2,min_free_bytes=0,max_sqlite_bytes=1024**3,max_output_bytes=1024**3,expected_variants=2,code_sha256={},runtime={})
   result=h.run(rule,h.DEADLINE);self.assertEqual(1,result['status_counts']['SCORED_EXACT_GRCH38_ALLELE_GENE']);self.assertEqual(1,result['status_counts']['UNSUPPORTED_COMPLEX_ALLELE'])
   import csv
   with gzip.open(Path(rule['output'])/'predictions.tsv.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
   scored=next(r for r in rows if r['status']=='SCORED_EXACT_GRCH38_ALLELE_GENE');self.assertEqual('0.60',scored['DS_AL']);self.assertEqual('-22',scored['DP_AL']);self.assertTrue(all(next(r for r in rows if r['status']=='UNSUPPORTED_COMPLEX_ALLELE')[field]=='' for field in h.SCORE_FIELDS))
   rule['output']=str(folder/'resumed');resumed=h.run(rule,h.DEADLINE);self.assertEqual(result['status_counts'],resumed['status_counts'])
   proof=dict(status='PASS_PUBLIC_OFFICIAL_SPLICEAI_STAGE_AND_CHECKPOINT_RESUME',patient_inputs=False,public_control=str(control),variants=2,scored=1,unsupported_explicit_missing=1,DS_AL=scored['DS_AL'],DP_AL=scored['DP_AL'],checkpoint_resume=True,producer_sha256=h.sha(ROOT/'scripts/run_spliceai_reference_context.py'))
   if os.environ.get('WGS_SPLICEAI_STAGE_CONTROL_RECEIPT'):Path(os.environ['WGS_SPLICEAI_STAGE_CONTROL_RECEIPT']).write_text(json.dumps(proof,indent=2)+'\n')
if __name__=='__main__':unittest.main()
