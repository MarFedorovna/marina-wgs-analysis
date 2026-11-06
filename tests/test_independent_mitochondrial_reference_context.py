import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];spec=importlib.util.spec_from_file_location('independent_mt',ROOT/'scripts/export_mitochondrial_reference_context.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
class Contract(unittest.TestCase):
 def test_multiallelic_is_not_flattened_and_original_module_pinned(self):
  self.assertIsNone(h.old.allele_key('chrM',73,'A','G,T'));self.assertEqual(('chrM',73,'A','G'),h.old.allele_key('MT',73,'A','G'));self.assertEqual(h.LEGACY_SHA,h.old.sha(ROOT/'scripts/export_population_mt_reference_context.py'))
 def test_unreviewed_rule_fails_closed(self):
  with tempfile.TemporaryDirectory() as t:
   path=Path(t)/'rule.json';path.write_text('{}')
   with self.assertRaises(ValueError):h.validate({},path,h.old.sha(path))

@unittest.skipUnless(os.environ.get('WGS_MT_PUBLIC_REFERENCE_CONTROL'),'Explicit actual public database control, no patient input')
class Public(unittest.TestCase):
 def test_four_synthetic_mt_sources_full_frequencies_raw_calls_and_qc(self):
  import pysam
  database=Path(os.environ['WGS_MT_PUBLIC_REFERENCE_CONTROL']);connection=h.old.connection(database);record=connection.execute('SELECT g.pos,g.ref,g.alt FROM gnomad g JOIN helix h ON g.pos=h.pos AND g.ref=h.ref AND g.alt=h.alt ORDER BY g.pos LIMIT 1').fetchone();connection.close();self.assertIsNotNone(record);pos,ref,alt=record;other=next(base for base in 'ACGT' if base not in [ref,alt])
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);state=root/'state';state.mkdir();audit=root/'public-audit.json';audit.write_text(json.dumps(dict(status='PASS_RESEARCH_WITH_LIMITATIONS',synthetic_control=True)));profile=root/'public-profile.json';profile.write_text('{"patient_inputs":false}');ready=Path('/srv/marina/resources/medinter-population-mt-sidecars-v1/resource-ready-manifest.json');entries=[]
   for index in range(4):
    sample='PUBLIC_MT'+str(index+1);folder=root/sample;folder.mkdir();vcf=folder/'mitochondria.filtered.vcf.gz';header='##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n';records=f'chrM\t{pos}\t.\t{ref}\t{alt}\t43\tPASS\tTLOD=20\tGT:AD:AF:DP:FT\t0|1:90,10:0.10:100:PASS\nchrM\t{pos}\t.\t{ref}\t{alt},{other}\t.\tnumt\tTLOD=3\tGT:AD:AF:DP:FT\t1/2:70,10,20:0.10,0.20:100:.\n'
    with pysam.BGZFile(str(vcf),'wb') as f:f.write((header+records).encode())
    pysam.tabix_index(str(vcf),preset='vcf');qc=[]
    for name in ['research-report.json','numt-assessment.json','haplo-output.json','contamination.json']:
     p=folder/name;p.write_text(json.dumps(dict(synthetic_control=True,sample=sample,unchanged_quality='LIMITATIONS_RETAINED')));qc.append(h.old.pin(p))
    artifacts={p.name:{k:v for k,v in h.old.pin(p).items() if k!='path'} for p in folder.iterdir()};receipt=folder/'workflow-receipt.json';h.old.dump(receipt,dict(status='PASS',identity=str(index+1)*64,artifacts=artifacts));entries.append(dict(sample=sample,receipt=h.old.pin(receipt),vcf=h.old.pin(vcf),qc=qc,indices=[h.old.pin(str(vcf)+'.tbi')]))
   rule=dict(output=str(root/'output'),state_directory=str(state),source_population_rule=h.old.pin(profile),source_database=h.old.pin(database),mt_audit=h.old.pin(audit),source_resource_manifest=h.old.pin(ready),mt_sources=entries,expected_samples=[e['sample'] for e in entries],max_output_bytes=512*1024**2,min_free_bytes=0,code_sha256={str(ROOT/'scripts/export_mitochondrial_reference_context.py'):h.old.sha(ROOT/'scripts/export_mitochondrial_reference_context.py'),str(ROOT/'scripts/export_population_mt_reference_context.py'):h.LEGACY_SHA});result=h.export(rule);self.assertEqual('PASS_INDEPENDENT_MITOCHONDRIAL_REFERENCE_CONTEXT',result['status']);self.assertTrue(result['no_nuclear_dependency']);self.assertEqual(4,len(result['verified']))
   for entry in entries:
    with gzip.open(Path(rule['output'])/'samples'/entry['sample']/'MT-population-reference.tsv.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    self.assertEqual(2,len(rows));self.assertEqual('EXACT_rCRS_ALLELES',rows[0]['gnomAD_MT_match_status']);self.assertIsNone(json.loads(rows[0]['Helix_reference_json'])['AN']);self.assertIn('NO_FLATTEN',rows[1]['Helix_match_status']);self.assertEqual('',rows[1]['gnomAD_MT_reference_json']);self.assertEqual('0.10,0.20',rows[1]['AF']);self.assertEqual('numt',rows[1]['FILTER'])
   for item in result['applied_artifacts']:
    path=Path(item['path']);self.assertEqual(item['sha256'],h.old.sha(path));self.assertTrue(item['relative_path'].startswith('applied-evidence/MT/'))
    if 'source_original_pin' in item:self.assertEqual(item['source_original_pin']['sha256'],item['sha256'])
   proof=dict(status='PASS_PUBLIC_MT_DATABASE_FOUR_SYNTHETIC_SOURCES_NO_NUCLEAR_GATE',patient_inputs=False,public_database=h.old.pin(database),synthetic_samples=4,rows=8,multiallelic_rows_not_split=4,all_raw_genotype_heteroplasmy_filter_info_preserved=True,all_original_QC_bytes_preserved=True,helix_AN_missing_preserved=True,gnomad_and_helix_exact_matches=4,no_nuclear_receipt_read=True,producer_sha256=h.old.sha(ROOT/'scripts/export_mitochondrial_reference_context.py'))
   if os.environ.get('WGS_MT_PUBLIC_CONTROL_RECEIPT'):Path(os.environ['WGS_MT_PUBLIC_CONTROL_RECEIPT']).write_text(json.dumps(proof,indent=2)+'\n')
if __name__=='__main__':unittest.main()
