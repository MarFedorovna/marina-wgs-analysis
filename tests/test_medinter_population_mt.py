"""Meaningful data/coordinate/denominator guards using wholly synthetic references."""
import base64
import gzip
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from medinter_population_mt_sources import allowed,check_format,crc32c,range_ok
from medinter_population_mt_common import IntervalLookup,PublicFasta,population_fields,population_info_projection,encode_population_info,decode_population_info
from medinter_population_mt_dossiers import tsv
from medinter_population_mt_coverage import validate as validate_coverage

class PublicPopulationTests(unittest.TestCase):
 def row(self,info,alt='C'):return ['chr1','11','.','A',alt,'.','AS_VQSR',info]
 def test_missing_denominator_is_NULL_never_zero(self):
  parsed=population_fields(self.row('AC=.;AF=.;nhomalt=.'))
  self.assertIsNone(parsed['AN']);self.assertEqual(parsed['AC'],'.');self.assertEqual(parsed['filter'],'AS_VQSR')
 def test_zero_called_denominator_keeps_missing_AF(self):
  parsed=population_fields(self.row('AC=0;AN=0;AF=.;nhomalt=0'))
  self.assertEqual(parsed['AN'],0);self.assertEqual(parsed['AF'],'.')
 def test_multiallelic_arrays_and_homozygote_bounds(self):
  self.assertEqual(population_fields(self.row('AC=1,2;AN=10;AF=0.1,0.2;nhomalt=0,1','C,G'))['AC'],'1,2')
  for bad in ['AC=1;AN=10;AF=0.1;nhomalt=1','AC=11;AN=10;AF=0.1;nhomalt=0','AC=1;AN=10;AF=NaN;nhomalt=0']:
   with self.subTest(bad=bad),self.assertRaises(ValueError):population_fields(self.row(bad))
 def test_no_autosomal_diploidy_assumption_for_sex_chromosome_counts(self):
  row=self.row('AC=1;AN=100;AF=0.01;nhomalt=1');row[0]='chrX'
  self.assertEqual(population_fields(row)['nhomalt'],'1')
 def test_projection_retains_every_population_and_quality_field(self):
  value='AC=1;AN=10;AF_afr=0.1;nhomalt_nfe=0;AS_QD=4;vep=large;CSQ=duplicate;non_ukb_AF=0.01;grpmax=afr'
  self.assertEqual(population_info_projection(value),'AC=1;AN=10;AF_afr=0.1;nhomalt_nfe=0;AS_QD=4;non_ukb_AF=0.01;grpmax=afr')
 def test_lossless_population_INFO_storage(self):
  original='AC=1;AN=100;AF=0.01;nhomalt=0;AS_QD=3;AF_afr=0.1;non_ukb_AN=20;FLAG'
  self.assertEqual(decode_population_info(encode_population_info(original)),original)
 def test_BED_to_one_based_position_boundaries(self):
  lookup=IntervalLookup([('chr1',10,20),('chr1',30,40)])
  self.assertFalse(lookup.contains('chr1',10));self.assertTrue(lookup.contains('chr1',11));self.assertTrue(lookup.contains('chr1',20));self.assertFalse(lookup.contains('chr1',21));self.assertFalse(lookup.contains('chr2',11))
 def test_public_fasta_across_line_boundaries(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'synthetic.fa';path.write_bytes(b'>chrM\nACGT\nTGCA\n');path.with_name(path.name+'.fai').write_text('chrM\t8\t6\t4\t5\n')
   ref=PublicFasta(path)
   try:self.assertEqual(ref.sequence('chrM',3,4),'GTTG')
   finally:ref.close()
 def test_CRC32C_official_standard_control(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'control';path.write_bytes(b'123456789');self.assertEqual(crc32c(path),base64.b64encode(bytes.fromhex('e3069283')).decode())
 def test_tabular_NCBI_preamble_and_gzip_integrity(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'synthetic.gz'
   with gzip.open(path,'wt') as stream:stream.write('##documentation preamble\n#VariationID\tSCV\n1\tSCV000001.1\n')
   check_format(path,'tsv_gzip')
   path.write_bytes(path.read_bytes()[:-5])
   with self.assertRaises((EOFError,OSError)):check_format(path,'tsv_gzip')
 def test_exact_resume_and_no_unapproved_or_query_endpoints(self):
  self.assertTrue(range_ok('bytes 10-99/100',10,100));self.assertFalse(range_ok('bytes 0-99/100',10,100))
  self.assertFalse(allowed('https://example.org/database'));self.assertFalse(allowed('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz?patient=anything'))
 def test_long_public_ClinVar_field_is_bounded_and_parsed(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'synthetic.tsv';path.write_text('#VariationID\tDescription\n1\t'+'x'*200000+'\n')
   self.assertEqual(len(next(tsv(path,'#VariationID\t'))['Description']),200000)
 def test_coverage_identity_and_missing_values(self):
  expected={'hgnc_id':'HGNC:CONTROL'};gene={'hgnc_id':'HGNC:CONTROL','reference_genome':'GRCh38','coverage':{'exome':[{'pos':10,'mean':None,'median':None,'over_10':None,'over_20':None,'over_30':None}],'genome':[]}}
  self.assertEqual(validate_coverage(gene,expected),{'exome':1,'genome':0});self.assertIsNone(gene['coverage']['exome'][0]['mean'])
  gene['hgnc_id']='HGNC:WRONG'
  with self.assertRaises(ValueError):validate_coverage(gene,expected)
if __name__=='__main__':unittest.main()
