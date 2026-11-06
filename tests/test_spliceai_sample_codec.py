import importlib.util
import json
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1];spec=importlib.util.spec_from_file_location('sample_safe_codec',ROOT/'scripts/spliceai_vcf_codec_v2.py');c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
class SampleScope(unittest.TestCase):
 def test_other_genotypes_membership_and_unknown_private_metadata_removed(self):
  fields=['sample','GT','FT','gene','transcript','spliceai_prediction_json','spliceai_status'];packet=dict(v=1,allele=['chr1',10,'A','G'],rows=[dict(i=1,values=['S1','0|1','PASS','GENE','TX.1','[]','OUTSIDE']),dict(i=2,values=['S2','1/1','PASS','OTHER','TX.2','[]','OUTSIDE'])],membership=dict(sample_states_json=json.dumps({'S1':'alternate','S2':'alternate'}),quality_passing_carriers='S1&S2',carrier_count=2),future_private_annotation={'S2':'private'},evidence={'SCV':[{'public_reference_only':True}]},gene_context_keys=['HGNC:1'],scope='cohort',sample=None)
  selected=c.select_sample(packet,fields,'S1');self.assertNotIn('membership',selected);self.assertNotIn('future_private_annotation',selected);self.assertNotIn('S2',json.dumps(selected));self.assertEqual('0|1',selected['rows'][0]['values'][1]);self.assertEqual(packet['evidence'],selected['evidence']);self.assertIn('S2',json.dumps(packet));self.assertEqual(selected,c.decode_packet(c.encode_packet(selected)));self.assertNotIn('S2',str(c.rebuild_readable(selected,fields)))
 def test_cohort_projection_apis_and_info_contract_unchanged(self):
  self.assertEqual(c.v1.INFO_IDS,c.INFO_IDS);self.assertIn('SpliceAI',c.INFO_IDS);self.assertEqual(13,len(c.SPLICE_FIELDS));schema=dict(v=1,scope='cohort');self.assertEqual(schema,c.decode_schema(c.encode_schema(schema)))
if __name__=='__main__':unittest.main()
