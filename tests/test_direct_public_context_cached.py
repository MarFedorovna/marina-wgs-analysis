import importlib.util
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1]
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
c=module('cached_adapter',ROOT/'scripts/direct_public_context_cached.py');fixtures=module('context_fixture',ROOT/'tests/test_direct_public_context.py')
class Cache(unittest.TestCase):
 def test_all_fields_exact_across_samples_gt_ft_and_ordinals(self):
  original,row=fixtures.fixture();cached=c.Context.__new__(c.Context);cached.__dict__.update(original.__dict__);cached.reset_projection_cache()
  for i in range(1,25):
   r=dict(row,sample='S'+str(i%4),GT='0|1' if i%2 else './.',FT='PASS' if i%3 else 'LowQual',transcript='ENST'+str(1+(i-1)//4)+'.1');expected=c.base.Context.augment(original,r,i);actual=cached.augment(r,i);self.assertEqual(expected,actual)
  self.assertEqual(18,cached.cache_hits);self.assertEqual(6,cached.cache_misses);original.close()
 def test_allele_and_clinvar_identity_change_invalidates(self):
  original,row=fixtures.fixture();cached=c.Context.__new__(c.Context);cached.__dict__.update(original.__dict__);cached.reset_projection_cache()
  for i,changes in enumerate([{},dict(clinvar_allele_id='88'),dict(pos='11'),dict(pos='10'),dict(gene='OTHER',hgnc_id='HGNC:999')],1):
   r=dict(row,**changes);self.assertEqual(c.base.Context.augment(original,r,i),cached.augment(r,i))
  self.assertEqual(5,cached.cache_misses);original.close()
if __name__=='__main__':unittest.main()
