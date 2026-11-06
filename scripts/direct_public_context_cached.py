#!/usr/bin/env python3
"""Bounded per-allele memoization of the preserved direct context projection."""
import hashlib
import importlib.util
from pathlib import Path
path=Path(__file__).resolve().parent/'direct_public_context_adapter.py'
if path.is_symlink() or hashlib.sha256(path.read_bytes()).hexdigest()!='ee6c57534d6e8c215de8f4063788a9c53ba7dfe3af5c2a94308ab6d08f6b760e':raise ValueError('Exact preserved context adapter required')
spec=importlib.util.spec_from_file_location('preserved_direct_context',path);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
APPENDED=base.APPENDED
SPLICE_FIELDS=base.SPLICE_FIELDS
class Context(base.Context):
 def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);self.reset_projection_cache()
 def reset_projection_cache(self):self._cache_allele=None;self._projection_cache={};self.cache_hits=0;self.cache_misses=0
 def augment(self,row,ordinal):
  allele=tuple(str(row.get(f,'')) for f in ['chrom','pos','ref','alt'])
  if allele!=self._cache_allele:self._projection_cache.clear();self._cache_allele=allele
  # These are every other row input read by preserved Context.augment/tissue
  # supplement/SpliceAI projection. GT/FT and sample are never inputs to them.
  identity=tuple(row.get(f) for f in ['hgnc_id','gene','transcript','clinvar_allele_id'])
  if identity in self._projection_cache:
   self.cache_hits+=1;result=dict(self._projection_cache[identity]);result['allref_base_master_row_number']=str(ordinal);return result
  self.cache_misses+=1;result=super().augment(row,ordinal)
  if len(self._projection_cache)<4096:self._projection_cache[identity]=dict(result)
  return result
