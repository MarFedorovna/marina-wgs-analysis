#!/usr/bin/env python3
"""Direct-packet codec with closed sample-scope metadata projection."""
import hashlib
import importlib.util
from pathlib import Path
path=Path(__file__).resolve().parent/'spliceai_vcf_codec.py'
if path.is_symlink() or hashlib.sha256(path.read_bytes()).hexdigest()!='2711491e5265e442638ddfedc2da0d3e890ea521bb94fda974b02f7d0c4552dd':raise ValueError('Exact preserved first SpliceAI codec required')
spec=importlib.util.spec_from_file_location('preserved_spliceai_codec',path);v1=importlib.util.module_from_spec(spec);spec.loader.exec_module(v1)
for name in ['CODEC','SPLICE_FIELDS','INFO_IDS','NAMESPACED_INFO_IDS','encode_packet','decode_packet','encode_schema','decode_schema','schema_header','decode_header_bundle','chunk_headers','public_bundle','rebuild_readable','headers','strip_overlay']:globals()[name]=getattr(v1,name)
SAMPLE_PACKET_KEYS={'v','allele','rows','evidence','gene_context_keys','scope','sample'}
def select_sample(packet,fields,sample):
 selected=v1.select_sample(packet,fields,sample)
 # Cohort membership contains other individuals' states. It remains available
 # in the joint source only; sample VCFs retain own row GT/FT and public sources.
 return {key:value for key,value in selected.items() if key in SAMPLE_PACKET_KEYS}
