#!/usr/bin/env python3
"""Versioned portable context codec extending the immutable reviewed base codec."""
import hashlib
import importlib.util
from pathlib import Path

BASE_SHA='18e7a2102aceb5b962892242763d7d331523d7669c93f16b12f9fe71a1eb238e'
path=Path(__file__).resolve().parent/'embed_vcf_evidence.py'
if path.is_symlink() or hashlib.sha256(path.read_bytes()).hexdigest()!=BASE_SHA:raise ValueError('Exact preserved base evidence codec required')
spec=importlib.util.spec_from_file_location('preserved_portable_evidence_codec',path);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
path=Path(__file__).resolve().parent/'spliceai_packet_context.py';spec=importlib.util.spec_from_file_location('spliceai_packet_projection',path);splice=importlib.util.module_from_spec(spec);spec.loader.exec_module(splice)
CODEC=base.CODEC
SPLICE_FIELDS=splice.SPLICE_FIELDS
INFO_IDS=[*base.INFO_IDS,'SpliceAI','WGS_SPLICEAI_STATUS','WGS_SPLICEAI_SOURCE']
NAMESPACED_INFO_IDS=INFO_IDS
encode_packet=base.encode_packet
decode_packet=base.decode_packet
encode_schema=base.encode_schema
decode_schema=base.decode_schema
schema_header=base.schema_header
decode_header_bundle=base.decode_header_bundle
chunk_headers=base.chunk_headers
public_bundle=base.public_bundle
select_sample=base.select_sample
def rebuild_readable(packet,fields):return dict(base.rebuild_readable(packet,fields),**splice.readable(packet,fields))
def headers(schema,coverage_token,public_token):
 return base.headers(schema,coverage_token,public_token)+[
  b'##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI1.3.1 five official models; GRCh38 GENCODE24 gene-union scope NOT VEP GENCODE49 transcript mapping; distance500 mask0; ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL; unsupported/outside-scope scores missing">\n',
  b'##INFO=<ID=WGS_SPLICEAI_STATUS,Number=.,Type=String,Description="Percent-encoded sample|gene|status|model-scope; missing scores never zero">\n',
  b'##INFO=<ID=WGS_SPLICEAI_SOURCE,Number=.,Type=String,Description="Percent-encoded full model/version/hash/FASTA/distance/mask source JSON; gene-level score not reassigned to unverified transcript">\n']
def strip_overlay(line):
 ending=b'\r\n' if line.endswith(b'\r\n') else b'\n' if line.endswith(b'\n') else b'';parts=(line[:-len(ending)] if ending else line).split(b'\t');parts[7]=b';'.join(p for p in parts[7].split(b';') if p.split(b'=',1)[0].decode() not in INFO_IDS) or b'.';return b'\t'.join(parts)+ending
