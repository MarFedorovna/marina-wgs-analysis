import copy
import csv
import gzip
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).parent))
import test_packet_views_fast as common
sys.path.insert(0,str(common.ROOT/'scripts'))
import export_packet_views_fast_v2 as m
import continue_packet_views_fast_v2 as controller
from test_packet_views_fast_continuation import ContinuationGuards


class V2Bounds(unittest.TestCase):
    def test_twentyfour_accepted_twentyfive_rejected(self):
        fixture=ContinuationGuards();fixture.setUp();self.addCleanup(fixture.doCleanups)
        fixture.rule.update(controller_sha256=controller.sha(controller.__file__),workers=24)
        controller.checked_rule(fixture.path,fixture.write())
        fixture.rule['workers']=25
        with self.assertRaisesRegex(ValueError,'bounded'):controller.checked_rule(fixture.path,fixture.write())
        with self.assertRaisesRegex(ValueError,'twenty-four'):m.export(dict(workers=25,deadline_utc=m.END.isoformat()))


def many_alleles():
    f=common.fixture(Path('.'));rows=[];packets=[];lines=[];baseproof=[];fullproof=[]
    for number in range(30):
        packet=copy.deepcopy(f['packets'][number%3]);pos=(number+1)*10;packet['allele'][1]=pos
        first=len(rows)+1;bd=common.proof.VectorDigest();fd=common.proof.VectorDigest()
        for item in packet['rows']:
            values=item['values'];values[f['fields'].index('pos')]=str(pos);values[f['fields'].index('variant_id')]=common.rank.variant_key(*packet['allele'])
            predictions=json.loads(values[f['fields'].index('spliceai_prediction_json')])
            for p in predictions:p['pos']=str(pos);p['variant_key']=common.rank.variant_key(*packet['allele'])
            values[f['fields'].index('spliceai_prediction_json')]=json.dumps(predictions)
            item['i']=len(rows)+1;rows.append(values);bd.row(values[:len(common.rank.FIELDS)]);fd.row(values)
        descriptor=dict(id=number,first_i=first,last_i=len(rows),row_count=len(packet['rows']),record_count=1)
        baseproof.append(dict(descriptor,canonical_row_sha256=bd.digest.hexdigest()));fullproof.append(dict(descriptor,canonical_row_sha256=fd.digest.hexdigest()))
        parts=f['lines'][number%3].rstrip('\n').split('\t');parts[1]=str(pos);info=m.info_dict(parts[7]);info['WGS_EVIDENCE']=common.codec.encode_packet(packet);parts[7]=m.info_text(info)
        lines.append('\t'.join(parts)+'\n');packets.append(packet)
    f.update(rows=rows,packets=packets,lines=lines,core_schema=common.schema(common.rank.FIELDS,baseproof),schema=common.schema(f['fields'],fullproof))
    return common.project_fixture(f)


@unittest.skipUnless(common.pysam and sys.platform=='linux','Linux 24-worker BGZF proof')
class WorkerInvariance(unittest.TestCase):
    def test_two_vs_twentyfour_workers_identical_decoded_tables_and_sample_packets(self):
        affinity=os.sched_getaffinity(0)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);f=many_alleles();folder=root/'fully';folder.mkdir();vcf=folder/'cohort.fully-annotated.vcf.gz'
            with common.pysam.BGZFile(str(vcf),'wb') as stream:stream.write((f['header']+''.join(f['lines'])).encode())
            common.pysam.tabix_index(str(vcf),preset='vcf')
            core=folder/'core-manifest.json';core.write_text(json.dumps(dict(status='DIRECT_PACKET_BASELINE_PASS',schema=f['core_schema'],no_legacy_receipt_created=True)))
            manifest=folder/'manifest.json';manifest.write_text(json.dumps(dict(schema_version=1,status='PASS_FULLY_ANNOTATED_VCF_DIRECT',schema=f['schema'],projection=f['projection'],samples=f['samples'],
                core_manifest=m.pin(core),enriched_vcf=m.pin(vcf),index=m.pin(str(vcf)+'.tbi'),core33_packets_byte_identical=True,direct_core_prefix_verified=True,packet_decode_verified=True,original_records_fields_preserved=True)))
            outputs=[];results=[]
            try:
                for workers in (2,24):
                    os.sched_setaffinity(0,affinity);output=root/('views-'+str(workers));outputs.append(output)
                    rule=dict(root_reviewed=True,synthetic=True,workers=workers,deadline_utc=m.END.isoformat(),output=str(output),fully_manifest=m.pin(manifest),exporter_sha256=m.sha(m.__file__),
                        approved_appended_fields=m.SPLICE_FIELDS,**{k:m.pin(common.ROOT/'scripts'/v) for k,v in [('proof_module','fast_vcf_evidence_candidate.py'),('codec_module','spliceai_allele_info_codec.py'),('view_policy','export_cohort_sample_views.py'),('ranking_module','annotation_stages.py')]})
                    results.append(m.export(rule))
            finally:os.sched_setaffinity(0,affinity)
            for relative in results[0]['artifacts']:
                if relative.startswith(('masters/','tables/')):
                    def content(path):return gzip.decompress(path.read_bytes()) if str(path).endswith('.gz') else path.read_bytes()
                    self.assertEqual(content(outputs[0]/relative),content(outputs[1]/relative),relative)
            for i,sample in enumerate(f['samples'],1):
                def packets(path):
                    with common.pysam.VariantFile(str(path)) as stream:return [common.allele_codec.decode_packet(r.info['WGS_EVIDENCE']) for r in stream.fetch()]
                relative='vcf/sample_%04d.annotated.vcf.gz'%i
                self.assertEqual(packets(outputs[0]/relative),packets(outputs[1]/relative))
            self.assertEqual(results[0]['derived_batch_merkle_root'],results[1]['derived_batch_merkle_root'])


if __name__=='__main__':unittest.main()
