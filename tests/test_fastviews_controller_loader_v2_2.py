"""Real controller -> dynamically loaded frozen exporter -> 24 fork workers.

Only the production path validator is replaced by a strict synthetic fixture
validator. The exporter, its importlib loader, pool, workers and output checks
execute unchanged; no patient files or production output paths are accessed.
"""
import csv
import datetime as dt
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(Path(__file__).parent))
import continue_packet_views_fast_v2_2 as controller
import test_packet_views_fast as common
import test_packet_views_fast_v2 as fixtures
import test_packet_views_parent_limits as guards
guards.controller=controller


class NewControllerParentGuards(guards.ParentLimits):
    """Reuse the three actual parent timer/memory tests against v2.2."""


@unittest.skipUnless(common.pysam and sys.platform=='linux','Real Linux 24-worker controller regression')
class ControllerPoolRegression(unittest.TestCase):
    def test_real_importlib_controller_and_twentyfour_worker_export(self):
        with tempfile.TemporaryDirectory(prefix='PUBLIC-controller-loader-proof-') as temp:
            root=Path(temp);f=fixtures.many_alleles();folder=root/'fully';folder.mkdir();vcf=folder/'cohort.fully-annotated.vcf.gz'
            with common.pysam.BGZFile(str(vcf),'wb') as stream:stream.write((f['header']+''.join(f['lines'])).encode())
            common.pysam.tabix_index(str(vcf),preset='vcf')
            pin=common.m.pin
            core=folder/'core-manifest.json';core.write_text(json.dumps(dict(status='DIRECT_PACKET_BASELINE_PASS',schema=f['core_schema'],no_legacy_receipt_created=True)))
            manifest=folder/'manifest.json';manifest.write_text(json.dumps(dict(schema_version=1,status='PASS_FULLY_ANNOTATED_VCF_DIRECT',schema=f['schema'],projection=f['projection'],samples=f['samples'],
                core_manifest=pin(core),enriched_vcf=pin(vcf),index=pin(str(vcf)+'.tbi'),core33_packets_byte_identical=True,direct_core_prefix_verified=True,packet_decode_verified=True,original_records_fields_preserved=True)))
            exporter=ROOT/'scripts/export_packet_views_fast_v2.py'
            self.assertEqual(pin(exporter)['sha256'],'078726c0d3ec9fc9ebf4f463c213913592092ce8254db18611befef49e630ba9')
            output=root/'cohort-sample-views-v1';state=root/'public-state'
            rule=dict(schema_version=1,root_reviewed=True,synthetic=False,workers=24,max_hours=.1,interval_seconds=10,deadline_utc=controller.END.isoformat(),
                fully_manifest_path=str(manifest),output=str(output),state_directory=str(state),controller_sha256=controller.sha(controller.__file__),
                exporter=pin(exporter),approved_appended_fields=common.m.SPLICE_FIELDS,authority_expected={'core.no_legacy_receipt_created':True},
                **{name:pin(ROOT/'scripts'/file) for name,file in [('proof_module','fast_vcf_evidence_candidate.py'),('codec_module','spliceai_allele_info_codec.py'),('view_policy','export_cohort_sample_views.py'),('ranking_module','annotation_stages.py')]})
            path=root/'public-rule.json';path.write_text(json.dumps(rule));digest=controller.sha(path)
            original_validator=controller.checked_rule
            def fixture_validator(actual_path,actual_digest):
                self.assertEqual(Path(actual_path),path);self.assertEqual(actual_digest,digest);self.assertEqual(controller.sha(path),digest)
                self.assertEqual(controller.sha(controller.__file__),rule['controller_sha256'])
                self.assertTrue(all(s.startswith('SYNTHETIC_') for s in f['samples']))
                self.assertTrue(root in output.parents and root in state.parents and root in manifest.parents)
                for name in ('exporter','proof_module','codec_module','view_policy','ranking_module'):
                    self.assertEqual(pin(rule[name]['path']),rule[name])
                return controller.load(path)
            controller.checked_rule=fixture_validator
            try:result=controller.continue_views(path,digest,execute=True,watch=False)
            finally:controller.checked_rule=original_validator
            self.assertEqual(result['status'],'PASS_VCF_FIRST_COHORT_SAMPLE_VIEWS_PRIMARY_UNCHANGED')
            loaded=sys.modules['reviewed_fastviews_exporter']
            self.assertEqual(loaded.process_batch.__module__,'reviewed_fastviews_exporter')
            self.assertIs(sys.modules[loaded.process_batch.__module__].process_batch,loaded.process_batch)
            completed=json.loads((output/'manifest.json').read_text());receipt=json.loads((output/'workflow-receipt.json').read_text())
            self.assertEqual(completed['execution_bounds']['workers'],24);self.assertEqual(completed['vcf_records'],30)
            self.assertEqual(completed['master_transcript_rows'],240);self.assertTrue(completed['sample_packet_scope_verified'])
            self.assertEqual(receipt['status'],'PASS');self.assertFalse((state/'owns.lock').exists())
            with gzip.open(output/'masters/cohort.master.tsv.gz','rt',newline='') as stream:
                table=list(csv.reader(stream,delimiter='\t'));self.assertEqual(table[1:],f['rows'])
            for index,sample in enumerate(f['samples'],1):
                with common.pysam.VariantFile(str(output/'vcf'/('sample_%04d.annotated.vcf.gz'%index))) as stream:
                    for record in stream.fetch():
                        packet=common.allele_codec.decode_packet(record.info['WGS_EVIDENCE'])
                        self.assertEqual(packet['sample'],sample)
                        self.assertTrue(all(r['values'][f['fields'].index('sample')]==sample for r in packet['rows']))


if __name__=='__main__':unittest.main()
