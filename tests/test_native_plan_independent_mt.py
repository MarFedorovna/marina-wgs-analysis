import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).parent))
from test_native_plan_sidecars import NativePlanSidecarTests as Harness,m


class IndependentMTPlanTests(unittest.TestCase):
    setUp=Harness.setUp;file=Harness.file;json=Harness.json;pin=Harness.pin
    sidecars=Harness.sidecars;primary_pins=Harness.primary_pins;inventory=Harness.inventory

    def mt(self):
        # The scientific exporter tests the actual public SQLite release. This
        # planner fixture substitutes only its known metadata release digest.
        ready=self.root/'resources/public-ready.json';self.file(ready,b'synthetic public release metadata')
        old_ready='465ee667cb78baddc3665d06adb15fb27ad54c47415ff14238fda5f0f42c3d85'
        self.helper=m.HELPER.replace(old_ready,self.pin(ready)['sha256'])
        audit=self.root/'private/current-MT-audit.json';self.json(audit,dict(status='PASS_RESEARCH_WITH_LIMITATIONS'))
        source_rule=self.root/'private/reviewed-original-population-rule.json';self.json(source_rule,dict(review_status='ROOT_REVIEWED_EXACT_LOCAL_ONLY_POPULATION_RULE',no_network=True))
        package=self.root/'private/mitochondrial-reference-context-candidate-v1/scripts'
        legacy=Path(__file__).resolve().parents[1]/'scripts/export_population_mt_reference_context.py'
        self.file(package/'export_population_mt_reference_context.py',legacy.read_bytes());self.file(package/'export_mitochondrial_reference_context.py',b'synthetic reviewed wrapper identity')
        folder=self.root/'outputs/mitochondrial-reference-context-v1';inputs=[];applied=[];artifacts={}
        def output(relative,raw,role,receipts,original=None):
            artifacts[relative]=self.file(folder/relative,raw)
            item=dict(path=str(folder/relative),relative_path='applied-evidence/MT/'+relative,**artifacts[relative],role=role,source_receipts=receipts)
            if original:item['source_original_pin']=original
            applied.append(item)
        for i in range(4):
            sample='SYNTHETIC_MT_'+str(i);source=self.root/'outputs/mitochondria-recovery-v2'/sample/'mt/calls'
            names=['mitochondria.filtered.vcf.gz','mitochondria.filtered.vcf.gz.tbi','research-report.json']
            outputs={name:self.file(source/name,b'unchanged synthetic original '+name.encode()) for name in names}
            self.json(source/'workflow-receipt.json',dict(status='PASS',identity=str(i+1)*64,artifacts=outputs))
            entry=dict(sample=sample,receipt=self.pin(source/'workflow-receipt.json'),vcf=self.pin(source/names[0]),qc=[self.pin(source/names[2])],indices=[self.pin(source/names[1])]);inputs.append(entry)
            output('samples/'+sample+'/MT-population-reference.tsv.gz',b'applied MT source context','MT_POPULATION_CONTEXT',[entry['receipt']])
            for original in [entry['vcf'],*entry['indices']]:output('samples/'+sample+'/'+Path(original['path']).name,Path(original['path']).read_bytes(),'MT_SOURCE_VCF',[entry['receipt']],original)
            output('samples/'+sample+'/MT-QC/00-research-report.json.txt',Path(entry['qc'][0]['path']).read_bytes(),'MT_SOURCE_QC',[entry['receipt']],entry['qc'][0])
        receipt_pins=[x['receipt'] for x in inputs]
        output('cohort/MT-population-reference-long.tsv.gz',b'four samples actual applied long table fixture','MT_POPULATION_CONTEXT',receipt_pins)
        output('index.html',b'local MT report','MT_POPULATION_CONTEXT',receipt_pins)
        meta=dict(schema_version=1,status='PASS_INDEPENDENT_MITOCHONDRIAL_REFERENCE_CONTEXT',source_population_rule=self.pin(source_rule),resource_ready_manifest=self.pin(ready),
            source_database=dict(path=str(self.root/'resources/NEVER_COPY.sqlite'),bytes=1,sha256='f069164654316043d8a96f47ad5b29413d7a28b0b9ba320648206110e42739b0'),mt_audit=self.pin(audit),
            source_mt_inputs=inputs,source_mt_receipts=receipt_pins,source_mt_vcfs=[x['vcf'] for x in inputs],samples=[x['sample'] for x in inputs],applied_artifacts=applied,
            producer_code_sha256={str(p):self.pin(p)['sha256'] for p in package.glob('*.py')},network_requests=0,no_nuclear_dependency=True,primary_and_original_MT_changed=False,no_automatic_ACMG=True,clinical_validation=False)
        self.json(folder/'manifest.json',meta);artifacts['manifest.json']={k:self.pin(folder/'manifest.json')[k] for k in ('bytes','sha256')}
        self.json(folder/'workflow-receipt.json',dict(schema_version=1,status='PASS',stage='mitochondrial-reference-context-v1',identity='7'*64,artifacts=artifacts))
        return folder

    def test_completed_independent_MT_context_originals_and_QC_selected_no_database(self):
        folder=self.mt()
        with patch.object(m,'HELPER',self.helper):plan,_=self.inventory()
        paths={e['RelativePath'] for e in plan['Files'] if '/mitochondrial-reference-context-v1/' in e['RelativePath']}
        self.assertEqual(sum('filtered.vcf.gz' in p for p in paths),8)
        self.assertEqual(sum('/MT-QC/' in p for p in paths),4)
        self.assertTrue(any('cohort/MT-population-reference-long.tsv.gz' in p for p in paths))
        self.assertFalse(any('resources/' in p for p in paths))

    def test_failure_and_missing_source_copy_refused(self):
        folder=self.mt();meta=json.loads((folder/'manifest.json').read_text());meta['applied_artifacts']=[e for e in meta['applied_artifacts'] if '/MT-QC/' not in e['relative_path']]
        self.json(folder/'manifest.json',meta);receipt=json.loads((folder/'workflow-receipt.json').read_text());receipt['artifacts']['manifest.json']={k:self.pin(folder/'manifest.json')[k] for k in ('bytes','sha256')};self.json(folder/'workflow-receipt.json',receipt)
        with patch.object(m,'HELPER',self.helper):
            with self.assertRaisesRegex(RuntimeError,'omitted'):self.inventory()


if __name__=='__main__':unittest.main()
