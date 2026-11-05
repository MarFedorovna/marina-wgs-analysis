import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SOURCE=Path(__file__).resolve().parents[1]/'scripts/continue_annotation_reference_v2.py'
spec=importlib.util.spec_from_file_location('conditional_annotation',SOURCE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class ConditionalTests(unittest.TestCase):
    def fixture(self,root):
        base=root/'base-config.json';candidate=root/'candidate-config.json';core=root/'core-manifest.json';sidecar=root/'resource-manifest-public-hcm-v2.json';binding=root/'binding.json';approval=root/'ROOT-approval.json'
        original=dict(annotation=dict(resource_manifest=str(core),cache_version=115),output_root='/public-synthetic-output')
        adapter=dict(name='gene_evidence',format='gene_evidence_tsv',path='/public-synthetic.tsv')
        cfg=copy.deepcopy(original);cfg['annotation'].update(local_adapters=[adapter],resource_manifest=str(sidecar))
        m.save(base,original);m.save(candidate,cfg)
        m.save(core,dict(status='PASS',config_sha256='a'*64,vep=dict(image='pinned-image'),resources={}))
        draft=dict(external_upstream_code={},approval_file=str(approval),annotation_config=dict(path=str(candidate),sha256=m.sha(candidate)))
        m.save(binding,draft)
        return dict(pins=[m.pin(base),m.pin(candidate)],package_pins=[],adapter_pins=[],
            original_annotation_config=str(base),candidate_annotation_config=str(candidate),approved_adapters=[adapter],
            sidecar_manifest=str(sidecar),draft_binding_path=str(binding),draft_binding_sha256=m.sha(binding),draft_binding=draft,
            base_core_manifest=str(core),expected_core_recipe_sha256='a'*64,expected_vep_image='pinned-image',
            resource_additions={'orphanet':dict(status='available_local_supplement',version='July2026')})

    def materializer(self,rule,extra_sidecar=False,extra_binding=False):
        value=m.expected_sidecar(m.load(rule['base_core_manifest']),rule['resource_additions'])
        if extra_sidecar:value['unexpected']='mutation'
        m.save(rule['sidecar_manifest'],value)
        final=copy.deepcopy(rule['draft_binding']);final['external_upstream_code'][rule['sidecar_manifest']]=m.sha(rule['sidecar_manifest'])
        if extra_binding:final['resources']='unapproved change'
        m.save(rule['draft_binding_path'],final)

    def test_synthetic_conditional_materialize_then_approval_of_exact_final_binding(self):
        with tempfile.TemporaryDirectory() as t:
            rule=self.fixture(Path(t));before=Path(rule['base_core_manifest']).read_bytes()
            with patch.object(m,'annotation_processes',return_value=[]):result=m.materialize_and_approve(rule,'f'*64,lambda:self.materializer(rule))
            self.assertEqual(result['approval']['binding_sha256'],m.sha(rule['draft_binding_path']))
            self.assertEqual(result['approval']['ROOT_reviewed_controller_rule_sha256'],'f'*64)
            self.assertEqual(Path(rule['base_core_manifest']).read_bytes(),before)

    def test_unexpected_cfg_or_pinned_source_mutation_refused(self):
        with tempfile.TemporaryDirectory() as t:
            rule=self.fixture(Path(t));path=Path(rule['candidate_annotation_config']);value=m.load(path);value['annotation']['cache_version']=116;m.save(path,value)
            with self.assertRaises(AssertionError):m.materialize_and_approve(rule,'f'*64,lambda:self.fail('Must not materialize'))
            self.assertFalse(Path(rule['draft_binding']['approval_file']).exists())

    def test_unexpected_sidecar_or_binding_digest_mutation_refused(self):
        for sidecar,binding in [(True,False),(False,True)]:
            with self.subTest(sidecar=sidecar,binding=binding),tempfile.TemporaryDirectory() as t:
                rule=self.fixture(Path(t))
                with self.assertRaises(AssertionError):m.materialize_and_approve(rule,'f'*64,lambda:self.materializer(rule,sidecar,binding))
                self.assertFalse(Path(rule['draft_binding']['approval_file']).exists())

    def test_other_annotation_job_prevents_approval(self):
        with tempfile.TemporaryDirectory() as t:
            rule=self.fixture(Path(t))
            with patch.object(m,'annotation_processes',return_value=[1234]):
                with self.assertRaises(AssertionError):m.materialize_and_approve(rule,'f'*64,lambda:self.materializer(rule))
            self.assertFalse(Path(rule['draft_binding']['approval_file']).exists())

    def test_unknown_existing_approval_never_overwritten(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'approval.json';path.write_text('{"binding_sha256":"wrong"}')
            before=path.read_bytes()
            with self.assertRaises(AssertionError):m.exact_approval(path,'a'*64,'b'*64)
            self.assertEqual(path.read_bytes(),before)

    def test_cfg_delta_does_not_allow_ignored_annotation_fields(self):
        base=dict(annotation=dict(resource_manifest='old',cache_version=115))
        candidate=copy.deepcopy(base);candidate['annotation'].update(resource_manifest='new',local_adapters=[],availability_overlay={'ignored':True})
        with self.assertRaises(AssertionError):m.check_cfg_delta(base,candidate,[],'new')

    def test_process_detection_is_identity_based_and_ignores_resource_installer(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            for pid,argv in [(1,['python3','/sealed/scripts/annotation_queue.py','run','--execute']),
                             (2,['python3','/installer/prepare_annotation_resources.py','--execute']),
                             (3,['python3','/sealed/scripts/workflow.py','annotate','--execute'])]:
                p=root/str(pid);p.mkdir();(p/'cmdline').write_bytes(b'\0'.join(s.encode() for s in argv)+b'\0')
            self.assertEqual(m.annotation_processes(root),[1,3])


if __name__=='__main__':unittest.main()
