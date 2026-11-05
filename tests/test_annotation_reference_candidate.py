import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location('annotation_reference_review',SCRIPTS/'prepare_annotation_reference_candidate.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class BindingCandidateTests(unittest.TestCase):
    def fixture(self,root):
        target=root/'review';target.mkdir();resource=root/'resources';resource.mkdir()
        base=resource/'resource-manifest.json';base.write_text(json.dumps(dict(status='PASS',config_sha256='a'*64,vep=dict(image='test@sha256:pinned'),resources={})))
        new=resource/'resource-manifest-public-hcm-v2.json';overlay=root/'overlay.json';overlay.write_text('{}')
        binding=target/'binding.json';binding.write_text(json.dumps(dict(approval_file=str(target/'approval.json'),external_upstream_code={},code_package=dict(path=str(root/'runtime/code-manifest.json'),sha256='b'*64))))
        plan=dict(base_manifest=str(base),new_manifest=str(new),new_binding=str(binding),overlay_path=str(overlay),overlay_sha256=m.sha(overlay),expected_core_recipe_sha256='a'*64,expected_vep_image='test@sha256:pinned',resource_additions={'orphanet':dict(status='available_local_supplement',version='July2026')})
        (target/'resource-overlay-plan.json').write_text(json.dumps(plan));return target,base,new,binding
    def test_new_sidecar_preserves_core_and_binds_through_existing_digest_guard(self):
        with tempfile.TemporaryDirectory() as t:
            target,base,new,binding=self.fixture(Path(t));original=base.read_bytes();verified=[]
            with patch.object(m,'module',return_value=types.SimpleNamespace(verify_bindings=lambda value:verified.append(value))):result=m.materialize(target)
            self.assertEqual(base.read_bytes(),original);self.assertEqual(m.load(binding)['external_upstream_code'][str(new)],m.sha(new))
            self.assertTrue(verified);self.assertEqual(result['status'],'FINAL_BINDING_REQUIRES_ROOT_REVIEW_NO_APPROVAL_CREATED');self.assertFalse((target/'approval.json').exists())
    def test_existing_approval_prevents_materialization_or_binding_mutation(self):
        with tempfile.TemporaryDirectory() as t:
            target,base,new,binding=self.fixture(Path(t));(target/'approval.json').write_text('approved');before=binding.read_bytes()
            with self.assertRaises(AssertionError):m.materialize(target)
            self.assertFalse(new.exists());self.assertEqual(binding.read_bytes(),before)


if __name__=='__main__':unittest.main()
