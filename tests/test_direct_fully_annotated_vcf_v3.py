"""Synthetic resume-contract checks; no real VCF, clinical IDs or databases."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import types
import unittest

spec = importlib.util.spec_from_file_location('resume_driver', Path(__file__).resolve().parents[1] / 'scripts/direct_fully_annotated_vcf_v3.py')
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)


class ResumeContractTests(unittest.TestCase):
    def fixture(self, root):
        eof=b'BGZF_EOF_SYNTHETIC';d.C=types.SimpleNamespace(BGZF_EOF=eof)
        body=root/'body.vcf.gz';body.write_bytes(b'unchanged synthetic compressed-body placeholder'+eof)
        job=dict(id=3,first_i=41,last_i=48,expected_rows=8,records=2)
        proof=dict(id=3,first_i=41,last_i=48,row_count=8,record_count=2,
            canonical_row_sha256='a'*64,body=d.pin(body))
        data={'batch-proof.json':proof,'full-vector-proof.json':dict(canonical_row_sha256='b'*64,plaintext_bytes=100),
            'annotation_counts.json':dict(records=2,transcript_sample_rows=8),
            'consequence_counts.json':dict(schema_version=1,definitions={},samples={},cohort={})}
        for name,value in data.items():(root/name).write_text(json.dumps(value))
        return job

    def test_complete_cache_preserves_body_and_vector_commitments(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);job=self.fixture(root)
            before={p.name:p.read_bytes() for p in root.iterdir()}
            result=d.complete_cache(root,job,1024)
            self.assertEqual(result['proof']['canonical_row_sha256'],'a'*64)
            self.assertEqual(result['full_proof']['canonical_row_sha256'],'b'*64)
            self.assertEqual({p.name:p.read_bytes() for p in root.iterdir()},before)

    def test_missing_commit_marker_requires_recompute_not_partial_reuse(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);job=self.fixture(root);(root/'full-vector-proof.json').unlink()
            self.assertIsNone(d.complete_cache(root,job,1024))

    def test_changed_body_wrong_ordinals_or_plaintext_ceiling_fail_closed(self):
        for change in ('body','ordinal','cap','eof'):
            with self.subTest(change=change),tempfile.TemporaryDirectory() as temp:
                root=Path(temp);job=self.fixture(root);limit=1024
                if change=='body':(root/'body.vcf.gz').write_bytes(b'changed')
                if change=='ordinal':job['first_i']=42
                if change=='cap':limit=99
                if change=='eof':
                    body=root/'body.vcf.gz';body.write_bytes(b'body without EOF')
                    p=json.loads((root/'batch-proof.json').read_text());p['body']=d.pin(body)
                    (root/'batch-proof.json').write_text(json.dumps(p))
                with self.assertRaises((ValueError,OSError)):
                    d.complete_cache(root,job,limit)


if __name__=='__main__':unittest.main()
