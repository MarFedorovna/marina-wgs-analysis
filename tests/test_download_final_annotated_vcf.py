import base64
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'transfer'))
import download_final_annotated_vcf as m


class FinalPairTests(unittest.TestCase):
    def plan(self):
        files=[dict(RelativePath='outputs/fully-annotated-vcf-v1/'+name,ServerPath=m.REMOTE+name,SourceReceipt='outputs/fully-annotated-vcf-v1/workflow-receipt.json',DownloadFileName=name,SizeBytes=size,Sha256=digest) for name,size,digest in [(m.NAME,8461702199,'a'*64),(m.NAME+'.tbi',1904867,'b'*64)]]
        manifest=dict(schema_version=1,status='PASS_FULLY_ANNOTATED_VCF_DIRECT',final_delivery_complete=True,original_records_fields_preserved=True,packet_decode_verified=True,direct_core_prefix_verified=True,core33_packets_byte_identical=True,
            enriched_vcf=dict(path=files[0]['ServerPath'],bytes=files[0]['SizeBytes'],sha256=files[0]['Sha256']),index=dict(path=files[1]['ServerPath'],bytes=files[1]['SizeBytes'],sha256=files[1]['Sha256']))
        def capture(name,value):
            raw=json.dumps(value).encode();return dict(path=m.REMOTE+name,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),content_base64=base64.b64encode(raw).decode())
        meta=capture('manifest.json',manifest)
        receipt=dict(schema_version=1,status='PASS',stage='fully-annotated-vcf-v1',identity='c'*64,artifacts={'manifest.json':{k:meta[k] for k in ('bytes','sha256')},**{e['DownloadFileName']:dict(bytes=e['SizeBytes'],sha256=e['Sha256']) for e in files}})
        return dict(schema_version=1,status='ROOT_REVIEWED_FINAL_ANNOTATED_VCF_DOWNLOAD',scientific_status='PASS_FULLY_ANNOTATED_VCF_DIRECT',deadline_utc=m.END.isoformat(),destination='NOT_OPENED',state_directory='NOT_CREATED',pause_file='NOT_OPENED',source_provenance=dict(manifest=meta,receipt=capture('workflow-receipt.json',receipt)),files=files)

    def test_dry_validates_actual_sized_pair_without_auth_or_state(self):
        p=self.plan();result=m.run(p,Path('AUTH_MUST_NOT_BE_READ'))
        self.assertEqual(result['status'],'PLAN_ONLY_NO_AUTH_KEY_NETWORK_OR_FILES');self.assertEqual(result['bytes'],8463607066)

    def test_wrong_pair_fake_pass_and_extended_deadline_refused(self):
        p=self.plan();p['files'][1]['ServerPath']+='wrong'
        with self.assertRaises(ValueError):m.run(p)
        p=self.plan();p['scientific_status']='COMPLETED_VEP_COMPONENT_NOT_FINAL'
        with self.assertRaises(ValueError):m.run(p)
        p=self.plan();p['source_provenance']['manifest']['content_base64']=base64.b64encode(b'{"status":"PASS"}').decode()
        with self.assertRaisesRegex(ValueError,'SHA/size'):m.run(p)
        p=self.plan();p['deadline_utc']='2027-01-01T00:00:00+00:00'
        with self.assertRaisesRegex(ValueError,'deadline'):m.run(p)

    def test_remaining_resume_bytes_and_unrelated_final_preserved(self):
        p=self.plan()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);entry=p['files'][0];part=root/(m.NAME+'.part');part.write_bytes(b'x'*100)
            identity=dict(RelativePath=entry['RelativePath'],RemotePath=entry['ServerPath'],ExpectedBytes=entry['SizeBytes'],ExpectedServerSha256=entry['Sha256'],SourceReceipt=entry['SourceReceipt'])
            part.with_name(part.name+'.source.json').write_text(json.dumps(identity))
            self.assertEqual(m.remaining_bytes(p,root),sum(e['SizeBytes'] for e in p['files'])-100)
            final=root/(m.NAME+'.tbi');final.write_bytes(b'unrelated')
            with self.assertRaisesRegex(ValueError,'refusing overwrite'):m.remaining_bytes(p,root)
            self.assertEqual(final.read_bytes(),b'unrelated')

    def test_space_reserve_checked_before_auth_or_state_creation(self):
        p=self.plan()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);before=set(root.iterdir());required=sum(e['SizeBytes'] for e in p['files'])
            with self.assertRaisesRegex(ValueError,'15GiB'):m.check_space(p,root,free_bytes=m.RESERVE+required-1)
            self.assertEqual(m.check_space(p,root,free_bytes=m.RESERVE+required),required)
            self.assertEqual(set(root.iterdir()),before)


if __name__=='__main__':unittest.main()
