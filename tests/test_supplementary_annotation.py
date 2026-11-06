"""Linux public/synthetic portable evidence and source-preservation controls."""
import argparse
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
import pysam

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('supplementary', ROOT/'scripts/annotate_supplementary_callsets.py')
q = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(q)
CODEC = Path(os.environ.get('SUPPLEMENTARY_CODEC', ROOT/'scripts/embed_vcf_evidence.py'))
c = q.module(CODEC, '_public_test_codec')


class SupplementaryTests(unittest.TestCase):
    def fixture(self, root, empty=False):
        samples = ['PUBLIC_A', 'PUBLIC_B', 'PUBLIC_C', 'PUBLIC_D']
        headers = ['##fileformat=VCFv4.2', '##contig=<ID=chrX,length=156040895>',
                   '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
                   '##FORMAT=<ID=FT,Number=1,Type=String,Description="FT">',
                   '##GATKCommandLine=<CommandLine="PRIVATE_PUBLIC_B_PATH">',
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)]
        source = root/'source.vcf.gz'; annotation = root/'annotation.vcf.gz'; table = root/'table.tsv.gz'
        rows = [] if empty else [
            ['chrX','3000001','source1','A','G','60','PASS','AC=3','GT:FT','1:PASS','0:PASS','0|1:PASS','./.:NoCall'],
            ['chrX','3001001','source2','AC','A','20','PASS','AC=1','GT:FT','0:PASS','0:PASS','0/1:LowGQ','./.:NoCall']]
        with gzip.open(source, 'wt') as stream:
            stream.write('\n'.join(headers)+'\n')
            for row in rows: stream.write('\t'.join(row)+'\n')
        with gzip.open(annotation, 'wt') as stream:
            stream.write('\n'.join(headers[:-1]+['##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP">',headers[-1]])+'\n')
            for row in rows:
                added = list(row); added[2] = 'new_vep_identifier'; added[7] += ';CSQ=G|missense_variant'
                stream.write('\t'.join(added)+'\n')
        fields = ['sample','chrom','pos','ref','alt','transcript','GT','FT','AlphaMissense_score','local_evidence_json']
        with gzip.open(table, 'wt', newline='') as stream:
            writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n'); writer.writeheader()
            for row in rows:
                for sample, raw in zip(samples, row[9:]):
                    gt, ft = raw.split(':')
                    writer.writerow(dict(sample=sample,chrom=row[0],pos=row[1],ref=row[3],alt=row[4],transcript='ENST_PUBLIC.1',
                                         GT=gt,FT=ft,AlphaMissense_score='.42',local_evidence_json='{"note":"public\tUnicode Ω\ncontrol"}'))
        return source, annotation, table, fields

    def test_full_persisted_roundtrip_and_original_fields(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name); src, ann, table, fields = self.fixture(root)
            output = root/'complete.vcf.gz'
            proof = q.embed(src, ann, table, output, c, lambda: None)
            self.assertEqual((proof['records'], proof['rows']), (2,8))
            with gzip.open(src,'rt') as a, gzip.open(output,'rt') as b:
                before=[r.rstrip().split('\t') for r in a if not r.startswith('#')]
                after=[r.rstrip().split('\t') for r in b if not r.startswith('#')]
            for original, final in zip(before,after):
                self.assertEqual(original[:7],final[:7]); self.assertEqual(original[8:],final[8:])
                self.assertTrue(final[7].startswith(original[7]+';'))
            with gzip.open(table,'rt') as a,gzip.open(proof['extracted_master']['path'],'rt') as b:
                self.assertEqual(list(csv.reader(a,delimiter='\t')),list(csv.reader(b,delimiter='\t')))
            with pysam.VariantFile(str(output)) as vcf:
                self.assertEqual(sum(1 for _ in vcf.fetch('chrX',3000000,3001002)),2)

    def test_sample_packet_and_header_privacy(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name); src, ann, table, fields = self.fixture(root)
            cohort=root/'cohort.vcf.gz'; q.embed(src,ann,table,cohort,c,lambda:None)
            sample=root/'sample.vcf.gz'
            with gzip.open(cohort,'rt') as stream,pysam.BGZFile(str(sample),'wb') as out:
                for line in stream:
                    if line.startswith('#CHROM'):line='\t'.join(line.rstrip().split('\t')[:10])+'\n'
                    elif not line.startswith('#'):line='\t'.join(line.rstrip().split('\t')[:10])+'\n'
                    out.write(line.encode())
            pysam.tabix_index(str(sample),preset='vcf')
            q.prune_sample(sample,c,fields,lambda:None,['PUBLIC_A','PUBLIC_B','PUBLIC_C','PUBLIC_D'])
            with gzip.open(sample,'rt') as stream:
                content=stream.read()
            self.assertNotIn('PUBLIC_B',content);self.assertNotIn('GATKCommandLine',content)
            for line in content.splitlines():
                if line.startswith('#'):continue
                token=next(t.split('=',1)[1] for t in line.split('\t')[7].split(';') if t.startswith(q.TAG+'='))
                packet=c.decode_packet(token)
                self.assertEqual(len(packet['rows']),1)
                self.assertTrue(all(r['values'][0]=='PUBLIC_A' for r in packet['rows']))

    def test_empty_input_is_explicit_empty_not_failure(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);src,ann,table,_=self.fixture(root,empty=True)
            proof=q.embed(src,ann,table,root/'empty.vcf.gz',c,lambda:None)
            self.assertEqual((proof['records'],proof['rows']),(0,0))

    def test_missing_master_row_fails(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);src,ann,table,_=self.fixture(root)
            with gzip.open(table,'wt') as stream:stream.write('sample\tchrom\tpos\tref\talt\n')
            with self.assertRaises(ValueError):q.embed(src,ann,table,root/'bad.vcf.gz',c,lambda:None)

    def test_budget_configuration_unchanged(self):
        base={'budget':{'planned_vm_usd':40},'annotation':{'forks':4},'resources':{}}
        modified=q.configure(base,'source','new-output')
        self.assertEqual(base['budget'],modified['budget']);self.assertEqual(base['annotation']['forks'],4)
        self.assertEqual(modified['annotation']['forks'],8)

    def test_deadline_extension_rejected(self):
        with self.assertRaises(ValueError):q.prepare(argparse.Namespace(deadline_utc='2026-11-06T07:00:00Z'))

    def test_upstream_path_escape_rejected(self):
        with self.assertRaises(ValueError):q.upstream_pin({'artifacts':{'source.vcf.gz':{'path':'/outside/source.vcf.gz','bytes':0,'sha256':'0'*64}}},Path('/reviewed'),'source.vcf.gz')

    def test_actual_list_artifact_caller_contract_and_canonical_profile_identity(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);caller=root/'caller';caller.mkdir()
            executor=root/'executor.py';executor.write_text('# public fixture\n')
            cfg=dict(executor=q.pin(executor),processing_model='PUBLIC_CONTROL_ONLY')
            reviewed=root/'ROOT-profile.json';q.write(reviewed,cfg)
            stored=caller/'reviewed-profile.json';q.write(stored,cfg)
            vcf=caller/'cohort.xy.soft-filtered.vcf.gz';vcf.write_bytes(b'public fixture')
            index=Path(str(vcf)+'.tbi');index.write_bytes(b'public index fixture')
            identity=q.hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest()
            artifacts=[q.pin(vcf),q.pin(index)]
            q.write(caller/'manifest.json',dict(status='PASS_XY_GERMLINE_CALLING',profile_sha256=identity,
                    producer_sha256=cfg['executor']['sha256'],source_profile=q.pin(stored),artifacts=artifacts))
            q.write(caller/'workflow-receipt.json',dict(status='PASS',stage='xy-germline-v1',
                    artifacts=artifacts+[q.pin(caller/'manifest.json')]))
            p=dict(upstream_output=str(caller),upstream_stage='xy-germline-v1',upstream_manifest_status='PASS_XY_GERMLINE_CALLING',
                   upstream_profile=q.pin(reviewed),upstream_profile_identity=identity,upstream_executor=cfg['executor'],
                   upstream_vcf_name=vcf.name)
            self.assertEqual(q.checked_upstream(p)['vcf'],q.pin(vcf))
            p['upstream_profile_identity']=q.sha(reviewed)
            with self.assertRaises(ValueError):q.checked_upstream(p)

    def test_postprune_manifest_and_receipt_pin_every_actual_artifact(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);src,ann,table,fields=self.fixture(root)
            cohort=root/'cohort.vcf.gz';q.embed(src,ann,table,cohort,c,lambda:None)
            views=root/'views';(views/'vcf').mkdir(parents=True)
            sample=views/'vcf/sample.vcf.gz'
            with gzip.open(cohort,'rt') as stream,pysam.BGZFile(str(sample),'wb') as out:
                for line in stream:
                    if line.startswith('#CHROM') or not line.startswith('#'):
                        line='\t'.join(line.rstrip().split('\t')[:10])+'\n'
                    out.write(line.encode())
            pysam.tabix_index(str(sample),preset='vcf')
            old=q.pin(sample)
            q.write(views/'manifest.json',dict(status='PASS',artifacts={'vcf/sample.vcf.gz':old}))
            q.prune_sample(sample,c,fields,lambda:None,['PUBLIC_A','PUBLIC_B','PUBLIC_C','PUBLIC_D'])
            self.assertNotEqual(old['sha256'],q.sha(sample))
            q.finalize_views(views,dict(public_control=True),lambda:None)
            manifest=q.load(views/'manifest.json');receipt=q.load(views/'workflow-receipt.json')
            for relative,item in manifest['artifacts'].items():
                actual=q.pin(views/relative)
                self.assertEqual(actual['bytes'],item['bytes']);self.assertEqual(actual['sha256'],item['sha256'])
                self.assertEqual(actual,receipt['artifacts'][relative])
            self.assertEqual(receipt['artifacts']['manifest.json'],q.pin(views/'manifest.json'))
            self.assertTrue(manifest['packet_sample_privacy_verified'])


if __name__=='__main__':unittest.main()
