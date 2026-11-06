import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
def imported(name,path):
    spec=importlib.util.spec_from_file_location(name,path);value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value
m=imported('vcf_views',ROOT/'scripts/export_cohort_sample_views_v3.py')
codec=imported('packet_codec',ROOT/'scripts/embed_vcf_evidence.py')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import test_cohort_sample_views as baseline
import test_cohort_sample_views_v2 as joined_test


def pin(path):
    path=Path(path);return dict(path=str(path),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def fully_fixture(fixture):
    joined_path=joined_test.joined(fixture.workflow);joined_test.applied(joined_path);joined=json.loads(joined_path.read_text())
    control_fields=joined['original_fields']+joined['appended_fields'];fields=control_fields+m.SPLICE_FIELDS;rows=[]
    with gzip.open(joined['enriched_master']['path'],'rt',newline='') as stream:
        reader=csv.reader(stream,delimiter='\t');assert next(reader)==control_fields
        for i,row in enumerate(reader,1):rows.append(dict(i=i,values=row+splice_values(dict(zip(control_fields,row)))))
    full=m.DigestText(fields)
    for row in rows:full.row(row['values'])
    schema=dict(v=1,codec=m.CODEC,scope='cohort',sample=None,original_fields=joined['original_fields'],appended_fields=joined['appended_fields']+m.SPLICE_FIELDS,
        full_fields=fields,row_count=len(rows),invariant_digest=joined['invariant_digest'],full_fields_digest=full.digest.hexdigest())
    folder=fixture.workflow.parent/'fully-annotated-vcf-v1';folder.mkdir()
    out=folder/'cohort.fully-annotated.vcf.gz'
    with gzip.open(fixture.stage6/'annotated.vcf.gz','rb') as inp,baseline.pysam.BGZFile(str(out),'wb') as dst:
        for line in inp:
            if line.startswith(b'#CHROM'):
                for header in codec.headers(schema,codec.encode_packet(dict(v=1,scope='PUBLIC_GENE_CONTEXT',genes={})),codec.encode_packet(dict(v=1,scope='PUBLIC_INSTALLED_REFERENCE_OBJECTS',sources=[]))):dst.write(header)
            if line.startswith(b'#'):dst.write(line);continue
            values=line.decode().rstrip('\r\n').split('\t');allele=[values[0],int(values[1]),values[3],values[4]]
            selected=[row for row in rows if int(row['values'][fields.index('pos')])==allele[1]]
            packet=dict(v=1,allele=allele,rows=selected,evidence=dict(population=[dict(dataset='exomes',source_version='gnomAD4.1.1',AF='0.001',AC='1',AN='1000',nhomalt='0',match_status='EXACT',population_FILTER='PASS')]),gene_context_keys=[],scope='cohort',sample=None)
            tags=dict(WGS_EVIDENCE=codec.encode_packet(packet),**codec.rebuild_readable(packet,fields));values[7]+=';'+m.info_text(tags)
            dst.write(('\t'.join(values)+'\n').encode())
    baseline.pysam.tabix_index(str(out),preset='vcf')
    manifest=dict(schema_version=1,status='PASS_FULLY_ANNOTATED_VCF',base_vcf=pin(fixture.stage6/'annotated.vcf.gz'),master_manifest=pin(joined_path),master=joined['enriched_master'],
        combined_receipt=pin(joined_path.parent/'workflow-receipt.json'),enriched_vcf=pin(out),index=pin(Path(str(out)+'.tbi')),schema=schema,
        original_fields=joined['original_fields'],appended_fields=joined['appended_fields']+m.SPLICE_FIELDS,full_fields=fields,row_count=len(rows),record_count=6,
        invariant_digest=joined['invariant_digest'],full_fields_digest=full.digest.hexdigest(),original_records_fields_preserved=True,packet_decode_verified=True,decoded_master_vectors_verified=True,indexed_seek_checks=6)
    splice=fixture.workflow.parent/'spliceai-reference-context-v1';splice.mkdir()
    prediction=splice/'predictions.tsv.gz'
    with gzip.open(prediction,'wt') as stream:stream.write('synthetic model predictions for schema test only\n')
    (splice/'manifest.json').write_text(json.dumps(dict(schema_version=1,status='PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT',transcript_scope='GENCODE24_GENE_UNION_NOT_ENST49_MAPPING')))
    (splice/'workflow-receipt.json').write_text(json.dumps(dict(schema_version=1,status='PASS',stage='spliceai-reference-context-v1',identity='4'*64,
        artifacts={file.name:{k:pin(file)[k] for k in ('bytes','sha256')} for file in [splice/'manifest.json',prediction]})))
    manifest.update(master_role='ORIGINAL_PREFIX_CONTROL',control_prefix_fields=control_fields,supplement_fields=m.SPLICE_FIELDS,
        spliceai_manifest=pin(splice/'manifest.json'),spliceai_receipt=pin(splice/'workflow-receipt.json'),spliceai_predictions=pin(prediction))
    path=folder/'manifest.json';path.write_text(json.dumps(manifest))
    (folder/'workflow-receipt.json').write_text(json.dumps(dict(schema_version=1,status='PASS',stage='fully-annotated-vcf-v1',identity='3'*64,
        artifacts={file.name:{k:pin(file)[k] for k in ('bytes','sha256')} for file in [path,out,Path(str(out)+'.tbi')]})))
    return path,joined_path


def splice_values(row):
    return ['SCORED',row['gene'],'0.00','0.15','0.90','0.03','0','-12','4','7','[{"DS_AG":0,"DS_DG":0.9}]',
        '{"release":"1.3.1","annotation":"GENCODE24"}','GENCODE24_GENE_UNION_NOT_ENST49_MAPPING']


class VcfPacketPureTests(unittest.TestCase):
    def test_sample_selection_retains_source_ordinals_and_public_allele_evidence(self):
        packet=dict(v=1,allele=['chr1',1,'A','T'],rows=[dict(i=1,values=['A','0/1']),dict(i=2,values=['B','1/.']),dict(i=5,values=['A','0/1'])],
            evidence=dict(population=[dict(dataset='exomes',AF='.1')]),gene_context_keys=['HGNC:1'])
        selected=codec.select_sample(packet,['sample','GT'],'A')
        self.assertEqual([r['i'] for r in selected['rows']],[1,5]);self.assertEqual(selected['sample'],'A');self.assertEqual(selected['evidence'],packet['evidence'])
        self.assertEqual(codec.decode_packet(codec.encode_packet(selected)),selected)


@unittest.skipUnless(baseline.pysam and sys.platform=='linux','Linux pysam required for actual BGZF VCF-first integration')
class VcfPacketIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.f=baseline.ViewIntegrationTests();self.f.setUp();self.addCleanup(self.f.doCleanups)
        self.path,self.joined=fully_fixture(self.f)
        root=ROOT/'scripts';self.rule=dict(workflow_root=str(self.f.workflow),output=str(self.f.workflow.parent/'cohort-sample-views-v1'),state_directory=str(self.f.root/'private-state'),
            ranking_module=str(root/'annotation_stages.py'),ranking_sha256=pin(root/'annotation_stages.py')['sha256'],
            v1_engine=str(root/'export_cohort_sample_views.py'),v1_engine_sha256=pin(root/'export_cohort_sample_views.py')['sha256'],
            v2_helper=str(root/'export_cohort_sample_views_v2.py'),v2_helper_sha256=pin(root/'export_cohort_sample_views_v2.py')['sha256'],
            codec_module=str(root/'embed_vcf_evidence.py'),codec_sha256=pin(root/'embed_vcf_evidence.py')['sha256'],
            fully_vcf_manifest=str(self.path),fully_manifest_status='PASS_FULLY_ANNOTATED_VCF',control_joined_manifest=str(self.joined),packet_append_fields=m.SPLICE_FIELDS,
            spliceai_manifest=str(self.f.workflow.parent/'spliceai-reference-context-v1/manifest.json'),spliceai_manifest_status='PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT')

    def test_all_views_from_vcf_packets_and_person_only_payloads(self):
        result=m.export(self.rule);output=Path(self.rule['output'])
        self.assertEqual(result['report_input_kind'],'FULLY_ANNOTATED_VCF_PACKETS_V1');self.assertTrue(result['sample_packet_scope_verified'])
        self.assertFalse(result['joined_control']['used_as_report_input'])
        self.assertEqual(pin(output/'vcf/cohort.fully-annotated.vcf.gz')['sha256'],pin(self.path.parent/'cohort.fully-annotated.vcf.gz')['sha256'])
        original=self.f.read_table(Path(json.loads(self.joined.read_text())['enriched_master']['path']))
        expected=[dict(row,**dict(zip(m.SPLICE_FIELDS,splice_values(row)))) for row in original]
        self.assertEqual(self.f.read_table(output/'masters/cohort.master.tsv.gz'),expected)
        for sample,label in result['sample_file_labels'].items():
            path=output/'vcf'/(label+'.annotated.vcf.gz');schema=None;ordinals=[]
            with gzip.open(path,'rt') as stream:
                for line in stream:
                    if line.startswith('##WGS_EVIDENCE_SCHEMA='):schema=m.header_schema(line,codec)
                    elif not line.startswith('#'):
                        values=line.rstrip('\r\n').split('\t');info=m.info_dict(values[7]);packet=codec.decode_packet(info['WGS_EVIDENCE'])
                        self.assertEqual(packet['scope'],'sample');self.assertEqual(packet['sample'],sample)
                        self.assertTrue(all(row['values'][schema['full_fields'].index('sample')]==sample for row in packet['rows']))
                        for other in self.f.samples:
                            if other!=sample:self.assertNotIn(other,info['WGS_TX_SCORES'])
                        ordinals.extend(row['i'] for row in packet['rows'])
                        gt=values[9].split(':')[0];self.assertIn('1',gt)
            self.assertEqual(schema['scope'],'sample');self.assertEqual(schema['sample'],sample)
            self.assertEqual(len(ordinals),schema['row_count']);self.assertEqual(ordinals,sorted(set(ordinals)))
        top=self.f.read_table(output/'tables/cohort.top20.tsv')
        self.assertEqual(len({r['allele_id'] for r in top}),len(top));self.assertTrue(all(r['ACMG_class']=='NOT_ASSIGNED' for r in top))
        self.assertTrue(all(r['spliceai_DS_AG']=='0.00' and r['spliceai_DP_AL']=='-12' for r in top))
        self.assertTrue((output/'applied-evidence/population/cohort/nuclear-population-long.tsv.gz').is_file())

    def test_changed_packet_cannot_pass_even_with_updated_file_hashes(self):
        value=json.loads(self.path.read_text())
        fields=value['schema']['full_fields']
        value=self.rewrite(lambda packet:packet['rows'][0]['values'].__setitem__(fields.index('GT'),'0/0'))
        with self.assertRaisesRegex(ValueError,'complete original/control row'):
            m.extract_master(value,codec,self.f.root/'from-vcf.tsv.gz',Path(value['master']['path']),baseline.m.Guard())

    def rewrite(self,change):
        value=json.loads(self.path.read_text());path=Path(value['enriched_vcf']['path']);temp=path.with_name('changed.vcf.gz')
        with gzip.open(path,'rt') as inp,baseline.pysam.BGZFile(str(temp),'wb') as out:
            changed=False
            for line in inp:
                if line.startswith('#'):out.write(line.encode());continue
                parts=line.rstrip('\r\n').split('\t');info=m.info_dict(parts[7]);packet=codec.decode_packet(info['WGS_EVIDENCE'])
                if not changed:change(packet);changed=True
                info['WGS_EVIDENCE']=codec.encode_packet(packet);parts[7]=m.info_text(info);out.write(('\t'.join(parts)+'\n').encode())
        temp.replace(path);value['enriched_vcf']=pin(path);return value

    def test_global_ordinals_reconstructed_not_assumed_vcf_packet_order(self):
        value=self.rewrite(lambda packet:packet['rows'].reverse());target=self.f.root/'reordered.from-vcf.tsv.gz'
        result=m.extract_master(value,codec,target,Path(value['master']['path']),baseline.m.Guard())
        self.assertEqual(result['rows'],value['row_count']);self.assertEqual(self.f.read_table(target),self.f.read_table(Path(value['master']['path'])))

    def test_no_annotation_source_record_fails_closed_without_silent_drop(self):
        value=self.rewrite(lambda packet:packet.__setitem__('rows',[]))
        with self.assertRaisesRegex(ValueError,'compatibility precondition'):
            m.extract_master(value,codec,self.f.root/'empty.from-vcf.tsv.gz',Path(value['master']['path']),baseline.m.Guard())


if __name__=='__main__':unittest.main()
