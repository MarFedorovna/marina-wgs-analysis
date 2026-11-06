import copy
import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import export_packet_views_fast as m
import fast_vcf_evidence_candidate as proof
import spliceai_vcf_codec as codec
import spliceai_allele_info_codec as allele_codec
import annotation_stages as rank
import export_cohort_sample_views as policy
try:import pysam
except ImportError:pysam=None


def schema(fields,batches):
    return dict(v=1,codec='json-utf8+zlib+base64url-nopad',scope='cohort',sample=None,original_fields=rank.FIELDS,
        appended_fields=fields[len(rank.FIELDS):],full_fields=fields,row_count=sum(p['row_count'] for p in batches),
        ordinal_in_vcf_order=True,vector_digest_contract=proof.DIGEST_CONTRACT,batches=batches,
        ordered_batch_merkle_root=proof.merkle_digest(fields,batches),producer_kind=m.PRODUCER)


def fixture(folder):
    samples=['SYNTHETIC_A','SYNTHETIC_B','SYNTHETIC_C','SYNTHETIC_D'];fields=rank.FIELDS+m.SPLICE_FIELDS
    rows=[];packets=[];batches=[];prefix=[];lines=[]
    for j,(gts,effect,score) in enumerate([(['0/1','0/1','0/0','./.'],'stop_gained','0.00'),
        (['1/.','0/0','0/1','0/0'],'missense_variant','0.80'),(['0/1']*4,'intron_variant','0.90')]):
        pos=10*(j+1);allele=['chr1',pos,'A','T'];begin=len(rows)+1;vectors=[];base=proof.VectorDigest();full=proof.VectorDigest()
        predictions=[dict(variant_key=rank.variant_key(*allele),chrom='chr1',pos=str(pos),ref='A',alt='T',gene='GENE',
            status='SCORED_EXACT_GRCH38_ALLELE_GENE',DS_AG='0.00',DS_AL=score,DS_DG='0.00',DS_DL='0.00',DP_AG='0',DP_AL='-22',DP_DG='0',DP_DL='0',
            raw_annotation='official-shaped synthetic test',source_json='{"model":"synthetic-only"}',transcript_scope='GENCODE24_GENE_UNION_NOT_ENST49_MAPPING')]
        fts=['PASS','PASS','LowDP' if j==1 else 'PASS','PASS']
        for tx in ['TX.1','TX.2']:
            for sample,gt,ft in zip(samples,gts,fts):
                row={f:'' for f in rank.FIELDS};row.update(sample=sample,variant_id=rank.variant_key(*allele),chrom='chr1',pos=str(pos),ref='A',alt='T',
                    GT=gt,FT=ft,DP='20',AD='10,10',GQ='50',site_filter='PASS',consequence=effect,gene='GENE',transcript=tx,
                    candidate_eligible='1' if policy.call_flags(gt,ft,'PASS')['qc'] and effect!='intron_variant' else '0',
                    consequence_rank='5' if tx=='TX.1' else '10')
                row['local_evidence_json']='{"raw":"synthetic TAB\tand LF\nremain lossless"}'
                basevalues=[row[f] for f in rank.FIELDS]
                extra=['SCORED_EXACT_GRCH38_ALLELE_GENE','GENE','0.00',score,'0.00','0.00','0','-22','0','0',json.dumps(predictions),'{}','GENCODE24_GENE_UNION_NOT_ENST49_MAPPING']
                values=basevalues+extra;base.row(basevalues);full.row(values);rows.append(values);vectors.append(dict(i=len(rows),values=values))
        descriptor=dict(id=j,first_i=begin,last_i=len(rows),row_count=len(vectors),record_count=1)
        batches.append(dict(descriptor,canonical_row_sha256=full.digest.hexdigest()));prefix.append(dict(descriptor,canonical_row_sha256=base.digest.hexdigest()))
        packet=dict(v=1,allele=allele,scope='cohort',sample=None,rows=vectors,evidence={},gene_context_keys=[],membership={'sample_states_json':json.dumps(samples)},unknown_private={'all_GT':gts});packets.append(packet)
        ac=sum(gt.count('1') for gt in gts);an=sum(c in '01' for gt in gts for c in gt)
        tags='AC=%d;AN=%d;AF=%s;WGS_EVIDENCE=%s'%(ac,an,str(ac/an),codec.encode_packet(packet))
        line='chr1\t%d\t.\tA\tT\t90\tPASS\t%s\tGT:FT:DP:AD:GQ\t%s\n'%(pos,tags,'\t'.join(gt+':'+ft+':20:10,10:50' for gt,ft in zip(gts,fts)))
        lines.append(line)
    original=schema(rank.FIELDS,prefix);extended=schema(fields,batches)
    header='##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
    for tag,num,kind in [('AC','A','Integer'),('AN','1','Integer'),('AF','A','Float')]:header+='##INFO=<ID=%s,Number=%s,Type=%s,Description="source">\n'%(tag,num,kind)
    for tag,num,kind in [('GT','1','String'),('FT','1','String'),('DP','1','Integer'),('AD','R','Integer'),('GQ','1','Integer')]:header+='##FORMAT=<ID=%s,Number=%s,Type=%s,Description="source">\n'%(tag,num,kind)
    header+=''.join(line.decode() for line in codec.headers(extended,codec.encode_packet(dict(v=1,scope='PUBLIC_GENE_CONTEXT',genes={})),codec.encode_packet(dict(v=1,scope='PUBLIC_INSTALLED_REFERENCE_OBJECTS',sources=[]))))
    header+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n'
    return dict(samples=samples,fields=fields,rows=rows,packets=packets,lines=lines,header=header,core_schema=original,schema=extended)


def project_fixture(f):
    source_json=json.dumps(dict(version='1.3.1',distance=500,mask=0,annotation_scope=allele_codec.SCOPE,assembly='GRCh38',
        model_sha256=['a'*64]*5,annotation_sha256='b'*64,reference_sha256='c'*64),sort_keys=True,separators=(',',':'))
    expected=[];lines=[];packets=[]
    for line,packet in zip(f['lines'],f['packets']):
        packet=copy.deepcopy(packet);old=packet['rows'][0]['values'];preds=json.loads(old[f['fields'].index('spliceai_prediction_json')])
        for prediction in preds:
            prediction['variant_key']=allele_codec.projection.key(*packet['allele'])
            prediction['source_json']=source_json
            prediction['raw_annotation']='|'.join([prediction['alt'],prediction['gene'],*[prediction[field] for field in allele_codec.SCORE_FIELDS]])
        for row in packet['rows']:
            row['values']=row['values'][:len(rank.FIELDS)];base=dict(zip(rank.FIELDS,row['values']))
            extra=allele_codec.enrich_row(base,preds,source_json);expected.append(row['values']+[extra[field] for field in m.SPLICE_FIELDS])
        parts=line.rstrip('\n').split('\t');info=m.info_dict(parts[7]);info['WGS_EVIDENCE']=allele_codec.encode_packet(packet)
        info.update(allele_codec.overlay_info(packet['allele'],preds,source_json));parts[7]=m.info_text(info)
        lines.append('\t'.join(parts)+'\n');packets.append(packet)
    physical=copy.deepcopy(f['core_schema']);projection=dict(v=1,projection_contract=allele_codec.CONTRACT,projection_fields=m.SPLICE_FIELDS,scope=allele_codec.SCOPE,model=json.loads(source_json))
    headers=[]
    for line in f['header'].splitlines(keepends=True):
        if line.startswith('##WGS_EVIDENCE_SCHEMA=') or line.startswith(('##INFO=<ID=SpliceAI,','##INFO=<ID=WGS_SPLICEAI_')):continue
        if line.startswith('#CHROM'):
            headers.append(allele_codec.schema_header(physical).decode());headers.extend(v.decode() for v in allele_codec.overlay_headers(projection))
        headers.append(line)
    f.update(rows=expected,packets=packets,lines=lines,header=''.join(headers),schema=physical,projection=projection)
    return f


class FastPureTests(unittest.TestCase):
    def test_missing_not_zero_and_unassessed_not_benign(self):
        self.assertIsNone(m.scored_max([dict(status='UNSUPPORTED',DS_AG='')]))
        self.assertEqual(m.scored_max([dict(status='SCORED_EXACT_GRCH38_ALLELE_GENE',DS_AG='0.00')]),0)
        with self.assertRaises(ValueError):m.scored_max([dict(status='SCORED_EXACT_GRCH38_ALLELE_GENE',DS_AG='nan')])

    def test_merkle_and_prefix_mutation_rejected(self):
        f=fixture(Path('.'));m.validate_schema(f['schema'],f['core_schema'],proof)
        bad=copy.deepcopy(f['schema']);bad['batches'][0]['canonical_row_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'commitment'):m.validate_schema(bad,f['core_schema'],proof)
        bad=copy.deepcopy(f['schema']);bad['full_fields'][0]='WRONG'
        with self.assertRaisesRegex(ValueError,'prefix'):m.validate_schema(bad,f['core_schema'],proof)

    def test_expired_never_creates_output(self):
        with self.assertRaisesRegex(ValueError,'deadline'):m.export(dict(deadline_utc='2000-01-01T00:00:00+00:00'))

    def test_canonical_vector_quote_roundtrip(self):
        fields=['a','b'];values=['embedded\tTAB','embedded\nLF']
        self.assertEqual(next(csv.reader(io.StringIO(proof.canonical_row_bytes(values).decode()),delimiter='\t')),values)


@unittest.skipUnless(pysam and sys.platform=='linux','Linux pysam BGZF integration')
class FastIntegration(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name);self.f=fixture(self.root)
        if getattr(self,'projected',False):self.f=project_fixture(self.f)
        folder=self.root/'fully-annotated-vcf-v1';folder.mkdir();source=folder/'cohort.fully-annotated.vcf.gz'
        with pysam.BGZFile(str(source),'wb') as stream:stream.write((self.f['header']+''.join(self.f['lines'])).encode())
        pysam.tabix_index(str(source),preset='vcf')
        core=folder/'core-manifest.json';core.write_text(json.dumps(dict(status='DIRECT_PACKET_BASELINE_PASS',schema=self.f['core_schema'],no_legacy_receipt_created=True)))
        manifest=folder/'manifest.json';value=dict(status='PASS_FULLY_ANNOTATED_VCF_DIRECT',schema_version=1,
            schema=self.f['schema'],samples=self.f['samples'],enriched_vcf=m.pin(source),index=m.pin(str(source)+'.tbi'),core_manifest=m.pin(core),
            direct_core_prefix_verified=True,packet_decode_verified=True,original_records_fields_preserved=True)
        if self.f.get('projection'):value.update(projection=self.f['projection'],core33_packets_byte_identical=True)
        manifest.write_text(json.dumps(value))
        self.rule=dict(root_reviewed=True,synthetic=True,workers=2,deadline_utc=m.END.isoformat(),output=str(self.root/'cohort-sample-views-v1'),fully_manifest=m.pin(manifest),
            exporter_sha256=m.sha(m.__file__),approved_appended_fields=m.SPLICE_FIELDS,**{name:m.pin(ROOT/'scripts'/file) for name,file in [
                ('proof_module','fast_vcf_evidence_candidate.py'),('codec_module','spliceai_vcf_codec.py'),('view_policy','export_cohort_sample_views.py'),('ranking_module','annotation_stages.py')]})
        if getattr(self,'projected',False):self.rule['codec_module']=m.pin(ROOT/'scripts/spliceai_allele_info_codec.py')

    def test_one_pass_lossless_masters_privacy_rank_and_splice_appendix(self):
        result=m.export(self.rule);out=Path(self.rule['output'])
        with gzip.open(out/'masters/cohort.master.tsv.gz','rt',newline='') as stream:
            rows=list(csv.reader(stream,delimiter='\t'));self.assertEqual(rows[0],self.f['fields']);self.assertEqual(rows[1:],self.f['rows'])
        with (out/'tables/cohort.top20.tsv').open() as stream:
            top=list(csv.DictReader(stream,delimiter='\t'));self.assertEqual(len(top),1);self.assertEqual(top[0]['pos'],'10');self.assertEqual(top[0]['pct_of4'],'50')
        with (out/'tables/spliceai/cohort.top20.tsv').open() as stream:
            top=list(csv.DictReader(stream,delimiter='\t'));self.assertEqual([r['pos'] for r in top],['30','10']);self.assertEqual(top[0]['candidate_eligible'],'0')
        for i,sample in enumerate(self.f['samples'],1):
            with pysam.VariantFile(str(out/'vcf'/('sample_'+str(i).zfill(4)+'.annotated.vcf.gz'))) as stream:
                self.assertEqual(list(stream.header.samples),[sample])
                for record in stream.fetch():
                    packet=codec.decode_packet(record.info['WGS_EVIDENCE']);self.assertEqual(packet['sample'],sample)
                    self.assertNotIn('membership',packet);self.assertNotIn('unknown_private',packet)
                    self.assertTrue(all(r['values'][self.f['fields'].index('sample')]==sample for r in packet['rows']))
                    self.assertEqual(record.info['AN'],sum(a is not None for a in record.samples[sample]['GT']))
        self.assertTrue(result['no_legacy_receipt_created']);self.assertFalse((out/'stages').exists())

    def test_packet_ordinal_mutation_rejected_without_promoted_output(self):
        # A valid enclosing compressed file SHA cannot conceal a changed vector proof.
        manifest_path=Path(self.rule['fully_manifest']['path']);manifest=m.load(manifest_path);source=Path(manifest['enriched_vcf']['path'])
        packet=copy.deepcopy(self.f['packets'][0]);packet['rows'][0]['i']=2
        line=self.f['lines'][0].replace(codec.encode_packet(self.f['packets'][0]),codec.encode_packet(packet))
        with pysam.BGZFile(str(source),'wb') as stream:stream.write((self.f['header']+line+''.join(self.f['lines'][1:])).encode())
        manifest['enriched_vcf']=m.pin(source);manifest_path.write_text(json.dumps(manifest));self.rule['fully_manifest']=m.pin(manifest_path)
        with self.assertRaisesRegex(ValueError,'ordinal'):m.export(self.rule)
        self.assertFalse(Path(self.rule['output']).exists())

    def test_source_identity_rejected_before_report_work(self):
        self.rule['authority_expected']={'core.source_pins':{'wrong':'source'}}
        with self.assertRaisesRegex(ValueError,'authority field missing'):m.export(self.rule)
        self.assertFalse(Path(self.rule['output']).exists())


@unittest.skipUnless(pysam and sys.platform=='linux','Linux pysam projection integration')
class AlleleInfoProjectionIntegration(FastIntegration):
    projected=True

    def test_projection_has_separate_derived_proof_and_no_physical_field_claim(self):
        result=m.export(self.rule)
        self.assertEqual(len(self.f['schema']['full_fields']),len(rank.FIELDS))
        self.assertEqual(result['logical_full_fields'],rank.FIELDS+m.SPLICE_FIELDS)
        self.assertNotEqual(result['derived_batch_merkle_root'],result['batch_merkle_root'])


if __name__=='__main__':unittest.main()
