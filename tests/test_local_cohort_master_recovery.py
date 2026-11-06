from pathlib import Path
import csv,gzip,hashlib,io,json,os,sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'transfer'))
import rebuild_local_cohort_master as r

@unittest.skipUnless(os.name=='nt','Actual Windows multiprocessing/checkpoint regression')
class LocalRecovery(unittest.TestCase):
    def fixture(self,folder):
        folder=folder/'source';folder.mkdir()
        c=r.preview.module('spliceai_allele_info_codec');proof=r.preview.module('fast_vcf_evidence_candidate')
        fields=['sample','chrom','pos','ref','alt','gene','GT','FT','DP','AD','GQ','site_filter']+[f'public_extra_{i}' for i in range(126)];self.assertEqual(len(fields),138);samples=['PUBLIC_1','PUBLIC_2','PUBLIC_3','PUBLIC_4'];logical=fields+c.SPLICE_FIELDS
        source=json.dumps(dict(version='1.3.1',distance=500,mask=0,annotation_scope=c.SCOPE,assembly='GRCh38',model_sha256=[str(i)*64 for i in range(1,6)],annotation_sha256='a'*64,reference_sha256='b'*64),sort_keys=True)
        prediction=dict(variant_key='chr11:47332146:T:G',chrom='chr11',pos='47332146',ref='T',alt='G',gene='MYBPC3',status='SCORED_EXACT_GRCH38_ALLELE_GENE',DS_AG='0.02',DS_AL='0.60',DS_DG='0.03',DS_DL='0.01',DP_AG='-164',DP_AL='-22',DP_DG='-50',DP_DL='-301',raw_annotation='G|MYBPC3|0.02|0.60|0.03|0.01|-164|-22|-50|-301',source_json=source,transcript_scope=c.SCOPE)
        source_proofs=[];derived=[];records=[];all_expected=[];ordinal=1
        for batchid,(allele,gene,predictions) in enumerate([(['chr1',20,'A','T'],'',[]),(['chr11',47332146,'T','G'],'MYBPC3',[prediction])]):
            rows=[];base_digest=proof.VectorDigest();logical_digest=proof.VectorDigest();first=ordinal
            for sid,gt in zip(samples,['0|1','./1','.','0/0']):
                row=dict.fromkeys(fields,'');row.update(sample=sid,chrom=allele[0],pos=str(allele[1]),ref=allele[2],alt=allele[3],gene=gene,GT=gt,FT='PASS' if gt=='0|1' else 'NoCall',DP='25',AD='12,13',GQ='40',site_filter='PASS',public_extra_0='Unicode α\tTAB\nLINE "quote"')
                values=[row[k] for k in fields];rows.append(dict(i=ordinal,values=values));base_digest.row(values);ordinal+=1
                # Reference implementation deliberately calls pure projection with the full row.
                expected=values+[c.enrich_row(row,predictions,source if predictions else '')[key] for key in c.SPLICE_FIELDS];all_expected.append(expected);logical_digest.row(expected)
            b=dict(id=batchid,first_i=first,last_i=ordinal-1,row_count=4,record_count=1,canonical_row_sha256=base_digest.digest.hexdigest());source_proofs.append(b);derived.append(dict(b,canonical_row_sha256=logical_digest.digest.hexdigest()))
            packet=dict(v=1,allele=allele,rows=rows,scope='cohort',sample=None);info=dict(WGS_EVIDENCE=c.encode_packet(packet),**c.overlay_info(allele,predictions,source if predictions else ''));text=';'.join(k+'='+v for k,v in info.items());records.append('\t'.join([allele[0],str(allele[1]),'.',allele[2],allele[3],'50','PASS',text,'GT:FT','0|1:PASS','./1:NoCall','.:NoCall','0/0:PASS'])+'\n')
        schema=dict(v=1,codec=c.CODEC,scope='cohort',sample=None,full_fields=fields,row_count=8,ordinal_in_vcf_order=True,batches=source_proofs,ordered_batch_merkle_root=proof.merkle_digest(fields,source_proofs))
        projection=dict(v=1,projection_contract=c.CONTRACT,projection_fields=c.SPLICE_FIELDS,scope=c.SCOPE,model={})
        vcf=folder/'public.vcf.gz';header=c.schema_header(schema).decode()+c.overlay_headers(projection)[-1].decode()+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n'
        with gzip.open(vcf,'wt',encoding='utf-8',newline='') as f:f.write(header+''.join(records))
        expected_source=dict(bytes=vcf.stat().st_size,sha256='d'*64);note=dict(Status='COMPLETED_SIZE_ONLY',SourceIdentity=dict(ExpectedBytes=vcf.stat().st_size,ExpectedServerSha256='d'*64),LocalMtimeNs=vcf.stat().st_mtime_ns);vcf.with_name(vcf.name+'.completed.json').write_text(json.dumps(note))
        view=dict(status='PASS',synthetic=True,report_input_kind='FULLY_ANNOTATED_VCF_PACKETS_V1',source_pins=dict(vcf=expected_source),logical_full_fields=logical,batch_merkle_root=schema['ordered_batch_merkle_root'],derived_batch_merkle_root=proof.merkle_digest(logical,derived),master_transcript_rows=8,vcf_records=2);mp=folder/'manifest.json';mp.write_text(json.dumps(view))
        expected=io.StringIO(newline='');w=csv.writer(expected,delimiter='\t',lineterminator='\n');w.writerow(logical);w.writerows(all_expected)
        return vcf,mp,expected.getvalue().encode()
    def test_real_two_worker_windows_full151field_semantic_equality_resume(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);vcf,manifest,expected=self.fixture(folder);out=folder/'reconstructed'
            result=r.recover(vcf,manifest,out,workers=2,max_hours=.1,synthetic=True);self.assertEqual(result['fields'],151);self.assertEqual(result['rows'],8);self.assertTrue(result['canonical_TSV_UTF8_LF_vectors_equal']);self.assertFalse(result['source_file_sha_recomputed'])
            final=out/'cohort.master.reconstructed.tsv.gz'
            with gzip.open(final,'rb') as f:self.assertEqual(f.read(),expected)
            stamps={p.name:p.stat().st_mtime_ns for p in (out/'verified-batches').glob('*.tsv.gz')};final.unlink();(out/'LOCAL-RECONSTRUCTION-RECEIPT.json').unlink()
            again=r.recover(vcf,manifest,out,workers=2,max_hours=.1,synthetic=True);self.assertEqual(again['derived_merkle_root'],result['derived_merkle_root']);self.assertEqual(stamps,{p.name:p.stat().st_mtime_ns for p in (out/'verified-batches').glob('*.tsv.gz')})
    def test_wrong_source_and_wrong_derived_merkle_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);vcf,mp,_=self.fixture(folder);v=json.loads(mp.read_text());v['source_pins']['vcf']['sha256']='f'*64;mp.write_text(json.dumps(v))
            with self.assertRaisesRegex(ValueError,'source identity'):r.recover(vcf,mp,folder/'wrong',2,.1,True)
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);vcf,mp,_=self.fixture(folder);v=json.loads(mp.read_text());v['derived_batch_merkle_root']='f'*64;mp.write_text(json.dumps(v))
            with self.assertRaisesRegex(ValueError,'vectors differ'):r.recover(vcf,mp,folder/'wrong',2,.1,True)
            self.assertFalse((folder/'wrong/cohort.master.reconstructed.tsv.gz').exists())

if __name__=='__main__':unittest.main()
