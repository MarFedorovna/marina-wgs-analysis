"""Frozen fast-reader compatibility with whole-allele mixed-cache regrouping."""
import copy
import gzip
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).parent))
import test_packet_views_fast as public
import test_packet_views_fast_v2 as fixtures
sys.path.insert(0,str(public.ROOT/'scripts'))
import export_packet_views_fast_v2 as frozen


class ResumedBatchContract(unittest.TestCase):
    def setUp(self):
        self.f=fixtures.many_alleles();self.core_fields=public.rank.FIELDS
        self.fields=self.core_fields+['allref_context_%02d'%i for i in range(33)]
        self.lines=[];self.packets=[]
        for line in self.f['lines']:
            columns=line.rstrip('\n').split('\t');info=frozen.info_dict(columns[7]);packet=public.allele_codec.decode_packet(info['WGS_EVIDENCE'])
            for row in packet['rows']:row['values']+=['context\tquoted\nunchanged']*33
            self.packets.append(packet);info['WGS_EVIDENCE']=public.allele_codec.encode_packet(packet);columns[7]=frozen.info_text(info);self.lines.append('\t'.join(columns)+'\n')

    def proofs(self,widths,fields):
        output=[];record=0;ordinal=1
        for number,width in enumerate(widths):
            digest=public.proof.VectorDigest();count=0
            for packet in self.packets[record:record+width]:
                for row in packet['rows']:digest.row(row['values'][:len(fields)]);count+=1
            output.append(dict(id=number,first_i=ordinal,last_i=ordinal+count-1,row_count=count,record_count=width,canonical_row_sha256=digest.digest.hexdigest()))
            record+=width;ordinal+=count
        self.assertEqual(record,len(self.packets));return public.schema(fields,output)

    def test_new_ordered_ids_and_mixed_whole_allele_boundaries_keep_all_vectors(self):
        old=self.proofs([10,10,10],self.fields)
        widths=[10]+[1]*20  # one unchanged cached group, two missing groups split.
        source=self.proofs(widths,self.fields);core=self.proofs(widths,self.core_fields)
        source.update(spliceai_storage=public.allele_codec.CONTRACT,projection_fields=frozen.SPLICE_FIELDS)
        self.assertEqual(frozen.validate_schema(source,core,public.proof),self.fields+frozen.SPLICE_FIELDS)
        self.assertNotEqual(old['ordered_batch_merkle_root'],source['ordered_batch_merkle_root'])
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'public.vcf.gz'
            with gzip.open(path,'wt',encoding='utf-8',newline='') as stream:stream.write('##synthetic-public-proof\n');stream.writelines(self.lines)
            jobs=list(frozen.source_jobs(path,source,core,[],Path(temp)/'unused-workers',hashlib.sha256(path.read_bytes()).hexdigest()))
        observed=[]
        for job in jobs:
            digest=public.proof.VectorDigest()
            for line in job['lines']:
                packet=public.allele_codec.decode_packet(frozen.info_dict(line.rstrip('\n').split('\t')[7])['WGS_EVIDENCE'])
                for row in packet['rows']:observed.append(row);digest.row(row['values'])
            self.assertEqual(digest.digest.hexdigest(),job['proof']['canonical_row_sha256'])
        self.assertEqual(observed,[row for packet in self.packets for row in packet['rows']])
        self.assertEqual([row['i'] for row in observed],list(range(1,241)))

    def test_rehashed_descriptor_cannot_hide_global_ordinal_gap(self):
        source=self.proofs([10]+[1]*20,self.fields);core=self.proofs([10]+[1]*20,self.core_fields)
        source['batches'][1]['first_i']+=1;core['batches'][1]['first_i']+=1
        for value in (source,core):value['ordered_batch_merkle_root']=public.proof.merkle_digest(value['full_fields'],value['batches'])
        source.update(spliceai_storage=public.allele_codec.CONTRACT,projection_fields=frozen.SPLICE_FIELDS)
        with self.assertRaisesRegex(ValueError,'Noncontiguous'):frozen.validate_schema(source,core,public.proof)


if __name__=='__main__':unittest.main()
