import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('xy_new',ROOT/'scripts/run_xy_germline.py');xy=importlib.util.module_from_spec(spec);spec.loader.exec_module(xy)

class XY(unittest.TestCase):
    def test_partition_union_and_boundaries(self):
        lock=xy.load(ROOT/'config/sex-chromosomes.resources.lock.json');parts=xy.intervals(lock)
        for chrom in ('chrX','chrY'):
            spans=[(a,b) for rows in parts.values() for c,a,b in rows if c==chrom]
            if chrom=='chrY':spans+=[tuple(x) for x in lock['par_one_based_closed'][chrom]]
            spans.sort();self.assertEqual(spans[0][0],1);self.assertEqual(spans[-1][1],lock['primary_contig_lengths'][chrom]);self.assertTrue(all(b+1==c for (a,b),(c,d) in zip(spans,spans[1:])))
    def fixture(self,d,gt='1\t0/1',info='.'):
        p=d/'v.vcf.gz'
        with gzip.open(p,'wt') as f:f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\th\td\nchrX\t10\t.\tA\tC\t50\tPASS\t'+info+'\tGT\t'+gt+'\n')
        return p
    def test_mixed_gt_actual_vectors_and_wrong_ploidy(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.fixture(Path(d));self.assertEqual(xy.vcf_check(p,['h','d'],{'h':1,'d':2},[('chrX',1,20)]),1)
            with self.assertRaises(ValueError):xy.vcf_check(p,['h','d'],{'h':2,'d':2},[('chrX',1,20)])
    def test_reference_block_boundary_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.fixture(Path(d),info='END=21')
            with self.assertRaises(ValueError):xy.vcf_check(p,['h','d'],{'h':1,'d':2},[('chrX',1,20)])
    def test_excluded_y_no_reference_invention(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.fixture(Path(d),'1\t.');self.assertEqual(xy.vcf_check(p,['h','d'],{'h':1},[('chrX',1,20)],['d']),1)
            p=self.fixture(Path(d),'1\t0/0')
            with self.assertRaises(ValueError):xy.vcf_check(p,['h','d'],{'h':1},[('chrX',1,20)],['d'])
    def test_final_merged_gt_and_duplicates(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.fixture(Path(d));samples=[{'sample':'h','X_NONPAR':1},{'sample':'d','X_NONPAR':2}]
            self.assertEqual(xy.check_final(p,['h','d'],samples,{'X_NONPAR':[('chrX',1,20)]}),1)
            with gzip.open(p,'at') as f:f.write('chrX\t10\t.\tA\tC\t50\tPASS\t.\tGT\t1\t0/1\n')
            with self.assertRaises(ValueError):xy.check_final(p,['h','d'],samples,{'X_NONPAR':[('chrX',1,20)]})
    def test_source_identity_no_large_hash(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bam';p.write_bytes(b'x');v={'path':str(p),'bytes':1,'mtime_ns':p.stat().st_mtime_ns,'sha256':'a'*64};self.assertEqual(xy.identity(v),p)
            v['bytes']=2
            with self.assertRaises(ValueError):xy.identity(v)
    def test_soft_filter_preserves_haploid_missing_and_all_gt(self):
        module=xy.module(ROOT/'scripts/calling_stages.py')
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);source=self.fixture(d,'1\t.');allrows=d/'all.vcf';selected=d/'selected.vcf';module.soft_filter_vcf(source,allrows,selected)
            row=[x for x in allrows.read_text().splitlines() if not x.startswith('#')][0].split('\t');self.assertEqual(row[9].split(':')[0],'1');self.assertEqual(row[10].split(':')[0],'.');self.assertIn('MissingAnnotation',row[6]);self.assertIn('NoCall',row[10])

if __name__=='__main__':unittest.main()
