"""Public/synthetic pileup and denominator checks, without patient inputs."""
import collections
import importlib.util
import json
import multiprocessing
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import datetime as dt
import unittest

spec=importlib.util.spec_from_file_location('xy_target_qc',Path(__file__).resolve().parents[1]/'scripts/xy_target_quantitative_qc_v2.py')
q=importlib.util.module_from_spec(spec);sys.modules[spec.name]=q;spec.loader.exec_module(q)


class XYTargetTests(unittest.TestCase):
    def test_public_targets_fail_closed_on_wrong_gc_par_overlap_and_oversize(self):
        class Ref:
            rows={'chr1':(10000000,), 'chrX':(156040895,), 'chrY':(57227415,)}
            def get(self,c,a,b):return ('ACGT'*126)[:b-a]
        header='chrom\tstart0\tend0\tgroup\tgc_fraction\n'
        good=['chr1\t3000000\t3000100\tAUTOSOME_MATCHED\t0.5','chrX\t3000000\t3000100\tX_NONPAR_EXON\t0.5','chrY\t3000000\t3000100\tY_NONPAR_EXON\t0.5']
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'targets.tsv';p.write_text(header+'\n'.join(good)+'\n')
            self.assertEqual(len(q.target_windows(p,Ref())),3)
            for bad in [good+['chrX\t3000050\t3000150\tX_NONPAR_EXON\t0.5'],
                        [good[0],good[1].replace('0.5','0.4'),good[2]],
                        [good[0],good[1].replace('3000000\t3000100','10001\t10101'),good[2]],
                        [good[0],good[1].replace('3000100','3000501'),good[2]]]:
                p.write_text(header+'\n'.join(bad)+'\n')
                with self.assertRaises(ValueError):q.target_windows(p,Ref())
    def test_pileup_ignores_markers_and_indel_sequence_but_counts_observations(self):
        value=q.parse_bases('A','^F.,+2tgG-1a*$','IIII',4)
        self.assertEqual(value,collections.Counter(A=2,G=1,**{'*':1}))
    def test_unknown_truncated_low_quality_and_inconsistent_pileup_fail(self):
        for bases,qual,dp in [('^','I',1),('.+4AG','I',1),('z','I',1),('..','I',2),('.','!',1)]:
            with self.subTest(bases=bases),self.assertRaises(ValueError):q.parse_bases('A',bases,qual,dp)
    def test_zero_depth_and_quality_filtered_het_hom_and_ambiguous_states(self):
        self.assertEqual(q.parse_bases('A','*','*',0),{})
        for counts,state in [(dict(A=6,G=6),'HET_BALANCED'),(dict(A=10),'HOM_REF_LIKE'),(dict(G=10),'HOM_ALT_LIKE'),(dict(A=8,G=2),'AMBIGUOUS_AB'),(dict(A=5,G=4),'LOW_DP'),(dict(A=8,G=8,T=2),'OFF_ALLELE_OR_NONBASE')]:
            self.assertEqual(q.classify(collections.Counter(counts),'A','G')['status'],state)
    def test_no_f_without_expected_he_or_sufficient_sites_and_spread(self):
        rows=[dict(pos=i*10000000+1,af=.5,status='HET_BALANCED') for i in range(12)]
        self.assertIsNone(q.heterozygosity(rows)['F_descriptive'])
        rows=rows*10
        self.assertIsNone(q.heterozygosity(rows,False)['F_descriptive'])
        result=q.heterozygosity(rows,True)
        self.assertEqual(result['F_descriptive'],-1)
        self.assertFalse(result['sex_assigned'])
    def test_par_boundaries_and_bed_union(self):
        self.assertFalse(q.in_par('chrX',10000));self.assertTrue(q.in_par('chrX',10001))
        self.assertTrue(q.in_par('chrY',57217415));self.assertFalse(q.in_par('chrY',57217416))
        windows=[dict(chrom='chrX',start=10,end=20),dict(chrom='chrX',start=18,end=30)]
        self.assertEqual(q.merge_regions(windows,[dict(chrom='chrX',pos=31)],['chrX']),[('chrX',10,31)])
    def test_wilson_zero_and_all_events_are_bounded(self):
        self.assertIsNone(q.wilson(0,0))
        self.assertEqual(q.wilson(0,10)[0],0)
        self.assertAlmostEqual(q.wilson(10,10)[1],1)

    @unittest.skipUnless(sys.platform=='linux' and shutil.which('samtools'),'Actual indexed public-BAM proof requires Linux samtools')
    def test_actual_four_public_bams_indexed_selection_quality_flags_and_zero_denominator(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);out=root/'output';out.mkdir();ref=root/'public.fa';sequence='ACGT'*125
            ref.write_text(''.join('>'+c+'\n'+sequence+'\n' for c in ('chr1','chrX','chrY')))
            subprocess.run(['samtools','faidx',str(ref)],check=True,capture_output=True)
            bed=root/'targets.bed';bed.write_text('chr1\t79\t149\nchrX\t79\t149\nchrY\t79\t149\n')
            samples=[]
            for sample in range(4):
                sam=root/f'PUBLIC_{sample+1}.sam';lines=['@HD\tVN:1.6\tSO:coordinate']+['@SQ\tSN:'+c+'\tLN:500' for c in ('chr1','chrX','chrY')]
                def read(name,chrom,flag=0,mq=60,alt=False,low=False):
                    bases=list(sequence[79:149]);quality=list('I'*70)
                    if alt:bases[20]='C'
                    if low:quality[20]='!'
                    lines.append('\t'.join(map(str,[name,flag,chrom,80,mq,'70M','*',0,0,''.join(bases),''.join(quality)])))
                if sample<3:
                    for i in range(12):read('auto'+str(i),'chr1')
                    for i in range(6 if sample==2 else 12):read('x'+str(i),'chrX',alt=sample==0 and i>=6)
                    for i,flag in enumerate((256,512,1024,2048)):read('excluded'+str(i),'chrX',flag=flag,alt=True)
                    read('lowmq','chrX',mq=20,alt=True);read('lowbq','chrX',alt=True,low=True)
                else:
                    for i in range(10):read('y'+str(i),'chrY')
                sam.write_text('\n'.join(lines)+'\n');bam=root/f'PUBLIC_{sample+1}.bam'
                subprocess.run(['samtools','sort','-o',str(bam),str(sam)],check=True,capture_output=True)
                subprocess.run(['samtools','index',str(bam)],check=True,capture_output=True)
                receipt=root/f'PUBLIC_{sample+1}.json';receipt.write_text('{}')
                samples.append(dict(sample_id='PUBLIC_'+str(sample+1),bam=q.stat_identity(bam),index=q.pin(str(bam)+'.bai'),receipt=q.pin(receipt)))
            target=root/'public-coding-exon-targets.tsv'
            target.write_text('chrom\tstart0\tend0\tgroup\tgc_fraction\n'+''.join(c+'\t79\t149\t'+g+'\t'+str((sequence[79:149].count('G')+sequence[79:149].count('C'))/70)+'\n' for c,g in [('chr1','AUTOSOME_MATCHED'),('chrX','X_NONPAR_EXON'),('chrY','Y_NONPAR_EXON')]))
            public_reference=q.Fasta(ref);public_targets=q.target_windows(target,public_reference);public_reference.close()
            cpu=sorted(os.sched_getaffinity(0));self.assertGreaterEqual(len(cpu),16)
            p=dict(samples=samples,output=str(out),deadline_utc=(dt.datetime.now(dt.timezone.utc)+dt.timedelta(seconds=90)).isoformat(),max_wall_seconds=90,cpu_ids=cpu[-16:],
                reference=q.stat_identity(ref),bed=q.pin(bed),expected_he_available=True,sites=[dict(chrom='chrX',pos=100,ref='T',alt='C',af=.5,source_id='PUBLIC_SNP')],
                windows=public_targets)
            with multiprocessing.get_context('fork').Pool(4) as pool:results=pool.map(q.measure,[(p,i) for i in range(4)])
            self.assertEqual([r['X_nonPAR']['n_called'] for r in results],[1,1,0,0])
            self.assertEqual(results[0]['X_nonPAR']['n_het'],1)
            with (out/'sample_0001.x-site-counts.tsv').open() as stream:rows=list(__import__('csv').DictReader(stream,delimiter='\t'))
            self.assertEqual((int(rows[0]['ref_count']),int(rows[0]['alt_count'])),(6,6))
            self.assertEqual(results[0]['survey']['chrY']['breadth_ge1'],0)
            self.assertGreater(results[3]['survey']['chrY']['breadth_ge1'],0)
            self.assertIsNone(results[3]['sex_chromosome_to_autosome_median_window_depth']['chrY'])
            self.assertTrue(all(r['X_nonPAR']['sufficiency_status']=='INSUFFICIENT_CALLABLE_SITES_OR_SPREAD' for r in results))
            self.assertTrue(all(r['X_nonPAR']['F_descriptive'] is None and not r['sex_assigned'] for r in results))


if __name__=='__main__':unittest.main()
