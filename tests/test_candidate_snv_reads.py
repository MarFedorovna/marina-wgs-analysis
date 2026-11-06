import importlib.util,json,tempfile,unittest,subprocess,sys,csv
from pathlib import Path
from types import SimpleNamespace
P=Path(__file__).resolve().parents[1]/'scripts/review_candidate_snv_reads.py'
spec=importlib.util.spec_from_file_location('snv_reads',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def read(seq='A',cigar=None,start=99,flag=0,mq=60,quals='default',name='template',rg='rg'):
 class R(SimpleNamespace):
  def get_tag(self,key):
   if self.rg is None:raise KeyError(key)
   return self.rg
 return R(query_sequence=seq,cigartuples=cigar if cigar is not None else [(0,len(seq))],reference_start=start,flag=flag,mapping_quality=mq,query_qualities=[30]*len(seq) if quals=='default' else quals,query_name=name,rg=rg)
T={'chrom':'chr1','pos':100,'ref':'A','alt':'C'}
class Tests(unittest.TestCase):
 def test_cigar_coordinates(self):
  r=read('TTACG',[(4,2),(0,1),(1,1),(7,1)],start=99)
  self.assertEqual(m.target_base(r,99),('A',30,None));self.assertEqual(m.target_base(r,100),('G',30,None))
  for op,reason in [(2,'DELETION_AT_LOCUS'),(3,'REF_SKIP_AT_LOCUS')]:
   self.assertEqual(m.target_base(read('AC',[(0,1),(op,3),(8,1)]),100)[2],reason)
   self.assertEqual(m.target_base(read('AC',[(0,1),(op,3),(8,1)]),103),('C',30,None))
  self.assertEqual(m.target_base(read('C',[(5,5),(6,2),(8,1),(5,3)]),99),('C',30,None))
  self.assertEqual(m.target_base(read('C',[(0,1),(9,1)]),99)[2],'UNSUPPORTED_CIGAR')
 def test_reverse_not_reverse_complemented(self):
  a=m.summarize_reads([read('C',flag=16)],T)
  self.assertEqual(a['counts']['reads_ALT_reverse'],1);self.assertEqual(a['counts']['templates_ALT'],1)
 def test_filters_and_missing_qualities(self):
  rr=[read('C',flag=bit) for bit in m.EXCLUDED_FLAGS]+[read('C',mq=255),read('C',mq=19),read('C',quals=None),read('C',quals=[255]),read('C',quals=[19])]
  x=m.summarize_reads(rr,T)
  self.assertEqual(x['counts'].get('passing_reads',0),0)
  self.assertEqual(x['exclusion_reasons']['MISSING_BASE_QUALITY'],2)
  self.assertEqual(x['exclusion_reasons']['MAPQ_UNKNOWN_255'],1)
  self.assertIn('NO_PASSING_BASE_COVERAGE_NOT_REFERENCE_GENOTYPE',x['descriptive_flags'])
 def test_template_consensus_and_readgroups(self):
  x=m.summarize_reads([read('C'),read('C',flag=16)],T)['counts']
  self.assertEqual(x['reads_ALT'],2);self.assertEqual(x['templates_ALT'],1);self.assertEqual(x['templates_ALT_both'],1)
  x=m.summarize_reads([read('C'),read('A')],T)['counts'];self.assertEqual(x['templates_ALT'],0);self.assertEqual(x['ambiguous_disagreeing_templates'],1)
  self.assertEqual(m.summarize_reads([read('C',rg='one'),read('C',rg='two')],T)['counts']['templates_ALT'],2)
  self.assertEqual(m.summarize_reads([read('C',name='*'),read('C',name=None)],T)['counts']['templates_ALT'],2)
  missing=m.summarize_reads([read('C',rg=None),read('C',rg=None)],T)['counts'];self.assertEqual(missing['templates_ALT'],2);self.assertEqual(missing['passing_reads_missing_RG_uncollapsed'],2)
 def test_no_coverage_no_genotype_and_cap(self):
  x=m.summarize_reads([],T);self.assertEqual(x['counts']['templates_ALT'],0);self.assertNotIn('GT',x)
  with self.assertRaisesRegex(ValueError,'cap exceeded'):m.summarize_reads([read(),read()],T,read_cap=1)
 def test_header_identity(self):
  self.assertEqual(m.header_identity(SimpleNamespace(header=SimpleNamespace(SAMheader=None)),'S')['status'],'HEADER_SM_UNAVAILABLE_EXPLICIT_MAPPING_ONLY')
  self.assertEqual(m.header_identity(SimpleNamespace(header=SimpleNamespace(SAMheader={'RG':[{'SM':'S','ID':'x'}]})),'S')['status'],'EXACT_HEADER_SM_MATCH')
  with self.assertRaisesRegex(ValueError,'SM differs'):m.header_identity(SimpleNamespace(header=SimpleNamespace(SAMheader={'RG':[{'SM':'OTHER'}]})),'S')
 def test_real_four_worker_cli_public_bam(self):
  import bamnostic
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);samples={('PUBLIC%d'%i):{'bam':bamnostic.example_bam,'index':bamnostic.example_bam+'.bai'} for i in range(4)}
   (root/'bams.json').write_text(json.dumps({'schema_version':1,'samples':samples}))
   with (root/'candidates.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['chrom','pos','ref','alt','hcm_phenotype_priority','sample_calls_json'],delimiter='\t');w.writeheader();w.writerow(dict(chrom='chr1',pos=100,ref='A',alt='C',hcm_phenotype_priority='1',sample_calls_json=json.dumps({s:{'GT':'./.','FT':'','DP':'','AD':'','GQ':''} for s in samples})))
   command=[sys.executable,str(P),'--candidates',str(root/'candidates.tsv'),'--bam-manifest',str(root/'bams.json'),'--output',str(root/'out'),'--workers','4']
   dry=subprocess.run(command,check=True,capture_output=True,text=True);self.assertIn('PLAN_ONLY_NO_BAM_OPEN',dry.stdout);self.assertFalse((root/'out').exists())
   subprocess.run(command+['--execute'],check=True,capture_output=True,text=True,timeout=30)
   result=json.loads((root/'out/summary.json').read_text());self.assertEqual(result['rows'],4);self.assertTrue(all(r['source_call']['GT']=='./.' and r['counts']['fetched_reads']==2 for r in result['results']))
   repeat=subprocess.run(command+['--execute'],capture_output=True,text=True);self.assertNotEqual(repeat.returncode,0)
 def test_actual_public_bam_worker_and_unsupported(self):
  import bamnostic
  entry={'bam':bamnostic.example_bam,'index':bamnostic.example_bam+'.bai'}
  t=dict(T,variant_id='public',gene='',source_sample_calls={'PUBLIC':{'GT':'./.','DP':'','AD':'','GQ':'','FT':''}},source_membership={})
  out=m.worker('PUBLIC',entry,[t,dict(t,alt='CT')],20,20,100000,30)
  self.assertEqual(len(out['rows']),2);self.assertEqual(out['rows'][0]['counts']['fetched_reads'],2)
  self.assertEqual(out['rows'][0]['counts']['reads_REF'],1);self.assertEqual(out['rows'][0]['counts']['reads_ALT'],0);self.assertEqual(out['rows'][0]['counts']['templates_REF'],1)
  self.assertEqual(out['rows'][0]['source_call']['GT'],'./.');self.assertEqual(out['rows'][1]['status'],'NOT_ASSESSED_UNSUPPORTED_NON_SNV')
  self.assertNotIn('query_name',json.dumps(out));self.assertNotIn('query_sequence',json.dumps(out))
if __name__=='__main__':unittest.main()
