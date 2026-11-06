from pathlib import Path
import copy,csv,datetime as dt,gzip,json,sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'));import rerank_phenotype_candidates as r
class Phenotype(unittest.TestCase):
 def fixture(self,root):
  views=root/'views';(views/'masters').mkdir(parents=True);(views/'tables').mkdir();samples=['PUBLIC_1','PUBLIC_2','PUBLIC_3','PUBLIC_4'];labels={s:f'sample_{i:04d}' for i,s in enumerate(samples,1)}
  fields=list(dict.fromkeys([*r.REQUIRED,'phenotype_context_rank','target_gene_context_rank','gene_validity_rank','af_max_observed','clinvar_significance','GT','FT','DP','GQ']))
  fields+=[f'public_original_{i}' for i in range(151-len(fields))];self.assertEqual(len(fields),151);artifacts={};rowcounts={};counts={};allrows={}
  for index,sample in enumerate(samples):
   rows=[]
   specs=[('A','WRONG','frameshift_variant','1')]+([('A','GOOD','missense_variant','1')] if index==0 else [])+[('B','DISPUTED','stop_gained','1' if index==0 else '0')]
   for name,gene,term,eligible in specs:
    v=dict.fromkeys(fields,'');v.update(sample=sample,chrom='chr1',pos='10' if name=='A' else '20',ref='A',alt='T',gene=gene,hgnc_id='',transcript='PUBLIC_TX_'+gene,candidate_eligible=eligible,quality_rank='1' if eligible=='1' else '0',sample_state='alternate',consequence=term,variant_id='1:'+('10' if name=='A' else '20')+':A:T',mane_select='',clinvar_review_rank='2' if gene=='WRONG' else '0',consequence_rank='3',population_rank='0',predictors_damaging='0',AlphaMissense_research_rank='0',phenotype_context_rank='0',target_gene_context_rank='0',gene_validity_rank='0',GT='0|1',FT='PASS',DP='30',GQ='99',clinvar_significance='Pathogenic' if gene=='DISPUTED' else '',af_max_observed='.2' if gene=='DISPUTED' else '',public_original_0='α\tquoted\nline')
    rows.append(v)
   p=views/'masters'/(labels[sample]+'.master.tsv.gz');r.table(p,fields,rows);artifacts[p.relative_to(views).as_posix()]={'bytes':p.stat().st_size,'sha256':'a'*64};rowcounts[sample]=len(rows);counts[sample]={'eligible_alleles':2 if index==0 else 1};allrows[sample]=rows
  membership=[]
  for pos,eligible,qc in [(10,samples,samples),(20,samples[:1],samples[:3])]:membership.append(dict(allele_id=f'chr1:{pos}:A:T',chrom='chr1',pos=str(pos),ref='A',alt='T',eligible_sample_count=str(len(eligible)),eligible_samples_json=json.dumps(eligible),quality_passing_carriers_json=json.dumps(qc),raw_observed_alt_samples_json=json.dumps(samples),carrier_count=str(len(qc)),pct_of4=str(len(qc)*25),sample_calls_json=json.dumps({s:{'GT':'0|1','FT':'PASS'} for s in samples})))
  matrix=views/'tables/allele-membership-matrix.tsv.gz';r.table(matrix,list(membership[0]),membership);artifacts[matrix.relative_to(views).as_posix()]={'bytes':matrix.stat().st_size,'sha256':'b'*64};counts.update(cohort={'eligible_alleles':2},shared_ge2={'eligible_alleles':2},all4={'eligible_alleles':1})
  m=dict(status='PASS',synthetic=False,report_input_kind='FULLY_ANNOTATED_VCF_PACKETS_V1',samples=samples,sample_file_labels=labels,logical_full_fields=fields,artifacts=artifacts,sample_master_rows=rowcounts,master_transcript_rows=36,counts=counts);(views/'manifest.json').write_text(json.dumps(m))
  profiles=root/'clinical.json';profiles.write_text(json.dumps(dict(schema_version=1,created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),assertion_source='PUBLIC synthetic user-report control',samples=[dict(sample_id=s,diagnoses=[dict(label='Hypertrophic cardiomyopathy',mondo_id=r.MONDO,confirmation='confirmed_by_user')],phenotypic_features=[dict(hpo_id=r.HPO,status='present',source='user_report')]) for s in samples])))
  hpo=root/'gene_to_phenotype.tsv';hpo.write_text('#ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\n1\tGOOD\tHP:0001639\tHypertrophic cardiomyopathy\n2\tWRONG\tHP:OTHER\tOther\n')
  cg=root/'clingen.csv';cg.write_text('#PUBLIC metadata\nGENE SYMBOL,GENE ID (HGNC),DISEASE LABEL,DISEASE ID (MONDO),CLASSIFICATION,MOI\nGOOD,,familial hypertrophic cardiomyopathy,MONDO:0005045,Definitive,AR\nGOOD,,hypertrophic cardiomyopathy,MONDO:0005045,Strong,AD\nWRONG,,dilated cardiomyopathy,MONDO:OTHER,Definitive,AD\nDISPUTED,,hypertrophic cardiomyopathy,MONDO:0005045,Disputed,AD\n')
  return views,profiles,hpo,cg,root/'new-output',allrows
 def test_actual_spawn_full_tx_rerank_preserves_source_and_all_scope_semantics(self):
  with tempfile.TemporaryDirectory() as d:
   views,profiles,hpo,cg,out,original=self.fixture(Path(d));before={p:r.stamp(p) for p in (views/'masters').glob('*.gz')};dry=r.run(views,profiles,hpo,cg,out);self.assertEqual(dry['sample_ALT_master_rows'],9);self.assertFalse(out.exists())
   result=r.run(views,profiles,hpo,cg,out,workers=2,execute=True,max_hours=.1);self.assertEqual(result['scanned_sample_ALT_rows'],9);self.assertEqual(result['source_cohort_all_sample_rows'],36);self.assertEqual(result['scopes']['cohort']['source_eligible_alleles'],2);self.assertEqual(result['scopes']['shared_ge2']['source_eligible_alleles'],2);self.assertEqual(result['scopes']['all4']['source_eligible_alleles'],1)
   workbook=r.load(out/'workbook-data.json');top=workbook['scopes']['cohort'];self.assertEqual(len(top),1);self.assertEqual(top[0]['gene'],'GOOD');self.assertEqual(top[0]['hcm_clingen_validity_rank'],3);self.assertEqual(top[0]['gene_validity_rank'],'0');self.assertEqual(top[0]['phenotype_context_rank'],'0');self.assertEqual(top[0]['af_max_observed'],'');self.assertIn('SOURCE_AF_MISSING_NOT_ZERO',top[0]['hcm_manual_review_flags_json']);self.assertEqual(top[0]['ACMG_class'],'NOT_ASSIGNED')
   aux=out/'transcript-evidence/sample_0001.eligible-transcripts.jsonl';lines=[json.loads(x) for x in aux.read_text(encoding='utf-8').splitlines()];self.assertEqual(len(lines),3);self.assertEqual({x['gene'] for x in lines},{'WRONG','GOOD','DISPUTED'});self.assertIn('HIGH_SOURCE_AF',next(x for x in lines if x['gene']=='DISPUTED')['hcm_manual_review_flags_json'])
   for p,old in before.items():self.assertEqual(r.stamp(p),old)
   self.assertEqual(workbook['scopes']['sample_0002'],[])
 def test_negative_wrong_disease_unknown_phenotype_and_exact_alleles(self):
  with tempfile.TemporaryDirectory() as d:
   _,_,hp,cg,_,_=self.fixture(Path(d));hpo,clingen=r.resources(hp,cg);self.assertNotIn('WRONG',clingen)
   base=dict(gene='DISPUTED',hgnc_id='',af_max_observed='',clinvar_significance='');x=r.context(base,hpo,clingen,True);self.assertEqual(x['hcm_phenotype_priority'],0);self.assertEqual(x['hcm_clingen_validity_rank'],0)
   good=dict(base,gene='GOOD');self.assertEqual(r.context(good,hpo,clingen,False)['hcm_phenotype_priority'],0)
   self.assertFalse(r.belongs('all4',{'PUBLIC_1'},['PUBLIC_1','PUBLIC_2','PUBLIC_3']));self.assertTrue(r.belongs('all4',{'PUBLIC_1'},['PUBLIC_1','PUBLIC_2','PUBLIC_3','PUBLIC_4']))
   self.assertNotEqual(r.allele(dict(chrom='chr1',pos='1',ref='A',alt='T')),r.allele(dict(chrom='chr1',pos='2',ref='A',alt='T')))
if __name__=='__main__':unittest.main()
