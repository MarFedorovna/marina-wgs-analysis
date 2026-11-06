import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

// Bundled artifact-tool dependency is supplied by the documented runtime loader.
// Patient input JSON and resulting workbooks stay outside the code repository.
const [input, output, previewDir, mtReviewPath, readReviewPath] = process.argv.slice(2);
if (!input || !output || !previewDir) throw new Error('Usage: node build_phenotype_workbook.mjs workbook-data.json output.xlsx preview-directory');
const data = JSON.parse(await fs.readFile(input, 'utf8'));
const wb = Workbook.create();
const overview = wb.worksheets.add('Обзор');
const names = {cohort:'Общий топ20',shared_ge2:'Общие от двух',all4:'Общие все четыре',sample_0001:'Образец 1',sample_0002:'Образец 2',sample_0003:'Образец 3',sample_0004:'Образец 4'};
const columns = [
  ['rank','Ранг','int',8],['gene','Ген','text',15],['allele_id','Вариант GRCh38','text',34],
  ['hgvs_c','HGVS c','text',45],['hgvs_p','HGVS p','text',40],['consequence','Последствие','text',36],
  ['hcm_priority_reason','Основание приоритета ГКМП','text',56],
  ['hcm_hpo_match','Совпадение HPO','int',18],['hcm_clingen_validity_rank','ClinGen контекстный ранг','int',20],
  ['carrier_count','QC носителей из 4','int',19],['pct_of4','Доля из 4','pct',15],
  ['quality_passing_carriers_json','QC носители','text',50],['eligible_samples_json','Образцы с кандидатом','text',50],
  ['sample','Представитель строки','text',20],['GT','GT представителя','text',18],['DP','DP представителя','int',18],['GQ','GQ представителя','number',18],
  ['af_max_observed','Макс. популяционная AF','number',22],['af_status','Статус AF','text',48],
  ['clinvar_significance','ClinVar','text',48],['clinvar_conflicts','Конфликты ClinVar','text',45],['clinvar_conditions','Состояния ClinVar','text',60],
  ['transcript','Транскрипт','text',29],['mane_select','MANE','text',24],
  ['SIFT','SIFT','text',35],['PolyPhen','PolyPhen','text',35],['FATHMM','FATHMM','text',35],['CADD_PHRED','CADD PHRED','number',18],
  ['AlphaMissense_score','AlphaMissense','number',19],['AlphaMissense_match_status','Статус AlphaMissense','text',55],
  ['spliceai_status','Статус SpliceAI','text',50],['spliceai_DS_AG','SpliceAI AG','number',17],['spliceai_DS_AL','SpliceAI AL','number',17],['spliceai_DS_DG','SpliceAI DG','number',17],['spliceai_DS_DL','SpliceAI DL','number',17],
  ['ACMG_class','Класс ACMG','text',22],['classification_status','Статус интерпретации','text',55],
];
function scalar(v,type) {
  if (v === undefined || v === null || v === '' || v === '.') return null;
  if (typeof v === 'boolean') return Number(v);
  if (type === 'text') {const s=typeof v==='object'?JSON.stringify(v):String(v); return s.startsWith('=')?"'"+s:s;}
  const n=Number(v);if (!Number.isFinite(n)) return String(v);
  return type==='pct'?n/100:n;
}
function formatSheet(sh, matrix, widths, tableName) {
  const n=matrix.length,k=matrix[0].length;
  const used=sh.getRangeByIndexes(0,0,n,k);used.values=matrix;
  used.format.font={name:'Arial',size:10};used.format.verticalAlignment='center';used.format.rowHeight=24;
  const head=sh.getRangeByIndexes(0,0,1,k);head.format.fill='#20364D';head.format.font={name:'Arial',size:10,bold:true,color:'#FFFFFF'};
  head.format.rowHeight=40;head.format.wrapText=true;head.format.horizontalAlignment='center';
  for(let c=0;c<k;c++) sh.getRangeByIndexes(0,c,n,1).format.columnWidth=widths[c]||24;
  sh.showGridLines=false;sh.freezePanes.freezeRows(1);sh.freezePanes.freezeColumns(Math.min(3,k));
  let letters='',cnum=k;while(cnum>0){cnum--;letters=String.fromCharCode(65+cnum%26)+letters;cnum=Math.floor(cnum/26);}
  if(n>1){const t=sh.tables.add(`A1:${letters}${n}`,true,tableName);t.showFilterButton=true;}
}
const overviewRows=[['Область','Число кандидатов','Строк в топе','Что означает'],
 ['Результат',null,null,'Исследовательская очередь для проверки, без автоматического класса ACMG.'],
 ['Фенотип',null,null,'ГКМП подтверждена пользователем у всех четырёх; HP:0001639.'],
 ['Сохранность',null,null,'Исходные VCF, master и прежние топы не перезаписываются.'],
 ['HPO',null,null,'Точное совпадение термина с геном; не семантическая оценка и не доказательство причинности.'],
 ['ClinGen',null,null,'Уровень связи гена с ГКМП. Наследование и механизм конкретного варианта требуют проверки.'],
 ['Общие варианты',null,null,'Совпадают CHROM, POS, REF, ALT; QC-присутствие и кандидатность по образцам различаются.'],
 ['Пропуски',null,null,'Пустая оценка или AF означает отсутствие данных, а не ноль.'],
 ['Полная аннотация',null,null,'Полные транскриптные доказательства находятся в TSV/JSONL рядом с отчётом и исходном VCF.']];
for(const [scope,title] of Object.entries(names)) {
  const raw=data.scopes[scope]||[];const rows=Array.isArray(raw)?raw:(raw.top20||[]);
  const sh=wb.worksheets.add(title);sh.tabColor='#365D82';
  const matrix=[columns.map(c=>c[1]),...rows.map(row=>columns.map(([key,,type])=>scalar(row[key],type)))];
  formatSheet(sh,matrix,columns.map(c=>c[3]),'Top_'+scope);
  for(let c=0;c<columns.length;c++) if(matrix.length>1 && columns[c][2]!=='text') {
    sh.getRangeByIndexes(1,c,matrix.length-1,1).setNumberFormat(columns[c][2]==='pct'?'0.0%':columns[c][2]==='int'?'#,##0':'0.########');
  }
  const counts=data.summary?.scopes?.[scope]||{};
  overviewRows.push([title,Number.isFinite(Number(counts.source_eligible_alleles))?Number(counts.source_eligible_alleles):null,rows.length,`Положительный контекст ГКМП: ${counts.hcm_context_positive_alleles??'не указан'}. Полный список и функциональные категории сохранены отдельно.`]);
}
if (mtReviewPath) {
  const mt=JSON.parse(await fs.readFile(mtReviewPath,'utf8'));
  if(mt.status!=='PASS_LOCAL_DESCRIPTIVE_MT_HCM_AUDIT') throw new Error('Closed descriptive MT review required');
  const samples=Object.keys(mt.alleles[0]?.sample_original_AF_HL_and_call_context||{});
  const head=['Вариант GRCh38','Прямые гены','Последствия','HGVS p','Прямые HPO-гены','Site PASS из 4','Site PASS образцы','Популяционный контекст',...samples.flatMap(s=>[`${s} AF`,`${s} FILTER`])];
  const rows=mt.alleles.map(a=>{
    const statuses=[...new Set(a.existing_frequency_contexts.map(x=>`${x.dataset}: ${x.match_status}`))].join('; ');
    return [a.allele_id,a.direct_genes.join(', '),a.direct_consequences.join(', '),[...new Set(a.direct_transcript_effects.map(e=>e.hgvs_p_raw).filter(Boolean))].join('; '),a.exact_HPO_HP0001639_direct_genes.join(', '),a.literal_site_PASS_count,a.literal_site_PASS_samples.join(', '),statuses,
      ...samples.flatMap(s=>{const c=a.sample_original_AF_HL_and_call_context[s];const fs=c.exact_ALT_fractions.map(x=>x.AF.exact_allele_raw).filter(x=>x!==null);return [fs.length===1?Number(fs[0]):fs.join('; ')||null,c.original_records.map(x=>x.FILTER).join('; ')||'NOT_OBSERVED'];})];
  });
  const sh=wb.worksheets.add('MT контекст');formatSheet(sh,[head,...rows],[34,28,40,50,30,18,55,70,...samples.flatMap(()=>[20,40])],'MT_Context');
  for(let c=8;c<head.length;c+=2)sh.getRangeByIndexes(1,c,rows.length,1).setNumberFormat('0.########');
  overviewRows.push(['MT контекст',mt.summary.unique_literal_PASS_alleles,null,'Отдельный описательный разбор site-PASS; AF сопоставлен точному ALT. HL сохраняется как пропуск. Ранг и ACMG не назначаются.']);
}
if (readReviewPath) {
  const review=JSON.parse(await fs.readFile(readReviewPath,'utf8'));
  if(review.status!=='PASS_LOCAL_INDEXED_SNV_READ_SUPPORT_REVIEW') throw new Error('Completed indexed read-support report required');
  const selected=review.results.filter(r=>JSON.parse(r.source_membership.eligible_samples_json).includes(r.sample));
  const headers=['Образец','Ген','Вариант GRCh38','Исходный GT','Исходный AD','Исходный DP','REF чтения','ALT чтения','REF шаблоны','ALT шаблоны','ALT вперёд','ALT обратно','Доля ALT чтений','Флаги просмотра'];
  const matrix=[headers,...selected.map(r=>{const c=r.counts;const total=c.reads_REF+c.reads_ALT+c.reads_OTHER;
    return [r.sample,r.gene,`${r.chrom}:${r.pos}:${r.ref}:${r.alt}`,r.source_call.GT,r.source_call.AD,scalar(r.source_call.DP,'int'),c.reads_REF,c.reads_ALT,c.templates_REF,c.templates_ALT,c.reads_ALT_forward,c.reads_ALT_reverse,total?c.reads_ALT/total:null,r.descriptive_flags.join('; ')||'Нет описательных флагов'];})];
  const sh=wb.worksheets.add('Чтения кандидатов');formatSheet(sh,matrix,[20,16,35,18,20,18,17,17,19,19,18,18,20,65],'Candidate_Read_Support');
  if(selected.length) sh.getRangeByIndexes(1,12,selected.length,1).setNumberFormat('0.0%');
  overviewRows.push(['Проверка чтений',selected.length,null,`MAPQ≥${review.parameters.mapq}, BaseQ≥${review.parameters.baseq}; повторяющиеся перекрывающиеся mates учтены по шаблонам. Полные четыре образца находятся в read-support TSV/JSON.`],
    ['Ограничение подсчёта',null,null,'Это отдельный фильтрованный подсчёт без BAQ и локальной сборки. Он не обязан совпадать с GATK DP/AD; описательные флаги не являются решением об артефакте или патогенности.']);
}
overviewRows.push(['Источники',null,null,'HPO 2026-09-01; ClinGen Gene–Disease Validity snapshot 2026-11-05; аннотации исходного завершённого набора.'],
 ['HPO документация',null,null,'https://obophenotype.github.io/human-phenotype-ontology/annotations/introduction/'],
 ['Термин ГКМП',null,null,'https://www.ncbi.nlm.nih.gov/medgen/2881'],
 ['Клиническая интерпретация',null,null,'https://www.clinicalgenome.org/tools/clingen-variant-classification-guidance/']);
formatSheet(overview,overviewRows,[27,24,20,125],'Report_Overview');
overview.getRangeByIndexes(1,3,overviewRows.length-1,1).format.wrapText=true;
overview.getRangeByIndexes(1,0,overviewRows.length-1,4).format.rowHeight=38;
overview.tabColor='#20364D';
wb.recalculate();
await fs.mkdir(path.dirname(output),{recursive:true});await fs.mkdir(previewDir,{recursive:true});
const inspect=await wb.inspect({kind:'table',range:"'Общий топ20'!A1:I6",include:'values',tableMaxRows:6,tableMaxCols:9,maxChars:3000});
await fs.writeFile(path.join(previewDir,'inspection.ndjson'),inspect.ndjson);
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},maxChars:2000});
await fs.writeFile(path.join(previewDir,'errors.ndjson'),errors.ndjson);
for(const name of ['Обзор',...Object.values(names),...(mtReviewPath?['MT контекст']:[]),...(readReviewPath?['Чтения кандидатов']:[])]) {
 const blob=await wb.render({sheetName:name,range:name==='Обзор'?'A1:D11':'A1:F8',scale:1,format:'png'});
 await fs.writeFile(path.join(previewDir,name+'.png'),new Uint8Array(await blob.arrayBuffer()));
}
await (await SpreadsheetFile.exportXlsx(wb)).save(output);
console.log(JSON.stringify({status:'EXPORTED_PHENOTYPE_REVIEW_WORKBOOK',sheets:8+Number(Boolean(mtReviewPath))+Number(Boolean(readReviewPath)),output}));
