# Паспорт Analysis 2

Статус: независимый новый анализ; фактическая завершённость определяется только PASS receipts и проверенными выходами.

## Предусмотренная цепочка

1. FASTQ intake и SHA-256.
2. Parabricks alignment/sort/markdup/QC.
3. GATK-compatible gVCF по каждому образцу.
4. GenomicsDBImport и GenotypeGVCFs.
5. Нормализация и split multiallelic records.
6. Сохранение полного VCF с FILTER/FT и отдельной quality subset.
7. Offline VEP и локальные базы/предикторы.
8. Пообразцовая и cohort-приоритизация.
9. Отдельный evidence ledger ACMG/AMP.
10. Отчёты и supplementary tables.

## Заявленные рабочие counts

- raw cohort: 864905;
- normalized/split: 907845;
- quality subset: 170822.

Эти значения являются ожидаемыми до тех пор, пока не пересчитаны из скопированных VCF и не связаны с PASS receipts.

## Отдельные ветки

- chrX/chrY с подтверждённой плоидностью;
- mtDNA с heteroplasmy/NUMT/circular-reference QC;
- coding-region coverage;
- CNV/SV/repeat expansions при наличии подходящих входов.

## Запрет смешивания

Старые VCF и прежние 29 вариантов не используются для фильтрации, top-10/top-20 или ACMG нового анализа. После фиксации независимого результата допускается отдельная concordance table.

