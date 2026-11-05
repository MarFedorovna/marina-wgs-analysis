# Одна аннотация: общая таблица и четыре индивидуальных представления

Этот отдельный локальный экспорт использует уже завершённую общую аннотацию.
Он не повторяет VEP/calling, не меняет stages6/7 или существующие outputs, не
отправляет варианты/генотипы/фенотипы API и не устанавливает новые базы. Четыре
sample IDs берутся только из private VCF header; в публичном коде их нет.

## Основной результат

`index.html` открывает простой набор:

1. Общая полная master TSV.gz, общий top20 и функциональные категории.
2. Для каждого из4образцов: annotated VCF.gz+TBI, master TSV.gz со всеми
   transcript rows наблюдаемых ALT, индивидуальный top20 и функциональные категории.

Полная cohort master — byte-identical копия stage7 `master_evidence.tsv.gz`.
Per-sample master сохраняет каждую исходную transcript строку при наличии ALT
в GT, включая filtered и partial `1/.`. Наблюдаемый ALT не означает доказанную
качественную/клиническую находку. Reference/no-call строки остаются в полном
cohort master и membership matrix; отсутствие записи в sample VCF не доказывает
callability или клинически отрицательный результат.

В таблицах сохранены позиция/REF/ALT, gene/HGNC, transcript/protein identity,
аминокислотная замена, HGVS, consequence, exon/protein position, genotype/FT/DP/AD/GQ,
ClinVar significance/review/conflicts/snapshot, ClinGen/gene context, source AF,
SIFT/PolyPhen/CADD/AlphaMissense и raw local evidence, включая unaligned FATHMM.
Предиктор без подтверждённого transcript match не получает выдуманного scalar score.
Cohort INFO AF и gnomAD AF находятся в разных полях. All-transcript appendix
сохраняет также полные VEP CSQ и все исходные source/provenance columns.

## Ранжирование и ACMG

Используется непосредственно `rank_key()` из exact SHA-pinned sealed
`annotation_stages.py`; candidate eligibility берётся из исходного
`candidate_eligible==1`, а не рассчитывается новой эвристикой. Для каждого образца
и каждой уникальной аллели выбирается лучший eligible transcript; cohort top20
выбирает лучший eligible sample/transcript этой аллели. Дубликаты транскриптов не
занимают дополнительные места. При числе eligible аллелей<20 список не дополняется.

Это исследовательская приоритизация с существующими HCM/phenotype context flags,
не диагноз и не автоматическая ACMG classification. В каждой summary строке
`ACMG_class=NOT_ASSIGNED`, criterion evidence/review — manual pending.
Отдельный `ACMG-manual-top20-worksheet.tsv.gz` содержит критерии, пустые strength,
source/rationale и explicit unassessed status. ClinVar external class остаётся
отдельным источником, а не собственным заключением. PVS1/PP3/BP4/LoF/GoF/NMD не
назначаются по consequence labels, числу программ или числу носителей.

## Membership и дополнительные приложения

Идентичность shared allele — ровно уже нормализованные CHROM/POS/REF/ALT из
одного source VCF. Разные аллели в одном гене не объединяются. Duplicate exact
VCF keys, multiallelic input, GT allele>1 или расхождения VCF↔master отклоняются.
Нормализация повторно не запускается.

`carrier_count` — число полных quality-passing ALT GT (FT=PASS и FILTER=PASS).
`pct_of4` =25×carrier_count:25/50/75/100; знаменатель всегда4. Наблюдаемые raw ALT
и raw percentage показаны отдельно. `assessed_count` включает только complete
quality-passing ALT/reference; missing/partial/filtered/unknown-quality и
nonsequence placeholders не считаются assessed reference. Явные списки unknown,
filtered, partial и unassessed могут пересекаться и не должны суммироваться.
`sample_calls_json` хранит GT/DP/AD/GQ/FT каждого образца.

Shared≥2/all4 top20 и functional tables — вспомогательные приложения. Их top20
требует≥2/4 quality-passing carriers и хотя бы одной исходной eligible строки;
индивидуальная eligibility каждого образца остаётся в matrix. Raw shared
functional membership может включать partial/filtered ALT и потому имеет отдельные
QC/candidate flags. `gene-overlap-SEPARATE-NOT-SHARED-ALLELES.tsv.gz` отдельно
показывает sample/gene overlaps и unique allele counts, допуская разные аллели.

Functional tables считают уникальные аллели, не transcript rows. Категория означает
наличие соответствующего эффекта хотя бы в одном исходном транскрипте; категория
может отличаться от representative transcript. Поэтому есть union consequences,
явный membership basis и ссылка на полную master. Категории пересекаются;
pLoF_candidates — исследовательская группа, не доказанная потеря функции.
Функциональные raw/QC/candidate counts являются allele-level membership и не
доказывают эффект каждого транскрипта/ткани или индивидуальную клиническую тяжесть.

## Sample VCF: INFO и техническая проверка

Сохраняются исходные site FILTER, FORMAT, phase и sample column, все CSQ и другие
аннотации. Оригинальные cohort INFO AC/AN/AF/NS/DP переносятся в `COHORT_*` с
исходными header declarations. Новые AC/AN/AF/NS вычисляются из called GT tokens
единственного образца; missing tokens исключаются из AN. Для `1/.` получается
AC=1,AN=1,AF=1, но `VIEW_QUALITY_STATE=PARTIAL_UNKNOWN`; это не уверенная оценка
частоты/генотипа. Остальные site INFO сохраняют cohort/source происхождение.

Сжатие — настоящий BGZF через проверенный pysam0.24.1, индекс — настоящий TBI.
Каждый готовый файл открывается с indexed fetch; проверяются sample header,
число записей и derived AC/AN/AF каждой записи. Zero-record sample VCF допустим
и сохраняет корректные header/BGZF/TBI; он не является отрицательным заключением.

## Запуск только после review

[Exporter](../scripts/export_cohort_sample_views.py) требует06_annotate и07_report
PASS и точное совпадение SHA256/bytes annotated VCF+TBI и master с receipts.
Input/receipt/ranking hashes повторно проверяются перед PASS. При расхождении
или deadline остаётся INCOMPLETE без PASS; исходный workflow не изменяется.

```bash
/isolated/pysam0241/bin/python scripts/export_cohort_sample_views.py --workflow-root /srv/marina/outputs/workflow-v1 --output /srv/marina/outputs/cohort-sample-views-v1 --ranking-module /srv/marina/private/annotation-runtime-v1/scripts/annotation_stages.py --ranking-sha256 EXACT_REVIEWED_SEALED_SHA
```

[Continuation](../scripts/continue_cohort_sample_views.py) по умолчанию только
показывает статус. Для выполнения нужен ROOT-reviewed exact rule SHA, обе code SHA,
sealed ranking SHA и проверенная isolated Python/pysam identity. Он ждёт06/07PASS,
держит собственный lock, повторно проверяет pins и вызывает только exporter.
Никаких VEP/container/controller start/stop/kill/restart нет. Изменять текущий
live annotation процесс ради этого экспорта не требуется.

```bash
/isolated/pysam0241/bin/python scripts/continue_cohort_sample_views.py --rule /private/ROOT-reviewed-views-rule.json --reviewed-rule-sha256 EXACT_REVIEWED_RULE_SHA --execute
```

Ограничения:≤2CPU affinity,4GiB address space, SQLite64MiB cache и disk-backed
working index; transcript evidence агрегируется по одной аллели, top queues≤20.
Sidecar сохраняет≥64GiB свободного диска плюс размер полного master appendix;
при нарушении останавливается. Wait≤8часов, с исходным абсолютным VM deadline
2026-11-06T03:46:48.046025UTC. Новый бюджет не создаётся. Raw TSV следует импортировать
как текст; исходные clinical/evidence strings не переписываются ради Excel.

Тесты содержат только synthetic IDs: реальный BGZF/TBI fetch, recalculated partial
GT counts/cohort INFO provenance, reference/missing/filtered distinction, exact
allele sharing vs same gene, four sample masters, all transcript preservation,
direct parity with sealed ranking, no padding, unknown ACMG, hash/GT/PASS gates,
ROOT review/code pins и остановка по deadline до чтения/создания outputs.
