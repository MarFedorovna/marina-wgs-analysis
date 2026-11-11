# Финальная сверка аннотаций в Excel

Исходный Excel создают скрипты [подготовки полных списков](../scripts/prepare_full_cvd_views.py) и [экспорта](../scripts/write_full_cvd_xlsx.py). Этап сверки принимает готовую книгу с листами `Все_прошедшие`, `X_Y_все`, `MT_все`, `Missense` и `Полнота_аннотаций`. Входные данные и все файлы сверки хранятся локально. Для запуска нужны Python 3.10+ и `openpyxl==3.1.5`.

```bash
python -m pip install openpyxl==3.1.5
```

## FATHMM на листе Missense

[Скрипт](../scripts/fill_fathmm_missense.py) выбирает все строки `Missense` с пустой оценкой, подтверждает аллели по QC VCF и записывает статус каждого запроса в TSV и на лист `FATHMM_Missense`. Используются код модуля OpenCRAVAT 2.4.1 и соответствующая база данных v2.3 из OpenCRAVAT Store. Передайте распакованную `FATHMM.db` с проверенной контрольной суммой и официальный архив кода, в котором находится `fathmm.py`.

```bash
python scripts/fill_fathmm_missense.py \
  --input /private/report.xlsx --qc-vcf /private/cohort.quality.vcf.gz \
  --database /private/FATHMM.db --code-archive /private/fathmm-2.4.1.code.zip \
  --output /private/report.fathmm.xlsx --audit-tsv /private/fathmm.audit.tsv
```

Для одного аллеля возможны оценки по нескольким транскриптам. Они сохраняются вместе с идентификаторами. Поле оценки для транскрипта из Excel заполняется при совпадении идентификатора Ensembl. Если модуль возвращает ошибку, скрипт проверяет числовые записи в локальной базе и отдельно отмечает происхождение восстановленной оценки.

## Сверка всех листов

[Сборщик](../scripts/reconcile_workbook.py) сопоставляет варианты по хромосоме, позиции, REF и ALT, с обработкой эквивалентных представлений малых indel при однозначном совпадении. Для транскриптных оценок учитываются ген и идентификатор транскрипта. Для полей качества, BAM и MT используется соответствующий образец. Поля DS и DP SpliceAI переносятся единым набором, когда он присутствует в источнике. Доля ALT и исходное поле HL в MT сохраняются раздельно. Если подходящего значения нет, поле остаётся пустым. Пропуски и источники видны на листах `Полнота_по_столбцам` и `Источники_дополнений`. Порядок строк и исходные ранги сохраняются.

Подготовьте локальный JSON-манифест. Относительные пути отсчитываются от папки манифеста. Для индексированного VCF ожидается соседний индекс `.tbi`. Для X/Y и MT используйте TSV.GZ, сформированные скриптом [подготовки веток](../scripts/prepare_mt_xy_annotation_inputs.py). Все ресурсы и результаты должны использовать согласованную сборку GRCh38 и rCRS для MT.

```json
{
  "annotated_vcf": "cohort.fully-annotated.vcf.gz",
  "autosome_qc_vcf": "cohort.quality.vcf.gz",
  "xy_qc_vcf": "XY.quality.vcf",
  "xy_annotations": "XY.quality.annotations.tsv.gz",
  "mt_annotations": "MT.site_PASS.annotations.tsv.gz",
  "hpo": "genes_to_phenotype.txt",
  "clingen": "clingen-gene-validity.csv",
  "opencravat": [
    {"label": "OpenCRAVAT chr1–22", "sqlite": "autosomes.sqlite"},
    {"label": "OpenCRAVAT X/Y", "sqlite": "xy.sqlite"},
    {"label": "OpenCRAVAT MT", "sqlite": "mt.sqlite"}
  ],
  "fathmm_database": "FATHMM.db",
  "fathmm_code": "fathmm-2.4.1.code.zip",
  "bam_support_tsv": ["cohort.read-support.tsv"],
  "refseq_gene_locations": {}
}
```

Можно передать несколько OpenCRAVAT SQLite для одной ветки в порядке предпочтения источников. `fathmm_database` и `fathmm_code` задаются парой, если требуется повторная локальная сверка FATHMM для X/Y. `bam_support_tsv` принимает выходы [скрипта поддержки SNV чтениями](../scripts/review_candidate_snv_reads.py). Необязательный `refseq_gene_locations` задаёт для символа гена массив из названия хромосомы и ссылки на проверенную запись NCBI Gene.

```bash
python scripts/reconcile_workbook.py run \
  --input /private/report.fathmm.xlsx --manifest /private/reconcile.local.json \
  --output /private/report.reconciled.xlsx --summary /private/completeness.tsv \
  --audit /private/reconciliation.audit.json --work-dir /private/reconcile-work

python scripts/reconcile_workbook.py verify \
  --input /private/report.fathmm.xlsx --output /private/report.reconciled.xlsx \
  --audit /private/reconciliation.audit.json
```

Проверка сверяет заполненные ячейки с журналом и сравнивает остальные значения с исходной книгой. JSON журнала, промежуточные файлы, итоговый Excel и TSV содержат сведения об исследовании и хранятся вне публичного репозитория. Счётчик заполненных ячеек включает повторения варианта на нескольких листах.
