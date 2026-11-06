# Marina WGS: воспроизводимый исследовательский пайплайн

Репозиторий содержит только исполняемые скрипты и это описание. Пациентские FASTQ/BAM/VCF, таблицы, отчёты, кэши баз, токены и тестовые файлы сюда не входят и не должны добавляться.

## Что делает пайплайн

Пайплайн рассчитан на четыре образца WGS человека (GRCh38) и разделён на независимые этапы. Каждый этап принимает явно заданные входы и пишет результат в отдельный каталог. Пути к данным и секреты задаются только в приватном конфигурационном файле или переменными окружения.

## Этапы от входа до отчёта

### 0. Контроль входов и воспроизводимость

1. Создать приватный `samples.tsv`: `sample_id`, `fastq_r1`, `fastq_r2`, ожидаемые контрольные суммы.
2. Проверить контрольные суммы, место на диске, версии инструментов, референс и доступ к базам.
3. Не помещать `samples.tsv`, FASTQ/BAM/VCF, токены и локальные конфиги в Git.

Основные скрипты: `configure_workflow.py`, `configure_pilot.py`, `preflight-server.sh`, `resource_status_snapshot.py`, `monitor_pipeline.py`.

### 1. Референс GRCh38

Скачать FASTA GRCh38 и проверить `FAI`, sequence dictionary, ALT/decoy и индексы BWA-MEM2. Для всех образцов используется один зафиксированный референс и его SHA-256.

Основные скрипты: `prepare_reference.py`, `reference_access_snapshot.py`, `probe_grch38_xy_reference.py`.

### 2. Выравнивание и BAM/QC

FASTQ выравниваются BWA-MEM2; при наличии совместимого GPU допускается эквивалентный Parabricks-этап. Затем выполняются сортировка, маркировка дубликатов, индексация BAM и базовая статистика покрытия/дубликатов. Результаты этого этапа не публикуются в репозитории.

Основные скрипты: `run_parabricks.py`, `workflow.py`, `calling_stages.py`, `gene_coverage.py`, `joint_scatter.py`.

### 3. Герминальный вызов вариантов

Для каждого образца строится gVCF (GATK HaplotypeCaller или эквивалент Parabricks), затем выполняются совместная генотипизация, нормализация REF/ALT и индексация. Для SNP/indel применяется один и тот же GRCh38 reference dictionary.

Основные скрипты: `calling_stages.py`, `workflow.py`, `parallel_annotation.py`.

### 4. Фильтрация качества

Сначала сохраняются исходные значения, затем формируется исследовательский PASS-набор. Использованные пороги:

- генотип: `DP >= 10`, `GQ >= 20`;
- сайт: `QUAL >= 30`;
- SNP: `QD >= 2`, `FS <= 60`, `SOR <= 3`, `MQ >= 40`, `MQRankSum >= -12.5`, `ReadPosRankSum >= -8`;
- indel: `QD >= 2`, `FS <= 200`, `ReadPosRankSum >= -20`;
- отсутствующие QC-аннотации не проходят молча: они помечаются для ручной проверки.

Это исследовательский фильтр, а не автоматическая клиническая классификация.

### 5. Первичная аннотация

Основной слой — Ensembl VEP для GRCh38 с зафиксированным cache. В VCF/master сохраняются исходный вариант, транскрипт, HGVS, CONSEQUENCE, ген, аллельные частоты и идентификаторы источников. Неизвестные или отсутствующие значения не заменяются догадками.

Основные скрипты: `annotation_stages.py`, `parallel_annotation.py`, `continue_fully_annotated_vcf.py`, `prepare_annotation_resources.py`, `prepare_all_resources.py`.

### 6. Источники аннотации и их назначение

- **ClinVar** — опубликованные клинические утверждения и уровень согласия; используется как evidence, не как безусловный диагноз.
- **gnomAD exomes/genomes** — популяционные AF/AC/AN и контроль редкости.
- **dbSNP** — идентификаторы вариантов.
- **dbNSFP** — агрегированные функциональные предикторы, когда доступен соответствующий public snapshot.
- **HGNC** — стабильные символы и идентификаторы генов.
- **MANE Select/Plus Clinical и RefSeq** — клинически согласованные транскрипты.
- **HPO** — связь признаков с генами; в HCM-контексте используется подтверждённый фенотип `HP:0001639`.
- **ClinGen Gene–Disease Validity, Dosage и VCEP** — сила gene–disease связи и gene-specific правила; слабые/спорные связи не получают положительный boost.
- **Orphanet, PubMed** — контекст заболевания и литература для ручной проверки.
- **UniProt, GTEx/pext** — белковый и тканевой контекст, не самостоятельное доказательство патогенности.
- **AlphaMissense, SIFT, PolyPhen, CADD, REVEL, FATHMM** — вычислительные предикторы; их баллы показываются отдельно и не заменяют ACMG evidence.
- **SpliceAI** — отдельный слой splice-предсказаний; версия, модель и статус должны быть записаны в manifest.
- **MITOMAP/gnomAD-MT/HelixMTdb** — только для mtDNA и только в MT-ветке.

HGMD, OMIM, GeneReviews, LOVD и VarSome не импортируются автоматически без соответствующей лицензии/API-доступа. Если источник недоступен, это явно указывается в manifest.

### 7. Исследовательский отбор кандидатов

Кандидат сохраняется, если он является допустимым альтернативным генотипом (`GT` с ALT, `FT=PASS`, сайт PASS) и выполняется хотя бы одно условие:

1. ClinVar содержит `Pathogenic`/`Likely pathogenic`; или
2. вариант функциональный и редкий: `max available gnomAD AF <= 0.01`; или
3. функциональный вариант не имеет AF — он сохраняется как `AF unavailable` для ручной проверки.

Функциональный набор включает loss-of-function (stop-gained, frameshift, splice donor/acceptor), start-lost, stop-lost, missense, in-frame, protein-altering и splice-region варианты. Причина включения, AF и исходные фильтры сохраняются в каждой строке.

### 8. Общие, индивидуальные и функциональные таблицы

Из одного общего master строятся:

- cohort: все четыре образца;
- `shared_ge2`: варианты минимум в двух образцах;
- `all4`: варианты во всех четырёх;
- отдельный scope каждого образца;
- top-20 для каждого scope;
- отдельные категории: `stop_gained`, `stop_lost`, `start_lost`, `frameshift`, `splice_donor_acceptor`, `splice_region`, `missense`, `synonymous`, `inframe_insertion`, `inframe_deletion` и прочие consequence-классы.

Top-20 — это исследовательское ранжирование, а не утверждение причинности. Для HCM применяется отдельный контекстный rerank по HPO/ClinGen; автоматический ACMG-класс не присваивается.

Основные скрипты: `continue_cohort_sample_views_v3.py`, `export_cohort_sample_views_v3.py`, `export_consequence_tables.py`, `rerank_phenotype_candidates.py`, `build_phenotype_workbook.mjs`.

### 9. X/Y и митохондриальная ветки

X/Y обрабатываются с отдельным QC покрытия и генотипов; XX/XY не объявляются клиническим кариотипом без независимого подтверждения. mtDNA анализируется отдельно: собственные правила гетероплазмии/AF, MITOMAP/gnomAD-MT-контекст и функциональные категории. Диплоидные autosomal-пороги нельзя бездумно переносить на MT.

Основные скрипты: `run_xy_germline.py`, `technical_xy_qc.py`, `sex_chromosome_plan.py`, `annotate_mitochondrial_consequences.py`, `audit_mitochondria.py`, `mt_phenotype_context_review.py`.

### 10. Поддержка чтениями и ручная интерпретация

Для выбранных SNV локально проверяются ALT-чтения, MAPQ/BaseQ, число независимых шаблонов и баланс цепей; затем выполняется визуальная проверка в IGV. Это техническая проверка, не замена ACMG/AMP. Финальный ACMG разбор делает исследователь по ClinGen/ACMG evidence codes (например, PVS1, PS, PM, PP, BS, BP) с сохранением ссылки на источник и ручного решения.

Основные скрипты: `review_candidate_snv_reads.py`, `research_review_contracts.py`, `validate_annotation_controls.py`.

## Минимальный запуск

Команды ниже — шаблон. Они требуют приватных путей и настроек; пациентские пути намеренно не указаны.

```bash
bash scripts/preflight-server.sh /path/to/work
python3 scripts/prepare_reference.py --lock /private/config/reference.lock.json --root /private/work/references --execute
python3 scripts/prepare_all_resources.py --root /private/work --execute
python3 scripts/workflow.py --help
python3 scripts/annotation_stages.py --help
python3 scripts/export_cohort_sample_views_v3.py --help
python3 scripts/review_candidate_snv_reads.py --help
```

На практике shell-файлы запускаются через `bash`, а Python — через `python3`. Перед полным запуском нужно проверить `--help` конкретной версии скрипта и записать версии всех инструментов.

## Ограничения

Скрипты формируют воспроизводимый исследовательский конвейер. Они не выдают медицинское заключение, не заменяют подтверждение варианта, родословную, валидацию лаборатории, ACMG-review или CNV/SV-специфический pipeline. В публичную копию не включаются результаты, идентификаторы пациентов, BAM/VCF, базы данных и секреты.
