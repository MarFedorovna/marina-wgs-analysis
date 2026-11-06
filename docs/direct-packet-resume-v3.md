# Продолжение прямой аннотации: v3

`scripts/direct_fully_annotated_vcf_v3.py` продолжает уже проверенный запуск v2 после превышения размера одного блока. VEP повторно не запускается. Итог остаётся один: дополненный общий VCF → таблицы, категории и top-20 общего набора и четырёх образцов.

Готовый блок используется только при наличии всех пяти файлов: `body.vcf.gz`, `batch-proof.json`, `full-vector-proof.json`, `annotation_counts.json`, `consequence_counts.json`. Проверяются SHA-256 тела, BGZF EOF, диапазон порядковых номеров и счётчики. Все сохранённые входные записи побайтно сверяются с зафиксированным завершённым VEP VCF. Только отсутствующие или незавершённые блоки делятся на задачи примерно по 5000 строк; одна аллель остаётся целиком.

Новый `resume-provenance.json` связывает старое проверенное правило, код, файлы блоков и новые номера блоков. Сохранённые тела и токены не переписываются. Поля, порядок строк, исходные записи и определения счётчиков сохраняются; создаются новые честные доказательства порядка. Квитанции старых стадий 06/07 не подделываются.

Подготовка требует приватных, уже проверенных путей и исходного плана четырёх образцов:

```bash
python3 scripts/direct_fully_annotated_vcf_v3.py prepare-resume \
  --parent-rule /PRIVATE/PARENT/ROOT-reviewed-rule.json \
  --parent-reviewed-rule-sha256 EXACT_PARENT_RULE_SHA256 \
  --old-work /PRIVATE/OLD_STATE/work-ID \
  --resume-public-proof /PRIVATE/CANDIDATE/public-resume-proof.json \
  --resume-cache-audit /PRIVATE/CANDIDATE/resume-cache-audit.json \
  --state-directory /PRIVATE/NEW_STATE \
  --rule /PRIVATE/CANDIDATE/rule-candidate.json
```

Запуск использует отдельное проверенное правило; без `--execute` команда только показывает план:

```bash
python3 scripts/direct_fully_annotated_vcf_v3.py run \
  --rule /PRIVATE/CANDIDATE/ROOT-reviewed-rule.json \
  --reviewed-rule-sha256 EXACT_RULE_SHA256 --execute
```

Скрипт сам не останавливает процессы. Сохраняются пределы: 24 процесса по одному CPU, общий бюджет памяти до 48 ГиБ, 64 МиБ текста на блок, 48 ГиБ общего текста и резерв диска 150 ГиБ. Исходный абсолютный срок VM не продлевается.

Публичная проверка намеренно вызвала превышение уменьшенного лимита, затем совместила шесть готовых блоков и две пересчитанные задачи. Все 4604 упорядоченных вектора, 33 дополнительных поля, исходные записи и счётчики совпали с исходным последовательным расчётом; старые готовые и частичные файлы остались неизменными.

`INTERMEDIATE_CORE33_PASS` означает готовность промежуточного VCF с 33 полями. Итоговый VCF публикуется после фактического SpliceAI overlay; таблицы получают его 13 полей из INFO того же итогового VCF. Наблюдение без чтения данных:

```bash
python3 scripts/monitor_pipeline.py --root /PROJECT_ROOT --once --human
```
