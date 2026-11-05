# Единый enriched master и прежнее ранжирование

Дополнение к [представлениям v1](cohort-sample-views-v1.md). Новые
[exporter](../scripts/export_cohort_sample_views_v2.py) и
[continuation](../scripts/continue_cohort_sample_views_v2.py) не меняют замороженные
v1-файлы, primary stages6/7, tissue/population controllers или VEP.

Основной вход — один завершённый
`outputs/combined-annotation-master-v1/master_evidence.enriched.tsv.gz`:
все исходные строки и колонки07master, затем namespaced `allref_*` reference context,
availability и provenance. Клиническая классификация и индивидуальная экспрессия
по reference context не назначаются. MT context остаётся отдельным приложением,
не добавляется как ядерная genotype строка.

Основные выходы сохраняются прежними: один cohort master/top20/functional types,
затем четыре sample master/VCF+TBI/top20/types. Новые context колонки присутствуют
во всех этих master/top20/representative functional таблицах. Общая master —
byte-identical копия enriched master; исходный07master остаётся неизменённым и
закреплён отдельным SHA. Sample VCF получают единожды сделанные общие VEP CSQ;
reference context находится в единой joined TSV, повторный VEP/API не выполняется.

Portable revision v2.1 также сохраняет declared `applied-evidence/population/...`
long tables и separate MT appendix indexes рядом с представлениями. Каждая копия
связана точным SHA/bytes с joined completion receipt **и** population producer
receipt; только известные applied result tables/QC, без raw public databases/cache.
Относительные `allref_*` evidence pointers относятся к export root и разрешаются
offline. Исходные receipts и таблицы не переписываются; HTML добавляет appendix links.

Перед экспортом проверяются:

- Точный root, completed manifest schema1/status`PASS_COMBINED_ANNOTATION_MASTER`
  и собственный PASS workflow receipt со stage`combined-annotation-master-v1`.
- Base/enriched bytes/SHA, исходный07receipt и exact06/07+tissue/population
  receipt pins. Собственная receipt должна закреплять manifest и enriched artifact.
- Полный порядок original header, уникальные только appended `allref_*` колонки
  без конфликтов и сохранённые original/rank eligibility invariant declarations.
- **Каждая строка**, каждое исходное значение и исходный порядок сравниваются с
  оригинальным07master. Новый SHA не может легализовать изменённый GT/FT/score/rank.
- `invariant_digest` — SHA256 reconstructed original TSV: `csv.writer` с delimiter
  TAB, lineterminator LF, стандартными quoting/doublequote, ordered original header
  и field vectors, UTF-8. Embedded TAB/LF обрабатываются как значения, не split/join.

Reader ограничивает один TSV field64MiB. Работа остаётся streaming,≤2CPU и4GiB
address space, с исходным deadline. После проверки используется отдельный freshly
loaded instance **exact SHA-pinned v1 engine**: меняется только входной master и
список отображаемых appended context колонок. Файл v1 engine и sealed rank_key
не переписываются. Candidate eligibility и rank code выполняются без изменений.

Новый controller ждёт06/07 **и completed enriched master**. Один06/07PASS не
разрешает baseline-only экспорт. ROOT-reviewed rule закрепляет новый code pair,
старый engine SHA, sealed ranking SHA и isolated pysam0.24.1 runtime. Default —
status only; выполнить можно только с exact reviewed-rule SHA. Wait≤8часов и
absolute2026-11-06T03:46:48.046025UTC, без нового бюджета или process actuation.

```bash
/isolated/pysam0241/bin/python scripts/continue_cohort_sample_views_v2.py --rule /private/ROOT-reviewed-joined-views-v2-rule.json --reviewed-rule-sha256 EXACT_RULE_SHA --execute
```

Выход сначала создаётся в скрытой staging directory. Только после всех source
rechecks и полной записи manifest он атомарно переводится через Linux
`renameat2(RENAME_NOREPLACE)` в canonical `outputs/cohort-sample-views-v1`.
Существующий output не заменяется. Незавершённый staging не принимается backup
allowlist как completed output. Backup schema adapter для joined provenance
должен рассматриваться отдельным review; старый parser не объявляет новый
неизвестный contract поддержанным автоматически.

Synthetic проверки: quote-aware invariant, source mutation при заново вычисленных
hashes, namespace rejection, exact upstream pins, baseline06/07 insufficient,
review/runtime/base-code/deadline pins; Linux pysam integration сравнивает counts,
rank/variant/sample/transcript порядок с original v1 и проверяет добавленные
колонки в cohort/sample masters/top20/functional tables. Пациентский экспорт не
объявляется выполненным до собственного нового PASS manifest.
