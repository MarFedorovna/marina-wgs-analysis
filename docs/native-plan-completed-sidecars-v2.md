# План локальной копии: завершённые таблицы и sample annotated VCF

Обновлён только [prepare_native_plan.py](../transfer/prepare_native_plan.py).
Активный native SFTP writer, managed SSH key, pin/known_hosts, локальная проверка
SIZE ONLY, исходные workflow receipts и завершённые результаты не изменяются.
Новые таблицы не считаются скачанными, пока отдельный downloader/importer их не
получил и не записал свою квитанцию.

## Что выбирается

Сохраняются4основных BAM+BAI, завершённые cohort VCF/indices, текущие финальные
`mitochondria-recovery-v2/.../calls` VCF/QC и существующие reports/provenance.
Добавлены `.tsv.gz`/`.csv.gz` и receipt-listed `.png`/`.jpg`/`.jpeg`/`.svg` plots,
необходимые для offline HTML. Изображения выбираются только из завершённых
разрешённых output receipts, не из ресурсов/датасетов. Четыре индивидуальных **annotated VCF+TBI** включаются;
индивидуальные **gVCF** и partition/checkpoint VCF остаются на сервере.

Новые sidecar результаты разрешены только из точных roots/schema:

| Root в outputs | Исходная квитанция и проверка |
|---|---|
| `consequence-tables-v1` | `manifest.json`, schema1/PASS/non-synthetic, exact reviewed exporter SHA, embedded completed06receipt и source transcript artifact pin |
| `tissue-protein-context-v1` | `manifest.json`, точный research-export status/keys, source06pin,06/07receipt SHA, exact public context-index SHA; только итоговая joined TSV |
| `population-reference-context-v1` | собственный `workflow-receipt.json`, schema1/PASS/stage/identity, population metadata status, source07pin и06/07receipt SHA |
| `cohort-sample-views-v1` | `manifest.json`, schema1/PASS/non-synthetic, exact sealed rank module SHA, пять source/receipt pins, четыре sample labels и восемь VCF/TBI artifacts |
| `combined-annotation-master-v1` | собственный schema1/PASS/stage workflow receipt, enriched master+manifest pins, исходный07master,06/07+tissue/pop receipts, namespaced columns и declared invariant/count |

В adapter v4 unified views принимаются только при
`report_input_kind=COMBINED_ANNOTATION_MASTER_V1`: exact joined manifest,
base/enriched/source pins, проверенные row/rank invariant flags и exact SHA
обоих reviewed v1-engine/v2-exporter files. Дополнительные source pins ограничены
двумя combined files, own receipt и двумя upstream sidecar receipts. Неизвестный
pin/resource отклоняется. Сам report exporter independently проверяет каждую
original строку; backup planner не повторяет этот проход и не переаннотирует данные.
Portable v2.1 views и combined master могут включать только declared
`applied-evidence/population/...` копии уже завершённых population result tables,
каждая с exact source_relative_path/role и двойным joined+population receipt
SHA/size binding. Это применённые таблицы доказательств, не публичные базы.
Поддерживаются их два небольших appendix indexes; неизвестный файл не разрешается.

Для каждой схемы есть отдельный allowlist конечных относительных путей. Произвольный
`manifest.json`/PASS внутри неизвестного dataset не допускается. Partial/failure
не выбираются. Traversal, symlink/parent symlink, resource/cache/raw пути,
неверный source/producer pin, digest format или изменённый size отклоняются.
Root-level completion manifest сохраняется как оригинал; нормализация artifacts
делается только в памяти. Primary workflow receipt не подменяется и не создаётся.

Public databases, reference FASTA/GTF/SQLite, raw API reference responses,
download caches и FASTQ не включаются. Указание reference path в provenance не
является командой копирования этого источника.

## Проверка и красивые папки

Expected SHA256 берётся из настоящей producer completion receipt; фактический
размер проверяется по stat. Большие BAM и таблицы **не перечитываются ради нового
SHA256**. Хешируются только небольшие receipts/manifests/code/context metadata.
Локальный SIZE ONLY режим продолжает хранить `ActualSha256=null`, отдельно от
исторических SHA-verified файлов. Равный размер не объявляется полным checksum PASS.

Server `RESULTS` помещает `.tsv.gz`/`.csv.gz` в `Reports`, sample VCF — в `VCF`.
Отдельный локальный `write_local_index.py` имеет свою category mapping; обновление
этой mapping следует рассматривать отдельно, чтобы compressed tables также
попали в локальный `Reports`. Данный patch не меняет importer/index/writer.

## Безопасное сравнение перед подключением наблюдателя

Без флагов выполнение остаётся dry run. `--compare` выполняет один pinned SSH
inventory и пишет **новый private plan**; не создаёт server hardlinks/RESULTS,
не меняет active plan и не запускает загрузку. Повторное использование существующего
output в compare режиме отклоняется. `--compare-with-plan` сравнивает только
ожидаемые path/size/SHA metadata и сохраняет counts added/removed/changed.

```bash
python3 transfer/prepare_native_plan.py --host CONFIRMED_SSH_ALIAS --ssh-config EXISTING_MANAGED_CONFIG --known-hosts EXISTING_PINNED_HOSTS --remote-root /srv/marina --workspace-root /home/ubuntu/workspace --jupyter-root /home/ubuntu --output /private/new-candidate-plan.json --compare --compare-with-plan /private/active-plan.json --deadline-utc 2026-11-06T03:46:48.046025+00:00
```

После ROOT review можно заменить только plan observer, загрузив новый HELPER в
новый bounded process. Старый уже работающий observer держит старую константу в
памяти. Native SFTP writer и его partial трогать не требуется. Watch≤8часов,
interval≥60секунд и исходный абсолютный VM deadline сохраняются.

Реальное read-only сравнение до завершения stages6/7:464готовых artifacts,
43receipts,4BAM,44current MTv2 artifacts,0gVCF;24added/269removed относительно старого
загруженного observer plan,0изменённых expected identities. Удаления из candidate
означают исключение старых/неподходящих источников из нового списка, а не удаление
файлов. Новые annotation sidecars на этот момент ещё ожидали завершения, поэтому
их будущая выборка подтверждается synthetic schema fixtures, а не заявлением о
готовой пациентской аннотации.

Тесты: четыре completed schemas, compressed Reports mapping, четыре sample VCF+TBI,
четыре BAM без повторного hashing, gVCF exclusion, pending/failure/unknown dataset,
resource smuggling, traversal, symlink, producer/context/source pin и size mismatch.
