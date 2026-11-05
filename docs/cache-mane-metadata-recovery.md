# Исправление проверки MANE в официальном VEP115 cache

Официальный GRCh38 cache был полностью получен и распакован, но прежний installer
ошибочно требовал строку `source_MANE` с версией 1.4. В фактическом `info.txt` этой
строки нет. Это software assertion failure, не повреждение архива и не доказательство
отсутствия всех runtime MANE labels. Отдельный MANE GRCh38 summary1.4 уже получен;
его нельзя выдавать за заявленную версию данных внутри cache.

`cache_metadata_records()` теперь проверяет точные TAB-separated source keys:
species/assembly, GRCh38.p14, GENCODE49, SIFT6.2.1, PolyPhen2.2.3, обе gnomAD4.1 AF
компоненты и dbSNP156. Missing/changed non-MANE metadata и duplicate keys дают отказ.
Если source_MANE отсутствует, `mane_cache.status` равен
`source_version_unreported_in_cache_metadata`, `version=null`. Если версия явно
сообщена, сохраняется она сама. Внешний MANE summary остаётся отдельным ресурсом.
Применимость transcript labels и сопоставление с summary проверяются публичными
VEP controls; этот repair их не объявляет выполненными.

Resume выполняется отдельным проверяемым helper с `--execute`:

```bash
python3 -u /srv/marina/private/reviewed-annotation-resource-recovery-v1/scripts/recover_annotation_extracted_stage.py \
  --config /srv/marina/private/annotation-installer-v3/config/annotation-resources.example.json \
  --resume-extracted-stage /srv/marina/resources/annotation/.partial-annotation-7691702c9fc446c0b77eecec34752ba2 \
  --archive-sha256 5871a5e34527cce76f7ac75399c33e306bf3be8083c63fe7d3459c56338f98de \
  --execute
```

До запуска проверяются: owned stage path, byte-identical base `resources.lock.json`,
существующий target без перезаписи, свободный резерв 150GiB и pinned VEP image.
Helper проверяет publisher BSD checksum/bytes и полный archive SHA, затем сравнивает
каждый распакованный файл с соответствующим tar member по SHA/размеру. Links/special
paths, лишние или отсутствующие files и изменённые small-source receipts отклоняются.
Создаётся настоящее validation receipt; только после всех проверок base bundle
продвигается через atomic no-replace. Обычный corrected final preparer повторно
проверяет base artifacts и переиспользует hardlinks/cache и готовые small sources.
Архив не скачивается заново и cache повторно не распаковывается.

Exact recipe bytes/SHA сохранены: final config
`0d65e5ca1c15cd16def5518e5ff528d2ab052da7a26fb6414046a5e77396fa76`.
Исходные 15 sealed runtime files, binding/rule и контроллеры не изменены. Исправление
не создаёт искусственный PASS и не запускает анализ пациентов. После реального core
manifest PASS прежний conditional controller отдельно выполняет public controls и
только затем этапы6/7.

Текущий repair имеет отдельные фактические PID/log/status. После ROOT review он
запущен как PID368794; журнал
`/srv/marina/private/annotation-resource-recovery-v1/recovery.log`, состояние
`/srv/marina/private/annotation-resource-recovery-v1/status.json`. Исторический
`NONRETRYABLE_FAILURE` старого supervisor не заменяет этот новый recovery status.
Ограничение — 30 минут и прежний deadline $40; automatic retry/подмена исходного
времени VM отсутствуют. Шесть synthetic tests и actual metadata/recipe smoke PASS.
