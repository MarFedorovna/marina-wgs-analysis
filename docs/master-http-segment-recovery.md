# Восстановление только недостающих блоков master

`transfer/recover_master_http_segments.py` — локальный приёмник отдельных ZIP-сегментов, сохранённых через обычную разрешённую HTTPS-загрузку. Он не подключается по SSH, не читает ключи, не загружает базы и не выполняет аннотацию. Старый SFTP-план не изменяется. Новый явно проверенный профиль разрешает только это восстановление до 2026-11-06T07:10:00Z.

```powershell
python transfer/recover_master_http_segments.py --profile "$env:HTTP_RECOVERY_PROFILE" --profile-sha256 "$env:HTTP_RECOVERY_PROFILE_SHA256"
python transfer/recover_master_http_segments.py --profile "$env:HTTP_RECOVERY_PROFILE" --profile-sha256 "$env:HTTP_RECOVERY_PROFILE_SHA256" --execute --watch
```

Без `--execute` проверяются только локальные метаданные и наличие ZIP. Исполняемый профиль должен иметь schema_version=1, status=ROOT_REVIEWED_MASTER_HTTP_RECOVERY_V1, producer_sha256, deadline_utc, interval_seconds (2–30), old_pids, download_directory и отдельный state_directory. Поля original_plan, captured_ledger и http_manifest — точные малые pins `{path,bytes,sha256}`. Профиль и реальные контрольные суммы данных остаются приватными.

Приёмник использует прежнюю блокировку dataset.lock, байт 0, и требует завершения старых процессов. Повторно использованный живым процессом PID также блокирует запуск; программа ничего не останавливает. Проверяются точный inode/file_id, полный размер и физическая Windows FileStandardInfo-аллокация собственного `.part`. Исходные expected/source/chunk_bytes и первоначальное множество недостающих 64 MiB блоков неизменны.

Серверный закрытый HTTP manifest имеет status PASS_CLOSED_MISSING_MASTER_HTTP_SEGMENTS, source/expected/file_id/initial_missing_indices и список segments. Каждый элемент содержит name, размер ZIP bytes, index, suboffset, абсолютный offset, payload_bytes и payload_crc32. Имена — `chunk-NNNN-sub-00.zip` / `-01.zip`. Source size/mtime должны быть неизменны, zip_crc_verified=true; большой исходный master не хешируется заново.

Каждый ZIP_STORED содержит ровно `block.bin` и `segment.json`. JSON имеет schema_version=1, status=CLOSED_MISSING_MASTER_SEGMENT, source, expected, file_id, index, suboffset, offset, bytes, initial_missing_indices и payload_crc32. Размер сегмента — 32 MiB или остаток. Только два допустимых suboffset: 0 и 32 MiB. Источник, границы и CRC должны точно совпасть с проверенным manifest и исходным снимком журнала.

До записи проверяется CRC всех членов ZIP. Затем данные записываются только в первоначально недостающий диапазон. Завершённые старые блоки не переписываются. После fsync сохраняется собственная долговечная квитанция 32 MiB. Старый 64 MiB bitmap получает отметку лишь когда готовы все его сегменты. При повторном запуске уже подтверждённые сегменты не записываются снова; неполный неподтверждённый сегмент можно повторить в том же разрешённом диапазоне.

Исходный статус ошибки и исходный master ledger копируются в новую папку состояния. Только после завершения всего bitmap и проверки остальных закрытых результатов собственный `.part` переименовывается в final без перезаписи. Создаётся исходная SIZE ONLY completion note с expected/source, mtime и actual_sha256=null. После проверки всех файлов набора прежний status.json получает COMPLETED_CLOSED_DATASET_SIZE_ONLY_HTTP_RECOVERY и whole_dataset_complete=true с новой provenance профиля. Ожидаемый серверный SHA сохраняется; новый полный SHA master не заявляется.

Резерв во время записи — 5 GiB, поскольку весь собственный файл уже физически выделен. По новому сроку или ошибке свои ZIP, журналы и незавершённые данные сохраняются, автоматического перезапуска нет. Ничего не удаляется. Публичные Windows-тесты проверяют границы, источник, CRC, inode, запрет записи старого completed-блока и фиксацию 64 MiB только после всех 32 MiB сегментов.
