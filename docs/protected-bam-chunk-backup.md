# Возобновляемая передача BAM через существующий защищённый Jupyter

Этот fallback использует только уже авторизованный Jupyter file download. Новые
listeners, public endpoints, auth tokens, cookies, ключи, ACL и способы обхода
аутентификации не создаются. Существующий native SFTP writer не останавливается
и его `.part`/source metadata не изменяются. Базы и FASTQ не передаются.

[Server exporter](../transfer/prepare_bam_chunk_export.py) принимает private config
с completed BAM/receipt, точным ожидаемым размером и существующей protected export
directory внутри Jupyter root. Он сохраняет исторический checksum из завершённой
квитанции, но **не вычисляет новый hash BAM или chunks**. Export ID получается из
метаданных, не содержимого. Manifest доставляется на PC отдельно по уже pinned
managed SSH; он не берётся из непроверенного download payload.

Один chunk —64MiB; последний может быть короче. Полный manifest перечисляет
непересекающиеся offsets, sizes и уникальные deterministic basenames. Exporter
сверяет размер/mtime/inode источника до и после чтения, пишет resumable partial,
fsync и atomic no-replace final. Повторный запуск пропускает завершённые chunks
правильного размера. По умолчанию materialize только два новых chunks; все chunks
не создаются до review фактической скорости. После одобренного proof можно увеличить
`--max-chunks`; проверяется место для всего выбранного нового объёма плюс15GiB.
Runtime ограничен30min и прежним абсолютным VM deadline, с сохранением partial.

```bash
python3 transfer/prepare_bam_chunk_export.py --config /private/chunk-export.json
python3 transfer/prepare_bam_chunk_export.py --config /private/chunk-export.json --execute --max-chunks 2
```

Из manifest берут `jupyter_relative_path` и используют **обычный Download в уже
авторизованном Jupyter file browser**, сохраняя нормальный origin/referer.
Прямой переход к `/files` может вернуть403 unknown origin; такой путь не повторяют
через подмену headers/cookies и не обходят защиту. Не публиковать
patient manifest/links в GitHub. Для проверки скорости сначала скачивают один
готовый64MiB chunk; это не считается завершённым BAM. До дальнейшего экспорта
сравнивают реальную скорость и остаток времени VM.

[Local assembler](../transfer/assemble_bam_chunks.py) работает без сети и credentials.
Он проверяет manifest на gaps/overlaps/duplicates/path traversal, точные chunk
sizes, отсутствие browser partial siblings и стабильность файла. Используются
только canonical filenames; повторные browser copies с суффиксами `(1)` не
подставляются молча вместо ожидаемых chunks.

```powershell
python transfer/assemble_bam_chunks.py --manifest .local/trusted-chunk-manifest.json --download-directory "$env:USERPROFILE\Downloads" --destination "$env:USERPROFILE\Downloads"
python transfer/assemble_bam_chunks.py --manifest .local/trusted-chunk-manifest.json --download-directory "$env:USERPROFILE\Downloads" --destination "$env:USERPROFILE\Downloads" --execute
```

Assembly идёт в отдельный `.chunk-assembly.part`, с OS file lock и atomic journal
после каждого fsynced chunk. После сбоя допускается убрать только незакоммиченный
tail собственной assembly с идентичным manifest. PC должен сохранить15GiB свободно
с учётом всего ещё не записанного итогового BAM. Chunks автоматически не удаляются.
Необходимое дополнительное место: все chunks плюс собранный BAM; это проверяется
до полной передачи и перед сборкой.

При native SFTP `.part` того же итогового файла assembler отказывается работать.
Существующий итоговый файл никогда не заменяется. После всех chunks и exact full
size final публикуется atomic hard link без overwrite; записывается completion
receipt **COMPLETED_SIZE_ONLY**. Существующий native importer может обработать
этот canonical completed download обычным size-only режимом. До полного final
имени/размера/receipt нельзя сообщать, что BAM сохранён на PC.

Size-only подтверждает длину и непрерывную сборку manifest offsets, но не
криптографическую идентичность содержимого. Исходные historical checksums
сохраняются отдельно без утверждения, что новый local BAM ими перепроверен.

Проверка: семь synthetic tests — resumable exact assembly, gap/duplicate/traversal
rejection, native writer preservation, browser partial rejection, no overwrite,
recovery собственного незакоммиченного tail и15GiB guard. Synthetic fixture
содержит только искусственные байты; patient data tests не читают.
