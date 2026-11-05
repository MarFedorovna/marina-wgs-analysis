# Одна команда для готовых небольших результатов

[Sync](../transfer/sync_completed_small_results.py) — отдельный native Windows
процесс для уже завершённых общих/индивидуальных таблиц, annotated VCF+TBI,
графиков, отчётов, текущих MTv2 results и gene coverage. Четыре BAM/BAI,
индивидуальные gVCF, raw reads и публичные базы исключены. Активный двух-BAM
writer, его `.part`, ключи, host pin, authentication и процессы не изменяются.

Используется существующий native Python с проверенным Paramiko5.1.0 и уже
существующая private SFTP configuration. [Пример профиля](../config/small-result-sync.example.json)
содержит только ссылки на этот config/plan и output paths; скопировать пример
следует в private ignored location и заполнить разрешённые пути. Секреты и
endpoint/key values не нужно переносить в профиль. Исполнитель читает существующий
auth config только в памяти и вызывает прежние `connect()`/`blocks()` с точным
pinned host key и RejectPolicy; ключ не копируется и trust не меняется.

Dry run — одна команда, без чтения auth config/ключа, сети, новых папок или импорта:

```powershell
& 'EXISTING_NATIVE_PYTHON.exe' 'transfer/sync_completed_small_results.py' --profile 'PRIVATE_SMALL_RESULT_PROFILE.json'
```

После ROOT review запуск/возобновление — та же команда с двумя флагами:

```powershell
& 'EXISTING_NATIVE_PYTHON.exe' 'transfer/sync_completed_small_results.py' --profile 'PRIVATE_SMALL_RESULT_PROFILE.json' --execute --watch
```

Профиль ссылается на live plan от reviewed `prepare_native_plan.py` v4. Этот plan
обновляет отдельный observer; sync сам не переаннотирует workflow и не объявляет
незавершённый sidecar готовым. Проверяются metadata SHA captured producer receipts
и связь каждого selected result с expected SHA/size. Индивидуальные annotated VCF
разрешены, gVCF исключены. Новые completed entries подхватываются каждый60секунд.
Отсутствующие пять обязательных sidecar roots остаются WAITING; временное отсутствие
pending downloads не означает окончательный полный backup.

Используется отдельный download subdirectory, собственный lock/status и отдельные
`.part`/`.source.json`. Продолжение требует прежние remote path/size/mtime и expected
producer SHA metadata; чужой partial/final не перезаписывается. SFTP читает
только недостающий contiguous suffix. Проверяется размер и remote size/mtime до
и после, затем происходит native atomic final rename и completion receipt.
Short chunks, изменившийся source, pause или deadline сохраняют partial без PASS.

Завершённые файлы передаются прежнему WSL `import_verified_download.py` с
`--verification size --completed-download --completion-evidence ssh_completed_response`.
Он использует общий короткий backup lock и существующий index/friendly folders.
Если importer занят, download сохраняется до следующей попытки; BAM writer не
останавливается. Другой native import watcher не видит отдельный вложенный
download directory. Для stable import используются private snapshots выбранных
plan entries и receipts; auth/key config в эти snapshots не включается.

Новые patient bytes **не хешируются**. `ActualSha256=null`; expected server SHA
остаётся provenance. Исторические SHA-verified записи не понижаются и не
объявляются заново проверенными. Уже корректно импортированные exact entries
пропускаются по существующей metadata/stat/index policy.

Профиль ограничивает один файл32GiB и все выбранные outputs64GiB, требует15GiB
свободного места. Превышение cap — явная остановка, а не тихий пропуск результата.
Предел можно отдельно рассмотреть до запуска. Wait/download≤8часов и исходный
absolute deadline2026-11-06T03:46:48.046025UTC сохраняются. Pause marker не удаляется
автоматически. Код не делает kill/restart, не меняет ports/ACL и не загружает базы.

Состояние находится в private `state_directory/status.json`; поля содержат counts,
received/expected bytes, missing completed roots и verification scope, без
endpoint/private-key/auth values. Ошибки транспорта не выводят сторонние
exception strings с конфигурацией. Это sync выбранных завершённых результатов,
не заявление о полной клинической интерпретации или SHA-verified full snapshot.

Тесты используют fake SFTP и synthetic receipt metadata: dry run без auth файла,
выбор четырёх sample VCF/TBI и исключение BAM/gVCF/resources, receipt tamper/caps,
resume ровно с оставшегося offset, changed remote mtime, short chunks,
pause/deadline до source access и explicit SIZE ONLY completion. Реальная
передача этим новым worker до ROOT review не запускается.
