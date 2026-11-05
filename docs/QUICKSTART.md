# Краткая инструкция

Все серверные команды выполняются из корня установленного кода. Данные проекта в примерах находятся в `/srv/marina`. Приватные планы запуска уже должны существовать: они связывают входы, версии источников, код и предел времени.

## Уже запущенный проект

Проверить состояние:

```bash
python3 scripts/monitor_pipeline.py --root /srv/marina --once --human
```

Наблюдать раз в минуту:

```bash
python3 scripts/monitor_pipeline.py --root /srv/marina --watch --interval 60 --human
```

Без `--human` выводится прежний JSON. Монитор читает небольшие квитанции и PID, не открывает BAM/VCF и не пересчитывает их SHA. Старая остановленная последовательная очередь показывается как историческая запись, если основной параллельный процесс уже работает. Завершение отмечается по квитанциям, без повторного аудита файлов.

`NOT_CREATED_OR_WAITING` означает, что итоговой квитанции ещё нет. Это не подтверждает запуск следующего контроллера. После готовности откройте `/srv/marina/outputs/cohort-sample-views-v1/index.html` или его сохранённую локальную копию.

## Рабочие скрипты и входы

Это один поток результатов. Контекстные таблицы присоединяются к общей аннотации; новый вызов вариантов для каждой таблицы не выполняется.

| Шаг | Скрипт запуска → исполнитель | Вход | Выход в `outputs/` |
|---|---|---|---|
| 1. Аннотация и отчёт | `parallel_annotation.py` → неизменённый `annotation_stages.py` | проверенный `05_filter`, GRCh38, VEP и локальные источники | `workflow-v1/stages/06_annotate/`, `07_report/` |
| 2. Ткань/белок | `continue_tissue_protein_context.py` → `export_tissue_protein_context.py` | завершённые 06/07 и подготовленные справочники | `tissue-protein-context-v1/` |
| 3. Популяции/ClinVar/MT | `continue_population_mt_reference_context.py` → `export_population_mt_reference_context.py` | завершённые 06/07, локальные источники, покрытие и MT | `population-reference-context-v1/` |
| 4. Единый master | `continue_combined_annotation_master.py` → `join_combined_annotation_master.py` | исходный master и завершённые шаги 2–3 | `combined-annotation-master-v1/` |
| 5. Общие и четыре пообразцовых набора | `continue_cohort_sample_views_v2.py` → `export_cohort_sample_views_v2.py` | 06/07 и единый обогащённый master | `cohort-sample-views-v1/` |
| 6. Дополнительный экспорт последствий | `run_consequence_export.sh` → `export_consequence_tables.py` | завершённые 06/07 | `consequence-tables-v1/` |
| 7. Локальная копия | `prepare_native_plan.py` + `native_sftp_download.py` | завершённые результаты на сервере | выбранная папка на компьютере |

Версия `v2` в имени экспортера — версия кода. Пользовательская выдача остаётся в одном каталоге `cohort-sample-views-v1/`.

## Как запускаются подготовленные этапы

Если процесс уже запущен, используйте монитор. Приведённые ниже команды нужны для подготовленного, ещё не запущенного этапа. Переменные указывают на существующие приватные файлы; их значения берутся из проверенного плана проекта и не сохраняются в GitHub. Используется именно проверенная копия кода, указанная в плане.

Основная аннотация с последующим отчётом:

```bash
python3 "$ANNOTATION_ENTRY" run --rule "$ANNOTATION_RULE" --reviewed-rule-sha256 "$ANNOTATION_RULE_SHA" --execute
```

`ANNOTATION_ENTRY` — путь к проверенной копии `parallel_annotation.py`. До пяти задач по хромосомам используют суммарно не более 20 CPU и 80 GiB. Кэш, референс, флаги VEP и источники остаются общими. Неполный результат прежнего процесса сохраняется, но не выдаётся за завершённый. Подробности: [параллельное выполнение](parallel-annotation.md).

Контекстные этапы и объединение запускаются по одному подготовленному правилу каждого шага:

```bash
python3 "$TISSUE_ENTRY" --rule "$TISSUE_RULE" --reviewed-rule-sha256 "$TISSUE_RULE_SHA" --execute
python3 "$POPULATION_ENTRY" run --rule "$POPULATION_RULE" --reviewed-rule-sha256 "$POPULATION_RULE_SHA" --execute
python3 "$MASTER_ENTRY" run --rule "$MASTER_RULE" --reviewed-rule-sha256 "$MASTER_RULE_SHA" --execute
python3 "$VIEWS_ENTRY" --rule "$VIEWS_RULE" --reviewed-rule-sha256 "$VIEWS_RULE_SHA" --execute
```

`*_ENTRY` соответствует скрипту запуска в таблице. Контроллеры ожидают нужные входы, не перезапускают основную аннотацию и не продлевают исходный срок работы. Без `--execute` эти команды показывают план или состояние. Не запускайте отдельные экспортеры параллельно с уже настроенными контроллерами.

## Файлы готовой выдачи

Пути ниже указаны относительно `cohort-sample-views-v1/`. Имена файлов образцов перечислены в `manifest.json`, поле `sample_file_labels`.

| Что нужно | Файл |
|---|---|
| Открыть всё | `index.html` |
| Полная общая таблица, все транскрипты | `masters/cohort.master.tsv.gz` |
| Общий топ-20 | `tables/cohort.top20.tsv` |
| Общая таблица последствий | `tables/cohort.functional.tsv.gz` |
| Категории и пути к отдельным файлам | `tables/functional-category-index.tsv` |
| Полная таблица образца | `masters/<метка>.master.tsv.gz` |
| Топ-20 образца | `tables/<метка>.top20.tsv` |
| Последствия образца | `tables/<метка>.functional.tsv.gz` |
| VCF образца | ссылки раздела образца в `index.html` |
| Фактические источники и контрольные суммы | `manifest.json` |

Недоступные оценки и источники отмечаются явно. Отсутствующий генотип не считается референсным. Категории считаются по уникальным аллелям с учётом транскриптных последствий; их суммы могут превышать число вариантов. MT остаётся отдельным приложением и не смешивается с ядерными строками master.

## Сохранить на компьютер

Планировщик `transfer/prepare_native_plan.py` выбирает завершённые результаты по квитанциям. Однократный просмотр нового плана из WSL, без изменения серверных ссылок и работающей загрузки:

```bash
python3 transfer/prepare_native_plan.py --host "$TRANSFER_HOST" --ssh-config "$TRANSFER_SSH_CONFIG" --known-hosts "$TRANSFER_KNOWN_HOSTS" --remote-root /srv/marina --workspace-root "$REMOTE_WORKSPACE" --jupyter-root "$REMOTE_JUPYTER_ROOT" --output "$NEW_NATIVE_PLAN" --compare --compare-with-plan "$ACTIVE_NATIVE_PLAN" --deadline-utc "$TRANSFER_DEADLINE_UTC"
```

Переменные относятся к существующему подключению и приватным файлам. `NEW_NATIVE_PLAN` — новый, ещё не существующий файл. Планировщик не скачивает файлы. [Подробности выбора готовых результатов](native-plan-completed-sidecars-v2.md).

Саму передачу выполняет `native_sftp_download.py` в Windows по **отдельной приватной конфигурации**. Если загрузка уже работает, второй экземпляр не запускают. Для подготовленного, ещё не запущенного профиля в PowerShell:

```powershell
python transfer/native_sftp_download.py --config "$env:NATIVE_SFTP_CONFIG" --execute
```

План метаданных не является конфигурацией SFTP. Новые таблицы и отчёты должны быть включены в проверенный список загрузки: работающий профиль BAM сам их не подхватывает. Ключи и параметры подключения не копируются в документацию или GitHub.

Этот способ сохраняет `.part` для продолжения и проверяет размер и неизменность источника. Завершение `COMPLETED_SIZE_ONLY` означает проверку размера, без нового SHA-256 содержимого. На компьютере открывайте `INDEX.html` соответствующей сохранённой выдачи; полная копия готова только после переноса всех выбранных файлов.

## Новый сервер

Начните с [установки](server-setup.md), [референса](reference-setup.md) и [подготовки источников](resource-bootstrap.md). Для нового анализа из BAM используются этапы `workflow.py call`, `joint`, `filter`, `annotate`, `report`; это отдельная подготовка, не команда возобновления текущего процесса. Завершённые этапы переиспользуются только при совпадении входов, кода и параметров. Предел времени скриптов не останавливает оплату облачной VM.

