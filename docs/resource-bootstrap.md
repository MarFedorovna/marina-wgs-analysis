# Одна серверная команда для ресурсов

Из корня проверенного проекта **на сервере**:

```bash
bash publicstages/00_prepare_all_resources.sh --execute
```

Без `--execute` команда только показывает план/состояние и не создаёт файлов или
загрузок. `--status` также только читает состояние. Пациентские файлы для этой
команды не нужны; она не запускает calling или аннотацию пациентов.

Если VEP уже скачивается под существующим супервизором, команда сообщает
`ALREADY_RUNNING_REUSING_READY_SMALLS_AND_ALPHA` и пути журналов. Она не останавливает
процесс, не создаёт второй VEP-download, не сбрасывает его части и не копирует кэш
на компьютер. Удерживаемая блокировка супервизора распознаётся и во время сетевой
паузы между повторами. Отдельная блокировка не допускает два одновременных запуска
единого bootstrap.

Автоматический набор: ClinVar и индекс, ClinGen gene validity/dosage, HGNC, HPO
gene/ontology, MANE summary, AlphaMissense canonical/isoforms, фиксированные VEP
115.2 Docker image и Ensembl 115 GRCh38 cache. Все существующие готовые файлы
переиспользуются; при подготовке/возобновлении сверяются версия, размер и SHA.
Частичный VEP-кэш докачивается через существующие проверенные 64 MiB range receipts;
готовность требует целого архива, publisher checksum и проверки распакованных файлов.

На чистом сервере обязательно передать исходное фактическое время создания VM:

```bash
bash publicstages/00_prepare_all_resources.sh --execute --created-utc ACTUAL_VM_CREATION_UTC
```

На текущем сервере timestamp читается из уже существующей private budget/supervisor
квитанции. Подмена этого времени запрещена. Лимит остаётся $40 при $2.39/час с
исходного создания VM, в рамках общего бюджета $50. Временные сетевые ошибки
повторяются с backoff 60–900 секунд до исходного абсолютного дедлайна; ошибки
контрольных сумм, версии, источника и доступа останавливают подготовку. Продление
лимита или гарантия завершения в бюджет не подразумеваются. Супервизор не прекращает
provider billing и не удаляет VM.

Эти $2.39/час и $40 относятся к текущей согласованной конфигурации. Для другой VM,
тарифа или бюджета сначала пересмотреть код/лимиты и фактический creation anchor;
команда откажется молча подменять существующий другой тариф или allowance.

Если небольшие ресурсы отсутствуют, их подготовка идёт в терминале перед запуском
отдельного серверного супервизора VEP. Для полностью отделённого холодного запуска
можно использовать одну строку с `nohup` и тем же фактическим timestamp:

```bash
nohup bash publicstages/00_prepare_all_resources.sh --execute --created-utc ACTUAL_VM_CREATION_UTC >>/srv/marina/logs/all-resources-bootstrap.log 2>&1 </dev/null &
```

## Что команда не обещает скачать

[resources-inventory.example.json](../config/resources-inventory.example.json)
разделяет источники на `auto`, `cache`, `manual`, `access_required`, `optional`.
Полный реестр требований также читается из
[evidence-sources.example.json](../config/evidence-sources.example.json).

SIFT, PolyPhen-2, gnomAD AF, rs/dbSNP и версии Ensembl проверяются по фактическому
`info.txt` после установки кэша. Full RefSeq/merged cache, полный dbSNP156, новый
gnomAD или все его AN/coverage/constraint не заявляются. MANE summary не является
полным RefSeq. [Последствия/SO](sources-and-consequences.md) рассчитываются VEP для
каждого транскрипта; это не отдельная база патогенности. LoF/GoF/NMD/PVS1 требуют
проверки механизма; LOFTEE пока optional/not installed, меток HC/LC нет.

OMIM, LOVD, GeneReviews, Orphanet, PubMed, VarSome и gene-specific ClinGen VCEP
остаются ручными источниками либо разрешёнными локальными экспортами. PMID — не
скачанный полный текст и не завершённый литературный обзор. dbNSFP4.9a из Analysis1
не переносится в новый анализ. FATHMM, REVEL, SpliceAI, CADD, Varity, MetaLR и
MaxEntScan остаются видимыми optional/not installed до выбора официального пакета,
проверки условий и реализации точного локального adapter. Эта команда не принимает
лицензии, не регистрируется в сервисах и не отправляет пациентские варианты наружу.

## Проверка состояния

```bash
bash publicstages/00_prepare_all_resources.sh --status
tail -n 25 /srv/marina/logs/annotation-databases.log
tail -n 10 /srv/marina/logs/annotation-supervisor.log
```

Повторный запуск безопасно переиспользует текущий процесс или готовые ресурсы.
Изменённый существующий рецепт/файл не перезаписывается. Завершение автоматического
набора не объявляет установленными все 69 источников и не заменяет экспертный
[допаудит интерпретации](medinter-gap-review.md).
