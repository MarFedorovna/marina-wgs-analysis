# REVEL: фактическая загрузка исходного архива, scores выключены

Показать состояние:

```bash
python3 /srv/marina/private/reviewed-revel-source-v1/scripts/download_revel_source.py
```

Получить/безопасно продолжить только исходный архив:

```bash
python3 /srv/marina/private/reviewed-revel-source-v1/scripts/download_revel_source.py --execute
```

Точная ссылка взята из [официальной REVEL downloads page](https://sites.google.com/site/revelgenomics/downloads)
на [Zenodo record 7072866](https://zenodo.org/records/7072866). Это v1.3, файл
`revel-v1.3_all_chromosomes.zip`: **667,102,707 B**, publisher MD5
`3ea2bc33e6b5455fc7e9899da863b5fe`. Legacy landing page пишет 526 MB; для ограничения
используется фактический Content-Length и размер linked Zenodo record.

Загрузка началась на VM 2026-11-05 17:51 UTC, PID338632, отдельный каталог
`/srv/marina/resources/revel-source-only-v1`; status.json и resource-receipt.json
различают DOWNLOADING и DOWNLOADED. Лог:
`/srv/marina/private/reviewed-revel-source-v1/download.log`. Source-only ready требует
полного размера, publisher MD5, собственного SHA256 и ZIP format check. Частичные
файлы не считаются готовым архивом. Повторный запуск не создаёт второго writer.

Это **не установленный/проверенный predictor**: архив не распаковывается, scores
пациентов не вычисляются, plugin и очередь аннотации не меняются. Условия авторов —
free for non-commercial use, для другого использования нужен отдельный контакт.
Назначение проекта как non-commercial здесь не подтверждено. Перед включением нужны
license-purpose review, adapter под VEP115, GRCh38 REF/ALT и transcript/protein
matching, обработка отсутствующего lifted GRCh38 position и публичные controls.
Исходная таблица содержит GRCh37 и lifted GRCh38; GRCh37 score нельзя выдавать за
проверенный GRCh38 allele без проверки.

Потолок архива 700 MB, максимум 30 минут/четыре попытки и прежний абсолютный
дедлайн $40. Семь уже полученных Orphadata/CSpec файлов вместе — 149,232,651 B;
с REVEL итоговый ожидаемый объём **816,335,358 B**, меньше 1 GB. Программа не
принимает новые условия, не регистрируется и не обращается с пациентскими вариантами
в API. Другие большие predictors описаны отдельными installation profiles и не
выдаются за установленные или благоприятные при отсутствии scores.
