# Read-only монитор одной командой

На текущем сервере изолированный проверяемый пакет:

```bash
python3 /srv/marina/private/reviewed-monitor-v1/scripts/monitor_pipeline.py --once
```

Наблюдение раз в минуту, максимум восемь часов (Ctrl-C завершает только монитор):

```bash
python3 /srv/marina/private/reviewed-monitor-v1/scripts/monitor_pipeline.py --watch --interval 60 --max-hours 8
```

При желании добавить `--report-dir /srv/marina/private/reviewed-monitor-v1/monitor-reports`:
создаются только отдельные `status.json` и `status.html` в этом каталоге. Без этой
опции монитор ничего не записывает. Другие status streams не изменяются.
JSON на stdout обновляется каждую итерацию. Watch также завершается при исходном
бюджетном дедлайне из квитанции супервизора; это не останавливает VM или billing.

Монитор читает небольшие квитанции и `/proc`: фактические PID и совпадение назначения
процесса, номер попытки загрузчика, VEP bytes по подтверждённым SHA range receipts,
текущий MT recovery v2 и его отдельный итоговый audit, coverage, очередь аннотации.
Логический sparse-размер и выделенные блоки показываются отдельно и не считаются
скачанными байтами. Гигабайтные SHA повторно не вычисляются. Старые ошибки MT v1 и
старые строки логов не определяют текущее состояние v2. Нет download, restart,
kill, вызовов Docker/GATK или анализа пациентских вариантов.

Локальная безопасная проверка читает только разрешённую квитанцию переноса:

```powershell
python scripts/monitor_pipeline.py --local --root . --once
```

Приватная конфигурация SFTP/ключи, BAM/VCF не читаются. Отсутствующие или временно
недоступные квитанции не превращаются в PASS. AlphaMissense controls выделены
отдельно; они не заменяют ожидаемые VEP controls и этапы 6/7.

## X/Y остаётся отдельным решением пользователя

Показать план (без запуска calling):

```bash
python3 /srv/marina/private/reviewed-monitor-v1/scripts/sex_chromosome_plan.py --config /srv/marina/private/workflow-v1.json
```

Планировщик и его lock находятся в отдельном пакете монитора; в исходном checkout
команда `python3 scripts/sex_chromosome_plan.py --config /srv/marina/private/workflow-v1.json`.
У планировщика нет executor и `--execute`. Он не устанавливает клинический пол и
не проводит соматический анализ. У четырёх образцов плоидность unknown; PAR
alignment contract не подтверждён. Технические X/Y ratios не снимают эти блокировки.
Пользователь отдельно решает следующий шаг после просмотра [плана X/Y](sex-chromosome-plan.md).

Проверки: `python -m unittest discover -s tests -p test_monitor_pipeline.py`.
