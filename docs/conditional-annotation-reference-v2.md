# Условное продолжение после ROOT review

Новый `continue_annotation_reference_v2.py` по умолчанию только готовит приватное
правило для review. Он не утверждает binding и не запускает очередь без явного
`--execute --reviewed-rule-sha256 EXACT_ROOT_REVIEWED_SHA`.

Rule фиксирует controller/materializer/dependency SHA, исходный и новый configs,
draft binding, overlay, public-control receipt, все четыре adapter files, исходный
15-file code package и external upstream sources, принятый stage5 receipt и прежний
deadline $40 / $2.39/h. Допустимы только уже проверенные 39-gene public panel adapters:
gene evidence, официальный CADD1.7 и raw API/REVEL evidence. Новые методические
предложения этим правилом не утверждаются.

Контроллер ждёт исходный VEP core PASS. Затем pinned materializer создаёт новый
sidecar manifest рядом с исходным (это необходимо неизменённому cache/ClinVar mount
контракту). Проверяется равенство `deepcopy(core) + EXACT reviewed resource_additions`,
неизменность исходного manifest и config delta только `annotation.local_adapters` /
`annotation.resource_manifest`. У final binding разрешено ровно одно отличие от
draft: SHA нового sidecar добавлен в существующий generic path→digest guard
`external_upstream_code`. Неизменённый sealed `verify_bindings` реально проверяет
этот путь; никаких новых игнорируемых config fields нет.

После этих проверок создаётся approval **точного final binding SHA**, с SHA самого
ROOT-reviewed deterministic rule. Только затем запускается прежний sealed
`annotation_queue.py`: 20 публичных VEP controls → stage6 → stage7. Никакие stages
3–5, VM/database jobs или calling не перезапускаются. Другая annotation queue или
workflow child, неизвестный approval, изменённый config/source/sidecar/binding
приводят к отказу. Контроллер не убивает чужие процессы и не ослабляет locks.

ROOT сначала отдельно проверяет, что старый controller действительно WAITING и не
имеет child, и останавливает только его; новый controller этого не делает.
Создание/выполнение правила разрешает ROOT после review кода и тестов, без отдельного
вопроса пользователю о внутренних pipeline approvals.

На VM подготовить правило:

```bash
python3 /srv/marina/private/reviewed-annotation-reference-controller-v1/scripts/continue_annotation_reference_v2.py
```

После ROOT review exact rule SHA и остановки только старого WAITING controller:

```bash
python3 /srv/marina/private/reviewed-annotation-reference-controller-v1/scripts/continue_annotation_reference_v2.py \
  --execute --reviewed-rule-sha256 EXACT_ROOT_REVIEWED_RULE_SHA
```

Показать только состояние:

```bash
python3 /srv/marina/private/reviewed-annotation-reference-controller-v1/scripts/continue_annotation_reference_v2.py --status
```

Мониторинг контроллера ограничен максимум восемью часами и прежним VM deadline.
Запущенная sealed queue дополнительно имеет собственный исходный budget guard.
Автоматического restart при ошибке нет; неизвестное частично выполненное действие
требует ROOT review. Consequence-export после stage7 пока остаётся отдельным
проверенным wrapper с gates обоих PASS receipts; он не добавлен в это правило.

Семь synthetic tests проверяют exact final approval, запрет неожиданного config /
source / sidecar / binding change, отказ при другом annotation job, невозможность
перезаписать чужой approval и process identity. Они не используют пациентов или сеть.
