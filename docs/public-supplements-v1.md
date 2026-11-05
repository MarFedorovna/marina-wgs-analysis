# Дополнительные публичные источники и полный охват требований

Новый пакет добавляет данные Orphadata Science и ClinGen CSpec отдельно от активного
VEP installer, MT и утверждённой очереди аннотации. Пациентские варианты не отправляются
в API. Ресурс не становится клиническим заключением и не включает автоматический ACMG.

На сервере показать фактические квитанции одной командой:

```bash
bash /srv/marina/private/reviewed-publicsupplements-v1/publicstages/01_prepare_public_supplements.sh
```

Подготовить только публичные источники (повторно использует готовые файлы):

```bash
bash /srv/marina/private/reviewed-publicsupplements-v1/publicstages/01_prepare_public_supplements.sh --execute
```

Установщик имеет отдельную блокировку, максимум четыре попытки на источник, максимум
30 минут на запуск, потолок источников 1 GB и исходный бюджетный дедлайн
2026-11-06T03:46:48.046025 UTC. Готовые источники проверяются по размеру и SHA256.
Частичные файлы сохраняются; resume требует неизменной HTTP identity и корректного
206 Content-Range. Если сервер не поддержал безопасный resume, частичный файл
сохраняется и требуется review; с нуля он молча не перезаписывается. HTTP 4xx,
кроме 408/429, требуют проверки запроса/доступа. У CSpec HEAD возвращал 400 при
работающем документированном GET, поэтому для HEAD 400/405 используется GET.

## Что действительно добавлено

Официальные страницы [Orphadata genes](https://sciences.orphadata.com/genes/),
[alignments](https://sciences.orphadata.com/alignments/) и
[phenotypes](https://sciences.orphadata.com/phenotypes/) перечисляют выпуск July 2026,
файлы от 23 June 2026 и CC BY 4.0. Точные ссылки взяты из HTML этих страниц,
не угаданы. Локальная квитанция фиксирует источник, release, транспортную identity,
время, размер, SHA256 и attribution. SHA здесь вычислен при получении; отдельный
publisher SHA256 не был предоставлен. OMIM ID в Orphadata не означает скачанный OMIM.

Фактически получены en_product6.xml (22,612,034 B), en_product1.xml (54,026,799 B),
en_product4.xml (47,893,357 B). XML root date соответствует 2026-06-23.

[CSpec Help](https://cspec.genome.network/cspec/ui/svi/help) документирует публичный
REST JSON. В текущем сервисе v1.3.14 `/SequenceVariantInterpretation` сообщает
247 сущностей; `/SequenceVariantInterpretation/id?detail=high` возвращает первую
страницу из 100. Страницы 2 и 3 запрошены с `pgSize=100`; adapter требует 247
уникальных IDs, равных entCount. Первая страница отдельно не считается полным реестром.
[ClinGen terms](https://www.clinicalgenome.org/docs/terms-of-use/) разрешают свободное
использование curated content под CC0; attribution/date сохраняются, linked-source
terms не отменяются. Registry snapshot — 2026-11-05, версии отдельных specs сохранены.

Реальная сборка `/srv/marina/private/publicsupplements-candidate-v1` содержит
8,500 gene–disease assertions, 52,456 cross-reference assertions и 116,640 уникальных
HPO assertions (24 дубликата по ключу учтены). SQLite: 34,263,040 B. Из всех 247 CSpec
сущностей отобраны **только current Released** с точным MONDO:0005045 и HGNC ID:
MYH7 v2.0; MYBPC3, TNNI3, TNNT2, TPM1, ACTC1, MYL2, MYL3 v1.0. Полные CriteriaCode
и RuleSet сохранены в восьми отдельных JSON для review. Это не утверждение, что у
всех остальных ГКМП-генов есть специфические правила.

## Offline adapter и граница интеграции

```bash
python3 /srv/marina/private/reviewed-publicsupplements-v1/scripts/public_supplement_adapter.py \
  --build-candidate /srv/marina/private/publicsupplements-candidate-v2
```

Команда читает только публичные источники и создаёт новый каталог; существующий
candidate не перезаписывает. `gene_evidence()` требует точные HGNC и MONDO ID.
В disease join допускается только источник MONDO с `E (Exact mapping...)` и
`Validated`; NTBT/BTNT и иные отношения сохранены как исходные утверждения.
Символьные aliases не угадываются. `annotate_variant_context()` — offline interface
для будущей интеграции: возвращает gene context/provenance, без variant pathogenicity
и без применения ACMG/PVS1. К пациентским данным этот adapter ещё не применялся.

Candidate имеет статус `CANDIDATE_PUBLIC_ADAPTER_NOT_BOUND`. ROOT review должен
проверить новые данные, код, controls и manifest перед отдельным изменением binding
этапов 6/7. Старое approval и sealed annotation runtime не менялись.

## Все требования и оставшиеся источники

[requirements-coverage-v1.json](../config/requirements-coverage-v1.json) содержит
все 69 источников, 61 строку матрицы и 66 source IDs: U01–19, M01–41, G42–44,
USER01–03. M26 — пустой неудавшийся конспект, требований из него не придумано.
Статический файл — inventory, без обещания завершённой аннотации. На сервере
добавляются только фактические небольшие квитанции:

```bash
python3 /srv/marina/private/reviewed-publicsupplements-v1/scripts/build_requirements_coverage.py \
  --server-root /srv/marina --output /srv/marina/private/publicsupplements-candidate-v1/requirements-coverage.json
```

Обязательны исследовательские stop-gain, stop-loss, start-loss, frameshift, splice,
позиции/HGVS/transcript.version, top10/top20 на образец и когорту, AlphaMissense с
точной идентичностью. LoF/GoF, NMD, механизм болезни и ACMG требуют evidence review;
consequence HIGH не становится PVS1. **KEGG enrichment остаётся withdrawn.**

[Predictor profiles](../config/predictor-installation-profiles-v1.json) сохраняют
official URLs, release/assembly, размеры, terms, access и adapter requirements для
REVEL, SpliceAI, FATHMM, CADD, VARITY, MetaLR, MaxEntScan, dbNSFP и LOFTEE. FATHMM
2.3 и FATHMM-XF — разные методы. VEP plugin release116 из audit — ссылка на
документацию, не установленный plugin для нашего VEP115. Отсутствующие оценки явно
MISSING, не benign/AF=0. dbNSFP gated institutional email и ND/noncommercial terms;
CADD 81G SNVs и большие supporting data не ставятся под видом small supplements.
Полные OMIM/HGMD/Mastermind, account-mediated SpliceAI и другие источники с доступом
требуют разрешённого локального экспорта; условия не принимаются автоматически.

Проверки: `python -m unittest discover -s tests -p 'test_public_supplement_adapter.py'`
и `test_requirements_coverage.py`. Проверяются exact IDs, Released state, полные
rules/provenance, duplicate refusal, broad≠exact, missing≠installed и withdrawn KEGG.
