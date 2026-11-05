# MedInter: независимые предложения дополнений v2

Проверка 2026-11-05. **Только предложения для решения пользователя. Ничего не скачано, не установлено, не зарегистрировано; pipeline не изменён; пациентские варианты внешним API не передавались.** Некоммерческая цель уже подтверждена и повторного подтверждения не требует. Приоритет MUST означает важность для обоснованной интерпретации, а не разрешение установки и не блокировку текущих исследовательских таблиц.

Проверены 44 SHA256 исходных MedInter/GATK по локальному inventory: 41 уникальный MedInter (20 PDF, 21 HTML) и 3 GATK. Просмотрены обзор/заголовки и тематический контекст всех извлечений, заново извлечены приоритетные страницы оригинальных PDF. Это не заявление о повторном визуальном просмотре каждого слайда. M26 содержит сообщение о неудаче конспекта; требований из него нет. Применимых AGENTS.md в workspace и проверенных родительских каталогах не найдено. Метод PDF skill применён для read-only чтения.

Старый medinter-gap-review.md и requirements-coverage-v1.json не являются актуальным свидетельством установки: ряд их MISSING/«не установлен» уже устарел. Фактические состояния ниже получены из текущего handoff и локальных configs/docs.

## Уже доступно или находится в работе — повторно не предлагать

ClinVar GRCh38 2026-10-25; ClinGen validity/dosage 2026-11-05; HPO 2026-09-01; HGNC 2026-10-30; MANE1.4; GENCODE49 и CDS mosdepth QC 4/4; VEP115.2 установлен, cache115 ещё загружается (84.7% по handoff 23:27), stage6 ожидает cache. AlphaMissense canonical/isoforms READY и публичные controls прошли; это ещё не controls реального VEP-run. Orphadata July2026 все три источника и bridge готовы. CSpec registry247/247 и восемь Released HCM gene-specific rules уже скачаны. GeneReviews официальный полный NLM архив готов.

SpliceAI1.3.1 с пятью weights и isolated runtime установлен, synthetic control прошёл; default GENCODE24 не равен нашему49, patient run не выполнялся. REVEL1.3 полный исходник и public HCM subset получены, transcript mapping проверяется. CADD1.7 public HCM39 CDS subset готов; полный 87GB файл не скачан. MyVariant public39 prefetch выполнен (metadata build20250624/dbNSFP4.8a/gnomAD2.1.1/ClinVar2025-05), raw arrays не подменяют новые core snapshots. Новый combined candidate binding ожидает root review. RefSeq cache115 только queued, не установлен. OMIM mim2gene готов; full OMIM/VarSome требуют отдельного доступа, LOVD старое MYH7 XML не принято как GRCh38 mapping.

Новая autosomal calling/filtering цепочка выполнена; все четыре MT технические проверки PASS_RESEARCH_WITH_LIMITATIONS, но NUMT boundaries и clinical LoD не валидированы. Один haplochecker flag не доказал контаминацию. X/Y без пола остаётся техническим QC; CNV/SV/repeats не выполнены.

## Короткий список решений

| ID | Приоритет | Что предлагается | Рекомендуемое решение |
|---|---|---|---|
| P01 | MUST | gnomAD: AC/AN, homozygotes, FILTER и покрытие локуса | Одобрить ограниченный план/методику, детали ниже |
| P02 | MUST | MITOMAP + HelixMTdb + gnomAD MT: специализированные частоты и публикации | Одобрить ограниченный план/методику, детали ниже |
| P03 | MUST | MT-specific manual ACMG worksheet и контекст ткани/гетероплазмии | Одобрить ограниченный план/методику, детали ниже |
| P04 | USEFUL | Значимость сердечного транскрипта/экзона: GTEx и pext | Одобрить ограниченный план/методику, детали ниже |
| P05 | MUST | NMD и новый стоп после frameshift: воспроизводимый transcript review | Одобрить ограниченный план/методику, детали ниже |
| P06 | MUST | Независимая локальная проверка HGVS: Mutalyzer или VariantValidator | Одобрить ограниченный план/методику, детали ниже |
| P07 | USEFUL | UniProt: домены, функциональные участки и альтернативные изоформы | Одобрить ограниченный план/методику, детали ниже |
| P08 | MUST | ClinVar SCV + первичная литература: локальное досье короткого списка | Одобрить ограниченный план/методику, детали ниже |
| P09 | MUST | Контроль сопоставления predictor scores и актуальных PP3/BP4 | Одобрить ограниченный план/методику, детали ниже |
| P10 | MUST | Локальная IGV-проверка приоритетных SNV/indel и MT | Одобрить ограниченный план/методику, детали ниже |
| P11 | MUST | Структурированные phenotype/pedigree inputs вместо ещё одной базы | Одобрить ограниченный план/методику, детали ниже |
| P12 | SEPARATE_BRANCH | CNV/SV и другие пропущенные классы: отдельный обоснованный этап | Отложить отдельным этапом |

## P01. gnomAD: AC/AN, homozygotes, FILTER и покрытие локуса (MUST)

Источник M10: `3.3 ACMG без формализма_ как на практике классифицировать варианты_watermark (1).pdf`, PDF с. 9–10. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M33: `Клинические базы данных и ресурсы для интерпретации генетических вариантов (1).html`, HTML [38:40], gnomAD. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** Есть AF из VEP-кэша; MyVariant содержит исторический gnomAD 2.1.1. Современных полных AC/AN/nhomalt/quality/ancestry и all-sites AN для исследуемых областей нет.
**Для нашего случая:** Проверка редкости приоритетных ядерных HCM SNV/indel; отсутствие записи не равно AC=0 при хорошем AN.
**Добавочная ценность:** Деноминатор и callable context позволяют отличить редкий вариант от плохо представленного участка. Частота когорты из четырёх образцов остаётся отдельной.

**Доступ:** Публичные release downloads. Предложен региональный публичный HCM39 extract, а не запросы по вариантам пациента или полный dump. Browser API schema нестабилен.
**Условия:** Историческая официальная gnomAD документация указывает ODbL; условия конкретного актуального файла сохранить до загрузки/распространения.
**Цена/аккаунт:** Публичный доступ; фиксированную цену egress не проверяли. Платные облачные ресурсы не согласованы. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Размер выбранного регионального extract и сетевого чтения неизвестен; сначала size/egress plan. Полные геномные файлы не предлагать в текущий бюджет.

**Решение:** Одобрить малый региональный план: release-specific AC/AN/quality + coverage/all-sites AN. Constraint добавлять только отдельным файлом v4.1.2: этот релиз исправляет rank/bin, не AF и не LOEUF values.

Первичные страницы: [1](https://gnomad.broadinstitute.org/downloads), [2](https://broadinstitute.github.io/gnomad_methods/_modules/gnomad/resources/grch38/gnomad.html), [3](https://discuss.gnomad.broadinstitute.org/t/stratified-coverage-counts-for-gnomad-v4/814/2), [4](https://gnomad.broadinstitute.org/news/2026-09-gnomad-v4-1-2/).

## P02. MITOMAP + HelixMTdb + gnomAD MT: специализированные частоты и публикации (MUST)

Источник M14: `4.2 Митохондриальная ДНК_ особенности анализа и интерпретации_watermark (2).pdf`, PDF с. 8, 12. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M35: `Митохондриальная ДНК особенности анализа и интерпретации (1).html`, HTML [18:15], таблица частот. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** MT calling/QC уже выполнен, но это не локальная база MITOMAP/HelixMTdb. Наличие nuclear gnomAD AF и HCM predictor subsets не закрывает MT.
**Для нашего случая:** MT-TI и прочие MT-кандидаты: rCRS NC_012920.1, homoplasmic/heteroplasmic counts, haplogroup context, published disease evidence.
**Добавочная ценность:** Добавляет mtDNA-specific знаменатели и литературу; патогенные утверждения MITOMAP отделять от популяционных записей. Helix и GenBank не считать доказанной здоровой контрольной когортой.

**Доступ:** Helix официально даёт публичную Download HelixMTdb ссылку. MITOMAP вернул 403 при read-only проверке, способ разрешённого bulk/download не подтверждён. gnomAD MT formats описаны авторами.
**Условия:** Helix просит цитирование и ссылку; полной лицензии скачиваемого файла пока не проверено. MITOMAP условия неизвестны: не объявлять CC0/unrestricted.
**Цена/аккаунт:** На Helix download странице платёж не указан; точная цена/стоимость передачи неизвестна. MITOMAP неизвестно. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Размеры трёх выбранных файлов неизвестны, без предварительного численного обещания.

**Решение:** Одобрить отдельный небольшой MT-reference план после проверки metadata/terms. Непонятный TOOLST из HTML не переименовывать в MitoTip и не устанавливать по догадке.

Первичные страницы: [1](https://www.mitomap.org/MITOMAP/Resources), [2](https://www.helix.com/mitochondrial-variant-database), [3](https://pmc.ncbi.nlm.nih.gov/articles/PMC8896463/).

## P03. MT-specific manual ACMG worksheet и контекст ткани/гетероплазмии (MUST)

Источник M14: `4.2 Митохондриальная ДНК_ особенности анализа и интерпретации_watermark (2).pdf`, PDF с. 8–11. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** PASS_RESEARCH_WITH_LIMITATIONS относится к техническому MT-run; manual nuclear ACMG worksheet не заменяет отдельные MT rules. Registry CSpec уже скачан: сначала проверить, есть ли нужный документ внутри, не дублировать загрузку.
**Для нашего случая:** Раздельная оценка mtDNA protein-coding/tRNA/rRNA, гаплогруппы, ткани, maternal history, гетероплазмии и NUMT flags.
**Добавочная ценность:** MT-TI не является missense protein: AlphaMissense/REVEL и ядерная NMD/splicing схема к нему неприменимы. Учебные 1.5% нельзя выдавать за подтверждённый LoD нашего workflow.

**Доступ:** Официальная статья и ClinGen affiliation/guidance публично читаемы; сначала использовать уже локальный registry, затем только недостающее разрешённое руководство.
**Условия:** Условия статьи и linked guidance сохранять отдельно; публичная доступность не означает права на весь контент MITOMAP.
**Цена/аккаунт:** Для ручной методики платный аккаунт не требуется; экспертиза не оценивалась как услуга. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Новых больших данных не нужно; размер недостающего guidance неизвестен.

**Решение:** Одобрить локальную MT evidence worksheet; вся классификация остаётся ручной. Один haplochecker flag не объявлять доказанной контаминацией.

Первичные страницы: [1](https://pmc.ncbi.nlm.nih.gov/articles/PMC7717623/), [2](https://clinicalgenome.org/affiliation/50027/).

## P04. Значимость сердечного транскрипта/экзона: GTEx и pext (USEFUL)

Источник M11: `3.4 Критерий PVS1_ когда loss-of-function действительно значим_watermark (1).pdf`, PDF с. 4, 9, 13. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** MANE1.4, GENCODE49 и все VEP transcripts уже есть/готовятся; тканевая экспрессия и включение экзона в локальном bundle отсутствуют.
**Для нашего случая:** Ручная проверка biologically relevant exon/transcript для MYBPC3 и других HCM pLoF/splice кандидатов.
**Добавочная ценность:** MANE/canonical сами по себе не доказывают экспрессию нужного экзона в сердце; pext является поддерживающим контекстом, не самостоятельным PVS1.

**Доступ:** GTEx open-access derived expression, gnomAD public gene model/pext downloads. Raw RNA-seq и donor-level данные защищены и для цели не нужны.
**Условия:** GTEx разрешает использование с acknowledgement; лицензию выбранного производного файла и условия gnomAD сохранить. Никакого доступа к protected donor data.
**Цена/аккаунт:** Open-access данные; отдельные крупные GTEx файлы бывают requester-pays. Нужны только малые summary, платный аккаунт не предполагается. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Размер summary/subset неизвестен; не загружать RNA-seq/QTL all-pairs.

**Решение:** Одобрить компактный HCM transcript/exon summary после shortlist. Проверять accession/version: GTEx/gnomAD gene models могут отличаться от GENCODE49; новый GTEx V11 использует GENCODE47.

Первичные страницы: [1](https://www.gtexportal.org/home/downloads/adult-gtex/overview), [2](https://www.gtexportal.org/home/faq), [3](https://gnomad.broadinstitute.org/news/2025-01-update-gnomad-browser-gene-model-hail-tables/), [4](https://gnomad.broadinstitute.org/news/2026-03-gnomad-v4-1-1/).

## P05. NMD и новый стоп после frameshift: воспроизводимый transcript review (MUST)

Источник M11: `3.4 Критерий PVS1_ когда loss-of-function действительно значим_watermark (1).pdf`, PDF с. 9. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M11: `3.4 Критерий PVS1_ когда loss-of-function действительно значим_watermark (1).pdf`, PDF с. 4. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** Есть coordinates, exon/CDS/protein lengths и pLoF flags; NMD остаётся unassessed. Текущий GENCODE49 уже достаточен как структурная основа; покупать новую базу не нужно.
**Для нашего случая:** Короткий список stop-gain/frameshift/start-lost/splice; только при disease mechanism и текущей gene-specific спецификации, допускающих PVS1.
**Добавочная ценность:** Выяснить положение фактического нового стопа, последнего exon junction, rescue isoforms/альтернативного старта и долю утраченного белка. Положение исходного frameshift не равно положению нового стопа.

**Доступ:** Локальная методика на существующих versioned FASTA/GTF/VEP и скачанных CSpec. LOFTEE не назван как требование в лекциях и не необходим для этой проверки.
**Условия:** Использовать уже зафиксированные лицензии локальных ресурсов; SVI guidance — отдельное цитирование.
**Цена/аккаунт:** Новая база/ключ не нужны; вычислительное время не измерено. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Новый dataset не нужен; объём таблицы зависит от shortlist.

**Решение:** Одобрить методический review без изменения текущего calling. Для MYH7 применимость PVS1 проверять по уже скачанному MYH7 v2.0, а не по generic HIGH.

Первичные страницы: [1](https://www.clinicalgenome.org/docs/variant-curation-sop/).

## P06. Независимая локальная проверка HGVS: Mutalyzer или VariantValidator (MUST)

Источник M12: `3.5 Как правильно записывать варианты_ номенклатура без ошибок и путаницы_watermark (1).pdf`, PDF с. 20–21. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M12: `3.5 Как правильно записывать варианты_ номенклатура без ошибок и путаницы_watermark (1).pdf`, PDF с. 20. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** VEP выдаёт HGVS, но независимый HGVS validator с подготовленными референсами не установлен. RefSeq cache115 только queued и не равен Mutalyzer reference cache.
**Для нашего случая:** Проверить shortlist indel/MNV и версии c./p./g.; различать VCF left-alignment и HGVS 3-prime rule.
**Добавочная ценность:** Независимая проверка nomenclature/reference mapping перед рукописью; не дублирует аннотацию всего VCF.

**Доступ:** Mutalyzer 3.1.1 официально описывает локальный CLI и file-based reference cache. Нужно заранее загрузить публичные references и проверять работоспособность с отключённой сетью: cache miss способен вызвать NCBI retrieval.
**Условия:** Mutalyzer source MIT по официальному репозиторию; права reference sequences отдельно. VariantValidator здесь альтернатива, его package/dependencies/terms ещё не проверены.
**Цена/аккаунт:** Mutalyzer CLI open source; reference egress/installation cost не измеряли. NCBI key в docs optional, не обязательная покупка. После одобрения: ключ не обязателен; для настройки retrieval требуется служебный email по документации. Его не запрашиваем и установку сейчас не выполняем.
**Размер:** Размер выбранного transcript/reference cache неизвестен; full chromosome cache заранее не обещается.

**Решение:** Одобрить локальный validator proof-of-concept на публичных controls и ограниченных reference IDs; не отправлять варианты пациента в hosted сервис.

Первичные страницы: [1](https://github.com/mutalyzer/mutalyzer), [2](https://mutalyzer.readthedocs.io/en/latest/), [3](https://mutalyzer.readthedocs.io/en/latest/usage.html).

## P07. UniProt: домены, функциональные участки и альтернативные изоформы (USEFUL)

Источник M09: `3.2 Главные базы данных для интерпретации вариантов_ что использовать и когда.pptx_watermark (1) (1).pdf`, PDF с. 3. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M07: `2.3 Референсный геном и базы данных_ как не утонуть в источниках и использовать их осмысленно.pptx_watermark (1).pdf`, PDF с. 20–21. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** Protein consequence/length уже есть. Локальный versioned UniProt feature/context extract для HCM белков не подтверждён; VEP protein ID не означает наличие curated domain evidence.
**Для нашего случая:** Оценка утраченного C-конца, функционального домена и механизма shortlisted variants.
**Добавочная ценность:** Добавляет curated feature coordinates/evidence и публикации. Само попадание в домен не назначает PM1; приоритет имеют уже скачанные HCM CSpec определения.

**Доступ:** Официальный protein resource; предложен public-gene/protein extract, без variant queries. API/download metadata и точные release records предстоит подготовить после одобрения.
**Условия:** Текущая help/license страница дала JS fallback; обнаружено расхождение с SIB manual. Лицензия конкретного файла требует проверки, CC-BY автоматически не приписывать.
**Цена/аккаунт:** Свободный доступ к просмотру; bulk/API price не подтверждён. Платная лицензия не согласована. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Размер выбранных HCM feature records неизвестен; full UniProt download не нужен.

**Решение:** Одобрить только небольшой reference-only extract после подтверждения terms; если CSpec уже содержит нужные regions, дополнительную базу отложить.

Первичные страницы: [1](https://www.uniprot.org/help/license), [2](https://www.uniprot.org/help/programmatic_access), [3](https://web.expasy.org/docs/userman.html).

## P08. ClinVar SCV + первичная литература: локальное досье короткого списка (MUST)

Источник M10: `3.3 ACMG без формализма_ как на практике классифицировать варианты_watermark (1).pdf`, PDF с. 11, 17, 19. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M33: `Клинические базы данных и ресурсы для интерпретации генетических вариантов (1).html`, HTML [41:20], VarSome/RAB27A, Mastermind links. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** Свежий ClinVar GRCh38 snapshot уже есть. Полный SCV evidence dossier и прочитанные variant-specific публикации для shortlist не подтверждены. GeneReviews уже скачан, его повторная покупка/загрузка не нужна.
**Для нашего случая:** Разобрать conflicts, disease, submitter, dates, SCV versions, assay validity, PS1/PM5 comparisons и независимость наблюдений.
**Добавочная ценность:** Сводный VCF не содержит весь объём SCV ObservedIn/assertion method/citations. Не считать ClinVar classification независимым PP5 и не дублировать одну публикацию несколько раз.

**Доступ:** NCBI публичные VCV/RCV XML содержат SCV и ObservedIn; сначала использовать доступные local data, затем разрешённые public gene-level citation/submission файлы или ручную очередь. Mastermind только optional shortcut, не предпосылка.
**Условия:** NCBI и права отдельных публикаций/submitter content сохранять. Статья/abstract не автоматически unrestricted full text. Условия Mastermind API/full export не подтверждены.
**Цена/аккаунт:** NCBI FTP без login; конкретные full-text статьи могут требовать доступа. Платный Mastermind/VarSome не одобрять по умолчанию; цены не проверены. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Размер SCV subset неизвестен. Полный ClinVar XML для первого shortlist не планировать без отдельной size/streaming оценки.

**Решение:** Одобрить ручное evidence dossier и компактные SCV/citation supplements; запросы пациентских вариантов внешним API остаются запрещены в рамках этого предложения.

Первичные страницы: [1](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/), [2](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/), [3](https://www.ncbi.nlm.nih.gov/clinvar/docs/submission_instructions/), [4](https://www.genomenon.com/mastermind).

## P09. Контроль сопоставления predictor scores и актуальных PP3/BP4 (MUST)

Источник M11: `3.4 Критерий PVS1_ когда loss-of-function действительно значим_watermark (1).pdf`, PDF с. 13. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M10: `3.3 ACMG без формализма_ как на практике классифицировать варианты_watermark (1).pdf`, PDF с. 18. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** SpliceAI, AlphaMissense, REVEL raw, CADD1.7 и MyVariant raw уже подготовлены. Остались применимость transcript/version/AA mapping и публичные controls реального VEP; новые массовые predictor пакеты не требуются.
**Для нашего случая:** Аннотация missense/splice HCM-кандидатов с exact allele, transcript identity и source-specific missing reasons.
**Добавочная ценность:** SpliceAI default GENCODE24 отличается от GENCODE49; REVEL ENST без версии; FATHMM legacy/XF/MKL и неодинаковые raw arrays не сводить к одному scalar. HCM39 CDS subset не содержит все noncoding/indel/MT оценки.

**Доступ:** Уже локальные runtime/weights/subsets; это предложение проверок, не повторной установки. MaxEntScan optional только при отдельной механистической необходимости.
**Условия:** Пользователь уже подтвердил noncommercial research; повторное подтверждение NC-цели не нужно. Условия каждого score источника сохраняются.
**Цена/аккаунт:** Новый платный API/full dbNSFP не нужен; контрольный runtime cost ещё не измерен. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Новая база не нужна; объём дополнительных controls неизвестен.

**Решение:** Одобрить validation gate после готовности VEP; predictor vote оставлять описательным, PP3/BP4 назначать только ручным применимым правилом.

Первичные страницы: [1](https://github.com/Illumina/SpliceAI), [2](https://sites.google.com/site/revelgenomics/downloads), [3](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/).

## P10. Локальная IGV-проверка приоритетных SNV/indel и MT (MUST)

Источник M06: `2.2 Что происходит после секвенирования_ биоинформатический анализ (пайплайн), метрики качества и что такое IGV.pptx_watermark (1).pdf`, PDF с. 10. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M06: `2.2 Что происходит после секвенирования_ биоинформатический анализ (пайплайн), метрики качества и что такое IGV.pptx_watermark (1).pdf`, PDF с. 38. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M19: `5.4 Верификация результатов_ когда и как перепроверять варианты_watermark (1).pdf`, PDF с. 3. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** BAM/BAI и QC есть, GENCODE49 mosdepth CDS QC прошёл. Это не доказательство ручного read-review каждого приоритетного варианта; наличие установленного IGV не подтверждалось.
**Для нашего случая:** Проверка strand/read-end/mapQ/baseQ/softclips/nearby indel/paralog context, MT NUMT flags и спорных allele balances.
**Добавочная ценность:** Порог DP/GQ/FT не заменяет визуальную поддержку. Итог review: reviewed_pass/concern/unassessed с причиной, а не автоматический clinical PASS.

**Доступ:** Официальный IGV Desktop локально открывает BAM/VCF/reference. Проверить имеющуюся установку; только если её нет — предложить portable desktop package. Внешний upload не нужен.
**Условия:** IGV MIT; Java/runtime условия отдельно. Release на официальной download странице 2.19.8, Java21+ для текущей 2.x линии.
**Цена/аккаунт:** Open source, платный аккаунт не нужен; размер package/runtime и время review не измерены. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Размер установки неизвестен; BAM уже существуют, повторный calling не требуется.

**Решение:** Одобрить offline review shortlist с журналом причин; identity/relatedness/contamination требуют отдельных проверок и не выводятся из изображения одного локуса.

Первичные страницы: [1](https://igv.org/doc/desktop/), [2](https://igv.org/doc/desktop/DownloadPage/index.html).

## P11. Структурированные phenotype/pedigree inputs вместо ещё одной базы (MUST)

Источник M08: `3.1 Как читать клинические данные пациента для генетического анализа_watermark (1).pdf`, PDF с. 5. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M10: `3.3 ACMG без формализма_ как на практике классифицировать варианты_watermark (1).pdf`, PDF с. 12. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M32: `Как читать клинические данные пациента для генетического анализа (1).html`, HTML: SAMS, Phenopacket. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** HPO, Orphadata, GeneReviews и ClinGen уже доступны. Не хватает индивидуальной клиники/возраста/пола/материала/семейных связей и критериев включения; database download этого не исправит.
**Для нашего случая:** HCM phenotype matching, differential context, inheritance/phase/de novo/segregation и tissue-specific MT interpretation.
**Добавочная ценность:** Локальная таблица HPO present/absent/unknown + onset + source, pedigree relationship, cardiac examinations, maternal tissue history. Четыре образца не считать автоматически семьёй или независимыми probands.

**Доступ:** Локальная форма; Phenopacket — необязательный обменный формат, SAMS upload не требуется. Сведения поступают только от пользователя/уполномоченного специалиста.
**Условия:** Новая data license не нужна; персональные сведения остаются локально и не попадают в proposal JSON/GitHub.
**Цена/аккаунт:** Новый аккаунт/платная база не нужны; клиническую консультацию не заказывали. После отдельного одобрения понадобятся только реально имеющиеся clinical/pedigree сведения, без регистрации и ключей.
**Размер:** Новая database не нужна; локальный шаблон ещё не создаётся этим review.

**Решение:** Одобрить перечень недостающих полей; сейчас не запрашивать документы по одному. До получения данных phenotype/de novo/segregation остаются unknown.

Первичные страницы: [1](https://phenopacket-schema.readthedocs.io/en/latest/).

## P12. CNV/SV и другие пропущенные классы: отдельный обоснованный этап (SEPARATE_BRANCH)

Источник M19: `5.4 Верификация результатов_ когда и как перепроверять варианты_watermark (1).pdf`, PDF с. 4. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.
Источник M06: `2.2 Что происходит после секвенирования_ биоинформатический анализ (пайплайн), метрики качества и что такое IGV.pptx_watermark (1).pdf`, PDF с. 35. Полный SHA — в JSON; локальные ссылки сохранены отдельно в private `.local/medinter-review-links-v2.json`.

**Пробел:** CNV/SV/repeat calling не выполнен; ClinGen dosage уже скачан, но это interpretation data, а не caller. Обычный nuclear SNV/indel и MT small-variant run не закрывают эту ветку.
**Для нашего случая:** Только если shortlist/фенотип/семейная картина оправдывает поиск CNV/SV или второго аллеля, либо нужно полноценное исследование иных классов.
**Добавочная ценность:** Закрывает иной класс событий; новый caller/QC/validation design нужен до выбора annotation databases. mtDNA deletion/mtCN также не доказываются готовым small-variant MT-run.

**Доступ:** Сначала отдельный design: BAM support, callable regions, benchmark, confirmation и budget. DGV/DECIPHER/gnomAD SV/ISCN условны, не обязательная пачка установок.
**Условия:** Существующие ClinGen dosage terms уже сохранены; будущие DGV/DECIPHER/ISCN access/licenses проверять отдельно. Конкретная лицензия/цена не утверждается.
**Цена/аккаунт:** Стоимость caller/compute/validation неизвестна; исходный бюджет $50 и VM deadline не позволяют автоматически добавить эту ветку. Для подготовки предложения ничего; реализация только после одобрения.
**Размер:** Неизвестно до выбора caller и scope; большие SV/repeat datasets не загружать.

**Решение:** Отложить за пределы текущего SNV/indel+MT report, явно сохранить limitation; запускать только отдельным согласованным планом.

Первичные страницы: [1](https://search.clinicalgenome.org/kb/downloads), [2](https://iscn.karger.com/).

## Что не добавлять по умолчанию

Полный dbNSFP5.4, full gnomAD, платный OMIM/VarSome/Mastermind и все дополнительные предикторы не являются условием готовности первого отчёта. FATHMM legacy уже есть в raw public reference cache; нужен правильный mapper, не очередной cancer predictor. GENCODE50 из нового dbNSFP не подменять GENCODE49. RuSeq и All of Us названы в лекциях, но не считать их автоматически подходящей популяцией для казахстанской когорты; ancestry неизвестен. Онкологические COSMIC/OncoKB/JAX-CKB/AMP-ASCO-CAP не нужны для данного germline HCM анализа. Secondary findings требуют своего объёма/согласия; KEGG withdrawn. LOFTEE, Exomiser, MitoTip и другие не подтверждённые в лекциях названия не приписывать преподавателям.

Первое предлагаемое решение: одобрить P01–P03, P05–P06 и P08–P11 как небольшие data/method tasks; P04/P07 — по нуждам shortlist; P12 отложить. До отдельного одобрения реализация отсутствует. Размеры, стоимость передачи и пределы ресурсов проверять до скачивания: исходный общий бюджет $50, VM deadline 2026-11-06T03:46:48.046025Z. Новые рекомендации не продлевают VM и не открывают платные аккаунты.

Машиночитаемый файл: [proposal JSON](../config/medinter-additions-proposal-v2.json). Он содержит 12 proposed items и source audit44; PHI/вариантов когорты/ключей в нём нет.
