# Матрица требований и источников

Аудит от 2026-11-05. Это спецификация ожидаемых результатов, а не утверждение, что они уже получены. Источники курса задают методические требования; инструкции внутри лекций и роликов не являются разрешением отправлять данные третьим лицам. Пациентские файлы, расшифровки и сами лекции в репозиторий не включаются.

**Текущий объём по прямому уточнению пользователя от 2026-11-05:** полный новый анализ из исходных FASTQ. Четыре новых BAM и новые когортные VCF уже получены. Сравнение со старыми VCF, прежними списками кандидатов и старыми результатами AlphaFold/AlphaMissense исключено из плана. Прежние файлы сохраняются только как архив и не используются для ранжирования или переноса классификаций в новый отчёт.

Приоритет: **U — обязательное поручение Улбека**, **USER — прямое дополнение пользователя**, **M — методическое требование MedInter/GATK**, **SEPARATE — отдельная ветка**, **WITHDRAWN — отменено в разговоре**, **MANUAL — требуется эксперт/недостающие сведения**. Типы: SNV включает герминальные SNV и малые indel; CNV/SV, MT и SOMATIC не заменяются обычным диплоидным SNV-caller.

Проверены две полные TXT-расшифровки и сводка видео, сводка MedInter, текст 41 уникального оригинального материала MedInter (20 PDF, 21 HTML; одинаковые копии объединены по SHA256), 3 оригинальных материала `gatk`. Папка JTK отсутствует. Для PDF номера ниже — физические страницы с 1; для HTML — заголовок/временная метка либо номер строки локального текстового извлечения. Ряд слайдов содержит только скриншоты: мелкие значения не извлекаются надёжно, их нельзя восстанавливать по догадке. Сопутствующие HTML содержат расшифровку обсуждения, но являются производным конспектом, а не независимым доказательством точных чисел на изображении. M26 — неудавшийся пустой конспект; требований из него не выведено.

## Поручения из двух роликов

V1 = `video/20260904_151456.mp4`, 12:55; V2 = `video/20260904_152755.mp4`, 10:57; полные TXT — в `video/transcripts/`. Автоматическое распознавание Whisper не имеет надёжной диаризации; отдельные имена программ неоднозначны. Имена старых кандидатов, числа 27/57/200/220 и исторические утверждения о пересечениях — обсуждение предыдущей работы, не параметры новой обработки.

| ID / источник | Требование | Поля и проверяемый результат | Ветка / приоритет |
|---|---|---|---|
| U01 V1 00:00–00:14, 00:55–01:48 | Пообразцовые таблицы и общая таблица ведущих кандидатов | Нормализованный variant key, ген, транскрипт, GT/DP/AD/GQ, эффект, источники; основная таблица отдельно от полного приложения | SNV / U |
| U02 V1 00:25–00:52; V2 00:28–01:28 | Пересечения между 2, 3 и 4 образцами | Матрица присутствия, число носителей, множество образцов; отдельно пересечения вариантов и генов; пропуски/низкое покрытие не считать REF | SNV / U |
| U03 V1 01:49–02:18, 10:08–10:44 | Анализ кардиогенов / ГКМП | Новый версионированный набор генов из актуального источника, HGNC ID, дата, disease/mechanism; индивидуальный клинический контекст указывается отдельно | SNV / U + MANUAL |
| U04 V1 03:59–04:45; V2 09:39–09:48 | ClinVar для конкретных вариантов, различать P и LP | Allele/Variation ID, significance, condition, review stars, conflict, evaluation date, release; отсутствие записи отдельно | SNV / U |
| U05 V1 04:51–05:58 | Сопоставить предикторы аминокислотных замен | SIFT, PolyPhen, FATHMM: score/prediction/version/missing; четвёртое название «Cart» может быть CADD, **не подтверждено**; concordance как описание, не 3–4 независимых ACMG-критерия | SNV / U |
| U06 V1 06:08–07:21; V2 09:49–10:14 | Отдельно stop-gain и stop-loss, положение изменения | HGVSg/c/p, strand, exon, CDS/protein position/length, relative position, NMD assessment, start-loss/frameshift/splice отдельно; stop-loss не автоматически LoF/PVS1 | SNV / U |
| U07 V1 07:21–08:01 | Enrichment/KEGG сначала предложено, затем снято | Не запускать по умолчанию; вернуть только при отдельном новом поручении | WITHDRAWN |
| U08 V1 08:02–08:43; уточнение пользователя 2026-11-05 | Историческое упоминание AlphaFold/AlphaMissense; текущую задачу задаёт USER02 | Прежние результаты не используются и не сравниваются. AlphaMissense применяется заново к новым вариантам по USER02 с версиями источников и проверкой идентичности. Моделирование AlphaFold не является обязательным этапом | Исторический источник / USER OVERRIDE |
| U09 V1 08:53–09:08, 09:31–09:42 | Сначала таблицы, затем гистограммы типов вариантов | Типы/число на образец, одинаковый набор фильтров, легенда, значения, знаменатель; SNV/indel и consequence не смешивать | SNV / U |
| U10 V1 10:08–10:44; V2 01:28–04:21 | ACMG для приоритетного короткого списка | По критерию: применён/не применён/неизвестно, сила, доказательство/PMID/дата, конфликт, экспертное заключение; не классифицировать автоматически тысячи сырых вызовов | SNV / U + MANUAL |
| U11 V2 01:28–04:21 | Платформа классификации, если доступна, иначе вручную | Название похоже на VarSome, **неуверенно**; независимая локальная evidence table обязательна, платная API не предпосылка | SNV / U + лицензия |
| U12 V1 10:54–12:24 | Частоты в других популяциях и повторяющиеся значимые кандидаты ГКМП | AF/AC/AN/Hom, ancestry/population и release, max/popmax с определением, cohort carrier count отдельно; 4 образца не доказывают биомаркер | SNV / U |
| U13 V2 00:00–00:24 | Указывать варианты и позиции, а не только гены | Build+CHROM+POS+REF+ALT; rsID при наличии; HGVS с версией транскрипта | SNV / U |
| U14 V2 04:29–06:33, 06:57–07:24, 08:42–09:14 | Клинические характеристики когорты и критерии включения | Предоставленные диагноз/анамнез/ЭКГ/эхо/биохимия и критерии; отсутствующие сведения = unknown; клиническое приложение и Materials/Methods; сейчас можно готовить результаты | MANUAL / U, ожидание данных |
| U15 V1 03:31–03:59; V2 07:26–07:55 | Рукопись Word по шаблону журнала | Introduction, Materials and Methods, Results, Discussion, tables/supplements; шаблон/журнал пока не определён | U + MANUAL |
| U16 V2 08:02–08:40 | Обсудить основные/повторяющиеся варианты по литературе | PMID/DOI, variant-specific evidence, disease/mechanism, ограничения; не обсуждать каждый сырой вызов | U + MANUAL |
| U17 V1 03:31–03:59; V2 10:21–10:30 | Совместное редактирование / материал для комментариев | Подготовленный документ/рабочая версия; загрузка в Google Drive и адресаты — отдельный шаг, без автоматической внешней отправки | U + MANUAL |
| U18 V1 03:05–03:21 | Навигация: основные гены, интерпретация, пересечения, связанные находки, сердечный скрининг | Индекс результатов/листы с описанием области; secondary findings зависят от согласия, не из самого упоминания | U + MANUAL |
| U19 V1 12:44–12:53 | Неясная завершающая реплика | Неразрешённая неоднозначность; новое требование не придумано | MANUAL |
| USER01 последнее прямое поручение | Top 10 и Top 20 на образец и когорту, отдельными таблицами | Прозрачный rank/score components/QC/source completeness; не дополнять пустоту слабыми кандидатами до квоты; ranking не pathogenicity classification | SNV / USER MUST |
| USER02 новое прямое поручение пользователя | AlphaMissense обязателен в таблицах и интерпретации | Score/class, источник/release, transcript.version и protein variant, статус точного сопоставления; для non-missense — not applicable, для отсутствующего/несопоставимого score — missing reason. Отдельный исследовательский показатель; не автоматическая ACMG-классификация. Требование относится к AlphaMissense и не меняет исторический смысл U08 или объём AlphaFold-моделирования | SNV / USER MUST |
| USER03 новые прямые поручения пользователя | Собрать понятные результаты на сервере и в папке проекта, затем один ZIP | Локально: 4 BAM+BAI, принятый общий VCF+index, QC/reports/code/settings; индивидуальные gVCF остаются на сервере. Новые backup-файлы проверяются по завершению/размеру с явной меткой size_only; ранее полученная SHA-история сохраняется. Базы/кэши и raw в архив не входят. ZIP собирается локально только после принятого cohort VCF | Операционный результат / USER MUST |

## Методическая матрица всех разделов

| ID / источник | Требование | Поля / результат | Ветка / приоритет |
|---|---|---|---|
| M01/M21, PDF 4–9 | Сначала определить цель, объём и ограничения анализа | Research/diagnostic scope, claimed WGS versus observed coverage, negative-result limitations; не гарантировать исключение болезни | M |
| M02/M37, PDF 7–27 | Эффект зависит от структуры гена, изоформы, ткани, псевдогенов и альтернативного сплайсинга | Transcript/version, CDS/exon/strand, tissue relevance, paralog/pseudogene ambiguity | SNV / M |
| M03/M41, PDF 4–31 | Сохранять типы вариантов и последствия раздельно | SNV/indel/CNV/SV/repeat; synonymous/missense/nonsense/splice/start/stop/frameshift; multi-consequence per transcript | SNV + SEPARATE |
| M04/M25, PDF 3–19 | Наследование и механизм болезни не выводятся из одного GT | AD/AR/X/MT/digenic, cis/trans/phase, de novo unknown, penetrance/expressivity/imprinting/anticipation/mosaicism | MANUAL / M |
| M05/M38, PDF 5–19 | Контроль исходных данных и границ технологии | Pairing, read quality, adapter/library information when known, depth **и breadth**, uneven coverage; не маркировать WES/WGS только по названию файла | M |
| M05 PDF 19–37; M20/M28 | Невидимые SNV-пайплайну классы и иной биоматериал | Ограничения CNV/SV/repeats/methylation/imprinting/LOH/mosaicism; MLPA/qPCR/CMA/FISH/Sanger/Hi-C/fragment analysis как возможное подтверждение, не автоматический заказ лаборатории | SEPARATE + MANUAL |
| M06/M27, PDF 3–17 | Прослеживаемая цепочка FASTQ→BAM→VCF→annotation | Checksums, full reference lock, tools/commands/container, read groups, duplicate flags, QC, caller/filter policy | M |
| M06/M27, PDF 18–44 | Фильтры по DP/GQ/AD/AB, частоте, эффекту и качеству поддержки | Raw annotations сохранены, explicit thresholds/config, pass/fail/reason, missing≠benign, FORMAT versus INFO DP; FILTER `.` = not filtered | SNV / M |
| M06 PDF 10,35–44; M19/M23 | Визуальная проверка в IGV | Build/reference match, strand/read position/mapQ/baseQ/duplicates/soft clips/nearby indel, manual-review state; отсутствие проверки не PASS | SNV / MANUAL |
| M07/M29, PDF 2–14; G42 3–12 | Зафиксировать точный референс | Assembly accession, FASTA+SHA, ALT/decoy/masking, contig dictionary/aliases; не смешивать hg19/hg38/T2T | M |
| M07 PDF 16–49 | Единые идентификаторы генов/транскриптов/белков | HGNC+aliases, NCBI Gene, RefSeq, Ensembl, MANE, LRG, UniProt; MANE Select/Plus Clinical и другие изоформы; TSL не clinical proof | SNV / M |
| M07 PDF 50–57; M12/M36 | Уникально и корректно обозначить вариант | Нормализация left-aligned VCF и allele split; HGVS 3′-правило отдельно; build, accession.version, c./g./p.; protein-only и rsID-only недостаточны | SNV / M |
| M12 PDF 20–23; M36 24:20–29:15 | Валидировать HGVS | VariantValidator/Mutalyzer; offline code/resources если подготовлены, иначе manual queue; не отправлять варианты пациента в веб автоматически | SNV / M + MANUAL |
| M08/M32, PDF 5–17 | Сопоставить клинику и фенотип | Анамнез/дебют/течение/исследования/родословная, HPO present/absent/unknown, OMIM/GeneReviews/ClinGen/HPO/Orphanet; SAMS показан как пример, без загрузки пациентов | MANUAL / M |
| M09/M33, PDF 3–33 | Клиническая запись не равна диагнозу | ClinVar/LOVD variant-condition assertions, stars/conflicts/submitter/date; phenotype and inheritance; VarSome/Franklin вспомогательные агрегаторы | SNV / M |
| M10/M34, PDF 3–28 | Все 5 классов ACMG/AMP и все категории доказательств | P/LP/VUS/LB/B; codes PVS1, PS1–4, PM1–6, PP1–5, BA1, BS1–4, BP1–7; состояние/сила/доказательство отдельно, unknown допустимо | SNV / M + MANUAL |
| M10 PDF 8–11; M34 05:04–07:35 | Частота в контексте болезни | gnomAD, RuSeq, All of Us (HTML); BA1/BS1 disease-specific, PM2 с актуальной ClinGen силой; отсутствие в БД не само доказательство P; ExAC исторический источник | SNV / M |
| M10 PDF 10–19; M34 07:35–21:35 | Сегрегация, de novo, phase, функциональные и вычислительные данные | PS2/PM6 требует parental evidence; PM3 cis/trans; PS3/BS3 assay validity; PS4 case enrichment не из 4 выбранных образцов; PP4 phenotype specificity | SNV / M + MANUAL |
| M10 PDF 27–50; M11/M31 | ClinGen ген-/болезнь-/критерий-специфичные рекомендации | Ссылка+версия VCEP/criteria specification, gene validity; не применять устаревшее PP5/BP6 как независимое ClinVar-доказательство; не считать коррелированные предикторы независимыми | SNV / M |
| M11 PDF 3–21; M31 | PVS1 не равен consequence HIGH | Established LoF mechanism, biologically relevant transcript/exon, NMD/last exon/last 50bp penultimate exon, lost protein fraction/domain, splice rescue/in-frame exon skipping, alternative start; strength reduction / not applicable / unknown | SNV / M + MANUAL |
| M11 PDF 13; M33 50:00; M34 18:09 | Контролировать применимость предикторов | SpliceAI для сплайсинга; SIFT/PolyPhen/FATHMM/REVEL/AlphaMissense/Varity; `MT` и `MUT Assessor` показаны такими сокращениями — не переименованы по догадке | SNV / M |
| M13; M30 00:00–конец | CNV/SV требуют своего caller/QC/классификации | Coordinates+uncertainty+SVTYPE/CN/end/size, support split/discordant/depth, genes/exons, HI/TS ClinGen, DECIPHER/DGV/gnomAD constraint, pLI/LOEUF/pHaplo/pTriplo; ISCN/HGVS по доступной точности | CNV/SV / SEPARATE |
| M14/M35, PDF 2–20; G42 9–10 | MT — особая референсная и интерпретационная ветка | rCRS NC_012920.1, circular/shifted-reference policy, NUMTs, heteroplasmy/DP/AF/detection limit, haplogroup/tissue, MITOMAP, MT-specific ACMG (humu.24107), MT genetic code; не применять ядерный NMD автоматически | MT / SEPARATE |
| M15/M40, PDF 25–55 | Не смешивать герминальную классификацию с опухолевой онкогенностью/ответом терапии | Tumour type/purity, matched normal, VAF, cis/trans, assay sensitivity; AMP/tier/actionability; COSMIC/My Cancer Genome/DoCM/PharmGKB/DGIdb/OncoKB/JAX-CKB/FDA; нет опухолевых входов — disabled | SOMATIC / SEPARATE + лицензия |
| M16/M39, PDF 3–13 | Случайные и вторичные находки — отдельная область согласия | Consent/age/scope, versioned ACMG SF list, disease+inheritance and P/LP; список SF не ограничивает весь ACMG-анализ 57 генами | MANUAL / SEPARATE |
| M17/M24, PDF 7–58 | Пенетрантность, возраст и экспрессивность | Gene–disease validity, phenotype/age of relatives, mosaicism; healthy carrier не автоматический benign; ClinVar P при слабой связи гена с болезнью не окончательный вывод | SNV / M + MANUAL |
| M18/M22, PDF 2–14 | Заключение объясняет интерпретацию и ограничения | Methods/reference/resources/dates/QC, genotype+classification+phenotype fit, alternatives, findings withheld/unknown, reanalysis triggers; VUS не клиническое решение | M + MANUAL |
| M19/M23, PDF 3–13 | План подтверждения и реклассификации | Technical-review queue, phase/segregation/functional evidence missing, alternative second allele; лабораторные действия только как рекомендации, не контакт/заказ | MANUAL / M |
| M20/M28, PDF 2–30 | Выбор метода и переанализ зависят от клиники/материала | Variant-class limitations, appropriate tissue, prior tests, updated resources/new phenotype; отрицательное NGS не исключает генетическую причину | M + MANUAL |
| G42 13–18; G44 32:00–99:00 | Управляемое распараллеливание и воспроизводимость | Dry-run, explicit resources, readgroup SM/RG/LB/PL provenance, idempotence, logs/resume; не выравнивать full-genome reads только на chr1 | M |
| G43 2–8 | Caller/BQSR/cohort policy и CPU/GPU версии | Варианты Parabricks/GATK/DeepVariant — альтернативы; BQSR/known-sites policy и совместимые версии фиксировать, совместимость не равна точному совпадению всех версий; small cohort ≠ автоматически VQSR | M |

### Дополнительные именованные ресурсы и границы уверенности

Финальная сверка имён исходных HTML дополняет перечень:

| Источник | Ресурс / требование | Данные / применение |
|---|---|---|
| M33 L20, M34 L30 | Sherloc | Уточнённая система классификации; ручное сопоставление правил, не второй независимый голос рядом с ACMG. [Авторская публикация](https://pubmed.ncbi.nlm.nih.gov/28492532/). |
| M33 L119/124/202 | Mastermind | Поиск публикаций по варианту, доступ/лицензия отдельно; local evidence worksheet с PMID/DOI. [Авторский продукт](https://www.genomenon.com/mastermind) указывает freemium и лицензируемые API/flat files. Бесплатный полный dump не предполагается. |
| M33 L136–137/271; M34 L130/311 | MetaLR, MaxEntScan | Предикторы из примеров Franklin/ACMG; версия/score/applicability/missing сохраняются, включение только с допустимым локальным источником. |
| M35 L74–79/126–127 | Helix, `TOOLST` | MT-популяционная таблица и показатель в конспекте. Точный продукт/release и расшифровка TOOLST не установлены; не подставлять MitoTip по догадке. |
| M24 L182–183 | CanVIG-UK, gene-specific BRCA1/BRCA2 guidance | Пример ген-специфичных порогов, не общий порог для ГКМП. [Официальная guidance](https://www.cangene-canvaruk.org/canvig-uk-guidance). |
| M32 L121/134–138/234/243 | LOINC, Phenopacket | Кодирование измерений и формат структурированного фенотипа, показанные в HPO/SAMS; не разрешение загружать пациента. |
| M36 L154; M29 L137 | CCDS; `LenVar` | CCDS — поле референсного coding sequence; LenVar — неясное имя в HTML при обсуждении LRG, вероятная ошибка конспекта; точную БД не выдумывать. |
| M40; G44 L7 | NCCN/AMP/ASCO/CAP; 1000 Genomes | Онкологические руководства — отдельная ветка; 1000 Genomes в GATK — источник учебных данных, не обязанность менять исходные 4 образца. |

## Полнота ресурсов и способ доступа

Машиночитаемый перечень находится в `config/evidence-sources.example.json`. Это **каталог требований**, не runtime-конфиг и не подтверждение загрузки. Для каждой локальной базы фактический lock должен содержать URL, дату/release, assembly/transcript version, SHA256, размер, условия использования. Незаполненные ресурсы дают `missing/unavailable/not_applicable`, а не выдуманные оценки. Ни один пациентский вариант автоматически не отправляется внешнему API.

Проверенные первичные страницы доступа (2026-11-05):

| Ресурс | API-free/offline путь и ограничение |
|---|---|
| ClinVar | Публичные FTP XML/VCF/TSV; VCF только precise-location summary и не весь набор представлений/неопределённых CNV. Нужны release+SHA, для conflicting submissions XML/TSV. [NCBI data guide](https://www.ncbi.nlm.nih.gov/guide/data-software/), [ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/). Не считать весь контент автоматически public-domain без условий submitters. |
| ClinGen | CSV gene-disease validity, dosage, VCEP pathogenicity и criteria registry; можно скачать публичные файлы и читать локально. [Официальный каталог](https://search.clinicalgenome.org/kb/downloads). |
| Ensembl/VEP | `--offline` с согласованными cache/FASTA, без database calls; SIFT/PolyPhen precomputed включены в human annotation infrastructure. Отдельные права каждого предиктора сохраняются. [VEP options](https://static.ensembl.org/info/docs/tools/vep/script/vep_options.html), [VEP FAQ](https://rc-metazoa.ensembl.org/info/docs/tools/vep/vep_faq.html). |
| AlphaMissense | Public precomputed predictions, CC-BY 4.0 с 2024-03-13; software Apache не заменяет лицензию scores. [DeepMind announcement](https://deepmind.google/blog/a-catalogue-of-genetic-mutations-to-help-pinpoint-the-cause-of-diseases/), [data README](https://storage.googleapis.com/dm_alphamissense/README.pdf). |
| SpliceAI | Downloadable precomputed annotations и local inference; source PolyForm Strict, models CC-BY-NC4.0, scores academic/not-for-profit; коммерческое применение требует других прав. Не включать автоматически без подтверждения применимости лицензии. [Illumina](https://github.com/illumina/spliceAI). |
| FATHMM | Локальная программа GPLv3; выбрать inherited-disease predictor/version, не cancer model; исходный HGMD training set не распространяется. [Авторские downloads](https://fathmm.biocompute.org.uk/downloads.html). |
| REVEL | Авторская страница с precomputed downloads; лицензию конкретного пакета сохранить, не выводить права из наличия URL. [Downloads](https://sites.google.com/site/revelgenomics/downloads). |
| OMIM | Web/manual; batch/API требует регистрации/ключа/разрешённых downloads, коммерческое использование/redistribution ограничены. Не скрейпить полный OMIM. [Авторская публикация](https://pmc.ncbi.nlm.nih.gov/articles/PMC4383985/), [NCBI legal statement](https://www.ncbi.nlm.nih.gov/sites/books/n/handbook/pdf/). |
| VarSome | Публичный web не означает бесплатную full API: подписка/токены и платный Standard. Необходимость API отсутствует: локальные upstream databases + evidence worksheet. [Официальные packages](https://docs.varsome.com/en/varsome-api-packages-explained). |
| gnomAD/HPO/HGNC/UniProt/LOVD/Orphanet/MITOMAP и иные | Публичный browser ≠ автоматически unrestricted full dump. Использовать официальные release/downloads и сохранить условия конкретного файла; в этой проверке полные актуальные условия части сайтов не удалось прочитать (JS/404/access), поэтому license status в примере `review_required`, без обхода или выдуманного CC0. |
| Franklin/DECIPHER/COSMIC/OncoKB/JAX-CKB/ISCN | Manual/licensed или access-dependent; автоматическое копирование/загрузка данных пациента отключено; ветки не обязательны для текущего germline HCM. |

Неподтверждённые названия: четвёртый predictor «Cart» (возможно CADD) и платформа во V2 (возможно VarSome). Сокращения `MT`, `MUT Assessor` сохраняются как в M33, без догадки о полном названии. CADD не назван достоверно в извлечённом тексте лекций. HGMD, dbNSFP, MutationTaster, GERP/phyloP/phastCons, MitoTip, CIViC, Exomiser и PanelApp **не подтверждены как требования в этих источниках**; их нельзя приписывать Улбеку/лекциям. Дополнительный выбор их разработчиком требует отдельного rationale и лицензии.

## Критерии приёмки и очередность

1. Завершённые alignment/QC с проверенными входами и checksum; оценить измеренное покрытие и ограничения каждого образца.
2. Зафиксировать germline calling/joint-genotyping/filtering policy, ресурсы и версии. chr1 technical pilot не является четырьмя полными клиническими результатами.
3. Локальная аннотация с полными source/missing flags; подготовить все U01–U13, USER01 и USER02 таблицы, отдельные LoF/stop sheets, overlap и графики. Stage 7.2 — приоритетные кандидаты и объяснимый ranking; AlphaMissense обязателен с явной применимостью и точностью сопоставления.
4. Экспертный evidence review ACMG/ClinGen, фенотип/семья и литература по короткому списку. Без клиники явно incomplete; без семейных данных не заявлять de novo/cis/trans.
5. Рукопись/приложения и очередь ручных/лицензионных действий. Не выдавать исследовательское ранжирование за диагноз.
6. CNV/SV/MT/somatic/secondary findings остаются отдельными disabled/conditional ветками, а не молчаливо выполненными задачами.

Локальные полные расшифровки и извлечения сохраняются в исходных папках и ignored `.local/`; ниже каталог даёт проверяемое соответствие source IDs оригиналам. Даты, числа и клинические примеры лекций не переносятся как наблюдения текущей когорты.

## Каталог оригиналов

<!-- Generated source catalogue follows; SHA256 identifies the source lecture, not patient data. -->

- **M01** — `med_inter\1.1 Когда и зачем назначать генетическое тестирование_ клинические сценарии и путь пациента_watermark (1).pdf`; 16 pages; SHA256 `640a4e0ecf92c2d31032889600cf12f9b6979aabb0444d569f5b325e7d80008f`.

- **M02** — `med_inter\1.2 NGS без магии_ ключевые принципы молекулярной биологии которые реально нужны врачу.pdf`; 27 pages; SHA256 `bf55886ec9324ab9e1fea963a73a6903155bb7fccc82b922d288d10b4f1d231a`.

- **M03** — `med_inter\1.3 Типы генетических вариантов_ как мутации влияют на структуру и функцию гена, генома, белка.pptx_watermark (1).pdf`; 34 pages; SHA256 `124939dc01277a95288ce8ad2ef98563dfcdd3f1b2f63949ab14b38670059f91`.

- **M04** — `med_inter\1.4 Как наследуются заболевания_ от простых моделей к сложной реальности.pptx_watermark (1).pdf`; 20 pages; SHA256 `64c3ab969816507bc793b4a62778dd8c05bf8d9ef00d8fd218e87c49428fbb13`.

- **M05** — `med_inter\2.1 NGS в реальной диагностике_ возможности_ограничения и когда нужны другие методы_watermark (1).pdf`; 37 pages; SHA256 `c5a123832c5eea2e772b477054229355d48c177aa8bda93822672f95b32affe9`.

- **M06** — `med_inter\2.2 Что происходит после секвенирования_ биоинформатический анализ (пайплайн), метрики качества и что такое IGV.pptx_watermark (1).pdf`; 45 pages; SHA256 `0f3d5966bd96b2e220e9bbe14214c99c68689a1cb065c681b1e3d2964c140df0`.

- **M07** — `med_inter\2.3 Референсный геном и базы данных_ как не утонуть в источниках и использовать их осмысленно.pptx_watermark (1).pdf`; 58 pages; SHA256 `499728490945d7965858640bb5883fc7fec8b3039f90fef84785d315736252e2`.

- **M08** — `med_inter\3.1 Как читать клинические данные пациента для генетического анализа_watermark (1).pdf`; 19 pages; SHA256 `80fced947459ea63b7e9f6aaab0e0dfc72f22a8a05a5c0a78aa73338d8e4c7a6`.

- **M09** — `med_inter\3.2 Главные базы данных для интерпретации вариантов_ что использовать и когда.pptx_watermark (1) (1).pdf`; 34 pages; SHA256 `a6bd1c729dbc770a9d44e621e01da0954a5c29368b77d6e74fb842c7e17b9954`.

- **M10** — `med_inter\3.3 ACMG без формализма_ как на практике классифицировать варианты_watermark (1).pdf`; 52 pages; SHA256 `254624cc0b66f1246ad85efcf8aaefa15dfa002707e988b16810dc8116345501`.

- **M11** — `med_inter\3.4 Критерий PVS1_ когда loss-of-function действительно значим_watermark (1).pdf`; 23 pages; SHA256 `002699d45e0c65830728efed81a62e5bbc763cbde0092c8e0f0cae4d8568e93c`.

- **M12** — `med_inter\3.5 Как правильно записывать варианты_ номенклатура без ошибок и путаницы_watermark (1).pdf`; 23 pages; SHA256 `1f07242faf3549274c0b75d57b9b332ec29181eda1d4d22228edcb6689ae0d03`.

- **M13** — `med_inter\4.1.2 CNV и структурные перестройки_ как их находить и интерпретировать. Номенклатура_watermark.pdf`; 14 pages; SHA256 `e9aa99a6c8c1e400c790d4b2563228255892a4f74097fe572ad1f7fd22f18983`.

- **M14** — `med_inter\4.2 Митохондриальная ДНК_ особенности анализа и интерпретации_watermark (2).pdf`; 20 pages; SHA256 `756b6ee30375fd4888a78d3937d3a5d35540437fa76d160bf001244cc9ec7d18`.

- **M15** — `med_inter\4.3 Соматические мутации в онкологии_ от секвенирования до клинического решения.pptx_watermark (1) (1).pdf`; 61 pages; SHA256 `e043cd7d5230c4a44edb5c8b8ed930a4124b3b61945cb4736f27b63f802e8fb0`.

- **M16** — `med_inter\5.1 Случайные находки_ что делать, если нашли не то, что искали_watermark (1).pdf`; 16 pages; SHA256 `e24a2af37361b4840f9e7601027016ea5655514ef50574fbb02510a13580f76a`.

- **M17** — `med_inter\5.2 почему один и тот же вариант ведёт себя по_разному_ пенетрантность.pdf`; 58 pages; SHA256 `82035f03cb6d2d4b6793d8f481f22fe0bc3c45361581a3c66e48db6d36ffd1a5`.

- **M18** — `med_inter\5.3 Как читать генетическое заключение и использовать его в практике.pptx_watermark.pdf`; 16 pages; SHA256 `50a46f80e1959ce5c937f138f0c4216b2c045a6e7ad8124acdb92dcba5a202c1`.

- **M19** — `med_inter\5.4 Верификация результатов_ когда и как перепроверять варианты_watermark (1).pdf`; 13 pages; SHA256 `9f306b294cd34b8a0275647ff360fc6e916b984e61f7f10429f7268573ae71c7`.

- **M20** — `med_inter\5.5 Как выбрать правильный метод тестирования под задачу пациента.pptx_watermark (1).pdf`; 31 pages; SHA256 `a975ad9a0f9bb763ef61a656585659549eca9ead6290e5106be407bea1077614`.

- **M21** — `med_inter\NGS в наследственных заболеваниях цели анализа и путь пациента (4).html`; HTML pages; SHA256 `7071505f00c6e6cbdafbaf8b1379b13bf05aa23b7e82c0ad4665435601980ede`.

- **M22** — `med_inter\Session 10_28_2026 (1).html`; HTML pages; SHA256 `2c49c4a85921912b1c77bf7239d87ddb1059a471ae8ce2a86f4c145741e2e9ce`.

- **M23** — `med_inter\Session 10_28_2026 (2).html`; HTML pages; SHA256 `0fba292ac3ccf5b493a19519df69b0bd9dd033530d90c98a653e0f0f29802a0e`.

- **M24** — `med_inter\Session 10_28_2026.html`; HTML pages; SHA256 `cca37d14e989b2f0df8bec0679230e3a37ab6b4b7b47d3cac02701e42a1fe7f7`.

- **M25** — `med_inter\Summary 2026-10-24 13-05 (2).html`; HTML pages; SHA256 `f8694b9b4db68e50edb22ef89026cb4dcc979711a6adbd9ba80f4be91333bab9`.

- **M26** — `med_inter\Summary 2026-10-28 15-39.html`; HTML pages; SHA256 `7ed9dccd0e5f6980a36f6c285144a62ed11829a6eb095b343ed556bda8796dce`.

- **M27** — `med_inter\Биоинформатический анализ данных секвенирования нового поколения пайплайн, контр (1).html`; HTML pages; SHA256 `9c81133625b7fe77ea02fa3d6c82785a9501d3568c1b6f3f9cf678d170c419bf`.

- **M28** — `med_inter\Выбор метода генетической диагностики фенотип, биоматериал и ограничения тестов (1).html`; HTML pages; SHA256 `041ec372a32bdf7208296797132650eec39dc3d0ccaa7c6418230e4bd49f6e0b`.

- **M29** — `med_inter\Занятие 25.10.2026.html`; HTML pages; SHA256 `533ef5d9dcd21a3d8ef6fb1138d4cd2d33a3255c5aea602cf3311e7860db8dd1`.

- **M30** — `med_inter\Занятие 26.10.2026 (1).html`; HTML pages; SHA256 `4e594c6d904ab551d13297211d36c5f28a0d57b49568b60e497c3708cc08f116`.

- **M31** — `med_inter\Интерпретация критерия PVS1 когда потеря функции действительно значима (1).html`; HTML pages; SHA256 `de456f765ae9e0b27ff5aba0375c1dfb5ef5c52fd5386c1cfce5267ca15b5758`.

- **M32** — `med_inter\Как читать клинические данные пациента для генетического анализа (1).html`; HTML pages; SHA256 `2db4e8fda607af7b769a8190e1ff238a832166c9fc3c90d78d5a1fc7366f1392`.

- **M33** — `med_inter\Клинические базы данных и ресурсы для интерпретации генетических вариантов (1).html`; HTML pages; SHA256 `361cc1c57aeb9168fad0a92bb66301b976b1f5f4f68544a073787e0a66cb7b6a`.

- **M34** — `med_inter\Критерии ACMG практическая классификация генетических вариантов (1).html`; HTML pages; SHA256 `a12af730eeb72bed9572cd3ccb1af7abd4a97acc91b5426f055b9781bc90668e`.

- **M35** — `med_inter\Митохондриальная ДНК особенности анализа и интерпретации (1).html`; HTML pages; SHA256 `a274bbe19683c8c7e1ca9b85dbdcb389f309d2589307a1637c4c12ddb9e4dbd8`.

- **M36** — `med_inter\Номенклатура генетических вариантов как записывать изменения по HGVS (2).html`; HTML pages; SHA256 `00b3929beb9306110b459e78bfb6522186950e6327ea0bb6ce81c52c8ab7d8fd`.

- **M37** — `med_inter\Общие принципы молекулярной биологии в основе NGS (1).html`; HTML pages; SHA256 `911f23f52ca11b4e49aa2a7d4b1f56dde9e0f4e228f483b866d076a1e825b483`.

- **M38** — `med_inter\Секвенирование в диагностике возможности, ограничения и выбор метода (1).html`; HTML pages; SHA256 `6662ff61e8c6829a44dbec877526733a689adfd080fcabf46457a9801d32f72a`.

- **M39** — `med_inter\Случайные и вторичные находки что делать, если нашли не то, что искали (1).html`; HTML pages; SHA256 `c930151a99eb79124769fcb5477d29e0b50dc03c53439a61a86fd67bbd92502d`.

- **M40** — `med_inter\Соматические мутации в онкологии от секвенирования до клинического решения (1).html`; HTML pages; SHA256 `a2e9345ab97497e4e080abde56792ebb66d173c39325c98b215dd2ea8d2a7c77`.

- **M41** — `med_inter\Типы мутаций и генетических вариантов от изменений генома до вариантов в отдельн (1).html`; HTML pages; SHA256 `9d3a4dc56ed03b499d875e1c74bd44517901f245786579d4927b4060348de5a8`.

- **G42** — `gatk\13. Анализ Y и MT хромосом 1.pdf`; 18 pages; SHA256 `faeb4efb33b2708761022ea94e76094835a0909adb746f3e7d88eabc174bcc93`.

- **G43** — `gatk\14. Анализ Y и MT хромосом P2.pdf`; 8 pages; SHA256 `80a9719ba4a11fed2b0bc800cb690de910a27bef79b3d72a61645d5bbdf0d81d`.

- **G44** — `gatk\Геном человека, референсные сборки и особенности митохондриальной и Y-хромосомы.html`; HTML pages; SHA256 `52d61bee5956fe45e48f1dea471eeaa93dd5b59ce316bf6c2f31091bd0066cd4`.
