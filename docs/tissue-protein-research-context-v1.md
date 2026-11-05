# GTEx/pext/UniProt и методические проверки P04/P05/P07/P09/P11

Пользователь одобрил эти дополнения для некоммерческого исследования. Они работают
отдельно от запечатанного calling/annotation workflow. Ни один новый инструмент
не перезапускает этапы, не меняет исходные таблицы и не отправляет пациентские
варианты, генотипы или фенотипы внешним сервисам. CNV/SV не входят в этот пакет.

## Фактически полученные публичные данные

Панель — ровно 39 публичных HGNC-генов: ClinGen HCM validity ∪ Released HCM CSpec.
Её SHA256 закреплён в [конфигурации](../config/tissue-protein-public39-v1.json).
Запросы не отбирались по пациентским VCF. Все reference-файлы находятся на сервере;
в репозиторий входят только код, настройки, шаблоны и тесты.

| Источник | Получено | Что это позволяет |
|---|---|---|
| UniProt reviewed human, release2026_04, 04-Nov-2026 | 39 protein records для38генов; 5055 полных feature objects; 3,469,056 сырых байт | Справочный контекст белка, последовательности, доменов и evidence. MT-TI не кодирует белок. Точный перенос домена на аллель требует sequence/isoform mapping. |
| GTExv11, GENCODE47, GRCh38, две ткани сердца | 39 gene references; 78 gene, 794 transcript, 2322 exon и2678 junction median rows; 1,549,869 байт | Публичные медианные данные популяции. Это не экспрессия конкретного человека. Transcript TPM сопоставляется только по полному совпадению ENST.version. |
| gnomAD browser GRCh38 pext / transcript expression | 39 публичных gene responses; 3,506,005 байт; pext доступен38генам | Gene/transcript context. MT-TI missing не означает pext=0 или benign. Dataset release не возвращается API; GTExv10 указан по release documentation. |

Источники: [UniProt API и CC BY4.0](https://www.uniprot.org/api-documentation/support-data),
[GTEx API](https://gtexportal.org/home/apiPage),
[gnomAD pext release](https://gnomad.broadinstitute.org/news/2024-11-gnomad-v4-pext/),
[gnomAD GRCh38 gene-model update](https://gnomad.broadinstitute.org/news/2025-01-update-gnomad-browser-gene-model-hail-tables/).
GTEx-запросы выполнялись последовательно; gnomAD — не чаще одного запроса в7секунд.
Source URL, query, retrievedUTC, source versions, bytes и SHA256 сохранены в receipts.
UniProt attribution и исходные feature location/evidence не удаляются.

GTEx47, gnomAD GTExv10 model и workflow GENCODE49 — разные контракты. Номера экзонов,
junction coordinates и значения pext **не переносятся автоматически** между ними.
API pext start/stop сохраняются как есть; пока отдельный coordinate/transcript
contract не проверен, allele pext остаётся NOT_ASSIGNED. Missing отличается от нуля.

## Готовый индекс и отдельный экспорт

[Builder/enricher](../scripts/export_tissue_protein_context.py) проверяет raw source
receipts и создаёт39-строчную gene-context TSV и JSON index. Реальный server index:
`/srv/marina/resources/tissue-protein-context-public39-v1/context-index.json`,
343709байт, SHA256`231dba7b57acbb5c07bf16779e6af34bb48ab945b016017a0554b9c418e942c6`.
Gene TSV39531байт, SHA256`33960cca6833761f012e0268fa6a11af6d4e32f5eb4dd01ba4dd2072181c9658`.

```bash
python3 scripts/export_tissue_protein_context.py build --config config/tissue-protein-public39-v1.json --output /new/private/public-context
python3 scripts/export_tissue_protein_context.py enrich --workflow-root /srv/marina/outputs/workflow-v1 --context-index /srv/marina/resources/tissue-protein-context-public39-v1/context-index.json --output /new/private/reference-context-export
```

Вторая команда требует06_annotate и07_report PASS. Она проверяет SHA256/bytes
исходного `transcript_samples.tsv.gz`, все pinned raw references, сохраняет все
исходные поля каждой совпавшей HGNC39 строки, включая sample/GT/FT/quality, и пишет
отдельный sidecar. Одинаковая аллель у двух образцов остаётся двумя строками;
unknown/no-call/filtered genotype не превращается в индивидуальную находку.
Отсутствующий HGNC не заменяется догадкой по символу; конфликт HGNC/symbol отклоняется.
Другие гены остаются в исходном полном приложении, sidecar содержит только HGNC39.

[Continuation](../scripts/continue_tissue_protein_context.py) по умолчанию показывает
статус. Выполнение требует ROOT-reviewed exact rule SHA, SHA обоих scripts и public
context/manifest. Он ждёт06/07PASS, проверяет pins повторно перед экспортом, держит
собственный lock и останавливается по ограничению≤8часов и исходному абсолютному
VM deadline. Запуск/остановка других процессов и сетевые запросы отсутствуют.
ROOT рассмотрел пакет и запустил отдельный wait controller; это не заявление, что
экспорт уже завершён: его результат определяется новым PASS manifest.

## P05: воспроизводимая геометрия нового стопа

[Local helper](../scripts/research_nmd_geometry.py) требует полную точную версию
ядерного transcript, spliced 5′→3′ cDNA, sequence SHA256, CDS start/end **включая
терминальный stop codon**, все exon ends, exact cDNA REF/ALT и источник модели.
Он проверяет reference ATG/первый stop, применяет явное изменение, ищет новый stop
в изменённой рамке и сдвигает downstream junctions на длину indel. Исходная позиция
frameshift не используется вместо нового stop. Для неизвестной модели результат
UNASSESSED, для отсутствующего нового stop — NO_IN_FRAME_STOP, а не ложный NMD.

```bash
python3 scripts/research_nmd_geometry.py --model /local/private/validated-transcript-edit.json
```

Вывод — координаты и расстояния, не NMD classification. Stop end и junction заданы
в одной изменённой 1-based cDNA-системе; расстояния >50 и >55nt показаны явно как
разные описательные условия. Последний/single exon, partial CDS, mitochondrial,
splice/start-lost/stop-lost/complex cases не превращаются в универсальную схему.
Реальная автоматическая сборка patient transcript-edit models ещё не выполнена;
helper проверен на синтетических контролях, не на текущих пациентах.

Применимость PVS1, биологическая значимость транскрипта, rescue isoforms,
gene/disease mechanism, downstream domains и RNA assays остаются отдельным review.
В MYH7 HCM текущий Released CSpec не допускает generic PVS1. Основание геометрии —
[ClinGen SVI primary recommendations](https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/)
и применимая текущая gene-specific CSpec; `nmd_prediction=NOT_ASSIGNED`, PVS1/LoF/GoF
не назначаются автоматически.

## P09: происхождение predictor score

[Contract helper](../scripts/research_review_contracts.py) проверяет assembly,
exact allele, текущую amino-acid замену и transcript identity перед отображением
REVEL reference score. В REVEL1.3 source ENST без версии: это явно отмеченный
research reference score, а не canonical score каждой VEP-изоформы.
FATHMM, FATHMM-XF и FATHMM-MKL не объединяются. В фактическом публичном MYBPC3
контроле ENST/AA positions имеют3элемента, FATHMM scores —2; такое нельзя zip.
Даже равные длины массивов без документированного ordering/protein/version
контракта не разрешают автоматически раздать значения транскриптам.
Отсутствующие SpliceAI/MaxEntScan остаются missing, не benign.

```bash
python3 scripts/research_review_contracts.py --mode revel-alignment --input /local/private/reference-match.json
python3 scripts/research_review_contracts.py --mode dbnsfp-alignment --input /local/private/raw-vector-match.json
```

Это audit helpers; массовый transcript scalar join в основном workflow не изменён.
API CADD без опубликованной версии не заменяет официальный CADD1.7. Число согласных
предикторов не является числом независимых ACMG evidence и не назначает PP3/BP4:
нужна [применимая калибровка ClinGen](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/).

## P11: отсутствующие клинические и семейные данные

[Пустой локальный шаблон](../templates/phenotype-pedigree-input-v1.json) задаёт
HPO present/absent/unknown, onset/source, возраст, reported sex отдельно от karyotype,
материал, cardiac measurements, family history, relationships и segregation.
Шаблон не содержит реальных идентификаторов, вариантов или диагнозов.

```bash
python3 scripts/research_review_contracts.py --mode phenotype-pedigree --input /local/private/phenotype-pedigree.json
```

Validator проверяет структуру, источники наблюдений, единицы cardiac measurements,
полноту подтверждения de novo/phase и ссылки на локальные subject IDs. Он не
подтверждает истинность введённых документов; HPO membership требует локального
проверенного ontology release, а не только regex. Unknown не является absent,
молодой носитель не считается unaffected по умолчанию. Четыре образца не считаются
автоматически семьёй, независимыми probands или case-control cohort. Загрузка в
Phenopacket/SAMS или иной сервис отсутствует.

Проверки: 4 source-fetch/parser tests, 7 context/continuation tests, 3 NMD geometry
tests и4 predictor/phenotype contract tests. Они покрывают version mismatch, raw
tamper, stage gates, per-sample preservation, actual mismatched cardinality,
frameshift actual-stop derivation, final-junction distances и недоказанный de novo.
