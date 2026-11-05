# Публичная HCM-панель и новые evidence adapters

Пользователь подтвердил **некоммерческое исследование**. Новая приватная approval
metadata хранится отдельно: исторические download receipts с прежним unknown-purpose
не переписываются. Остаются техническая валидация и отдельный review binding; purpose
больше не является блокировкой использования допустимых NC-ресурсов.

## Область и фактические данные

Новая панель v2 — объединение 38 точных HGNC genes из 40 ClinGen validity curations
с MONDO:0005045 и восьми current Released HCM CSpec gene specifications: **39 genes**.
TNNT2 добавлен по CSpec; в текущем validity CSV для него есть DCM/ARVC, но нет точного
HCM MONDO. Эти классификации не повышаются до HCM validity. Исходный panel38 не
перезаписан; TNNT2 получен отдельными incremental API/CADD jobs. У 38 генов есть
поддерживаемые nuclear CDS GENCODE49; MT-TI — отдельная MT-область, не нулевое покрытие.

Orphadata/CSpec bridge v2: 78 строк, 39 HGNC genes, все восемь released HCM specs.
`public-hcm-gene-evidence.tsv` содержит компактные source/version/SHA/IDs и пути
к полным rules. Он поддерживается существующим `gene_evidence_tsv`; публичный
синтетический контроль реальным sealed runtime подтвердил попадание новых источников
в `gene_context_json`, без изменения gene-validity rank или ACMG/PVS1/LoF/GoF.
Текущий loader связывает gene context по canonical SYMBOL, HGNC source provenance
сохранён явно; per-variant HGNC контроль входит в review нового binding.

**CADD GRCh38 v1.7** получен только для фиксированных публичных CDS, не по VCF.
[Официальный download](https://cadd.bihealth.org/download) ссылается на большой
87,473,403,655-byte BGZF и 2,761,840-byte TBI. Полный BGZF не скачивается. Изолированный
pysam 0.24.1 установлен из [официального PyPI](https://pypi.org/project/pysam/0.24.1/)
по publisher SHA256. HTTPS byte ranges проходят через локальный ограничитель с
реальным учётом upstream bytes; index проверен по publisher MD5. Для baseline CDS
получено 730,833 SNV rows / 135,135,232 HTTP bytes; для TNNT2 — 4,161 rows /
3,276,800 bytes. Каждый REF проверен по подготовленной GRCh38 FASTA; source header
v1.7 и три альтернативы публичного контрольного основания прошли контроль. Canonical
`CADD_RAW`/`CADD_PHRED` TSV совместим с существующим exact-allele local adapter.
Это не all-genome, indel, noncoding или MT CADD.

**MyVariant public reference cache** использует только gene names из фиксированных
публичных HGNC-панелей, без patient coordinates/ALT/VCF/phenotype. Основные 38 gene
queries дали 568,603 records, 775,525,942 raw bytes; TNNT2 — ещё 2,428 records.
[Документация MyVariant](https://docs.myvariant.info/en/latest/doc/data.html) отдельно
предупреждает о hg19 default и source-specific terms. Запросы явно `assembly=hg38`;
API `_id`, structured hg38 mapping, REF/ALT и GRCh38 FASTA проверяются совместно.
Raw records сохраняются независимо от результатов проверки; несовместимые записи не
попадают в нормализованный exact-allele adapter. Metadata: API build 20250624,
dbNSFP **4.8a**, не полный текущий dbNSFP5.4. API CADD subscore version не установлен;
он не заменяет официальный CADD1.7 и новый ClinVar/gnomAD/VEP115.

FATHMM legacy, FATHMM-XF, FATHMM-MKL, REVEL, CADD, MetaLR, VARITY, SIFT/PolyPhen raw
vectors сохраняются раздельно с версиями и page SHA. У публичного MYBPC3 control
FATHMM vectors имеют две оценки, Ensembl transcript array — три entries: их нельзя
слепо zip/flatten или раздать каждой VEP-изоформе. Поэтому raw-vector TSV добавляет
их в `local_evidence_json`; он не создаёт ложный canonical FATHMM scalar. SpliceAI /
MaxEntScan missing остаётся missing; отдельный локальный SpliceAI model runtime —
другой ресурс с собственной GENCODE24 областью и controls.

**REVEL v1.3**: publisher MD5 и SHA256 полного 667,102,707-byte архива прошли проверку.
Единственный CSV внутри — 6,496,023,640 uncompressed bytes. Он потоково фильтруется
по новой public CDS39-панели; целый CSV не распаковывается. Используются только
lifted GRCh38 coordinate, точные REF/ALT и проверка reference base. Unlifted и
REF-mismatch не превращаются в scores. Source ENST IDs не имеют версий; amino-acid
ref/alt и IDs сохраняются. Raw `REVEL_reference_score` не выдаётся за автоматически
подтверждённый score каждой versioned VEP115-изоформы: нужен exact allele + текущая
AA substitution + transcript-base match с явным статусом отсутствия source version.

## Одна команда мониторинга

```bash
python3 /srv/marina/private/reviewed-public-reference-finalize-v1/scripts/public_reference_status.py --once
```

Или watch, максимум 8 часов и не позже исходного VM deadline:

```bash
python3 /srv/marina/private/reviewed-public-reference-finalize-v1/scripts/public_reference_status.py --watch --interval 60 --max-hours 8
```

Новый read-only monitor включает основной VEP/MT/coverage/annotation и все публичные
sources, включая GeneReviews, SpliceAI, RefSeq queue, OMIM identifiers и LOVD.
Он не читает секреты, не хеширует гигабайты, не перезапускает процессы и не меняет
существующие streams. Optional `--report-dir` должен указывать отдельный каталог
`monitor-reports` нового пакета. `--local --root .` читает только прежнюю allowlisted
квитанцию переноса, без SSH/SFTP configuration.

## Новый ROOT review candidate

`merge_public_hcm_candidates.py` объединяет совместимые CADD source TSV с отказом
при conflicting values, добавляет gene bridge и только готовые raw API/REVEL
adapters. Результат — `config-overlay.json` с append-only local_adapters и точными
source/version/bytes/SHA. Он **не применяет** overlay и не меняет approval старой
очереди. ROOT должен создать и проверить новый config/resource-manifest/binding,
сохранив старые locked core resources. Gene context и raw model scores не являются
автоматической clinical classification; KEGG enrichment остаётся withdrawn.
