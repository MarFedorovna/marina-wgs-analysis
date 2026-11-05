# Одобренные методики P03/P06/P10

Версия 1, 2026-11-05 UTC. Пользователь одобрил P01–P11. Этот пакет реализует
отдельные методические инструменты P03/P06/P10; текущие sealed runtime, calling,
annotation controllers и задания stages6/7 не изменяются. Новые пустые шаблоны
не означают, что пациентские варианты просмотрены или классифицированы.

## MT evidence worksheet

[Шаблон MT](../templates/mt-acmg-evidence-v1.json) фиксирует rCRS, тип MT-гена,
материал, долю альтернативного аллеля, haplogroup, NUMT/contamination flags,
популяционные знаменатели и независимые публикации. Каждый критерий добавляется
отдельной записью с версией guidance, силой, основанием и противоречиями. Исходное
состояние — unknown, без автоматического ACMG class.

Основной первичный источник — [MT specifications, PMID32906214](https://pmc.ncbi.nlm.nih.gov/articles/PMC7717623/).
Перед применением проверяют актуальную [ClinGen guidance](https://clinicalgenome.org/affiliation/50027/)
и уже скачанный CSpec registry. Не переносить nuclear NMD/splicing-схему на mtDNA,
protein missense predictors на tRNA/rRNA или общие учебные пороги на конкретный
тест. Технический MT PASS не подтверждает клиническую чувствительность. Материнские
ткани/assay limits, происхождение haplogroup counts и независимость probands могут
оставаться unknown. Haplochecker warning не равен доказанной контаминации.

## Offline HGVS proof

[Setup](../scripts/setup_offline_hgvs_public_controls.py) устанавливает Mutalyzer3.1.1
в отдельный venv из официального PyPI. [Hash lock](../config/offline-hgvs-mutalyzer3.1.1-py310.lock.txt)
содержит 23 точные версии и SHA256 source artifacts; pip25.3, setuptools80.9.0,
wheel0.45.1 фиксированы отдельно. Часть upstream packages распространяется только
как sdist; нужен C/C++ build toolchain. Это изолированный Python3.10 Linux lock,
не дополнение к SpliceAI/TensorFlow/VEP среде.

Официальные [Mutalyzer docs](https://mutalyzer.readthedocs.io/en/latest/usage.html)
и [source MIT](https://github.com/mutalyzer/mutalyzer) проверены. Публичные FASTA/GFF3
fixtures закреплены commit171f36ada65b228c443b359925305bd85e9fa40e; их SHA сохраняется
в source manifest. Использованы NM_003002.4 и NC_012920.1 и два публичных versioned
NCBI HCM references, без patient inputs, hosted HGVS запросов, email, API key или аккаунта.

[Runner](../scripts/run_offline_hgvs.py) применяет тот же local retrieve_raw boundary,
что и upstream fixture tests, плюс запрет socket/DNS. Проверка дала **6/6 PASS**:
публичный MT SNV, transcript deletion, неправильная REF-база отклонена, отсутствующий
reference явно отклонён, синтетические замены start-codon для двух публичных HCM
references корректно нормализованы. Guard self-test PASS, попыток сети normalizer=0.

Из существующей MANE1.4 summary подтверждены MYBPC3 NM_000256.3 ↔ ENST00000545968.6
и MYH7 NM_000257.4 ↔ ENST00000355349.4. [Offline concordance check](../scripts/verify_hgvs_public_reference_concordance.py)
нашёл точные ENST versions, gene и strand в фактическом GENCODE49: **2/2 PASS**.
Это проверка идентификаторов, не genomic-to-transcript conversion каждого варианта.
[Публичные контрольные квитанции](../config/offline-hgvs-public-proof-v1.json) сохраняют
точные reference URLs/digests и границы выполненной проверки.

Это доказательство локальной работоспособности HGVS normalizer. Оно **не сертифицирует**
актуальные HCM transcripts MANE1.4/GENCODE49/VEP115 и не доказывает genomic↔transcript
mapping для каждого кандидата. Для применения к shortlist нужно добавить только
публичные versioned references, проверить их связь с текущим MANE и cache115,
дополнить source manifest и прогнать matching controls. При отсутствии reference
возвращается unavailable, не скачивание по пациентскому вводу. RefSeq VEP cache
и Mutalyzer reference fixtures — разные форматы.

Из нового пустого каталога:

```bash
python3 scripts/setup_offline_hgvs_public_controls.py --root /isolated/offline-hgvs --mane-summary "$PUBLIC_MANE14_SUMMARY"
/isolated/offline-hgvs/venv/bin/python scripts/run_offline_hgvs.py --root /isolated/offline-hgvs --public-controls
python3 scripts/verify_hgvs_public_reference_concordance.py --hgvs-root /isolated/offline-hgvs --gencode-gtf "$PUBLIC_GENCODE49_GTF"
```

Рабочий расширенный server proof расположен отдельно в resources/offline-hgvs-public-controls-v4;
v1/v2 были неуспешными setup attempts (старый pip без --report; отсутствие sdist
при первоначальном wheels-only режиме). Они не считаются успешными установками.

## IGV manual review

[Шаблон IGV](../templates/igv-read-review-v1.json) сохраняет reference/contig match,
индекс BAM, версию transcript track и настройки отображения. Проверяются независимые
фрагменты, обе цепи, read ends, baseQ/mapQ, duplicates, read groups, softclips,
nearby indel, повторы/паралоги, ориентация пар и MT NUMT context. Для indel показывают
основания и вставки/делеции; сохраняют влияние downsampling и фильтров. Видимые
read counts могут отличаться от FORMAT/AD из-за алгоритма и display filters.

[IGV guide](https://igv.org/doc/desktop/UserGuide/tracks/alignments/viewing_alignments_basics/)
описывает BAM/BAI, настройки качества и отображения. [IGV Desktop](https://igv.org/doc/desktop/)
распространяется под MIT. На headless server IGV/java не обнаружены в PATH;
GUI не устанавливался, patient screenshot/session не создавались. Полные BAM/BAI
и FASTA можно просматривать локально после их получения; внешняя загрузка не нужна.

Итог — reviewed_support/reviewed_concern/inconclusive/unassessed с причиной.
Он не присваивает pathogenicity, не подтверждает sample identity по одному локусу
и не заменяет необходимую orthogonal validation.
