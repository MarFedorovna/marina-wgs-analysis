# Количественная проверка X/Y

`scripts/xy_sparse_quantitative_qc.py` читает выбранные области четырёх индексированных BAM. Итоговый VCF, плоидность и другие расчёты не меняются. Это технические измерения; пол и кариотип не назначаются.

На X вне PAR используется независимая публичная панель распространённых двухаллельных SNP из проверенного VEP115 cache: глобальная AF 1000 Genomes phase3 0,1–0,9, точное совпадение REF с GRCh38. Частоты четырёх исследуемых образцов для отбора панели не используются. Для каждого сайта сохраняются REF/ALT/прочие основания, глубина, баланс аллелей и причина исключения.

Фильтры: MQ≥30, BQ/BAQ≥20; исключены unmapped, secondary, supplementary, duplicate и QC-fail. Перекрытие пар учитывается один раз. Классифицируются сайты с ACGT DP 10–100 и долей прочих наблюдений ≤5%; HET требует AB 0,25–0,75 и ≥3 наблюдений каждой аллели. Промежуточный AB остаётся неопределённым. Это классификация счётчиков, а не полноценный вызов вариантов.

`het_fraction` использует только классифицированные сайты. `F_descriptive = 1 − Hobs/Hexp` выводится лишь при наличии внешней AF, ≥100 классифицированных сайтов и ≥10 покрытых 10-Мб интервалов. Это прокси по глобальной популяции: происхождение, LD и отбор покрытых сайтов могут смещать результат. Интервал Wilson — приближение по числу сайтов, не диагностическая уверенность. При недостатке данных вывод остаётся неопределённым.

Для X, Y и аутосом измеряются распределённые окна по 1000 ACGT-оснований: средняя глубина и доля оснований с глубиной ≥1/5/10. Из X/Y исключены PAR; окна с N не включаются. Это охват выбранных окон, не всей Y. Высокий MQ не заменяет маску уникальной картируемости. Нулевой сигнал в разреженной панели не доказывает отсутствие Y; неравномерность аутосомных окон показывается отдельно.

Приватная подготовка не читает выравнивания BAM и не вычисляет их новые хеши:

```bash
python3 scripts/xy_sparse_quantitative_qc.py prepare \
  --config /PRIVATE/technical-xy-qc-config.json \
  --panel /PUBLIC_PANEL/x_nonpar_common_snps.tsv \
  --panel-manifest /PUBLIC_PANEL/validated-contract.json \
  --output /NEW_QC_OUTPUT --profile /PRIVATE/candidate-profile.json \
  --deadline-utc EXACT_SEPARATE_QC_DEADLINE_UTC
```

Проверенное правило запускается отдельно. Без `--execute` остаётся только просмотр плана:

```bash
python3 scripts/xy_sparse_quantitative_qc.py run \
  --profile /PRIVATE/ROOT-reviewed-profile.json \
  --reviewed-profile-sha256 EXACT_SHA256 --execute
```

Пределы: четыре образца одновременно, до 16 CPU суммарно, память ≤24 ГиБ, не более 40 минут и отдельного абсолютного срока. Чтение выбирается через индекс `samtools view -M -L`, затем ограничивается теми же позициями в pileup. Полный BAM не сканируется и не хешируется. Результаты: `quantitative-xy-profile.json`, четыре таблицы `sample_000N.x-site-counts.tsv`, метрики окон, журналы и квитанция.

Источники: [samtools pileup и фильтры](https://www.htslib.org/doc/samtools-mpileup.html), [индексированный multi-region iterator](https://www.htslib.org/doc/samtools-view.html), [GRCh38 PAR](https://useast.ensembl.org/info/genome/genebuild/human_PARS.html), [ограничения гетерозиготности/F](https://www.cog-genomics.org/plink/2.0/basic_stats).
