# X/Y: только проверяемый план, без запуска calling

```bash
python3 scripts/sex_chromosome_plan.py --config /srv/marina/private/workflow-v1.json
```

Команда читает существующие манифесты/квитанции и показывает блокировки. Она не
меняет BAM, FASTA, MT или завершённые аутосомные результаты. В программе **нет
`--execute`** и нет вызова Docker/GATK: даже при заполненных данных выдаются только
шаблоны команд для следующего отдельного review. Germline X/Y не является
соматическим анализом.

У всех четырёх текущих образцов пол/кариотип не подтверждён. Отношения X/Y из
технического QC не устанавливают XX/XY и не задают плоидность. Поэтому исходные
`karyotype`, `x_nonpar_ploidy`, `y_nonpar_ploidy`, `x_par_ploidy` равны `unknown`.
Для просмотра таблицы/интервалов можно создать новый приватный план:

```bash
python3 scripts/sex_chromosome_plan.py --config /srv/marina/private/workflow-v1.json \
  --write-plan /srv/marina/private/sex-chromosomes-review-v1
```

Это создаёт только `plan.json`, `ploidy-review.tsv`, PAR/nonPAR `.intervals` и BED.
Существующий каталог не перезаписывается. В TSV нельзя заполнять значения по
догадке: нужны явные copy-number/model evidence и его источник. XX/XY здесь —
поддерживаемые модели обработки, не автоматически установленный клинический
кариотип. Мозаицизм/анеуплоидия и другие модели требуют отдельного плана.

## PAR и референс

Официальные [GRC GRCh38 PAR](https://www.ncbi.nlm.nih.gov/grc/human) закреплены в
[lock](../config/sex-chromosomes.resources.lock.json): 1-based closed X PAR1
10001–2781479, PAR2 155701383–156030895; Y PAR1 10001–2781479, PAR2
56887903–57217415. BED явно переводится в 0-based half-open. NonPAR — точное
дополнение на X/Y. Y-PAR не выдаётся второй независимой записью.

Простой diploid chrX-only calling на BAM, выровненном на неподтверждённый/unmasked
X+Y reference, **не считается проверенным PAR handling**. Broad описывает
[Y-PAR hard masking в analysis-set reference](https://raw.githubusercontent.com/broadinstitute/gatk-docs/master/gatk3-dictionary/Reference_Genome_Components.md).
Нужно проверить реальную последовательность и то, что этот mask присутствовал при
исходном выравнивании; название GRCh38 само этого не доказывает.

```bash
python3 scripts/sex_chromosome_plan.py --config /srv/marina/private/workflow-v1.json --preflight
```

`--preflight` только читает FASTA, проверяет SHA/FAI и считает N/ACGT в Y-PAR;
варианты не вызывает. Сам нынешний mask ещё не доказывает историю BAM. Для снятия
PAR-блокировки дополнительно нужен отдельно проверенный
[alignment contract](../config/sex-chromosome-alignment-contract.example.json),
связанный с точным FASTA и SHA исходных alignment receipts каждого образца.
Пример не является PASS. Если требуются новый masked reference и remapping,
сначала отдельно подготовить и проверить их из исходных ридов; не редактировать
старую FASTA/BAM и не просто убрать Y-PAR из интервалов.

## Что пока не выполнено

Положительные плоидности должны быть заданы отдельно для каждой части; для
поддерживаемых XX/XY требуется явное согласованное заполнение, без diploid default.
Y=0 означает **skip/not assessed, без GT**, а не `0/0`, «генотип референсный» или
доказательство отсутствия Y. Unknown, низкое покрытие и zero-copy различаются.

Шаблоны используют существующий immutable GATK4.3 image, HaplotypeCaller с явным
`--sample-ploidy`, CPU4/16GiB без GPU. Они содержат `SET_*` и не исполняются.
Будущему executor нужно заново проверить исходный абсолютный бюджет $40, BAM
identity, выходные каталоги, gVCF/плоидность и смешанное joint-genotyping. Здесь
не реализованы X/Y calling, cohort merge, анеуплоидия/CNV или клинический вывод.
MT Mutect2 и текущие серверные процессы этой программой не затрагиваются.
