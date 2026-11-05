# Параллельная CNV/SV feasibility

По новому поручению P12 имеет статус **PARALLEL_RESEARCH_FEASIBILITY**. Выполнены
официальный public-tool smoke proof и read-only fitcheck новых полных BAM в двух
фиксированных публичных gene windows. **Patient SV/CNV calling не запускался.**
Основная адаптация имеющейся annotation сохраняет приоритет. Прежние VCF/кандидаты
не использовались, BAM subset не создавался.

Сильное обогащение CDS при низком широком геномном покрытии означает, что нельзя
объявлять данные равномерными WGS по имени файла. Это также не устанавливает
конкретную capture-платформу. Для CNV важна стабильность относительного покрытия
сопоставимых targets; для SV — прочтения у обеих breakpoints и их mates. Высокое
покрытие CDS не устанавливает пригодность для genome-wide SV и не позволяет
исключать структурные события. Индивидуальные результаты fitcheck сохраняются
в private receipt, а не в публичном репозитории.

## Проверенное ПО и публичный пример

[DELLY2.7.0](https://github.com/dellytools/delly/releases/tag/v2.7.0) — maintained
официальный release 2026-10-05, BSD-3-Clause. Linux-amd64 binary9871504 bytes
проверен по publisher SHA256 из GitHub asset metadata. Public example закреплён
commit3d505adc17afccb4443cf97dd01f78afc1d5a267; FASTA/BAM/index/license сохранены
с собственными SHA256. Общая передача **11,458,448 bytes**, в пределах50MB.
Версия подтверждена исполняемым файлом. `delly sr` на upstream short-read example
завершился успешно, output BCF читается и содержит одну запись. Это **CLI fixture
smoke PASS**, не golden-truth accuracy benchmark и не измерение чувствительности.

[Официальный README](https://github.com/dellytools/delly) описывает indexed,
duplicate-marked BAM, short-read `sr` и отдельный whole-genome read-depth `cnv`.
Последний не поддерживает exome/targeted данные. Поэтому `delly cnv` не запускался
и не предлагается для автоматического запуска при таком профиле покрытия.

Исполняемые команды реализации:

```bash
python3 scripts/cnv_sv_public_proof_and_fitcheck.py --public-proof-root /isolated/new-delly-public-proof
python3 scripts/cnv_sv_public_proof_and_fitcheck.py --fitcheck-config "$EXISTING_NEW_BAM_COVERAGE_CONFIG" --private-output /private/new-readonly-fitcheck
```

[Script](../scripts/cnv_sv_public_proof_and_fitcheck.py) содержит фиксированные
официальные URLs, checksum, deadline и explicit limits. Public proof использует
два CPU/4GiB и timeout180s для малой команды. [Публичная квитанция](../config/delly270-public-proof-v1.json)
не содержит patient data. [План](../config/cnv-sv-feasibility-plan-v1.json) фиксирует
оставшиеся ограничения и факт отсутствия patient caller.

## Read-only пригодность полных BAM

Проверены четыре source alignment receipts, BAM quickcheck и неизменность размера/
mtime исходных файлов. Из полных indexed BAM прочитаны regions MYBPC3 и MYH7
(MANE1.4 gene body +50kb каждого фланка), с агрегацией read flags/mapQ/softclips/SA,
positive template lengths и depth/breadth. Общее время fitcheck около3.5s.

`samtools view` читает региональные records из полного BAM. `samtools depth -aa -s
-q 0 -Q 20 -G 3844` считает локальный профиль с фиксированными filters; это отдельная
метрика с собственными ограничениями, не новая оценка всего генома. Не создавались
mate-subset BAM и reference background для caller.

Результаты подтверждают выраженную неоднородность этих публичных gene+flank окон.
Это основание оставить whole-genome sensitivity и необозримые breakpoint regions
неоценёнными. Отсутствие найденного события нельзя считать отрицательным SV-тестом.

## Следующий допустимый этап

Для запуска patient caller требуется отдельный meaningful full-BAM или проверенный
target-aware дизайн и root review. На mate-subset BAM SV calling запрещён:
ограниченная выборка искажает background/insert-size context и теряет часть
supplementary/remote breakpoint support. Подготовленная ранее идея такого запуска
отклонена до исполнения. Автоматического расширения на четыре образца нет.

[CNVkit](https://cnvkit.readthedocs.io/en/master/quickstart.html) допускает germline
targeted данные и flat reference, но реальный capture BED и matching normal
reference нужны для оценки систематических biases. GENCODE49 CDS BED не является
capture-design BED; четыре исследуемых образца не становятся normal panel
автоматически. Пока подходящий дизайн не установлен, targeted CNV calling не
выполняется. MT deletions, repeat expansions, methylation и genome-wide CNV остаются
непроверенными отдельными задачами. VM не продлевается; исходный deadline и бюджет
сохраняются, stage6/7 имеет приоритет.
