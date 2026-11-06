# Отдельный germline calling X/Y

Новый `scripts/run_xy_germline.py` не изменяет аутосомный или MT-анализ. Исполнитель требует проверенный приватный профиль и точный SHA. Без `--execute` выполняется только проверка входов и плана.

```bash
python3 scripts/run_xy_germline.py --profile "$XY_PROFILE" --profile-sha256 "$XY_PROFILE_SHA256"
python3 scripts/run_xy_germline.py --profile "$XY_PROFILE" --profile-sha256 "$XY_PROFILE_SHA256" --execute
```

Это отдельная явно разрешённая пользователем модель обработки по техническому QC, а не клинически установленный пол или кариотип. Четыре образца задаются в приватном профиле: X-nonPAR плоидность 1/1/2/2, X-PAR 2/2/2/2, Y-nonPAR 1/1/skip/skip. Пропущенный Y остаётся неоценённым с отсутствующим GT, а не референсным 0/0. PAR выдаётся только в координатах X после доказанного Y-PAR masking при исходном выравнивании.

Профиль schema_version=1/status=ROOT_REVIEWED_XY_GERMLINE_V1 содержит processing_model=USER_AUTHORIZED_INFERRED_FROM_QC_NOT_CLINICAL_KARYOTYPE; samples (sample, X_NONPAR, X_PAR, Y_NONPAR, bam, bai, alignment_receipt); reference/reference_fai/reference_dict; qc_receipt; reference_mask_proof; alignment_contract; par_lock; executor; filter_module; bcftools; samtools; immutable gatk_image; public_control; output; deadline_utc; min_free_gib. Малые pins имеют path/bytes/sha256; BAM и FASTA — также mtime_ns, без повторного хеширования гигабайт. Исторический BAM SHA сверяется с точной исходной alignment receipt. Contract связывает все receipt SHA, точный reference SHA, официальный PAR GFF SHA, mask_present_at_original_alignment=true/par_policy=Y_PAR_HARD_MASKED_AT_ALIGNMENT/par_report_coordinates=X_only и разрешённую модель обработки.

Ресурсы фиксированы: resources={call_workers:6,call_cpus:4,call_memory_gib:16,max_cpu:24,max_memory_gib:120}. Десять отдельных GATK4.3 HaplotypeCaller дают gVCF с явной плоидностью. Для каждой из трёх частей CombineGVCFs и GenotypeGVCFs используют плоидность исходных генотипов без принудительной общей diploid-настройки. GATK запускается без сети, в закреплённом digest CPU image, с UID/GID владельца. Команды записываются только в приватный журнал.

Оригинальные joint.raw.vcf.gz и индексы сохраняются. Производный результат проходит `bcftools norm -m -any --multi-overlaps . -f REF -c e`: удалённый альтернативный аллель становится missing, а не REF. Поэтому производные частичные GT сохраняются и помечаются существующим консервативным soft_filter_vcf; все исходные GT остаются в raw-артефактах. Отсутствие DP/GQ/аннотаций не считается хорошим качеством. AD/PL проецирует bcftools; сумма AD не переопределяется. Итоговый cohort.xy.soft-filtered.vcf.gz содержит все четыре колонки, индекс TBI и отдельные X/Y scope receipts. Счётчик record_count относится к разложенным biallelic-строкам.

Перед пациентскими входами обязателен реальный отдельный синтетический контроль на точном image и bcftools: HaplotypeCaller для haploid ALT и diploid HET, смешанные CombineGVCFs/GenotypeGVCFs, затем реальная multi-ALT decomposition с missing вместо удалённого ALT. Контроль не читает BAM проекта:

```bash
python3 scripts/run_xy_germline.py --public-control --gatk-image "$GATK_CPU_IMAGE_DIGEST" --control-output "$NEW_PUBLIC_CONTROL_DIR" --execute
```

Используются `samtools` и новый проверенный `bcftools` из PATH; proof связывает SHA обоих бинарных файлов. Само наличие синтетического VCF не заменяет выполненный контроль инструментов.

Абсолютный срок — не позднее 2026-11-06T06:20:00Z и не более двух часов от запуска. Это отдельно разрешённое окно X/Y, не изменение старого срока анализа. Контроллер имеет flock и общий таймер; subprocess ограничены тем же сроком. При ошибке сохраняются уникальные attempt-каталоги. Переиспользуются только завершённые задачи с точным profile identity и проверенными артефактами; чужие выходы не заменяются. После всех частей повторно проверяются финальные GT/плоидность, порядок, число строк, отсутствие повторных аллелей и границы частей. PASS calling не означает завершённую VEP-аннотацию или клиническую интерпретацию; последующая аннотация выполняется отдельно.

Основания: [GATK о смешанной плоидности](https://gatk.broadinstitute.org/hc/en-us/articles/360035889691-Does-GATK-work-on-non-diploid-organisms), [bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm). Мозаицизм, CNV/анеуплоидии, X-transposed/повторные области и недостаточное покрытие этим запуском не диагностируются.

## Воспроизводимая установка bcftools 1.22

Этот запуск использует отдельный префикс, не заменяя системный bcftools. Официальный [release 1.22](https://github.com/samtools/bcftools/releases/tag/1.22) содержит HTSlib 1.22 в архиве. Закреплённый архив — 8 176 878 байт, SHA-256 `f2ab9e2f605b1203a7e9cbfb0a3eb7689322297f8c34b45dc5237fe57d98489f` (зафиксирован при фактической загрузке; это не подпись издателя).

```bash
PREFIX="$HOME/.local/tools/bcftools-1.22"
mkdir -p "$PREFIX/source"
cd "$PREFIX/source"
curl --fail --location --output bcftools-1.22.tar.bz2 \
  https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2
printf '%s  %s\n' \
  f2ab9e2f605b1203a7e9cbfb0a3eb7689322297f8c34b45dc5237fe57d98489f \
  bcftools-1.22.tar.bz2 | sha256sum --check -
tar -xjf bcftools-1.22.tar.bz2
cd bcftools-1.22
./configure --prefix="$PREFIX/installed" --disable-plugins \
  --disable-libcurl --disable-bz2 --disable-lzma --disable-ref-cache
make -j8
make install
"$PREFIX/installed/bin/bcftools" --version
sha256sum "$PREFIX/installed/bin/bcftools"
```

Нужны C compiler, make, zlib development headers и стандартные POSIX tools. Сборка ограничивалась 8 CPU и 600 секундами; фактическая успешная сборка заняла 17 секунд. Отключённые curl/ref-cache/cloud-возможности не нужны для локального VCF/FASTA. Системные программы не заменялись. Фактически собранный бинарный файл — 9 394 608 байт, SHA-256 `612798aeab1ac5a4697f68972b3fc1a895f6b55dc4648115c6076ccf70d450b3`; размер и SHA другой сборки могут отличаться, поэтому новый бинарник заново фиксируется в профиле и проходит контроль.

Перед синтетическим контролем добавьте проверенный новый `installed/bin` в PATH. Затем профиль calling должен содержать **точно тот же** путь, размер и SHA bcftools, которые записаны в успешном public-control.json. Сборка сама по себе не разрешает запуск на образцах и не заменяет проверку плоидности/генотипов.
