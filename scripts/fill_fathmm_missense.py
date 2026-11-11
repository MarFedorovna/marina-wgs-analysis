"""Audit and populate missing FATHMM scores for the complete Missense sheet."""
import argparse
import csv
import gzip
import json
import sqlite3
from collections import Counter
from contextlib import closing
from copy import copy
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from fathmm_lookup import load_annotator, score_variant
from reconcile_workbook_vcf import key as variant_key


AUDIT_SHEET = "FATHMM_Missense"
SUMMARY_SHEET = "Полнота_аннотаций"
VARIANT_FIELDS = ("Хромосома", "Позиция", "REF", "ALT", "Ген")
SCORE_FIELDS = ("FATHMM", "FATHMM transcript")
SOURCE_STATUS = {
    "NUMERIC_SCORE": "Есть числовая оценка",
    "NO_RECORD": "Запись в базе не найдена",
    "NO_NUMERIC_SCORE": "Запись есть, числовой оценки нет",
    "AMBIGUOUS": "Неоднозначная запись в базе",
}
MODULE_STATUS = {
    "SCORE_RETURNED": "Оценка возвращена",
    "NO_SCORE": "Оценка не возвращена",
    "ERROR": "Ошибка модуля",
    "MULTIPLE_ROWS": "Несколько строк в базе",
}


def select_missing(input_path, vcf):
    workbook = load_workbook(input_path, read_only=True, data_only=True)
    source = iter(workbook["Missense"].values)
    header = next(source)
    expected = ("Ранг", "Образцы-носители", "Ген", "Хромосома", "Позиция", "REF", "ALT", "HGVS-c", "HGVS-p")
    if header[:len(expected)] != expected or header[11] != "Транскрипт":
        raise ValueError("Unexpected Missense columns")
    score_column = header.index("FATHMM")
    rows = list(source)
    workbook.close()
    by_key = {}
    for row in rows:
        if not all(row[i] not in (None, "", ".") for i in range(3, 7)):
            continue
        key = variant_key(*row[3:7])
        if key in by_key:
            raise ValueError(f"Duplicate allele in Missense: {key}")
        by_key[key] = row
    missing = {key: row[:12] for key, row in by_key.items() if row[score_column] in (None, "")}
    confirmed = set()
    with gzip.open(vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alts, _ = line.split("\t", 5)
            for alt in alts.split(","):
                key = variant_key(chrom, pos, ref, alt)
                if key in missing:
                    confirmed.add(key)
            if len(confirmed) == len(missing):
                break
    if len(confirmed) != len(missing):
        raise ValueError(f"{len(missing) - len(confirmed)} alleles could not be verified in the VCF")
    return len(by_key), list(missing.values())


def scan_database(input_path, vcf, database, code_archive, tsv):
    total, candidates = select_missing(input_path, vcf)
    print(f"Verified VCF keys: {len(candidates)}; unique Missense alleles: {total}", flush=True)
    results = []
    with closing(sqlite3.connect(database.resolve().as_uri() + "?mode=ro&immutable=1", uri=True)) as connection:
        cursor = connection.cursor()
        annotator = load_annotator(cursor, code_archive)
        for row in candidates:
            results.append(score_variant(row, annotator, cursor))
            if len(results) % 250 == 0:
                print(f"Checked {len(results)}/{len(candidates)}", flush=True)
    counts = dict(Counter(r["source_status"] for r in results))
    module = dict(Counter(r["module_status"] for r in results))
    scored = [r for r in results if r["scores"] is not None]
    recovered = sum(r["module_status"] == "ERROR" for r in scored)
    disambiguated = sum(r["module_status"] == "MULTIPLE_ROWS" for r in scored)
    transcript_matches = sum(r["workbook_transcript_score"] is not None for r in scored)
    stats = dict(total=total, already_filled=total-len(results), checked=len(results), added=len(scored),
                 recovered=recovered, disambiguated=disambiguated,
                 transcript_matches=transcript_matches, module=module, source=counts)
    with tsv.open("w", encoding="utf-8-sig", newline="") as handle:
        if results:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(results[0]))
            writer.writeheader()
            writer.writerows(results)
        else:
            handle.write("status\nNO_MISSING_VARIANTS\n")
    print(json.dumps(stats, ensure_ascii=False), flush=True)
    return stats, results


def add_audit(workbook, results):
    sheet = workbook.create_sheet(AUDIT_SHEET)
    sheet.freeze_panes = "A2"
    sheet.sheet_view.showGridLines = False
    sheet.append([
        "Ранг", "Образцы", "Ген", "Хромосома", "Позиция GRCh38", "REF", "ALT", "HGVS-p",
        "Статус в базе", "Статус модуля", "FATHMM по транскриптам", "Транскрипты FATHMM",
        "Транскрипт в книге", "FATHMM для транскрипта в книге", "Примечание к проверке", "Источник оценки",
    ])
    for r in results:
        sheet.append([
            r["rank"], r["samples"], r["gene"], r["chrom"], r["pos"], r["ref"], r["alt"], r["hgvs_p"],
            SOURCE_STATUS[r["source_status"]], MODULE_STATUS[r["module_status"]], r["scores"], r["transcripts"],
            r["workbook_transcript"], r["workbook_transcript_score"], r["error"], r["score_source"],
        ])
    header_fill = PatternFill("solid", fgColor="23425D")
    header_font = Font(color="FFFFFF", bold=True)
    wrap = Alignment(vertical="center", wrap_text=True)
    for cell in sheet[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = wrap
    sheet.row_dimensions[1].height = 35
    alternating_fill = PatternFill("solid", fgColor="F0F4F8")
    score_fill = PatternFill("solid", fgColor="DFEFE5")
    error_fill = PatternFill("solid", fgColor="FCE9DA")
    for i, result in enumerate(results, 2):
        cells = sheet[i]
        for cell in cells:
            cell.alignment = wrap
            if i % 2:
                cell.fill = alternating_fill
        if result["scores"] is not None:
            cells[8].fill = score_fill
        if result["module_status"] == "ERROR":
            cells[9].fill = error_fill
        cells[13].number_format = "0.00"
        lines = max((len(result["transcripts"] or "") + 39) // 40, (len(result["scores"] or "") + 22) // 23)
        sheet.row_dimensions[i].height = min(170, max(29, 14 * lines + 7))
    for col, width in enumerate([8, 20, 14, 13, 17, 6, 6, 20, 35, 23, 27, 48, 25, 25, 49, 43], 1):
        sheet.column_dimensions[get_column_letter(col)].width = width
    sheet.auto_filter.ref = f"A1:P{sheet.max_row}"
    notes = sheet.max_row + 2
    sheet.cell(notes, 1, "Источник")
    sheet.cell(notes, 2, "OpenCRAVAT FATHMM 2.4.1, база FATHMM v2.3").hyperlink = "https://store.opencravat.org/modules/fathmm/2.4.1/fathmm.yml"
    sheet.cell(notes + 1, 2, "Сопоставление транскриптов по стабильному идентификатору Ensembl без суффикса версии.")
    sheet.cell(notes + 2, 2, "Столбец K содержит опубликованные числовые оценки по транскриптам. Столбец N заполнен при совпадении транскрипта.")


def update_summary(workbook, stats):
    sheet = workbook[SUMMARY_SHEET]
    total = stats["total"]
    filled = stats["already_filled"] + stats["added"]
    new = stats["added"]
    checked = stats["checked"]
    source = stats["source"]
    for values in (
        ["Missense с пропуском", "FATHMM v2.3", "вариант", checked, new, checked - new,
         "Проверены все пропуски",
         f"Нет записи в базе: {source.get('NO_RECORD', 0)}. Запись без числовой оценки: "
         f"{source.get('NO_NUMERIC_SCORE', 0)}. Неоднозначные записи: {source.get('AMBIGUOUS', 0)}."],
        ["Весь лист Missense", "FATHMM v2.3", "вариант", total, filled, total - filled,
         "Итог после проверки", f"Ранее заполнено: {stats['already_filled']}. Добавлено: {new}. "
         f"Извлечено из базы после ошибки модуля: {stats['recovered']}. "
         f"Выбрано по транскрипту при нескольких строках базы: {stats['disambiguated']}. Лист FATHMM_Missense."],
    ):
        sheet.append(values)
        for cell in sheet[sheet.max_row]:
            cell._style = copy(sheet.cell(4, cell.column)._style)
            cell.alignment = Alignment(vertical="center", wrap_text=True)
        sheet.row_dimensions[sheet.max_row].height = 86
    if sheet.auto_filter.ref:
        sheet.auto_filter.ref = f"A1:H{sheet.max_row}"


def update_variant_sheets(workbook, results):
    scores = {(r["chrom"], r["pos"], r["ref"], r["alt"], r["gene"]): r for r in results if r["scores"] is not None}
    counts = {}
    for sheet in workbook:
        header = {cell.value: cell.column for cell in sheet[1]}
        keys = ["Хромосома", "Позиция", "REF", "ALT", "Ген", "FATHMM", "FATHMM transcript"]
        if not all(k in header for k in keys):
            continue
        seen = set()
        for row in sheet.iter_rows(min_row=2, max_col=max(header[k] for k in keys[:5])):
            key = tuple(row[header[k] - 1].value for k in keys[:5])
            result = scores.get(key)
            if result is None:
                continue
            cell = sheet.cell(row[0].row, header["FATHMM"])
            raw = result["scores"]
            score = raw if ";" in raw else float(raw)
            if cell.value not in (None, "", score, raw):
                raise ValueError(f"Unexpected existing FATHMM score in {sheet.title}!{cell.coordinate}")
            cell.value = score
            if result["transcripts"]:
                sheet.cell(row[0].row, header["FATHMM transcript"], result["transcripts"])
            seen.add(key)
        if seen:
            counts[sheet.title] = len(seen)
    if any(counts.get(name, 0) != len(scores) for name in ("Missense", "Все_прошедшие")):
        raise ValueError("FATHMM updates did not cover both Missense and the primary autosomal sheet")
    return len(scores), counts


def write_workbook(input_path, output, stats, results):
    workbook = load_workbook(input_path)
    if AUDIT_SHEET in workbook.sheetnames:
        raise ValueError("The input workbook already contains a FATHMM audit")
    count, affected_sheets = update_variant_sheets(workbook, results)
    assert count == stats["added"]
    add_audit(workbook, results)
    update_summary(workbook, stats)
    partial = output.with_suffix(".writing.xlsx")
    workbook.save(partial)
    workbook.close()
    partial.replace(output)
    print(f"Workbook saved. Updated scores across {len(affected_sheets)} existing sheets.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("input", "qc-vcf", "database", "code-archive", "output", "audit-tsv"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    inputs = (args.input, args.qc_vcf, args.database, args.code_archive)
    for path in inputs:
        if not path.is_file():
            parser.error(f"Input file not found: {path}")
    for path in (args.output, args.audit_tsv):
        if path.exists() or path.resolve() in {x.resolve() for x in inputs}:
            parser.error(f"Output path already exists: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
    if args.output.resolve() == args.audit_tsv.resolve():
        parser.error("Workbook and audit TSV paths must differ")
    stats, results = scan_database(args.input, args.qc_vcf, args.database, args.code_archive, args.audit_tsv)
    write_workbook(args.input, args.output, stats, results)


if __name__ == "__main__":
    main()
