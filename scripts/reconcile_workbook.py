"""Fill supported annotation gaps in an existing research workbook and verify the result."""
import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from copy import copy
from itertools import zip_longest
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from reconcile_workbook_sources import (
    ALLELE_FIELDS, GENE_FIELDS, MODEL_GROUPS, SPLICE,
    clean, collect_evidence, json_value, key, number, stable,
)

COVERAGE = "Полнота_по_столбцам"
PROVENANCE = "Источники_дополнений"
COORDS = ("Хромосома", "Позиция", "REF", "ALT")
MAIN_SHEETS = {"Все_прошедшие", "X_Y_все", "MT_все"}
CLINVAR = ("ClinVar", "ClinVar review", "ClinVar заболевание", "ClinVar ID")
GNOMAD = ("gnomAD4 AF", "gnomAD4 AC", "gnomAD4 AN")
GTEX = ("GTEx gene", "GTEx tissue")
MITOMAP = ("MITOMAP заболевание", "MITOMAP статус", "MITOMAP score", "MITOMAP PMID")
GENE_GROUP = ("ClinGen заболевание", "ClinGen классификация")
SCORE_COLUMNS = {field for group in MODEL_GROUPS for field in group} | set(SPLICE)
COPIED = set(ALLELE_FIELDS) | set(GENE_FIELDS) | SCORE_COLUMNS | {
    "gnomAD AF max", "ССЗ-контекст ClinGen", "Категории по всем транскриптам",
}
ANNOTATED_HEADERS = {"Образец метрик", "AF вызова"}
CVD_TERMS = re.compile(
    r"cardio|cardiac|heart|aort|vascular|arterial|artery|arrhyth|long.qt|short.qt|brugada|"
    r"ventricul|atrial|hypercholesterolemia|lipoprotein|marfan|sudden cardiac", re.I
)
BAM_COUNTS = {
    "REF шаблоны": "templates_REF", "ALT шаблоны": "templates_ALT", "Другие шаблоны": "templates_OTHER",
    "ALT чтения +": "reads_ALT_forward", "ALT чтения -": "reads_ALT_reverse",
}


def identical(left, right):
    if left == right:
        return True
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return abs(left - right) <= 1e-14 * max(1, abs(left), abs(right))
    return False


def compatible(field, current, candidate):
    if identical(current, candidate):
        return True
    if field.endswith(" transcript"):
        return {stable(x.strip()) for x in str(current).split(";")} == {
            stable(x.strip()) for x in str(candidate).split(";")
        }
    if field in {"ClinVar", "ClinVar review"}:
        normalize = lambda value: re.sub(r"[\s_]+", " ", str(value)).strip().casefold()
        return normalize(current) == normalize(candidate)
    return False


def distinct(values):
    data = {}
    for v in values:
        if clean(v) is not None:
            data[json.dumps(v, sort_keys=True)] = v
    return list(data.values())


def one(values):
    choices = distinct(values)
    return choices[0] if len(choices) == 1 else None


def ensembl_csq(row):
    return {
        "SYMBOL": row.get("gene"), "Feature": row.get("transcript"), "Feature_type": "Transcript",
        "Consequence": row.get("consequence"), "HGVSc": row.get("hgvs_c"), "HGVSp": row.get("hgvs_p"),
        "MANE_SELECT": row.get("mane_select"), "CANONICAL": row.get("canonical"),
        "SIFT": row.get("SIFT"), "PolyPhen": row.get("PolyPhen"),
    }


def choose_effect(annotations, gene, transcript):
    candidates = [
        x for x in annotations if x.get("Feature_type") == "Transcript" and clean(x.get("SYMBOL")) and clean(x.get("Feature"))
    ]
    if gene:
        candidates = [x for x in candidates if x["SYMBOL"] == gene]
    elif len({x["SYMBOL"] for x in candidates}) != 1:
        return []
    if transcript:
        exact = [x for x in candidates if x["Feature"] == transcript]
        if exact:
            return exact
        return [x for x in candidates if stable(x["Feature"]) == stable(transcript) and
                ("." not in str(x["Feature"]) or "." not in str(transcript))]
    if len({x["Feature"] for x in candidates}) == 1:
        return candidates
    for field in ("MANE_SELECT", "CANONICAL"):
        filtered = [x for x in candidates if clean(x.get(field)) and str(x[field]).upper() not in {"NO", "FALSE", "0"}]
        if len({x["Feature"] for x in filtered}) == 1:
            return filtered
    return []


def hgvs(value):
    return str(value).split(":", 1)[-1] if clean(value) is not None else None


def rsids(values):
    matches = {rs for v in values if v for rs in re.findall(r"\brs[0-9]+\b", str(v))}
    return "; ".join(sorted(matches, key=lambda x: int(x[2:]))) if matches else None


def parsed_prediction(value):
    if not isinstance(value, str):
        return None, None
    match = re.fullmatch(r"([^()]+)\(([-+0-9.eE]+)\)", value)
    return (number(match[2]), match[1]) if match else (None, None)


def get_model_rows(rec, selected):
    s = stable(selected)
    sift = [r for r in json_value(rec.get("sift__all"), []) if len(r) >= 3 and stable(r[0]) == s]
    revel = [r for r in json_value(rec.get("revel__all"), []) if len(r) >= 2 and stable(r[0]) == s]
    pp = [r for r in json_value(rec.get("polyphen2__all"), []) if len(r) >= 8 and stable(r[0]) == s and clean(r[7])]
    result = {}
    if len(sift) == 1 and number(sift[0][1]) is not None:
        result["SIFT"] = {"SIFT": number(sift[0][1]), "SIFT prediction": sift[0][2], "SIFT transcript": sift[0][0]}
    if len(revel) == 1 and number(revel[0][1]) is not None:
        result["REVEL"] = {"REVEL": number(revel[0][1]), "REVEL transcript": revel[0][0]}
    if len(pp) == 1:
        result["PolyPhen"] = {"PolyPhen": pp[0][7], "PolyPhen transcript": pp[0][0]}
    return result


def workbook_allele_key(row):
    return key(*(row[field] for field in COORDS))


class Proposal:
    def __init__(self, row):
        self.values = dict(row)
        self.changes = {}

    def put(self, field, value, source):
        value = clean(value)
        if field in self.values and clean(self.values[field]) is None and value is not None:
            self.values[field] = value
            self.changes[field] = (value, source)

    def group(self, values, source):
        values = {f: v for f, v in values.items() if f in self.values and clean(v) is not None}
        present = {f: self.values[f] for f in values if clean(self.values[f]) is not None}
        if any(not compatible(f, old, values[f]) for f, old in present.items()):
            return
        for field, value in values.items():
            self.put(field, value, source)


class Evidence:
    def __init__(self, sources, paths, bam_sources, refseq_locations):
        self.refseq_locations = refseq_locations
        with sources.open(encoding="utf-8") as handle:
            raw = json.load(handle)
        for name, table in raw.items():
            setattr(self, name, {(tuple(k) if isinstance(k, list) else k): v for k, v in table})
        self.xy_fathmm_by_gene = defaultdict(list)
        for k, score in self.xy_fathmm.items():
            self.xy_fathmm_by_gene[k[:-1]].append(score)
        self.hpo = defaultdict(set)
        self.entrez = defaultdict(set)
        with paths["hpo"].open(encoding="utf-8", newline="") as handle:
            for r in csv.DictReader(handle, delimiter="\t"):
                self.hpo[r["gene_symbol"]].add(r["hpo_name"])
                self.entrez[r["gene_symbol"]].add(r["ncbi_gene_id"])
        self.clingen = defaultdict(list)
        with paths["clingen"].open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            for header in reader:
                if header and header[0] == "GENE SYMBOL":
                    break
            for fields in reader:
                r = dict(zip(header, fields))
                if re.fullmatch(r"HGNC:\d+", r.get("GENE ID (HGNC)", "")) and clean(r.get("CLASSIFICATION")):
                    self.clingen[r["GENE SYMBOL"]].append(r)
        self.gene_chrom = defaultdict(set)
        for k, entry in self.auto.items():
            for csq in entry.get("csq", []):
                if csq.get("Feature_type") == "Transcript" and csq.get("SYMBOL"):
                    self.gene_chrom[csq["SYMBOL"]].add("chr" + k[0])
        for dataset in (self.xy, self.mt):
            for k, entries in dataset.items():
                for r in entries:
                    if clean(r.get("gene")) and clean(r.get("transcript")):
                        self.gene_chrom[r["gene"]].add("chr" + ("M" if k[0] == "MT" else k[0]))
        self.bam = defaultdict(list)
        for bam_file in bam_sources:
            with bam_file.open(encoding="utf-8", newline="") as handle:
                for r in csv.DictReader(handle, delimiter="\t"):
                    k = key(*(r[x] for x in ("chrom", "pos", "ref", "alt")))
                    self.bam[k + (r["sample"], r["gene"])].append(r)

    def propose(self, row, inherited):
        k = workbook_allele_key(row)
        p = Proposal(row)
        auto = self.auto.get(k, {})
        xy_rows = self.xy.get(k, [])
        mt_rows = self.mt.get(k, [])
        jobs = self.oc.get(k, [])
        call = self.qc.get(k, {})
        for field, value in inherited.items():
            p.put(field, value, "Та же запись в исходной книге")

        if k[0] == "MT":
            self.mt_proposals(p, mt_rows, jobs)
            return p

        source_csq = auto.get("csq", []) + [ensembl_csq(x) for x in xy_rows]
        for x in xy_rows:
            item = json_value(x.get("vep_all_csq_json"), {})
            source_csq.extend(y for y in ([item] if isinstance(item, dict) else item) if isinstance(y, dict))
        effects = choose_effect(source_csq, clean(p.values.get("Ген")), clean(p.values.get("Транскрипт")))
        effect_source = "VEP 115, тот же аллель и транскрипт"
        if effects:
            symbol = one(x.get("SYMBOL") for x in effects)
            selected = one(x.get("Feature") for x in effects)
            p.put("Ген", symbol, effect_source)
            p.put("Транскрипт", selected, effect_source)
            p.put("Consequence выбранного транскрипта", one(x.get("Consequence") for x in effects), effect_source)
            p.put("HGVS-c", one(hgvs(x.get("HGVSc")) for x in effects), effect_source)
            p.put("HGVS-p", one(hgvs(x.get("HGVSp")) for x in effects), effect_source)

        categories = {term for x in source_csq for term in str(x.get("Consequence") or "").split("&") if term}
        p.put("Категории по всем транскриптам", "; ".join(sorted(categories)), "VEP 115, все эффекты аллеля")
        transcript = clean(p.values.get("Транскрипт"))
        gene = clean(p.values.get("Ген"))
        sample = clean(p.values.get("Образец метрик"))
        self.qc_proposals(p, call, sample)
        self.opencravat_proposals(p, jobs, gene, transcript)
        self.gene_proposals(p, gene)

        for entry in effects:
            value, label = parsed_prediction(entry.get("SIFT"))
            if value is not None:
                p.group({"SIFT": value, "SIFT prediction": label, "SIFT transcript": entry.get("Feature")}, effect_source)
            if clean(entry.get("PolyPhen")):
                p.group({"PolyPhen": entry["PolyPhen"], "PolyPhen transcript": entry.get("Feature")}, effect_source)

        wgs = [x for x in auto.get("transcript", []) if transcript and x["transcript"] == transcript]
        wgs.extend(x for x in xy_rows if transcript and x.get("transcript") == transcript)
        if wgs:
            p.put("CADD PHRED", one(number(x.get("CADD_PHRED")) for x in wgs), "Аннотированный VCF/TSV, точный транскрипт")
            fathmm = one(number(x.get("FATHMM")) for x in wgs)
            if fathmm is not None and "missense_variant" in categories:
                p.group({"FATHMM": fathmm, "FATHMM transcript": transcript}, "Исходная аннотация FATHMM, точный транскрипт")
            alpha = one(number(x.get("AlphaMissense_score")) for x in wgs if str(x.get("AlphaMissense_match_status", "")).startswith("exact_"))
            if alpha is not None:
                p.group({"AlphaMissense": alpha, "AlphaMissense transcript": transcript}, "AlphaMissense, точный транскрипт из VCF/TSV")

        xy_fathmm = self.xy_fathmm.get(k + (gene, stable(transcript)))
        if not xy_fathmm:
            alternatives = self.xy_fathmm_by_gene.get(k + (gene,), [])
            if alternatives and len({(x["scores"], x["transcripts"]) for x in alternatives}) == 1:
                xy_fathmm = alternatives[0]
        if xy_fathmm:
            self.fathmm_proposal(p, xy_fathmm)

        af = auto.get("frequencies", []) + [number(v) for x in auto.get("csq", []) for f, v in x.items()
                                         if f.startswith(("gnomADe_", "gnomADg_"))]
        af += [number(x.get("af_max_observed")) for x in xy_rows]
        af += [number(v) for _, rec in jobs for f, v in rec.items() if f == "gnomad4__af" or f.startswith("gnomad4__af_")]
        p.put("gnomAD AF max", max((v for v in af if v is not None and 0 <= v <= 1), default=None), "Максимум доступных AF того же аллеля")
        p.put("dbSNP", rsids([auto.get("id"), call.get("id"), *[x.get("Existing_variation") for x in source_csq]]), "VCF/VEP, rsID того же аллеля")
        for x in xy_rows:
            p.group(dict(zip(CLINVAR, (x.get(f) for f in ("clinvar_significance", "clinvar_review", "clinvar_conditions", "clinvar_ids")))),
                    "ClinVar из локальной аннотации X/Y")
        cv = auto.get("clinvar", {})
        if call.get("n_alts") == 1:
            p.group({"ClinVar": cv.get("ClinVar_CLNSIG"), "ClinVar review": cv.get("ClinVar_CLNREVSTAT"),
                     "ClinVar заболевание": cv.get("ClinVar_CLNDN")}, "ClinVar из аннотированного VCF")
        vectors = [tuple(number(v, integer=i >= 4) for i, v in enumerate(x[2:10])) for x in auto.get("splice", [])
                   if len(x) == 10 and x[0] == k[3] and gene and x[1] == gene]
        if len(set(vectors)) == 1 and all(v is not None for v in vectors[0]) and all(0 <= v <= 1 for v in vectors[0][:4]):
            p.group(dict(zip(SPLICE, vectors[0])), "SpliceAI из аннотированного VCF, точный аллель и ген")
        return p

    @staticmethod
    def qc_proposals(p, record, sample):
        if not record:
            return
        source = "QC VCF, точный аллель и образец"
        p.put("QUAL", number(record["qual"]), source)
        p.put("Фильтр сайта", record["filter"], source)
        for name, genotype in record["samples"].items():
            p.put("GT " + name, genotype.get("GT"), source)
        if not sample or sample not in record["samples"]:
            return
        gen = record["samples"][sample]
        for field, value in {"GT": gen.get("GT"), "AD": gen.get("AD"), "DP": number(gen.get("DP"), True),
                             "GQ": number(gen.get("GQ"), True), "Фильтр генотипа": gen.get("FT")}.items():
            p.put(field, value, source)
        depths = [number(x, True) for x in str(gen.get("AD") or "").split(",")]
        index = record["allele_number"]
        if len(depths) == record["n_alts"] + 1 and all(v is not None and v >= 0 for v in depths) and sum(depths) > 0:
            p.put("Доля ALT в AD", depths[index] / sum(depths), "Рассчитано из AD того же аллеля и образца")

    def opencravat_proposals(self, p, jobs, gene, transcript):
        excluded = set(CLINVAR + GNOMAD + GTEX + MITOMAP)
        for label, rec in jobs:
            p.put("Статус OpenCRAVAT", "Есть строка результата", label)
            for field, src in ALLELE_FIELDS.items():
                if field not in excluded:
                    p.put(field, rec.get(src), label)
            for names in (CLINVAR, GNOMAD, GTEX):
                p.group({f: rec.get(ALLELE_FIELDS[f]) for f in names}, label)
            p.put("dbSNP", rsids([rec.get("clinvar__dbsnp_id")]), label + ", ClinVar dbSNP")
            exact_models = get_model_rows(rec, transcript) if transcript else {}
            for values in exact_models.values():
                p.group(values, label + ", совпавший транскрипт")
            for names, fields in (
                (MODEL_GROUPS[0], ("sift__score","sift__prediction","sift__transcript")),
                (MODEL_GROUPS[2], ("revel__score","revel__transcript")),
                (MODEL_GROUPS[3], ("alphamissense__am_pathogenicity","alphamissense__transcript_id")),
                (MODEL_GROUPS[4], ("fathmm__fathmm_score","fathmm__ens_tid")),
            ):
                if clean(rec.get(fields[0])) is not None:
                    p.group(dict(zip(names, (rec.get(f) for f in fields))), label)
            if clean(rec.get("polyphen2__hvar_pred")):
                prediction = rec["polyphen2__hvar_pred"]
                matches = {r[0] for r in json_value(rec.get("polyphen2__all"), []) if len(r) >= 8 and r[7] == prediction and r[0]}
                p.group({"PolyPhen": prediction, "PolyPhen transcript": ";".join(sorted(matches)) if matches else None}, label)
            splice = {f: rec.get("spliceai__" + f.split()[1].lower()) for f in SPLICE}
            if all(clean(v) is not None for v in splice.values()):
                p.group(splice, label)
        if gene:
            for label, rec in self.genes.get(gene, []):
                for field, src in GENE_FIELDS.items():
                    if field not in GENE_GROUP:
                        p.put(field, rec.get(src), label + ", тот же ген")
                p.group({f: rec.get(GENE_FIELDS[f]) for f in GENE_GROUP}, label + ", тот же ген")

    @staticmethod
    def fathmm_proposal(p, score):
        raw = score["scores"]
        p.group({"FATHMM": number(raw) if ";" not in raw else raw, "FATHMM transcript": score["transcripts"]},
                "FATHMM v2.3, точный аллель X/Y")

    def gene_proposals(self, p, gene):
        if not gene:
            return
        terms = ";".join(sorted(term for term in self.hpo.get(gene, ()) if clean(term)))
        if len(terms) <= 32767:
            p.put("HPO", terms, "HPO 2026-09-01, точный символ гена")
        p.put("NCBI Gene ID", one(self.entrez.get(gene, ())), "HPO 2026-09-01, единственный NCBI Gene ID")
        records = self.clingen.get(gene, ())
        pairs = sorted({(r["DISEASE LABEL"], r["CLASSIFICATION"]) for r in records if clean(r.get("DISEASE LABEL"))})
        if pairs:
            p.group({"ClinGen заболевание": ";".join(x[0] for x in pairs),
                     "ClinGen классификация": ";".join(x[1] for x in pairs)},
                    "ClinGen Gene Validity 2026-11-05, точный символ гена")
        context = sorted({f"{r['DISEASE LABEL']} ({r['CLASSIFICATION']}; {r['DISEASE ID (MONDO)']})" for r in records
                          if r["CLASSIFICATION"].lower() in {"definitive", "strong", "moderate"} and
                          (CVD_TERMS.search(r["DISEASE LABEL"]) or CVD_TERMS.search(r["GCEP"]))})
        p.put("ССЗ-контекст ClinGen", "; ".join(context), "ClinGen, исходное правило отбора ССЗ")

    def technical_proposals(self, name, row):
        p = Proposal(row)
        if name == "ГКМП_поддержка_BAM":
            records = self.bam.get(workbook_allele_key(row) + (row.get("Образец"), row.get("Ген")), [])
            flags = []
            for r in records:
                counts = json_value(r["counts"], {})
                if r["status"] == "ASSESSED_SNV_READ_SUPPORT" and all(
                    counts.get(source) == row[field] for field, source in BAM_COUNTS.items()
                ):
                    labels = json_value(r["descriptive_flags"], None)
                    if isinstance(labels, list):
                        flags.append("; ".join(labels) if labels else "Без описательных флагов")
            p.put("Признаки", one(flags), "BAM read-support, тот же аллель, образец и счётчики")
        if name == "ГКМП_покрытие":
            p.put("Хромосома", one(self.gene_chrom.get(row.get("Ген"), ())), "VEP 115, единственная хромосома гена")
            location = self.refseq_locations.get(row.get("Ген"))
            if location:
                p.put("Хромосома", location[0], "NCBI Gene, GRCh38.p14, " + location[1])
        return p

    @staticmethod
    def mt_proposals(p, variants, jobs):
        sample = p.values.get("Образец")
        transcript = clean(p.values.get("Транскрипт"))
        selected = [r for r in variants if r.get("sample") == sample and
                    (not transcript or r.get("transcript") == transcript)]
        source = "Аннотация MT rCRS, точный аллель и образец"
        for field, raw in {"HGVS-c": "hgvs_c", "HGVS-p": "hgvs_p", "GT": "GT", "DP": "DP",
                           "AD": "AD", "AF вызова": "AF", "HL исходный": "HL", "FILTER": "FILTER",
                           "gnomAD MT AF hom": "gnomAD_MT_AF_hom", "gnomAD MT AF het": "gnomAD_MT_AF_het",
                           "Helix MT AF hom": "Helix_AF_hom", "Helix MT AF het": "Helix_AF_het"}.items():
            value = one(r.get(raw) for r in selected)
            if field in {"DP", "AF вызова", "HL исходный", "gnomAD MT AF hom", "gnomAD MT AF het", "Helix MT AF hom", "Helix MT AF het"}:
                value = number(value, integer=field == "DP")
            p.put(field, value, source)
        for label, rec in jobs:
            p.group({f: rec.get(ALLELE_FIELDS[f]) for f in CLINVAR}, label)
            p.group({f: rec.get(ALLELE_FIELDS[f]) for f in MITOMAP}, label)
        for r in selected:
            p.group(dict(zip(CLINVAR[:3], (r.get(f) for f in ("ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN")))), source)


def row_type(sheet, header):
    return set(COORDS).issubset(header) and any(field in header for field in ANNOTATED_HEADERS)


def existing_annotations(workbook):
    inherited = {}
    for name in MAIN_SHEETS:
        sheet = workbook[name]
        values = iter(sheet.values)
        header = next(values)
        for data in values:
            row = dict(zip(header, data))
            k = workbook_allele_key(row) + (row.get("Ген"), row.get("Транскрипт"))
            inherited[k] = {field: row[field] for field in COPIED if field in row and clean(row[field]) is not None}
    return inherited


def reason(field, sheet):
    if sheet == "ГКМП_покрытие" and field != "Хромосома":
        return "У гена нет измерения в исходной CDS-выборке. Статус указан в строке."
    if sheet.startswith("FATHMM_") and field in {"Ошибка", "Примечание к проверке"}:
        return "Примечание указывается, если оно относится к результату строки."
    if field in SCORE_COLUMNS:
        return "Модель не вернула оценку для данного аллеля или транскрипта либо эффект вне области модели."
    if field in CLINVAR or field in MITOMAP or field in {"dbSNP", "OMIM ID"}:
        return "Нет подтверждённого значения нужного поля в сверенных записях."
    if field in GNOMAD or field.startswith(("gnomAD", "Helix")):
        return "Нет сопоставимого поля частоты или числа аллелей в сверенных ресурсах."
    if field in GENE_FIELDS or field in GTEX or field == "ССЗ-контекст ClinGen":
        return "В сверенных источниках нет значения для этого гена или ассоциации аллеля."
    if field in {"Ранг", "Ранг в списке ССЗ"}:
        return "На листах без ранжирования исходный ранг отсутствует."
    if field.startswith("GT ") or field in {"GT", "DP", "GQ", "AD", "HL исходный", "Фильтр генотипа"}:
        return "В исходной записи образца это поле отсутствует или генотип не вызван."
    if field in {"Ген", "Транскрипт", "HGVS-c", "HGVS-p", "Consequence выбранного транскрипта"}:
        return "Не подтверждено однозначное значение для этого аллеля и транскрипта."
    return "Нет подтверждённого значения в доступных источниках или поле неприменимо к строке."


def add_table(book, name, header, rows, widths):
    sheet = book.create_sheet(name)
    sheet.append(header)
    for row in rows:
        sheet.append(row)
    sheet.freeze_panes = "C2"
    sheet.sheet_view.showGridLines = False
    for c in sheet[1]:
        c.fill = PatternFill("solid", fgColor="23425D")
        c.font = Font(color="FFFFFF", bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="center")
    sheet.row_dimensions[1].height = 33
    sheet.auto_filter.ref = sheet.dimensions
    for i, width in enumerate(widths, 1):
        sheet.column_dimensions[get_column_letter(i)].width = width
    return sheet


def build(input_path, output, summary, changes, evidence):
    book = load_workbook(input_path)
    inherited = existing_annotations(book)
    edits = defaultdict(dict)
    provenance = Counter()
    counts = []
    for sheet in book:
        header = tuple(c.value for c in sheet[1])
        annotation_table = row_type(sheet, header)
        has_coords = {"REF", "ALT"}.issubset(header)
        before = Counter()
        after = Counter()
        total = 0
        for cells in sheet.iter_rows(min_row=2):
            row = dict(zip(header, (c.value for c in cells)))
            if has_coords and (clean(row["REF"]) is None or clean(row["ALT"]) is None):
                continue
            if not any(clean(v) is not None for v in row.values()):
                continue
            total += 1
            before.update(f for f in header if clean(row[f]) is not None)
            if annotation_table:
                k = workbook_allele_key(row) + (row.get("Ген"), row.get("Транскрипт"))
                proposal = evidence.propose(row, inherited.get(k, {}))
            else:
                proposal = evidence.technical_proposals(sheet.title, row)
            for col, field in enumerate(header):
                if field in proposal.changes:
                    value, source = proposal.changes[field]
                    cell = cells[col]
                    cell.value = value
                    edits[sheet.title][cell.coordinate] = value
                    provenance[(source, field)] += 1
            after.update(f for f in header if clean(proposal.values[f]) is not None)
        for field in header:
            counts.append([sheet.title, field, total, before[field], after[field], after[field] - before[field],
                           total - after[field], reason(field, sheet.title) if total > after[field] else None])
        print(f"{sheet.title.encode('ascii', 'backslashreplace').decode()}: {total} rows, {len(edits[sheet.title])} fills", flush=True)

    coverage_header = ["Лист", "Столбец", "Строк", "Заполнено до", "Заполнено после", "Добавлено", "Пустых после", "Причина оставшихся пропусков"]
    coverage = add_table(book, COVERAGE, coverage_header, counts, [29, 38, 11, 15, 17, 13, 14, 90])
    for i in range(2, coverage.max_row + 1):
        coverage.cell(i, 8).alignment = Alignment(wrap_text=True, vertical="center")
        coverage.row_dimensions[i].height = 30
    provenance_rows = [[source, field, count] for (source, field), count in sorted(provenance.items())]
    source_sheet = add_table(book, PROVENANCE, ["Источник", "Столбец", "Заполнено ячеек"], provenance_rows, [67, 39, 19])
    source_sheet.cell(source_sheet.max_row + 2, 1, "Ранг и критерии рейтинга отражают исходный отбор. После дополнения аннотаций порядок вариантов не пересчитывался.")
    sheet = book["Полнота_аннотаций"]
    sheet.append(["Все таблицы вариантов", "Локальная сверка", "ячейки", None, sum(provenance.values()), None,
                  "Заполнено подтверждёнными данными", "Подробные счётчики по каждому столбцу на листе Полнота_по_столбцам."])
    for cell in sheet[sheet.max_row]:
        cell._style = copy(sheet.cell(4, cell.column)._style)
        cell.alignment = Alignment(vertical="center", wrap_text=True)
    sheet.row_dimensions[sheet.max_row].height = 64
    if sheet.auto_filter.ref:
        sheet.auto_filter.ref = f"A1:H{sheet.max_row}"
    partial = output.with_suffix(".writing.xlsx")
    book.save(partial)
    book.close()
    partial.replace(output)
    with summary.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(coverage_header)
        writer.writerows(counts)
    changes.write_text(json.dumps({"edits": edits, "counts": counts, "provenance": provenance_rows}, ensure_ascii=False),
                       encoding="utf-8")
    print(f"Saved {sum(provenance.values())} confirmed additions.", flush=True)


def verify(input_path, output, changes):
    audit = json.loads(changes.read_text(encoding="utf-8"))
    before = load_workbook(input_path, read_only=True)
    after = load_workbook(output, read_only=True)
    assert after.sheetnames == before.sheetnames + [COVERAGE, PROVENANCE]
    checked = 0
    for name in before.sheetnames:
        edits = audit["edits"].get(name, {})
        for row_no, (old, new) in enumerate(zip_longest(before[name].values, after[name].values), 1):
            if name == "Полнота_аннотаций" and old is None:
                assert new and new[1] == "Локальная сверка"
                continue
            assert old is not None and new is not None, (name, row_no)
            for col, (a, b) in enumerate(zip_longest(old, new), 1):
                expected = a
                if edits and (cell := f"{get_column_letter(col)}{row_no}") in edits:
                    assert clean(a) is None
                    expected = edits[cell]
                    checked += 1
                assert identical(expected, b) or clean(expected) is None and clean(b) is None, (name, row_no, col)
    counts = list(after[COVERAGE].values)[1:]
    assert len(counts) == len(audit["counts"])
    assert checked == sum(len(edits) for edits in audit["edits"].values()) == sum(x[5] for x in counts)
    for a, b in zip(counts, audit["counts"], strict=True):
        assert tuple(b) == a
    before.close()
    after.close()
    print(f"Verified {checked} filled cells and all unchanged workbook values.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run_parser = commands.add_parser("run", help="extract evidence, update a copy, and verify it")
    run_parser.add_argument("--manifest", type=Path, required=True)
    run_parser.add_argument("--summary", type=Path, required=True)
    run_parser.add_argument("--work-dir", type=Path, required=True)
    verify_parser = commands.add_parser("verify", help="verify an existing workbook against its audit")
    for subparser in (run_parser, verify_parser):
        subparser.add_argument("--input", type=Path, required=True)
        subparser.add_argument("--output", type=Path, required=True)
        subparser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "verify":
        verify(args.input, args.output, args.audit)
        return

    config_path = args.manifest.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    path_keys = ("annotated_vcf", "autosome_qc_vcf", "xy_qc_vcf", "xy_annotations",
                 "mt_annotations", "hpo", "clingen", "fathmm_database", "fathmm_code")
    paths = {name: config_path.parent / config[name] for name in path_keys if config.get(name)}
    jobs = [(job["label"], config_path.parent / job["sqlite"]) for job in config.get("opencravat", ())]
    bam_sources = [config_path.parent / item for item in config.get("bam_support_tsv", ())]
    args.work_dir.mkdir(parents=True, exist_ok=True)
    sources = args.work_dir / "reconciliation_sources.json"
    extract = args.work_dir / "reconciliation_vcf.jsonl"
    outputs = (args.output, args.summary, args.audit, sources, extract)
    if len({path.resolve() for path in (*outputs, args.input)}) != len(outputs) + 1:
        parser.error("Input workbook and output paths must be distinct")
    for path in outputs:
        if path.exists():
            parser.error(f"Output path already exists: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
    collect_evidence(args.input, paths, jobs, extract, sources)
    evidence = Evidence(sources, paths, bam_sources, config.get("refseq_gene_locations", {}))
    build(args.input, args.output, args.summary, args.audit, evidence)
    verify(args.input, args.output, args.audit)


if __name__ == "__main__":
    main()
