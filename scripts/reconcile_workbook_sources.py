import csv
import gzip
import json
import math
import sqlite3
from collections import defaultdict
from contextlib import closing
from urllib.parse import unquote

from reconcile_workbook_vcf import collect_targets, key, write_extract
from fathmm_lookup import load_annotator, score_variant


ALLELE_FIELDS = {
    "CADD PHRED": "cadd__phred", "dbSNP": "dbsnp__rsid", "OMIM ID": "omim__omim_id",
    "gnomAD4 AF": "gnomad4__af", "gnomAD4 AC": "gnomad4__ac", "gnomAD4 AN": "gnomad4__an",
    "GTEx gene": "gtex__gtex_gene", "GTEx tissue": "gtex__gtex_tissue",
    "ClinVar": "clinvar__sig", "ClinVar review": "clinvar__rev_stat",
    "ClinVar заболевание": "clinvar__disease_names", "ClinVar ID": "clinvar__id",
    "MITOMAP заболевание": "mitomap__disease", "MITOMAP статус": "mitomap__status",
    "MITOMAP score": "mitomap__score", "MITOMAP PMID": "mitomap__pubmed",
}
GENE_FIELDS = {
    "ClinGen заболевание": "clingen__disease", "ClinGen классификация": "clingen__classification",
    "HPO": "hpo__term", "UniProt": "uniprot__acc", "NCBI Gene ID": "ncbigene__entrez",
    "PubMed count": "pubmed__n", "LoFtool": "loftool__loftool_score", "pHI": "phi__phi", "Prec": "prec__prec",
}
SPLICE = tuple(f"SpliceAI {kind}_{effect}" for kind in ("DS", "DP") for effect in ("AG", "AL", "DG", "DL"))
MODEL_GROUPS = (
    ("SIFT", "SIFT prediction", "SIFT transcript"), ("PolyPhen", "PolyPhen transcript"),
    ("REVEL", "REVEL transcript"), ("AlphaMissense", "AlphaMissense transcript"),
    ("FATHMM", "FATHMM transcript"),
)
VARIANT_SELECT = {
    "base__uid", "base__chrom", "base__pos", "base__ref_base", "base__alt_base", "base__hugo",
    "base__transcript", "base__all_mappings", "original_input__chrom", "original_input__pos",
    "original_input__ref_base", "original_input__alt_base", "clinvar__dbsnp_id",
    "sift__all", "sift__score", "sift__transcript", "sift__prediction",
    "polyphen2__all", "polyphen2__hvar_pred", "revel__all", "revel__score", "revel__transcript",
    "alphamissense__am_pathogenicity", "alphamissense__transcript_id",
    "fathmm__fathmm_score", "fathmm__ens_tid",
} | set(ALLELE_FIELDS.values()) | {"spliceai__" + name.split()[1].lower() for name in SPLICE} | {
    f"gnomad4__af_{p}" for p in ("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "rem")
}


def clean(value):
    return None if value is None or isinstance(value, str) and value.strip() in {"", ".", "NA", "N/A", "./.", ".|.", "[]", "{}"} else value


def number(value, integer=False):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(result) or integer and result != int(result):
        return None
    return int(result) if integer else result


def json_value(value, default):
    try:
        return json.loads(value) if value else default
    except (TypeError, ValueError):
        return default


def stable(value):
    return str(value or "").split(".", 1)[0]


def trimmed(k):
    chrom, pos, ref, alt = k
    ref, alt = ref.replace("-", ""), alt.replace("-", "")
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return chrom, pos, ref or "-", alt or "-"


def entries(connection, table, wanted):
    fields = [row[1] for row in connection.execute(f"PRAGMA table_info({table})") if row[1] in wanted]
    return connection.execute(f"SELECT {','.join(fields)} FROM {table}")


def load_opencravat(targets, jobs):
    normalized = defaultdict(set)
    for k in targets:
        normalized[trimmed(k)].add(k)
    variants = defaultdict(list)
    genes = defaultdict(list)
    for label, path in jobs:
        if not path.is_file():
            raise FileNotFoundError(path)
        exact = defaultdict(list)
        with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)) as connection:
            connection.row_factory = sqlite3.Row
            for sqlite_row in entries(connection, "variant", VARIANT_SELECT):
                row = dict(sqlite_row)
                original = tuple(row.get(f"original_input__{field}") for field in ("chrom", "pos", "ref_base", "alt_base"))
                if all(v is not None for v in original):
                    matches = {key(*original)} & targets
                else:
                    matches = set()
                if not matches:
                    normalized_key = trimmed(key(*(row[f"base__{field}"] for field in ("chrom", "pos", "ref_base", "alt_base"))))
                    matches = normalized.get(normalized_key, set())
                    if len(matches) != 1 or (all(v is not None for v in original) and trimmed(key(*original)) != normalized_key):
                        continue
                if matches:
                    exact[next(iter(matches))].append(dict(row))
            for k, rows in exact.items():
                unique = {json.dumps({f:v for f,v in row.items() if f != "base__uid"}, sort_keys=True) for row in rows}
                if len(unique) == 1:
                    variants[k].append([label, rows[0]])
            for row in entries(connection, "gene", {"base__hugo", *GENE_FIELDS.values()}):
                symbol = clean(row["base__hugo"])
                if symbol:
                    genes[symbol].append([label, dict(row)])
        print(f"{label}: {sum(len(v)==1 for v in exact.values())} single matching alleles", flush=True)
    return variants, genes


def parse_vcf_calls(path, targets):
    use = gzip.open if str(path).endswith(".gz") else open
    records = {}
    with use(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
            if line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            for alt_num, alt in enumerate(f[4].split(","), 1):
                k = key(f[0], f[1], f[3], alt)
                if k in targets:
                    values = {s: dict(zip(f[8].split(":"), v.split(":"))) for s, v in zip(samples, f[9:])}
                    rec = {"id": f[2], "qual": f[5], "filter": f[6], "samples": values, "allele_number": alt_num,
                           "n_alts": len(f[4].split(","))}
                    if k in records and records[k] != rec:
                        raise ValueError(f"Conflicting source VCF records: {path.name}, {k}")
                    records[k] = rec
    return records


def annotations_tsv(path, targets):
    grouped = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            alt = r["alt"]
            if "," in alt and r.get("matched_alt"):
                alternates = alt.split(",")
                index = int(r["allele_number"]) - 1
                if not 0 <= index < len(alternates) or alternates[index] != r["matched_alt"]:
                    continue
                alt = alternates[index]
                for field in ("AF", "HL"):
                    values = r[field].split(",")
                    r[field] = values[index] if len(values) == len(alternates) else ""
                depths = r["AD"].split(",")
                r["AD"] = ",".join((depths[0], depths[index + 1])) if len(depths) == len(alternates) + 1 else ""
            k = key(r["chrom"], r["pos"], r["ref"], alt)
            if k in targets:
                grouped[k].append(r)
    return grouped


def extract_autosomal_evidence(extract_path):
    result = {}
    with extract_path.open(encoding="utf-8") as handle:
        schema = json.loads(next(handle))["header"]["csq"]
        keep = ("SYMBOL", "Feature", "Feature_type", "Consequence", "HGVSc", "HGVSp", "SIFT", "PolyPhen",
                "Existing_variation", "MANE_SELECT", "CANONICAL", "ALLELE_NUM")
        for line in handle:
            k, row = json.loads(line)
            info = dict(item.split("=", 1) for item in row["info"].split(";") if "=" in item)
            csq = []
            for record in info.get("CSQ", "").split(","):
                values = dict(zip(schema, record.split("|")))
                if values.get("ALLELE_NUM") != str(row["allele_number"]):
                    continue
                csq.append({field: values[field] for field in (*keep, *[f for f in schema if f.startswith(("gnomADg_", "gnomADe_"))])
                            if clean(values.get(field)) is not None})
            transcript = []
            for entry in info.get("WGS_TX_SCORES", "").split(","):
                pieces = unquote(entry).split("|")
                if len(pieces) < 4 or not pieces[1]:
                    continue
                transcript.append({"sample": pieces[0], "transcript": pieces[1],
                    **dict(x.split(":", 1) for x in pieces[3:] if ":" in x)})
            frequencies = [number(info.get(f"WGS_GNOMAD_{cohort}_AF")) for cohort in ("EXOMES", "GENOMES")]
            splice = [x.split("|") for x in unquote(info.get("SpliceAI", "")).split(",")]
            clinvar = {name: unquote(info.get(name, "")) for name in
                       ("ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN", "ClinVar_ALLELEID")}
            result[tuple(k)] = {"csq": csq, "transcript": transcript, "frequencies": frequencies,
                               "splice": splice, "clinvar": clinvar, "id": row["id"]}
    return result


def local_xy_fathmm(xy_annotations, database, code_archive):
    result = {}
    with closing(sqlite3.connect(database.resolve().as_uri() + "?mode=ro&immutable=1", uri=True)) as conn:
        cursor = conn.cursor()
        available = {row[0] for row in cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        annotator = load_annotator(cursor, code_archive)
        for k, variants in xy_annotations.items():
            chrom, pos, ref, alt = k
            if "chr"+chrom not in available or len(ref) != 1 or len(alt) != 1 or alt not in "ACGT":
                continue
            missense = [r for r in variants if "missense_variant" in r.get("consequence", "").split("&")]
            for gene, transcript in sorted({(r["gene"], r["transcript"]) for r in missense if clean(r["gene"])}):
                row = (1, None, gene, "chr"+chrom, pos, ref, alt, None, None, None, None, transcript)
                scored = score_variant(row, annotator, cursor)
                if scored["scores"] is not None:
                    result[k+(gene, stable(transcript))] = scored
    return result


def collect_evidence(workbook, paths, jobs, extract_path, output_path):
    auto, xy, mt = collect_targets(workbook)
    write_extract(paths["annotated_vcf"], extract_path, (auto, xy, mt))
    oc, genes = load_opencravat(auto | xy | mt, jobs)
    auto_qc = parse_vcf_calls(paths["autosome_qc_vcf"], auto)
    xy_qc = parse_vcf_calls(paths["xy_qc_vcf"], xy)
    if set(auto_qc) != auto or set(xy_qc) != xy:
        raise ValueError("Incomplete variant matching against QC VCF")
    xa = annotations_tsv(paths["xy_annotations"], xy)
    ma = annotations_tsv(paths["mt_annotations"], mt)
    af = extract_autosomal_evidence(extract_path)
    database, code = paths.get("fathmm_database"), paths.get("fathmm_code")
    if (database is None) != (code is None):
        raise ValueError("Supply both FATHMM database and module code, or neither")
    yf = local_xy_fathmm(xa, database, code) if database else {}
    tables = {"oc": oc, "genes": genes, "qc": {**auto_qc, **xy_qc},
              "xy": xa, "mt": ma, "auto": af, "xy_fathmm": yf}
    with output_path.open("w", encoding="utf-8") as out:
        json.dump({name: list(table.items()) for name, table in tables.items()},
                  out, ensure_ascii=False, separators=(",", ":"))
    print(json.dumps({name: len(table) for name, table in tables.items()}), flush=True)
