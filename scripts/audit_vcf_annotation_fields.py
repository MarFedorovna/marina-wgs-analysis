#!/usr/bin/env python3
"""Report declared and populated annotation tags in a VCF without changing it."""
import argparse
import csv
import gzip
from pathlib import Path

INFO_FIELDS = ("CSQ", "SpliceAI", "SpliceAI_pred", "FATHMM", "FATHMM_pred", "CADD", "CADD_PHRED",
               "REVEL", "REVEL_score", "AlphaMissense", "AlphaMissense_pred", "SIFT", "PolyPhen",
               "CLNSIG", "CLNREVSTAT", "gnomADg_AF", "gnomADe_AF", "AF", "RS")
FORMAT_FIELDS = ("GT", "DP", "GQ", "PS", "PID", "PGT", "FT")


def opener(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith((".gz", ".bgz")) else open(path, encoding="utf-8")


def present(value):
    return value not in (None, "", ".") and any(part != "." for part in value.split(","))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    declared = {"INFO": set(), "FORMAT": set()}
    counts = {kind: {tag: 0 for tag in tags} for kind, tags in (("INFO", INFO_FIELDS), ("FORMAT", FORMAT_FIELDS))}
    total = 0
    with opener(args.vcf) as stream:
        for line in stream:
            if line.startswith("##INFO=<ID="):
                declared["INFO"].add(line.split("##INFO=<ID=", 1)[1].split(",", 1)[0])
            elif line.startswith("##FORMAT=<ID="):
                declared["FORMAT"].add(line.split("##FORMAT=<ID=", 1)[1].split(",", 1)[0])
            elif not line.startswith("#"):
                cols = line.rstrip().split(chr(9))
                if len(cols) < 8:
                    continue
                total += 1
                info = {}
                for item in cols[7].split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info[key] = value
                for tag in INFO_FIELDS:
                    if present(info.get(tag)):
                        counts["INFO"][tag] += 1
                if len(cols) > 9:
                    fmt = cols[8].split(":")
                    values = [sample.split(":") for sample in cols[9:]]
                    for tag in FORMAT_FIELDS:
                        if tag in fmt:
                            ix = fmt.index(tag)
                            if any(ix < len(sample) and present(sample[ix]) for sample in values):
                                counts["FORMAT"][tag] += 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter=chr(9))
        w.writerow(("namespace", "field", "declared_in_header", "records_with_value", "records_total"))
        for kind in ("INFO", "FORMAT"):
            for tag in counts[kind]:
                w.writerow((kind, tag, str(tag in declared[kind]).lower(), counts[kind][tag], total))
    print("records=%d; report=%s" % (total, args.out))


if __name__ == "__main__":
    main()
