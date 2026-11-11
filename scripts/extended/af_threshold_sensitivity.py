"""Counts of shortlisted variants under progressively stricter maximum-AF thresholds.

Uses the final workbook. Variants without a frequency are kept at every threshold.
The 0.004% threshold is the maximum credible population AF for HCM (Whiffin et al., 2017).
"""
import argparse
import json

from ext_common import load_sheet

parser = argparse.ArgumentParser()
parser.add_argument("--workbook", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()


def locus(row):
    return tuple(row[k] for k in ("Хромосома", "Позиция", "REF", "ALT"))


overall = load_sheet(args.workbook, "Все_прошедшие")
candidates = {locus(r) for r in load_sheet(args.workbook, "ГКМП_кандидаты")}
protein = set()
for sheet in ("Missense", "Frameshift", "Inframe", "Stop_gained", "Stop_lost", "Start_lost"):
    protein |= {locus(r) for r in load_sheet(args.workbook, sheet)}
plof = {locus(r) for r in load_sheet(args.workbook, "pLoF_candidates")}
result = []
for label, limit in (("<=1%", 0.01), ("<=0.1%", 0.001), ("<=0.01%", 0.0001), ("<=0.004%", 0.00004), ("no frequency", None)):
    kept = [r for r in overall if r["gnomAD AF max"] is None or (limit is not None and r["gnomAD AF max"] <= limit)]
    keys = {locus(r) for r in kept}
    result.append({
        "threshold": label, "total": len(kept),
        "without_af": sum(r["gnomAD AF max"] is None for r in kept),
        "protein_altering": len(keys & protein), "plof": len(keys & plof),
        "clingen_cvd_genes": sum(bool(r["ССЗ-контекст ClinGen"]) for r in kept),
        "two_or_more_carriers": sum(r["Число носителей"] >= 2 for r in kept),
        "hcm_candidates": len(keys & candidates),
    })
    print(result[-1])
json.dump(result, open(args.out, "w"), ensure_ascii=False, indent=1)
