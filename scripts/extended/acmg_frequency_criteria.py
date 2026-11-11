"""Frequency and computational ACMG/AMP criteria for the HCM gene-screen candidates.

Thresholds: BA1 at FAF >= 0.1%, BS1 at FAF >= 0.02%, PM2_Supporting at FAF < 0.004%
(ClinGen MYH7 specifications, Kelly et al., 2018; applied by analogy to other dominant
genes, not to genes with a recessive association). REVEL: BP4_Supporting <= 0.290,
BP4_Moderate <= 0.183, PP3_Supporting >= 0.644 (Pejaver et al., 2022).

FAF95 is approximated: the allele count is reconstructed from the gnomAD v4.1 exome AF
and an approximate allele number per group. Check the values on the gnomAD site before
publication. Criteria that need family, segregation or functional data are not assessed.
"""
import argparse
import json

import pandas as pd
from scipy.stats import chi2

from ext_common import load_sheet

ALLELE_NUMBER = {"AFR": 33480, "AMR": 44724, "EAS": 44896, "NFE": 1180062, "SAS": 86258}
RECESSIVE = {"TRIM63"}

parser = argparse.ArgumentParser()
parser.add_argument("--workbook", required=True)
parser.add_argument("--table", required=True, help="var.tsv from extract_vcf_table.py")
parser.add_argument("--out", required=True)
args = parser.parse_args()

cols = ["chrom", "pos", "ref", "alt"] + [f"gnomADe_{g}_AF" for g in ALLELE_NUMBER]
v = pd.read_csv(args.table, sep="\t", usecols=cols, dtype=str)
v["key"] = v.chrom + ":" + v.pos + ":" + v.ref + ":" + v.alt
v = v.set_index("key")
result = []
for r in load_sheet(args.workbook, "ГКМП_кандидаты"):
    z = v.loc[f'{r["Хромосома"]}:{r["Позиция"]}:{r["REF"]}:{r["ALT"]}']
    best = None
    for group, an in ALLELE_NUMBER.items():
        af = float(z[f"gnomADe_{group}_AF"] or 0) if pd.notna(z[f"gnomADe_{group}_AF"]) else 0.0
        ac = int(round(af * an))
        faf = chi2.ppf(0.05, 2 * ac) / 2 / an if ac else 0.0
        if best is None or faf > best[1]:
            best = (group, faf, af)
    group, faf, af = best
    criteria = []
    if r["Ген"] not in RECESSIVE:
        criteria.append("BA1" if faf >= 0.001 else "BS1" if faf >= 0.0002 else "PM2_Supporting" if faf < 0.00004 else None)
    if r["Consequence выбранного транскрипта"] == "missense_variant" and r["REVEL"] is not None:
        revel = r["REVEL"]
        criteria.append("BP4_Moderate" if revel <= 0.183 else "BP4_Supporting" if revel <= 0.290 else "PP3_Supporting" if revel >= 0.644 else None)
    criteria = [c for c in criteria if c]
    outcome = "B" if "BA1" in criteria else "LB (provisional)" if "BS1" in criteria else "VUS"
    result.append({"gene": r["Ген"], "protein": r["HGVS-p"], "group": group, "af_pct": round(af * 100, 4),
                   "faf95_pct": round(float(faf) * 100, 4), "revel": r["REVEL"], "criteria": criteria, "class": outcome})
    print(result[-1])
json.dump(result, open(args.out, "w"), ensure_ascii=False, indent=1)
