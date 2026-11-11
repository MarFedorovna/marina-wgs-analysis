"""All alleles in the nuclear HCM genes, by position relative to exons and by rarity.

Inputs are the two outputs of extract_vcf_table.py. A carrier needs FILTER = PASS at
the site, an alternate genotype, FT = PASS, DP >= 10 and GQ >= 20. The position class
is taken from the MANE Select (or canonical) transcript of the listed gene.
"""
import argparse
import json
import re

import numpy as np
import pandas as pd

from ext_common import HCM_GENES, labels, sample_names

parser = argparse.ArgumentParser()
parser.add_argument("--table", required=True)
parser.add_argument("--genes", required=True)
parser.add_argument("--out-prefix", required=True)
args = parser.parse_args()

S = sample_names(args.table)
L = labels(S)
v = pd.read_csv(args.table, sep="\t", dtype=str, low_memory=False)
v = v[v.genes.fillna("").apply(lambda g: bool(set(g.split(",")) & set(HCM_GENES)))].copy()
v["key"] = v.chrom + ":" + v.pos + ":" + v.ref + ":" + v.alt
pop_cols = [c for c in v.columns if re.match(r"gnomAD[eg]_(AFR|AMR|EAS|NFE|SAS|MID|FIN|ASJ)_AF", c)] + ["gnomADe_AF", "gnomADg_AF"]
v["af_max"] = v[pop_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)


def carriers(row):
    found = []
    for s in S:
        gt = row[f"{s}_gt"].replace("|", "/")
        if "1" in gt and row[f"{s}_ft"] == "PASS" and float(row[f"{s}_dp"] or 0) >= 10 and float(row[f"{s}_gq"] or 0) >= 20:
            found.append(L[s] + ("(hom)" if gt == "1/1" else ""))
    return "+".join(found)


v["carriers"] = v.apply(carriers, axis=1)
blocks = {}
for line in open(args.genes):
    j = json.loads(line)
    blocks[f'{j["chrom"]}:{j["pos"]}:{j["ref"]}:{j["alt"]}'] = j["csq"]


def pick(key):
    c = [b for b in blocks.get(key, []) if b["SYMBOL"] in HCM_GENES and b["Feature_type"] == "Transcript"]
    best = [b for b in c if b["MANE_SELECT"]] or [b for b in c if b["CANONICAL"] == "YES"] or c
    return best[0] if best else {}


def position_class(block):
    q = block.get("Consequence", "")
    if any(t in q for t in ("missense", "frameshift", "stop_gained", "stop_lost", "start_lost", "inframe", "protein_altering")):
        return "coding, protein-altering"
    if "splice_donor_variant" in q or "splice_acceptor_variant" in q:
        return "canonical splice site"
    if "synonymous" in q or "stop_retained" in q:
        return "synonymous"
    if "splice" in q:
        return "intron, splice region"
    if "intron" in q:
        offsets = [int(z) for z in re.findall(r"[+-](\d+)", block.get("HGVSc", "").split("c.")[-1])]
        distance = min(offsets) if offsets else None
        if distance is not None and distance <= 50:
            return "intron, <=50 bp from exon"
        return "intron, 51-300 bp" if distance is not None and distance <= 300 else "intron, >300 bp"
    if "UTR" in q:
        return "UTR"
    return "outside gene (+-5 kb)" if "stream" in q else (q or "unannotated")


rows = []
for _, r in v.iterrows():
    b = pick(r.key)
    rows.append({"key": r.key, "gene": b.get("SYMBOL") or r.genes, "filter": r["filter"], "carriers": r.carriers,
                 "af_max": r.af_max, "class": position_class(b), "hgvsc": b.get("HGVSc", "").split(":")[-1],
                 "clinvar": b.get("ClinVar_CLNSIG", "") or b.get("CLIN_SIG", "")})
t = pd.DataFrame(rows)
q = t[(t["filter"] == "PASS") & (t.carriers != "")].copy()
q["rarity"] = np.where(q.af_max.isna(), "absent", np.where(q.af_max <= 0.001, "<=0.1%", np.where(q.af_max <= 0.01, "<=1%", ">1%")))
print("alleles:", len(t), "| site PASS:", int((t["filter"] == "PASS").sum()), "| with qualifying carrier:", len(q))
print(pd.crosstab(q["class"], q["rarity"], margins=True).to_string())
rare = q[(q.rarity != ">1%") & ~q["class"].isin(["coding, protein-altering", "synonymous"])]
print("rare noncoding:", len(rare))
q.to_csv(args.out_prefix + ".carriers.tsv", sep="\t", index=False)
rare.to_csv(args.out_prefix + ".rare_noncoding.tsv", sep="\t", index=False)
