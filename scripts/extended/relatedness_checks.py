"""Pairwise relatedness, cross-sample contamination check and runs of homozygosity.

Input is the table from extract_vcf_table.py. Sites: autosomal biallelic PASS SNPs
with gnomAD genome AF 5-95%. Relatedness uses sites where every sample has depth
>= 25, GQ >= 30 and a clean allele balance. Long stretches without opposite
homozygotes (IBS0) are compared with a block-permutation expectation, because
exome markers are sparse and such stretches also arise by chance.
"""
import argparse
import itertools
import json

import numpy as np
import pandas as pd

from ext_common import AUTOSOMES, allele_depths, genotype, labels, sample_names

parser = argparse.ArgumentParser()
parser.add_argument("--table", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--permutations", type=int, default=200)
args = parser.parse_args()

S = sample_names(args.table)
L = labels(S)
use = ["chrom", "pos", "ref", "alt", "nalt", "filter", "gnomADg_AF"] + [f"{s}_{f}" for s in S for f in ("gt", "dp", "gq", "ad")]
text = {f"{s}_{f}": str for s in S for f in ("gt", "ad")} | {"chrom": str, "filter": str, "ref": str, "alt": str}
d = pd.read_csv(args.table, sep="\t", usecols=use, dtype=text, low_memory=False)
af = pd.to_numeric(d.gnomADg_AF, errors="coerce")
snp = d.chrom.isin(AUTOSOMES) & (d.nalt == 1) & (d.ref.str.len() == 1) & (d.alt.str.len() == 1) & (d["filter"] == "PASS")
common = snp & (af > 0.05) & (af < 0.95)
G, ALT, N = {}, {}, {}
for s in S:
    G[s] = genotype(d[f"{s}_gt"])
    ALT[s], N[s] = allele_depths(d[f"{s}_ad"])
result = {"contamination": {}, "pairs": {}, "roh": {}}

# 1. Contamination: ALT reads in homozygous-reference calls, split by the partner genotype.
loose = common.copy()
for s in S:
    loose &= (d[f"{s}_dp"] >= 20) & (d[f"{s}_gq"] >= 30) & G[s].notna()
for s in S:
    ref_calls = loose & (G[s] == 0)
    row = {"alt_fraction_in_hom_ref": round(float(ALT[s][ref_calls].sum() / N[s][ref_calls].sum()), 5)}
    for o in S:
        if o != s:
            sub = ref_calls & (G[o] == 2)
            row[f"partner_{L[o]}_hom_alt"] = round(float(ALT[s][sub].sum() / N[s][sub].sum()), 5)
    het = loose & (G[s] == 1)
    row["het_alt_fraction"] = round(float(ALT[s][het].sum() / N[s][het].sum()), 4)
    result["contamination"][L[s]] = row

# 2. Relatedness on clean sites.
clean = common.copy()
for s in S:
    ab = ALT[s] / N[s]
    ok = ((G[s] == 1) & (ab >= 0.35) & (ab <= 0.65)) | ((G[s] == 0) & (ab <= 0.02)) | ((G[s] == 2) & (ab >= 0.98))
    clean &= (N[s] >= 25) & (d[f"{s}_gq"] >= 30) & ok
x = d.loc[clean, ["chrom", "pos"]].copy()
for s in S:
    x[s] = G[s][clean].astype(int)
result["clean_sites"] = int(len(x))
positions = {c: g.pos.values for c, g in x.groupby("chrom", sort=False)}


def long_runs(flags):
    """Stretches >= 5 Mb with >= 100 consecutive sites and no IBS0 site."""
    total, count, longest = 0, 0, 0
    for c, pos in positions.items():
        idx = np.flatnonzero(flags[c])
        bounds = np.concatenate(([-1], idx, [len(pos)]))
        for u, v in zip(bounds[:-1], bounds[1:]):
            if v - u - 1 >= 100 and pos[v - 1] - pos[u + 1] >= 5_000_000:
                span = int(pos[v - 1] - pos[u + 1])
                total, count, longest = total + span, count + 1, max(longest, span)
    return total, count, longest


rng = np.random.default_rng(7)
for a, b in itertools.combinations(S, 2):
    ga, gb = x[a].values, x[b].values
    het_a, het_b = (ga == 1).sum(), (gb == 1).sum()
    both_het = ((ga == 1) & (gb == 1)).sum()
    ibs0 = ((ga == 0) & (gb == 2)) | ((ga == 2) & (gb == 0))
    low = min(het_a, het_b)
    king = (both_het - 2 * ibs0.sum()) / (2 * low) + 0.5 - (het_a + het_b) / (4 * low)
    x["_ibs0"] = ibs0
    flags = {c: g["_ibs0"].values for c, g in x.groupby("chrom", sort=False)}
    total, count, longest = long_runs(flags)
    null = []
    for _ in range(args.permutations):
        shuffled = {}
        for c, z in flags.items():
            blocks = [z[i:i + 20] for i in range(0, len(z), 20)]
            shuffled[c] = np.concatenate([blocks[i] for i in rng.permutation(len(blocks))])
        null.append(long_runs(shuffled)[0])
    result["pairs"][f"{L[a]}-{L[b]}"] = {
        "ibs0_rate": round(float(ibs0.mean()), 4), "king_robust": round(float(king), 4),
        "runs": count, "runs_total_mb": round(total / 1e6, 1), "longest_mb": round(longest / 1e6, 1),
        "null_total_mb_mean": round(float(np.mean(null)) / 1e6, 1),
        "null_total_mb_p95": round(float(np.percentile(null, 95)) / 1e6, 1),
    }
    print(L[a], L[b], result["pairs"][f"{L[a]}-{L[b]}"])

# 3. Runs of homozygosity: >= 1.5 Mb, >= 50 markers, gaps < 1 Mb, at most one heterozygote per 100 markers.
called = snp.copy()
for s in S:
    called &= (d[f"{s}_dp"] >= 10) & (d[f"{s}_gq"] >= 20) & G[s].notna()
y = d.loc[called, ["chrom", "pos"]].copy()
for s in S:
    y[s] = G[s][called].astype(int)
for s in S:
    segments = []
    for c, grp in y.groupby("chrom", sort=False):
        pos, het = grp.pos.values, grp[s].values == 1
        i, n = 0, len(pos)
        while i < n:
            if het[i]:
                i += 1
                continue
            j, n_het, last = i, 0, i
            while j + 1 < n and pos[j + 1] - pos[j] < 1_000_000:
                if het[j + 1]:
                    if n_het + 1 > max(1, (j + 1 - i) // 100):
                        break
                    n_het += 1
                else:
                    last = j + 1
                j += 1
            if pos[last] - pos[i] >= 1_500_000 and last - i + 1 >= 50:
                segments.append(int(pos[last] - pos[i]))
            i = last + 1
    result["roh"][L[s]] = {"segments": len(segments), "total_mb": round(sum(segments) / 1e6, 1),
                           "longest_mb": round(max(segments, default=0) / 1e6, 1),
                           "segments_over_5mb": sum(v >= 5e6 for v in segments)}
print(result["contamination"], result["roh"], result["clean_sites"])
json.dump(result, open(args.out, "w"), ensure_ascii=False, indent=1)
