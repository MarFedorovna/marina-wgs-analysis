"""Supervised ancestry proportions from gnomAD group allele frequencies.

Each genotype is modelled as binomial with an allele frequency that is a mixture of
gnomAD genome group frequencies; proportions are fitted by maximum likelihood on the
simplex. To avoid ascertainment through the sample itself, sites are chosen where at
least one of the other samples carries the allele. One SNP per 100 kb; bootstrap 95% CI.
This is an approximation and does not replace PCA with a reference panel.
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from ext_common import AUTOSOMES, genotype, labels, sample_names

parser = argparse.ArgumentParser()
parser.add_argument("--table", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--groups", default="NFE,EAS,SAS,AFR")
parser.add_argument("--bootstrap", type=int, default=200)
args = parser.parse_args()

S = sample_names(args.table)
L = labels(S)
groups = args.groups.split(",")
cols = [f"gnomADg_{g}_AF" for g in groups]
use = ["chrom", "pos", "ref", "alt", "nalt", "filter", "gnomADg_AF"] + cols + [f"{s}_{f}" for s in S for f in ("gt", "dp", "gq")]
d = pd.read_csv(args.table, sep="\t", usecols=use, dtype={f"{s}_gt": str for s in S} | {"chrom": str, "filter": str, "ref": str, "alt": str}, low_memory=False)
keep = d.chrom.isin(AUTOSOMES) & (d.nalt == 1) & (d.ref.str.len() == 1) & (d.alt.str.len() == 1) & (d["filter"] == "PASS")
for s in S:
    d[s] = genotype(d[f"{s}_gt"])
    keep &= (d[f"{s}_dp"] >= 10) & (d[f"{s}_gq"] >= 20) & d[s].notna()
x = d[keep].copy()
for c in cols + ["gnomADg_AF"]:
    x[c] = pd.to_numeric(x[c], errors="coerce")


def fit(g, P):
    k = P.shape[1]

    def nll(q):
        p = P @ q
        return -(g * np.log(p) + (2 - g) * np.log(1 - p)).sum()

    def grad(q):
        p = P @ q
        return -((g / p - (2 - g) / (1 - p))[:, None] * P).sum(0)

    constraint = {"type": "eq", "fun": lambda q: q.sum() - 1, "jac": lambda q: np.ones(k)}
    return minimize(nll, np.ones(k) / k, jac=grad, method="SLSQP", bounds=[(0, 1)] * k,
                    constraints=[constraint], options={"maxiter": 300, "ftol": 1e-10}).x


rng = np.random.default_rng(1)
result = {}
for s in S:
    others = [o for o in S if o != s]
    sel = (x[others].sum(1) > 0) & x[cols].notna().all(1) & (x.gnomADg_AF > 0.05) & (x.gnomADg_AF < 0.95)
    y = x[sel]
    y = y[~(y.chrom + ":" + (y.pos // 100000).astype(str)).duplicated()]
    P = np.clip(y[cols].values, 1e-3, 1 - 1e-3)
    g = y[s].values.astype(float)
    q = fit(g, P)
    boot = np.array([fit(g[i], P[i]) for i in (rng.integers(0, len(g), len(g)) for _ in range(args.bootstrap))])
    low, high = np.percentile(boot, [2.5, 97.5], axis=0)
    result[L[s]] = {"n_snps": int(len(g)), "estimate": dict(zip(groups, q.round(3).tolist())),
                    "ci95": {p: [round(float(low[i]), 3), round(float(high[i]), 3)] for i, p in enumerate(groups)}}
    print(L[s], result[L[s]])
json.dump(result, open(args.out, "w"), ensure_ascii=False, indent=1)
