"""Exon read-depth screen for deletions and duplications in a gene list.

Mean depth is computed over the coding part of each exon of the MANE Select transcript
(cdot JSON built from the Ensembl GTF, https://github.com/SACGF/cdot). Reads that are
duplicates, secondary, supplementary, QC-fail or have MAPQ < 20 are ignored. Depth is
divided by the sample median over autosomal exons and compared with the median of the
other samples. chrX depth is rescaled by the X copy number of each sample.

With a handful of samples the single-exon noise is high, so only homozygous losses,
whole-gene shifts and runs of >= 3 adjacent exons are interpretable.
"""
import argparse
import gzip
import json

import numpy as np
import pandas as pd
import pysam

parser = argparse.ArgumentParser()
parser.add_argument("--bam", action="append", required=True, help="LABEL=path, repeat per sample")
parser.add_argument("--x-copies", required=True, help="comma-separated X copy number per sample, in --bam order")
parser.add_argument("--cdot", required=True)
parser.add_argument("--genes", required=True, help="text file, one gene symbol per line")
parser.add_argument("--out-prefix", required=True)
args = parser.parse_args()

bams = dict(item.split("=", 1) for item in args.bam)
S = list(bams)
x_copies = dict(zip(S, map(int, args.x_copies.split(","))))
genes = [g.strip() for g in open(args.genes) if g.strip()]
contigs = {f"NC_0000{i:02d}": f"chr{i}" for i in range(1, 23)} | {"NC_000023": "chrX"}
exons = []
for tid, t in json.load(gzip.open(args.cdot))["transcripts"].items():
    build = t["genome_builds"].get("GRCh38")
    if t.get("gene_name") in genes and build and "MANE_Select" in (build.get("tag") or ""):
        chrom = contigs.get(build["contig"].split(".")[0])
        if chrom and build.get("cds_start") is not None:
            for e in sorted(build["exons"]):
                start, end = max(e[0], build["cds_start"]), min(e[1], build["cds_end"])
                if end > start:
                    exons.append((t["gene_name"], e[2] + 1, chrom, start, end))
table = pd.DataFrame(exons, columns=["gene", "exon", "chrom", "start", "end"])
for s, path in bams.items():
    bam = pysam.AlignmentFile(path)
    depth = []
    for _, chrom, start, end in table[["gene", "chrom", "start", "end"]].itertuples(index=False):
        bases = 0
        for r in bam.fetch(chrom, start, end):
            if r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail or r.mapping_quality < 20:
                continue
            bases += r.get_overlap(start, end) or 0
        depth.append(bases / (end - start))
    table[s] = depth
is_x = table.chrom == "chrX"
norm = table[S] / table.loc[~is_x, S].median()
for s in S:
    norm.loc[is_x, s] *= 2 / x_copies[s]
usable = (norm.median(axis=1) >= 0.3) & (table[S].median(axis=1) >= 15)
for s in S:
    others = [o for o in S if o != s]
    table[f"log2_{s}"] = np.log2((norm[s] + 1e-3) / (norm[others].median(axis=1) + 1e-3))
table["usable"] = usable
summary = {"exons": int(len(table)), "usable": int(usable.sum()), "noise": {}, "hom_loss": {}, "runs": [], "flagged": {}}
ok = table[usable].sort_values(["gene", "exon"])
for s in S:
    values = ok.loc[ok.chrom != "chrX", f"log2_{s}"]
    sd = 1.4826 * np.median(np.abs(values - values.median()))
    z = (table[f"log2_{s}"] - values.median()) / sd
    hit = usable & ((table[f"log2_{s}"] <= -0.6) | (table[f"log2_{s}"] >= 0.45)) & (z.abs() >= 3)
    summary["noise"][s] = round(float(sd), 3)
    summary["flagged"][s] = int(hit.sum())
    others = table[[o for o in S if o != s]].median(axis=1)
    lost = table[(table[s] < 0.05 * others) & (others >= 15)]
    summary["hom_loss"][s] = [f"{g}:{e}" for g, e in zip(lost.gene, lost.exon)]
    for gene, grp in ok.groupby("gene"):
        v, ex = grp[f"log2_{s}"].values, grp.exon.values
        for i in range(len(v) - 2):
            w = v[i:i + 3]
            if ex[i + 2] - ex[i] == 2 and ((w.mean() <= -0.6 and (w <= -0.35).all()) or (w.mean() >= 0.45 and (w >= 0.3).all())):
                summary["runs"].append({"sample": s, "gene": gene, "exons": f"{ex[i]}-{ex[i + 2]}", "mean_log2": round(float(w.mean()), 2)})
gene_mean = ok.groupby("gene")[[f"log2_{s}" for s in S]].mean().round(2)
gene_mean["usable_exons"] = ok.groupby("gene").size()
print(json.dumps(summary, indent=1))
table.to_csv(args.out_prefix + ".exons.tsv", sep="\t", index=False)
gene_mean.to_csv(args.out_prefix + ".genes.tsv", sep="\t")
json.dump(summary, open(args.out_prefix + ".summary.json", "w"), indent=1)
