"""Flatten a VEP-annotated joint VCF into a per-allele table.

Output: var.tsv (genotypes, depths, gnomAD group frequencies, gene symbols) and
genes.jsonl (all VEP consequence blocks for alleles in the genes of interest).
Run per group of contigs when the run time per call is limited, then concatenate.
"""
import argparse
import itertools
import json

import pysam

from ext_common import HCM_GENES

parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True)
parser.add_argument("--contigs", required=True, help="comma-separated, e.g. chr1,chr2")
parser.add_argument("--out-table", required=True)
parser.add_argument("--out-genes", required=True)
parser.add_argument("--header", action="store_true", help="write the header line")
args = parser.parse_args()

vcf = pysam.VariantFile(args.vcf)
fields = vcf.header.info["CSQ"].description.split("Format: ")[1].split("|")
ix = {k: i for i, k in enumerate(fields)}
pops = [k for k in fields if k.startswith("gnomADg_")] + [k for k in fields if k.startswith("gnomADe_")]
samples = list(vcf.header.samples)
genes = set(HCM_GENES)
with open(args.out_table, "w") as out, open(args.out_genes, "w") as gene_out:
    if args.header:
        cols = ["chrom", "pos", "ref", "alt", "nalt", "filter", "qual"]
        cols += [f"{s}_{f}" for s in samples for f in ("gt", "dp", "gq", "ft", "ad")]
        out.write("\t".join(cols + pops + ["genes", "worst"]) + "\n")
    records = itertools.chain.from_iterable(vcf.fetch(c) for c in args.contigs.split(","))
    for rec in records:
        csq = [c.split("|") for c in rec.info.get("CSQ", ())]
        for number, alt in enumerate(rec.alts or (), 1):
            mine = [c for c in csq if c[ix["ALLELE_NUM"]] == str(number)] or (csq if len(rec.alts) == 1 else [])

            def first(key):
                return next((c[ix[key]] for c in mine if c[ix[key]]), "")

            row = [rec.chrom, rec.pos, rec.ref, alt, len(rec.alts), ";".join(rec.filter.keys()) or ".", rec.qual]
            for s in samples:
                call = rec.samples[s]
                gt = "/".join("." if a is None else "1" if a == number else "0" if a == 0 else "2" for a in call["GT"])
                ad = call.get("AD")
                row += [gt, call.get("DP"), call.get("GQ"), call.get("FT"),
                        ",".join(map(str, ad)) if ad and None not in ad else ""]
            symbols = sorted({c[ix["SYMBOL"]] for c in mine if c[ix["SYMBOL"]]})
            row += [first(k) for k in pops] + [",".join(symbols), first("Consequence")]
            out.write("\t".join("" if v is None else str(v) for v in row) + "\n")
            if genes & set(symbols):
                blocks = [dict(zip(fields, c)) for c in mine if c[ix["SYMBOL"]] in genes]
                gene_out.write(json.dumps({"chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": alt, "csq": blocks}) + "\n")
