#!/usr/bin/env python3
"""Read-backed phasing for a four-sample VCF; plan-only unless --execute."""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def run(args, capture=False):
    p = subprocess.run(args, check=True, text=True, stdout=subprocess.PIPE if capture else None)
    return p.stdout.strip() if capture else ""


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require(ok, msg):
    if not ok:
        raise SystemExit(msg)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--bam-manifest", type=Path, required=True, help="Private TSV: sample_id<TAB>bam_path")
    ap.add_argument("--expected-sha256", required=True, help="From the verified transfer receipt")
    ap.add_argument("--output", type=Path, required=True, help="New output directory; must not exist")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--execute", action="store_true", help="Run only after all preflight checks pass")
    a = ap.parse_args()
    vcf, ref, manifest, out = [p.resolve() for p in (a.vcf, a.reference, a.bam_manifest, a.output)]
    for tool in ("bcftools", "samtools", "tabix", "whatshap"):
        require(shutil.which(tool), "Missing required executable: " + tool)
    require(vcf.is_file() and vcf.stat().st_size > 0, "Missing or empty VCF")
    actual = sha256(vcf)
    require(actual == a.expected_sha256.lower(), "VCF SHA-256 mismatch; refusing to run")
    require(ref.is_file() and Path(str(ref) + ".fai").is_file(), "Reference FASTA or .fai missing")
    require(Path(str(vcf) + ".tbi").is_file() or Path(str(vcf) + ".csi").is_file(), "VCF index missing")
    require(not out.exists(), "Output directory already exists; refusing to overwrite")
    samples = run(["bcftools", "query", "-l", str(vcf)], True).splitlines()
    require(len(samples) == 4 and len(set(samples)) == 4, "VCF must contain exactly four distinct samples")
    bams = {}
    for n, line in enumerate(manifest.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#") or line.lower().startswith("sample_id	"):
            continue
        fields = line.split("	")
        require(len(fields) == 2, "Manifest line %d must have sample_id<TAB>bam_path" % n)
        sample, path = fields[0], Path(fields[1]).expanduser().resolve()
        require(sample not in bams and path.is_file(), "Duplicate sample or missing BAM in manifest line %d" % n)
        require(path.suffix.lower() == ".bam", "Expected BAM for " + sample)
        require(Path(str(path) + ".bai").is_file() or path.with_suffix(".bai").is_file() or Path(str(path) + ".csi").is_file(), "BAM index missing for " + sample)
        run(["samtools", "quickcheck", "-v", str(path)])
        header = run(["samtools", "view", "-H", str(path)], True)
        sms, rids = set(), set()
        for row in header.splitlines():
            if row.startswith("@RG	"):
                tags = dict(x.split(":", 1) for x in row.split("	")[1:] if ":" in x)
                if tags.get("SM"): sms.add(tags["SM"])
                if tags.get("ID"): rids.add(tags["ID"])
        require(sms == {sample} and rids, "BAM @RG SM/ID does not match sample " + sample)
        bams[sample] = path
    require(set(bams) == set(samples), "BAM manifest sample IDs must exactly match VCF sample IDs")
    plan = {"status": "PREFLIGHT_PASS_PLAN_ONLY", "vcf_sha256": actual, "samples": samples,
            "workers": a.workers, "output": str(out)}
    require(1 <= a.workers <= 4, "workers must be 1..4")
    if not a.execute:
        print(json.dumps(plan, indent=2))
        return
    out.mkdir(parents=True)
    def phase(sample):
        d = out / sample
        d.mkdir()
        one = d / "input.vcf.gz"
        phased = d / (sample + ".phased.vcf.gz")
        run(["bcftools", "view", "-s", sample, "-Oz", "-o", str(one), str(vcf)])
        run(["bcftools", "index", "-t", str(one)])
        run(["whatshap", "phase", "--reference", str(ref), "--indels", "--tag", "PS",
             "--sample", sample, "-o", str(phased), str(one), str(bams[sample])])
        run(["bcftools", "index", "-t", "-f", str(phased)])
        run(["whatshap", "stats", "--tsv=" + str(d / "stats.tsv"),
             "--block-list=" + str(d / "blocks.tsv"), str(phased)])
        return {"sample": sample, "phased_vcf": str(phased), "stats": str(d / "stats.tsv"), "blocks": str(d / "blocks.tsv")}
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        results = list(pool.map(phase, samples))
    report = {"status": "READ_BACKED_PHASED_VCF_PER_SAMPLE", "source_vcf_sha256": actual,
              "samples": results, "workers": a.workers,
              "note": "Unphased variants remain unphased when reads do not resolve them; no parental origin or ACMG class is inferred."}
    (out / "phasing-run.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
