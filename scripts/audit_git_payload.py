#!/usr/bin/env python3
"""Refuse a Git payload containing genomic/patient files or likely secrets."""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

BLOCKED_SUFFIXES = (
    ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".bai", ".cram", ".crai",
    ".sam", ".vcf", ".vcf.gz", ".g.vcf.gz", ".bcf", ".tbi", ".csi", ".ped",
    ".sqlite", ".db", ".sql", ".dump", ".pem", ".key", ".ppk",
)
BLOCKED_PARTS = {"data", "raw", "raw_data", "input", "inputs", "legacy", "results", "output", "outputs", "reports", "manifests", "backups", ".ssh"}
SECRET_RE = re.compile(r"(?i)(BEGIN (RSA|OPENSSH|EC) PRIVATE KEY|api[_-]?key\s*[:=]|secret\s*[:=]|token\s*[:=])")
SAMPLE_RE = re.compile(r"(?i)\bSD0(?:55|57|60|67)(?:_[A-Z0-9]+)?\b")


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True, encoding="utf-8", errors="replace")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--staged", action="store_true")
    args = parser.parse_args()
    names = git("diff", "--cached", "--name-only", "--diff-filter=ACMR") if args.staged else git("ls-files")
    failures = []
    for raw in names.splitlines():
        rel = Path(raw)
        lower = rel.as_posix().lower()
        if any(part.lower() in BLOCKED_PARTS for part in rel.parts) or any(lower.endswith(s) for s in BLOCKED_SUFFIXES):
            failures.append(f"blocked path: {raw}")
            continue
        try:
            data = Path(raw).read_bytes()
        except OSError:
            continue
        if len(data) > 10 * 1024 * 1024:
            failures.append(f"unexpected tracked file >10 MiB: {raw}")
            continue
        text = data.decode("utf-8", errors="ignore")
        if SECRET_RE.search(text):
            failures.append(f"possible secret in: {raw}")
        if SAMPLE_RE.search(text):
            failures.append(f"sample identifier in public payload: {raw}")
    if failures:
        print("\n".join(failures))
        return 2
    print(f"PASS: audited {len(names.splitlines())} paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

