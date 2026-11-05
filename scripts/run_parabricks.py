#!/usr/bin/env python3
"""Single-GPU FASTQ -> marked BAM -> interval/ploidy-specific gVCF queue.

Default: validate real inputs and print commands. Only --execute runs workloads.
Private manifests, provenance, logs, and biological outputs must stay outside git.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import uuid

TAG = "nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1"
DIGEST = re.compile(r"nvcr\.io/nvidia/clara/clara-parabricks@sha256:[0-9a-f]{64}\Z")
IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,99}\Z")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def contained_file(root, relative):
    require(isinstance(relative, str) and relative, "Missing input relative path.")
    require(not Path(relative).is_absolute(), "Manifest inputs must be relative to their configured root.")
    path = (root / relative).resolve(strict=True)
    require(path.is_relative_to(root) and path.is_file(), f"Input escapes its root or is not a file: {relative}")
    return path


def checked_file(root, relative, expected):
    require(re.fullmatch(r"[0-9a-f]{64}", expected or ""), f"Provide SHA256 for {relative}.")
    path = contained_file(root, relative)
    require(sha256(path) == expected, f"SHA256 mismatch: {path}")
    return path


def read_tsv(path, columns):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames and set(columns).issubset(reader.fieldnames), f"Missing columns in {path}: {columns}")
        rows = list(reader)
    require(rows, f"Empty manifest: {path}")
    require(all(None not in row and all(v is not None and '\x00' not in v for v in row.values()) for row in rows),
            f"Malformed TSV in {path}")
    return rows


def receipt_valid(directory, expected):
    receipt = directory / "receipt.json"
    if not receipt.is_file():
        return False
    saved = json.loads(receipt.read_text(encoding="utf-8"))
    return saved.get("fingerprint") == expected and bool(saved.get("outputs")) and all(
        sha256(contained_file(directory.resolve(), name)) == digest for name, digest in saved["outputs"].items())


class Queue:
    def __init__(self, config_path, execute):
        self.cfg = json.loads(Path(config_path).read_text(encoding="utf-8-sig"))
        self.execute = execute
        self.commands = []
        self.docker_prefix = ["docker"]
        cfg = self.cfg
        require(cfg["image"] == TAG and DIGEST.fullmatch(cfg["image_digest"]),
                "Set image to pinned 4.8.0-1 and image_digest to its pulled RepoDigests[0].")
        for name in ("reads_root", "reference_root", "output_root", "scratch_root"):
            require(Path(cfg[name]).is_absolute(), f"{name} must be absolute.")
        self.reads = Path(cfg["reads_root"]).resolve(strict=True)
        self.refs = Path(cfg["reference_root"]).resolve(strict=True)
        self.out = Path(cfg["output_root"]).resolve()
        self.scratch = Path(cfg["scratch_root"]).resolve()
        roots = (self.reads, self.refs, self.out, self.scratch)
        require(all(',' not in str(path) for path in roots), "Docker bind-mount roots must not contain commas.")
        require(all(not a.is_relative_to(b) and not b.is_relative_to(a)
                    for i, a in enumerate(roots) for b in roots[i + 1:]), "Input, reference, output and scratch roots must not overlap.")
        require(isinstance(cfg["bqsr_enabled"], bool), "Set bqsr_enabled explicitly.")
        require(cfg.get("processing_assumptions"), "Document processing_assumptions.")
        require(1 <= int(cfg["cpu_limit"]) <= 30, "Use a 1..30 CPU limit on this host.")
        require(32 <= int(cfg["container_memory_gib"]) <= 180, "Use 32..180 GiB container memory on this host.")
        require(1 <= int(cfg["sort_memory_gb"]) < int(cfg["container_memory_gib"]), "Invalid sort memory limit.")
        require(0 < float(cfg["max_run_hours"]) <= 16 and
                float(cfg["max_run_hours"]) * float(cfg["hourly_usd"]) <= 40, "Processing allocation exceeds $40; allow for setup/transfer.")
        require(re.fullmatch(r"[0-9]+", str(cfg["gpu_device"])), "Choose one numeric physical GPU device; MIG is unsupported.")
        reference = cfg["reference"]
        require(all(reference.get(k) and not reference[k].startswith("SET_") for k in ("assembly_label", "source")),
                "Document the exact reference identity and source.")
        self.fasta = checked_file(self.refs, reference["fasta"], reference["sha256"])
        require(self.fasta.suffix in (".fa", ".fasta", ".fna"), "Use an uncompressed FASTA.")
        for suffix in (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"):
            require(Path(str(self.fasta) + suffix).is_file(), f"Missing BWA-MEM/reference index: {self.fasta}{suffix}")
        require(self.fasta.with_suffix(".dict").is_file(), "Missing sequence dictionary (.dict).")
        known = cfg.get("known_sites", [])
        require(bool(known) == cfg["bqsr_enabled"], "BQSR requires explicit known-sites; leave known_sites empty when BQSR is disabled.")
        self.known = [checked_file(self.refs, r["vcf"], r["sha256"]) for r in known]
        for path in self.known:
            require(any(Path(str(path) + s).is_file() for s in (".tbi", ".csi", ".idx")), f"Missing known-sites index: {path}")
        self.samples = read_tsv(cfg["samples_manifest"],
                               ("sample_id", "r1", "r2", "read_group_id", "library", "library_origin", "platform", "platform_unit"))
        for field in ("sample_id", "read_group_id"):
            require(len({s[field] for s in self.samples}) == len(self.samples), f"Duplicate {field}; one FASTQ pair per sample is supported.")
        hashes = {}
        for sample in self.samples:
            for key in ("sample_id", "read_group_id", "library"):
                require(IDENTIFIER.fullmatch(sample[key]) and not sample[key].startswith("SET_"), f"Invalid {key}.")
            require(sample["library_origin"] in ("confirmed", "processing_assumption"), "Set library_origin explicitly.")
            for key in ("platform", "platform_unit"):
                require(not sample[key].startswith("SET_") and not any(c in sample[key] for c in "\t\r\n\\"),
                        f"Invalid {key}; leave unknown values empty.")
            sample["platform_unit_original"] = sample["platform_unit"] or "unknown"
            sample["platform_unit_effective"] = sample["platform_unit"] or "analysis_unit_" + sample["read_group_id"]
            sample["platform_unit_origin"] = "manifest_supplied" if sample["platform_unit"] else "processing_assumption"
            sample["processing_assumptions"] = [] if sample["platform_unit"] else [
                "Original platform unit is unknown; effective PU is a technical input-group label, not a confirmed physical flowcell."
            ]
            pair = [contained_file(self.reads, sample[r]) for r in ("r1", "r2")]
            require(pair[0] != pair[1], "R1 and R2 must be different files.")
            for key, path in zip(("r1", "r2"), pair):
                sample[key] = path.relative_to(self.reads).as_posix()
            hashes.update({str(p): sha256(p) for p in pair})
        self.runner_sha256 = sha256(Path(__file__).resolve())
        self.base = fingerprint({"config": cfg, "samples": self.samples, "reads_sha256": hashes,
                                 "runner_sha256": self.runner_sha256})
        self.deadline = time.monotonic() + float(cfg["max_run_hours"]) * 3600

    def run(self, command, log, capture=False):
        command = list(map(str, command))
        self.commands.append(command)
        print(json.dumps({"command": command}, ensure_ascii=False), flush=True)
        if not self.execute:
            return ""
        remaining = self.deadline - time.monotonic()
        require(remaining > 0, "Processing deadline reached. VM billing continues until termination.")
        with Path(log).open("ab") as stream:
            try:
                result = subprocess.run(command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE if capture else stream,
                                        stderr=stream, check=True, timeout=remaining)
                if capture:
                    stream.write(result.stdout)
                    return result.stdout.decode("utf-8")
            except (subprocess.TimeoutExpired, KeyboardInterrupt):
                if command[:len(self.docker_prefix) + 1] == self.docker_prefix + ["run"] and "--name" in command:
                    subprocess.run(self.docker_prefix + ["stop", "--time", "60", command[command.index("--name") + 1]], check=False, timeout=90)
                raise
        return ""

    def docker(self, token, args):
        cfg = self.cfg
        return [*self.docker_prefix, "run", "--rm", "--name", "marina-" + token, "--network", "none",
                "--gpus", f"device={cfg['gpu_device']}", "--cpus", str(cfg["cpu_limit"]),
                "--memory", f"{cfg['container_memory_gib']}g", "--user", f"{os.getuid()}:{os.getgid()}",
                "--mount", f"type=bind,src={self.reads},dst=/inputs,readonly",
                "--mount", f"type=bind,src={self.refs},dst=/reference,readonly",
                "--mount", f"type=bind,src={self.out},dst=/outputs",
                "--mount", f"type=bind,src={self.scratch},dst=/scratch",
                cfg["image_digest"], "pbrun", *args]

    def stage(self, sample, kind, call=None):
        sid = sample["sample_id"]
        label = kind if call is None else "call-" + call["callset"]
        target = self.out / sid / label
        stamp = fingerprint({"base": self.base, "stage": label, "call": call})
        if target.exists():
            require(receipt_valid(target, stamp), f"Existing stage {target} changed or is incomplete; refusing to overwrite.")
            print(f"Verified completed stage: {sid}/{label}")
            return
        token = uuid.uuid4().hex[:10]
        partial = target.with_name(target.name + ".partial-" + token)
        scratch = self.scratch / sid / (label + "-" + token)
        output = "/outputs/" + partial.relative_to(self.out).as_posix()
        tmp = "/scratch/" + scratch.relative_to(self.scratch).as_posix()
        log = partial / "stage.log"
        if self.execute:
            partial.mkdir(parents=True, exist_ok=False)
            scratch.mkdir(parents=True, exist_ok=False)
            require(shutil.disk_usage(self.out).free >= int(self.cfg["free_disk_gib_min"]) * 1024**3,
                    "Less than configured free disk reserve; stage not started.")
        ref = "/reference/" + self.fasta.relative_to(self.refs).as_posix()
        if kind == "align":
            rg = f"@RG\\tID:{sample['read_group_id']}\\tSM:{sid}\\tLB:{sample['library']}"
            if sample["platform"]:
                rg += "\\tPL:" + sample["platform"]
            rg += "\\tPU:" + sample["platform_unit_effective"]
            args = ["fq2bam", "--ref", ref, "--in-fq", "/inputs/" + sample["r1"], "/inputs/" + sample["r2"], rg,
                    "--out-bam", output + "/marked.bam", "--out-duplicate-metrics", output + "/duplicate-metrics.txt",
                    "--memory-limit", str(self.cfg["sort_memory_gb"]), "--bwa-primary-cpus", "16", "--bwa-cpu-thread-pool", "1",
                    "--bwa-options=-K 10000000", "--num-gpus", "1", "--tmp-dir", tmp,
                    "--monitor-usage", "--logfile", output + "/parabricks.log"]
            if self.cfg["bqsr_enabled"]:
                for site in self.known:
                    args += ["--knownSites", "/reference/" + site.relative_to(self.refs).as_posix()]
                args += ["--out-recal-file", output + "/bqsr.table"]
        else:
            aligned = self.out / sid / "align"
            if self.execute:
                require(receipt_valid(aligned, fingerprint({"base": self.base, "stage": "align", "call": None})), "Run/verify alignment first.")
            args = ["haplotypecaller", "--ref", ref, "--in-bam", "/outputs/" + (aligned / "marked.bam").relative_to(self.out).as_posix(),
                    "--gvcf", "--out-variants", output + "/variants.g.vcf.gz", "--num-gpus", "1", "--ploidy", call["ploidy"],
                    "--interval-file", "/reference/" + call["intervals"], "--num-htvc-threads", "8", "--tmp-dir", tmp,
                    "--logfile", output + "/parabricks.log"]
            if self.cfg["bqsr_enabled"]:
                args += ["--in-recal-file", "/outputs/" + (aligned / "bqsr.table").relative_to(self.out).as_posix()]
        start = dt.datetime.now(dt.timezone.utc).isoformat()
        self.commands = []
        processing_metadata = {"original_platform_unit": sample["platform_unit_original"],
                               "effective_platform_unit": sample["platform_unit_effective"],
                               "platform_unit_origin": sample["platform_unit_origin"],
                               "processing_assumptions": self.cfg["processing_assumptions"] + sample["processing_assumptions"]}
        print(json.dumps({"sample": sid, "processing_metadata": processing_metadata}, ensure_ascii=False), flush=True)
        self.run(self.docker(token, args), log)
        if kind == "align":
            self.run(["samtools", "quickcheck", "-v", partial / "marked.bam"], log)
            self.run(["samtools", "index", "-@", "4", partial / "marked.bam"], log)
            for qc in ("flagstat", "idxstats", "stats"):
                self.run(["samtools", qc, partial / "marked.bam"], partial / (qc + ".txt"))
        else:
            gvcf = partial / "variants.g.vcf.gz"
            self.run(["bcftools", "index", "--force", "--tbi", gvcf], log)
            samples = self.run(["bcftools", "query", "-l", gvcf], log, capture=True)
            require(not self.execute or samples.splitlines() == [sid], "Unexpected gVCF sample name.")
            self.run(["bcftools", "stats", gvcf], partial / "vcf-stats.txt")
        if self.execute:
            if kind == "align":
                require((partial / "duplicate-metrics.txt").is_file(), "Missing duplicate metrics.")
                require(not self.cfg["bqsr_enabled"] or (partial / "bqsr.table").is_file(), "Missing BQSR report.")
            receipt = {"fingerprint": stamp, "sample": sid, "image_tag": TAG, "config": self.cfg, "commands": self.commands,
                       "runner_sha256": self.runner_sha256,
                       "original_platform": sample["platform"] or "unknown", "library_origin": sample["library_origin"],
                       **processing_metadata,
                       "started_utc": start, "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                       "outputs": {p.name: sha256(p) for p in partial.iterdir() if p.is_file()}}
            (partial / "receipt.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            partial.rename(target)
            print(f"Completed and hashed: {target}. Back up this stage off the ephemeral VM.")

    def process(self, stage, requested):
        selected = [s for s in self.samples if not requested or s["sample_id"] in requested]
        require(selected and (not requested or set(requested) == {s['sample_id'] for s in selected}), "Unknown sample selection.")
        calls = []
        if stage in ("call", "all"):
            calls = read_tsv(self.cfg["calls_manifest"], ("sample_id", "callset", "intervals", "intervals_sha256", "ploidy"))
            require(len({(c["sample_id"], c["callset"]) for c in calls}) == len(calls), "Duplicate callset rows.")
            for call in calls:
                require(call["sample_id"] in {s["sample_id"] for s in self.samples}, "Unknown sample in calls manifest.")
                require(IDENTIFIER.fullmatch(call["callset"]) and not call["callset"].startswith("SET_"), "Set callset names explicitly.")
                require(call["ploidy"] in ("1", "2"), "Specify ploidy 1 or 2 explicitly; no sex is inferred.")
                path = checked_file(self.refs, call["intervals"], call["intervals_sha256"])
                call["intervals"] = path.relative_to(self.refs).as_posix()
            require(all(any(c["sample_id"] == s["sample_id"] for c in calls) for s in selected), "Missing call intervals/ploidy for a selected sample.")
        if self.execute:
            require(sys.platform == "linux", "Execute only on the Linux GPU server.")
            for tool in ("docker", "samtools", "bcftools"):
                require(shutil.which(tool), f"Missing {tool}; bootstrap first.")
            docker_access = subprocess.run(["docker", "info"], stdin=subprocess.DEVNULL,
                                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=30)
            if docker_access.returncode:
                require(shutil.which("sudo"), "Docker is not accessible; sudo -n docker is unavailable.")
                self.docker_prefix = ["sudo", "-n", "docker"]
                subprocess.run(self.docker_prefix + ["info"], check=True, stdin=subprocess.DEVNULL,
                               stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=30)
                print("Using sudo -n docker; container outputs retain the invoking user's UID:GID.")
            self.out.mkdir(parents=True, exist_ok=True)
            self.scratch.mkdir(parents=True, exist_ok=True)
            import fcntl
            self.lock = (self.out / ".gpu-queue.lock").open("a")
            try:
                fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise ValueError("Another GPU queue is running for this output root.") from None
            ids = [subprocess.check_output(self.docker_prefix + ["image", "inspect", name, "--format", "{{.Id}}"], text=True).strip()
                   for name in (TAG, self.cfg["image_digest"])]
            require(ids[0] == ids[1], "Digest does not match the locally pulled 4.8.0-1 image.")
        for sample in selected:
            if stage in ("align", "all"):
                self.stage(sample, "align")
            if stage in ("call", "all"):
                for call in calls:
                    if call["sample_id"] == sample["sample_id"]:
                        self.stage(sample, "call", call)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--stage", choices=("align", "call", "all"), default="align")
    parser.add_argument("--sample", action="append", help="Repeat for selected samples; omitted means all manifest samples.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--execute", action="store_true", help="Launch on the Linux GPU server.")
    mode.add_argument("--dry-run", action="store_true", help="Validate real inputs and print commands; this is the default.")
    args = parser.parse_args()
    try:
        os.umask(0o077)
        Queue(args.config, args.execute).process(args.stage, args.sample)
        if not args.execute:
            print("DRY RUN: no containers, uploads, outputs, or analyses were started.")
    except (ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
