#!/usr/bin/env python3
"""Prepare a pinned Broad GRCh38 bundle; local plan unless --execute."""
import argparse
import base64
from collections import OrderedDict
import ctypes
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import urllib.parse
import urllib.request
import uuid


FASTA = "Homo_sapiens_assembly38.fasta"
DICT = "Homo_sapiens_assembly38.dict"
BUNDLE = "broad-hg38-v0-altaware-g1575676516681666"
INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
ANCHORS = {"chr1": 248956422, "chr2": 242193529, "chrX": 156040895,
           "chrY": 57227415, "chrM": 16569}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path, algorithm="sha256"):
    result = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result


def source_url(obj):
    return ("https://storage.googleapis.com/storage/v1/b/" + obj["bucket"] +
            "/o/" + urllib.parse.quote(obj["object"], safe="") +
            "?alt=media&generation=" + obj["generation"])


def download(obj, directory):
    target = directory / obj["destination"]
    partial = directory / (obj["destination"] + ".download")
    require(not target.exists() and not partial.exists(), "Refusing download overwrite")
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    total = 0
    started = time.monotonic()
    last_progress = 0
    print(f"Downloading {obj['destination']}: {obj['bytes']} bytes, generation {obj['generation']}", flush=True)
    with urllib.request.urlopen(source_url(obj), timeout=120) as source, partial.open("xb") as out:
        require(source.status == 200, "Unexpected download HTTP status")
        while block := source.read(8 * 1024 * 1024):
            total += len(block)
            require(total <= obj["bytes"], "Download larger than locked object")
            md5.update(block)
            sha256.update(block)
            out.write(block)
            if total - last_progress >= 128 * 1024 * 1024:
                print(f"Download progress {obj['destination']}: {total}/{obj['bytes']} bytes", flush=True)
                last_progress = total
        out.flush()
        os.fsync(out.fileno())
    require(total == obj["bytes"], "Downloaded object size differs from lock")
    require(base64.b64encode(md5.digest()).decode("ascii") == obj["md5_base64"],
            "Downloaded object MD5 differs from publisher metadata")
    partial.rename(target)
    print(f"Downloaded and hash-verified {obj['destination']} in {time.monotonic() - started:.1f}s", flush=True)
    return {**obj, "url": source_url(obj), "sha256": sha256.hexdigest()}


def load_fai(path):
    records = OrderedDict()
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split("\t")
        require(len(fields) == 5 and fields[0] not in records, "Malformed/duplicate FAI row")
        records[fields[0]] = tuple(int(value) for value in fields[1:])
    return records


def load_dictionary(path):
    records = OrderedDict()
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:])
        require(fields["SN"] not in records and re.fullmatch(r"[0-9a-fA-F]{32}", fields["M5"]),
                "Malformed/duplicate dictionary row")
        records[fields["SN"]] = (int(fields["LN"]), fields["M5"].lower())
    return records


def validate_fasta(directory):
    """Check every FASTA base against source dictionary M5 and every FAI offset."""
    fai = load_fai(directory / (FASTA + ".fai"))
    dictionary = load_dictionary(directory / DICT)
    require(list(fai) == list(dictionary), "FAI and dictionary contig order differ")
    require(len(fai) == 3366, "Unexpected Broad hg38/v0 contig count")
    for name, length in ANCHORS.items():
        require(name in fai and fai[name][0] == length, "Unexpected GRCh38 anchor " + name)
    observed = []
    name = None
    md5 = None
    length = 0
    first_offset = first_bases = first_width = None
    prior_bases = prior_width = None
    probe = None
    probe_tail = b""

    def finish():
        if name is None:
            return
        require(length == dictionary[name][0] and md5.hexdigest() == dictionary[name][1],
                "FASTA length or sequence MD5 differs for " + name)
        require((length, first_offset, first_bases, first_width) == fai[name],
                "FASTA layout differs from FAI for " + name)
        observed.append(name)

    with (directory / FASTA).open("rb") as stream:
        while True:
            offset = stream.tell()
            raw = stream.readline()
            if not raw:
                break
            if raw.startswith(b">"):
                finish()
                name = raw[1:].split()[0].decode("ascii")
                require(name in dictionary and name not in observed, "Unexpected/duplicate FASTA contig")
                md5, length = hashlib.md5(), 0
                first_offset = first_bases = first_width = None
                prior_bases = prior_width = None
                continue
            require(name is not None, "Sequence before FASTA header")
            sequence = raw.rstrip(b"\r\n").upper()
            require(sequence and re.fullmatch(rb"[ACGTRYSWKMBDHVN]+", sequence), "Invalid FASTA bases")
            if first_offset is None:
                first_offset, first_bases, first_width = offset, len(sequence), len(raw)
            if prior_bases is not None:
                require(prior_bases == first_bases and prior_width == first_width,
                        "Short/internal variable-width FASTA sequence line")
            require(len(sequence) <= first_bases, "FASTA line longer than first line")
            md5.update(sequence)
            length += len(sequence)
            prior_bases, prior_width = len(sequence), len(raw)
            if name == "chr1" and probe is None:
                probe_window = probe_tail + sequence
                probe_match = re.search(rb"[ACGT]{100}", probe_window)
                if probe_match:
                    probe = probe_match.group(0)
                probe_tail = probe_window[-99:]
        finish()
    require(observed == list(fai), "FASTA and FAI contig order differ")
    require(sum(row[0] for row in fai.values()) == 3217346917, "Unexpected total reference bases")
    require(probe is not None, "No public-reference smoke-test sequence found")
    alt_count = alt_mapped = alt_unmapped = 0
    for line in (directory / (FASTA + ".alt")).read_text(encoding="ascii").splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        require(len(fields) >= 6 and fields[0] in fai, "ALT sidecar query contig absent from locked FASTA")
        flag = int(fields[1])
        if flag & 4:
            require(fields[2] == "*" and fields[3] == "0", "Malformed unmapped ALT SAM row")
            alt_unmapped += 1
        else:
            require(fields[2] in fai and 1 <= int(fields[3]) <= fai[fields[2]][0],
                    "Mapped ALT target contig/coordinate absent from locked FASTA")
            alt_mapped += 1
        alt_count += 1
    return fai, probe, {"contigs": len(fai), "total_bases": 3217346917,
                        "alt_sidecar_rows": alt_count, "alt_mapped_rows": alt_mapped,
                        "alt_unmapped_rows": alt_unmapped, "sequence_md5_all_contigs": "PASS"}


def promote_without_overwrite(stage, final):
    """Ubuntu/Linux atomic directory rename with kernel-enforced no-replace."""
    libc = ctypes.CDLL(None, use_errno=True)
    require(hasattr(libc, "renameat2"), "Atomic no-replace rename unavailable")
    rename = libc.renameat2
    rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    rename.restype = ctypes.c_int
    if rename(-100, os.fsencode(stage), -100, os.fsencode(final), 1) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), str(final))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--root", default="/srv/marina/references")
    parser.add_argument("--index-mode", choices=("prebuilt", "rebuild"), default="prebuilt")
    parser.add_argument("--bwa-binary", default="/usr/bin/bwa")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    lock_path = Path(args.lock).resolve(strict=True)
    lock_bytes = lock_path.read_bytes()
    lock = json.loads(lock_bytes)
    require(lock.get("schema_version") == 1 and lock.get("bwa_version") == "0.7.17", "Unexpected lock")
    require(re.fullmatch(r"nvcr\.io/nvidia/clara/clara-parabricks@sha256:[0-9a-f]{64}", lock["image"]),
            "Use an immutable official Parabricks image reference")
    expected_names = {FASTA, FASTA + ".fai", DICT, FASTA + ".alt"} | {FASTA + s for s in INDEX_SUFFIXES}
    require({obj["destination"] for obj in lock["objects"]} == expected_names and len(lock["objects"]) == 9,
            "Unexpected reference source list")
    for obj in lock["objects"]:
        require(obj["bucket"] == "gcp-public-data--broad-references" and obj["object"].startswith("hg38/v0/"),
                "Unexpected public reference source")
        require(obj["generation"].isdigit() and obj["bytes"] > 0, "Unpinned source metadata")
    root = Path(args.root)
    require(root.is_absolute() and ".." not in root.parts and len(root.parts) >= 4,
            "Use a dedicated absolute reference root")
    final = root / (BUNDLE + "-" + args.index_mode)
    selected = [o for o in lock["objects"] if args.index_mode == "prebuilt" or
                o["destination"] not in {FASTA + s for s in INDEX_SUFFIXES}]
    print(json.dumps({"mode": "execute" if args.execute else "dry-run", "destination": str(final),
                      "download_bytes": sum(o["bytes"] for o in selected),
                      "source_objects": [source_url(o) for o in selected],
                      "image": lock["image"], "index_mode": args.index_mode,
                      "bwa_index": "published Broad .fasta.64.* indexes, verified with BWA 0.7.17" if args.index_mode == "prebuilt" else "bwa 0.7.17 index -a bwtsw",
                      "alt_policy": lock["alt_policy"]}, indent=2), flush=True)
    if not args.execute:
        return
    require(os.name == "posix", "Execution requires the Linux analysis environment")
    bwa_path = Path(args.bwa_binary).resolve(strict=True)
    require(bwa_path.is_file() and os.access(bwa_path, os.X_OK), "Install official Ubuntu BWA 0.7.17 first")
    print("Checking host BWA 0.7.17 and package identity before any reference download", flush=True)
    version = subprocess.run([str(bwa_path)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
    require(re.search(r"Version:\s*0\.7\.17(?:\W|$)", version.stdout),
            "Expected host BWA 0.7.17; do not substitute another version silently")
    package_version = subprocess.run(["dpkg-query", "-W", "-f=${Package} ${Version} ${Architecture}\\n", "bwa"],
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=15)
    require(package_version.returncode == 0 and package_version.stdout.startswith("bwa 0.7.17"),
            "Expected official Ubuntu BWA package metadata")
    print("Host BWA version/package checks passed", flush=True)
    root.mkdir(parents=True, mode=0o700, exist_ok=True)
    require(not root.is_symlink(), "Reference root must not be a symlink")
    with (root / ".prepare-reference.lock").open("a") as lock_file:
        fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if final.exists():
            receipt = json.loads((final / "reference-receipt.json").read_text())
            require(receipt["lock_sha256"] == hashlib.sha256(lock_bytes).hexdigest(),
                    "Existing reference uses another lock; refusing overwrite")
            require(receipt["index_mode"] == args.index_mode, "Existing reference uses another index mode")
            for artifact in receipt["artifacts"]:
                p = final / artifact["path"]
                require(p.is_file() and p.stat().st_size == artifact["bytes"] and
                        digest(p).hexdigest() == artifact["sha256"], "Existing reference changed")
            print("Existing completed reference passed all recorded SHA256 checks; skipped.")
            return
        require(shutil.disk_usage(root).free >= 25 * 2**30, "Need at least 25 GiB free for reference preparation")
        stage = root / (".partial-reference-" + uuid.uuid4().hex)
        stage.mkdir(mode=0o700)
        (stage / "source-lock.json").write_bytes(lock_bytes)
        (stage / "bwa-version.txt").write_text(version.stdout, encoding="utf-8")
        (stage / "bwa-package.txt").write_text(package_version.stdout, encoding="utf-8")
        sources = [download(obj, stage) for obj in selected]
        print("Validating full FASTA sequence MD5 values, FAI layout, dictionary and ALT mapping", flush=True)
        fai, probe, validation = validate_fasta(stage)
        print("Full FASTA/contig validation passed: " + json.dumps(validation), flush=True)
        if args.index_mode == "rebuild":
            with (stage / "bwa-index.log").open("w") as log:
                subprocess.run([str(bwa_path), "index", "-a", "bwtsw", str(stage / FASTA)],
                               check=True, stdout=log, stderr=subprocess.STDOUT, timeout=7200)
        for suffix in INDEX_SUFFIXES:
            require((stage / (FASTA + suffix)).is_file() and (stage / (FASTA + suffix)).stat().st_size > 0,
                    "Missing newly generated BWA index " + suffix)
        ann = (stage / (FASTA + ".ann")).read_text(encoding="ascii").splitlines()
        header = ann[0].split()
        require(int(header[0]) == 3217346917 and int(header[1]) == len(fai), "BWA ANN header mismatch")
        require([ann[i].split()[1] for i in range(1, len(ann), 2)] == list(fai), "BWA ANN contig order mismatch")
        require([int(ann[i].split()[1]) for i in range(2, len(ann), 2)] == [r[0] for r in fai.values()],
                "BWA ANN contig lengths mismatch")
        (stage / "reference-probe.fastq").write_bytes(b"@reference_probe\n" + probe + b"\n+\n" + b"I" * len(probe) + b"\n")
        print("Running public-reference BWA index smoke test (timeout 300s)", flush=True)
        with (stage / "index-smoke-test.sam").open("w") as output, (stage / "index-smoke-test.log").open("w") as log:
            subprocess.run([str(bwa_path), "mem", "-t", "2", str(stage / FASTA),
                            str(stage / "reference-probe.fastq")], check=True, stdout=output, stderr=log, timeout=300)
        sam_records = [line.split("\t") for line in (stage / "index-smoke-test.sam").read_text().splitlines()
                       if line and not line.startswith("@")]
        require(any(not int(row[1]) & 4 and row[2] in fai for row in sam_records), "Fresh index failed smoke alignment")
        print("Index smoke test passed; calculating final artifact SHA256 values", flush=True)
        artifacts = [{"path": p.name, "bytes": p.stat().st_size, "sha256": digest(p).hexdigest()}
                     for p in sorted(stage.iterdir()) if p.is_file()]
        receipt = {"schema_version": 1, "assembly_label": lock["assembly_label"],
                   "lock_sha256": hashlib.sha256(lock_bytes).hexdigest(), "image": lock["image"],
                   "bwa_binary": str(bwa_path), "bwa_binary_sha256": digest(bwa_path).hexdigest(),
                   "bwa_version": "0.7.17", "bwa_package": package_version.stdout.strip(), "sources": sources,
                   "index_mode": args.index_mode,
                   "validation": validation, "index_smoke_test": "PASS", "artifacts": artifacts,
                   "completed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   "old_vcf_reference_identity": "Not asserted; new uniform reanalysis bundle"}
        (stage / "reference-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        require(not final.exists(), "Destination appeared during preparation; refusing overwrite")
        promote_without_overwrite(stage, final)
        print(json.dumps({"completed": str(final), "fasta_sha256": digest(final / FASTA).hexdigest(),
                          "validation": validation}, indent=2), flush=True)


if __name__ == "__main__":
    main()
