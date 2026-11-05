# GPU server setup and staged analysis

The chosen target is Brev / Lambda `gpu_1x_a100_sxm4`: one A100 40 GB, 30 vCPU,
200 GiB RAM, 512 GiB local SSD. The console quote was $2.39/hour including disk
on 2026-11-05. It is marked pre-release, has **no stop/start**, and Brev says an
exhausted credit balance causes deletion of the VM and all its data. Copy verified
outputs back after each sample. The scripts do not rent, terminate, or top up a VM.

Docker is the container engine; NVIDIA Container Toolkit makes its GPU available.
For four samples, a resumable single-GPU queue is sufficient orchestration. Kubernetes
and a second scheduler are not needed. GPU stages run one at a time; input/QC/backup
steps may be scheduled around them. Four simultaneous full pipelines on the same GPU
are deliberately avoided.

## Install after SSH is available

These commands run **on the Ubuntu GPU VM**, from a checkout of this repository.
Use a normal sudo-capable provider account; do not paste SSH keys into the repo.

```bash
sudo bash scripts/bootstrap-server.sh /srv/marina
sudo bash scripts/preflight-server.sh /srv/marina
```

Bootstrap preserves an existing provider Docker installation, uses the official
Docker/NVIDIA apt repositories when installation is needed, and pins a newly installed
NVIDIA Container Toolkit to `1.20.1-1`. It never replaces the GPU driver. If GPU runtime
configuration requires restarting Docker, active containers make it stop for review.
An unavailable pinned package produces an error rather than silently falling back.
Ubuntu 22.04, 24.04 and 26.04 x86_64 are supported by this script.

Parabricks is pinned to `nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1`. Bootstrap
records package versions and image inspection under `/srv/marina/logs`. Put the
actual pulled image's `RepoDigests[0]` in the private run configuration. The runner
uses this immutable digest and checks it matches the locally pulled version tag.
Provider Docker and host QC packages are recorded as installed, not silently described
as reproducibly pinned versions. `DOCKER_VERSION` can pin a fresh Docker installation.

Run the Python queue as the normal VM user, **without `sudo python3`**. If that user
cannot access the Docker daemon directly, the runner uses `sudo -n docker` for daemon
commands and timeout cleanup. Containers still receive the original user's numeric
UID:GID, so output files remain readable/writable by that user. A provider account
without either direct Docker access or non-interactive sudo fails before analysis.
Bootstrap does not add users to the root-equivalent `docker` group or change SSH access.

The A100 is on NVIDIA's tested list. Current Parabricks 4.8 documentation requires a
CUDA-13.3-compatible driver, Docker >=20.10 and NVIDIA Container Toolkit. Checking
`nvidia-smi` inside a container confirms GPU exposure, not a completed analysis or a
guarantee every CUDA workload will run. The first pilot confirms real compatibility.
Preflight checks at least 24 available CPU threads, 100 GB total host RAM, 36,000 MiB
GPU memory, 150 GiB free disk and absence of a MIG device.

## Private data/config contract

The transfer scripts use this layout beneath their `RemoteRoot` (use `/srv/marina`
to match these examples):

```text
/srv/marina/
  inputs/raw/          # Eight original R1/R2 FASTQ.gz; no duplicate copy
  inputs/legacy/       # Previous VCF and report files, preserved as originals
  inputs/provenance/   # Local audit/plan, outside git
  references/          # Exact FASTA, indexes, known-sites if enabled, interval files
  private/             # run.json, samples.tsv, calls.tsv; never commit
  outputs/             # Atomic completed stage directories and checksummed receipts
  scratch/             # Per-stage temporary data
  logs/                # Installation/preflight logs
```

Copy `config/run.example.json`, `config/samples.example.tsv` and
`config/calls.example.tsv` into `/srv/marina/private/` as `run.json`, `samples.tsv`
and `calls.tsv`. For local preparation, use the repository's ignored `.local/`
directory. Placeholder values are rejected. Real identifiers and biological results
belong in private manifests, not committed templates.

`samples.tsv` has eight TSV columns:

| Column | Meaning |
|---|---|
| `sample_id` | Explicit chosen sample label used as BAM `SM` and gVCF sample name |
| `r1`, `r2` | Relative filenames beneath `reads_root` |
| `read_group_id` | Unique technical group ID |
| `library` | Explicit library grouping ID |
| `library_origin` | `confirmed` or `processing_assumption` |
| `platform`, `platform_unit` | Confirmed values, or empty when unknown |

The approved processing assumption for this data is **one available FASTQ pair = one
analytical library**. Identifiers such as `analysis_library_01` express that grouping;
they do not assert a confirmed wet-lab library. The templates label this origin.
Unknown platform metadata stays unknown; no `ILLUMINA` tag is invented.

`calls.tsv` has five TSV columns: `sample_id`, `callset`, `intervals`,
`intervals_sha256`, `ploidy`. The interval file is relative to `reference_root`.
Use separate non-overlapping callsets when ploidy differs. Ploidy must explicitly be
`1` or `2`; the code does not guess sex, PAR boundaries or mitochondrial heteroplasmy.
An autosomal-only pilot is permissible when labelled as such. Do not describe it as
an all-chromosome final callset. Sex chromosome/MT processing requires its own policy.

## Reference and BQSR choices

There is no reference FASTA or known-sites bundle in the supplied data. Choose one
fixed official GRCh38 bundle for the **new uniform reanalysis**, record its source,
assembly label, download date and SHA256, and use it for every sample. The old VCF
annotation label `GRCh38.105` does not identify the exact alignment FASTA. A new bundle
does not silently reproduce the old laboratory reference.

Before running analysis, place an uncompressed `.fa`/`.fasta`/`.fna` under
`reference_root`, together with `.fai`, its `.dict`, and BWA-MEM `.amb`, `.ann`, `.bwt`,
`.pac`, `.sa` indexes. These are **BWA-MEM indexes**, not BWA-MEM2 indexes. Record the
FASTA SHA256 in `run.json`; the runner verifies it. Interval files also require their
real SHA256 in `calls.tsv`. Exact input read hashes are computed into the run fingerprint.
Index identity/contig agreement should additionally be checked during the pilot.

`bqsr_enabled` defaults to **false**. With the available unknown platform provenance,
this is an explicit conservative processing choice, not a claim BQSR was done.
To enable it, review the choice, set `bqsr_enabled=true`, and provide same-reference
known-sites entries `{"vcf":"relative.vcf.gz","sha256":"..."}` with `.tbi`/`.csi`/`.idx`.
The runner will refuse enabled BQSR without known-sites. When enabled, `fq2bam`
generates the report; HaplotypeCaller applies it via `--in-recal-file`. The stored
`marked.bam` is the marked BAM, not a promised permanently BQSR-rewritten BAM.

## Pilot and queue

Validate config, input hashes and planned commands first. These commands do not launch
containers unless `--execute` is present. They target Linux Python >=3.9.

```bash
python3 scripts/run_parabricks.py --config /srv/marina/private/run.json --stage all --sample SAMPLE_01 --dry-run
tmux new -s marina
python3 scripts/run_parabricks.py --config /srv/marina/private/run.json --stage align --sample SAMPLE_01 --execute
python3 scripts/run_parabricks.py --config /srv/marina/private/run.json --stage call --sample SAMPLE_01 --execute
```

Inspect duplicate/flagstat/idxstats/stats and gVCF statistics, time, resource use and
read/sample identity. Verify the downloaded outputs and resolve differences before
expanding to the four samples. A low-depth dataset remains low-depth after rerunning.
No clinical interpretation or normal/negative conclusion is implied by this pipeline.

After the pilot, a command without `--sample` queues all manifest samples:

```bash
python3 scripts/run_parabricks.py --config /srv/marina/private/run.json --stage all --execute
```

Successful stages are atomically renamed from `.partial-*` directories only after
basic output checks and SHA256 receipts. Repeating the exact config skips a stage
only if its fingerprint and output checksums still match. Changed/incomplete final
directories cause a stop instead of overwrite. Failed partial directories remain for
inspection; review their disk use before rerunning. Do not edit completed stage logs.
One process lock prevents duplicate queues for the same output root. Always use that
same root; changing it is not a method to run another concurrent GPU queue.

Default limits are 24 CPU threads/160 GiB container RAM, 100 GB sorting limit, and one
GPU. `max_run_hours` bounds a single runner invocation; timed-out Docker workloads are
stopped. It **does not stop cloud billing or enforce an account spending cap**, and it
does not include VM setup/upload time. Repeated invocations each get a new deadline.
The default 12 hours represents $28.68 of processing, not a measured runtime estimate.
Retain the rest of the $50 for installation, transfer, verification and contingency.

This code ends at marked BAM plus per-sample/per-ploidy gVCF and technical statistics.
Cohort joint genotyping, harmonized filters, annotation, CNV/SV/MT analysis and clinical
interpretation are subsequent explicit stages. Do not label the gVCFs as final filtered
clinical VCFs. No SQL dump is executed or imported by these scripts.

## Official documentation consulted

- [Parabricks 4.8 requirements](https://docs.nvidia.com/clara/parabricks/get-started/installation-requirements)
- [fq2bam options](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/fq2bam)
- [HaplotypeCaller options and mismatches](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/haplotypecaller)
- [NVIDIA Container Toolkit installation](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)
- [Docker Engine installation on Ubuntu](https://docs.docker.com/engine/install/ubuntu/)
