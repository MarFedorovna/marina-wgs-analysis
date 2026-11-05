# Broad GRCh38 reference preparation

Prepared and validated on the selected VM on 2026-11-05. The complete public reference and matching index set were downloaded, verified and atomically promoted in **193.6 seconds**, measured from staged source-lock creation to the completed receipt. This is an observed preparation time, not a promise for future connections or patient analysis.

## Exact resource and purpose

Use the Broad Institute public `hg38/v0/Homo_sapiens_assembly38.fasta` from `gs://gcp-public-data--broad-references`. This is a **new uniform reanalysis reference**, not a claim to reproduce the previous laboratory's FASTA or annotated VCFs. The Broad resource bundle documentation identifies this bucket for GRCh38 WGS short-variant workflows:

- https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
- https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/resources_hg38.json

The second official Broad repository explicitly lists the matching `.fasta.64.alt/.amb/.ann/.bwt/.pac/.sa` resources. The script's default is to download that **complete .64 set**, rename each destination prefix to `Homo_sapiens_assembly38.fasta.*`, and validate it before use. File bytes are unchanged. It does not mix the differently sized non-.64 BWT/SA objects with this set. The supplied BWA indexes are classic BWA-MEM indexes, not BWA-MEM2 indexes.

The source lock records each object name, GCS generation, size, publisher MD5 and CRC32C. Downloads request the locked generation over HTTPS and verify exact byte count and MD5. A new local SHA256 is computed for each downloaded/generated artifact. Publisher SHA256 values were not available and are not invented.

## Locked download sizes

| Destination | Bytes |
| --- | ---: |
| `Homo_sapiens_assembly38.fasta` | 3,249,912,778 |
| `Homo_sapiens_assembly38.fasta.fai` | 160,928 |
| `Homo_sapiens_assembly38.dict` | 581,712 |
| `Homo_sapiens_assembly38.fasta.alt` | 487,553 |
| `Homo_sapiens_assembly38.fasta.amb` | 20,199 |
| `Homo_sapiens_assembly38.fasta.ann` | 455,474 |
| `Homo_sapiens_assembly38.fasta.bwt` | 3,217,347,004 |
| `Homo_sapiens_assembly38.fasta.pac` | 804,336,731 |
| `Homo_sapiens_assembly38.fasta.sa` | 1,608,673,512 |
| **Total** | **8,881,975,891 (8.27 GiB)** |

FASTA object generation: `1575676516681666`. Publisher FASTA MD5 in base64: `f/E0lT3MqMiZdFO7uAtrXg==`. The complete lock is `config/broad-hg38-v0.sources.lock.json`. The downloaded FASTA SHA256 is `93157a161863464c9435062fd67c173fdaf99cb8b32f1455018361387ffa5564`.

Small public FAI/dictionary/ALT payloads were read locally to check their publisher hashes and metadata. They describe **3,366 contigs and 3,217,346,917 bases**. GRCh38 anchors: chr1 248,956,422; chr2 242,193,529; chrX 156,040,895; chrY 57,227,415; chrM 16,569. The ALT sidecar has 3,177 non-header SAM rows: 792 mapped and 2,385 unmapped. The unmapped rows legitimately use SAM flag 4, target `*` and position 0; they are preserved unchanged. These are public-reference characteristics, not patient metadata.

## Validation and promotion

Preparation runs under a reference-root process lock in a unique `.partial-reference-*` directory. It scans the complete FASTA, verifies all contig sequence MD5 values against the frozen dictionary, and checks FAI lengths, offsets, wrapping and contig order. It validates all ALT query contigs and mapped target contigs, accepts the source's legitimate unmapped SAM rows, and checks BWA ANN sequence order/lengths. A public-reference synthetic read then tests that the pinned BWA binary can load the index and produce a mapped alignment. This is an index smoke check, not patient analysis or an accuracy benchmark.

The immutable container is:

```text
nvcr.io/nvidia/clara/clara-parabricks@sha256:44206e2b522e87370c1107c774f52c5859af9b8f3408aa7232de4bbe3c24953c
```

The stock Parabricks image does not provide a normal `bwa` CLI. Reference index validation instead uses the host's official Ubuntu BWA **0.7.17** package at `/usr/bin/bwa`, consistent with the aligner's documented baseline. The host binary version and `dpkg-query` package identity are checked **before any reference download**, and the binary SHA256/package version are recorded in the final receipt. The [project setup wrapper](../scripts/setup_project.sh) installs the official Ubuntu BWA package before this stage. For standalone preparation, install it with `sudo apt-get update` and `sudo apt-get install --no-install-recommends bwa`; confirm it reports 0.7.17. The reference script itself does not install it or silently substitute another version. The pinned Parabricks image above remains the downstream analysis image; HaplotypeCaller baseline remains GATK 4.3.0.0. Published dictionaries/FAI are verified rather than regenerated with a different GATK version.

After success, all resource hashes, source identities, downstream image identity, host BWA identity, validation summary and command logs are retained. Linux `renameat2(RENAME_NOREPLACE)` atomically promotes the staged directory, refusing any existing target. Repeat execution verifies every recorded artifact SHA256 and skips only an unchanged completed bundle. Errors preserve the uniquely named partial directory for inspection and never replace a completed bundle. The script does not automatically remove failed partials. It prints object start/end and download progress, announces full FASTA validation, and limits the index smoke alignment to 300 seconds.

ALT policy is explicit: the Broad published ALT map is present as `.fasta.alt` for ALT-aware BWA handling. It is not a primary-assembly-only reference and does not strip ALT/HLA/decoy contigs or change any sequence bases. Call intervals and sex/ploidy policy remain separate explicit decisions; this reference-preparation script does not infer them. BQSR and known-sites files remain unconfigured and disabled by the existing run configuration.

## Commands

Run from the repository checkout. The default dry run reads the lock and prints a plan without a network connection, Docker invocation or reference-root creation:

```bash
python3 scripts/prepare_reference.py \
  --lock config/broad-hg38-v0.sources.lock.json \
  --root /srv/marina/references
```

Only after approving the source/ALT policy and checking the installed host BWA package:

```bash
python3 scripts/prepare_reference.py \
  --lock config/broad-hg38-v0.sources.lock.json \
  --root /srv/marina/references \
  --execute
```

Default completed directory:

```text
/srv/marina/references/broad-hg38-v0-altaware-g1575676516681666-prebuilt/
```

Keep the existing runner's `reference_root=/srv/marina/references` and set its private `reference.fasta` to that directory's relative FASTA path. Populate `reference.sha256` from the completed receipt, and record the exact assembly label/source. No private run configuration is edited by this script.

## Timing estimate and fallback

For the default 8.88 GB download, ideal transfer time is approximately 12 minutes at 100 Mbit/s or 2.4 minutes at 500 Mbit/s, before overhead. Full sequence/hash/index validation adds several minutes depending on disk/CPU. **Budget 10–30 minutes on a healthy cloud connection**, not a guaranteed completion time. Slow network, unavailable generations, a binary mismatch or an index-validation error stops progress. This adds public reference network traffic from GCS to the VM and does not resend patient data to GCS. The smoke check exercises host BWA loading of the indexes; it does not establish that a full Parabricks patient analysis has passed.

Alternative `--index-mode rebuild --execute` downloads only the FASTA/FAI/dictionary/ALT objects (3.25 GB), then builds the five BWA indexes with the verified host `bwa 0.7.17 index -a bwtsw`. This avoids relying on a prebuilt binary index but consumes more CPU time. **Allow roughly 30–90 extra minutes** for whole-human BWA indexing; actual time on the selected VM has not been measured. A 7,200-second indexing timeout is enforced. Both modes require at least 25 GiB free for preparation; the later analysis runner independently enforces its larger scratch-space reserve.

Official indexing/reference guidance:

- https://github.com/lh3/bwa/releases/tag/v0.7.17
- https://packages.ubuntu.com/jammy/bwa
- https://docs.nvidia.com/clara/parabricks/tool-reference/tools/fq2bam

## Current status

Public reference preparation completed successfully using the prebuilt mode. All nine source objects passed publisher MD5 and exact byte-size checks; the complete FASTA passed per-contig sequence MD5, FAI/dictionary/order/layout checks; ALT/ANN validation and host BWA index smoke alignment passed. The final directory appeared only after these checks and contains `reference-receipt.json`. Ubuntu BWA `0.7.17-6 amd64` was installed without upgrading other packages. No patient FASTQ alignment, variant calling or private run-configuration edit was performed by this preparation. The source FASTA and reference sidecars were not modified.
