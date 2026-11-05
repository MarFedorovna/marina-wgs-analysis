# Remaining full-reference alignment

Run from the reviewed checkout on the GPU VM. After the first sample's pilot alignment has completed and the full input transfer is verified, prepare a separate queue:

```bash
python3 scripts/configure_remaining_align.py --execute
```

This builder only writes private configuration. It reads the canonical eight-column `/srv/marina/private/samples.tsv`, the pilot configuration and completed alignment receipt, and the complete transfer verification marker. It verifies that the first sample's metadata matches the canonical manifest, its recorded alignment outputs retain their SHA256 values, and every canonical FASTQ matches the expected-source upload manifest's size and SHA256. The manifest's own SHA256 must match the verified transfer plan identity.

The builder excludes the first sample and creates `/srv/marina/private/remaining-align-samples.tsv` and `/srv/marina/private/remaining-align.json`. Existing identical files can be reused; changed files are preserved and cause an error. New private files use mode 0600 and exclusive creation. Omitting `--execute` validates and prints a plan without writing configuration; this still reads and hashes the inputs.

With the private queue prepared, start alignment:

```bash
python3 scripts/run_parabricks.py --config /srv/marina/private/remaining-align.json --stage align --execute
```

The queue runs the remaining canonical samples sequentially on the single GPU. In this four-sample workflow, that means three new alignments. It uses the same complete reference, immutable Parabricks image and processing settings as the pilot, with BQSR disabled. The first completed BAM remains in its original verified pilot directory and is not copied or realigned by this queue.

Each new alignment is retained under `/srv/marina/outputs/remaining-align/<sample_id>/align/`, with marked BAM, index, duplicate metrics, samtools flagstat/idxstats/stats, command logs and a hash-bound receipt. These are private outputs, not GitHub artifacts. The alignment stage does not generate new gVCFs or variants and does not set a sex-chromosome, MT, CNV or clinical interpretation policy.

Repeating the runner command reuses a completed alignment only when its fingerprint and output checksums still agree. Changed or incomplete final outputs cause a stop. Keep the same output root so the process lock protects this queue; do not start a parallel GPU queue under another root. Inspect failed partial directories before retrying.

The generated configuration gives each runner invocation a four-hour timeout. This stops that invocation's workloads at its deadline; cloud billing continues afterward. The timeout is not a promised completion time or an account spending cap. Retain verified completed outputs before terminating the VM.
