# Local indexed SNV read-support review

`review_candidate_snv_reads.py` adds a separate research read-support report to an existing candidate table. It reads indexed local BAM intervals only. It does not call variants, rerun annotation, modify existing results, infer clinical classifications, or send information to an external service.

## Isolated dependency

Use Python on Windows and a separate environment. Bamnostic 1.3 is a pure Python BAM reader. The explicit BAI path is required by this workflow; index discovery is not relied upon.

```powershell
python -m venv .venv-read-review
& '.venv-read-review/Scripts/python.exe' -m pip install 'bamnostic==1.3'
& '.venv-read-review/Scripts/python.exe' scripts/review_candidate_snv_reads.py --help
```

The public [Bamnostic API](https://bamnostic.readthedocs.io/en/stable/bamnostic.html) documents `AlignmentFile(path, 'rb', index_filename=..., require_index=True)`, half-open indexed fetch, and BAM sequence orientation. The [1.3 implementation](https://bamnostic.readthedocs.io/en/stable/_modules/bamnostic/core.html) exposes `cigartuples` as `(operation, length)` and `query_qualities` as numerical Phred scores without an ASCII offset. The script treats 255 as unavailable quality.

## Inputs and command

The candidate TSV, optionally gzip compressed, must contain `chrom`, 1-based `pos`, `ref`, `alt`, `hcm_phenotype_priority`, and `sample_calls_json`. Only rows whose HCM priority is exactly `1` are selected, up to 500 unique exact alleles. No top20 restriction is applied. The original sample-call JSON, including GT/FT/DP/AD/GQ and source flags where present, is retained in the new report.

The private BAM manifest explicitly maps the same four sample identifiers as the candidate source calls:

```json
{
  "schema_version": 1,
  "samples": {
    "sample_a": {"bam": "local/sample_a.bam", "index": "local/sample_a.bam.bai"},
    "sample_b": {"bam": "local/sample_b.bam", "index": "local/sample_b.bam.bai"},
    "sample_c": {"bam": "local/sample_c.bam", "index": "local/sample_c.bam.bai"},
    "sample_d": {"bam": "local/sample_d.bam", "index": "local/sample_d.bam.bai"}
  }
}
```

Keep this manifest and generated reports private. BAM `@RG` SM values, when available, must exactly match the explicit sample identifier. A header without SM remains explicitly marked as unavailable, rather than being interpreted as identity confirmation. This check does not establish biological identity. Contig names and bounds must match the BAM dictionary exactly; there is no automatic `chr` aliasing.

```powershell
$reviewPython = '.venv-read-review/Scripts/python.exe'
& $reviewPython scripts/review_candidate_snv_reads.py `
  --candidates "$env:HCM_CANDIDATES" `
  --bam-manifest "$env:READ_REVIEW_BAM_MANIFEST" `
  --output "$env:READ_REVIEW_OUTPUT" `
  --workers 4 --mapq 20 --baseq 20
```

This first command validates the candidate list and file metadata without opening BAM bodies or creating output. Add `--execute` after reviewing the plan. Output must be a new absent directory. A failed partial report is preserved and cannot be reused as a completed report; choose a new output directory for a reviewed repeat.

Execution uses at most four processes, 500 targets, and 100,000 fetched reads per target per sample. `--max-hours` defaults to two and accepts at most four; workers check their deadline between targets and fetched reads. Reaching the read cap or deadline fails the complete report rather than silently truncating counts. Every mapped BAM and BAI is checked for unchanged size and modification time before/after use; no BAM/BAI payload hash is computed.

## Counting rules

Only single-base distinct A/C/G/T REF and ALT alleles are assessed. Indels and other alleles remain as `NOT_ASSESSED_UNSUPPORTED_NON_SNV` rows for every sample. Missing coverage never becomes a reference genotype or a benign assessment.

Each target fetch is `[POS-1, POS)`. CIGAR traversal supports M, I, D, N, S, H, P, `=`, and X; obsolete B and unknown operations are excluded with explicit reasons. A deletion or reference skip over the locus does not supply an aligned base. Stored BAM SEQ is used in reference orientation, including reverse-strand reads, without another reverse complement.

Reads are excluded for unmapped, secondary, QC-fail, duplicate, or supplementary flags; MAPQ below the threshold or 255; missing/low base quality; or an ambiguous base. Reasons are counted. Flag reasons can overlap and therefore are not an additive partition. Passing read counts split REF/ALT/other by forward/reverse alignment orientation.

Template consensus groups passing bases by `(RG, query_name)` only in memory. Concordant overlapping mates count once; differing passing bases mark an ambiguous template excluded from template allele counts. Unknown query names remain separate observations and are explicitly counted. Templates with both read orientations are reported as `both`, rather than double counted. Missing RG likewise keeps each read separate, even for identical known query names, and is counted as `passing_reads_missing_RG_uncollapsed`: exact mate deduplication is unavailable in those observations. Raw names and sequences are never written.

Descriptive flags include fewer than three ALT templates, ALT reads on only one strand, or no passing base coverage. The ALT-template flag is also emitted for reference/non-carrier rows and requires interpretation with the unchanged source eligibility/call flags. They are not automatic artifact decisions, quality thresholds for diagnosis, ACMG evidence, or negative clinical findings. These counts need not equal GATK DP/AD: this is a targeted read-level check with specified filters, no BAQ, no local reassembly and no variant caller likelihood model.

## Outputs and verification

- `sample_0001.read-support.tsv` through `sample_0004.read-support.tsv`, in private manifest order: all selected targets and explicit source sample calls.
- `cohort.read-support.tsv`: combined sample rows.
- `summary.json`: metrics, original calls, source metadata, producer code SHA, bounds and limitations.
- `targets.igv.bed`: local 0-based half-open target coordinates, generic labels; opening IGV is a separate manual step.
- `status.json`: PASS only after all sample results and unchanged-source checks succeed.

```powershell
& $reviewPython -m unittest discover -s tests -p test_candidate_snv_reads.py -v
```

Tests cover clipping/insertion/deletion/reference skip and `=/X`, reverse orientation, quality/flag filters, missing quality, template agreement/disagreement, unknown names, no reference-genotype inference, read-cap failure, header mismatch, unsupported alleles, and an actual four-process Windows CLI integration using the dependency's public example BAM/BAI. No patient production is executed by these tests.
