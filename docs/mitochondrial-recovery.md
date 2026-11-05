The research MT caller uses the pinned GATK 4.3 normal/8000-base-shifted design. Original autosomal BAMs, autosomal calls and earlier MT checkpoints remain unchanged. `scripts/recover_mitochondria.py` writes a separate output namespace and validates source BAM/receipt identity before reusing completed MT selection/alignment checkpoints.

An actual run exposed `IllegalStateException: Padded span must contain active span` in `FilterAlignmentArtifacts.makeAssemblyRegionFromVariantReads`. The [pinned upstream source](https://github.com/broadinstitute/gatk/blob/4.3.0.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/realignmentfilter/FilterAlignmentArtifacts.java) constructs an assembly window ending at the supporting read end plus 50 without clipping to contig length. The same expression remains in the 4.6.2.0 source. This tool also compares the numeric BWA-image reference ID with the `-R` dictionary index; MT-only and whole-reference dictionaries must not be mixed for that comparison.

Recovery converts the normal MT BAM through textual SAM using the complete, verified whole-reference sequence dictionary. It checks that every alignment record is unchanged, then uses the matching whole-reference FASTA and BWA image for the additional artifact assessment. This changes numeric BAM reference IDs only; coordinates, reads, read groups and genotypes are retained.

Only candidates entirely within normal-reference positions 576–16024 enter this artifact assessment. Shifted/lifted circular-boundary candidates remain in the complete filtered VCF with `NUMTArtifactCheckUnavailable` and `NUMT_ASSESSMENT=UNAVAILABLE_CIRCULAR_BOUNDARY`; they do not enter the separate PASS subset. All existing filters are retained, and every raw allele row and FORMAT genotype/DP/AD/AF is checked for preservation. An assessment failure stops the recovery; it is not converted into an assessed result.

The recovery explicitly passes `--dont-skip-filtered-variants true`. Upstream FAA otherwise skips already filtered candidates, which the strict row-preservation assertion detected in the first real recovery attempt. That failed attempt and its raw intermediates remain unchanged. Recovery version 2 uses a fresh output namespace.

Synthetic recovery validation recovered the known interior and control-region alleles, kept both rows and their DP/AD/AF, assessed the interior allele without a spurious alignment artifact, and explicitly flagged the boundary allele as unavailable. Output artifact SHA256 verification passed. The unit fixture additionally checks position 16560, prior-filter retention and unchanged genotype fields. The end-boundary unit test is a merge-policy test, not evidence of successful GATK boundary realignment.

`scripts/validate_mito_recovery.py` additionally runs actual FAA against the synthetic interior allele after assigning an existing `weak_evidence` filter. The filtered row was assessed, its original filter and genotype fields were retained, and the interior allele had no spurious alignment artifact. Control and fabricated end-boundary rows remained explicitly unavailable. This regression uses only the earlier generated synthetic fixture; the fabricated end-boundary row tests handling policy, not variant-caller sensitivity.

Example, using private user-supplied inputs:

```bash
python3 scripts/recover_mitochondria.py \
  --workflow-config /srv/project/private/workflow.json \
  --resources /srv/project/resources/mitochondria/resource-manifest.json \
  --source-root /srv/project/outputs/mitochondria-original \
  --output-root /srv/project/outputs/mitochondria-recovery-v1 \
  --execute
```

Without `--execute` the command validates inputs and emits a plan. The source and destination must differ. Execution is sequential with four-CPU/16-GiB container limits and the existing global VM-budget deadline. The source whole-genome BAMs are read only. Re-execution validates completed recovery checkpoint identity and hashes rather than overwriting them.

After every sample checkpoint completes, perform the independent technical audit:

```bash
python3 scripts/audit_mitochondria.py \
  --output-root /srv/project/outputs/mitochondria-recovery-v1 \
  --samples-manifest /srv/project/private/samples.tsv \
  --report /srv/project/private/mitochondria-technical-audit.json
```

The audit verifies all completed artifact sizes/SHA256, sample headers, indexed chrM retrieval, exact raw/filtered allele keys and GT/DP/AD/AF, NUMT status counts and the separate site-PASS subset. It does not modify analysis outputs. Record counts are not counts of pathogenic variants. Use a new report filename for a repeated audit. Reports and sample manifests contain private data and belong outside Git.

These are research small-variant/heteroplasmy outputs. Clinical sensitivity, low-level heteroplasmy limit of detection, mtDNA deletions, tissue specificity and complete exclusion of polymorphic NUMTs are not validated. A completed technical receipt does not establish clinical callability or pathogenicity.
