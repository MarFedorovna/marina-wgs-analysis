# Direct portable reference context

The new direct-packet producer calls the preserved scientific annotation loop and appends context in its output sink. This adapter does not rerun VEP, change original fields, assign ACMG criteria, change eligibility/ranking, or claim completion of the legacy 06/07 stages. Baseline field-vector and quality-counter proofs belong to the separately reviewed direct producer.

`direct_public_context_adapter.py` builds a server-only indexed derivative of the approved public39 source release `465ee667cb78baddc3665d06adb15fb27ad54c47415ff14238fda5f0f42c3d85`. The actual index contains 220733 population rows, 243349 exome-AN positions, 109389 exact nuclear allele-map rows, 248946 versioned SCV rows, 207404 citations and 38 nuclear gene-coverage profiles. It is 1210716160 bytes, SHA256 `93cdbfd3d077952f10c21164fc3b6c87c4592ccc41122f28b994a21ce4ccf885`. Preparation is `READY_PUBLIC_DIRECT_CONTEXT_INDEX_NOT_APPLIED`, not an applied patient result.

All data comes from existing verified local public references. The adapter opens the derivative SQLite read-only with immutable mode and a 16 MiB page cache. It caches the current allele across its transcript/sample rows and compact coverage summaries by HGNC; it does not repeatedly inflate full gene-bin arrays or scan the unindexed original submissions table per row. All reference databases remain on the server.

The API is:

```python
context = Context(index_manifest, spliceai_manifest=None, verify_index=True,
                  expected_source_vcf=None)
extra_strings = context.augment(original_row_dict, one_based_global_ordinal)
evidence = context.allele_evidence(chrom, pos, ref, alt)
bundles = context.header_bundles()
context.close()
```

`APPENDED` is the exact ordered list of 33 `allref_*` fields preserved from the prior combined-master design. When an actual completed SpliceAI manifest is configured, `SPLICE_FIELDS` adds the exact thirteen fields documented in `spliceai-reference-context.md`. Every returned value is a string. The input row dictionary is not changed. Direct statuses explicitly describe direct HGNC/allele application; they do not assert comparison with a nonexistent legacy master.

Tissue context uses the preserved `supplement` function with exact HGNC and symbol checks and explicit GTEx/pext transcript-version match statuses. Allele-specific pext/domain mappings remain unverified. Population joins require exact GRCh38 CHROM/POS/REF/ALT. Full ancestry/quality INFO is retained separately from readable global counts. Exome AN is position context for v4.1 exomes only; missing is distinct from an observed AN=0. Gene coverage remains HGNC-labelled gene-bin context, never exact locus depth.

ClinVar SCV links require the source's explicit exact-allele AlleleID/VariationID map. A reported core AlleleID is never treated as a VariationID. Coordinate/snapshot conflicts produce a review flag and no scalar row assignment. Full source SCV assertions and citations remain available with the validated map in the allele packet for manual review. They do not establish independent cases, article review, PP5, or an automatic clinical class.

`allele_evidence` returns complete applied public `population`, `exome_AN`, `allele_map`, `SCV` and `citations` rows. `header_bundles` returns `public_reference` and `gene_context`: full compact GTEx/pext/UniProt source objects with hashes and full public gene-coverage objects. Both are embedded using the bounded chunk codec, without needing external database files for decoding the final VCF. The total inflated supplemental header cap remains 64 MiB. The gene bundle's `genes` values are native public gene-coverage objects; their positions are retained as returned by the source.

SpliceAI configuration requires workflow receipt PASS/stage `spliceai-reference-context-v1`, matching manifest/prediction artifact hashes, five model hashes, D500/mask0, GENCODE24 gene scope and per-row matching model/reference provenance. A producer supplies `expected_source_vcf` from its accepted original input pin to reject scores from another input. Missing/outside-panel rows have empty scores/source JSON and explicit status. A complete true zero is retained. Exact gene results populate scalar row fields; different/ambiguous gene context remains in the full prediction JSON and is not assigned to VEP transcripts.

`spliceai_vcf_codec.py` extends the immutable base codec SHA `18e7a2102aceb5b962892242763d7d331523d7669c93f16b12f9fe71a1eb238e`. It exports the same packet/schema/chunk APIs, plus `headers`, `rebuild_readable`, `select_sample`, `strip_overlay`, `INFO_IDS` and `NAMESPACED_INFO_IDS`. New IDs are standard `SpliceAI` and `WGS_SPLICEAI_STATUS`, `WGS_SPLICEAI_SOURCE`. The standard tag retains `ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`. Existing typed gnomAD fields and original CSQ/INFO/FORMAT remain intact. Per-sample packets/readable tags are rebuilt from only that sample's complete rows; global headers contain public reference objects, not other samples' genotypes.

New direct authority must name its actual core manifest, original-prefix batch proofs and extended-field batch proofs/Merkle contract. Merkle proofs must not be relabelled as the older single full-TSV SHA. Final VCF completion requires verified baseline prefix, packet decoding, original VCF record preservation, indexed access and pinned applied-source provenance. Mitochondrial genotypes/QC remain in the separate MT appendix and VCF.
