# Final VCF: core33 plus self-contained SpliceAI INFO

The reviewed direct producer may run the heavy original annotation loop and its 33 public-context fields while the independent local SpliceAI calculation is still running. Its output is explicitly intermediate: `/srv/marina/outputs/direct-core-context-v1`, status `INTERMEDIATE_CORE33_PASS`, `final_delivery_complete=false`.

`apply_spliceai_allele_overlay.py` makes the final `/srv/marina/outputs/fully-annotated-vcf-v1/cohort.fully-annotated.vcf.gz` and Tabix index after both producers PASS. It adds real allele-level SpliceAI evidence without decoding or re-encoding the 48.5 million core row vectors. Every original VCF record, INFO/CSQ/FORMAT/genotype string and `WGS_EVIDENCE` token remains byte-identical after removing only the three added INFO entries. The original physical schema, original-prefix proof batches and core33 Merkle root remain unchanged. The final manifest is `PASS_FULLY_ANNOTATED_VCF_DIRECT`; its workflow receipt is PASS/stage `fully-annotated-vcf-v1`. No legacy stage06/07 receipt is created.

## Physical INFO and logical table projection

The three added INFO IDs are `SpliceAI`, `WGS_SPLICEAI_STATUS` and `WGS_SPLICEAI_EVIDENCE`. The standard SpliceAI entries are:

```text
ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
```

DS values are actual official-model delta scores. DP values are signed genomic offsets. Distance is 500 bp and mask is 0. The gene model is the official bundled GRCh38 GENCODE24 gene union; it is not reassigned to VEP GENCODE49 transcripts. The approved public panel has 39 genes, 38 applicable nuclear gene models, and 2326 selected source alleles. Outside-panel variants and unsupported alleles do not receive invented zeros. A calculated `0.00` remains a real score.

`WGS_SPLICEAI_EVIDENCE` uses UTF-8 JSON, zlib, and unpadded URL-safe Base64, with a 64 MiB inflated payload bound. Its complete object is:

```json
{"v":1,"allele":["chr11",47332146,"T","G"],"status":"SCORED",
 "scope":"GENCODE24_GENE_UNION_NOT_ENST49_MAPPING",
 "predictions":["complete exact-allele, per-gene prediction dictionaries"],
 "source_json":"complete model/version/FASTA/annotation-hash/distance/mask JSON"}
```

Each prediction contains its exact allele, gene, status, DS/DP strings, standard raw annotation, full source JSON and model scope. Full five-model hashes, annotation hash and reference hash are inside the VCF. Whole databases or model-weight files are not embedded. If an allele was selected but unsupported, full evidence retains the reason and scores remain missing. Outside scope is explicit `WGS_SPLICEAI_STATUS=OUTSIDE_APPROVED_PUBLIC39_MODEL_GENE_SPANS_NOT_ZERO`, `SpliceAI=.` and `WGS_SPLICEAI_EVIDENCE=.`.

`WGS_SPLICEAI_SCHEMA` contains the encoded projection metadata. The exact contract is `ALLELE_INFO_V1_ROW_PROJECTION_ON_READ`. The final manifest's `schema` is the unchanged physical original-fields-plus-33 schema. Its separate `projection` object matches the decoded header and declares the thirteen logical SpliceAI fields. Consumers verify physical batch proofs first, then derive those fields from the same VCF and record separate derived-vector proofs. A Merkle root is never relabelled as a single full-TSV digest.

The pinned codec `spliceai_allele_info_codec.py` exposes:

```python
predictions, source_json = projection_from_info(info_dict, allele)
extra = enrich_row(original_row, predictions, source_json)
```

This requires no external database or prediction TSV. Scalar row scores require an unambiguous exact model-gene match. The full applicable model predictions remain available when gene/transcript mapping is unverified. Main ranking, eligibility and ACMG/classification fields remain unchanged. The separate SpliceAI research list can use actual maximum DS among quality-passing ALT carriers without an automatic pathogenicity threshold.

Sample packets retain only that sample's complete logical rows. The closed projection drops cohort `membership` and unknown private metadata. Public allele evidence and public reference headers may remain. Readable INFO is rebuilt from the selected logical rows, including full SpliceAI evidence. MT genotypes and QC remain separate and are not mixed into nuclear rows.

## Controller and verification

Preparation only:

```sh
PYTHON scripts/continue_spliceai_allele_overlay.py prepare \
  --base-rule BASE33_ROOT_REVIEWED_RULE --spliceai-rule SPLICEAI_ROOT_REVIEWED_RULE \
  --python /srv/marina/private/cadd-public-panel-v1/venv/bin/python \
  --rule NEW_CANDIDATE_RULE --output /srv/marina/outputs/fully-annotated-vcf-v1 \
  --state-directory /srv/marina/private/spliceai-allele-overlay-state-v1
```

After ROOT review, the exact copied rule must have review status `ROOT_REVIEWED_EXACT_SPLICEAI_ALLELE_OVERLAY_RULE`:

```sh
PYTHON scripts/continue_spliceai_allele_overlay.py run \
  --rule ROOT-reviewed-rule.json --reviewed-rule-sha256 EXACT_REVIEWED_SHA --execute
```

The controller waits for actual core33 and SpliceAI PASS receipts, locks its private state and never starts/stops/retries those producers. It pins their reviewed code, accepted original stage05 VCF identity, exact model/annotation/public-control provenance, and the existing pysam 0.24.1 runtime. Limits are 2 CPUs, 4 GiB address space, 48 GiB output, 150 GiB remaining disk, at most eight hours and the original absolute deadline `2026-11-06T03:46:48.046025Z`. A hard controller alarm and periodic checksum/record/disk checks preserve the deadline. New output promotion is atomic and refuses replacement.

An independent final BGZF reread verifies the original data-record digest after stripping only the added tags, exact old core token digest, original header sequence, every actual model prediction and every source allele. Random indexed access is tested on each represented contig. The public control proof covers missing/outside/unsupported states, real official-model DS/DP, phased/missing/filtered synthetic genotypes, unchanged physical tokens/schema and rejection of a different accepted input. It uses no patient genotypes or external patient API queries.

The complete master and all cohort/sample tables must be exported by the reviewed fast views controller with this projection contract. The older `embed_vcf_evidence.py decode` reconstructs only the physical core33 vectors and must not be described as the complete SpliceAI-enriched master decoder. Use the exact reviewed views rule, for example `continue_packet_views_fast.py --rule RULE --rule-sha256 SHA --execute --watch` (or the subsequent explicitly reviewed version of that controller).
