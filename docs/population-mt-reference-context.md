# Separate local population, SCV and MT reference context

`export_population_mt_reference_context.py` reads completed primary stage 06 and 07 receipts and verifies the SHA256/size of the stage-07 `master_evidence.tsv.gz`. It writes a new output tree, preserving original per-sample column values, GT/FT/DP/AD/GQ, sample state, site filters and primary ranking fields. It never rewrites the primary master, VCFs, report or ranking.

The exact public-resource ready manifest is pinned to SHA256 `465ee667cb78baddc3665d06adb15fb27ad54c47415ff14238fda5f0f42c3d85`. Each of its derived reference artifacts is checked against its own size/SHA256 before use. SQLite connections are read-only and query-only. Raw public gene-coverage responses are separately pinned and checked against the public HGNC identity/build. No networking libraries, external API calls or source-database mutations are used by the exporter.

## Output tables

The output is a new sibling of the primary workflow: `/srv/marina/outputs/population-reference-context-v1`. It contains:

- Per-sample and cohort copies of selected master rows with supplemental link/status columns. Selection covers the fixed public39 gene context and exact local public allele matches; unrelated rows are not expanded into another genome-wide report.
- A cohort long population table with separate gnomAD v4.1.1 exome/genome records, FILTER, AC/AN/AF/nhomalt and the complete retained ancestry/quality INFO. Missing source records remain explicitly missing, with blank values rather than AF=0.
- A separate exome all-sites AN v4.1 position-context table. It is not a genome AN or depth-coverage substitute.
- A public gene-coverage long table for matched HGNC genes. It is gene-bin context, with v4-exome/v3-genome provenance; it is never joined to a patient's nearest locus or used as proof of locus callability.
- A validated ClinVar AlleleID↔VariationID map, SCV/version and citation long tables, and per-sample SCV links. `clinvar_ids` is never reinterpreted as AlleleID. An explicitly supplied AlleleID must agree with the exact GRCh38 chromosome/POS/REF/ALT map; conflicts produce no SCV join. If the upstream AlleleID is absent, only the exact source-allele map can supply the separately labelled reference link. Assertions are not PP5, independent cases or an automatic classification.
- Per-sample MT reference tables and a cohort MT long table from the four reviewed final `mitochondria.filtered.vcf.gz` files. Each original CHROM/POS/REF/ALT/QUAL/FILTER/INFO/FORMAT/sample string and GT/AD/AF/DP/FT is preserved. Multiallelic records are not split or flattened. Helix's absent site AN remains absent. Reviewed NUMT, contamination, research-report and haplocheck source QC files are copied as raw provenance, without interpreting one flag as proven contamination or claiming clinical LoD.

The final `workflow-receipt.json` has a PASS identity and artifact map compatible with completed-output discovery. Earlier partial output directories are retained on failure; final promotion is exclusive. The final output cap is 2 GiB, the controller is limited to eight hours, and the original absolute VM deadline remains `2026-11-06T03:46:48.046025Z`.

The existing native-download plan helper currently excludes `.tsv.gz` suffixes. Root must add these new completed table outputs through a separately reviewed transfer allowlist update before claiming they are included in the local backup. This delta does not modify that existing helper or its active refresher.

## Review and execution

`continue_population_mt_reference_context.py prepare` creates a private rule candidate. It pins the unchanged primary 15-file runtime/config, accepted stage-05 receipt, public resources, coverage responses and all four actual MT-v2 receipts/final VCFs/QC files. It reads no patient variant data for an external request and starts no controller.

The candidate's `PREPARED_FOR_ROOT_REVIEW_NOT_RUNNING` status is rejected for execution. After root review, an exact reviewed rule must carry `ROOT_REVIEWED_EXACT_LOCAL_ONLY_POPULATION_RULE` and be passed with its SHA256:

```sh
python3 continue_population_mt_reference_context.py run \
  --rule PRIVATE_ROOT_REVIEWED_RULE.json \
  --reviewed-rule-sha256 REVIEWED_SHA256 --execute
```

The controller waits only for primary 06/07 PASS; it cannot start, restart or change the primary workflow. It rechecks reviewed pins at point of use, then invokes the local exporter. State and rule files stay private, outside the completed output tree. The rule candidate contains actual sample paths/pins and must never be published with the generic source code.

Synthetic controls cover exact allele matching and ID-type separation, unchanged GT/FT/sample quality/rank values, missing population AF versus an actual source AF=0, unsplit MT multiallelic calls, original MT filters and FORMAT arrays, master tamper rejection, both primary PASS gates, root-review denial and the unextendable absolute deadline. They do not represent an executed export on real cases.
