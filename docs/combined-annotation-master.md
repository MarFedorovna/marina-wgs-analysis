# One combined master for cohort and sample views

`join_combined_annotation_master.py` creates one new enriched master after the primary annotation/report and both actually applied reference enrichments complete. It preserves every original stage-07 master row, field value, field order and row order. GT, FT, sample quality, source classification, candidate eligibility and ranking remain unchanged. It appends 33 explicitly named `allref_*` fields containing verified tissue context, actual population/AN/SCV summaries, source versions, review flags and compact auxiliary evidence keys.

The source is the completed `07_report/master_evidence.tsv.gz`, verified against its PASS receipt. The tissue producer uses stage-06 transcript rows, so the join does **not** rely on a shared row number. It hashes the full original-field vector in exact ordered TSV form and checks row multiplicity. Every public-gene tissue row must correspond to a stage-07 original row; original-field differences or orphan supplemental rows stop the export. Duplicate identical original rows retain their multiplicity. Population sample rows receive the same full-row check and must agree with the exact chromosome/POS/REF/ALT key.

Population summaries come from the completed applied long tables, not a resource registry. Exome/genome values remain separate and name gnomAD v4.1.1. Missing AF is missing, never zero. Exome all-sites AN v4.1 is checked against its separate applied position-context table; it is not genome AN or depth. Gene coverage is summarized only as gene-bin context with v4-exome/v3-genome provenance, never a nearest-locus depth claim. SCV references are checked against the explicit exact-allele AlleleID↔VariationID map, and only the row's validated Variation IDs contribute summaries. Numerical ID equality alone is never a join rule.

GTEx heart values, transcript-version match statuses, pext context and UniProt gene-protein context come from the matching actually enriched tissue row. Unproven allele/exon/domain assignments remain unassigned. New review flags identify transcript-version context differences, population source filters or differing SCV assertion texts. They do not rerank, reclassify, assign PP5, resolve condition-specific disagreements or establish independent observations.

The master is `/srv/marina/outputs/combined-annotation-master-v1/master_evidence.enriched.tsv.gz`. Applied nuclear long tables and the separate MT appendix are copied byte-for-byte under `applied-evidence/population/`, with upstream receipt/size/SHA bindings. No public database caches are copied. MT rows are **not** added to the nuclear four-sample matrix. `auxiliary-index.json`, `applied-evidence/MT-appendix-index.json` and `index.html` provide one entrypoint and evidence keys. The downstream cohort/sample views adapter reads this exact enriched master so its master, top20 and type tables retain all `allref_*` columns.

## Manifest contract

`manifest.json` uses `schema_version=1`, `status=PASS_COMBINED_ANNOTATION_MASTER` and contains:

- `base_master` and `enriched_master`: exact path, bytes and SHA256.
- `original_fields`, `appended_fields`, `row_count` and `invariant_digest`.
- `original_rows_fields_order_preserved=true` and `ranking_eligibility_unchanged=true` after checks, not as an assumption.
- Exact `stage_receipts` for 06/07 and two `upstream_receipts`: tissue `manifest.json` and population `workflow-receipt.json`.
- `producer_code_sha256`, applied row counts, zero orphan-context counts, auxiliary artifacts with `source_upstream_receipt`/`source_relative_path`, and `MT_rows_added_to_nuclear_master=0`.

The invariant digest is SHA256 over the original ordered header and field vectors reconstructed using Python `csv.writer(delimiter='\t', lineterminator='\n')`, default TSV quoting and UTF-8. Embedded TAB/newline values are preserved through standard quoting. The output's original-field projection must have the same digest. The downstream views adapter independently compares each original row/vector and this digest. A PASS `workflow-receipt.json` with a 64-hex identity pins the enriched master, manifest and applied auxiliary results.

## Continuation and bounds

`continue_combined_annotation_master.py prepare` writes a private candidate and starts nothing. It pins the primary 15-file runtime/config, accepted stage-05 receipt, producer code, public panel, tissue context index and exact public-resource ready manifest. Execution requires a separately ROOT-reviewed rule and its exact SHA256. It waits for 06/07 PASS and both successful enrichment controller states; a failed/stopped dependency stops this continuation without retrying or changing primary annotations or views.

The helper runs with at most two process-affinity CPUs and 4 GiB address space, a disk-backed SQLite index with 64 MiB cache, no mmap, an 8 GiB private-index limit, a 6 GiB completed-output limit, and a 150 GiB VM analysis reserve. It obeys the unchanged `2026-11-06T03:46:48.046025Z` absolute deadline and an eight-hour maximum. Hidden staging is exclusively promoted only after original-row, source-pin and context-consumption invariants pass. Neither producer rules nor the primary runtime/views are modified.

Tests use synthetic data: reordered stage6/7 rows, duplicate rows, embedded TAB/LF, missing versus real AF, unchanged GT/FT/rank, SCV ID-type separation, actual auxiliary/MT separation, four dependency gates, root-review denial, and strict CPU/memory/deadline limits. They are not a claim that a real combined-case export has already run.
