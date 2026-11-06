# Independent mitochondrial result backup

The completed `mitochondrial-reference-context-v1` output is an optional,
independent result producer. It uses the four actual MT-v2 call receipts and the
previously reviewed local MT reference release. It has no stage07 or nuclear
annotation dependency, and does not delay the cohort/sample VCF delivery.

The strict native planner checks its exact producer pair, reference metadata,
technical audit, four source receipts and source VCF/QC/index identities. It
requires every declared original QC/VCF/index copy and all four reference tables.
Only completed result tables, original filtered MT VCF/index copies, QC, HTML and
provenance enter the plan. The reference SQLite database is never copied.

The revision2 small-result worker can import this completed producer while its
only nuclear completion prerequisite remains `cohort-sample-views-v1`. Missing
or unsupported population frequency is not interpreted as zero. This result
adds research context; it does not assign a clinical class or detection limit.

Local capacity is evaluated against all selected results and resumable partials
while preserving the 15GiB reserve. Insufficient space blocks the new plan; it
does not silently omit requested tables or delete existing files.
