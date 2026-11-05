# Standalone consequence tables v1

`scripts/export_consequence_tables.py` reads the completed stage-6
`transcript_samples.tsv.gz`. It does not replace stage 6/7, existing top-10/20
queues, the HTML report or the sealed Word report. This is a separate research
export, not an automatic clinical interpretation. No remote services are used.

The export contains `index.html`, `index.tsv`, `manifest.json`, and a separate
gzip TSV for each sample and for the cohort in every consequence category.
Core categories remain present with a header and count zero when absent:
stop gained/lost, start lost, frameshift, splice donor/acceptor/region, missense,
synonymous, inframe insertion/deletion. Every additional observed SO term is
exported too, including intronic, UTR, intergenic and unknown future terms.
Empty consequence values have a distinct `missing_consequence` group.

`pLoF_candidates` overlaps with the consequence categories. It includes stop
gained, frameshift, donor/acceptor consequences or an explicit source candidate
flag. This label does not establish LoF, PVS1 or pathogenicity. Stop loss does
not establish GoF. Mechanism review defaults to `UNASSESSED`; existing explicit
`mechanism_review_status` values and raw `mechanism_status` are copied exactly.
No additional clinical evidence is invented.

Counts are unique source `variant_id` alleles per scope/category. Multiple
transcripts and overlapping categories must not be summed as independent
variants. Each source variant ID is validated against its genomic coordinates.
Sample tables require at least one numeric ALT genotype token greater than zero.
`1/.` is an observed partial ALT call and remains marked unknown and partial;
`./.` does not enter sample tables. Cohort tables retain all source alleles,
with observed ALT, reference, unknown and invalid genotype sample sets kept
separate. Absent source records are not inferred to be reference samples.
Site FILTER and FT are copied separately with PASS/FILTERED/UNKNOWN quality
levels; candidate eligibility is copied separately. Low quality is not benign.

All original transcript rows and columns, including quality, ClinVar, predictors
and provenance, are copied exactly into `evidence/transcript_effects.tsv.gz`.
Three export identity columns are prepended. Each category row points to that
file and gives `export_allele_id` plus an optional sample ID; match those
columns to recover *all* underlying transcript effects. Aggregate quality and
source statuses cover all source transcript rows for that allele/scope, not
only the category-matching transcript. JSON cells preserve genotype strings.
Import gzip TSV as text (including GT); Excel's automatic date/formula import
must not be used for exact raw evidence. Free display labels in the TSV index
are formula escaped; HTML labels are escaped. Raw evidence semantics are exact.

SQLite uses a 32 MiB page cache and disk temp storage. Alleles, memberships and
deduplication live on disk. Python retains one allele's distinct call states
(normally one per sample), never all transcript rows. Outputs are ordered by
stable first source occurrence IDs, with categories in lexical order. Existing
output directories are rejected. A PASS manifest appears only at completion;
failed directories are marked INCOMPLETE and should not be used. SHA256 checks
before and after extraction detect source changes. Plan disk capacity for the
SQLite working database, complete raw evidence, and overlapping category tables.

Run meaningful synthetic tests without patient data:

```sh
python3 -m unittest discover -s tests -p test_export_consequence_tables.py -v
```

To export a reviewed source directly (new output only):

```sh
python3 scripts/export_consequence_tables.py \
  --source /private/completed/transcript_samples.tsv.gz \
  --receipt /private/completed/workflow-receipt.json \
  --output /private/separate-consequence-export-v1
```

The optional wrapper requires BOTH stage 6 and stage 7 PASS, checks the source
against the completed stage-6 SHA256, records both receipts, and writes outside
the workflow tree. It never restarts, edits, or signals running annotation.
By default it exits NOT_READY immediately. A bounded continuation can be
started after review with `--wait-seconds` (0–43200). No watcher is installed
or started by this package. Example prepared command for the server:

```sh
bash /srv/marina/private/reviewed-consequence-export-v1/scripts/run_consequence_export.sh \
  --workflow-root /srv/marina/outputs/workflow-v1 \
  --output /srv/marina/private/consequence-export-v1 \
  --wait-seconds 0
```

Patient outputs and source identities stay private. Publish code and synthetic
test descriptions only. Fixture counts are synthetic and cannot be reported as
patient results. Report integration is a separately reviewed sidecar link or
appendix; it must not overwrite sealed reports.
