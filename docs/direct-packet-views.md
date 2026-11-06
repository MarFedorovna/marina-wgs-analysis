# VCF-first cohort and four sample views

The active candidate reads one completed self-contained cohort VCF. Its stored
row vectors retain the exact frozen annotation fields plus the 33 public-context
fields. Full SpliceAI predictions, model/reference hashes, scope and status live
in that same VCF's allele INFO. The reviewed pure projection adds 13 logical
SpliceAI columns when generating tables; it does not claim those columns were
present in the stored core33 row packets.

The new direct producer has its own authority. It never creates synthetic legacy
stage06/stage07 receipts or requires the historical joined CSV as report input.
Both ordered canonical-TSV batch commitments are verified: the complete stored
vectors and their exact frozen original-field prefix. Derived logical vectors
have a separately reported commitment.

## One reviewed continuation command

The rule is private, contains exact code/source/control identities, and must be
reviewed before launch. The default command reports readiness only:

```sh
python3 scripts/continue_packet_views_fast.py \
  --rule "$FAST_VIEWS_RULE" --rule-sha256 "$FAST_VIEWS_RULE_SHA256"
```

After the exact rule and source package are reviewed, the same command with
`--execute --watch` waits for `PASS_FULLY_ANNOTATED_VCF_DIRECT`, then produces
the independent report views. It controls only its own lock and worker pool.
It neither stops nor starts VEP, calling, database installation or other queues.
Polling is bounded to eight hours and the original absolute VM deadline.

## Outputs and interpretation

`cohort-sample-views-v1` contains one full cohort VCF/TBI, four individual
VCFs/TBIs, the cohort master and four sample masters, top20 tables, functional
consequence tables/categories, exact-allele memberships and a manual ACMG
worksheet. All transcript vectors and raw quality fields survive. Sample
masters/VCFs retain observed ALT calls, including filtered and partial calls,
with explicit flags. Missing calls are never converted to reference.

Individual VCF packets contain only that individual's rows. Cohort membership
objects and unrecognized private packet keys are removed. Public reference
headers and public allele/model evidence may remain. Original cohort AC/AN/AF
are labelled `COHORT_*`; subset AC/AN/AF are calculated from the selected GT.
Every emitted sample byte stream is hashed, reread from BGZF, and checked again
while its compressed blocks are assembled. TBI queries independently check
sample identity, ALT presence and derived counts.

The primary top20 uses the unchanged frozen `rank_key` and source eligibility,
one entry per exact normalized allele, without padding. Functional groups count
distinct alleles and may overlap. Shared alleles and gene overlap are different
auxiliary views. No automatic ACMG class, PVS1/PP3 criterion, NMD result, LoF/GoF
mechanism or clinical burden is inferred.

SpliceAI has a separate prediction-only scored table and cohort/sample top20 by
actual maximum delta among complete PASS ALT calls. This can include intronic
alleles excluded from the original functional ranking. No clinical threshold or
new main-ranking rule is introduced. Missing/unsupported scores remain missing;
actual zero remains zero. GENCODE24 gene-union predictions are not reassigned to
unverified GENCODE49 transcript versions.

## Bounds and backup

At most eight workers share an eight-CPU affinity. Each worker is bounded to
1.5GiB address space. A source batch is limited to 64MiB text; queued source
batches total at most 128MiB, plus one current batch. Top20 candidates are merged
incrementally; complete transcript vectors are written once into ordered gzip
shards, without a giant SQL sort or repeated CSV scans. Parent merging keeps
gzip member order and writes each final table header once.

The strict native planner recognizes the new direct authority and copies the
canonical final VCF once. Intermediate core/joined matrices, reference database
caches and gVCFs stay on the server. BAM expected hashes retain their historical
meaning; local size-only verification remains explicitly labelled.

Before publishing a plan containing final results, supply the existing local
destination and independent small-result staging paths using
`--local-results-root` and `--local-results-staging`. The metadata-only capacity
inventory subtracts existing completed files and resumable partial bytes. If all
required results cannot fit above the unchanged 15GiB reserve, the new plan is
not published. No requested tables are silently dropped, no existing files are
deleted, and the active transfer worker is not restarted.
