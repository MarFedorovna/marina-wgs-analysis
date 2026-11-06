# Mitochondrial functional tables

`scripts/export_mt_functional_tables.py` organizes an actual completed
[MT consequence annotation](mt-consequence-annotation.md) into cohort, shared and
individual views. It performs no calling, reannotation, public API request or
database download.

## Input authority and sharing

The exporter requires `PASS_MT_CONSEQUENCE_ANNOTATION` and the actual
`mitochondrial-consequence-annotation-v1` PASS receipt, bound to the reviewed source
profile and producer. It decodes each sample's persisted
`WGS_MT_CONSEQUENCE_EVIDENCE` packets, checks every field vector against that sample's
`MT-functional-all-transcripts.tsv.gz`, and verifies ordered-vector digests. The
complete original transcript rows are retained in the output masters.

Sharing uses **exact chrM/POS/REF/ALT identity**, not gene overlap. Source ALT vectors
seed the observed alleles. `ALLELE_NUM` projects transcript consequences onto the
corresponding ALT; unresolved multiallelic effects remain in the complete master
and are not guessed onto either allele. Different transcript effects for one allele
are unioned without counting another allele or another sample.

There are two explicitly separate levels:

| Level | Membership |
|---|---|
| `observed` | The allele occurs in a source Mutect2 call record |
| `site_PASS` | At least one corresponding source record has literal `FILTER=PASS` |

Observed call-record membership is not a nuclear GT carrier classification.
Original `GT=0/0` does not erase an observed MT record with an AF/HL estimate. No
nuclear GQ threshold, diploid rule or `candidate_eligible` field is introduced.
Missing sample records are labelled `NOT_OBSERVED_NOT_ASSESSED_NO_REFERENCE_INFERENCE`;
they are never converted to reference genotypes.

Tables retain the original ALT vector and index, GT, AF vector, HL, AD vector, DP,
QUAL, FILTER, FORMAT and raw sample string, plus frequency context and the full
transcript appendix. This is descriptive evidence, without LoF/GoF/PVS1, automatic
ACMG, clinical burden or limit-of-detection conclusions.

## Output layout

The new output root is `mitochondrial-functional-views-v1`:

```text
index.html
manifest.json
workflow-receipt.json
masters/cohort.master.tsv.gz
masters/sample_0001.master.tsv.gz ... sample_0004.master.tsv.gz
tables/functional-category-index.tsv
tables/{cohort,shared_ge2,all4,sample_0001...sample_0004}/
  {observed,site_PASS}/all.functional.tsv.gz
  {observed,site_PASS}/<consequence>-<stable-short-hash>.tsv.gz
```

`shared_ge2` requires at least two members of the selected level; `all4` requires
all four. Sample labels are mapped to original sample IDs in the manifest. Every
observed VEP consequence receives a category file. Fixed stop-gain, stop-loss,
start-loss, frameshift, donor/acceptor splice, splice-region, missense, synonymous,
in-frame insertion/deletion and protein-altering categories also exist when their
observed count is zero. Zero means zero observations in that source scope, rather
than missing data or proven biological absence. Categories can overlap.

`PASS_MT_FUNCTIONAL_VIEWS` plus the actual PASS receipt with stage
`mitochondrial-functional-views-v1` identifies completion. Every declared output has
its own size/SHA256 pin. Progress is written separately to `mt-functional-status.json`
beside the reviewed execution profile.

## Commands and limits

Variables below refer to existing reviewed private inputs and a new output. Expired profiles are historical; a new run needs a newly reviewed deadline within the code limit.

```sh
python scripts/export_mt_functional_tables.py prepare \
  --helper scripts/annotate_supplementary_callsets.py \
  --codec scripts/embed_vcf_evidence.py \
  --input-profile "$REVIEWED_MT_PROFILE" \
  --output "$NEW_MT_FUNCTIONAL_OUTPUT" \
  --profile "$NEW_MT_FUNCTIONAL_PROFILE_CANDIDATE" \
  --deadline-utc "$REVIEWED_DEADLINE_UTC"
```

A separately reviewed copy sets `root_reviewed=true` and
`status=ROOT_REVIEWED_MT_FUNCTIONAL_VIEWS`. Without `--execute`, execution is plan-only:

```sh
python scripts/export_mt_functional_tables.py run \
  --helper scripts/annotate_supplementary_callsets.py \
  --profile "$REVIEWED_PROFILE" --reviewed-profile-sha256 "$REVIEWED_PROFILE_SHA256" --execute
```

The bounded exporter uses two CPUs, 4 GiB address space, a 512 MiB output cap and
150 GiB disk reserve. This run's code refuses deadlines later than 2026-11-06
06:25 UTC; its reviewed production profile uses 06:20, with an independent outer
process-group stop by 06:22. Existing outputs are not overwritten.

The public test suite includes four actual synthetic MT VCF/BGZF/TBI pipelines,
genotype-versus-heteroplasmy semantics, exact multiallelic sharing, missing-not-reference,
different transcript effects on one allele, unresolved ALT projection and explicit
zero categories. No patient export is performed by the tests.

## Suggested README entry

MT annotation is delivered in two separate stages: [mitochondrial consequence
annotation](mt-consequence-annotation.md) preserves original Mutect2 AF/HL and QC while
adding offline VEP/ClinVar and existing population context inside each sample VCF;
[MT functional views](mt-functional-tables.md) reconstruct all transcript evidence
from those VCFs and provide individual, cohort, shared-at-least-two and all-four
tables for every functional consequence. Observed records and literal site-PASS
records remain separate; missing is not reference, and no clinical class is assigned.
