# Mitochondrial phenotype context review

This local tool reviews the exact alleles in the completed mitochondrial functional export's `cohort/site_PASS/all.functional.tsv.gz`. PASS describes the source caller's site filter. It does not establish disease causality.

```sh
python scripts/mt_phenotype_context_review.py \
  --results-root "$RESULTS_ROOT" \
  --hpo-genes "$HPO_GENES_TO_PHENOTYPE_TSV" \
  --output "$NEW_MT_CONTEXT_OUTPUT"
```

The results root must contain `Reports/mitochondrial-functional-views-v1` and `Reports/mitochondrial-consequence-annotation-v1`, with their completed manifests, workflow receipts, functional table, and per-sample `master.from-vcf.tsv.gz` files. Sample names come from the source manifest. The HPO input is a local `genes_to_phenotype` TSV containing `gene_symbol`, `hpo_id`, and `hpo_name`. Only exact `HP:0001639` gene associations are used. No network calls, new patient-payload hashes, or source changes occur.

The tool joins exact mitochondrial position, REF, and matched ALT. Direct transcript consequences determine the gene context. Genes present only in upstream/downstream annotations are listed separately and cannot supply a phenotype match. A gene-level HPO association does not classify an individual variant. The tool preserves existing ClinVar assertions as source fields without adjudicating them.

All source GT, AD, AF, HL, DP, filters, and sample vectors are retained and crosschecked against the consequence masters. Exact ALT fractions use the original ALT index. For example, ALT `G,T` with AF `0.945,0.037` gives T an AF of `0.037`; it does not take the maximum vector value. Missing HL or population frequency stays missing, and an unobserved call is not a reference genotype. No diploid model is imposed on mitochondrial calls.

Existing gnomAD MT and HelixMTdb context, match statuses, raw reference JSON, and separate homoplasmic/heteroplasmic population AF fields are retained. Multiallelic no-join statuses remain unavailable. A source AF above 0.01 is reported as a descriptive study-frequency flag only, with its dataset and field name; it is not an ACMG threshold or a benign classification.

The new output directory contains `review.json` with every selected allele and full context, `alleles.tsv` with one row per exact allele, `sample-allele-fractions.tsv` with raw and exact-ALT sample fractions, and `summary.json` with counts. Existing output directories are refused. No clinical rank, ACMG criteria, LoF/GoF claim, or causal HCM label is assigned. The mitochondrial review remains separate from nuclear variant ranking.

Run the synthetic controls with:

```sh
python -m unittest discover -s tests -p test_mt_phenotype_context_review.py
```
