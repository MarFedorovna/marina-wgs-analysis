# Approved interpretation additions

The user approved MedInter proposal P01–P11, including GTEx/pext and UniProt, on
2026-11-06 (Asia/Qyzylorda). This records authorization and the implementation
contract; it is not a completion receipt or a statement that every resource has
already annotated patient variants.

Primary analysis continues with the reviewed local resource bundle. Newly fetched
supplements are prepared independently and must not change a running annotation
config, its frozen code, its input VCF, or existing stage receipts.

Parallel work:

- P01/P02/P08: bounded public population and mitochondrial references, ClinVar
  submission evidence and citation metadata. Keep exact release, assembly,
  allele-count denominators, coverage scope and source limitations.
- P04/P07: public heart-tissue GTEx summaries, gnomAD pext and reviewed UniProt
  features for the fixed public gene panel. Preserve transcript versions; GTEx
  gene models do not automatically match the analysis annotation release.
- P03/P05/P06/P09/P10/P11: mitochondrial evidence worksheet, NMD review method,
  independent offline HGVS validation, predictor alignment checks, IGV review
  checklist and structured phenotype/pedigree inputs.

External requests use fixed public reference genes or general source downloads.
Patient variants, genotypes and clinical records are not submitted to hosted APIs.
Reference datasets and installation caches stay on the server. Public repository
updates contain generic code, configuration templates and documentation only.

Resource receipts must distinguish downloaded source, tested runtime, prepared
local adapter and evidence actually applied to a patient result. Missing
annotations remain missing; manual ACMG/PVS1/mechanism judgments are not inferred
from a consequence label or a predictor vote.

The later instruction to work on CNV/SV in parallel authorizes a separate
feasibility assessment and preparation of a suitable research workflow. It does
not establish uniform whole-genome coverage or validate a four-sample CNV model.
The primary SNV/indel annotation does not wait for that branch. Any eventual
CNV/SV results need their own method, QC, limitations and receipts.

All jobs retain the existing server budget deadline; this approval does not
authorize a new server, a paid API subscription or extension of the budget.
