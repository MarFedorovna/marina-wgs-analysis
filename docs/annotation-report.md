# Offline annotation and evidence report

The Russian [source and consequence map](sources-and-consequences.md) explains the user-supplied Analysis 1/2 lists, actual versus planned resource use, SNV/indel and transcript consequences, and the distinction between pLoF candidates and established LoF/GoF mechanisms.

Workflow stages `annotate` and `report` use `scripts/annotation_stages.py` through the shared stage runtime. Annotation reads **all** normalized, decomposed soft-filtered records from `05_filter/cohort.filtered.vcf.gz`, not just its quality subset. The official VEP 115.2 image runs with `--offline`, the matched 115 cache, the fixed analysis FASTA and Docker networking disabled. The [database installer](databases.md) must have a verified final resource manifest before patient annotation can begin.

Real sequence alleles are annotated across transcripts. Spanning-deletion `*`, symbolic/breakend/nonsequence placeholders are retained in the annotated VCF and master table without being counted as SNVs, independent sequence alleles or small-variant candidates. A one-to-one ordered merge preserves the input sample GT/FT. A skipped real allele or changed genotype causes a stop. Missing/failed FT and partial GT remain unknown/filtered; they never become confidently reference samples.

Outputs include an indexed annotated VCF; `transcript_samples.tsv.gz`; variant and gene intersections; raw/QC/candidate counts; and exact resource availability. All transcript/sample rows remain in the master table. Per-sample supplementary tables retain observed ALT calls, including filtered calls. The complete evidence includes GT/DP/AD/GQ/QUAL, AD-derived VAF, transcript/HGVS/exon/protein fields, historical cache ClinVar separately from the newer custom ClinVar snapshot, review/conflicts/condition identifiers, population values and missingness, gene-disease records and local adapter assertions.

## Context and research prioritization

Set a structured `annotation.disease_context` to rank ClinGen evidence for a particular controlled disease/inheritance context, for example `{"mondo_ids":["MONDO:REPLACE_WITH_REVIEWED_ID"],"inheritance":["AD"]}`. A free-text cohort description is recorded but does not make a Definitive association for an unrelated disease into HCM evidence. Optional `target_genes` membership is a separate ranking field; provide `target_gene_list_provenance` with source/version/date. No old 27/57/200-gene list is inferred. Optional local HPO terms may be common to the cohort or provided in `hpo_terms_by_sample`; matches are exact-term descriptions, not automated diagnoses.

The main research queue requires a quality-passing ALT plus either an external ClinVar P/LP assertion for review, or a functional consequence with observed AF below the configured research threshold or missing population data requiring review. Missing AF is never labelled zero or rare. Defaults are maximum observed AF 0.01 and research CADD PHRED threshold 20. Both are research prioritization settings, not ACMG criteria or disease-specific validated frequency cutoffs.

[AlphaMissense canonical and isoform scores](alphamissense.md) are additionally mandatory by the user's explicit requirement. They have separate model-score/label/source/version/matching fields and a separate research-ranking component. Exact transcript-version and protein-change compatibility is required; the source's GENCODE v32 prediction is not silently applied to every newer transcript. A model label is not an ACMG class, and AlphaMissense does not become the uncertain fourth tool in the original vote.

Top10/top20 per sample and cohort are deterministic, without padding. Their columns show context, gene validity, external review/conflicts, consequence, population and predictor components separately. Low-quality observed ALT with a ClinVar P/LP assertion goes to a separate urgent manual-review table. It is not ranked as a validated candidate. Full tables retain the other records and transcript effects.

Protein beginning/middle/end means thirds of the stated transcript protein length. It does not define a functional domain or NMD. Stop-gain, stop-loss, start-loss, frameshift, donor/acceptor and other splice, synonymous and missense fields are separate. A pLoF candidate is not automatic LoF, PVS1, GoF or pathogenicity. The SIFT/PolyPhen/FATHMM/CADD 3-of-4 vote is only a correlated research description; missing predictions remain incomplete, and the fourth predictor identity is explicitly unconfirmed unless configured otherwise.

The counts use unique normalized alleles, not transcript rows. `consequence_counts.json` separates exclusive primary consequences, exclusive SNV/indel type histograms and overlapping consequence flags, at raw/QC/candidate levels per sample and cohort. SVG/PNG figures use pinned Matplotlib 3.10.8 through `tools.report_python`, default `/srv/marina/tools/report-venv/bin/python`. The report HTML links these figures and all tables. A local Word export is created with `scripts/export_report_docx.py`; no automatic external document upload occurs.

## Local database adapters

The [source registry](../config/evidence-sources.example.json) is an audited requirements catalogue, not proof of installation. Every resource gets an actual available/manual/unavailable/separate-branch status. Current dbNSFP registration/subscription, HGMD Professional and OMIM access are not assumed. A public legacy FATHMM backend and separately licensed CADD bulk data require their own prepared local backend or exports; downloading a source file alone does not establish prediction coverage.

Authorized local exports can be attached through `annotation.local_adapters` in a new reviewed private workflow configuration, without editing the completed public resource bundle. Each entry needs `name`, `path` (absolute), `bytes`, `sha256`, `version`, `status: "available"`, `local_adapter: true`, `source_url`, `license_status`, and `format`. Actual file sizes/hashes are checked. Changed inputs require a new reviewed output root rather than overwriting completed stages. Include adapter paths in the outer workflow's input identity when enabling them.

For `format: "variant_evidence_tsv"`, use a normalized biallelic TSV with:

```text
chrom  pos  ref  alt  transcript  FATHMM  CADD_PHRED  source  source_version  source_url  assertion_date
```

These names represent TAB-separated columns. `transcript` may be `*` for an allele-level assertion; use an exact transcript/version for transcript-specific predictions. Optional `gnomAD_AF` and `gnomAD_<ancestry>_AF` fields carry separately versioned AF; optional other columns are preserved verbatim in `local_evidence_json`. Conflicting duplicate predictor values do not receive a fabricated consensus. This adapter does not silently expand a protein-only prediction into every transcript or permit unnormalized coordinate matches.

The `gene_evidence` local TSV can hold `gene`, disease, inheritance, mechanism and source/curation metadata. HPO gene adapters use `gene` or `gene_symbol` plus `hpo_id`; pathway adapters use `gene` and `pathway`. Do not treat a gene annotation or dosage score as a variant mechanism. Review the original assertion and keep disease/inheritance contexts separate.

## Manual evidence and scope

`ACMG_manual_worksheet.tsv.gz` starts with unassessed criteria, not automatically assigned evidence. An assessor must enter source, date, rationale and strength using current gene-specific ClinGen guidance. Legacy PP5/BP6 labels are not automatically applicable. [ACMG/AMP 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/) and [ClinGen PVS1 guidance](https://pubmed.ncbi.nlm.nih.gov/30192042/) are framework references, not evidence that a particular criterion applies.

Pathway/enrichment analysis was withdrawn in the primary video and is disabled by default. Explicitly requested local descriptions remain descriptive. Enrichment additionally requires an explicit assay-appropriate background containing all selected genes; missing/invalid backgrounds return an explicit status. All tested pathways and BH-adjusted values are retained without a disease conclusion.

The scope and limitations report records reference/tools/callsets, actual filtering policy, stage receipt hashes, resource versions, and missing phenotype/coverage/inheritance evidence. The small-variant workflow does not establish CNV/SV/repeat, MT heteroplasmy, sex-chromosome or mosaicism completeness. Missing records/low coverage do not establish reference or a negative clinical result. No de novo, compound heterozygous, GoF, ACMG class or disease conclusion is inferred automatically.
