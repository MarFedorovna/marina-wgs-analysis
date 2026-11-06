# Independent mitochondrial population reference appendix

`export_mitochondrial_reference_context.py` applies the already prepared public MT references to the four successful original mitochondrial-v2 call sets. It wraps the preserved `mt_exports` function from exporter SHA `4fc954c6cacce92a199fa477fcb75a3c662ecb1fcd0b0dce0ef65b6e13301832`. It does not require nuclear stage06, stage07, a cohort TSV, SpliceAI or the direct nuclear producer. No database is downloaded, and no primary or MT-v2 input is changed.

Inputs are copied from the previously ROOT-reviewed MT source profile: all four actual PASS call receipts, final filtered VCFs and their receipt-bound indices, original research/NUMT/haplocheck/contamination QC, the MT-v2 technical audit and the exact public resource release. The MT SQLite database is SHA `f069164654316043d8a96f47ad5b29413d7a28b0b9ba320648206110e42739b0`. Only this 81.9 MB public reference database is opened read-only; the nuclear databases are not opened or exported.

The output is `/srv/marina/outputs/mitochondrial-reference-context-v1`:

- `cohort/MT-population-reference-long.tsv.gz` contains the separate four-sample MT long table.
- `samples/<sample>/MT-population-reference.tsv.gz` contains that sample's complete original calls and exact-allele population context.
- `samples/<sample>/MT-QC/` contains byte-preserved original QC, labelled by its original source pins.
- Original filtered MT VCFs and their indices are copied without changes into each sample directory.
- `index.html`, `manifest.json` and `workflow-receipt.json` provide the entry point, full source provenance and completion authority.

All original CHROM/POS/REF/ALT/QUAL/FILTER/INFO/FORMAT/sample strings and raw GT/AD/AF/DP/FT are compared row-by-row with the original VCF, in original order. Raw filtered calls remain included. Heteroplasmy is not replaced with population AF. gnomAD-MT and HelixMTdb source dictionaries retain their full fields. A missing population record is not AF=0. Helix's unknown AN stays null. Multiallelic or symbolic alleles are not split, flattened, or assigned a biallelic frequency. Original NUMT, haplocheck and contamination limitations remain intact; no new clinical classification or limit-of-detection claim is made.

The manifest is `PASS_INDEPENDENT_MITOCHONDRIAL_REFERENCE_CONTEXT`. The workflow receipt is PASS/stage `mitochondrial-reference-context-v1`. Producer identity is an exact two-path SHA256 map for the new wrapper and preserved exporter. `source_mt_inputs` includes the four original VCF/receipt/QC/index pins; `applied_artifacts` binds every derived/copied output to source receipts, byte sizes and SHA256. Copied VCF/QC/index entries also carry their exact `source_original_pin`. All portable applied paths begin `applied-evidence/MT/`; no reference database appears in the transfer artifact list.

Preparation is non-executing:

```sh
python3 scripts/export_mitochondrial_reference_context.py prepare \
  --population-rule ORIGINAL_ROOT_REVIEWED_POPULATION_RULE \
  --rule NEW_MT_CANDIDATE_RULE \
  --output /srv/marina/outputs/mitochondrial-reference-context-v1 \
  --state-directory /srv/marina/private/mitochondrial-reference-context-state-v1
```

After exact ROOT review, status must be `ROOT_REVIEWED_EXACT_INDEPENDENT_MT_REFERENCE_RULE`:

```sh
python3 scripts/export_mitochondrial_reference_context.py run \
  --rule ROOT-reviewed-rule.json --reviewed-rule-sha256 EXACT_SHA --execute
```

Limits are one CPU, 1 GiB address space, 512 MiB output, 150 GiB free disk, one hour maximum and the unchanged absolute VM deadline `2026-11-06T03:46:48.046025Z`. A private process lock and deadline alarm bound execution. Partial failures are never promoted or retried upstream. The only small state file is `/srv/marina/private/mitochondrial-reference-context-state-v1/status.json`.

Public controls use the actual approved public MT database and four wholly synthetic single-sample VCFs/QC fixtures. They verify exact gnomAD/Helix matches, preserved unknown Helix AN, unsplit multiallelic calls, phased/missing/raw heteroplasmy fields, unchanged filters/QC bytes, receipt-bound portable paths and the absence of any nuclear-stage dependency. Production execution remains a separate reviewed action.
