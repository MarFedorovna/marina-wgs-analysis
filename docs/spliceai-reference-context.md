# Local SpliceAI gene-context predictions

This supplement calculates actual SpliceAI 1.3.1 scores locally from the five official model weights, the accepted GRCh38 FASTA and official bundled GENCODE v24 gene-union annotation. It does not assign these scores to VEP GENCODE49 transcripts or replace the original VEP output. Model code and weights retain their publisher licenses; the project has explicit user confirmation of noncommercial research.

The fixed approved HCM panel contains 39 genes, including TNNT2 and MT-TI. Thirty-eight nuclear gene models are present in the bundled annotation. MT-TI is not applicable to this nuclear model and remains missing with scope status. The read-only audit found 2326 accepted source alleles in these gene spans: 1818 SNVs, 489 indels and 19 complex/unsupported alleles. Outside-panel alleles are not predicted and must not receive invented zeros. The old MyVariant public cache has no populated SpliceAI fields in its 572031 scanned public records.

Scientific parameters are distance 500 bp and mask 0 (unmasked raw model delta scores). The official result is `ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`. DS values are the official two-decimal scores; DP values are signed genomic offsets from the variant position. Exact REF/ALT, gene, model/annotation hashes, assembly, distance, mask, runtime and scope are retained. Unsupported complex, symbolic/non-ACGT, overlength, reference-mismatching and unavailable-window cases have explicit status and empty scores. A real `0.00` is retained only when a complete score vector was calculated.

The fixed public MYBPC3 control `chr11:47332146 T>G` produces `G|MYBPC3|0.02|0.60|0.03|0.01|-164|-22|-50|-301`. A direct Keras inference wrapper feeds the official `get_delta_scores` algorithm. It was compared with original Keras `predict` on this public control: zero maximum probability difference and identical DS/DP. It avoids repeated prediction-dispatch overhead and uses no changed model weights or inference formulas. The measured median was 2.402 seconds per allele with two CPU threads, approximately 93 minutes for the selected scope before additional safeguards. No GPU download or patient API request is necessary.

`scripts/continue_spliceai_reference_context.py prepare` pins the accepted stage05 VCF and receipt, exact public39 panel, model/annotation hashes, existing runtime wheel manifest, reference index, original binding, read-only audit and public control. It writes a candidate without execution. Root reviews the exact rule and changes its status to `ROOT_REVIEWED_EXACT_SPLICEAI_REFERENCE_RULE`; only then may the existing isolated runtime execute:

```sh
/srv/marina/resources/spliceai-research-runtime-v1/venv/bin/python scripts/continue_spliceai_reference_context.py run --rule ROOT-reviewed-rule.json --reviewed-rule-sha256 EXACT_REVIEWED_SHA --execute
```

The source FASTA SHA comes from the accepted original reference binding; size/mtime and index are pinned, and every scored allele has an actual local REF alignment check. No BAM is rehashed. The existing runtime and primary workflow are not modified. A private lock and identity-bound SQLite checkpoint prevent duplicate writers and allow exact-source resumption. Limits are 2 CPU threads, 8 GiB address space, 1 GiB checkpoint/output each, 150 GiB disk reserve, 3.5 hours maximum and the unchanged hard VM deadline `2026-11-06T03:46:48.046025Z`. Progress contains counts and timing only.

Success is `/srv/marina/outputs/spliceai-reference-context-v1/{predictions.tsv.gz,manifest.json,workflow-receipt.json}` with manifest `PASS_LOCAL_SPLICEAI_REFERENCE_CONTEXT`, workflow `PASS` and stage `spliceai-reference-context-v1`. The state is `/srv/marina/private/spliceai-reference-context-state-v1/status.json`. Placeholder-only results cannot receive PASS.

The new VCF codec appends these exact thirteen fields after the original and 33 public-context fields:

```text
spliceai_status spliceai_gene spliceai_DS_AG spliceai_DS_AL spliceai_DS_DG spliceai_DS_DL spliceai_DP_AG spliceai_DP_AL spliceai_DP_DG spliceai_DP_DL spliceai_prediction_json spliceai_source_json spliceai_transcript_scope
```

The join uses exact GRCh38 alleles and gene symbols. Scalar row scores require an unambiguous same-gene model result. The complete applicable allele/gene predictions remain in `spliceai_prediction_json`; mapping to a different gene or to an unverified transcript stays explicit. The original master fields, original ranking and classification remain unchanged. The cohort VCF includes standard SpliceAI INFO with real scores and offsets, while per-sample packets retain only that sample's rows. This new dependency is separate from the preserved earlier VCF codec and must be reviewed before deployment.
