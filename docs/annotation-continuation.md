# Detached research annotation and report

This continuation starts from a completed, hash-verified stage 5. It does not call
variants again or change stages 3–5. The accepted VCF scope is inherited; an
autosomal callset does not establish coverage of X, Y, mitochondrial variants or
CNV/SV. The separate mitochondrial workflow has its own resources and receipts.

Prepare an isolated code package on the server after reviewing its sources:

```bash
python3 scripts/annotation_queue.py prepare \
  --base-config /srv/marina/private/workflow-v1.json \
  --runtime-root /srv/marina/private/annotation-runtime-v1 \
  --upstream-code-root /home/ubuntu/workspace/marina-wgs-analysis \
  --created-utc ACTUAL_VM_CREATION_TIME_WITH_TIMEZONE
```

The creation timestamp must come from the actual billed instance, not the time of
a restart. The allowance is capped at $40 at $2.39/hour within the $50 total
budget. Stopping this process does not stop VM billing.

Preparation checks stage-5 artifacts and the analysis FASTA, copies an explicit
code allowlist, and seals a separate private annotation configuration. Existing
sealed artifacts may only be reused with identical bytes. The private bindings
record the original configuration, stage-5 receipt, FASTA, immutable code package,
upstream code, AlphaMissense manifest and expected VEP resource recipe. They are
not committed to GitHub. An internal review approval file must contain the exact
`binding_sha256` printed by preparation; a changed binding requires new review.

Start or resume the prepared queue with one command:

```bash
nohup python3 -u /srv/marina/private/annotation-runtime-v1/scripts/annotation_queue.py run --execute \
  >>/srv/marina/logs/annotation-queue.log 2>&1 </dev/null &
```

An exclusive lock prevents duplicate queues. The queue waits for the exact code
review approval, the verified offline VEP 115 GRCh38 resource manifest, and the
required verified AlphaMissense canonical/isoform supplement. It checks the
unchanged VM allowance while waiting and running. It then:

1. Runs twenty public ClinVar variants, including MYH7/MYBPC3 and consequence
   controls, with synthetic ALT/reference/unknown/filtered genotypes through the
   actual offline VEP module. Exact ClinVar IDs, SIFT/PolyPhen population, genotype
   distinctions, unique-allele counts and ranking size limits must pass.
2. Runs stage 6 against all soft-labelled stage-5 records. Resource availability,
   transcript identity mismatches and missing evidence stay explicit. FATHMM or
   other unavailable optional sources are not represented as loaded.
3. Runs stage 7 with master/Supplementary tables, per-sample and cohort top10/top20,
   urgent low-quality ClinVar P/LP review tables, figures, HTML and local DOCX.

Inspect progress without starting computation:

```bash
python3 /srv/marina/private/annotation-runtime-v1/scripts/annotation_queue.py run --status
tail -n 30 /srv/marina/logs/annotation-queue.log
```

Detailed child logs are under `/srv/marina/logs/annotation-continuation/`.
`RESEARCH_REPORT_COMPLETE` means the technical outputs passed their receipts.
Model labels, external ClinVar classifications and research ranks remain separate
from clinical classification. ACMG criteria require documented manual assessment;
no diagnosis, GoF assignment or automatic pathogenicity class is produced.

The queue stops its own active child process group and recorded VEP containers at
the allowance deadline. It does not stop unrelated jobs, remove failed partials,
or terminate the VM. Repeating the command reuses only valid completed receipts;
changed code/configuration or upstream artifacts fail closed.
