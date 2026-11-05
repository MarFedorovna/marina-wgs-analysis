# Operations and budget

## Before starting a paid VM

Prepare scripts and local manifests before provisioning. Confirm the selected machine, provider, current rate, usable account credits, storage behavior, deployment region and the mechanism for ending billing. Read-only inspection or an environment-selection form does not create an instance.

The Lambda A100 offering selected on 2026-11-05 quoted $2.39/hour with its fixed 512 GiB SSD included. At that quote:

| Running time | Quoted compute/storage cost |
| --- | ---: |
| 6 hours | $14.34 |
| 10 hours | $23.90 |
| 16 hours | $38.24 |

These are arithmetic scenarios, not processing-time promises or measured costs. Allow for provisioning/setup, uploads, downloads, provider rounding, taxes and any additional transfer charges. A $50 task budget does not establish that the provider account has $50 usable credit.

This offering cannot be stopped and resumed. Treat it as one continuous work session, preserve input files separately, and download completed results while work continues. Do not rely on exhausting credits to end a run: the console warns that doing so deletes the instance and all its data.

## First sample

Run one complete sample before scheduling the remaining samples. Record elapsed time per stage, peak system RAM, sampled GPU memory usage, free disk space and output sizes. Project the remaining cost using the measured first-sample behavior and leave time for verification and downloads. Restrict concurrency until those measurements support an increase.

A comparison with a CPU workflow must use the same reference, analysis intervals, read-group interpretation, resource versions and corresponding algorithms. To isolate caller differences, start from the same processed BAM. Agreement is a reproducibility metric, not a substitute for a truth dataset or technical validation.

## Retaining results

Keep both the original callset and explicitly filtered outputs. Preserve indexes, QC summaries, duplicate metrics, recalibration reports when used, command logs, software/image versions, and reference/resource checksums. Keep these artifacts outside GitHub if they contain sample-level information.

Download outputs after each completed sample and verify hashes on the receiving machine. A successful transfer process alone is insufficient evidence of a verified copy. Confirm that the final retained outputs can be read before terminating the only remote copy.

## Ending billing

After verifying retained results, terminate the VM through Brev/provider controls and check the resulting state. Closing SSH, the browser, Jupyter, Docker or a terminal does not end cloud billing. Delete retained cloud data only after an explicit decision and verified copies; VM termination can remove its local SSD.

## Scope

The automated technical workflow prepares aligned reads, variant calls and QC. It does not establish a diagnosis or automatically classify candidates under ACMG/AMP. Low coverage remains a limitation of the data even when compute processing succeeds.
