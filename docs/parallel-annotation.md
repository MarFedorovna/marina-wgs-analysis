# Parallel execution of the primary annotation

`scripts/parallel_annotation.py` wraps an unchanged, checksum-verified annotation
runtime. It intercepts only the frozen offline VEP command. It partitions the
sequence-allele input into autosome blocks, runs at most five simultaneous jobs
with four CPU cores and 16 GiB each, and merges them in the original input order.
The existing frozen functions then retain nonsequence records, add AlphaMissense,
make tables, and create the report. This remains one primary annotation followed
by its usual tables and report.

The wrapper retains the immutable VEP image, cache, reference, custom sources,
four VEP forks, scientific flags, and `--network=none`. The total allocation is
20 CPU cores and 80 GiB, with no GPU requirement. Speedup depends on chromosome
sizes and serial work already completed; five concurrent jobs do not promise a
fivefold improvement.

## Public proof before private execution

Run `prove` against an existing twenty-variant public-control receipt. It verifies
the public input checksum and synthetic sample names, runs serial annotation and
two-worker scatter with the same frozen configuration and sources, and compares
every final VCF data line, including complete INFO/CSQ and FORMAT/sample fields.
Only execution time and input/output command paths in explicit provenance header
lines are normalized. Scientific headers must match. Transcript/sample tables,
intersections, counts, and resource metadata must match, and both bgzip/TBI outputs
must support indexed access. At least two chromosomes are required. This proof is
public-only and does not stop a live private annotation.

`prepare` binds a candidate rule to the exact wrapper SHA256, public proof,
annotation binding, original budget deadline, and current controller/queue PID,
start time, and command arguments. The existing frozen package remains unchanged.
The rule is private: it includes deployment paths and process identities.

## Explicit user handoff

The generic server command is:

```sh
python3 /private/reviewed/scripts/parallel_annotation.py run --rule /private/reviewed/rule.json --reviewed-rule-sha256 EXACT_REVIEWED_SHA256 --execute
```

Without `--execute`, the command verifies the plan and leaves the serial job
running. The exact rule and proof must be reviewed before the executable command
is handed to the user. Preparing or reviewing this command does not authorize an
agent to launch it on the user's behalf.

If verified annotation is already complete, the command exits without stopping
the existing report queue. If verified annotation and report are complete, it
returns the report path. Otherwise the user launch stops only the exact pinned
controller and its queue, waits for clean queue exit and the annotation lock,
and refuses execution if another annotation/report process remains. It preserves
the partial serial directory. An incomplete serial VEP file is not a completed
checkpoint and is not silently adopted. Completed upstream stages are verified
and reused without rerunning calling or filtering.

The wrapper writes stages 06 and 07 to their original primary paths so existing
downstream watchers can find them. The new stage identity includes this wrapper,
the unchanged frozen helper hashes, execution parameters, the exact rule, and
the public proof. `parallel-annotation-producer.json` records every shard and the
new producer; it does not claim that the old sealed runtime performed scattering.
Watchers that additionally pin the old annotation producer must be explicitly
reviewed against this new producer before their rules are changed.

The original VM creation time, rate, and spending cap remain binding. This
controller does not extend the deadline or stop VM billing. Failed partial
directories remain for inspection; there is no automatic replacement of completed
stages, automatic diagnosis, canonical FATHMM inference, or ACMG/PVS1 assignment.
