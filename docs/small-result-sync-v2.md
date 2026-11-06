# Small-result sync for the VCF-first delivery

Revision2 requires the completed canonical `cohort-sample-views-v1` producer.
The four historical CSV/context producers are not completion prerequisites for
the direct VCF-first pipeline. Independently completed MT reference results and
existing completed QC/MT files may also be copied; they do not delay nuclear views.

```powershell
python transfer/sync_completed_small_results_v2.py --profile "$env:SMALL_RESULTS_PROFILE"
```

The default performs local inventory only. It does not open authentication
configuration, keys, a network connection or an importer, and creates no transfer
state. After review, append `--execute --watch` to sync completed outputs with
the existing pinned-host SFTP and size-only importer.

The private revision2 profile reuses the authentication reference, destination,
staging partials, unchanged absolute deadline and 15GiB reserve. It has a new own
state/lock directory. Before an operational handoff, the previous owned worker
must have exited and released its lock; this script never stops another worker.
No BAM duplicate, gVCF or reference database is selected. Historical expected
server hashes are recorded; patient file bytes are not newly hashed locally.
