# Bounded native SFTP parallel probe

`transfer/native_sftp_parallel_probe.py` reads exactly four MiB: one serial MiB, followed by three MiB split over four or eight independent SFTP connections. It reuses the existing strict host-key-pinned connection helper and managed key in memory, opens the source read-only, and checks size and mtime before and after the parallel ranges. It does not create or replace patient files or calculate patient data hashes. Connection setup is excluded from the parallel payload speed and reported separately.

```powershell
& '.local/native-ssh-venv/Scripts/python.exe' transfer/native_sftp_parallel_probe.py --config PRIVATE-CONFIG.json --receipt PRIVATE-PROBE-RECEIPT.json --connections 4 --execute
```

Without `--execute`, the CLI performs no network or key read. Keep endpoint configurations and probe receipts outside version control.

The actual 5 November 2026 run read 4,194,304 bytes total. The one-MiB serial sample took 48.0 seconds (0.02083 MiB/s); four simultaneous 0.75-MiB ranges took 47.109 seconds (aggregate 0.06368 MiB/s). Although this is 3.06 times the short serial sample, it is below the existing useful downloader's approximately 0.14 MiB/s long-run rate. The existing writer also completed 16 MiB between 18:19:31.74 and 18:20:49.64 UTC, approximately 0.205 MiB/s, while this probe was active. These measurements do not support replacing the running downloader with parallel small ranges. All probe connections closed, and the existing writer and validated partial remained active.

A larger per-connection range would require a separate, explicitly bounded benchmark to determine useful sustained throughput; the four-MiB probe alone does not establish that result. No accelerated downloader was started based on these measurements.

Four unit tests cover exact total/disjoint read ranges for both supported connection counts, invalid concurrency rejection before connection, source-size mismatch rejection before any read, and truncated baseline rejection before parallel connections. Tests use synthetic in-memory reads and no network or SSH key.
