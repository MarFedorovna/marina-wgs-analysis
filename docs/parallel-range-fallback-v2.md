# Transfer-only fallback v2

The new standalone transfer base accepts the explicitly authorized local download
ceiling `2026-11-06T04:50:00Z`. Earlier transfer sources and scientific analysis
controllers keep their original deadline. Preparing this fallback does not stop,
restart or modify the running downloader.

It reuses the same captured final VCF/TBI plan, native trust, sixteen range files,
source identities, remote mtime and 15GiB reserve. Only a fresh reviewed plan's
deadline changes. Existing completed files are never overwritten.

Before assembly, all sixteen range receipts must match the exact source/layout.
The persistent assembly marker is written before appending. On restart, the
prefix must remain within the original prefix/total bounds; every already
appended byte is compared with the retained ordered ranges without computing a
new patient-file hash. The rest is appended from its actual offset, never
truncated. Missing receipts, changed layout or a differing appended tail fail
closed and preserve every file.

Earlier range-v1 code did not persist this marker. A larger legacy prefix is
therefore rejected unless ROOT provides reviewed evidence authorizing that
assembly phase. The fallback does not infer ownership from file size alone.
