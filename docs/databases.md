# Install annotation databases directly on the server

Run from this reviewed repository checkout **on the Linux x86_64 server**, with Python 3.9+, GNU coreutils and working Docker access (directly or through non-interactive sudo):

```bash
bash publicstages/01_prepare_databases.sh
```

This standalone command has no patient/BAM/workflow-config prerequisites. It first stages the public ClinVar/ClinGen/HPO/HGNC/MANE files, then downloads the exact **Ensembl VEP 115 GRCh38 indexed cache** directly to the server, verifies it, and prepares the final resource bundle. It pulls only the pinned official VEP 115.2 container. It does not change the GPU driver, run patient analysis, or send patient variants to any API.

For a background run that survives a disconnected terminal:

```bash
mkdir -p /srv/marina/logs
nohup bash publicstages/01_prepare_databases.sh \
  > /srv/marina/logs/annotation-databases.log 2>&1 < /dev/null &
echo $!
tail -f /srv/marina/logs/annotation-databases.log
```

Check that the existing process has finished before starting another foreground command. A download lock prevents simultaneous writers. The transfer uses 16 bounded HTTPS ranges by default (`MARINA_DOWNLOAD_WORKERS=8` selects fewer connections). Every range requires HTTP 206, exact Content-Range/length and unchanged ETag/Last-Modified. Completed 64 MiB ranges have SHA256 receipts. Restart the **same command** to verify and resume owned partial ranges; no completed database is overwritten. The direct-server download directory is `/srv/marina/resources/.downloads/`. No workstation/OneDrive staging is used by this command.

The cache archive is exactly 25,272,957,721 bytes (23.54 GiB). The publisher's BSD rotating checksum is 01402 with 24,680,623 blocks. This is a weak, noncryptographic publisher checksum, not a publisher SHA256. Preparation verifies it plus exact size and records a newly computed SHA256 for transfer/reuse integrity. Cache extraction rejects links, special files, duplicate/traversal paths, excessive expanded size and insufficient disk reserve. At least 150 GiB free is retained for later work, in addition to the expansion requirement. Failed partials remain for diagnosis. A resource download invocation has a six-hour cache-transfer deadline; restart is explicit and does not extend the VM budget automatically.

The frozen base recipe is `config/annotation-base-resources.example.json`. The final recipe, `config/annotation-resources.example.json`, adds HPO/HGNC/MANE through a versioned v2 bundle. It reuses verified immutable base files and cache with hardlinks instead of a second 25 GB download or extraction. A standalone small-source cache allows the individual public databases to finish before the large VEP cache. Their individual `resource.json` receipts are not a declaration that the full annotation stage is ready.

**READY** means this exact file exists with `status: PASS`, and all its listed artifact hashes still verify:

```text
/srv/marina/resources/annotation/vep115.2-cache115-clinvar20261025-clingen20261105-publicextras-v2/resource-manifest.json
```

The manifest records source URL, version/snapshot, retrieval time, byte count, SHA256, publisher checks, cache info, container digest and actual availability. A `.partial-*` directory, download progress or installed image alone is **not READY**. Small-source receipts live under `annotation/.public-source-cache/`. All these public databases stay on the server; GitHub stores only their installer, configuration and documentation.

Current pinned resources are ClinVar GRCh38 2026-10-25; ClinGen gene validity/dosage snapshots 2026-11-05; HPO ontology and gene associations v2026-09-01 (published SHA256); HGNC GCS generation 1793369077513436 from 2026-10-30 (publisher MD5); MANE GRCh38 v1.4 summary; and VEP/cache 115.2/115. The dated ClinGen endpoints are dynamic: a different file date causes a stop, requiring an explicitly reviewed new snapshot recipe. gnomAD AF, SIFT/PolyPhen and transcript/MANE component versions come from the verified cache's `info.txt`; they are not silently relabelled as a newer release.

The [audited source registry](../config/evidence-sources.example.json) lists additional manual, licensed, historical and separate analysis branches. Availability is reported explicitly. A current dbNSFP bulk download requires registration or a commercial subscription; HGMD Professional and unrestricted OMIM dumps are not assumed available. Author FATHMM v2.3 has a public legacy precomputed SQL database and GPLv3 tool, but requires an isolated legacy local backend/export; no database dump is executed into the user's databases. CADD whole-genome v1.7 scores are about 81.5 GiB and have separate use terms; neither their full download nor a completed four-predictor consensus is implied by the core bundle. Authorized local exact-allele exports can be attached with explicit version/SHA/availability. Missing scores never become benign predictions.

The selected VM continues charging while downloading or after this process exits. Direct Ensembl network throughput can vary greatly; progress is observed, not a promised completion time. Keep the account budget and free space under review and retain verified patient results separately before terminating the VM.

Official sources: [VEP download](https://sep2025.archive.ensembl.org/info/docs/tools/vep/script/vep_download.html), [VEP 115 cache contents](https://sep2025.archive.ensembl.org/info/docs/tools/vep/script/vep_cache.html), [ClinVar downloads](https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/), [ClinGen downloads](https://search.clinicalgenome.org/kb/downloads), [HPO releases](https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01), [HGNC downloads](https://www.genenames.org/download/), [MANE v1.4](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/).
