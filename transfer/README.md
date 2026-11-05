# Verified data transfer

These scripts transfer inputs separately from the GitHub repository. They require Windows PowerShell 5.1 or newer, Windows OpenSSH (`ssh` and `sftp`), and Python 3 on the Linux destination. They install no local software.

`New-TransferManifest.ps1` hashes FASTQ files in `raw_data/` and existing files in `data/`. It excludes filenames beginning with `Copy` or `Копия`. Add `-IncludeProvenance` to include Markdown reports directly under the data root. Review the private manifest before uploading. SQL databases are not discovered or selected automatically; add an explicitly approved source as a manifest record if one exists.

```powershell
./transfer/New-TransferManifest.ps1 `
  -DataRoot 'C:/path/to/private-data' `
  -OutputPath './.local/transfer-manifest.json' `
  -IncludeProvenance

# Prints the plan without contacting the server or sending data.
./transfer/Send-VerifiedData.ps1 `
  -ManifestPath './.local/transfer-manifest.json' `
  -HostName 'confirmed-server.example' -UserName 'ubuntu' `
  -RemoteRoot '/home/ubuntu/project-data'
```

Use the host, SSH port and user returned by the actual deployed instance. Independently verify and trust its SSH host key first. The transfer requires an authenticated SSH agent or `-IdentityFile`; it never disables host key checks or stores credentials. Add `-Execute` only when the destination and authorized data have been reviewed.

Uploads first enter a dedicated `.upload-staging/` directory. New files use `put`; partial staging files use `reput` to resume after interruption. A staging file that is already complete and matches its SHA256 is finalized without uploading it again. On the server, `verify_upload.py` verifies the expected manifest, file sizes and SHA256 digests before moving completed files into `inputs/raw`, `inputs/legacy`, and optionally `inputs/provenance`. A repeat run verifies existing completed files and skips them. It refuses to overwrite a final destination with differing content. If corruption is reported, rerun with `-RestartStaging` to overwrite only unverified staging files.

The manifest, SFTP batch files and local state contain real paths and identifiers. Keep them in `.local/`, which must remain ignored by Git. Source inputs are never moved or deleted. Keep independent copies of all outputs before terminating an ephemeral server.
