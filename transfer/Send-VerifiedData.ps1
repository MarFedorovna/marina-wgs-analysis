[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$ManifestPath,
    [Parameter(Mandatory = $true)][string]$HostName,
    [Parameter(Mandatory = $true)][string]$UserName,
    [Parameter(Mandatory = $true)][string]$RemoteRoot,
    [ValidateRange(1, 65535)][int]$Port = 22,
    [string]$IdentityFile,
    [switch]$Execute,
    [switch]$RestartStaging
)

$ErrorActionPreference = 'Stop'
if ($HostName -notmatch '^[A-Za-z0-9][A-Za-z0-9.-]*$') { throw 'Use a hostname, IPv4 address or SSH config alias.' }
if ($UserName -notmatch '^[A-Za-z_][A-Za-z0-9_-]*$') { throw 'Invalid SSH user name.' }
if ($RemoteRoot -notmatch '^/[A-Za-z0-9_./-]+$' -or $RemoteRoot -match '(^|/)\.\.(/|$)') {
    throw 'RemoteRoot must be a dedicated absolute POSIX directory without spaces or parent traversal.'
}
$RemoteRoot = $RemoteRoot.TrimEnd('/')
if ($RemoteRoot -in @('', '/', '/home', '/root', '/tmp', '/srv', '/mnt', '/data', '/workspace')) {
    throw 'Use a dedicated project subdirectory, not a shared root directory.'
}
$manifestFile = (Resolve-Path -LiteralPath $ManifestPath).Path
$manifest = Get-Content -LiteralPath $manifestFile -Raw -Encoding UTF8 | ConvertFrom-Json
if ($manifest.SchemaVersion -ne 1 -or @($manifest.Files).Count -eq 0) { throw 'Unsupported or empty manifest.' }
$remoteFiles = [System.Collections.Generic.List[object]]::new()
$seen = @{}
foreach ($entry in $manifest.Files) {
    $relative = [string]$entry.RelativePath
    if (-not $relative -or $relative.StartsWith('/') -or $relative.Contains('\') -or $relative -match '(^|/)\.\.(/|$)' -or $relative -match '[\x00-\x1f]') {
        throw 'Unsafe manifest destination path.'
    }
    if ($seen.ContainsKey($relative)) { throw 'Duplicate manifest destination path.' }
    $seen[$relative] = $true
    if ([string]$entry.Sha256 -cnotmatch '^[0-9a-f]{64}$') { throw 'Invalid SHA256 in manifest.' }
    $file = Get-Item -LiteralPath $entry.SourcePath
    if ($file.PSIsContainer -or [int64]$file.Length -ne [int64]$entry.SizeBytes) { throw 'Input size changed; rebuild the private manifest.' }
    $remoteFiles.Add([ordered]@{ RelativePath = $relative; SizeBytes = [int64]$entry.SizeBytes; Sha256 = [string]$entry.Sha256 })
}
$bytes = ($remoteFiles | ForEach-Object { [int64]$_.SizeBytes } | Measure-Object -Sum).Sum
Write-Output ("Plan: {0} files, {1:N2} GiB -> {2}@{3}:{4}." -f $remoteFiles.Count, ($bytes / 1GB), $UserName, $HostName, $RemoteRoot)
if (-not $Execute) {
    Write-Output 'Dry run only. Verify the SSH host key independently, authenticate with your own key, then add -Execute.'
    return
}
foreach ($command in @('ssh', 'sftp')) { $null = Get-Command $command -ErrorAction Stop }
$sshOptions = @('-o', 'BatchMode=yes', '-o', 'StrictHostKeyChecking=yes', '-o', 'ConnectTimeout=20')
if ($IdentityFile) { $sshOptions += @('-i', (Resolve-Path -LiteralPath $IdentityFile).Path) }
$target = "${UserName}@${HostName}"
function Quote-Posix([string]$Value) { return "'" + $Value.Replace("'", "'\" + "'" + "'") + "'" }
function Quote-Sftp([string]$Value) {
    if ($Value -match '[\x00-\x1f]') { throw 'Control characters cannot be passed to SFTP.' }
    return '"' + $Value.Replace('\', '\\').Replace('"', '\"') + '"'
}
function Invoke-Remote([string]$Command) {
    $result = & ssh @sshOptions -p $Port $target $Command
    if ($LASTEXITCODE -ne 0) { throw "Remote operation failed (SSH exit $LASTEXITCODE)." }
    return ($result -join "`n")
}
function Invoke-SftpBatch([string[]]$Lines, [string]$BatchPath) {
    [IO.File]::WriteAllText($BatchPath, (($Lines -join "`n") + "`n"), [Text.UTF8Encoding]::new($false))
    & sftp @sshOptions -P $Port -b $BatchPath $target
    if ($LASTEXITCODE -ne 0) { throw "SFTP failed (exit $LASTEXITCODE). Re-running resumes staging uploads." }
}

# Rehash immediately before sending, so a stale manifest never silently changes the upload.
foreach ($entry in $manifest.Files) {
    if ((Get-FileHash -LiteralPath $entry.SourcePath -Algorithm SHA256).Hash.ToLowerInvariant() -cne $entry.Sha256) {
        throw 'Input content changed; rebuild the private manifest.'
    }
}
$privateDirectory = Join-Path (Split-Path -Parent $manifestFile) 'transfer-state'
$null = New-Item -ItemType Directory -Force -Path $privateDirectory
$remotePayload = [ordered]@{ SchemaVersion = 1; Files = @($remoteFiles.ToArray()) } | ConvertTo-Json -Depth 6 -Compress
$payloadPath = Join-Path $privateDirectory 'remote-manifest.json'
[IO.File]::WriteAllText($payloadPath, $remotePayload, [Text.UTF8Encoding]::new($false))
$planId = (Get-FileHash -LiteralPath $payloadPath -Algorithm SHA256).Hash.ToLowerInvariant()
$helperLocal = Join-Path $PSScriptRoot 'verify_upload.py'
$helperId = (Get-FileHash -LiteralPath $helperLocal -Algorithm SHA256).Hash.ToLowerInvariant()
$helperRemote = "$RemoteRoot/.transfer/verify-$helperId.py"
$manifestRemote = "$RemoteRoot/.transfer/$planId.json"
$stagingRoot = "$RemoteRoot/.upload-staging/$planId"
$null = Invoke-Remote ("umask 077; command -v python3 >/dev/null && command -v sha256sum >/dev/null && mkdir -p -- " + (Quote-Posix "$RemoteRoot/.transfer") + ' ' + (Quote-Posix $stagingRoot))
$setupBatch = @(
    ('put ' + (Quote-Sftp ($helperLocal.Replace('\', '/'))) + ' ' + (Quote-Sftp $helperRemote)),
    ('put ' + (Quote-Sftp ($payloadPath.Replace('\', '/'))) + ' ' + (Quote-Sftp $manifestRemote))
)
Invoke-SftpBatch $setupBatch (Join-Path $privateDirectory 'setup.sftp')
$helperCommand = 'python3 ' + (Quote-Posix $helperRemote) + ' --root ' + (Quote-Posix $RemoteRoot) + ' --manifest ' + (Quote-Posix $manifestRemote) + ' --plan-id ' + $planId
$plan = (Invoke-Remote ($helperCommand + ' --mode plan')) | ConvertFrom-Json
$pending = @($plan.Pending)
$resume = @($plan.Resume)
if ($pending.Count -gt 0) {
    $remainingBytes = [int64]0
    foreach ($index in $pending) { $remainingBytes += [int64]$manifest.Files[[int]$index].SizeBytes }
    $availableText = Invoke-Remote ('df -Pk -- ' + (Quote-Posix $RemoteRoot) + " | tail -n 1 | awk '{print `$4}'")
    $availableBytes = [int64]$availableText.Trim() * 1024
    if ($availableBytes -lt ($remainingBytes + 10GB)) { throw 'Insufficient free space for pending uploads plus 10 GiB reserve.' }
    $uploadBatch = [System.Collections.Generic.List[string]]::new()
    foreach ($index in $pending) {
        $source = [string]$manifest.Files[[int]$index].SourcePath
        $staging = "$stagingRoot/$('{0:D4}' -f [int]$index).part"
        $verb = if ($RestartStaging -or $resume -notcontains [int]$index) { 'put' } else { 'reput' }
        $uploadBatch.Add($verb + ' ' + (Quote-Sftp ($source.Replace('\', '/'))) + ' ' + (Quote-Sftp $staging))
    }
    Invoke-SftpBatch @($uploadBatch.ToArray()) (Join-Path $privateDirectory 'upload.sftp')
}
$summary = (Invoke-Remote ($helperCommand + ' --mode finalize')) | ConvertFrom-Json
if (@($summary.Failed).Count -gt 0 -or [int]$summary.Complete -ne @($manifest.Files).Count) {
    throw 'Verification incomplete. Keep originals; rerun with -RestartStaging if a staged file is corrupt.'
}
Write-Output ("Verified {0} files with SHA256 and promoted them to their final destinations." -f $summary.Complete)
