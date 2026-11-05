[CmdletBinding()]
param(
    [Parameter(Mandatory)][string]$HostName,
    [Parameter(Mandatory)][string]$KnownHosts,
    [Parameter(Mandatory)][string]$Destination,
    [string]$StatusFile,
    [string]$RemoteRoot = '/srv/project-data',
    [string]$RemoteCode,
    [string]$SshHostName,
    [ValidateRange(1, 65535)][int]$SshPort = 22,
    [string]$HostKeyAlias,
    [string]$SshConfig = '/root/.brev/ssh_config',
    [string]$Distribution = 'Ubuntu-22.04',
    [switch]$Execute,
    [switch]$Watch,
    [switch]$AllowOneDrive,
    [ValidateRange(60, 3600)][int]$Interval = 60,
    [ValidateRange(0.01, 8)][double]$MaxHours = 8
)
$ErrorActionPreference = 'Stop'
# No credentials or real endpoints are embedded. Existing WSL SSH trust is reused.
function Convert-WindowsPathForWsl([string]$Value) {
    if ($Value -match '^[A-Za-z]:[\\/]') { return $Value.Replace('\', '/') }
    return $Value
}
$scriptPath = Join-Path $PSScriptRoot 'pull_verified_results.py'
$wslScript = & wsl.exe -d $Distribution -- wslpath -u (Convert-WindowsPathForWsl $scriptPath)
if ($LASTEXITCODE -ne 0) { throw 'Cannot resolve script in WSL.' }
$wslScript = "$wslScript".Trim()
$arguments = @('-d', $Distribution, '--', 'python3', $wslScript,
    '--host', $HostName, '--ssh-config', (Convert-WindowsPathForWsl $SshConfig),
    '--known-hosts', (Convert-WindowsPathForWsl $KnownHosts),
    '--remote-root', $RemoteRoot, '--destination', (Convert-WindowsPathForWsl $Destination),
    '--interval', "$Interval", '--max-hours', "$MaxHours")
if ($RemoteCode) { $arguments += @('--remote-code', $RemoteCode) }
if ($SshHostName) { $arguments += @('--ssh-hostname', $SshHostName, '--ssh-port', "$SshPort", '--host-key-alias', $HostKeyAlias) }
if ($StatusFile) { $arguments += @('--status-file', (Convert-WindowsPathForWsl $StatusFile)) }
if ($Execute) { $arguments += '--execute' }
if ($Watch) { $arguments += '--watch' }
if ($AllowOneDrive) { $arguments += '--allow-onedrive' }
& wsl.exe @arguments
if ($LASTEXITCODE -ne 0) { throw "Verified backup exited with code $LASTEXITCODE. Existing completed backups remain intact." }
