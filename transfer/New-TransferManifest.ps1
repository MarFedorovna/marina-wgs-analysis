[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$DataRoot,
    [Parameter(Mandatory = $true)][string]$OutputPath,
    [switch]$IncludeProvenance,
    # Build the localized Copy prefix without a source-encoding dependency in Windows PowerShell 5.1.
    [string[]]$ExcludeName = @('Copy*', (-join @([char]0x041A, [char]0x043E, [char]0x043F, [char]0x0438, [char]0x044F, '*')))
)

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path -LiteralPath $DataRoot).Path.TrimEnd('\', '/')
$records = [System.Collections.Generic.List[object]]::new()
$categories = @(
    @{ Source = 'raw_data'; Destination = 'inputs/raw'; Filter = '*.fastq.gz' },
    @{ Source = 'data'; Destination = 'inputs/legacy'; Filter = '*' }
)
foreach ($category in $categories) {
    $directory = Join-Path $root $category.Source
    if (-not (Test-Path -LiteralPath $directory -PathType Container)) { continue }
    foreach ($file in (Get-ChildItem -LiteralPath $directory -File -Filter $category.Filter | Sort-Object Name)) {
        $excluded = $false
        foreach ($pattern in $ExcludeName) { if ($file.Name -like $pattern) { $excluded = $true; break } }
        if ($excluded) { continue }
        $records.Add([pscustomobject]@{
            SourcePath = $file.FullName
            RelativePath = "$($category.Destination)/$($file.Name)"
            SizeBytes = [int64]$file.Length
            Sha256 = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
            LastWriteTimeUtc = $file.LastWriteTimeUtc.ToString('o')
        })
    }
}
if ($IncludeProvenance) {
    foreach ($file in (Get-ChildItem -LiteralPath $root -File -Filter '*.md' | Sort-Object Name)) {
        $records.Add([pscustomobject]@{
            SourcePath = $file.FullName
            RelativePath = "inputs/provenance/$($file.Name)"
            SizeBytes = [int64]$file.Length
            Sha256 = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
            LastWriteTimeUtc = $file.LastWriteTimeUtc.ToString('o')
        })
    }
}
if ($records.Count -eq 0) { throw 'No matching input files were found.' }
if (($records.RelativePath | Select-Object -Unique).Count -ne $records.Count) { throw 'Duplicate remote relative paths.' }
$manifest = [ordered]@{
    SchemaVersion = 1
    CreatedUtc = [DateTime]::UtcNow.ToString('o')
    Files = @($records.ToArray())
}
$destination = [IO.Path]::GetFullPath($OutputPath)
$null = New-Item -ItemType Directory -Force -Path ([IO.Path]::GetDirectoryName($destination))
[IO.File]::WriteAllText($destination, ($manifest | ConvertTo-Json -Depth 6), [Text.UTF8Encoding]::new($false))
$bytes = ($records | Measure-Object -Property SizeBytes -Sum).Sum
Write-Output ("Prepared private manifest: {0} files, {1:N2} GiB. No data was uploaded." -f $records.Count, ($bytes / 1GB))
