param(
    [string]$OutputDir = "",
    [switch]$SkipBuild
)

$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
$Burnham = Join-Path $Root "build/manual/burnham.exe"

if (-not $SkipBuild) {
    & (Join-Path $PSScriptRoot "build-clang.ps1")
}

if (-not (Test-Path $Burnham)) {
    throw "burnham.exe was not found at $Burnham; run scripts/build-clang.ps1 first or omit -SkipBuild"
}

if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    $OutputDir = Join-Path $Root "build/package"
}

$cmake = Get-Content (Join-Path $Root "CMakeLists.txt") -Raw
$version = "0.1.0"
if ($cmake -match "project\(Burnham VERSION ([^\s\)]+)") {
    $version = $Matches[1]
}

$packageName = "burnham-$version-windows-x64"
$packageRoot = Join-Path $OutputDir $packageName
$zipPath = Join-Path $OutputDir "$packageName.zip"
$manifestPath = Join-Path $OutputDir "$packageName.manifest.json"

Remove-Item -Recurse -Force $packageRoot -ErrorAction SilentlyContinue
Remove-Item -Force $zipPath, $manifestPath -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path (Join-Path $packageRoot "bin") -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $packageRoot "docs") -Force | Out-Null

Copy-Item $Burnham (Join-Path $packageRoot "bin/burnham.exe") -Force
Copy-Item (Join-Path $Root "README.md") (Join-Path $packageRoot "README.md") -Force
Copy-Item (Join-Path $Root "docs/*.md") (Join-Path $packageRoot "docs") -Force

$hash = Get-FileHash -Algorithm SHA256 (Join-Path $packageRoot "bin/burnham.exe")
[System.IO.File]::WriteAllText((Join-Path $packageRoot "SHA256SUMS.txt"), "$($hash.Hash)  bin/burnham.exe`n", [System.Text.Encoding]::ASCII)

Compress-Archive -Path (Join-Path $packageRoot "*") -DestinationPath $zipPath -Force

$files = Get-ChildItem $packageRoot -Recurse -File | ForEach-Object {
    $_.FullName.Substring($packageRoot.Length + 1).Replace("\", "/")
}
$manifest = [ordered]@{
    schema_version = 1
    package_name = $packageName
    version = $version
    generated_at_utc = (Get-Date).ToUniversalTime().ToString("o")
    package_root = $packageRoot
    zip_path = $zipPath
    zip_bytes = (Get-Item $zipPath).Length
    files = @($files)
}
$json = $manifest | ConvertTo-Json -Depth 5
[System.IO.File]::WriteAllText($manifestPath, $json + "`n", [System.Text.Encoding]::ASCII)

Write-Host "Created package $zipPath"
Write-Host "Wrote manifest $manifestPath"