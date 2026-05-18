$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
$BuildDir = Join-Path $Root "build/manual"
$Compiler = "C:\Program Files\LLVM\bin\clang++.exe"

if (-not (Test-Path $Compiler)) {
    throw "clang++ was not found at $Compiler"
}

New-Item -ItemType Directory -Path $BuildDir -Force | Out-Null

$Sources = @(
    "src/chain_index.cpp",
    "src/commands.cpp",
    "src/cigar.cpp",
    "src/fs.cpp",
    "src/hts.cpp",
    "src/sample_index.cpp",
    "src/sam.cpp",
    "src/text.cpp",
    "src/vcf.cpp"
)

Push-Location $Root
try {
    & $Compiler -std=c++20 -I "include" @Sources "src/main.cpp" -o "build/manual/burnham.exe"
    & $Compiler -std=c++20 -I "include" @Sources "tests/unit_tests.cpp" -o "build/manual/burnham_tests.exe"
} finally {
    Pop-Location
}

Write-Host "Built build/manual/burnham.exe and build/manual/burnham_tests.exe"
