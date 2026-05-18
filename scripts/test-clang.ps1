$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
Push-Location $Root
try {
    & "$PSScriptRoot/build-clang.ps1"
    & "build/manual/burnham_tests.exe"
} finally {
    Pop-Location
}
