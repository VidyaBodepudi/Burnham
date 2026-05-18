param(
    [int]$BenchRecordCount = 120,
    [int]$BenchRepeatCount = 1,
    [switch]$SkipBench,
    [switch]$SkipPackage
)

$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
$Phase5Dir = Join-Path $Root "build/phase5"
Remove-Item -Recurse -Force $Phase5Dir -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $Phase5Dir -Force | Out-Null

$script:Steps = @()

function Invoke-Phase5Step {
    param(
        [string]$Name,
        [scriptblock]$Action
    )

    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    try {
        & $Action
        $timer.Stop()
        $script:Steps += [ordered]@{
            name = $Name
            status = "passed"
            milliseconds = [Math]::Round($timer.Elapsed.TotalMilliseconds, 3)
        }
    } catch {
        $timer.Stop()
        $script:Steps += [ordered]@{
            name = $Name
            status = "failed"
            milliseconds = [Math]::Round($timer.Elapsed.TotalMilliseconds, 3)
            error = $_.Exception.Message
        }
        throw
    }
}

Invoke-Phase5Step -Name "unit-tests" -Action {
    & (Join-Path $PSScriptRoot "test-clang.ps1")
}

$benchmarkPath = Join-Path $Phase5Dir "phase4-smoke-benchmark.json"
if (-not $SkipBench) {
    Invoke-Phase5Step -Name "benchmark-smoke" -Action {
        & (Join-Path $PSScriptRoot "bench-phase4.ps1") -RecordCount $BenchRecordCount -RepeatCount $BenchRepeatCount -WarmupCount 0 -OutputPath $benchmarkPath -SkipBuild
    }
    Invoke-Phase5Step -Name "benchmark-compare-smoke" -Action {
        & (Join-Path $PSScriptRoot "compare-phase4.ps1") -BaselinePath $benchmarkPath -CurrentPath $benchmarkPath -RegressionThresholdPercent 0 -FailOnRegression
    }
}

$packageOutputDir = Join-Path $Phase5Dir "package"
if (-not $SkipPackage) {
    Invoke-Phase5Step -Name "package-smoke" -Action {
        & (Join-Path $PSScriptRoot "package-release.ps1") -OutputDir $packageOutputDir -SkipBuild
        $zip = Get-ChildItem $packageOutputDir -Filter "*.zip" | Select-Object -First 1
        $manifest = Get-ChildItem $packageOutputDir -Filter "*.manifest.json" | Select-Object -First 1
        if ($null -eq $zip -or $null -eq $manifest) {
            throw "package zip or manifest was not created"
        }
    }
}

$summary = [ordered]@{
    schema_version = 1
    phase = "phase5-verification-ci-packaging"
    generated_at_utc = (Get-Date).ToUniversalTime().ToString("o")
    bench_record_count = $BenchRecordCount
    bench_repeat_count = $BenchRepeatCount
    benchmark_path = if ($SkipBench) { $null } else { $benchmarkPath }
    package_output_dir = if ($SkipPackage) { $null } else { $packageOutputDir }
    steps = $script:Steps
}

$summaryPath = Join-Path $Phase5Dir "verification.json"
$json = $summary | ConvertTo-Json -Depth 6
[System.IO.File]::WriteAllText($summaryPath, $json + "`n", [System.Text.Encoding]::ASCII)

Write-Host "Phase 5 verification passed"
Write-Host "Wrote verification summary to $summaryPath"