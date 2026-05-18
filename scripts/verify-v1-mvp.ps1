param(
    [int]$BenchRecordCount = 80,
    [int]$BenchRepeatCount = 1,
    [switch]$SkipPhase5
)

$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
$MvpDir = Join-Path $Root "build/v1-mvp"
$Burnham = Join-Path $Root "build/manual/burnham.exe"
Remove-Item -Recurse -Force $MvpDir -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $MvpDir -Force | Out-Null

$script:Steps = @()

function Invoke-MvpStep {
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

if (-not $SkipPhase5) {
    Invoke-MvpStep -Name "phase5-verification" -Action {
        & (Join-Path $PSScriptRoot "verify-phase5.ps1") -BenchRecordCount $BenchRecordCount -BenchRepeatCount $BenchRepeatCount
    }
} else {
    Invoke-MvpStep -Name "direct-clang-build" -Action {
        & (Join-Path $PSScriptRoot "build-clang.ps1")
    }
}

Invoke-MvpStep -Name "doctor-htslib-status" -Action {
    if (-not (Test-Path $Burnham)) {
        throw "burnham.exe was not found at $Burnham"
    }
    $doctorOutput = & $Burnham doctor 2>&1
    $doctorText = $doctorOutput -join "`n"
    [System.IO.File]::WriteAllText((Join-Path $MvpDir "doctor.txt"), $doctorText + "`n", [System.Text.Encoding]::ASCII)
    if ($LASTEXITCODE -ne 0) {
        throw "doctor failed with exit code $LASTEXITCODE"
    }
    if ($doctorText -notmatch "htslib:") {
        throw "doctor did not report htslib status"
    }
}

Invoke-MvpStep -Name "vcpkg-manifest" -Action {
    $manifestPath = Join-Path $Root "vcpkg.json"
    if (-not (Test-Path $manifestPath)) {
        throw "vcpkg.json is missing"
    }
    $manifest = Get-Content $manifestPath -Raw | ConvertFrom-Json
    if ($manifest.dependencies -notcontains "htslib") {
        throw "vcpkg.json does not declare htslib"
    }
}

Invoke-MvpStep -Name "vcpkg-cmake-preset" -Action {
    $presetPath = Join-Path $Root "CMakePresets.json"
    $presets = Get-Content $presetPath -Raw | ConvertFrom-Json
    $presetNames = @($presets.configurePresets | ForEach-Object { $_.name })
    if ($presetNames -notcontains "vcpkg-debug") {
        throw "CMakePresets.json does not include vcpkg-debug"
    }
    $cmakeLists = Get-Content (Join-Path $Root "CMakeLists.txt") -Raw
    if ($cmakeLists -notmatch "BURNHAM_ENABLE_HTSLIB") {
        throw "CMakeLists.txt does not expose BURNHAM_ENABLE_HTSLIB"
    }
}

Invoke-MvpStep -Name "binary-transform-guard" -Action {
    $fakeBam = Join-Path $MvpDir "fake.bam"
    $outSam = Join-Path $MvpDir "out.sam"
    [System.IO.File]::WriteAllText($fakeBam, "not a bam`n", [System.Text.Encoding]::ASCII)
    $previousPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    try {
        $guardOutput = & $Burnham sort --input $fakeBam --output $outSam --force 2>&1
        $exitCode = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $previousPreference
    }
    $guardText = $guardOutput -join "`n"
    [System.IO.File]::WriteAllText((Join-Path $MvpDir "binary-transform-guard.txt"), $guardText + "`n", [System.Text.Encoding]::ASCII)
    if ($exitCode -eq 0) {
        throw "sort unexpectedly accepted BAM input"
    }
    if ($guardText -notmatch "text SAM only") {
        throw "binary transform guard did not explain the text SAM boundary"
    }
}

$cmakeCommand = Get-Command cmake -ErrorAction SilentlyContinue
$vcpkgCommand = Get-Command vcpkg -ErrorAction SilentlyContinue
$doctorText = if (Test-Path (Join-Path $MvpDir "doctor.txt")) { Get-Content (Join-Path $MvpDir "doctor.txt") -Raw } else { "" }

$summary = [ordered]@{
    schema_version = 1
    milestone = "v1-mvp-production-io-boundary"
    generated_at_utc = (Get-Date).ToUniversalTime().ToString("o")
    cmake_available = $null -ne $cmakeCommand
    vcpkg_available = $null -ne $vcpkgCommand
    htslib_doctor_output = $doctorText.Trim()
    steps = $script:Steps
    note = "This verifies the v1 MVP text foundation plus optional htslib/vcpkg integration boundary. Full BAM/CRAM transformation ports require htslib-enabled CMake/vcpkg builds and follow-on command ports."
}

$summaryPath = Join-Path $MvpDir "verification.json"
$json = $summary | ConvertTo-Json -Depth 6
[System.IO.File]::WriteAllText($summaryPath, $json + "`n", [System.Text.Encoding]::ASCII)

Write-Host "V1 MVP verification passed"
Write-Host "Wrote verification summary to $summaryPath"