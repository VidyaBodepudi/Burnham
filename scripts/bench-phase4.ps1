param(
    [int]$RecordCount = 2000,
    [int]$RepeatCount = 3,
    [int]$WarmupCount = 1,
    [string]$OutputPath = "",
    [switch]$SkipBuild
)

$ErrorActionPreference = "Stop"

if ($RecordCount -lt 10) {
    throw "RecordCount must be at least 10"
}
if ($RepeatCount -lt 1) {
    throw "RepeatCount must be at least 1"
}
if ($WarmupCount -lt 0) {
    throw "WarmupCount must be zero or greater"
}

$Root = Split-Path -Parent $PSScriptRoot
$BuildDir = Join-Path $Root "build/manual"
$BenchDir = Join-Path $Root "build/bench/phase4"
$Burnham = Join-Path $BuildDir "burnham.exe"

if (-not $SkipBuild) {
    & (Join-Path $PSScriptRoot "build-clang.ps1")
}

if (-not (Test-Path $Burnham)) {
    throw "burnham.exe was not found at $Burnham; run scripts/build-clang.ps1 first or omit -SkipBuild"
}

Remove-Item -Recurse -Force $BenchDir -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $BenchDir -Force | Out-Null

if ([string]::IsNullOrWhiteSpace($OutputPath)) {
    $OutputPath = Join-Path $BenchDir "phase4-benchmark.json"
}

function Write-AsciiFile {
    param(
        [string]$Path,
        [string]$Content
    )
    [System.IO.File]::WriteAllText($Path, $Content, [System.Text.Encoding]::ASCII)
}

function Get-FileSizeOrZero {
    param([string]$Path)
    if ([string]::IsNullOrWhiteSpace($Path) -or -not (Test-Path $Path)) {
        return 0
    }
    return (Get-Item $Path).Length
}

function Get-DirectorySizeOrZero {
    param([string]$Path)
    if ([string]::IsNullOrWhiteSpace($Path) -or -not (Test-Path $Path)) {
        return 0
    }
    $total = 0L
    Get-ChildItem -Path $Path -Recurse -File | ForEach-Object { $total += $_.Length }
    return $total
}

function Get-Median {
    param([double[]]$Values)
    $sorted = @($Values | Sort-Object)
    if ($sorted.Count -eq 0) {
        return 0.0
    }
    $middle = [int]($sorted.Count / 2)
    if (($sorted.Count % 2) -eq 1) {
        return [double]$sorted[$middle]
    }
    return ([double]$sorted[$middle - 1] + [double]$sorted[$middle]) / 2.0
}

function New-SamFixture {
    param(
        [string]$Path,
        [int]$RequestedRecords
    )

    $builder = [System.Text.StringBuilder]::new()
    [void]$builder.Append("@HD`tVN:1.6`tSO:unknown`n")
    [void]$builder.Append("@RG`tID:rg1`tLB:lib1`tSM:sample1`n")
    [void]$builder.Append("@SQ`tSN:chr1`tLN:2000000`n")

    $pairTemplates = [Math]::Max(1, [int]($RequestedRecords / 4))
    $singleRecords = [Math]::Max(1, $RequestedRecords - ($pairTemplates * 2))

    for ($i = 0; $i -lt $pairTemplates; ++$i) {
        $position1 = 1 + (($i % 500) * 10)
        $position2 = $position1 + 100
        $quality = if (($i % 2) -eq 0) { "IIIII" } else { "!!!!!" }
        [void]$builder.Append("pair$i`t65`tchr1`t$position1`t60`t5M`t=`t$position2`t105`tAAAAA`t$quality`tRG:Z:rg1`n")
        [void]$builder.Append("pair$i`t129`tchr1`t$position2`t60`t5M`t=`t$position1`t-105`tCCCCC`t$quality`tRG:Z:rg1`n")
    }

    for ($i = 0; $i -lt $singleRecords; ++$i) {
        $position = 1 + (($i % 500) * 10)
        $quality = if (($i % 2) -eq 0) { "IIIII" } else { "!!!!!" }
        [void]$builder.Append("single$i`t0`tchr1`t$position`t50`t5M`t*`t0`t0`tGGGGG`t$quality`tRG:Z:rg1`n")
    }

    Write-AsciiFile -Path $Path -Content $builder.ToString()
}

function New-BedFixture {
    param(
        [string]$Path,
        [int]$RequestedIntervals
    )

    $builder = [System.Text.StringBuilder]::new()
    for ($i = 0; $i -lt $RequestedIntervals; ++$i) {
        $start = ($i % 10000) * 10
        $end = $start + 5
        [void]$builder.Append("chr1`t$start`t$end`tinterval$i`n")
    }
    Write-AsciiFile -Path $Path -Content $builder.ToString()
}

function Invoke-BenchmarkOnce {
    param([string[]]$Arguments)

    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    & $Burnham @Arguments *> $null
    $exitCode = $LASTEXITCODE
    $timer.Stop()

    if ($exitCode -ne 0) {
        throw "burnham $($Arguments -join ' ') failed with exit code $exitCode"
    }

    return [Math]::Round($timer.Elapsed.TotalMilliseconds, 3)
}

function Invoke-BenchmarkCommand {
    param(
        [string]$Name,
        [string[]]$Arguments,
        [string]$PrimaryInputPath = "",
        [string]$PrimaryOutputPath = ""
    )

    for ($i = 0; $i -lt $WarmupCount; ++$i) {
        [void](Invoke-BenchmarkOnce -Arguments $Arguments)
    }

    $samples = @()
    for ($i = 0; $i -lt $RepeatCount; ++$i) {
        $samples += Invoke-BenchmarkOnce -Arguments $Arguments
    }

    $min = [double]($samples | Measure-Object -Minimum).Minimum
    $max = [double]($samples | Measure-Object -Maximum).Maximum
    $mean = [double]($samples | Measure-Object -Average).Average
    $median = Get-Median -Values $samples
    $throughput = if ($median -gt 0.0) { [Math]::Round($RecordCount / ($median / 1000.0), 3) } else { 0.0 }

    return [ordered]@{
        name = $Name
        command = "burnham " + ($Arguments -join " ")
        repeat_count = $RepeatCount
        warmup_count = $WarmupCount
        samples_milliseconds = @($samples | ForEach-Object { [Math]::Round([double]$_, 3) })
        min_milliseconds = [Math]::Round($min, 3)
        median_milliseconds = [Math]::Round($median, 3)
        mean_milliseconds = [Math]::Round($mean, 3)
        max_milliseconds = [Math]::Round($max, 3)
        records_per_second = $throughput
        input_bytes = Get-FileSizeOrZero -Path $PrimaryInputPath
        output_bytes = Get-FileSizeOrZero -Path $PrimaryOutputPath
        exit_code = 0
    }
}

$SamPath = Join-Path $BenchDir "input.sam"
$BedPath = Join-Path $BenchDir "input.bed"
$ChainPath = Join-Path $BenchDir "simple.chain"
$IndexPath = Join-Path $BenchDir "simple.bci"

New-SamFixture -Path $SamPath -RequestedRecords $RecordCount
New-BedFixture -Path $BedPath -RequestedIntervals $RecordCount
Write-AsciiFile -Path $ChainPath -Content "chain 1 chrLift 3000000 + 1000 2001000 chr1 2000000 + 0 2000000 1`n2000000`n"

$sortedCoordinate = Join-Path $BenchDir "sorted-coordinate.sam"
$sortedQueryname = Join-Path $BenchDir "sorted-queryname.sam"
$alignmentSummary = Join-Path $BenchDir "alignment-summary.txt"
$markedSam = Join-Path $BenchDir "marked.sam"
$liftedBed = Join-Path $BenchDir "lifted.bed"
$liftedSam = Join-Path $BenchDir "lifted.sam"

$results = @()
$results += Invoke-BenchmarkCommand -Name "validate-sam-summary" -Arguments @("validate-sam", "--input", $SamPath, "--summary") -PrimaryInputPath $SamPath
$results += Invoke-BenchmarkCommand -Name "sort-coordinate" -Arguments @("sort", "--input", $SamPath, "--output", $sortedCoordinate, "--order", "coordinate", "--force") -PrimaryInputPath $SamPath -PrimaryOutputPath $sortedCoordinate
$results += Invoke-BenchmarkCommand -Name "sort-queryname" -Arguments @("sort", "--input", $SamPath, "--output", $sortedQueryname, "--order", "queryname", "--force") -PrimaryInputPath $SamPath -PrimaryOutputPath $sortedQueryname
$results += Invoke-BenchmarkCommand -Name "alignment-summary" -Arguments @("alignment-summary", "--input", $SamPath, "--output", $alignmentSummary, "--force") -PrimaryInputPath $SamPath -PrimaryOutputPath $alignmentSummary
$results += Invoke-BenchmarkCommand -Name "mark-dup" -Arguments @("mark-dup", "--input", $SamPath, "--output", $markedSam, "--metrics-file", (Join-Path $BenchDir "mark-dup.metrics.json"), "--duplicate-type-tags", "--force") -PrimaryInputPath $SamPath -PrimaryOutputPath $markedSam
$results += Invoke-BenchmarkCommand -Name "index-chain" -Arguments @("index-chain", "--chain", $ChainPath, "--output", $IndexPath, "--source-side", "query", "--force") -PrimaryInputPath $ChainPath -PrimaryOutputPath $IndexPath
$results += Invoke-BenchmarkCommand -Name "lift-bed" -Arguments @("lift-bed", "--index", $IndexPath, "--input", $BedPath, "--output", $liftedBed, "--metrics-file", (Join-Path $BenchDir "lift-bed.metrics.json"), "--force") -PrimaryInputPath $BedPath -PrimaryOutputPath $liftedBed
$results += Invoke-BenchmarkCommand -Name "lift-sam-chain" -Arguments @("lift-sam", "--index", $IndexPath, "--input", $SamPath, "--output", $liftedSam, "--metrics-file", (Join-Path $BenchDir "lift-sam.metrics.json"), "--force") -PrimaryInputPath $SamPath -PrimaryOutputPath $liftedSam

$summary = [ordered]@{
    schema_version = 2
    phase = "phase4-text-foundation"
    generated_at_utc = (Get-Date).ToUniversalTime().ToString("o")
    record_count = $RecordCount
    repeat_count = $RepeatCount
    warmup_count = $WarmupCount
    benchmark_dir = $BenchDir
    benchmark_dir_bytes = Get-DirectorySizeOrZero -Path $BenchDir
    note = "This benchmark measures the current text-format Burnham foundation only; htslib/BAM/CRAM performance must be measured after binary IO integration."
    commands = $results
}

$json = $summary | ConvertTo-Json -Depth 8
Write-AsciiFile -Path $OutputPath -Content ($json + "`n")

Write-Host "Wrote Phase 4 benchmark results to $OutputPath"
Write-Host "Benchmark fixture directory: $BenchDir"
Write-Host "Recorded $($results.Count) commands with $RepeatCount measured repeat(s) each"