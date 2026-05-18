param(
    [Parameter(Mandatory = $true)]
    [string]$BaselinePath,
    [Parameter(Mandatory = $true)]
    [string]$CurrentPath,
    [double]$RegressionThresholdPercent = 15.0,
    [switch]$FailOnRegression
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path $BaselinePath)) {
    throw "Baseline benchmark JSON not found: $BaselinePath"
}
if (-not (Test-Path $CurrentPath)) {
    throw "Current benchmark JSON not found: $CurrentPath"
}
if ($RegressionThresholdPercent -lt 0.0) {
    throw "RegressionThresholdPercent must be zero or greater"
}

function Get-CommandMilliseconds {
    param($Command)
    if ($null -ne $Command.median_milliseconds) {
        return [double]$Command.median_milliseconds
    }
    if ($null -ne $Command.milliseconds) {
        return [double]$Command.milliseconds
    }
    throw "Command '$($Command.name)' does not contain median_milliseconds or milliseconds"
}

$baseline = Get-Content $BaselinePath -Raw | ConvertFrom-Json
$current = Get-Content $CurrentPath -Raw | ConvertFrom-Json

$baselineByName = @{}
foreach ($command in $baseline.commands) {
    $baselineByName[$command.name] = $command
}

$rows = @()
$regressions = @()
foreach ($command in $current.commands) {
    if (-not $baselineByName.ContainsKey($command.name)) {
        $rows += [pscustomobject]@{
            Name = $command.name
            BaselineMs = $null
            CurrentMs = Get-CommandMilliseconds -Command $command
            DeltaPercent = $null
            Status = "new"
        }
        continue
    }

    $baselineMs = Get-CommandMilliseconds -Command $baselineByName[$command.name]
    $currentMs = Get-CommandMilliseconds -Command $command
    $deltaPercent = if ($baselineMs -gt 0.0) { (($currentMs - $baselineMs) / $baselineMs) * 100.0 } else { 0.0 }
    $status = if ($deltaPercent -gt $RegressionThresholdPercent) { "regression" } elseif ($deltaPercent -lt (-1.0 * $RegressionThresholdPercent)) { "improved" } else { "ok" }

    $row = [pscustomobject]@{
        Name = $command.name
        BaselineMs = [Math]::Round($baselineMs, 3)
        CurrentMs = [Math]::Round($currentMs, 3)
        DeltaPercent = [Math]::Round($deltaPercent, 2)
        Status = $status
    }
    $rows += $row
    if ($status -eq "regression") {
        $regressions += $row
    }
}

$rows | Format-Table -AutoSize

if ($regressions.Count -gt 0) {
    Write-Host "Detected $($regressions.Count) regression(s) above $RegressionThresholdPercent percent"
    if ($FailOnRegression) {
        exit 1
    }
} else {
    Write-Host "No regressions above $RegressionThresholdPercent percent"
}