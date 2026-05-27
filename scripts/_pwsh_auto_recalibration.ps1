#requires -Version 5
# §8.12 — Auto-recalibration scheduled-task driver.
#
# Reads configs/auto_recalibration.yaml; if any trigger_files mtime is newer
# than the .last_run marker, runs the cascade. Each cascade step is logged
# to data/calibration/cascade/<step_id>_<timestamp>.log.
#
# Register as a Windows Scheduled Task:
#   schtasks /Create /TN "MAMMAL recalibration" /TR `
#     "powershell.exe -ExecutionPolicy Bypass -File C:\...\scripts\_pwsh_auto_recalibration.ps1" `
#     /SC DAILY /ST 02:00 /F

[CmdletBinding()]
param(
    [string]$ConfigPath = (Join-Path (Split-Path $PSScriptRoot -Parent) "configs\auto_recalibration.yaml"),
    [switch]$Force,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path $PSScriptRoot -Parent
Push-Location $RepoRoot

function Write-Log {
    param([string]$msg, [string]$level = "INFO")
    $ts = (Get-Date).ToString("yyyy-MM-dd HH:mm:ss")
    Write-Host "[$ts][$level] $msg"
}

# --- 1. Locate Python env (mammal_env conda env) ----------------------------
$EnvPython = Join-Path $env:USERPROFILE ".conda\envs\mammal_env\python.exe"
if (-not (Test-Path $EnvPython)) {
    $EnvPython = (Get-Command python).Source
    Write-Log "mammal_env not found; falling back to system python: $EnvPython" "WARN"
}

# --- 2. Parse YAML (PowerShell 5 has no native YAML; use python) ------------
$cfgJson = & $EnvPython -c @"
import json, yaml, sys
with open(r'$ConfigPath','r',encoding='utf-8') as f:
    cfg = yaml.safe_load(f)
print(json.dumps(cfg))
"@
$cfg = $cfgJson | ConvertFrom-Json

# --- 3. Check trigger files vs last_run marker ------------------------------
$markerPath = Join-Path $RepoRoot $cfg.last_run_marker
$lastRun = $null
if (Test-Path $markerPath) {
    $lastRun = (Get-Item $markerPath).LastWriteTime
}

$triggered = $false
foreach ($tf in $cfg.trigger_files) {
    $p = Join-Path $RepoRoot $tf
    if (Test-Path $p) {
        $mt = (Get-Item $p).LastWriteTime
        if ($null -eq $lastRun -or $mt -gt $lastRun) {
            Write-Log "Trigger file newer: $tf (modified $mt)"
            $triggered = $true
        }
    }
}

if (-not $triggered -and -not $Force) {
    Write-Log "No trigger files newer than last run; exiting." "INFO"
    Pop-Location
    return
}

# --- 4. Run §8.16 calibrator QC FIRST as the gating step ---------------------
Write-Log "Running §8.16 calibrator QC gate"
if (-not $DryRun) {
    & $EnvPython "scripts/38_v5_calibrator_qc.py" | Out-Null
}

# Parse QC JSON for Tier-A REFIT_NEEDED — block if found
$qcDir = Join-Path $RepoRoot "data\calibration\qc"
$blocked = $false
if ($cfg.gating_policy.block_on_refit_needed) {
    foreach ($tier_a in $cfg.gating_policy.tier_a_targets) {
        $qcFile = Join-Path $qcDir "$tier_a.json"
        if (Test-Path $qcFile) {
            $qc = Get-Content $qcFile -Raw | ConvertFrom-Json
            if ($qc.status -eq "REFIT_NEEDED") {
                Write-Log "Tier-A target $tier_a needs REFIT (Delta-rho $($qc.delta_rho)) -- BLOCKED" "ERROR"
                $blocked = $true
            }
        }
    }
}

if ($blocked) {
    Write-Log "Cascade BLOCKED by §8.16 gate. Refit Tier-A calibrators or pass -Force." "ERROR"
    Pop-Location
    exit 1
}

# --- 5. Run cascade steps in order -------------------------------------------
$cascadeDir = Join-Path $RepoRoot "data\calibration\cascade"
if (-not (Test-Path $cascadeDir)) { New-Item -ItemType Directory -Path $cascadeDir -Force | Out-Null }
$ts = (Get-Date).ToString("yyyyMMdd_HHmmss")

foreach ($step in $cfg.cascade) {
    Write-Log "Running cascade step: $($step.id) — $($step.description)"
    $logFile = Join-Path $cascadeDir "$($step.id)_$ts.log"
    if ($DryRun) {
        Write-Log "  DRY-RUN: would execute: $($step.command)"
        continue
    }
    try {
        Invoke-Expression "$($step.command) *> '$logFile'"
        Write-Log "  $($step.id) OK — logged to $($logFile)"
    } catch {
        Write-Log "  $($step.id) FAILED: $_" "ERROR"
        if ($step.blocking -eq $true) {
            Write-Log "  Blocking step failed — aborting cascade." "ERROR"
            Pop-Location
            exit 2
        }
    }
}

# --- 6. Update marker ---------------------------------------------------------
$marker = @{
    last_run_utc = (Get-Date).ToUniversalTime().ToString("o")
    cascade_steps_run = $cfg.cascade | ForEach-Object { $_.id }
}
$marker | ConvertTo-Json -Depth 5 | Set-Content -Encoding utf8 $markerPath

Write-Log "Cascade complete. Marker updated at $markerPath" "INFO"
Pop-Location
