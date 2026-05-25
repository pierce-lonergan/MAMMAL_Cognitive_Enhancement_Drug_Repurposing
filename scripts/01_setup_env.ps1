# Stage 1 - Environment setup for MAMMAL drug repurposing pipeline.
#
# Creates conda environment 'mammal_env' on Windows 11, installs PyTorch nightly
# with CUDA 12.8 (required for Blackwell sm_120 / RTX 5070), installs MAMMAL and
# this package in editable mode, then runs validation checks.
#
# Usage:  .\scripts\01_setup_env.ps1
# Re-run is safe (idempotent: skips existing env unless -Recreate is passed).

param(
    [switch]$Recreate,
    [switch]$SkipModelDownload,
    [string]$EnvName = "mammal_env",
    [string]$PythonVersion = "3.10"
)

$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path -Parent $PSScriptRoot

Write-Host ""
Write-Host "===========================================================" -ForegroundColor Cyan
Write-Host " MAMMAL drug repurposing - environment setup" -ForegroundColor Cyan
Write-Host "===========================================================" -ForegroundColor Cyan
Write-Host " Env name      : $EnvName"
Write-Host " Python        : $PythonVersion"
Write-Host " Repo root     : $RepoRoot"
Write-Host " Recreate env  : $Recreate"
Write-Host ""

# --- 1. Locate conda --------------------------------------------------------
# Prefer the _conda.exe binary at the install root (bypasses the Scripts\conda.exe
# Windows launcher bug that breaks when the user profile path contains spaces).
$candidatePaths = @(
    "$env:USERPROFILE\miniconda3\_conda.exe",
    "$env:USERPROFILE\anaconda3\_conda.exe",
    "C:\miniconda3\_conda.exe",
    "C:\anaconda3\_conda.exe",
    "C:\ProgramData\miniconda3\_conda.exe",
    "C:\ProgramData\anaconda3\_conda.exe"
)
$conda = $null
foreach ($p in $candidatePaths) {
    if (Test-Path $p) { $conda = $p; break }
}
if (-not $conda) {
    $onPath = Get-Command conda -ErrorAction SilentlyContinue
    if ($onPath) { $conda = $onPath.Source }
}
if (-not $conda) {
    Write-Error "conda not found. Install Miniconda first (https://docs.conda.io/en/latest/miniconda.html) then rerun."
    exit 1
}
Write-Host "[1/7] conda found at: $conda" -ForegroundColor Green

# --- 2. Create / verify the conda env ---------------------------------------
$envExists = (conda env list) -match "^$EnvName\s"
if ($envExists -and $Recreate) {
    Write-Host "[2/7] Removing existing env '$EnvName' (-Recreate passed)..." -ForegroundColor Yellow
    conda env remove -n $EnvName -y
    $envExists = $false
}
if (-not $envExists) {
    Write-Host "[2/7] Creating conda env '$EnvName' with Python $PythonVersion..." -ForegroundColor Green
    conda create -n $EnvName "python=$PythonVersion" -y
    if ($LASTEXITCODE -ne 0) { Write-Error "conda env creation failed"; exit 1 }
} else {
    Write-Host "[2/7] Env '$EnvName' already exists (use -Recreate to rebuild)." -ForegroundColor Green
}

# --- 3. Resolve env python path ---------------------------------------------
# Don't rely on `conda activate` (PowerShell init quirks). Call the env's python directly.
# conda may place envs in either <prefix>/envs/ or ~/.conda/envs/ depending on
# install location and write permissions.
$condaPrefix = Split-Path -Parent (Split-Path -Parent $conda)
$candidateEnvs = @(
    (Join-Path $condaPrefix "envs\$EnvName\python.exe"),
    (Join-Path $env:USERPROFILE ".conda\envs\$EnvName\python.exe")
)
$envPython = $null
foreach ($p in $candidateEnvs) {
    if (Test-Path $p) { $envPython = $p; break }
}
if (-not $envPython) {
    Write-Error "Env python not found. Tried: $($candidateEnvs -join ', ')"
    exit 1
}
Write-Host "[3/7] Env python: $envPython" -ForegroundColor Green

# --- 4. Install PyTorch nightly (Blackwell sm_120 requirement) --------------
Write-Host "[4/7] Installing PyTorch nightly with CUDA 12.8 (Blackwell sm_120)..." -ForegroundColor Green
& $envPython -m pip install --pre `
    --index-url https://download.pytorch.org/whl/nightly/cu128 `
    torch torchvision
if ($LASTEXITCODE -ne 0) { Write-Error "PyTorch nightly install failed"; exit 1 }

# --- 5. Install MAMMAL + this package ---------------------------------------
Write-Host "[5/7] Installing biomed-multi-alignment..." -ForegroundColor Green
# Without [examples] extra — tiledbsoma (cellxgene-census dep) fails wheel build on Windows.
# We don't need cellxgene tooling for DTI/BBBP/ClinTox heads; the example task code
# (mammal/examples/dti_bindingdb_kd/) ships with the base package.
& $envPython -m pip install "biomed-multi-alignment>=0.2.4"
if ($LASTEXITCODE -ne 0) { Write-Error "biomed-multi-alignment install failed"; exit 1 }

Write-Host "[5/7] Installing mammal-repurposing in editable mode..." -ForegroundColor Green
Push-Location $RepoRoot
try {
    & $envPython -m pip install -e ".[dev,notebook]"
    if ($LASTEXITCODE -ne 0) { Write-Error "editable install failed"; exit 1 }
} finally {
    Pop-Location
}

# --- 6. Validate CUDA + GPU --------------------------------------------------
Write-Host "[6/7] Validating CUDA + GPU..." -ForegroundColor Green
$cudaCheck = @'
import torch, sys
print(f"PyTorch         : {torch.__version__}")
print(f"CUDA available  : {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"Device          : {torch.cuda.get_device_name(0)}")
    print(f"Compute cap.    : {torch.cuda.get_device_capability(0)}")
    print(f"CUDA runtime    : {torch.version.cuda}")
else:
    print("WARNING: CUDA not available - inference will run on CPU (slow).", file=sys.stderr)
    sys.exit(2)
'@
& $envPython -c $cudaCheck
if ($LASTEXITCODE -ne 0) {
    Write-Warning "CUDA validation reported issues (exit code $LASTEXITCODE). Inference may fall back to CPU."
}

# --- 7. Trigger MAMMAL DTI head download (optional) -------------------------
if ($SkipModelDownload) {
    Write-Host "[7/7] Skipping model weight download (-SkipModelDownload passed)." -ForegroundColor Yellow
} else {
    Write-Host "[7/7] Downloading MAMMAL DTI head (~1.8 GB, first-run only)..." -ForegroundColor Green
    $modelCheck = @'
from mammal.model import Mammal
m = Mammal.from_pretrained("ibm/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd")
print(f"Loaded MAMMAL DTI model OK. Parameter count: {sum(p.numel() for p in m.parameters()):,}")
'@
    & $envPython -c $modelCheck
    if ($LASTEXITCODE -ne 0) { Write-Error "MAMMAL model load failed"; exit 1 }
}

Write-Host ""
Write-Host "===========================================================" -ForegroundColor Cyan
Write-Host " Setup complete." -ForegroundColor Cyan
Write-Host "===========================================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Next steps:" -ForegroundColor White
Write-Host "  conda activate $EnvName"
Write-Host "  python scripts/02_fetch_targets.py"
Write-Host "  python scripts/03_fetch_compounds.py"
Write-Host "  python scripts/04_score_dti.py"
Write-Host "  python scripts/05_sanity_check.py"
Write-Host ""
