#!/bin/bash
# WSL2 — install PyG + TxGNN in a SEPARATE venv, isolated from mammal_env.
#
# Why separate venv: PyG's CUDA extensions (pyg-lib, torch-scatter, etc) are
# compiled against specific PyTorch tags. The matching wheel index is at
# https://data.pyg.org/whl/torch-2.7.0+cu128.html. Our mammal_env uses PyTorch
# 2.12.0.dev (required for Blackwell sm_120 cuequivariance kernels in Boltz).
# 2.12.0.dev is too far ahead of any PyG wheel set — pyg-lib's libpyg.so
# segfaults (Bus error) when loaded against 2.12 due to ABI changes in
# torch::Tensor's storage handling.
#
# Solution: txgnn_env venv with PyTorch 2.7.0+cu128 (the version PyG targets).
# TxGNN doesn't need cuequivariance. Run TxGNN inference via:
#     wsl -d Ubuntu -u root -- bash -c "source /root/txgnn_env/bin/activate && python scripts/23_v3_cluster_c.py"

set -euo pipefail
VENV="/root/txgnn_env"

echo "=== Step 1: apt prerequisites (idempotent) ==="
SUDO=""
if [ "$(id -u)" -ne 0 ]; then SUDO="sudo"; fi
$SUDO apt-get update -qq
$SUDO apt-get install -y --no-install-recommends python3 python3-venv python3-dev build-essential

echo ""
echo "=== Step 2: Fresh venv at $VENV ==="
if [ ! -d "$VENV" ]; then
    python3 -m venv "$VENV"
fi
source "$VENV/bin/activate"
python -m pip install --upgrade --quiet pip

echo ""
echo "=== Step 3: PyTorch 2.7.0 nightly cu128 (PyG-compatible) ==="
pip install --quiet --pre torch==2.7.0+cu128 torchvision \
    --index-url https://download.pytorch.org/whl/nightly/cu128 \
    || pip install --quiet --pre torch torchvision --index-url https://download.pytorch.org/whl/nightly/cu128

echo ""
echo "=== Step 4: PyG core wheels (matching torch 2.7.0+cu128) ==="
pip install --quiet \
    pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
    -f https://data.pyg.org/whl/torch-2.7.0+cu128.html

pip install --quiet torch-geometric igraph pandas pyarrow

echo ""
echo "=== Step 5: TxGNN ==="
pip install --quiet txgnn 2>/dev/null \
    || pip install --quiet git+https://github.com/mims-harvard/TxGNN.git

echo ""
echo "=== Step 6: Verify ==="
python -c "
import torch
print(f'PyTorch: {torch.__version__}; CUDA: {torch.cuda.is_available()}')
import torch_geometric
print(f'torch_geometric: {torch_geometric.__version__}')
import igraph
print(f'igraph: {igraph.__version__}')
try:
    from txgnn import TxGNN, TxData
    print('txgnn: import OK')
except Exception as e:
    print(f'txgnn: IMPORT FAILED ({type(e).__name__}: {e})')
"

echo ""
echo "PrimeKG download script: scripts/_wsl2_download_primekg.sh"
echo "Cluster C orchestrator: scripts/23_v3_cluster_c.py (run inside this venv)"
