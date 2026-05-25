#!/bin/bash
# Quick status check for the WSL2 Boltz overnight sweep.
# Usage from Windows: wsl -d Ubuntu -u root -- bash "/mnt/c/.../scripts/_wsl2_sweep_status.sh"
# Or from inside Ubuntu: bash _wsl2_sweep_status.sh

PROJ="/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing"
LOG="/tmp/wsl2_boltz_sweep.log"

echo "=== Process status ==="
if pgrep -f "_wsl2_boltz_full_sweep.py" >/dev/null; then
    echo "RUNNING (PID: $(pgrep -f _wsl2_boltz_full_sweep.py))"
else
    echo "NOT RUNNING — sweep may have finished or crashed"
fi

echo ""
echo "=== GPU activity ==="
nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.free --format=csv 2>&1 | head -2

echo ""
echo "=== Sweep log tail (last 8 lines) ==="
if [ -f "$LOG" ]; then
    tail -8 "$LOG"
else
    echo "no log at $LOG"
fi

echo ""
echo "=== Parquet row count ==="
source /root/mammal_env/bin/activate 2>/dev/null
python3 -c "
import pandas as pd
from pathlib import Path
p = Path('$PROJ/data/results/v2/boltzina_affinity.parquet')
if p.exists():
    df = pd.read_parquet(p)
    print(f'Rows: {len(df)}')
    print(f'Targets: {df[\"target_uniprot\"].nunique()}')
    print(f'Compounds: {df[\"compound_name\"].nunique()}')
    print(f'NaN: {df[\"affinity_pred_value\"].isna().sum()}')
    print(f'Most recent scored_at: {df[\"scored_at\"].max()}')
else:
    print('no parquet yet')
"
