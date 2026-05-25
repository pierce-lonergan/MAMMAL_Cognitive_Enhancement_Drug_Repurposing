"""Quick inspection of the top hits at the 3 failed positive-control targets.

Helps decide whether the sanity-gate failure is genuine signal loss vs a
"ChEMBL legitimately found stronger binders than the canonical drugs" case.
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
df = pd.read_parquet(ROOT / "data/results/dti_scores.parquet")

print(f"Total: {len(df)} pairs, {df['target_uniprot'].nunique()} targets, "
      f"{df['compound_name'].nunique()} compounds")
print(f"NaN predicted_pkd: {df['predicted_pkd'].isna().sum()}")
print()

FAILED = [
    ("Q01959", "SLC6A3 (DAT) — expected: methylphenidate, modafinil, d-amphetamine"),
    ("P23975", "SLC6A2 (NET) — expected: atomoxetine, methylphenidate"),
    ("P36544", "CHRNA7 — expected: galantamine, encenicline (PAMs)"),
]

for uniprot, label in FAILED:
    print(f"## {label} [{uniprot}]")
    top = df[df["target_uniprot"] == uniprot].nlargest(15, "predicted_pkd")[
        ["compound_name", "predicted_pkd"]
    ]
    print(top.to_string(index=False))
    print()

# Check pKd distribution width per target
print("## pKd std per target (low std = compressed; signal hard to extract)")
stats = df.groupby("target_uniprot")["predicted_pkd"].agg(["min", "median", "max", "std"])
stats["range"] = stats["max"] - stats["min"]
stats = stats.sort_values("range", ascending=False)
print(stats.to_string())
