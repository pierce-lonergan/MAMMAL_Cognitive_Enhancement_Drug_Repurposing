# Boltz-2 Allosteric Audit — CHRNA7 Rescue Gate

**Allosteric targets passing gate**: 2/2

Re-running the v1 allosteric audit with Boltz-2 `affinity_probability_binary` as the rank-scoring metric instead of MAMMAL pKd. Gate: ≥2 allosteric ligands in top 25% per target.

## CHRNA7 (P36544)

**Allosteric rescue gate**: ✅ (2/3 ligands in top 25%)

| Allosteric ligand | rank | percentile | Boltz score |
|---|---|---|---|
| tc-5619 | 5 | 100% ✅ | 0.766 |
| encenicline | 4 | 80% ✅ | 0.683 |
| galantamine | 2 | 40% ❌ | 0.498 |

## PDE4D (Q08499)

**Allosteric rescue gate**: ✅ (1/1 ligands in top 25%)

| Allosteric ligand | rank | percentile | Boltz score |
|---|---|---|---|
| bpn14770 | 3 | 100% ✅ | 0.963 |

Orthosteric reference:

| Orthosteric ligand | rank | percentile | Boltz score |
|---|---|---|---|
| rolipram | 2 | 67% | 0.907 |

## Positive-control retention

- pitolisant @ Q9Y5N1: rank 2, pct 100% ✅

**Positive-control retention PASS**: all checked controls in top 10 at their cognate target.
