# Hypothesis Audit v1 — V5 Pre-Commitment Validation

Falsifiable re-test of every pre-committed claim from V3/V4/V5 design docs against the live production artifacts. **Brutal honesty mode** — no claim is grandfathered.

**Summary**: PASS=19 | DEGRADE=3 | FAIL=0 | INSUFFICIENT_DATA=0 (of 22 hypotheses)

## Verdicts

| ID | Claim | Status | Measured | Expected |
|---|---|---|---|---|
| H1 | Tanimoto ρ beats MAMMAL ρ at every audited cognition target | **PASS** | 7 wins, 0 losses | ≥7 wins, 0 losses |
| H2 | SLC6A3 (DAT) isotonic post-cal ρ ∈ [+0.45, +0.65] (fit time) | **PASS** | ρ=+0.619 | [+0.45, +0.65] |
| H3 | SLC6A2 (NET) isotonic post-cal ρ ∈ [+0.30, +0.55] (fit time) | **PASS** | ρ=+0.396 | [+0.30, +0.55] |
| H4 | ≥5 of 7 positive controls in top-20% at expected target (via calibrated_pkd) | **DEGRADE** | 3 of 7 targets PASS | ≥5 / 7 |
| H5 | Negative controls (peripheral-only) average BELOW library median rrf_score | **DEGRADE** | neg_mean=0.476; lib_median=0.062; n_neg_in_lib=11 | < 0.062 |
| H6 | Top-25 PASS-only set spans ≥5 distinct mechanism classes | **PASS** | 15 classes: ['anti_epileptic', 'cardiac_beta_blocker', 'chem | ≥5 |
| H7 | §8.0b-zn assigns expected status + tier-1 hits for 8 reference compounds | **PASS** | 8 of 8 correct | ≥7 of 8 |
| H8 | §7.5 pocket DB validation: 13/13 gates pass | **PASS** | 13/13 | 13/13 |
| H9 | §8.15 disagreement signal surfaces ≥1 novel_scaffold AND ≥1 activity_cliff | **PASS** | novel_scaffold=2067, activity_cliff=2113 | ≥1 each |
| H10 | CHRNA7 rescue: TC-5619 + encenicline in top-25 by Boltzina affinity | **PASS** | found={'tc-5619': 1, 'encenicline': 5}; in_top25=2/2 | 2/2 in top-25 |
| H11 | SLC6A3 (Tier A) calibrator drift |Δρ| ≤ 0.20 between fit and audit | **DEGRADE** | Δρ = -0.190 | |Δρ| ≤ 0.10 (PASS) / ≤ 0.20 (DEGRADE) /  |
| H12 | Adding §8.7 MoA ranker preserves v6 top-3 in v7 | **PASS** | v6=['d-amphetamine', 'methylphenidate', 'bupropion']; v7=['d | set equality |
| H13 | BALM adapter availability probe is well-formed | **PASS** | {'available': False, 'reason': 'BALM repo not found for subp | dict with 'available' boolean |
| H14 | OT Genetics fetcher includes ≥4 canonical cognition GWAS | **PASS** | 5 studies registered | ≥4 cognition GWAS |
| H15 | Cluster D PyMC stub path produces σ(θ) ∈ (0, 1) per target | **PASS** | w_pipeline={'P22303': 0.7721691574281991, 'Q01959': 0.372124 | all w ∈ (0, 1) |
| H16 | TxGNN per-disease API: 5 anchors + availability + weighted mean | **PASS** | 5 anchors, avail=False, wmean=0.70 | None |
| H17 | V7.1 PBPK PET anchors produce finite bounded residuals | **PASS** | residuals=[-0.191, -0.54, -0.65] | None |
| H18 | V7.2 PRISMA: 12 classes, all peak_g ≤ Roberts ceiling 0.50 | **PASS** | n_classes=12, n_violations=0 | None |
| H19 | V8.1 LINCS: 3 neural cell lines weighted 1.0; probe well-formed | **PASS** | NPC/NEU/SHSY5Y all weight=1.0 | None |
| H20 | V8.1b JUMP-CP: 13 sources × 3 embedding types registered | **PASS** | sources=13, embeddings=3 | None |
| H21 | V6.B.4: 15 ref compounds + 4 gates + Gate 1 fires on g > 0.50 | **PASS** | n_compounds=15, gate1='FAIL' | None |
| H22 | V7.3: Cluster D multiplicative gate monotonic in relevance_post | **PASS** | g(high)=0.135, g(low)=0.015 | None |

## Detail per hypothesis

### H1 — PASS

**Claim**: Tanimoto ρ beats MAMMAL ρ at every audited cognition target

**Expected**: ≥7 wins, 0 losses
**Measured**: 7 wins, 0 losses

_Note_: Parsed from reports/tanimoto_baseline_v1.md table

### H2 — PASS

**Claim**: SLC6A3 (DAT) isotonic post-cal ρ ∈ [+0.45, +0.65] (fit time)

**Expected**: [+0.45, +0.65]
**Measured**: ρ=+0.619

_Note_: From router_decisions.csv at fit time (NOT audit). H11 tests calibrator drift since fit.

### H3 — PASS

**Claim**: SLC6A2 (NET) isotonic post-cal ρ ∈ [+0.30, +0.55] (fit time)

**Expected**: [+0.30, +0.55]
**Measured**: ρ=+0.396

### H4 — DEGRADE

**Claim**: ≥5 of 7 positive controls in top-20% at expected target (via calibrated_pkd)

**Expected**: ≥5 / 7
**Measured**: 3 of 7 targets PASS

_Note_: DTI source: dti_scores_calibrated.parquet

_Raw_:
```json
{
  "passes": {
    "P22303": [
      "donepezil"
    ],
    "Q01959": [
      "methylphenidate",
      "d-amphetamine"
    ],
    "P21728": [
      "aripiprazole"
    ]
  },
  "fails": {
    "P22303": [
      "rivastigmine(94.0%)",
      "galantamine(49.0%)"
    ],
    "Q01959": [
      "modafinil(46.3%)"
    ],
    "P23975": [
      "atomoxetine(44.6%)",
      "methylphenidate(44.0%)"
    ],
    "Q9Y5N1": [
      "pitolisant(88.6%)"
    ],
    "P36544": [
      "encenicline(71.8%)",
      "galantamine(99.7%)",
      "tc-5619(72.8%)"
    ],
    "Q08499": [
      "rolipram(83.2%)",
      "bpn14770(99.0%)"
    ]
  }
}
```

### H5 — DEGRADE

**Claim**: Negative controls (peripheral-only) average BELOW library median rrf_score

**Expected**: < 0.062
**Measured**: neg_mean=0.476; lib_median=0.062; n_neg_in_lib=11

_Note_: From final_ranking_v7_moa.parquet

_Raw_:
```json
{
  "neg_compound_scores": {
    "atenolol": 1.013,
    "enalapril": 1.011,
    "ibuprofen": 1.007,
    "metformin": 0.939,
    "ranitidine": 0.899,
    "naproxen": 0.101,
    "cetirizine": 0.07,
    "loratadine": 0.064,
    "warfarin": 0.059,
    "simvastatin": 0.038,
    "omeprazole": 0.037
  }
}
```

### H6 — PASS

**Claim**: Top-25 PASS-only set spans ≥5 distinct mechanism classes

**Expected**: ≥5
**Measured**: 15 classes: ['anti_epileptic', 'cardiac_beta_blocker', 'chembl_binder', 'cholinergic', 'dopaminergic', 'glutamatergic_ampa', 'glutamatergic_nmda', 'isr_modulator', 'maoi', 'noradrenergic', 'nri_ndri', 'nsaid', 'parkinsons_dopamine', 'phosphodiesterase', 'sigma']

### H7 — PASS

**Claim**: §8.0b-zn assigns expected status + tier-1 hits for 8 reference compounds

**Expected**: ≥7 of 8
**Measured**: 8 of 8 correct

_Raw_:
```json
{
  "misclassified": []
}
```

### H8 — PASS

**Claim**: §7.5 pocket DB validation: 13/13 gates pass

**Expected**: 13/13
**Measured**: 13/13

### H9 — PASS

**Claim**: §8.15 disagreement signal surfaces ≥1 novel_scaffold AND ≥1 activity_cliff

**Expected**: ≥1 each
**Measured**: novel_scaffold=2067, activity_cliff=2113

### H10 — PASS

**Claim**: CHRNA7 rescue: TC-5619 + encenicline in top-25 by Boltzina affinity

**Expected**: 2/2 in top-25
**Measured**: found={'tc-5619': 1, 'encenicline': 5}; in_top25=2/2

### H11 — DEGRADE

**Claim**: SLC6A3 (Tier A) calibrator drift |Δρ| ≤ 0.20 between fit and audit

**Expected**: |Δρ| ≤ 0.10 (PASS) / ≤ 0.20 (DEGRADE) / > 0.20 (FAIL)
**Measured**: Δρ = -0.190

_Note_: Audit n=10 vs fit n=23

### H12 — PASS

**Claim**: Adding §8.7 MoA ranker preserves v6 top-3 in v7

**Expected**: set equality
**Measured**: v6=['d-amphetamine', 'methylphenidate', 'bupropion']; v7=['d-amphetamine', 'methylphenidate', 'bupropion']

### H13 — PASS

**Claim**: BALM adapter availability probe is well-formed

**Expected**: dict with 'available' boolean
**Measured**: {'available': False, 'reason': 'BALM repo not found for subprocess fallback. Set

### H14 — PASS

**Claim**: OT Genetics fetcher includes ≥4 canonical cognition GWAS

**Expected**: ≥4 cognition GWAS
**Measured**: 5 studies registered

### H15 — PASS

**Claim**: Cluster D PyMC stub path produces σ(θ) ∈ (0, 1) per target

**Expected**: all w ∈ (0, 1)
**Measured**: w_pipeline={'P22303': 0.7721691574281991, 'Q01959': 0.3721242776511587, 'P36544': 0.33236988665395295}

### H16 — PASS

**Claim**: TxGNN per-disease API: 5 anchors + availability + weighted mean

**Measured**: 5 anchors, avail=False, wmean=0.70

### H17 — PASS

**Claim**: V7.1 PBPK PET anchors produce finite bounded residuals

**Measured**: residuals=[-0.191, -0.54, -0.65]

### H18 — PASS

**Claim**: V7.2 PRISMA: 12 classes, all peak_g ≤ Roberts ceiling 0.50

**Measured**: n_classes=12, n_violations=0

### H19 — PASS

**Claim**: V8.1 LINCS: 3 neural cell lines weighted 1.0; probe well-formed

**Measured**: NPC/NEU/SHSY5Y all weight=1.0

### H20 — PASS

**Claim**: V8.1b JUMP-CP: 13 sources × 3 embedding types registered

**Measured**: sources=13, embeddings=3

### H21 — PASS

**Claim**: V6.B.4: 15 ref compounds + 4 gates + Gate 1 fires on g > 0.50

**Measured**: n_compounds=15, gate1='FAIL'

### H22 — PASS

**Claim**: V7.3: Cluster D multiplicative gate monotonic in relevance_post

**Measured**: g(high)=0.135, g(low)=0.015

---

Generated by `scripts/41_v5_hypothesis_audit.py`. JSON ledger at `data\results\v2\hypothesis_audit_v1.json`.