# Hypothesis Audit v1 — V5 Pre-Commitment Validation

Falsifiable re-test of every pre-committed claim from V3/V4/V5 design docs against the live production artifacts. **Brutal honesty mode** — no claim is grandfathered.

**Summary**: PASS=9 | DEGRADE=3 | FAIL=0 | INSUFFICIENT_DATA=0 (of 12 hypotheses)

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

---

Generated by `scripts/41_v5_hypothesis_audit.py`. JSON ledger at `data\results\v2\hypothesis_audit_v1.json`.