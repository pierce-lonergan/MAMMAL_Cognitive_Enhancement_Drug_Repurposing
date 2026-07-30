# B1 — Does the drug x EXPERIENCE reframe survive its own evidence?

Pre-registered test of the hypothesis that durable cognitive gain is a drug x experience interaction rather than a drug property. Reproduced by `scripts/124_paired_experience_contrast.py` from `data/raw/paired_experience_ledger.csv`.

**Pre-registered criteria (fixed before computing):** SUCCESS = post-washout effect higher in PAIRED than UNPAIRED studies, significant by n-weighted permutation, **and** surviving the removal of every study with n < 25. KILL = not significant, or carried entirely by n < 25 studies.

## VERDICT: **KILL**

| test | n-weighted contrast (paired − unpaired) | permutation p | studies |
|---|---|---|---|
| full set | -0.071 | 0.2692 | 15 |
| n >= 25 only | -0.076 | 0.3345 | 6 |

## Why the second criterion decides it

- Largest n among **any** positive-direction study: **21**.
- Paired studies with n >= 100: **5**, of which positive-direction: **0**.

This is the small-study signature, not a mechanism. The paired literature is bimodal: a handful of small positives (n = 8, 9, 21) and a set of very large nulls (n = 593, 1047, 5907). The unpaired literature is uniformly null at every size. A hypothesis that only holds below n = 25 and reverses above n = 100 is indistinguishable from publication bias plus regression to the mean.

## The assembled evidence

| study | compound | paired | population | off-drug | direction | effect | metric | n | verification |
|---|---|---|---|---|---|---|---|---|---|
| focus_affinity_effects | fluoxetine | Y | stroke_patients | N | no_effect | 0.03 | SMD_cochrane_pooled_motor | 5907 | UNVERIFIED_IN_PAYLOAD |
| dcycloserine_ipd_meta | d_cycloserine | Y | anxiety_patients | Y | no_effect_at_followup | 0.19 | hedges_g_at_followup_CI_crosses_ze | 1047 | UNVERIFIED_IN_PAYLOAD |
| dars_dopamine_rehab | levodopa | Y | stroke_patients | N | negative | 0.78 | odds_ratio_favouring_placebo | 593 | UNVERIFIED_IN_PAYLOAD |
| levodopa_patching_amblyopia | levodopa | Y | amblyopia_patients | Y | no_effect |  | nan | 139 | UNVERIFIED_IN_PAYLOAD |
| vortioxetine_cognitive_training | vortioxetine | Y | depression_patients | N | no_effect_at_endpoint | 0.21 | hedges_g_decayed_from_0.57 | 100 | UNVERIFIED_IN_PAYLOAD |
| walker_batson_dexamphetamine | dexamphetamine | Y | stroke_patients | N | positive_fragile |  | nan | 21 | UNVERIFIED_IN_PAYLOAD |
| chamoun_2017 | donepezil | Y | healthy_adults | Y | positive_preliminary |  | within_arm_only_no_between_group_c | 9 | CONFIRMED |
| rokem_silver_2013 | donepezil | Y | healthy_adults | Y | positive | 0.78 | between_group_d_derived_from_repor | 8 | CONFIRMED |
| shellshear_2015 | levodopa | Y | healthy_adults | Y | positive_conditional |  | effect_size_unobtainable_paywall | ? | CONFIRMED_EFFECT_UNOBTAINABLE |
| ssri_amblyopia_meta | ssri_various | Y | amblyopia_patients | N | positive_subclinical | 0.09 | logMAR_improvement | ? | UNVERIFIED_IN_PAYLOAD |
| gilleen_2014 | modafinil | Y | healthy_adults | Y | no_effect_level | 0.84 | AUC_learning_RATE_not_level | ? | UNVERIFIED_IN_PAYLOAD |
| rucker_psilocybin_2022 | psilocybin | N | healthy_adults | Y | no_effect | 0 | null_on_designed_endpoint | 89 | CONFIRMED |
| esketamine_washout | esketamine | N | patients | Y | no_effect |  | nan | ? | UNVERIFIED_IN_PAYLOAD |
| ache_inhibitor_washout | donepezil_galantamine | N | patients | Y | no_effect |  | nan | ? | UNVERIFIED_IN_PAYLOAD |
| stimulant_discontinuation | methylphenidate_guanfacine | N | patients | Y | no_effect |  | nan | ? | UNVERIFIED_IN_PAYLOAD |

## Methodological honesty

1. **Direction, not pooled effect size.** The assembled effects are in incommensurable units (SMD, logMAR, odds ratios, learning-curve AUC). Pooling them numerically would be a category error, so the test scores DIRECTION, which is both commensurable and what the hypothesis actually predicts. Effect values are carried in the table for inspection, not summed.
2. **Unknown n gets weight 1**, the minimum — it refuses to let an unsized study dominate.
3. **Verification status is on every row.** Rows marked `UNVERIFIED_IN_PAYLOAD` had their verification block truncated in the research payload. The verdict does not depend on their exact values, only on their direction, which is triangulated across three mechanistically independent programmes (serotonergic, dopaminergic, glutamatergic).
4. **Patient-population rows are retained but flagged.** They test whether drug x experience produces durable gain *anywhere*, which is the weaker and more favourable version of the hypothesis. Even that version fails at scale.

## What this licenses, and what it forbids

The kill criterion **fired, as pre-registered and as predicted**. Consequences, stated before the fact and honoured now:

- **No further predictor work on compound-level durability is licensed.** The hypothesis that would have justified it has no support at the level of evidence available.
- **The deliverable is the negative result** plus the assay-indexed window screen (B2).
- **The one legitimate positive framing** is narrower than the programme was built around: drug x training may durably improve a *single trained skill* (Rokem & Silver, n = 8, no transfer), and even there the placebo arm also retained its learning.
- **What would overturn this:** a pre-registered, between-group, placebo-controlled drug + identical-training trial in healthy adults, n >= 30/arm (see B4 for the exact N*), with the retained off-drug LEVEL as primary endpoint — and the drug arm winning. Plus transfer off the trained task, which nothing in the current literature demonstrates.

---

Generated by `scripts/124_paired_experience_contrast.py`.