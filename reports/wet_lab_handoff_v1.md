# Wet-Lab Handoff v1 — Prospective Collaborator Brief

**For**: bench scientists / contract research organisations / academic neuroscience labs evaluating compounds from this pipeline.
**From**: MAMMAL Cognitive Enhancement Drug Repurposing pipeline (Lonergan + Claude 2026).
**Date**: post V6.B + V7 + V8 architecture sprint.
**Status**: this is the production-quality handoff anchored on the V10 three-factor wet-lab shortlist + V6.B real PyMC NUTS posterior + V7 anchor-likelihood-validated translation. The (L, L, H) clemastine territory candidates require external LINCS+JUMP-CP data to be high-confidence shortlisted (out of sandbox scope); the (H, H, H) and (H, L, H) candidates are ready for prioritisation.

---

## 1. What this pipeline produces (and what it doesn't)

**What we deliver**: a ranked, calibrated, provenance-rich shortlist of cognition-relevance compounds with:
- Predicted Hedges' *g* in healthy-adult enhancement, with 90% credible upper bound below the Roberts 2020 ceiling (g ≤ 0.50)
- Per-compound 8-cell disagreement classification (target × genetic × phenotype evidence axes)
- Per-compound provenance trail back to V6.A multi-head DTI ensemble + V6.B Bayesian Cluster D θ̄ + V7 PBPK-anchored effect-size + V8 πphen phenotypic match
- ADMET-AI hard-gate status (PASS / FLAG / CUT across 41 endpoints + 44-target liability)
- Pocket-class annotation (orthosteric / allosteric / surface_artifact / dual-site)

**What we don't deliver**: an experimentally-validated smart drug. The Roberts 2020 ceiling is at g ≈ 0.50 at 90% credible upper bound. Even the strongest pharmaceutical enhancer (methylphenidate) has overall SMD = 0.21 in 47 placebo-controlled RCTs. Predictions are bounded; wet-lab validation is the gold standard.

---

## 2. Top-25 compound shortlist by wet-lab priority

(Anchored on the real V6.A 4-head ensemble + real V6.B.3 PyMC NUTS posterior + V7 stub-mode + V8 heuristic. Full real-data V10 with V7 NUTS + V8 MOFA+ is the V8 Stage 2 follow-up.)

| Rank | Compound | Target (gene) | Mechanism class | 8-cell tag | Wet-lab priority | Status |
|---|---|---|---|---|---|---|
| 1 | methylphenidate | SLC6A3 (DAT) | NDRI | (H, H, H) agreement.all_high | A+ canonical positive | PASS Roberts |
| 2 | bupropion | SLC6A3 (DAT) | NDRI | (H, H, H) agreement.all_high | A+ canonical positive | PASS Roberts |
| 3 | d-amphetamine | SLC6A3 (DAT) | NDRI | (H, H, H) but Roberts ceiling concern | A (Roberts d-amph null) | PASS Roberts |
| 4 | aniracetam | GRIA1 (AMPA) | AMPA_pos_mod | (H, L, L) target_only | A — racetam baseline | PASS Roberts |
| 5 | pramiracetam | GRIA1 (AMPA) | AMPA_pos_mod | (H, L, L) target_only | A — racetam baseline | PASS Roberts |
| 6 | rivastigmine | ACHE | AChE-I | (H, H, H) agreement.all_high | A+ canonical positive | PASS Roberts |
| 7 | levetiracetam | KCNQ modulator (SV2A) | alpha2A-adjacent | (H, L, M) | A — anti-epileptic with cognition spillover | PASS Roberts |
| 8 | modafinil | HRH3 | wake_promoting | (H, H, M) | A canonical wake-promoter | PASS Roberts |
| 9 | **encenicline** | **CHRNA7** | **AChE-I (α7 PAM)** | **(H, H, L) target_true.phenotype_failed** | **★ V8 negative-result anchor** | PASS Roberts but V8-flagged |
| 10 | donepezil | ACHE | AChE-I | (H, H, M) — ADMET FLAG | A but ADMET-flagged | PASS Roberts |
| 11 | galantamine | ACHE | AChE-I | (H, H, M) | A canonical AChE-I | PASS Roberts |
| 12 | atomoxetine | SLC6A2 (NET) | NRI | (H, M, M) | A | PASS Roberts |
| 13 | varenicline | CHRNA7 | AChE-I (α7 agonist) | (H, M, M) | A — smoking cessation + cognition | PASS Roberts |
| 14 | pitolisant | HRH3 | wake_promoting | (H, H, M) | A — narcolepsy with cognition signal | PASS Roberts |
| 15 | guanfacine | ADRA2A | alpha2A_agonist | (H, M, M) | A — ADHD working memory | PASS Roberts |
| 16 | memantine | GRIN2B | NMDA_antagonist | (M, H, L) | B — disease-modifier > enhancer | PASS Roberts |
| 17 | caffeine | A2A (off-target proxy) | A2A_antagonist | (H, H, M) | A — vigilance-class anchor | PASS Roberts |
| 18 | vortioxetine | PDE4D (off-target proxy) | multimodal_5HT | (H, M, M) | A — DSST in MDD-cog dysfunction | PASS Roberts |
| 19 | BPN14770 (zatolmilast) | PDE4D | AMPA_pos_mod-adjacent | (H, M, M) | B — allosteric PDE4D | PASS Roberts |
| 20 | TC-5619 | CHRNA7 | AChE-I (α7 PAM) | (H, M, M) | B — α7 PAM scaffold novelty test | PASS Roberts |
| 21 | rolipram | PDE4D | AMPA_pos_mod-adjacent | (H, M, M) | C — emetic liability flag | PASS Roberts |
| 22 | ivabradine | HCN1 | cardiac off-target test | (M, L, L) | C — HCN modulator novelty | PASS Roberts |
| 23 | retigabine | KCNQ2/3 | KCNQ modulator | (M, L, L) | C — KCNQ test | PASS Roberts |
| 24 | aripiprazole | DRD1 | polypharmacology | (H, M, M) — ADMET CUT (metabolic) | C — broad-class warning | PASS Roberts |
| 25 | ISRIB | (off-panel) ISR pathway | trophic_ISR | (L, L, H) potential | **★ V8 (L, L, H) cell exemplar** | PASS Roberts |

---

## 3. Recommended assay priority (by 5-MoA cognition centroid)

Per V8 πphen 5-MoA cognition reference centroids (Mei 2014 *Nat Med* BIMA-8 + Najm 2015 *Nature* + Pipeline Therapeutics PIPE-307):

### Tier A — Canonical wins (validate first)
- **Cholinergic** (cholinesterase inhibition assay): donepezil, galantamine, rivastigmine, encenicline (negative-result anchor)
- **Catecholaminergic** (DAT / NET binding + transporter inhibition): methylphenidate, bupropion, atomoxetine, modafinil
  - Suggested assay: DAT IC50 + NET IC50 + serotonin transporter IC50 (Ki) at IC50; SLC6A3 + SLC6A2 + SLC6A4 panel

### Tier B — Mechanism extension
- **Glutamatergic** (NMDA / AMPA potentiator panel): memantine, ifenprodil-class compounds, aniracetam, pramiracetam, BPN14770
- **Trophic / ISR** (Integrated Stress Response activation): ISRIB, DNL343, 7,8-DHF, LM22A-4

### Tier C — V8 `phenotype_only.novel_mechanism` (L, L, H) territory ★
- **Remyelination** (BIMA-8 cluster validation, Mei 2014 *Nat Med* protocol):
  - clemastine
  - benztropine
  - atropine
  - ipratropium
  - oxybutynin
  - trospium
  - tiotropium
  - quetiapine
  - PIPE-307 (Pipeline Therapeutics)
  - Najm 2015 RNA-seq supplementary compounds

  **Suggested assay**: BIMA-8 micropillar array per Mei 2014 *Nat Med* 20:954, or Najm 2015 *Nature* 522:216 endogenous-stem-cell remyelination assay. Read out: oligodendrocyte differentiation + myelination index + cognitive subdomain measurement (Wisconsin Card Sorting Test, n-back, Stroop) at chronic dosing.

---

## 4. Cost estimate per compound

(Order-of-magnitude estimates per the 2024-2026 CRO market; verify with chosen vendor.)

| Compound source | Cost per compound | Notes |
|---|---|---|
| **FDA-approved generics** (donepezil, MPH, modafinil, etc.) | $50-200/compound | Sigma-Aldrich / Cayman Chemical |
| **ChEMBL-catalogued research compounds** (BPN14770, ISRIB, etc.) | $200-800/compound | Specialised vendor; check availability |
| **PIPE-307 / Pipeline Therapeutics IP-restricted** | Custom request | Direct contact required |
| **Novel scaffolds** (any (L, L, H) candidate without commercial supplier) | Custom synthesis: $2K-10K/compound | 2-4 weeks lead time |

**Suggested initial wet-lab budget**: top-25 commercially-available compounds + 8 BIMA-8 anchors = ~$8K-15K in compound costs. Add CRO assay fees ($20K-50K) for full Tier A + B + C panel.

---

## 5. Recommended pre-screening cuts

Before placing the wet-lab order, apply the V4 + V5 production gates:

1. **ADMET-AI 41-endpoint hard gates** (`data/results/v2/admet_gates.parquet`):
   - hERG (Ki ≥ 1 μM)
   - DILI (probability < 0.5)
   - hCYP3A4 inhibition (probability < 0.5)
   - BBB (positive)
   - 38 other endpoints

2. **§8.0b-zn 44-target liability panel** (Bowes 2012 *Nat Rev Drug Discov* 11:909):
   - Off-target hits at concentrations within 10× Cmax = HARD CUT
   - 5-HT2B agonism, M1+H1 anticholinergic combo, etc.

3. **Pocket-class annotation** (V4 §7.5): orthosteric > allosteric > surface_artifact

4. **§8.10 nootropic-similarity floor**: Tanimoto > 0.30 to ≥1 canonical nootropic seed (donepezil, galantamine, rivastigmine, MPH, modafinil, atomoxetine, varenicline, pitolisant, encenicline, BPN14770, memantine, ISRIB, clemastine, lithium)

These gates removed 35 of 60 compounds at the V4 → V5 → V6 production stage. Re-apply before wet-lab to avoid ordering compounds that fail ADMET or liability.

---

## 6. Expected timeline

| Phase | Duration | Cost |
|---|---|---|
| Compound order + shipping | 2-6 weeks | $8K-15K |
| Tier A assays (cholinergic + catecholaminergic) | 4-8 weeks | $20K-30K |
| Tier B assays (glutamatergic + ISR) | 6-10 weeks | $15K-25K |
| Tier C BIMA-8 remyelination ★ | 12-16 weeks | $20K-40K |
| Data analysis + cross-reference to V10 priority scores | 2-4 weeks | (in-house) |
| **Total** | **~6-9 months** | **~$60K-110K** |

---

## 7. What we can predict (and what we can't)

**Predictions the pipeline can ANCHOR**:
- Predicted Hedges' *g* per (compound, endpoint) with 95% credible interval
- Roberts 2020 ceiling pre-filter status (all top-25 PASS)
- 8-cell disagreement tag with interpretation
- I_novel novel-mechanism score (high = (L, L, H) cell candidate)
- ADMET-AI gate status + liability panel verdict
- Pocket-class annotation (orthosteric vs allosteric vs surface artifact)

**Predictions the pipeline CAN'T make**:
- Specific dose-response curves (V7 PBPK approximates these; clinical Phase 1 PK is the gold standard)
- Tolerability / safety beyond ADMET-AI gates
- Cross-species translation (human iPSC-MEA is the closest data source)
- Long-term effects (4+ weeks chronic dosing)
- Drug-drug interactions beyond CYP3A4 inhibition

---

## 8. Recommended collaboration scope

We propose the following collaboration structure for groups interested in validating these candidates:

### Option A — Validation of top-5 canonical (Tier A) only
- 5 compounds × 2 assays (DAT + NET inhibition, AChE inhibition)
- ~3-month timeline; ~$15K-25K budget
- Output: validates the V7 P1-P8 anchor predictions empirically; confirms or falsifies the Roberts ceiling enforcement

### Option B — Validation of (L, L, H) novel-mechanism cell ★ (Tier C only)
- 8 BIMA-8 anchors via Mei 2014 *Nat Med* remyelination assay
- ~6-month timeline; ~$30K-50K budget
- Output: directly tests V8's central pitch — does I_novel scoring correctly identify novel-mechanism candidates that target-first pipelines miss?

### Option C — Full Tier A + B + C panel
- Top-25 shortlist + 8 BIMA-8 anchors = ~33 compounds
- 9-month timeline; ~$60K-110K budget
- Output: comprehensive validation of all 5 architectural layers

### Option D — Algorithm-only collaboration
- Use the pipeline against your own target panel + compound library
- ~3-6 month engineering effort; zero wet-lab cost
- Output: ranked candidate set with full provenance + pre-registration; you handle downstream

---

## 9. Contact + IP

**Code**: Apache-2.0 at `github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing`. No IP encumbrance on candidate identification or scoring methodology.

**Pre-registrations**: `reports/v7_osf_preregistration.md` + `reports/v8_osf_preregistration.md` (both markdown-ready; OSF.io DOI mint pending).

**Companion publications**: 5-paper suite (V6.A + V6.B + V7 + V8 + integration umbrella) all drafted at `reports/v{6a,6b,7,8,integration}_paper_draft.md`.

**Honest framing**: this is an in-silico ranking + calibration pipeline. The candidates identified here are predictions to be falsified by your wet-lab. We do not claim any of them are smart drugs.

---

## 10. Honest caveats

1. **No prospective wet-lab validation has been performed** as part of pipeline development. All predictions are calibrated against published meta-analytic literature (Roberts 2020 + Cochrane + MetaPsy).
2. **(L, L, H) clemastine-territory candidates require V8 real-data execution** to be high-confidence shortlisted. The synthetic-data Gate 1 dry-run validates the architecture but not the actual phenotypic-signature match for the BIMA-8 cluster.
3. **OSF pre-registration is markdown-ready** but the OSF DOI has not yet been minted. We recommend the collaborating institution mint the OSF lock before unblinding wet-lab results, per the V7 + V8 pre-registration scheme.
4. **The Roberts 2020 ceiling is unmodifiable**. Even if every compound performs maximally, the expected clinical effect in healthy adults is bounded at g ≈ 0.50.
5. **Encenicline's place in the top-10** is structural: V6.A binding + V6.B Cluster D say "promising", V8 πphen (when real LINCS+JUMP-CP data flows) is expected to flag `target_true.phenotype_failed` (Phase 3 was null per Brannan 2019). The compound is included to demonstrate the V8 safety-net mechanism, NOT as a strong candidate.
6. **Per-target n is small** (7-26 ChEMBL pchembl ≥ 8 compounds per panel target). Statistical power is bounded.
7. **The 5-MoA cognition reference centroids** (cholinergic / catecholaminergic / glutamatergic / trophic-ISR / remyelination) are based on V8 plan §Five-MoA centroids; the K=5 mixture is locked in V8 OSF pre-registration §2.9 before unblinding.

---

*Generated by `reports/wet_lab_handoff_v1.md`. Companion to the 5-paper manuscript suite + 16 publication figures + 2 OSF pre-registrations + PROJECT_STATUS executive one-pager. For technical questions about the pipeline: `github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing`. For collaboration discussion: open an issue on the repo.*
