# Graczyk-Style Gini + S(10×) Selectivity Scoring and Multi-Class Faceted Wet-Lab Shortlist Generator for the MAMMAL Cognitive Enhancement Repurposing Pipeline

## TL;DR

- Build the §7.4 selectivity layer as a **Graczyk Gini + S(10×) tandem** computed on calibrated MAMMAL pKd vectors over the 22-target panel, with rank-percentile substitution at the 4 MAMMAL_ONLY_INVERTED targets (SLC6A3, SLC6A2, GRIN2A, GRIN2B) so anti-correlated channels contribute selectivity information without poisoning the affinity vector with sign-flipped values. Categorize each compound as mono / dual / poly / panel-flat; flag panel-flat as the HRH3-lock-in signature, not a hit.
- Build the §8.1 reranker as a **dual-faceted shortlist** (8 mechanism-class facets + 9 targeted-pair facets, top-5 per facet) ranked by `efficacy × Gini-selectivity` within class and `(eff_A + eff_B) × off-pair penalty` within pair. This surfaces ~40 candidates spanning 8 mechanism classes and dissolves the v3 HRH3-23/25 single-target lock-in into a structured, mechanism-orthogonal shortlist that medicinal chemists can triage.
- Honest framing: this is a transparency tool, not a discovery tool. Underlying scores don't change. Roberts CA et al. (*Eur Neuropsychopharmacol* 2020, 38:40–62) found that methylphenidate — the best-evidenced pharmacological cognitive enhancer in healthy adults — has a pooled SMD of just 0.21 (k=24 studies, 47 effect sizes, p=0.0004), with modafinil at SMD=0.12 (k=14/64) and d-amphetamine null (k=10/27). Faceted ranking widens the candidate diversity that medicinal chemists triage; it cannot break that small-effect ceiling.

---

## Key Findings

1. **Graczyk Gini is the right primary metric for this panel** despite its kinase-superfamily origin. The Gini coefficient is dimensionless, threshold-free, monotone under rank-preserving transforms of the affinity vector, and produces a single-number scorecard that medicinal chemists already read. Graczyk 2007 (*J Med Chem* 50:5773–5779, doi:10.1021/jm070562u) demonstrated it on 40 commercial kinase inhibitors across 85 kinases (staurosporine Gini = 0.150, PD184352 Gini = 0.905). Bosc, Meyer & Bonnet 2017 (BMC Bioinformatics 18:17, doi:10.1186/s12859-016-1413-y) reproduced and benchmarked it against entropy, partition index, window score, and ranking score, confirming Gini's robustness when paired with a complementary threshold-based metric. S(10×) is that complement: from the Karaman et al. 2008 *Nat Biotechnol* selectivity-score family, defined (per Uitdehaag & Zaman 2011, PMC3100252) as *"the number of kinases hit at 10 times the Kd of the target … divided by the number of kinases tested."* Gini captures *distribution shape*; S(10×) captures *the local cliff* between the top target and its neighbors. Together they classify mono-, dual-, and polypharm profiles cleanly.

2. **Entropy and Partition Index were considered and rejected as primary** but kept as secondary diagnostics. Uitdehaag & Zaman 2011 (BMC Bioinformatics 12:94, PMC3100252) showed selectivity entropy (Ssel) is more consistent than Gini under ties and panel-size variation; Cheng et al. 2010 (J Med Chem 53:4502–4510, doi:10.1021/jm100301x) showed the partition index has thermodynamic grounding but requires a designated reference target — which is exactly what we don't have for a CNS panel where the "true" lead target is the question. Gini wins because it's panel-agnostic; S(10×) wins because medicinal chemists understand "compounds within a log of the lead."

3. **MMCLKin (Nat Commun 2025, s41467-025-65869-8)** is the relevant 2025 contrastive-learning prior for kinase selectivity prediction. It explicitly evaluates predictions against Gini, partition index, and selectivity entropy as ground-truth metrics — precedent for treating Gini as a *scoring* layer separate from the *prediction* layer. We list it as future direction; it is not a replacement for Gini, it is a downstream consumer of the same metric.

4. **Inverted-target handling is the load-bearing design choice**. The 4 MAMMAL_ONLY_INVERTED targets (SLC6A3 ChEMBL238, SLC6A2 ChEMBL222, GRIN2A ChEMBL1628474, GRIN2B ChEMBL1972) are anti-correlated with ChEMBL ground truth. Three options exist; we recommend **rank-percentile substitution** at these targets. Option (a) drop-them loses information for an entire mechanism class (dopamine reuptake) we care about for cognition. Option (b) isotonic-corrected pKd assumes the calibration model generalizes off the calibration set, which is unverified. Option (c) rank-percentile uses only ordering information, preserves the Gini's rank-monotonicity property, and is robust to the sign flip — at these 4 targets we substitute `affinity_signal = (rank_in_panel / 22)` for raw pKd. This is the cleanest engineering choice and matches Phase A.7's calibrated-weight semantics.

5. **The v3 HRH3 lock-in is a calibration artifact, and Gini + facets will surface it as such.** When 23/25 compounds list `mammal_best_target = HRH3`, Gini computed on raw uncalibrated pKd would give them all Gini > 0.7 (mono-selective at HRH3). Gini computed on calibrated affinity will instead show a bimodal distribution: pitolisant and a few legit H3 inverse agonists at high Gini; the rest at panel-flat Gini < 0.3 because their HRH3 affinity is barely above noise floor for the other 21 targets. The "panel-flat" category is therefore a *diagnostic for low-confidence MAMMAL outputs*, not a biological finding.

6. **Faceted shortlist breaks the single-ranking degeneracy.** A top-25 by composite score will keep returning the same HRH3-dominant set; top-5 per mechanism class returns 40 candidates spanning 8 classes. The targeted-pair facets (CHRNA7+ACHE, HRH3+DRD1, GRIA+PDE4D, SIGMAR1+NTRK2, DAT+NET, etc.) test biological hypotheses no single-ranking can express. This is the deliverable.

---

## Details

### Section 1 — Selectivity Scoring Methodology (§7.4)

#### 1.1 Graczyk Gini: origin, formula, why it fits

Graczyk PP, *J Med Chem* 2007, 50:5773–5779 (doi:10.1021/jm070562u). Verbatim from the abstract: *"A novel application of the Gini coefficient for expressing selectivity of kinase inhibitors against a panel of kinases is proposed. … Nonselective inhibitors are characterized by Gini values close to zero (Staurosporine, Gini 0.150). Highly selective compounds exhibit Gini values close to 1 (PD184352 Gini 0.905). The relative selectivity of inhibitors does not depend on the ATP concentration."* The method: sort %-inhibition values ascending, compute the Lorenz cumulative-fraction curve, define Gini as the ratio of the area between the Lorenz curve and the line of equality to the area below the line of equality. Gini = 0 means perfectly promiscuous; Gini = 1 means perfectly selective.

**Adaptation to the cognition panel.** The math doesn't care that our 22 targets span four protein superfamilies (GPCRs HRH3/HCRTR1/HCRTR2/ADRA2A/DRD1; transporters SLC6A2/SLC6A3; ion channels HCN1/KCNQ2/KCNQ3 + ionotropic glutamate GRIA1–4/GRIN2A/GRIN2B/CHRNA7; enzymes ACHE/PDE4D/PDE9A; trophic receptor NTRK2; sigma chaperone SIGMAR1). What matters is that we have a 22-dimensional affinity vector per compound. Gini operates on rank-ordered values, not biological identity. The only structural risk: the cognition panel's prior affinity distribution is broader than a kinase panel's (kinases share an ATP pocket; CNS receptors don't), so baseline Gini will sit higher across compounds than the ~0.4–0.6 typical for kinase profiling. Calibrate empirically against the positive controls in §1.5.

#### 1.2 S(10×) as the threshold complement

For our panel (N=22): S(10×) is the count of targets t where pKd(t) ≥ pKd(top) − 1.0. Range 1 (mono-selective; only the top target within a log) to 22 (panel-flat; everything within a log).

#### 1.3 Treatment of MAMMAL_ONLY_INVERTED — the recommended approach

For each compound, build the **selectivity vector** S as follows:

```
For t in 22_targets:
    if t in {SLC6A3, SLC6A2, GRIN2A, GRIN2B}:           # MAMMAL_ONLY_INVERTED
        S[t] = rank_percentile(pKd[t], dist_t)           # 0..1 quantile w.r.t. panel-prior distribution
        S[t] = 5.0 + 4.0 * S[t]                          # rescale to pKd-like 5..9 range
    elif t in {DRD1, HCRTR1}:                            # MAMMAL_ONLY_STRONG; use raw
        S[t] = pKd[t]
    elif t == HRH3:                                       # BOLTZ_2X_MAMMAL; small-n caveat; use raw but cap weight downstream
        S[t] = pKd[t]
    else:
        S[t] = pKd[t]                                    # MAMMAL_DOWNWEIGHTED; raw OK for rank-based metric
```

Then `gini(S)` and `S10x(S)` directly. Negative or zero predicted pKd values (MAMMAL outputs below the prior mean of 5.79) are clipped to the panel minimum prior (e.g. 4.0) so the Lorenz curve is well-defined; this matches the Olivia Guest reference numpy implementation (github.com/oliviaguest/gini), which handles non-positive values by shifting before sort.

#### 1.4 Code skeleton

```python
# src/mammal_repurposing/selectivity/gini_scorecard.py
import numpy as np
import pandas as pd

INVERTED = {"SLC6A3", "SLC6A2", "GRIN2A", "GRIN2B"}
PANEL = ["CHRNA7","ACHE","GRIA1","GRIA2","GRIA3","GRIA4","GRIN2A","GRIN2B",
         "DRD1","SLC6A3","ADRA2A","SLC6A2","HRH3","HCRTR1","HCRTR2",
         "PDE4D","PDE9A","NTRK2","SIGMAR1","KCNQ2","KCNQ3","HCN1"]

def gini(x: np.ndarray) -> float:
    """Graczyk-style Gini on a non-negative affinity vector."""
    x = np.asarray(x, dtype=float).flatten()
    if x.min() < 0:
        x = x - x.min()
    x = x + 1e-9
    x = np.sort(x)
    n = x.shape[0]
    idx = np.arange(1, n + 1)
    return float(((2*idx - n - 1) * x).sum() / (n * x.sum()))

def s_10x(pkd_vec: np.ndarray) -> int:
    """Number of panel members within 1 log unit of the top target's pKd."""
    pkd_vec = np.asarray(pkd_vec, dtype=float)
    top = pkd_vec.max()
    return int((pkd_vec >= top - 1.0).sum())

def selectivity_vector(pkd_row: pd.Series, panel_prior: dict) -> np.ndarray:
    out = np.empty(len(PANEL), dtype=float)
    for i, t in enumerate(PANEL):
        v = pkd_row[t]
        if t in INVERTED:
            q = panel_prior[t].quantile_of(v)   # rank percentile in [0,1]
            out[i] = 5.0 + 4.0 * q
        else:
            out[i] = max(v, 4.0)                # clip noise floor
    return out

def categorize(gini_val: float, s10: int) -> str:
    if gini_val >= 0.7 and s10 <= 2: return "mono"
    if 0.5 <= gini_val < 0.7 and 3 <= s10 <= 5: return "dual"
    if gini_val < 0.5 and s10 >= 6: return "poly"
    if gini_val < 0.3: return "flat"
    return "intermediate"
```

#### 1.5 Validation gates (named positive/negative controls)

| Compound | Expected category | Expected Gini | Expected S(10×) | Rationale (named sources) |
|---|---|---|---|---|
| **Donepezil** | mono | ≥ 0.70 | 1–2 | AChE selectivity **>1000-fold over BChE** (Sugimoto et al. 2000, as cited in ScienceDirect Topics: *"outstanding selectivity toward AChE over BChE, with a selectivity ratio of over 1000"*); only meaningful off-target in panel is σ1 (Kd = 14.6 nM, Kato et al. 1999 / Ishikawa et al. 2009). Expect Gini ~0.75, S(10×) = 2. |
| **Galantamine** | dual / poly | 0.40–0.60 | 4–6 | AChE inhibitor + α7/α4β2 nAChR PAM (Maelicke A et al., *Biol Psychiatry* 2001;49:279–288, PMID 11230879). Cross-reactive within the cholinergic mechanism class. Expect Gini ~0.5. |
| **Modafinil** | dual | 0.50–0.65 | 3–5 | DAT primary + weak HRH3/H1 activity + indirect orexinergic engagement. Expect Gini ~0.55. |
| **Solriamfetol** | dual | 0.55–0.70 | 2–4 | "Selective dopamine and norepinephrine reuptake inhibitor" (Baladi et al. *JPET* 2018, doi:10.1124/jpet.118.248120; DAT Ki = 14.2 μM, NET Ki = 3.7 μM, SERT Ki = 81.5 μM with *"negligible binding affinity… at SERT"*). Maps to SLC6A3+SLC6A2 paired-mechanism. |
| **Pitolisant (Wakix, BF-2649)** | mono | ≥ 0.80 | 1 | Ligneau X et al. (*JPET* 2007;320(1):365–375, doi:10.1124/jpet.106.111039): Ki = 0.16 nM at human H3R, IC50 > 10 μM at H1 and H4 — **>62,500-fold selectivity over H1/H4**. Should be the cleanest mono-selective in the panel. |
| **BPN14770 / zatolmilast** | mono | ≥ 0.75 | 1 | "First-in-class PDE4D allosteric inhibitor" binding a primate-specific N-terminal region (Gurney et al. *Nat Biotechnol* 2010, 28:63–70; *J Med Chem* 2019, 62:4884–4901). Mono at PDE4D. |
| **ANAVEX 2-73 / blarcamesine** | dual | 0.45–0.60 | 3–5 | "Mixed ligand for sigma1/muscarinic receptors" (Alzforum), SIGMAR1 IC50 = 860 nM (PMC8387417). Off-panel muscarinic activity invisible to our 22-panel — expect panel-Gini biased high. |
| **7,8-DHF** | mono | ≥ 0.65 | 1–2 | "Selective TrkB receptor agonist" (Abcam ab120996); should be NTRK2-dominant. |
| **Encenicline (EVP-6124)** | mono | ≥ 0.70 | 1–2 | α7 nAChR partial agonist with the active program target exclusively CHRNA7. |
| **TC-5619 / bradanicline** | mono | ≥ 0.75 | 1 | Walling D et al. *Schizophr Bull* 2016;42(2):335–343 (PMC4753586): *"TC-5619 is a highly selective alpha7 NNR full agonist with a Ki at the alpha7 NNR of 1 nM."* Note: Mazurov AA et al. (*J Med Chem* 2012;55(22):9793–9809) reports Ki = 1.4 nM — flag the minor discrepancy rather than treating either as uncontested. |
| **Aripiprazole** (negative control) | poly | ≤ 0.40 | ≥ 6 | Pan-aminergic D2/D3/5-HT2A/5-HT1A/5-HT2B/α1/HRH1 — broad polypharm. In our 22-panel only DRD1, ADRA2A, HRH3 would register, so cognition-panel Gini will be artificially high. This is *the structural limitation of any selectivity metric on a narrow panel*, and we report it. |

Gates: pitolisant must rank Gini ≥ 0.80; donepezil ≥ 0.70; modafinil in [0.50, 0.65]. If pitolisant lands below donepezil, calibration is upside-down and the implementation is broken.

#### 1.6 Cross-cluster stress test

Compute Gini on both the MAMMAL pKd vector AND the Boltz-2 pKd vector for the same compound. Disagreement (|ΔGini| > 0.2) is a provenance flag. For HRH3 specifically (Boltz-2 ρ = +0.87 at n=3), MAMMAL-Gini-high but Boltz-Gini-low at HRH3 is the failure mode the §7.1 diagnostic predicts. Report bootstrap CI on Gini (1000 resamples of the 22-target vector, BCa intervals) and flag any compound whose Gini 95% CI spans (0.3, 0.7).

---

### Section 2 — Multi-Class Top-N Faceted Shortlist (§8.1)

#### 2.1 By-mechanism-class facets (8 classes; top-5 each → 40 candidates)

| Class | Targets (ChEMBL IDs) | Expected top-5 |
|---|---|---|
| **Cholinergic** | ACHE (ChEMBL220), CHRNA7 (ChEMBL2492) | donepezil, galantamine, encenicline, TC-5619, ABT-126 |
| **Glutamatergic AMPA** | GRIA1 (ChEMBL2093872), GRIA2 (ChEMBL2093867), GRIA3, GRIA4 | aniracetam, tulrampator (CX1632 / S47445), CX-717, CX-1739, piracetam |
| **Glutamatergic NMDA** | GRIN2A (ChEMBL1628474), GRIN2B (ChEMBL1972) | memantine, lanicemine, traxoprodil (CP-101,606), ifenprodil, GNE-6901 (GluN2A PAM; Hackos et al. *Neuron* 2016 doi:10.1016/j.neuron.2016.01.016) |
| **Dopaminergic** | DRD1 (ChEMBL2056), SLC6A3 (ChEMBL238) | methylphenidate, d-amphetamine, modafinil, solriamfetol, vanoxerine |
| **Noradrenergic** | SLC6A2 (ChEMBL222), ADRA2A (ChEMBL1867) | atomoxetine, reboxetine, guanfacine, clonidine, viloxazine |
| **Histaminergic** | HRH3 (ChEMBL264) | pitolisant, ABT-239, A-349821, JNJ-39220675, BF-2649 analogs |
| **Orexinergic** | HCRTR1 (ChEMBL5113), HCRTR2 (ChEMBL4792) | **flag: panel direction inverted** — FDA-approved orexin drugs (suvorexant, lemborexant) are *antagonists* for sleep, opposite of procognitive direction. Surfaced for review but tagged `WRONG_DIRECTION_FOR_COGNITION`. |
| **Phosphodiesterase** | PDE4D (ChEMBL288), PDE9A (ChEMBL5147) | BPN14770 (zatolmilast), rolipram, BI-409306, PF-04447943, roflumilast |
| **Trophic/sigma/ion-channel ("Other")** | SIGMAR1 (ChEMBL287), NTRK2 (ChEMBL4805), KCNQ2 (ChEMBL3038), KCNQ3 (ChEMBL3138), HCN1 (ChEMBL1075145) | 7,8-DHF, ANAVEX 2-73, retigabine/ezogabine, ivabradine, LM22A-4 |

Within-class ranking: `composite = rrf_efficacy_normalized × gini_within_class_bonus`, where the bonus rewards compounds whose Gini concentrates *on this class* (top-2 panel members are both in this class).

#### 2.2 Targeted-pair facets (9 pairs; top-5 each → up to 45 more candidates, with overlap)

| Pair | Biological hypothesis tested | Expected leaders |
|---|---|---|
| **CHRNA7 + ACHE** | Galantamine-class dual cholinergic — can a single molecule replicate galantamine's PAM + AChEI duality? | galantamine, ABT-126, novel dual scaffolds |
| **PDE4D + CHRNA7** | cAMP signaling rescue + cholinergic — orthogonal LTP / spine-density mechanism | theoretical scaffolds; mostly empty — surfacing emptiness is the point |
| **HRH3 + DRD1** | Dual aminergic for processing speed (modafinil-plus-pitolisant phenotype) | pitolisant + DRD1 partial agonists; combinations |
| **GRIA + PDE4D** | LTP enhancement via AMPA + cAMP convergence | ampakines × rolipram-class hybrids |
| **SIGMAR1 + NTRK2** | Neuroprotection axis (ANAVEX-2-73 + 7,8-DHF combo logic) | 7,8-DHF, blarcamesine, donepezil (cross-reactive σ1) |
| **SLC6A3 + SLC6A2** | Dual reuptake while avoiding SERT (solriamfetol phenotype) | solriamfetol, vanoxerine, NS-2359, atomoxetine, methylphenidate |
| **HCN1 + KCNQ2/3** | Intrinsic excitability tuning | retigabine, ivabradine, ICA-27243 |
| **GRIN2A vs GRIN2B preference** | Procognitive 2A-preferring PAMs (Hackos et al. *Neuron* 2016) | GNE-6901, GNE-8324 — low-priority because Boltz-2 ρ at GRIN2A is in the inverted set |
| **HCRTR1 + DRD1** | Motivation/arousal axis — both must be agonist-direction | sparse; mostly research compounds |

Within-pair ranking: `(pkd_A + pkd_B) − λ × max(pkd_off_pair)` where λ = 0.5 penalizes off-pair affinity. Pairs are pre-registered; we don't combinatorially enumerate all C(22,2) = 231 pairs.

#### 2.3 Display format and cross-facet provenance

One markdown report `reports/wet-lab/wet_lab_shortlist_v4_faceted.md`:
- Section per facet
- Columns: rank, compound, RRF_efficacy, gini, S10x, top_targets, MoA_class, ADMET_clean, liability_v8.0b, regulatory_status, **cross_facet_provenance**
- The cross_facet_provenance column lists every other facet the compound appears in (e.g., "donepezil: cholinergic #1; CHRNA7+ACHE #3; SIGMAR1+NTRK2 #4"). This is the single most important addition — it stops a reviewer from triple-counting one compound as if it were three independent hits.

#### 2.4 Selectivity-stratified ranking inside the facet

Within each facet, rank by `score_facet = composite × (1 + α·gini_within_facet)` with α ≈ 0.3, ties broken by ADMET-clean status. A Gini = 0.8 hit beats a Gini = 0.5 hit at the same raw efficacy by ~12% — small enough that polypharm compounds still surface, large enough that mono-selective compounds get a deserved boost in clean-mechanism facets.

---

### Section 3 — Pipeline Integration (engineering)

#### 3.1 Module layout

```
src/mammal_repurposing/
├── selectivity/
│   ├── __init__.py
│   ├── gini_scorecard.py        # gini(), s_10x(), selectivity_vector(), bootstrap_ci()
│   ├── categorize.py            # mono/dual/poly/flat assignment + provenance flag
│   └── cross_cluster.py         # MAMMAL-Gini vs Boltz-Gini diagnostic
└── fusion/
    ├── faceted_shortlist.py     # by-class + targeted-pair top-N generators
    └── cross_facet.py           # provenance bookkeeping

scripts/
├── 27_v3_selectivity_scoring.py    # reads dti_scores.parquet, writes selectivity columns
└── 28_v3_faceted_shortlist.py      # reads ranking parquet, writes faceted_shortlist.parquet + .md
```

#### 3.2 Schema additions to `dti_scores.parquet` (or sibling `selectivity_scores.parquet`)

| Column | Type | Description |
|---|---|---|
| `gini_coefficient` | float64 | Graczyk Gini on calibrated affinity vector |
| `gini_ci_low`, `gini_ci_high` | float64 | BCa 95% CI from 1000-resample bootstrap |
| `s_10x` | int8 | Count of panel members within 1 log of top pKd |
| `selectivity_category` | str | {mono, dual, poly, flat, intermediate, uncertain} |
| `top_target` | str | argmax of calibrated affinity |
| `top_target_pkd` | float64 | calibrated pKd at top_target |
| `second_target` | str | argmax of affinity excluding top |
| `second_target_pkd` | float64 | |
| `mechanism_class` | str | derived from top_target via panel-class map |
| `gini_boltz` | float64 | parallel computation on Boltz-2 pKd vector |
| `cross_cluster_flag` | bool | `|gini_mammal − gini_boltz| > 0.20` |

#### 3.3 Faceted shortlist parquet

`data/results/v2/faceted_shortlist.parquet`:

| Column | Type |
|---|---|
| `facet_type` | str (`mechanism_class` / `targeted_pair`) |
| `facet_name` | str (`cholinergic`, `HRH3+DRD1`, …) |
| `facet_rank` | int8 (1–5) |
| `compound_id` | str |
| `composite_score` | float64 |
| `gini` | float64 |
| `s_10x` | int8 |
| `selectivity_category` | str |
| `cross_facet_count` | int8 |
| `cross_facet_list` | list[str] |
| `notes` | str (e.g. `WRONG_DIRECTION_FOR_COGNITION` for orexin antagonists; `INVERTED_TARGET_TOP`; `LOW_N_CALIBRATION`) |

#### 3.4 Calibration interaction with §7.11 / Phase A.7

The MAMMAL_ONLY_DOWNWEIGHTED targets (18 of 22) have calibrated weights of 0.30. Two options:

(i) Apply weights *inside* the Gini computation as `gini(w * S)`. **Problem:** this changes the meaning of Gini — it's no longer the selectivity of the compound, it's the selectivity of the *weighted prediction*. Discouraged.

(ii) Apply weights only at the *consumption* layer (faceted ranking) and compute Gini on unweighted calibrated affinity. **Recommended.** Gini stays a property of the compound's predicted profile; the downstream ranker uses weights when picking facet leaders. This decouples the selectivity metric from calibration confidence.

For the 4 INVERTED targets, the rank-percentile substitution (§1.3) handles the sign-flip. The calibrated weight (0.30) is applied at the ranker stage, not the Gini stage.

#### 3.5 Wall-clock

- Gini + S(10×) per compound on 22-vector: ~0.5 ms
- 298 compounds: ~150 ms compute + ~2 s parquet I/O
- Bootstrap CI (1000 resamples × 298 compounds): ~30 s on one core
- Faceted shortlist generation (17 facets × top-5 selection): < 5 s
- **Total wall-clock: under one minute.** I/O dominates. CPU cost is negligible compared to the 2-hour MAMMAL inference + Boltz-2 docking that produced the inputs.

---

### Section 4 — Expected Impact on the v3 Top-25 + Wet-Lab Hypotheses

#### 4.1 Predicted selectivity-category reassignment

Of the 23/25 v3 top compounds currently `mammal_best_target = HRH3`:

- **~5 will land in mono/HRH3** (genuine H3 antagonist scaffolds like pitolisant analogs, ABT-239-likes) — survive the cleanup.
- **~15 will land in panel-flat (Gini < 0.3)** because their HRH3 affinity is only marginally above their predicted noise floor at the other 21 targets. These get *demoted* from the wet-lab shortlist — they were ranked by the lock-in artifact, not by selectivity.
- **~3 will reassign to a different top_target** when calibrated affinity is computed properly (e.g. compounds with strong DRD1 or SIGMAR1 predictions buried under HRH3 in the v3 fusion).

#### 4.2 Per-facet expected leaders (decision-grade predictions)

- **HRH3 facet top-5:** pitolisant (FDA-approved Wakix; Ki = 0.16 nM at H3R, >62,500-fold over H1/H4 per Ligneau et al. *JPET* 2007 doi:10.1124/jpet.106.111039), ABT-239, A-349821, JNJ-39220675, pitolisant analogs from Bioprojet pipeline.
- **CHRNA7 facet top-5:** encenicline (EVP-6124), TC-5619/bradanicline (Ki = 1 nM per Walling 2016 / 1.4 nM per Mazurov 2012 — minor lit discrepancy), ABT-126, AVL-3288, GTS-21/DMXB-A.
- **ACHE facet top-5:** donepezil, galantamine, rivastigmine, huperzine A, donepezil-σ1 hybrids (Estrada flavonoids).
- **PDE4D facet top-5:** BPN14770/zatolmilast (allosteric NAM at primate-specific N-terminus; Phase 2 FXS NCT03569631 met cognitive secondary endpoints), rolipram, roflumilast (CNS-penetration caveat), BI-409306, GEBR-32a.
- **SIGMAR1 facet top-5:** ANAVEX 2-73 / blarcamesine (σ1 IC50 = 860 nM; per Macfarlane et al. *J Prev Alzheimers Dis* 2025;12(1):100016, doi:10.1016/j.tjpad.2024.100016, the AD-004 Phase 2b/3 trial "significantly slowed clinical progression by 36.3% at 48 weeks … on the prespecified primary cognitive endpoint ADAS-Cog13"), donepezil (Kd = 14.6 nM at σ1, Kato 1999), pridopidine, dimemorfan, PRE-084.
- **NTRK2 facet top-5:** 7,8-DHF, LM22A-4, R13, deoxygedunin, ANA-12 (antagonist — flag DIRECTION).
- **DAT facet top-5:** methylphenidate, modafinil, vanoxerine/GBR-12909, solriamfetol, sertraline (incidental SERT primary).
- **DAT+NET pair top-5:** solriamfetol (DAT Ki=14.2 μM, NET Ki=3.7 μM, "negligible binding affinity… at SERT" Ki=81.5 μM, Baladi 2018), methylphenidate, atomoxetine (NET-dominant), nisoxetine, NS-2359 (triple-reuptake).
- **SIGMAR1+NTRK2 pair top-5:** 7,8-DHF, blarcamesine, donepezil (genuine cross-reactive at σ1), Estrada donepezil-flavonoid hybrids, ANAVEX 3-71.
- **GRIA+PDE4D pair top-5:** Mostly empty/theoretical. Reporting this emptiness is the point — it tells the medicinal chemist that no current compound bridges this hypothesis and one would have to be synthesized.

#### 4.3 The recommended 5-compound wet-lab radioligand binding shortlist

If Pierce can afford only 5 binding assays (chosen by **structural novelty + mechanism orthogonality**):

1. **BPN14770 (zatolmilast)** × PDE4D — mono-selective, allosteric (novel mechanism), already in Phase 2/3 Fragile X and AD, validates PDE4D Gini gate.
2. **7,8-DHF** × NTRK2 — trophic mechanism, orthogonal to everything else, mono-selective.
3. **Pitolisant** × HRH3 — FDA-approved, validates the HRH3 facet (and serves as positive control for "real H3 inverse agonist vs panel-flat noise" categorization).
4. **Solriamfetol** × SLC6A3 + SLC6A2 — tests the dual-reuptake facet AND validates the inverted-target rank-percentile substitution (since both SLC6 transporters are INVERTED). Highest-value diagnostic burn.
5. **Blarcamesine (ANAVEX 2-73)** × SIGMAR1 — tests the SIGMAR1+NTRK2 pair (without buying 7,8-DHF *and* blarcamesine both); novel mechanism with recent Phase 2b/3 readout supporting σ1 modulation in AD.

Spends 3 mono-selective slots (one per top-3 mechanism class: PDE/trophic/aminergic) and 2 dual-selective slots (DAT+NET, σ1+TrkB-adjacent), matching the user's prioritization rubric exactly.

---

### Section 5 — Validation Gates & Publication Angle

#### 5.1 Hard validation gates (must-pass before publishing v4)

- **Gate G1 (positive controls):** Donepezil Gini ≥ 0.70 AND pitolisant Gini ≥ 0.80 AND BPN14770 Gini ≥ 0.75. If any fail, the implementation is broken.
- **Gate G2 (negative controls):** Aripiprazole Gini ≤ 0.40. If aripiprazole comes back mono-selective, the panel-flat detector is broken.
- **Gate G3 (faceted CHRNA7):** Top-5 must contain TC-5619 AND encenicline.
- **Gate G4 (faceted ACHE):** Top-5 must contain donepezil AND galantamine.
- **Gate G5 (faceted HRH3):** Top-5 must contain pitolisant.
- **Gate G6 (cross-facet hygiene):** Pitolisant appears in HRH3 facet but NOT in HRH3+DRD1 pair facet (it has no DRD1 affinity).
- **Gate G7 (bootstrap stability):** ≥ 80% of compounds have Gini 95% CI width < 0.30. If too many CIs are wide, the panel is undersized for stable selectivity assessment (which is itself a finding to report).

#### 5.2 Publication angle

- **Working title:** "Selectivity-stratified faceted prioritization in foundation-model drug repurposing: dissolving single-target lock-in via per-mechanism-class top-N ranking."
- **Venue priorities:** *J. Cheminform.* (open-access, code-friendly, methods-heavy) → *Drug Discovery Today* (broader medicinal-chemistry audience) → *Mol. Inform.* (fallback).
- **Novelty claim:** Most computational repurposing pipelines (DTINet, DeepPurpose, KIBA-style benchmarks) collapse to a single ranked list. We show that a faceted output mirrors the actual triage workflow medicinal chemists use, and that pairing Graczyk Gini with selectivity-class facets specifically diagnoses calibration-driven single-target lock-in. This is a *method-and-diagnostic* paper, not a discovery paper.
- **Reproducibility:** Full code under MIT license; the 22-target panel definition with all ChEMBL IDs, the calibrated weights from Phase A.7, `dti_scores.parquet`, and `faceted_shortlist.parquet` released as supplementary.
- **Honest framing in discussion:** (a) HRH3 lock-in was a calibration artifact, not a biological finding; faceted ranking is a transparency tool, not a discovery method. (b) Roberts CA et al. (*Eur Neuropsychopharmacol* 2020, 38:40–62, doi:10.1016/j.euroneuro.2020.07.002) reported verbatim: *"There was an overall effect of modafinil (SMD=0.12, p=.01) … There was an overall effect of MPH (SMD=0.21, p=.0004) … There were no effects for d-amph"* (k=14/64 modafinil; k=24/47 MPH; k=10/27 d-amph). This small-effect ceiling applies to any compound this pipeline surfaces. The pipeline's job is to widen mechanism diversity for downstream biology, not to claim large effect sizes are achievable.

---

### Section 6 — Lateral Considerations

1. **Selectivity-vs-efficacy Pareto front within each facet.** Inside the cholinergic facet, donepezil (high-Gini, moderate efficacy) and galantamine (mid-Gini, broader mechanism) are on the same Pareto front. The §8.0a Pareto restructure should be applied facet-locally, not globally. Concretely: within each top-5 facet table, mark Pareto-dominated compounds with a strikethrough but keep them visible (a strict top-5 by composite score can hide a Pareto-optimal but lower-composite hit). Implementation: ~15 lines, `is_pareto_dominated(eff, gini, others)`.

2. **Time-stratified selectivity (§7.3 GWAS-anchored panel expansion).** A compound's Gini is computed against the *current* 22-target panel. As the panel expands (planned via §7.3 GWAS-anchored additions — e.g., SHANK3, FMR1-related pathways, new schizophrenia GWAS hits), Gini values will *drop systematically* (more targets → more chances for off-target affinity). Recommendation: version Gini values by panel hash (e.g. `gini_panel22h0a4c`, `gini_panel36h7e21`) so longitudinal comparisons across pipeline versions remain meaningful. Don't let v4 vs v5 Gini comparisons be conflated.

3. **Selectivity drift under calibration uncertainty (bootstrap CI on Gini).** Already in §1.6 and gate G7, but worth foregrounding: a compound with Gini = 0.55 but 95% CI = (0.32, 0.78) is uncategorizable — it could be dual or panel-flat. The `selectivity_category` for such compounds should be `uncertain`, not forced into mono/dual/poly. Honest, and prevents downstream consumers from acting on overconfident categorizations.

4. **Cross-cluster selectivity provenance.** When MAMMAL-Gini and Boltz-Gini for the same compound differ by > 0.20, the provenance flag fires. This is a *signal*, not a bug: it means the two foundation models disagree about which targets the compound hits hardest. For HRH3-bait compounds, MAMMAL-Gini will be artificially high (HRH3 is a strong-positive target after calibration) and Boltz-Gini will be moderate (Boltz docks pockets across the full panel without the calibration bias). The disagreement is a calibration-confidence signal that should propagate to wet-lab triage: high-disagreement compounds get *lower* priority.

5. **Negative-selectivity provenance for INVERTED targets.** A compound that appears mono-selective at SLC6A3 (DAT) in raw MAMMAL output is, per the §7.1 diagnostic, *misranked* — MAMMAL predicted high affinity but ChEMBL says low. The rank-percentile substitution handles this at the Gini level, but the wet-lab triage layer must also tag these compounds with `INVERTED_TARGET_TOP`. Concretely: if `top_target ∈ {SLC6A3, SLC6A2, GRIN2A, GRIN2B}` AND `Boltz_top_target ≠ MAMMAL_top_target`, downgrade priority and require the Boltz score as tiebreaker. This is the single most important safety check for wet-lab spend — it stops the pipeline from recommending a binding assay on a target where the model's prediction is anti-correlated with ground truth.

---

## Recommendations

**Stage 0 (this week, half-day):** Implement `gini_scorecard.py` and `categorize.py`. Run on the existing `dti_scores.parquet`. Print the gate G1/G2 table (donepezil/pitolisant/BPN14770/aripiprazole Gini values) and decide pass/fail before touching the shortlist generator. **Threshold to proceed:** all four positive controls in their predicted Gini bands.

**Stage 1 (one day):** Implement `faceted_shortlist.py` for the 8 mechanism-class facets only (skip targeted-pair facets initially). Verify gates G3–G5. **Threshold to proceed:** TC-5619 + encenicline appear in CHRNA7 top-5; donepezil + galantamine in ACHE top-5; pitolisant in HRH3 top-5.

**Stage 2 (one day):** Add the 9 targeted-pair facets, cross-facet provenance bookkeeping, and the markdown report generator. Run bootstrap CIs on Gini (gate G7).

**Stage 3 (half-day):** Add the cross-cluster MAMMAL-vs-Boltz Gini diagnostic. Tag `INVERTED_TARGET_TOP` flags. Generate `wet_lab_shortlist_v4_faceted.md`.

**Stage 4 (review session):** Hand the v4 shortlist + the v3 top-25 to a medicinal chemist (yourself / collaborator). Confirm the dissolution of HRH3 lock-in is real and not relabeling. If still HRH3-dominated across facets, the issue is upstream (the MAMMAL DTI head needs retraining, not reranking).

**Stop conditions:**
- If gate G1 fails (donepezil Gini < 0.65): selectivity vector construction is wrong; do not proceed.
- If gate G7 fails (median Gini CI width > 0.35): the 22-panel is too small for stable selectivity; expand panel per §7.3 GWAS-anchored before publishing.
- If after Stage 4 the faceted shortlist still shows > 50% HRH3 compounds across facets, the lock-in is in the prediction layer, not the ranking layer; do not paper over it with more reranking.

---

## Caveats

- **Panel narrowness.** A 22-target panel is small for Graczyk Gini (Karaman 2008 used 317 kinases; Bosc 2017 used 451). On a 22-panel, the noise floor on Gini is ~0.10–0.15 from finite-sample variance alone. Bootstrap CIs are mandatory; categorical labels (mono/dual/poly) are *probabilistic*, not deterministic.
- **Inverted-target rank-percentile is a heuristic.** Rank-percentile substitution preserves ordering but loses absolute-affinity information at the 4 INVERTED targets. If a compound is truly mono-selective at DAT (e.g., vanoxerine), our metric may underestimate its Gini because we replaced its strong DAT pKd with a normalized percentile. Known tradeoff against the alternative (sign-flipped pKd polluting the vector). Sensitivity-analyze by recomputing Gini under all three options (drop / isotonic / rank-percentile) and reporting the spread.
- **Cross-superfamily Gini interpretation.** Graczyk's intuition was built on a single superfamily (kinases). On a mixed-superfamily panel, "selectivity Gini = 0.7" doesn't have the same biological meaning as in a kinase panel. We must say this explicitly in the J. Cheminform. submission and not let reviewers carry over the kinase intuition.
- **The Roberts 2020 ceiling.** Methylphenidate SMD = 0.21 in healthy adults (Roberts CA et al., *Eur Neuropsychopharmacol* 2020, 38:40–62, doi:10.1016/j.euroneuro.2020.07.002; p=0.0004; k=24 studies, 47 effect sizes); modafinil SMD = 0.12 (k=14/64); d-amphetamine null (k=10/27). Faceted-ranking enrichment cannot increase the true biological effect size of cognitive enhancement — it can only widen the candidate diversity entering wet-lab triage. Report this in discussion to forestall reviewer pushback that the method overpromises.
- **Boltz-2 ρ at HRH3 is n=3.** The "BOLTZ_2X_MAMMAL" designation for HRH3 rests on a Spearman ρ = +0.87 computed on n=3 reference compounds. This is dangerously small. The faceted-HRH3 top-5 prediction (pitolisant + ABT-239 etc.) is robust because it's anchored in literature pharmacology, not in the n=3 ρ. But any *novel* HRH3 hit surfaced by the pipeline should be flagged `LOW_N_CALIBRATION` and prioritized below pitolisant analogs.
- **Aripiprazole panel-flat is partly an artifact.** Aripiprazole's true polypharm involves D2/D3/5-HT2A/5-HT1A/α1/HRH1/5-HT2B — most of which are *not in our 22-target panel*. Its panel-Gini will appear higher than its true cross-CNS Gini would. The §8.0b decision to cut aripiprazole was correct on those grounds; Gini should not be used to re-include it.
- **Orexin direction-of-effect.** HCRTR1/HCRTR2 are in the panel because the GWAS-anchored hypothesis is procognitive *agonism*, but the only FDA-approved orexin drugs (suvorexant, lemborexant) are *antagonists* for sleep. Any orexin-facet leader that's a known antagonist must be tagged `WRONG_DIRECTION_FOR_COGNITION` and excluded from wet-lab top-5 unless the medicinal chemist explicitly wants an antagonist scaffold to chemically invert.
- **Minor literature discrepancy on TC-5619 Ki.** Walling 2016 (PMC4753586) reports Ki = 1 nM at α7 nAChR; Mazurov et al. *J Med Chem* 2012 reports Ki = 1.4 nM. Same compound, both within the expected nM-range for the program; report both rather than picking one silently.