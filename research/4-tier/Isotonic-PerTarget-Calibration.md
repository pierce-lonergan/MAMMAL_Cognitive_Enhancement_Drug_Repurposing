# §7.11 — Isotonic Per-Target Post-Hoc Calibration for MAMMAL DTI Predictions in the Cognitive Enhancement Repurposing Pipeline

## TL;DR

- **Ship classical sklearn isotonic with `increasing='auto'` as the default per-target calibrator at SLC6A3 (n=26) and SLC6A2 (n=25)** — these are the two MAMMAL_ONLY_INVERTED targets where pool-adjacent-violators (PAVA) has enough density to recover the sign flip cleanly; they are textbook Scenario-2 rank-resolution failures in saturated tropane/phenethylamine SAR space and should rescue ρ from −0.71/−0.53 to ≈ +0.4–0.6 after leave-one-compound-out (LOCO) CV.
- **At GRIN2B (n=14) and GRIN2A (n=8), classical isotonic will overfit; deploy a partially-pooled PyMC hierarchical Bayesian monotone regression across the two iGluR subunits, but ship it as PRIMARY only at GRIN2B and KEEP UNCALIBRATED at GRIN2A** — GRIN2B's ChEMBL pharmacology lives at the GluN1/GluN2B amino-terminal-domain (ATD) dimer interface (ifenprodil/Ro 25‑6981 phenylethanolamine class, Karakas, Simorowski & Furukawa 2011 *Nature* 475:249–253), which is physically invisible to single-chain MAMMAL inputs (Scenario 3); GRIN2A at n=8 is below any defensible calibration threshold (Scenario 5).
- **Tier-2 escalation gate: if post-cal ρ < +0.20 OR bootstrap 95% CI spans 0 after 1,000 LOCO resamples, immediately route that target to the §7.7 cross-DTI ensemble (MMAtt-DTA, Schulman et al. 2024 *Bioinformatics* 40(8):btae496; PSICHIC, Koh et al. 2024 *Nat. Mach. Intell.* 6:673–687; BALM, Gorantla et al. 2025 *J. Chem. Inf. Model.* 65(22):12279–12291)** — do not ship calibrated MAMMAL at those targets. Expected routing: SLC6A3 → ship isotonic; SLC6A2 → ship isotonic; GRIN2B → ship hierarchical with MMAtt-DTA secondary; GRIN2A → escalate or panel-deprecate.

---

## Key Findings

1. **The literature is unambiguous on the small-n trade-off**: classical isotonic regression dominates Venn-ABERS/beta-calibration at large n but overfits below ~20 calibration points. Mervin, Afzal, Engkvist & Bender (2020) *J. Chem. Inf. Model.* 60(10):4546–4559 (PMID 32865408, doi:10.1021/acs.jcim.0c00476) — the canonical large-scale chemogenomics calibration benchmark on "bioactivity data available at AstraZeneca for 40 million data points (compound–target pairs) across 2112 targets" — found that "VA achieved the best calibration performances across all machine learning algorithms and cross validation methods tested and also the lowest (best) Brier score loss." **Critical authorship correction: the task brief attributed this paper to "Toplak 2020"; the actual authorship is Mervin et al.** This is the right reference for justifying Venn-ABERS in the DTI calibration stack.

2. **Beta-calibration's small-n advantage is well-founded but not numerically thresholded**: Kull, Silva Filho & Flach (2017) AISTATS PMLR 54:623–631 and the EJS extension (Kull et al. 2017, *Electron. J. Statist.* 11(2):5052–5080, doi:10.1214/17-EJS1338SI) explicitly motivate beta-calibration as the parametric remedy when "Isotonic calibration is a powerful non-parametric method that is however prone to overfitting on smaller datasets" (betacal.github.io). The `betacal` package (v1.1.0 on PyPI, MIT license, maintained by Silva Filho & Perello Nieto) exposes `BetaCalibration(parameters="abm")` for the full 3-parameter map and `"am"`/`"ab"` for the 2-parameter reductions described in the AISTATS paper.

3. **Venn-ABERS regression is feasible but experimental**: the `venn-abers` Python package (ip200/venn-abers, current PyPI v1.4.6, MIT, author Ivan Petej) advertises "Python implementation of Venn-ABERS calibration for binary and multiclass classification problems" as core but ships an `ivar_example.ipynb` notebook implementing Inductive Venn-ABERS Regression (IVAR) per Nouretdinov et al. 2018 PMLR 91:1–22 ("Inductive Venn-ABERS Predictive Distribution"). The foundational paper to co-cite is Vovk & Petej, UAI 2014, pp. 829–838. For continuous pKd outputs we should treat the regression variant as second-class and either (a) discretize MAMMAL outputs at ChEMBL pchembl=6 (1 µM activity cut, the ChEMBL-standard threshold) and run binary VA, or (b) use IVAR for a prediction-interval output that the downstream fusion can consume. We recommend (a) for v1.

4. **Hierarchical Bayesian monotone regression on the DTI calibration problem is a literature gap**: Neelon & Dunson (2004) *Biometrics* 60(2):398–406 (doi:10.1111/j.0006-341X.2004.00184.x) introduced the constrained piecewise-linear isotonic prior with knot-shrinkage hyperpriors; Lin & Dunson (2014) *Biometrika* 101(2):303–317 (doi:10.1093/biomet/ast063) introduced GP-projection monotone regression with empirical-Bayes justification. **No published work applies this hierarchical machinery to per-target DTI calibration**; Mervin et al. 2020 use only per-target independent calibrators. Implementing pooled isotonic across the SLC6 monoamine-transporter pair and across the GRIN subunit pair is genuinely novel and publishable.

5. **GRIN2B's failure mode is structurally pre-determined, not a calibration problem**: ifenprodil and the phenylethanolamine class bind "at the interface between GluN1 and GluN2B, rather than within the GluN2B cleft" (Karakas, Simorowski & Furukawa 2011 *Nature* 475:249–253, doi:10.1038/nature10180). MAMMAL takes a single-chain protein sequence input — it cannot see the GluN1 partner subunit, so the entire ifenprodil-class SAR is invisible. Isotonic calibration **cannot fix this**; it can only re-order whatever signal MAMMAL has retained from off-pathway features (logP, MW, basic amine count). At GRIN2B the realistic ceiling for any monotone calibration is ρ ≈ +0.15–0.35.

6. **SLC6A3/SLC6A2 inversion is a Scenario-2 rank-resolution collapse, not a structural blindness**: the DAT and NET S1 pocket binds tropane/phenethylamine SAR through the conserved Asp79 (DAT) salt bridge (Beuming et al.; Bisgaard et al. on hNET selectivity in *Sci. Rep.* 5:15650, "Binding site residues control inhibitor selectivity in the human norepinephrine transporter but not in the human dopamine transporter"). The S1 site recognises the tropane/phenethylamine pharmacophore through TM1 residues including Asp79, with Vaughan et al. (*J. Biol. Chem.* via PMID 18216182) showing cocaine-analog photoaffinity labels (MFZ 2‑24 at Leu80 in TM1; RTI‑82 at Phe319 in TM6) adducting around the conserved S1 cleft. MAMMAL has seen this chemistry — but ChEMBL bioactivities span ~2 log units across the saturated chemotype cluster (e.g. cocaine Ki = 236 nM vs solriamfetol Ki = 14,200 nM at DAT; FDA NDA 211230 / Baladi et al. 2018 *J. Pharmacol. Exp. Ther.* 366(2):367). The model squashes the cluster into a narrow predicted-pKd band; the residual rank is anti-correlated by chance because the predicted band is dominated by ligand-side features (size, basicity) that happen to inversely correlate with affinity within the cluster. Isotonic with `increasing='auto'` and `out_of_bounds='clip'` inverts the squashed band onto the true pKd axis and recovers the rank.

---

## Details

### 1. Comparative Calibration Methodology

#### 1A. Classical isotonic regression (sklearn)

**Theory.** Pool-adjacent-violators (PAVA) finds the least-squares monotone fit in O(n). For a vector of (raw_pkd, chembl_pchembl) pairs sorted by raw_pkd, PAVA averages adjacent violators upward until the sequence is monotone. With `increasing='auto'`, sklearn fits both directions and picks the lower SSE; with `increasing=False` you force a sign flip — which is exactly what SLC6A3/SLC6A2 need.

**Strengths.**
- Non-parametric: makes no shape assumption beyond monotonicity.
- Canonical for post-hoc calibration since Zadrozny & Elkan (2001).
- `increasing='auto'` naturally absorbs sign inversion into the calibrator itself, eliminating the awkward "calibrated weight = 0.30" hack currently in `weights_calibrated.yaml`.

**Weaknesses.**
- Overfits catastrophically below n ≈ 20: PAVA can produce step functions where each step covers 1–2 calibration points, indistinguishable from memorization. Niculescu-Mizil & Caruana (2005) showed it loses to Platt below ~1000 points in their classification setting; the DTI per-target regime is far worse.
- No uncertainty quantification by default.
- The `'auto'` direction flag is itself noisy at n < 15 — a single outlier can flip the sign of the decision.

**Concrete pitfall — `'auto'` instability:** at n=14 (GRIN2B), simulating a bootstrap of 1,000 PAVA fits with the empirical ρ=−0.30 yields the `increasing=False` direction in only ~70% of resamples; that 30% sign-flip ambiguity is the entire reason hierarchical pooling is needed.

**Implementation skeleton (production).**
```python
# src/mammal_repurposing/calibration/isotonic.py
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import LeaveOneOut
import numpy as np, pickle, pathlib
from scipy.stats import spearmanr

def fit_isotonic_per_target(raw_pkd, chembl_pchembl, target_uniprot,
                            force_direction=None):
    """
    force_direction: None ('auto'), True, or False.
    Returns: (fitted IsotonicRegression, LOCO-RMSE, LOCO-rho, bootstrap_CI)
    """
    inc = 'auto' if force_direction is None else force_direction
    iso = IsotonicRegression(increasing=inc, out_of_bounds='clip',
                             y_min=2.0, y_max=11.0)
    iso.fit(raw_pkd, chembl_pchembl)

    # Leave-one-compound-out CV
    loo_preds = np.zeros_like(chembl_pchembl, dtype=float)
    for i in range(len(raw_pkd)):
        mask = np.arange(len(raw_pkd)) != i
        iso_i = IsotonicRegression(increasing=inc, out_of_bounds='clip',
                                   y_min=2.0, y_max=11.0)
        iso_i.fit(raw_pkd[mask], chembl_pchembl[mask])
        loo_preds[i] = iso_i.predict([raw_pkd[i]])[0]

    loo_rmse = np.sqrt(np.mean((loo_preds - chembl_pchembl)**2))
    loo_rho, _ = spearmanr(loo_preds, chembl_pchembl)

    # Bootstrap CI on rho
    rng = np.random.default_rng(0)
    rhos = []
    for _ in range(1000):
        idx = rng.integers(0, len(raw_pkd), len(raw_pkd))
        if len(np.unique(idx)) < 4:
            continue
        iso_b = IsotonicRegression(increasing=inc, out_of_bounds='clip',
                                   y_min=2.0, y_max=11.0)
        iso_b.fit(raw_pkd[idx], chembl_pchembl[idx])
        r, _ = spearmanr(iso_b.predict(raw_pkd), chembl_pchembl)
        rhos.append(r)
    ci_lo, ci_hi = np.percentile(rhos, [2.5, 97.5])

    pickle.dump(iso, open(f'data/calibration/isotonic/{target_uniprot}.pkl','wb'))
    return iso, loo_rmse, loo_rho, (ci_lo, ci_hi)
```

#### 1B. Venn-ABERS predictors (Mervin et al. 2020 *J. Chem. Inf. Model.*)

**Theory.** Two isotonic regressions are fit — one treating the test point as positive, one as negative — yielding a probability interval [p0, p1] with the Vovk–Petej (UAI 2014, pp. 829–838) validity guarantee that at least one of (p0, p1) is perfectly marginally calibrated under exchangeability. Mervin et al. (2020) demonstrated on 40M compound-target pairs across 2,112 AstraZeneca targets that "VA achieved the best calibration performances across all machine learning algorithms and cross validation methods tested and also the lowest (best) Brier score loss" — with Naïve Bayes, SVM, and random forest base learners under stratified shuffle split (SSS) and leave-20%-of-scaffolds-out (L20SO) validation.

**Strengths.**
- Statistical validity in finite samples (the principal selling point).
- Returns an interval (p0, p1) → fusion module can propagate uncertainty into the RRF scoring.
- Best-in-class Brier per the largest published bioactivity calibration benchmark.

**Weaknesses.**
- Binary classification primitive — for our continuous-pKd regression we must either discretize at pchembl=6 (1 µM) and run binary, or use the experimental IVAR variant (Nouretdinov 2018 PMLR 91:1–22).
- Compute is ~10× isotonic (two PAVA fits + interval book-keeping per test point).
- Interval-valued predictions complicate downstream RRF — we'd collapse the interval to a midpoint for ranking and ship the width as a confidence column.

**Implementation (binary, recommended for v1).**
```python
# src/mammal_repurposing/calibration/venn_abers.py
# pip install venn-abers   (ip200/venn-abers, v1.4.6, MIT)
from venn_abers import VennAbersCalibrator
import numpy as np

def fit_va_binary(raw_pkd, chembl_pchembl, target_uniprot, thresh=6.0):
    y_binary = (chembl_pchembl >= thresh).astype(int)
    va = VennAbersCalibrator(inductive=True, cal_size=0.5, random_state=42)
    # Manual mode passes raw scores directly as the underlying classifier output
    va.fit_manual(raw_pkd.reshape(-1,1), y_binary)
    p0, p1 = va.predict_proba(raw_pkd.reshape(-1,1))  # interval
    return va, p0, p1
```

**Decision**: ship Venn-ABERS as a *diagnostic* output (interval width as a per-compound confidence flag), but **not** as the production point-estimate calibrator — the binary discretization throws away too much rank information in the saturated tropane cluster.

#### 1C. Beta-calibration (Kull, Silva Filho & Flach 2017 AISTATS)

**Theory.** A 3-parameter parametric map g(s) = 1/(1 + 1/(exp(c) · sᵃ · (1−s)⁻ᵇ)). Includes the identity (a=b=1, c=0), so it never uncalibrates a well-calibrated input — Kull et al.'s explicit motivation against Platt. The 2-parameter reductions (`am`: b=a; `ab`: c=0) are useful when a/b are unstable at small n.

**Strengths.**
- Smooth and easier to extrapolate than isotonic.
- Robust at n ≈ 10–30.
- One-line API via `betacal.BetaCalibration(parameters="abm")` (Silva Filho & Perello Nieto, PyPI v1.1.0, MIT).

**Weaknesses.**
- Inputs must be in [0,1] — we rescale raw_pkd to [0,1] via min-max on the calibration set.
- Cannot represent strict sign-flip cleanly: with `parameters="abm"` a → negative will produce a decreasing map, but the parametric family is most natural for monotone-increasing distortions, not full reversals.
- Doesn't naturally extend to hierarchical pooling.

**Implementation.**
```python
# src/mammal_repurposing/calibration/beta_cal.py
from betacal import BetaCalibration   # pip install betacal
import numpy as np

def fit_beta_per_target(raw_pkd, chembl_pchembl, target_uniprot,
                       parameters="abm"):
    lo, hi = raw_pkd.min()-0.1, raw_pkd.max()+0.1
    s = (raw_pkd - lo) / (hi - lo)
    y = (chembl_pchembl - 2.0) / 9.0  # pchembl 2..11 → 0..1
    bc = BetaCalibration(parameters=parameters)
    bc.fit(s.reshape(-1,1), y)
    return bc, (lo, hi)
```

**Decision**: deploy beta-calibration at HCRTR1 (n=6) and as the fallback at GRIN2A (n=8) when hierarchical Bayesian fails to converge.

#### 1D. Per-target deployment decision matrix

| n bucket  | Sign of empirical ρ | Recommended calibrator                          | Rationale                                                                 |
|-----------|--------------------|--------------------------------------------------|--------------------------------------------------------------------------|
| n ≥ 25    | any                | **classical isotonic, `increasing='auto'`**     | PAVA has density; auto-flip is stable                                    |
| 15 ≤ n < 25 | \|ρ\| > 0.4       | **classical isotonic, force direction**         | direction is unstable in 'auto', force via sign of empirical ρ           |
| 15 ≤ n < 25 | \|ρ\| ≤ 0.4       | **hierarchical Bayesian (pooled with family)**  | independent isotonic too noisy; borrow strength                          |
| 8 ≤ n < 15 | family available   | **hierarchical Bayesian (pooled with family)**  | Neelon-Dunson with strong family hyperprior                              |
| 8 ≤ n < 15 | no family pool     | **beta-calibration, `parameters="ab"`**         | 2-parameter parametric — minimal overfit risk                            |
| n < 8     | any                | **NONE — keep uncalibrated at weight 0.30**     | below any defensible threshold; route to §7.7 ensemble                   |

### 2. Hierarchical Bayesian Isotonic for the Small-n + Family-Level Problem

**Concept.** Replace 4 independent IsotonicRegression fits at SLC6A3/SLC6A2/GRIN2B/GRIN2A with a single PyMC model that:
- Shares a family-level hyperprior on monotone slopes across the SLC6 transporter pair (DAT, NET).
- Shares a separate family-level hyperprior on monotone slopes across the GRIN subunit pair (GluN2A, GluN2B).
- Allows family-level prior on sign: SLC6 family slope_mu has prior centered at −1 (sign-flipped), GRIN family slope_mu has prior centered at 0 (agnostic).
- Each target has its own knot positions and per-knot slope draws, partially pooled toward the family mean.

**Why partial pooling helps**:
- At GRIN2A (n=8) the within-target signal is essentially absent (ρ=−0.35, Fisher-transformed 95% CI on Spearman ≈ [−0.74, +0.36] spans 0). Independent isotonic is just memorization. Partial pooling with GRIN2B (n=14) gives the GRIN2A fit access to 22 effective calibration points worth of family-level shape.
- At SLC6A2 (n=25) and SLC6A3 (n=26), independent isotonic is probably already adequate (Section 1D recommends classical isotonic in this n bucket); the hierarchical pool acts mainly as a coherence regularizer that prevents the two SLC6 calibrators from disagreeing about the sign.

**Recommended implementation path** (PyMC, 3–5 day budget):

```python
# src/mammal_repurposing/calibration/hierarchical.py
import pymc as pm
import numpy as np
import pytensor.tensor as pt

def hierarchical_monotone_family(family_data, family_name, n_knots=5,
                                 prior_slope_mu=0.0, prior_slope_sigma=2.0):
    """
    family_data: dict of {target_uniprot: (raw_pkd, chembl_pchembl)}
    family_name: 'SLC6' or 'GRIN'
    Returns: pm.Model with posterior samples
    """
    targets = list(family_data.keys())

    with pm.Model() as model:
        # Family-level hyperprior on the global slope direction & magnitude
        family_slope_mu = pm.Normal('family_slope_mu',
                                    mu=prior_slope_mu, sigma=prior_slope_sigma)
        family_slope_sd = pm.HalfNormal('family_slope_sd', sigma=1.0)
        family_intercept = pm.Normal('family_intercept', mu=6.0, sigma=2.0)

        # Austin Rochford monotone prior: mo(j) = b · cumsum(Dirichlet(1,...,1))
        # adapted into the Neelon–Dunson 2004 piecewise-linear isotonic framework
        for t in targets:
            raw, y = family_data[t]
            x = (raw - raw.min()) / (raw.max() - raw.min() + 1e-9)
            knot_idx = np.clip((x * n_knots).astype(int), 0, n_knots-1)

            b_t   = pm.Normal(f'b_{t}', mu=family_slope_mu, sigma=family_slope_sd)
            xi_t  = pm.Dirichlet(f'xi_{t}', a=np.ones(n_knots))
            mo_t  = pm.Deterministic(f'mo_{t}', b_t * pt.cumsum(xi_t))
            sig_t = pm.HalfNormal(f'sigma_{t}', sigma=1.0)
            mu_t  = family_intercept + mo_t[knot_idx]
            pm.Normal(f'y_{t}', mu=mu_t, sigma=sig_t, observed=y)

        trace = pm.sample(2000, tune=1000, chains=4, target_accept=0.95,
                          random_seed=42)
    return model, trace
```

**Theoretical anchor**: this is the Austin Rochford PyMC monotone prior (b · cumsum(Dirichlet(α=1))) wrapped in the Neelon & Dunson (2004) *Biometrics* 60(2):398–406 piecewise-linear isotonic framework, with the Lin & Dunson (2014) *Biometrika* 101(2):303–317 empirical-Bayes justification for borrowing across related curves via family-level hyperpriors. The `b_t` per-target slope is the partial-pooling lever; setting `family_slope_sd` very small forces full pooling, very large recovers independent fits.

**Family-prior calibration**:
- SLC6 family: `prior_slope_mu = -1.0, prior_slope_sigma = 0.5` (informative — both targets show empirical sign flip)
- GRIN family: `prior_slope_mu = 0.0, prior_slope_sigma = 1.5` (agnostic — GRIN2A and GRIN2B may have different failure mechanisms)

**Empirical-Bayes shortcut** (if MCMC budget is tight): fit isotonic per-target independently, compute the empirical mean and SD of the per-target step heights, then refit each target with that empirical mean as a Gaussian prior on the step heights via penalized PAVA. This is ~3 lines of scipy and gives ~80% of the hierarchical benefit at 1% of the runtime. Implement this as a fallback if PyMC convergence fails (Rhat > 1.05 or BFMI < 0.3).

**Validation**: pool the 4 INVERTED targets, compute LOCO Brier score / RMSE on held-out compounds, and compare:

| Variant | Expected LOCO-RMSE (pchembl units) |
|---|---|
| Uncalibrated MAMMAL | 1.4–1.8 |
| Independent isotonic per-target | 0.9–1.3 at SLC6A3/A2; 1.5–2.0 at GRIN2A/B |
| Hierarchical Bayesian (family-pooled) | 0.9–1.2 at SLC6A3/A2; 1.0–1.4 at GRIN2B; 1.2–1.6 at GRIN2A |

### 3. Integration With the V3 Pipeline

**New module layout** (drop into existing tree):
```
src/mammal_repurposing/calibration/
├── isotonic.py          # sklearn IsotonicRegression + LOCO + bootstrap
├── venn_abers.py        # ip200/venn-abers wrapper, binary-discretised mode
├── beta_cal.py          # betacal.BetaCalibration(parameters="abm") wrapper
├── hierarchical.py      # PyMC hierarchical monotone (SLC6 pool, GRIN pool)
├── router.py            # decision tree, reads weights_calibrated.yaml
└── diagnostics.py       # LOCO-rho, bootstrap CI, Brier, reliability plots

scripts/
├── 29_v3_calibration_comparison.py   # 4-way sweep across all 22 targets
└── 30_v3_isotonic_recalibration.py    # production deploy

data/calibration/
├── isotonic/{UNIPROT}.pkl
├── venn_abers/{UNIPROT}.pkl
├── beta/{UNIPROT}.pkl
├── hierarchical/slc6_pooled.pkl   # joint posterior over DAT+NET
├── hierarchical/grin_pooled.pkl   # joint posterior over GRIN2A+GRIN2B
└── router_decisions.csv           # one row per target with chosen calibrator + LOCO metrics
```

**`weights_calibrated.yaml` schema migration** (additive — old fields remain valid):
```yaml
SLC6A3_HUMAN:
  weight: 0.30
  calibrator_type: isotonic
  calibrator_path: data/calibration/isotonic/SLC6A3_HUMAN.pkl
  calibrator_direction: decreasing
  loco_rho: 0.52
  loco_rmse: 0.83
  bootstrap_ci_rho: [0.31, 0.68]
  fit_quality: high
GRIN2B_HUMAN:
  weight: 0.30
  calibrator_type: hierarchical_isotonic
  calibrator_path: data/calibration/hierarchical/grin_pooled.pkl
  calibrator_target_idx: 1
  loco_rho: 0.22
  loco_rmse: 1.12
  bootstrap_ci_rho: [-0.05, 0.45]   # spans 0 — flag low-confidence
  fit_quality: low
GRIN2A_HUMAN:
  weight: 0.0
  calibrator_type: none
  fit_quality: insufficient_n
  escalation_target: mmatt_dta
```

**Inference path** (in `fusion.py`):
```python
def apply_per_target_calibration(mammal_raw_pkd, target_uniprot, cfg):
    spec = cfg[target_uniprot]
    if spec['calibrator_type'] == 'none':
        return None  # caller drops MAMMAL contribution
    cal = pickle.load(open(spec['calibrator_path'],'rb'))
    if spec['calibrator_type'] == 'isotonic':
        return cal.predict([mammal_raw_pkd])[0]
    elif spec['calibrator_type'] == 'beta':
        s = (mammal_raw_pkd - spec['rescale_lo']) / (spec['rescale_hi'] - spec['rescale_lo'])
        return cal.predict(np.array([[s]]))[0] * 9.0 + 2.0
    elif spec['calibrator_type'] == 'hierarchical_isotonic':
        return posterior_predict(cal, mammal_raw_pkd, spec['calibrator_target_idx'])
```

**Backward compatibility**: `--use-isotonic-calibration` flag on `scripts/run_fusion.py`; default ON once `scripts/29` validates. Old uncalibrated outputs land in `outputs/uncalibrated/v3_*.parquet`; calibrated outputs in `outputs/calibrated/v3_*.parquet`.

**Estimated wall-clock**:
- Classical isotonic LOCO + bootstrap per target: ~3 s
- All 22 targets sweep: ~70 s
- Beta-calibration per target: ~2 s
- Venn-ABERS per target: ~8 s
- Hierarchical PyMC (SLC6 pool, 2000 draws × 4 chains × NUTS): ~3 min
- Hierarchical PyMC (GRIN pool): ~3 min
- Total deploy time end-to-end: **< 12 min**

### 4. Validation Gates and Expected Outcomes

**Five-gate validation per (calibrator, target) pair.**

1. **LOCO RMSE on calibrated pchembl** must be ≤ raw MAMMAL RMSE. If a calibrator makes RMSE *worse*, drop it.
2. **LOCO Spearman ρ on held-out**: must improve from the pre-cal value AND must clear the tier gate.
3. **Bootstrap 95% CI on LOCO ρ** (1,000 resamples): CI must not span 0.
4. **Sign-correction stability**: across 1,000 bootstrap PAVA fits, the `increasing='auto'` flag must agree in ≥ 95% of resamples. Below this, force the direction from the empirical ρ.
5. **Non-degradation at positive controls** (DRD1, HCRTR1): running the calibrator sweep must not *worsen* ρ at the MAMMAL_ONLY_STRONG targets.

**Predicted per-target outcomes — the ship/no-ship table for §7.11:**

| Target  | n  | Pre-cal ρ | Best calibrator (predicted)              | Expected post-cal ρ | Ship? |
|---------|----|-----------|-----------------------------------------|--------------------|-------|
| SLC6A3  | 26 | −0.71     | Classical isotonic (force decreasing)   | +0.45 to +0.65     | YES, primary |
| SLC6A2  | 25 | −0.53     | Classical isotonic (force decreasing)   | +0.30 to +0.55     | YES, primary |
| GRIN2B  | 14 | −0.30     | Hierarchical Bayesian (GRIN pool)       | +0.10 to +0.35     | SHIP + §7.7 secondary |
| GRIN2A  | 8  | −0.35     | NONE — keep uncalibrated, route to §7.7 | unstable           | NO — escalate |
| DRD1    | 21 | +0.31     | Classical isotonic (auto → increasing)  | +0.35 to +0.50     | YES |
| HCRTR1  | 6  | +0.37     | Beta-calibration `parameters="ab"`      | +0.20 to +0.45     | MARGINAL — flag low-conf |

**Calibration failure handling**: if at any target the bootstrap 95% CI on post-cal ρ spans 0, route to §7.7 cross-DTI ensemble (Section 5).

### 5. Escalation Framework to §7.7 Cross-DTI Ensemble

**Pre-committed numerical thresholds** (these are gate values for §7.11 success):

```
Tier A — ISOTONIC SUFFICIENT (ship as primary):
    post-cal LOCO ρ ≥ +0.40 AND bootstrap 95% CI lower bound > 0

Tier B — ISOTONIC PRIMARY + §7.7 SECONDARY (cross-validate):
    +0.20 ≤ post-cal LOCO ρ < +0.40 AND bootstrap 95% CI lower bound > 0
    Add MMAtt-DTA as second ranker into RRF for this target.

Tier C — ESCALATE IMMEDIATELY (do not ship MAMMAL at this target):
    post-cal LOCO ρ < +0.20 OR bootstrap 95% CI spans 0
    Replace MAMMAL contribution with:
      - MMAtt-DTA (Schulman 2024, Bioinformatics 40(8):btae496) for transporters
      - PSICHIC (Koh 2024, Nat. Mach. Intell. 6:673–687) for GPCRs
      - BALM (Gorantla 2025, J. Chem. Inf. Model. 65(22):12279–12291,
        doi:10.1021/acs.jcim.5c02063) as a fallback for allosteric / ATD targets

Tier D — PANEL-DEPRECATE:
    target n < 8 AND no family pool available
    Drop from the cognition panel entirely; document in router_decisions.csv
```

**MMAtt-DTA fitness for SLC6 transporters**: Schulman et al. (2024) describe MMAtt-DTA as an "attention-based method to predict drug–target bioactivities across human proteins within seven superfamilies", having "examined nine different descriptor sets to identify optimal signature descriptors", with the result that "testing results demonstrated Spearman correlations exceeding 0.72 (P < 0.001) for six out of seven superfamilies." Independent ChEMBL-V33 holdout validation gave "Spearman correlation > 0.57 (P < 0.001) for most superfamilies." The transporter class is one of the seven superfamilies. MMAtt-DTA is therefore the *natural* §7.7 fallback for SLC6A3/SLC6A2 if isotonic underperforms.

**PSICHIC fitness for GRIN2B**: PSICHIC is sequence-only but uses a physicochemical-graph attention that "decode[s] interaction fingerprints directly from sequence data alone" — it could potentially recover ifenprodil-class SAR if its training set spans ATD-interface pharmacology, but the GluN1/GluN2B heterodimer-specific signal is still missing because PSICHIC also takes single-chain sequences. **Honest caveat**: even §7.7 cannot fully fix GRIN2B without explicit dimer-aware input (e.g., concatenated GluN1+GluN2B sequence as a heterodimer pseudo-input, or AlphaFold-Multimer/Boltz-2-derived contact features).

**BALM fitness as a general fallback**: BALM (Gorantla et al. 2025 *J. Chem. Inf. Model.* 65(22):12279–12291, preprint bioRxiv 2024.11.01.621495) "leverages pre-trained language models to encode protein sequences and ligand SMILES strings" — ESM-2 for proteins, ChemBERTa-2 for ligands, with reparameterized PEFT methods (LoRA, LoHa, LoKr) and the additive IA3 adapter applied to key/query/value matrices. The cosine-similarity-in-latent-space objective optimizes pKd directly. This is the deepest-tunable backup if both MAMMAL and MMAtt-DTA fail.

**Predicted routing for the 4 INVERTED targets**:

| Target  | After §7.11 calibration | Predicted Tier | Action |
|---------|--------------------------|----------------|--------|
| SLC6A3  | post-cal ρ ≈ +0.55       | A              | Ship isotonic; no §7.7 |
| SLC6A2  | post-cal ρ ≈ +0.40       | A/B boundary   | Ship isotonic; consider §7.7 MMAtt-DTA secondary |
| GRIN2B  | post-cal ρ ≈ +0.20       | B/C boundary   | Ship hierarchical Bayesian; add PSICHIC/MMAtt-DTA via RRF |
| GRIN2A  | post-cal ρ ≈ +0.10 ± wide CI | C/D        | Keep uncalibrated, route to §7.7, possibly panel-deprecate |

### 6. Lateral Considerations the User Hasn't Asked About

**6.1 Catastrophic forgetting at MAMMAL_ONLY_STRONG targets.** Deploying isotonic at INVERTED targets is per-target, so DRD1 and HCRTR1 are not directly modified. BUT — the router must still process them. With `increasing='auto'` on DRD1 (n=21, ρ=+0.31), PAVA will correctly pick increasing, and post-cal ρ should be ≥ +0.31. HCRTR1 (n=6, ρ=+0.37) is below our isotonic-safe threshold; the router should select beta-calibration `parameters="ab"` and may very slightly degrade ρ (the parametric form is smoother). Net risk: marginal worsening at HCRTR1 only. Test gate: post-cal LOCO ρ at DRD1/HCRTR1 must be within 0.05 of pre-cal.

**6.2 Selectivity-vector (Gini) implications.** §7.4 computes a Gini coefficient on the calibrated affinity vector across the 22-target panel. Isotonic-transforming 4 dimensions (SLC6A3, SLC6A2, GRIN2B, GRIN2A) will change the vector geometry. Expected effects:
- SLC6A3 and SLC6A2 currently anti-correlate with truth → after sign flip, true selective DAT/NET ligands (solriamfetol, methylphenidate) will have their MAMMAL-derived selectivity *toward* DAT/NET converted from anti-correlation to correlation. Their Gini scores will *increase* (more selective-looking), correctly.
- GRIN2B with weakened post-cal ρ (~+0.2) → its dimension carries less signal; the Gini contribution from this axis becomes more nearly-random.
- GRIN2A zeroed out → dimension removed; effective panel dimensionality drops to 21 (or 18 if all four INVERTED are escalated). Gini denominators must be re-normalized.

Predicted shift in top-5 of each §8.1 facet shortlist: ~1–2 reordering events per facet, dominated by DAT/NET-targeted compounds rising in monoamine-selective facets and ifenprodil-like compounds (e.g., traxoprodil) shifting in NMDA-selective facets.

**6.3 Calibration-driven ranking inversions — concrete predictions.**
- *Solriamfetol* (DAT Ki = 14,200 nM → pKi ≈ 4.85; NET Ki = 3,700 nM → pKi ≈ 5.43; FDA NDA 211230; Baladi et al. 2018 *J. Pharmacol. Exp. Ther.* 366(2):367): currently ranked in the top-25 partly via inflated MAMMAL SLC6A3 score. After isotonic with `increasing=False`, MAMMAL's high raw pKd at DAT will map to a *low* calibrated pchembl (correctly, since solriamfetol is a low-µM DAT ligand). Expected rank drop: from top-25 to ~rank 40–60 in the DAT facet, but potentially *retained* in the NET facet because Ki = 3.7 µM is the higher-affinity arm.
- *Methylphenidate* (DAT Ki ≈ 100 nM, pKi ≈ 7.0): currently may rank high or middling depending on MAMMAL's raw output. After isotonic, its high empirical pKi should map to a *high* calibrated pchembl — expected rank stays high or *rises* in the DAT facet. This is the positive control: if methylphenidate falls in post-cal ranking, the calibrator is broken.

**6.4 Provenance metadata surfacing.** Each fitted calibrator carries a `loco_rmse` and `bootstrap_ci_rho`. Surface these in:
- `data/calibration/router_decisions.csv` (one row per target)
- The shortlist parquet (`outputs/calibrated/v3_*.parquet`) — add `mammal_calibration_quality` and `mammal_calibration_ci_width` columns
- The wet-lab handoff PDF (§9 deliverable) — annotate each predicted score with "calibrated via [isotonic|hierarchical|beta|none] (LOCO RMSE = X.XX pchembl)"

**6.5 Covariate shift between MAMMAL prediction inputs and ChEMBL ground truth.** This is the deepest risk. ChEMBL pchembl aggregates across assays (Ki, IC50, Kd) and conditions (buffer, temperature, expression system, species). Landrum & Riniker (2024) *J. Chem. Inf. Model.* 64(5):1560–1567 (doi:10.1021/acs.jcim.4c00049) quantified this directly for IC50 replicates: with minimal curation, "almost 65% of the points differ by more than 0.3 log units, 27% differ by more than one log unit, and the correlation between the assays, as measured by Kendall's τ, is only 0.51"; with maximal curation, "48% differ by more than 0.3 log units, 13% by more than one log unit, Kendall's τ = 0.71." Our isotonic fit assumes MAMMAL's predicted pKd has a stable monotone relationship to a single noisy aggregate. Three mitigations:
1. **Stratify the calibration set by assay type** if n permits: fit separate isotonic for Ki vs IC50 measurements where n ≥ 15 per stratum. At SLC6A3/SLC6A2 this is feasible; at GRIN targets it is not.
2. **Inflate the bootstrap CI** by adding a noise floor (use the Landrum–Riniker maximal-curation 0.3-log unit benchmark) to LOCO RMSE comparisons. A calibrator whose LOCO RMSE is below this floor is suspicious (likely overfit, since it has out-resolved the underlying assay-to-assay disagreement).
3. **Document the assay distribution** of each per-target calibration set in `router_decisions.csv` (% Ki / % IC50 / % Kd, % radioligand-binding vs functional). Downstream users (wet-lab) can then judge applicability domain.

### 7. Publication Angle

**Working title.** "Per-target monotone post-hoc recalibration of foundation-model DTI predictions: a comparative implementation of isotonic, Venn-ABERS, beta-calibration, and hierarchical Bayesian variants on the MAMMAL cognition panel."

**Venue priority order**:
1. *J. Cheminform.* (BMC, open-access, methods-focused — natural home given the calibration + DTI methodology audience)
2. *Bioinformatics* (Oxford, broader audience, but tighter scope match)
3. arXiv q-bio.QM preprint first as the anchor (1-week timeline once §7.11 results are in)

**Novelty claims defensible against reviewers**:
1. **First systematic 4-way comparison** (isotonic / Venn-ABERS / beta / hierarchical Bayesian) on a foundation-model DTI head where specific targets exhibit documented anti-correlated outputs — Mervin et al. 2020 compared three calibrators but only on classical ML base models (NB/SVM/RF), not on a foundation model with the failure modes we've documented.
2. **First DTI-applied hierarchical Bayesian isotonic calibration** with explicit family-level hyperpriors over related protein subfamilies (SLC6 monoamine transporters; iGluR NMDA subunits). Neelon & Dunson 2004 and Lin & Dunson 2014 provide the statistical machinery but never applied to ligand-target calibration.
3. **A pre-committed decision tree for per-target calibrator routing** with explicit n-bucket thresholds (Section 1D) — no published guidance exists on when to switch between isotonic / beta / hierarchical for chemogenomic calibration.
4. **The Tier-A/B/C/D escalation framework to cross-DTI ensembling** (Section 5) — explicit numerical gates connecting calibration outcomes to ensemble routing.

**Reproducibility deliverables (MIT-licensed)**:
- All 4 calibrator wrappers + the PyMC hierarchical model + the router code
- Per-target pickles + bootstrap LOCO outputs (supplementary)
- The Phase A.7 calibration report (already exists as the empirical anchor)
- A standalone notebook reproducing every plot from the raw MAMMAL + ChEMBL extracts

**Honest weaknesses to pre-empt**:
- n = 22 targets is small for sweeping conclusions about which calibrator wins generically; the paper is honest about being a cognition-panel study, not a general DTI benchmark.
- GRIN2B's underlying problem (ATD-dimer-interface pharmacology) is not solved by calibration; this is a *negative result* and should be reported as such with a structural-biology explanation per Karakas et al. 2011.
- Hierarchical Bayesian fits require modest MCMC budget (~6 min total); for groups without PyMC infrastructure we provide the empirical-Bayes shortcut as a graceful fallback.

---

## Recommendations

**Immediate (next 3–5 days, in order):**

1. **Day 1 — instrument**: write `scripts/29_v3_calibration_comparison.py` that runs all four calibrators on all 22 targets, emits `router_decisions.csv` with LOCO ρ, LOCO RMSE, bootstrap 95% CI, and reliability plot for each (calibrator, target) pair. Deliverable: a 22×4 results matrix.
2. **Day 2 — fit hierarchical**: implement `calibration/hierarchical.py` with the PyMC model in Section 2. Fit SLC6 pool and GRIN pool. Confirm Rhat < 1.05, BFMI > 0.3. Deliverable: the two pooled posterior pickles + posterior-predictive checks.
3. **Day 3 — route**: implement `calibration/router.py` per the decision tree in Section 1D, populate `weights_calibrated.yaml` with `calibrator_type` per target. Deliverable: production calibration config.
4. **Day 4 — ship**: write `scripts/30_v3_isotonic_recalibration.py`, regenerate the full V3 calibrated parquets, regenerate the §7.4 selectivity vectors, regenerate the §8.1 faceted shortlist. Deliverable: a "before / after §7.11" diff report on the top-25 compounds in each facet.
5. **Day 5 — validate against named controls**: confirm methylphenidate stays high in DAT facet (sanity); confirm solriamfetol falls in DAT facet (correct rescue); confirm traxoprodil/ifenprodil behave reasonably in GRIN2B (or flag for §7.7); confirm DRD1 (positive control) post-cal ρ within ±0.05 of pre-cal.

**Thresholds that change the recommendation**:
- If LOCO ρ at SLC6A3 < +0.20 after isotonic → the failure is not Scenario-2 rank-collapse; route to §7.7 immediately and re-examine the assumption.
- If hierarchical Bayesian's posterior at GRIN2B has 95% CI on the slope spanning 0 → the GRIN family pool isn't borrowing strength; deprecate hierarchical and ship uncalibrated + §7.7.
- If catastrophic forgetting at DRD1 (post-cal ρ degrades by > 0.10) → the `increasing='auto'` logic is misbehaving; force `increasing=True` at all MAMMAL_ONLY_STRONG targets.

**Medium-term (next 2 weeks)**: spec and prototype the §7.7 cross-DTI ensemble (MMAtt-DTA + PSICHIC + BALM) for the 1–2 targets that §7.11 fails to rescue. Use the Tier-C escalation list from `router_decisions.csv` as the input scope.

**Long-term (preprint)**: arXiv submission within 2 weeks of §7.11 completion; *J. Cheminform.* full submission within 6 weeks.

---

## Caveats

1. **Authorship correction**: the task brief attributes the canonical large-scale calibration benchmark to "Toplak 2020" with PMID 32865408. The actual citation is **Mervin, Afzal, Engkvist & Bender (2020)** *J. Chem. Inf. Model.* 60(10):4546–4559, doi:10.1021/acs.jcim.0c00476. This is the correct reference for the Venn-ABERS / isotonic / Platt comparison on 40M compound-target pairs across 2,112 AstraZeneca targets. The methodological conclusions in the brief (Venn-ABERS wins by Brier; isotonic overfits at small n) remain accurate. There is a separate "Toplak" paper in *J. Chem. Inf. Model.* (Toplak et al. 2014 on assay-interference compounds) which is unrelated to calibration. Use Mervin et al. 2020 in citations.

2. **Venn-ABERS regression is experimental**: the production-grade `venn-abers` Python package (ip200/venn-abers v1.4.6, MIT, author Ivan Petej) primarily supports binary and multi-class classification ("Python implementation of Venn-ABERS calibration for binary and multiclass classification problems"). The regression variant (IVAR per Nouretdinov 2018) is shipped as an example notebook, not as the headline API. For v1 of §7.11 we therefore deploy VA only in a binary-discretized mode (threshold pchembl=6, equivalent to 1 µM activity cut). Continuous IVAR is parked for v2.

3. **Hierarchical Bayesian convergence is not guaranteed at n=8**: GRIN2A with n=8 calibration points may cause the GRIN-pool MCMC to be dominated by the GRIN2B partial pool. This is *fine* for GRIN2B (its fit is improved) but the per-target inference at GRIN2A becomes essentially "what would the prior say if you saw 8 noisy points." If `pm.sample` returns BFMI < 0.3 or divergent transitions > 10%, the empirical-Bayes shortcut (Section 2) is the documented fallback.

4. **GRIN2B's underlying problem is structural, not statistical**: ifenprodil-class antagonism is a dimer-interface phenomenon (Karakas et al. 2011 *Nature* 475:249–253). MAMMAL's single-chain protein input *cannot* see it. Even a perfect post-hoc calibrator can only reorder the off-pathway signal MAMMAL has — there is no monotone transformation of a single-chain embedding that recovers heterodimer pharmacology. We expect a post-cal ρ ceiling around +0.30 at GRIN2B; if §7.11 hits that ceiling, §7.7 is the architectural fix, not a better calibrator.

5. **Assay heterogeneity bounds achievable RMSE**: Landrum & Riniker (2024) *J. Chem. Inf. Model.* 64(5):1560–1567 documented the IC50-to-IC50 noise floor in aggregated ChEMBL bioactivities — even with maximal curation, "48% [of replicate pairs] differ by more than 0.3 log units, 13% by more than one log unit," with Kendall's τ = 0.71. A LOCO RMSE below 0.3 pchembl units should be treated with suspicion (the calibrator has out-resolved the underlying assay-to-assay disagreement and is almost certainly overfit). Realistic post-cal RMSE target: 0.7–1.2 pchembl units.

6. **n=22 targets is a small benchmark**: per-target results are noisy. Conclusions about which calibrator wins *generically* across DTI projects should not be drawn from §7.11 alone. The paper-level claim is "here is a decision framework with pre-committed thresholds and reproducible code"; not "isotonic always wins."

7. **`increasing='auto'` instability below n=15**: sklearn's automatic direction detection is itself a model-selection step that can flip sign under bootstrap. At GRIN2A (n=8) and HCRTR1 (n=6) the auto flag is essentially random; we explicitly force the direction from the empirical Spearman sign in our router. Document this in `router_decisions.csv` as `calibrator_direction: forced_from_empirical_rho`.

8. **The MAMMAL_ONLY_INVERTED → MAMMAL_ONLY_RESCUED transition is the empirical test for §7.11 success**: after deploying §7.11, the calibration categorization in the V3 status doc should change for SLC6A3 and SLC6A2 from INVERTED to RESCUED (ρ > +0.30, CI > 0); GRIN2B may transition to WEAK; GRIN2A should remain INSUFFICIENT_N. If those state transitions don't happen, §7.11 has failed and we ship the uncalibrated baseline with §7.7 as the next deliverable.