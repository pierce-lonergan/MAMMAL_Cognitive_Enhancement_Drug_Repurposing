# MH3 — Per-Cell-Line Random Effect on V8 Joint Posterior

## Deep Research Note: Defending the U2OS-to-Brain Transfer Claim

**Status**: research direction synthesis, prelude to V8 OSF §7 amendment
**Companion to**: `GAPS_AND_RESEARCH_DIRECTIONS.md` (MH3 entry), `v8_osf_preregistration.md` §7, V8 Methods draft
**Scope**: methodology + implementation + empirical anchoring + pre-registered validation criteria
**Effort estimate (revised from 1 week)**: 2 weeks for honest version; 3 days for cosmetic version

---

## 1. The actual problem, stated honestly

The V8 πphen pipeline trains chemCPA on LINCS L1000 (multiple cell lines, transcriptomic) and consumes JUMP-CP `cpg0016` (U2OS osteosarcoma, morphological) as the dominant phenotypic anchor. The (L, L, H) discovery score is then asserted to triage compounds for **brain-relevant cognitive effects** (clemastine territory: remyelination, neuronal plasticity, synaptic function).

The skeptical reviewer's one-line objection is: *"You're inferring CNS pharmacology from a bone tumor cell line."* And U2OS is, in fact, derived from a moderately differentiated osteosarcoma of the tibia — an epithelial-like mesenchymal cancer cell, nothing remotely neuronal.

That reviewer is not wrong to ask. The question is whether the V8 framework can give a **principled, quantitative** answer rather than a wave-of-the-hand "well, mitochondria are mitochondria."

The current V8 joint posterior uses MOFA+ ARD for per-modality variance attribution. ARD tells you "factor *k* loads on modality *m* with weight *w*." It does **not** decompose perturbation effects into a transferable cross-context component plus a cell-line-idiosyncratic component. Without that decomposition, the model has no formal expression for "what fraction of the U2OS signal we should trust to generalize."

MH3 fixes that. The fix is a textbook hierarchical Bayesian construct — a per-cell-line random effect on the perturbation latent — but the **value** is not the technique itself. The value is that the resulting posterior gives you (a) a quantitative ICC for cell-line-attributable variance, (b) a posterior contrast between the cross-context and cell-specific perturbation effects, and (c) calibrated predictive intervals that explicitly include cell-line transfer uncertainty. All three together let you write a defensible Methods + Discussion paragraph instead of a defensive one.

---

## 2. Why MOFA+ ARD alone is structurally insufficient

A quick disambiguation, because this is easy to muddle.

MOFA+ (the V8 backbone) is a Bayesian Group Factor Analysis model that infers latent factors *Z* such that each modality *m*'s data matrix *Y^m* ≈ *Z* *W^m* + ε, with ARD priors on *W^m* shrinking factors that don't load on modality *m* to zero. It supports a "groups" axis — you can declare that samples partition into groups (cell line A vs cell line B vs cell line C) and MOFA+ will then estimate per-(factor, group) variance to identify factors that are differentially active across groups. MOFA+ provides a probabilistic framework for comprehensive integration of structured single-cell data with multiple sample groups and data modalities, using sparsity priors and hierarchical variance regularization.

That sounds like it should already solve MH3. It doesn't, for three reasons:

**Reason 1 — Groups model variance, not transferability.** The MOFA+ group axis decomposes *which factors are active in which group*. It does not produce a per-compound posterior on a transferable effect vs. a cell-specific effect. If a compound modulates factor 7, MOFA+ tells you "factor 7 is active in U2OS but not A549" — which is informative for QC but doesn't give you a posterior on `β_compound_transferable`.

**Reason 2 — ARD shrinks weights, not random effects.** ARD is a hierarchical prior on the factor loading matrix *W*. It performs feature selection at the factor-modality interface, not partial pooling of compound effects across cellular contexts. A factor can have a large loading in U2OS and a small one in iPSC neurons without that fact ever entering the V8 (L, L, H) discovery score.

**Reason 3 — The V8 phenotype anchor is U2OS-only.** Because `cpg0016` is U2OS-exclusively, there is no other Cell-Painting cell line in the joint inference for MOFA+ to actually contrast against. The only multi-cell-line modality is LINCS L1000 (and partly the mouse iPSC-MEA data). The transfer claim, by construction, asks the model to assert that signal extracted from a single-cell-line phenotypic modality is interpretable in another tissue context. That requires a model term that **explicitly represents cellular context as a source of variance in the compound-level effect**, not a model that just absorbs context into the latent factors.

This is why MH3 lives as a separate item from V8's existing MOFA+ ARD machinery. It's not a redundant variance term — it's a different mathematical object operating on a different parameter (the compound effect, not the factor loading).

---

## 3. Mathematical specification

The cleanest construction, which threads through MH3 + MH7 (species) coherently, is a partially crossed random-effects model on the per-compound latent perturbation embedding. Let:

- *c* ∈ {1, …, *C*} index compounds
- *l* ∈ {1, …, *L*} index cell lines (U2OS, A549, iPSC-cortical, iPSC-mDA, …)
- *s* ∈ {human, mouse} index species
- *m* ∈ {morphology, transcriptomic, electrophysiology, …} index modalities
- *k* ∈ {1, …, *K*} index endpoints (the V8 (L, L, H) bands or the chemCPA-decoded transcriptomic distance)

The observed phenotypic signal for compound *c* in cell line *l* of species *s*, on modality *m*, endpoint *k*, replicate *r* is:

```
y_{c,l,s,m,k,r} = μ_m,k + β_{c,k} + α_{l,k} + γ_{s,k} + δ_{c,l,k} + ε_{c,l,s,m,k,r}
```

Where:

| Term | Interpretation | Prior |
|------|----------------|-------|
| `μ_m,k` | global intercept for (modality, endpoint) | weakly informative Normal |
| `β_{c,k}` | **transferable compound effect** — the quantity V8 (L, L, H) ultimately wants to score | `β_{c,k} ~ Normal(0, σ_β,k²)` with `σ_β,k ~ HalfNormal` |
| `α_{l,k}` | **per-cell-line random effect (MH3)** | `α_{l,k} ~ Normal(0, σ_α,k²)`, `σ_α,k ~ HalfNormal(0.5)` |
| `γ_{s,k}` | **per-species random effect (MH7)** | `γ_{s,k} ~ Normal(0, σ_γ,k²)`, `σ_γ,k ~ HalfNormal(0.3)` |
| `δ_{c,l,k}` | **compound × cell-line interaction** — the cell-specific portion of the compound effect that does not transfer | `δ_{c,l,k} ~ Normal(0, σ_δ,k²)` |
| `ε` | residual noise | `Normal(0, σ_ε,k²)` |

The interaction term `δ_{c,l,k}` is the methodologically important addition. Without it, MH3 reduces to a "different cell lines have different baselines" model, which is sociologically interesting but does not defend the transfer claim. The transfer claim is fundamentally a claim about whether the compound effect is dominated by the global β or the local δ. Modeling `δ` explicitly with its own variance component is what makes the transfer claim testable.

The **transferable score** used in V8 downstream (the input to the (L, L, H) gating) is then:

```
β̂_{c,k} = posterior mean of β_{c,k}
ŝ_{c,k} = posterior std of β_{c,k}    # honest predictive uncertainty
```

And the **transferability index** for compound *c* on endpoint *k* is:

```
T_{c,k} = E[ |β_{c,k}| / (|β_{c,k}| + std(δ_{c,·,k})) | data ]
```

Compounds with `T_{c,k}` close to 1 have effects dominated by the transferable component; compounds with `T_{c,k}` close to 0 are essentially U2OS-idiosyncratic. The (L, L, H) gate in V8 should be tightened with a `T_{c,k} > 0.6` filter (or whatever threshold the prior predictive check supports). This is the **single most important deliverable from MH3** for the V8 paper — it converts "U2OS-to-brain transfer" from a rhetorical claim into a per-compound score with a posterior distribution.

### Population-level variance partition

Across compounds, the analogous quantity is the **cell-line intraclass correlation coefficient (ICC)**:

```
ICC_cell,k = σ_α,k² / (σ_β,k² + σ_α,k² + σ_γ,k² + σ_δ,k² + σ_ε,k²)
```

The intraclass correlation coefficient (also called the variance partition coefficient) decomposes total variance into the proportion attributable to systematic differences between clusters versus within-cluster variation, and is calculated directly from the variance estimates of mixed-effects models.

The interpretation matters: **a small ICC_cell is good news for the transfer claim** (most of the variance is in the transferable β term, not in the cell-line term), but the more important quantity is the **interaction ICC**:

```
ICC_int,k = σ_δ,k² / (σ_β,k² + σ_δ,k²)
```

`ICC_int` directly measures the fraction of compound-effect variance that is cell-line-idiosyncratic. If `ICC_int < 0.2` on endpoint *k*, then U2OS-derived β estimates are a good proxy for the transferable effect on that endpoint. If `ICC_int > 0.5`, then U2OS is not a useful surrogate for that endpoint and the V8 paper should explicitly demote those endpoints in the (L, L, H) score.

Either result is **publishable**. The honest version of MH3 is a calibration exercise, not a victory lap.

---

## 4. PyMC implementation (with numpyro backend)

Given that V6.B.5 NUTS already runs on numpyro after the recent sprint, the V8 hierarchical extension shares the same infrastructure. The non-centered parameterization is mandatory here — with five variance components, the centered version will divergence-bomb. When there is insufficient data in a hierarchical model, the variables being inferred end up having correlation effects, making it difficult to sample; one solution is to obtain more data, but when this isn't possible, reparameterization by creating a non-centered model from the centered model is the standard remedy.

Reference implementation, drop-in for the V8 joint posterior:

```python
import pymc as pm
import pytensor.tensor as pt
import numpy as np

def build_v8_hierarchical_with_cell_random_effect(
    y,                  # shape (N,) observed signal
    compound_idx,       # shape (N,) compound index
    cell_idx,           # shape (N,) cell-line index
    species_idx,        # shape (N,) species index
    endpoint_idx,       # shape (N,) endpoint index
    n_compounds, n_cells, n_species, n_endpoints,
):
    coords = {
        "compound": np.arange(n_compounds),
        "cell":     np.arange(n_cells),
        "species":  np.arange(n_species),
        "endpoint": np.arange(n_endpoints),
        "obs":      np.arange(len(y)),
    }
    with pm.Model(coords=coords) as model:
        # ---- Global intercepts per endpoint ----
        mu_k = pm.Normal("mu_k", 0.0, 1.0, dims="endpoint")

        # ---- Hyperpriors on variance components (per endpoint) ----
        sigma_beta  = pm.HalfNormal("sigma_beta",  0.5, dims="endpoint")
        sigma_alpha = pm.HalfNormal("sigma_alpha", 0.5, dims="endpoint")
        sigma_gamma = pm.HalfNormal("sigma_gamma", 0.3, dims="endpoint")
        sigma_delta = pm.HalfNormal("sigma_delta", 0.5, dims="endpoint")
        sigma_eps   = pm.HalfNormal("sigma_eps",   1.0, dims="endpoint")

        # ---- Non-centered random effects ----
        # Transferable compound effect: β_{c,k} = sigma_beta_k * raw
        beta_raw  = pm.Normal("beta_raw",  0.0, 1.0, dims=("compound", "endpoint"))
        beta      = pm.Deterministic("beta",  beta_raw * sigma_beta, dims=("compound", "endpoint"))

        # Cell-line random effect: α_{l,k}
        alpha_raw = pm.Normal("alpha_raw", 0.0, 1.0, dims=("cell", "endpoint"))
        alpha     = pm.Deterministic("alpha", alpha_raw * sigma_alpha, dims=("cell", "endpoint"))

        # Species random effect: γ_{s,k}  (MH7 — co-fit with MH3, see §7 below)
        gamma_raw = pm.Normal("gamma_raw", 0.0, 1.0, dims=("species", "endpoint"))
        gamma     = pm.Deterministic("gamma", gamma_raw * sigma_gamma, dims=("species", "endpoint"))

        # Compound × cell interaction: δ_{c,l,k}
        delta_raw = pm.Normal("delta_raw", 0.0, 1.0, dims=("compound", "cell", "endpoint"))
        delta     = pm.Deterministic("delta", delta_raw * sigma_delta, dims=("compound", "cell", "endpoint"))

        # ---- Linear predictor ----
        eta = (
            mu_k[endpoint_idx]
            + beta[compound_idx, endpoint_idx]
            + alpha[cell_idx, endpoint_idx]
            + gamma[species_idx, endpoint_idx]
            + delta[compound_idx, cell_idx, endpoint_idx]
        )
        sigma_obs = sigma_eps[endpoint_idx]
        y_obs = pm.Normal("y_obs", mu=eta, sigma=sigma_obs, observed=y, dims="obs")

        # ---- Derived quantities for reporting ----
        icc_cell  = pm.Deterministic(
            "icc_cell",
            sigma_alpha**2 / (sigma_beta**2 + sigma_alpha**2 + sigma_gamma**2 + sigma_delta**2 + sigma_eps**2),
            dims="endpoint",
        )
        icc_inter = pm.Deterministic(
            "icc_inter",
            sigma_delta**2 / (sigma_beta**2 + sigma_delta**2),
            dims="endpoint",
        )

    return model


def fit(model):
    with model:
        idata = pm.sample(
            draws=2000, tune=2000, chains=4,
            target_accept=0.95,
            nuts_sampler="numpyro",     # JAX backend, in-process — same fix as V6.B.5 sprint
            random_seed=20260527,
            return_inferencedata=True,
        )
    return idata
```

**Compute budget**: With ~1000 compounds × ~3 cell-line buckets (U2OS, iPSC-neural, A549) × ~8 endpoints, the model has ~24,000 latent parameters in `δ` plus a few thousand in `β`, `α`, `γ`. On numpyro/JAX with a non-centered parameterization, 4 chains × 4000 draws should complete in 5–15 minutes on the RTX 5070, comfortably under the V6.B.5 8-second baseline scaled to this dimensionality.

**Divergence triage**: If divergences appear, in priority order: (1) tighten `sigma_alpha` prior to `HalfNormal(0.3)`, since the cell-line variance is bounded above by what cpg0000 actually shows (see §5); (2) consider modeling `α_{l,k}` and `γ_{s,k}` as a single combined random effect for the U2OS/iPSC-cortical/iPSC-mDA case where cell-line and species are partly aliased; (3) try the partial-centering trick from Betancourt's case study — center `β` and non-center the rest. Centering only some of the terms is sometimes necessary; in some situations centered parametrization is superior to non-centered.

---

## 5. Empirical anchoring — why this isn't just statistical theater

A skeptical reviewer will ask: "Fine, you have a posterior on `T_{c,k}`. Why should I believe the prior on `σ_α` and `σ_δ` is even in the right ballpark?" The answer is that you can fit those variance components directly from data that already exists in JUMP-CP itself.

### 5.1 cpg0000 — the gold-standard prior calibration dataset

The JUMP-CP pilot dataset `cpg0000-jump-pilot` is purpose-built for this exact problem. It contains 300+ compounds and 160+ genes (CRISPR knockout and overexpression) profiled in both A549 and U2OS cells at two timepoints, designed so that resulting phenotypes can be investigated with matching perturbations targeting the same gene. The dataset contains 51 plates spanning three perturbation types (CRISPR, ORF, Compound) across the two cell lines, with 20,959,860 segmented single cells.

This is the empirical bedrock for MH3. The workflow:

1. Pull the well-aggregated CellProfiler profiles for cpg0000 (compound subset).
2. Fit a stripped-down version of the §3 model with just `β`, `α`, `δ`, `ε` (no species — both lines are human) on matched (compound, A549, U2OS) pairs across the 300+ shared compounds.
3. Report empirical ICC_cell and ICC_int across the ~1700-dim CellProfiler feature space, then aggregated into the V8 endpoint categories (mitochondrial, ER, nuclear, cytoskeletal, …).
4. Use the **empirical posterior median of σ_α and σ_δ** from this fit as the **informative prior** for the full V8 model.

This is a closed loop. The V8 paper Methods can say: *"We calibrated the per-cell-line variance prior on cpg0000 cross-(A549, U2OS) compound pairs, finding empirical ICC_cell of [X.X, Y.Y] (95% CrI) on mitochondrial-organelle endpoints and [X.X, Y.Y] on cytoskeletal endpoints. We use these as informative priors for σ_α in the full V8 joint posterior."* That paragraph is bulletproof.

**Effort line item**: this is ~3 days of work (cpg0000 download + ETL is already scaffolded; the PyMC fit is the §4 model with two cells; result generation is automatic). It should be the **first** thing done in MH3 sprint, not the last.

### 5.2 Gorgogietas/Wilbertz 2025 — the published precedent for U2OS → neuron transfer

The strongest external evidence that U2OS Cell Painting profiles carry neuron-relevant signal comes from Gorgogietas et al. in Scientific Reports (2025), which used patient-derived SNCA triplication (SNCA-4x) and isogenic control (SNCA-corr) midbrain dopaminergic neurons, applied high-content imaging-based morphological profiling, screened 1,020 compounds, and identified top-scoring compounds that restored healthy profiles in SNCA-4x neurons, increasing tyrosine hydroxylase and decreasing α-synuclein levels.

The directly relevant validation step: they validated the effect of identified compounds on cellular phenotypes using JUMP consortium data — compound-induced morphological profiles in U2OS cells by Cell Painting — and observed that, in a different cell type and using a different staining assay, compounds with similar mode of action induced similar morphological profiles, supporting the notion that the identified compounds, while acting on different targets, might possibly act on related pathways leading to a similar phenotypic outcome.

This is the empirical existence proof. It is **not** a proof of universal transferability — it's a single-disease, single-class-of-mechanism (mitochondrial uncoupling) demonstration that some U2OS Cell Painting signals carry information that generalizes to mDA neurons. MH3 turns that anecdotal cross-validation into a per-compound, per-endpoint posterior. The V8 paper should cite Gorgogietas 2025 in the Introduction and in the Discussion as both motivation and external concordance.

### 5.3 Anderson et al. eLife 2025 — iPSC cortical neuron CP as the modality bridge

A second important external anchor: recent eLife work jointly profiled cell morphology and gene expression at single-cell resolution across 60,000 iPSC-derived cortical neurons at three developmental time points using Cell Painting and scRNA-seq, showing that iPSC-derived cortical neurons are a relevant model for a range of brain-related complex traits including schizophrenia and bipolar disorder, and that disease heritability can also be captured in the morphological feature space.

This is important because it gives V8 a published feature-space for iPSC cortical neuron CP that can be directly co-modeled with U2OS CP under the §3 hierarchical specification. If V8 obtains even a modest sample of iPSC cortical neuron CP data on a handful of reference compounds (BDNF stimulators, clemastine, donepezil), it can directly anchor `σ_α` and `σ_δ` for the U2OS-to-cortical-neuron transfer at non-zero data weight rather than relying solely on cpg0000 calibration extrapolation.

### 5.4 The Cross-Cell-Line Methodological Literature

This isn't a novel methodology — there's a parallel line of work in oncology drug screening that has been wrestling with the same identifiability problem for a decade. The closest prior art is Snijder et al. on domain-invariant features for mechanism of action prediction in a multi-cell-line drug screen, which used multi-task autoencoders with adaptive models to construct domain-invariant feature representations across cell lines, applied to two triple-negative breast cancer cell lines. The V8 MH3 approach is the explicit-Bayesian version of that domain-invariance argument — instead of learning a domain-invariant representation, you decompose the perturbation effect into a transferable component and a cell-specific residual.

Methodologically adjacent and worth citing in V8 Methods: Dirmeier & Beerenwinkel (2022) structured hierarchical models for probabilistic inference from perturbation screening data, which uses classical hierarchical models combined with Markov random fields to encode biological prior information for high-noise pan-cancer perturbation screens. The MH3 formulation is the simpler, design-focused cousin (random effects rather than MRF priors) but the philosophical lineage is the same: noisy multi-context perturbation data is exactly the regime where hierarchical Bayesian models earn their keep.

### 5.5 What about chemCPA?

chemCPA already has cell-type covariates in its conditioning. CPA allows context transfer — predicting the effect of a perturbation on unseen cell types or transferring perturbation effects from one context to another — and enables batch effect removal on a latent space and gene expression space. chemCPA is an encoder-decoder architecture for studying perturbational effects of unseen drugs, combining transfer learning with architecture surgery, demonstrating that training on existing bulk RNA HTS datasets can improve generalization performance.

But chemCPA's cell-type embedding is a **point estimate** in the decoder, not a posterior over a transferable vs. cell-specific decomposition of the compound effect. MH3 is doing for the V8 joint posterior what chemCPA does not do for itself: turn the cell-type contribution into a quantified uncertainty rather than an absorbed latent. If reviewers ask why V8 doesn't just use chemCPA's cell-type embedding directly, the answer is that chemCPA's design lets you *predict* an effect in an unseen cell line, but doesn't give you the per-compound transferability posterior that the V8 (L, L, H) discovery score requires for calibrated downstream gating.

### 5.6 CellPainTR and the batch-vs-biology question

The transferability question is entangled with batch effects. CellPainTR (2025) is a Transformer-based architecture designed to learn foundational representations of cellular morphology robust to batch effects, validated on JUMP, where it outperforms ComBat and Harmony in both batch integration and biological signal preservation. MH3 lives at a different level of the stack — it's a perturbation-effect decomposition, not a feature-extraction step. A clean V8.3 architecture should run cpg0016 through CellPainTR (or equivalent batch-aware embedder) **first**, then feed the de-batched embeddings into the MH3 hierarchical model. Otherwise the `α_{l,k}` random effect will mop up site-and-batch variance instead of cell-line biology, which would understate the transferable component (good news in the wrong way).

---

## 6. The defense strategy — how MH3 actually answers the reviewer

When the reviewer asks "How do you defend U2OS-to-brain transfer?", the V8 Methods + Discussion answers in four concrete sentences:

1. **"We decompose the compound effect into a transferable cross-context component β and a per-cell-line interaction δ."** (cite §3 model)

2. **"On the cpg0000 pilot dataset, where the same 300 compounds were profiled in both U2OS and A549, we estimate σ_α/σ_β = X.X and σ_δ/σ_β = Y.Y, giving an empirical ICC_int of Z.Z on the V8 endpoint feature set."** (cite §5.1 calibration)

3. **"For each compound in our (L, L, H) shortlist, we report the per-compound transferability index T_c with 95% credible intervals; compounds with T_c < 0.6 are flagged as U2OS-restricted in our published Tier 3 list."** (cite §3 derived quantities)

4. **"Our shortlist is concordant with the independently published Gorgogietas et al. 2025 mDA neuron findings, which identified mitochondrial-uncoupling-related cellular phenotypes via the same cross-cell-line validation strategy."** (cite §5.2)

That's the response. It's three lines of math plus one external concordance — which is what passes peer review at Nat Mach Intell.

---

## 7. Why MH3 and MH7 should be done jointly, not sequentially

The current `GAPS_AND_RESEARCH_DIRECTIONS.md` lists MH3 (per-cell-line) and MH7 (per-species, mouse-to-human iPSC-MEA) as separate items with separate ~1-week estimates. **Don't do them sequentially.** Three reasons:

**Identifiability collision.** When the dominant species split in the data is "mouse iPSC neurons vs. human U2OS," the per-cell-line and per-species random effects are partly aliased. If you fit MH3 first without `γ_{s,k}`, the `α_{l,k}` posterior will absorb most of the mouse-vs-human variance, then fitting MH7 on top will require refitting from scratch and the prior on `σ_α` will already be wrong.

**Shared infrastructure.** Both terms are non-centered random effects with identical prior structure. The §4 PyMC model adds MH7 as four lines of code on top of MH3. The 2-week-honest estimate for MH3 already includes MH7 if they're done together; doing them separately roughly doubles the timeline.

**Cleaner Discussion paragraph.** The V8 Discussion is much stronger with a paragraph that says "we explicitly modeled variance attributable to cell line (`σ_α`), to species (`σ_γ`), and to compound-cell interaction (`σ_δ`), with ICC values of [a, b, c]" than with two paragraphs covering them separately. Combined treatment also lets you make the honest claim that "the cell-line variance dominates the species variance for mitochondrial endpoints but is dominated by species variance for synaptic-electrophysiology endpoints" — which is exactly the kind of substantive biological inference that justifies the methodological complexity.

**Recommendation**: revise the sprint sequence so the MH3 effort line item explicitly bundles MH7. Total: ~2 weeks, not ~2 weeks then ~1 week.

---

## 8. Pre-registered validation gates (OSF amendment for §7)

The V8 OSF pre-registration currently mentions per-cell-line random effects at §7 as a future amendment. Concrete pass/fail criteria for the production fit, locked before unblinding:

| Gate | Criterion | Pass | Conditional Pass | Fail |
|------|-----------|------|------------------|------|
| G1 — convergence | R̂ ≤ 1.01, ESS_bulk ≥ 400 on all σ hyperparams | both met | one met | neither |
| G2 — divergence count | divergences / total draws | < 0.5% | 0.5–2% with prior tightening rationale | > 2% |
| G3 — cpg0000 calibration consistency | informative prior posterior overlap with cpg0000 marginal | 95% CrI overlap | 50% overlap | non-overlapping |
| G4 — transferability posterior identifiability | std(δ) / std(β) prior-to-posterior contraction | ≥ 30% | 10–30% | < 10% |
| G5 — Gorgogietas concordance | of top-20 V8 (L, L, H) hits with T_c ≥ 0.6, fraction that overlap with the Gorgogietas mDA hit list at compound-or-MoA level | ≥ 25% | 10–25% | < 10% |
| G6 — Mondrian conformal coverage (MH4 hand-shake) | empirical coverage of the conformal 80% interval over held-out compounds | 75–85% | 70–80% / 80–90% | <70% or >90% |

G3 is the most important gate. If the cpg0000-calibrated prior is far from the full-fit posterior, then either the cpg0000 pilot is unrepresentative of `cpg0016` (likely, given site-and-batch differences) or the model is misspecified. Either is a publishable finding but changes the V8 paper substantially.

G5 is the cheapest external concordance check and the most defensible against reviewer skepticism — a low-overhead, high-credibility validation that the model's transferability scores recover an independently published cross-cell-line concordance pattern.

---

## 9. Honest limitations

This section is for the V8 Discussion as a Limitations subsection. None of these are dealbreakers; all of them need to be stated explicitly.

**L1 — The cpg0000 prior is A549-vs-U2OS, not iPSC-neuron-vs-U2OS.** The empirical ICC_cell calibrated on cpg0000 is from two epithelial-like cancer lines. The variance between U2OS and iPSC cortical neurons is plausibly larger. The honest framing in the V8 paper is that cpg0000-calibrated σ_α is a **lower bound** on the true U2OS-to-neuron variance, and the production model treats it as such by using a HalfNormal prior with the cpg0000 posterior median as scale (not as a tight point estimate).

**L2 — Cell-Painting U2OS has no synaptic, electrophysiological, or arborization phenotypes by construction.** U2OS doesn't form synapses or fire action potentials. The transferability posterior can be high for cytoskeletal, mitochondrial, ER, and autophagy endpoints; it is structurally undefined for synaptic-specific endpoints. The V8 (L, L, H) gate should explicitly **exclude** synaptic-electrophysiology endpoints from the U2OS-anchored portion of the score and source those endpoints only from the iPSC-MEA modality.

**L3 — δ_{c,l,k} estimation requires multi-cell-line data per compound.** For compounds appearing only in `cpg0016` (U2OS-only), `δ_{c,U2OS,k}` and `β_{c,k}` are not separately identifiable. The posterior on `T_{c,k}` for those compounds will be wider, properly reflecting the lack of identifiability. The V8 (L, L, H) score should weight these compounds lower in the ranking, or alternatively use the population-level σ_δ posterior as the prior on `δ_{c,U2OS,k}` (partial pooling, which the model already does naturally).

**L4 — Roberts ceiling still applies.** Even a model that perfectly defends the transfer claim is bounded by the Roberts 2020 g ≈ 0.50 90% credible upper bound on cognitive enhancement effect sizes. MH3 doesn't loosen this ceiling; it tightens the calibration around it.

**L5 — Batch effects vs. cell-line effects.** Even with cpg0016 being U2OS-only, the 12 partner sites generate site-level batch variance that aliases with the cell-line effect in any analysis that doesn't first batch-correct. The MH3 model assumes `cpg0016` has been embedded by a batch-aware feature extractor (CellPainTR or equivalent) **before** entering the joint posterior. If the V8 pipeline currently feeds raw CellProfiler features straight in, the `σ_α` estimate will be inflated by site batch variance.

**L6 — Compound × cell-line interaction is an enormous tensor.** With C compounds × L cell lines × K endpoints, the `δ` tensor scales linearly in C·L·K. For 1000 compounds × 3 cells × 8 endpoints that's 24,000 parameters — fine. For 100,000 compounds × 10 cells × 30 endpoints that's 30,000,000 parameters — not fine. The honest paper framing is that MH3 is fit on the curated V8 anchor subset (~1000 compounds) not on the full cpg0016 catalog; the trained model's `σ_β`, `σ_α`, `σ_δ` hyperparameters are then used as fixed values for inference on the remaining compounds (empirical Bayes).

---

## 10. Sprint plan

**Total effort: 2 weeks (revised up from the 1 week in current `GAPS_AND_RESEARCH_DIRECTIONS.md`)**

| Day | Task | Deliverable |
|-----|------|-------------|
| 1 | cpg0000 ETL: download compound subset well-aggregated profiles for A549 + U2OS | `data/cache/cpg0000/profiles_compound_a549_u2os.parquet` |
| 2 | Fit two-cell-line stripped-down §3 model on cpg0000 | `reports/v8_mh3_cpg0000_calibration.md`, posterior parquet |
| 3 | Generate cpg0000 ICC table by endpoint category | calibration table for V8 Methods §X |
| 4–5 | Refactor V8 joint posterior to add α, γ, δ random effects (PyMC + numpyro) | `scripts/v8_joint_with_random_effects.py` |
| 6–7 | Run full V8 fit with informative cpg0000 priors on σ_α, σ_δ | inference data parquet, divergence triage notes |
| 8 | Compute per-compound T_c for all V8 (L, L, H) shortlist | `reports/v8_transferability_index_per_compound.csv` |
| 9 | Pre-registered gate evaluation (G1–G6) | `reports/v8_mh3_gate_evaluation.md` |
| 10 | Gorgogietas concordance check (G5) | overlap table + visualization |
| 11–12 | V8 Methods + Discussion writing for MH3 + MH7 paragraph | V8 manuscript updates |
| 13 | OSF §7 amendment (informative prior values, gate criteria) | `reports/v8_osf_amendment_mh3.md` |
| 14 | Sanity-check holdout: refit with one cell-line withheld, check β recovery | sensitivity analysis notes |

**Dependencies**: cpg0000 must be downloaded (S3 boto3 UNSIGNED is already validated per recently-resolved gap #3); numpyro must be installed (resolved per recently-resolved gap #1); RTX 5070 GPU should be available for the full V8 fit if it exceeds CPU-feasible time.

**De-risked**: the infrastructure for both data acquisition and inference is already in place. The only genuinely novel work is the model specification and the cpg0000 calibration, neither of which has hidden surprises.

---

## 11. What this still doesn't solve

Even with MH3 + MH7 done well, the V8 paper has open Discussion items that no amount of hierarchical modeling can address:

- **Prospective wet-lab validation** — `MH9` and the wet-lab handoff doc are the only paths to external validation. MH3 makes the in-silico shortlist more defensible; it doesn't make it experimentally verified.
- **The MH3 transferability score is a model-internal quantity.** A reviewer can still ask "but why should I believe the model's decomposition of variance is correct?" The cpg0000 calibration + Gorgogietas concordance are the two answers, but they're partial.
- **Brain-region specificity** — U2OS-to-cortical-neuron and U2OS-to-mDA-neuron and U2OS-to-cerebellar-Purkinje are all different transfer problems. MH3 lumps them, mostly because iPSC differentiation protocols for the latter two are not widely represented in public CP datasets. The honest framing is that V8 estimates a generic "U2OS-to-neuron-lineage" transferability, not a brain-region-specific one.

Stating all three of these in the Limitations is the price of admission for an honest publication. None of them downgrades MH3; they bound its claim.

---

## 12. Decision

**Score**: this is a **9/10** must-have research direction for the V8 paper. The current MH3 entry in `GAPS_AND_RESEARCH_DIRECTIONS.md` is correct to flag it as "the single most important methodological improvement"; the only revision needed is bumping the effort from 1 week to 2 weeks and explicitly bundling MH7.

**Recommendation**: promote MH3 + MH7 (combined) to the top of the methodological sprint queue, ahead of MH4 (Mondrian conformal). Reasoning: G6 above shows that Mondrian conformal calibration naturally consumes the MH3 posterior; doing MH4 before MH3 would require redo-ing the calibration on a posterior that doesn't include random effects. Sequencing MH3 → MH4 saves rework.

**One thing not to do**: do not skip the cpg0000 calibration step (§5.1, days 1–3 of sprint plan). It is the single most defensible empirical anchor in this whole exercise. Without it the V8 paper has a Bayesian model with priors a reviewer can dismiss as arbitrary; with it the model has priors a reviewer can fight only by attacking the JUMP-CP consortium itself.

---

*Document generated as a research-direction synthesis for `GAPS_AND_RESEARCH_DIRECTIONS.md` MH3. Companion to V8 OSF pre-reg §7 amendment, V8 manuscript Methods, and V8 sprint planning. Sources span Chandrasekaran et al. 2023 (JUMP-CP), Gorgogietas et al. 2025 (mDA neurons + JUMP validation), Anderson et al. 2025 (iPSC cortical CP), Argelaguet et al. 2020 (MOFA+), Hetzel et al. 2022 (chemCPA), Dirmeier & Beerenwinkel 2022 (SHM), Snijder et al. 2019 (multi-cell-line MoA), and the PyMC/Stan hierarchical-modeling literature.*
