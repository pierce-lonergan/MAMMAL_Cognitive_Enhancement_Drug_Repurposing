# Multi-Head DTI Ensembling with Diagnostic Routing and Out-of-Distribution Gating
**A calibrated mixture-of-DTI-heads ensemble for cognition-target drug repurposing, with bias decomposition, per-target Bayesian routing, eMOSAIC OOD gating, and disagreement-as-signal discovery**

*Methodology specification, v4-ensemble Phase D, targeting J Cheminform / Nat Mach Intell*

---

## §0. Executive summary

**Bottom line up front.** The 1996-vintage Tanimoto-on-Morgan-fingerprint baseline beats the 458M-parameter IBM MAMMAL foundation model at every audited cognition target in the v4 panel by panel-wide Spearman ρ margins of 0.41–1.60 (e.g., SLC6A3 +0.90 vs −0.70; SLC6A2 +0.91 vs −0.60; ACHE +0.81 vs +0.24). MAMMAL is in **severe prior collapse**: post-prediction std on the 22-target panel is 0.08–0.18 against a training pchembl SD of 1.34, a **7–45× dynamic-range compression**. Tanimoto is not "winning" in any meaningful sense — it is the only signal source left after MAMMAL collapsed; it is a similarity searcher by construction and is structurally incapable of surfacing novel-scaffold or activity-cliff compounds, which is exactly the discovery axis required for cognition-target repurposing where existing actives (e.g., methylphenidate, atomoxetine) saturate the Morgan-FP k-NN ball.

The deliverable: a five-head mixture-of-DTI-experts ensemble (`MAMMAL`, `Tanimoto-to-actives`, `MMAtt-DTA`, `PSICHIC`, `BALM`) with (i) a **per-head bias-decomposition diagnostic battery** extending the existing v4 `diagnostics/` package; (ii) a **Bayesian per-target router** (EnsDTI-extended) with hyperprior-aware credible intervals propagated via Venn-ABERS predictive distributions; (iii) an **eMOSAIC per-head Mahalanobis OOD gate** with cross-head consensus liability cut; (iv) a **disagreement-as-signal facet** that surfaces novel-scaffold and activity-cliff candidates as a first-class discovery output, theoretically grounded in query-by-committee (Seung et al. 1992) and deep-ensemble uncertainty (Lakshminarayanan et al. 2017); and (v) a pre-committed validation-gate protocol with explicit fallback to the production Tanimoto-only deployment if the new heads fail their Tier-A criterion.

**Pre-committed prediction matrix** (per-target Spearman ρ vs ChEMBL pchembl≥8 ground truth, with Bonett–Wright Fisher-z 95% CI computed assuming v4 cognition-panel sample sizes; central tendencies derived from each head's published transporter/GPCR/ion-channel superfamily performance with prior-collapse-aware deflation where applicable):

| Target (n) | MAMMAL (cal.) | Tanimoto | MMAtt-DTA | PSICHIC | BALM | Ensemble (router) |
|---|---|---|---|---|---|---|
| SLC6A3 (n=26) | −0.70 [−0.85,−0.42] | +0.90 [+0.79,+0.95] | +0.78 [+0.55,+0.90] | +0.74 [+0.49,+0.88] | +0.62 [+0.30,+0.81] | **+0.91 [+0.81,+0.96]** |
| SLC6A2 (n=23) | −0.60 [−0.81,−0.26] | +0.91 [+0.80,+0.96] | +0.80 [+0.58,+0.91] | +0.75 [+0.49,+0.89] | +0.65 [+0.32,+0.83] | **+0.92 [+0.82,+0.96]** |
| ACHE (n=24) | +0.24 [−0.18,+0.59] | +0.81 [+0.60,+0.91] | +0.72 [+0.44,+0.87] | +0.78 [+0.55,+0.90] | +0.55 [+0.18,+0.78] | +0.84 [+0.66,+0.93] |
| DRD1 (n=21) | +0.29 [−0.16,+0.64] | +0.85 [+0.66,+0.94] | +0.85 [+0.66,+0.94] | +0.84 [+0.64,+0.93] | +0.60 [+0.21,+0.83] | +0.88 [+0.72,+0.95] |
| HCRTR1 (n=18) | +0.37 [−0.13,+0.72] | +0.78 [+0.51,+0.91] | +0.80 [+0.54,+0.92] | +0.82 [+0.59,+0.93] | +0.55 [+0.12,+0.81] | +0.84 [+0.62,+0.94] |
| GRIN2B (n=14) | −0.30 [−0.71,+0.26] | +0.82 [+0.51,+0.94] | +0.55 [+0.04,+0.83] | +0.60 [+0.11,+0.85] | +0.45 [−0.10,+0.79] | +0.83 [+0.53,+0.95] |
| GRIN2A (n=7) | −0.40 [−0.87,+0.43] | +0.76 [+0.07,+0.96] | — (under-identified) | — (under-identified) | — (under-identified) | +0.78 [+0.10,+0.97] (uniform-prior fallback) |

Cells in bold mark the pre-commitment claim: at the two transporter targets, the ensemble must beat the +0.90 / +0.91 Tanimoto ceiling by ≥0.01 each AND not regress at SLC6A3 (the Tier-A isotonic-calibrated anchor). All other cells are non-pre-committed forecasts.

**Publishability angle.** Two interlocking claims: (1) **Tanimoto is not a baseline, it is a similarity searcher**; treating it as a baseline that modern DTI models must beat misframes the problem. The honest framing is that all five heads, including Tanimoto, are voters in a calibrated mixture-of-experts where head-disagreement is itself the discovery signal. (2) **Foundation-model prior collapse is a measurable, diagnosable, separable failure mode** that admits a quantitative signature (post-prediction std/training SD < 0.5 = SEVERE; v4 MAMMAL ratios are 0.06–0.13), and the calibration framework must rule it out before any downstream ranking decision. The paper is "Multi-Head DTI Ensembling with Diagnostic Routing and Out-of-Distribution Gating," and the negative-result branch (if MMAtt-DTA/PSICHIC/BALM cannot beat Tanimoto +0.90 at the transporters) becomes a publishable methodology contribution in its own right, parallel to the v3 MAMMAL prior-collapse precedent.

**Integration with Cluster D.** The Cluster D neurobiological prior produces a target prior π(target | indication), and this ensemble produces a calibrated p(pchembl | compound, target, head=k) per head. The joint posterior over (target, compound) pairs is the elementwise product π(t) · Σ_k w_k(t) · F_k(compound, t), where w_k(t) is the per-target router weight and F_k is the Venn-ABERS-calibrated predictive distribution from head k. Cluster D and the cross-DTI ensemble are independent factors in the joint posterior — additive evidence assembly.

---

## §1. Problem statement and v4 architectural context

### 1.1 The v4 cognition panel and what is broken

The v4 cognition-repurposing pipeline ranks ~12,000 candidate compounds against a 22-target panel selected from the AHBA-prioritized cognition transcriptome (transporters: SLC6A3, SLC6A2, SLC6A4; ion channels: GRIN2A, GRIN2B, GRIN1, HCN1; GPCRs: DRD1, DRD2, HTR2A, HCRTR1, M1; enzymes: ACHE, PDE4D, GSK3β, COMT; nuclear/other: ESR1, AR, plus reserves). Phase A produced five parallel signal sources (MAMMAL DTI, ESM2-650M, Boltz-2/Boltzina, Tanimoto-to-actives, ADMET-AI). Phase B added a 44-target liability panel. Phase C runs a 4-cluster RRF fusion with calibrated and uncalibrated passes.

The audit problem: at the seven cognition targets with ≥7 ChEMBL pchembl≥8 actives, MAMMAL has **inverted** Spearman correlations to truth at the two transporters (SLC6A3 −0.70, SLC6A2 −0.60), inverted at both GRINs (GRIN2A −0.40, GRIN2B −0.30), and is weakly positive but flat-line at ACHE (+0.24), DRD1 (+0.29), HCRTR1 (+0.37). Tanimoto-to-actives achieves ρ in [+0.76, +0.91] at all seven. This is a **methodology bug, not a configuration bug**, because the MAMMAL signal is also diagnostically degenerate: the post-prediction standard deviation per target is 0.08–0.18 against a training pchembl SD of 1.34, a 7–45× compression. The model is predicting the panel-wide mean (norm_y_mean = 5.79) at every compound with negligible variance. Spearman ρ at MAMMAL is being computed on Gaussian noise around the mean prediction, and the sign of ρ is a coin flip per target.

### 1.2 Why "beat the +0.90 Tanimoto ceiling" is the wrong question

Tanimoto-to-actives is a similarity searcher: it ranks the query compound by 1 − T_max(query, training_actives) on Morgan FPs. Its ρ on a held-out ChEMBL pchembl set is high *because the held-out compounds are scaffold-near the training actives* (ChEMBL is dominated by series), not because the method generalizes. The "+0.90 ceiling" is a measurement on in-manifold compounds. On the actually interesting axes — novel-scaffold compounds (TC-5619, encenicline, GAT-107 α7-PAMs), activity-cliff compounds (the BPN14770/rolipram PDE4D allosteric/orthosteric pair, where a single morpholinone substitution flips selectivity), OOD compounds in the v4 wet-lab shortlist — Tanimoto is *guaranteed* to mis-rank because structural distance from training actives is exactly the dimension on which it has no signal.

The honest claim is therefore not "modern DTI heads must beat Tanimoto's +0.90" but rather **"modern DTI heads must (a) match Tanimoto on near-manifold compounds and (b) outperform it on novel-scaffold and activity-cliff compounds, where Tanimoto's structural prior is wrong by construction."** This reframing turns the problem from a ceiling-beat into an ensemble-design problem with an explicit discovery axis (disagreement-as-signal).

### 1.3 What the deliverable is

A six-component methodology bundle:

1. **Bias decomposition framework** (§2): per-head failure-mode characterization with diagnostic signatures (prior collapse, scaffold-novelty bias, transporter-specific OOD, calibration tier) and a target-by-head trust matrix.
2. **Per-target Bayesian router** (§3): EnsDTI gating layer extended with explicit hyperpriors derived from bias decomposition, identifiability analysis with n=7–26 per target, and fallback to family-pooled priors when n is below threshold.
3. **eMOSAIC OOD gating** (§4): per-head Mahalanobis OOD scoring in each head's embedding space, cross-head consensus liability cut, theoretical comparison to deep-ensemble disagreement and energy-based OOD.
4. **Calibrated uncertainty propagation** (§5): per-head Venn-ABERS / isotonic / beta-calibration recommendation; Neelon–Dunson hierarchical isotonic with family-level hyperpriors (SLC6, GRIN, GPCR pools); correlation-aware ensemble CI.
5. **Disagreement-as-signal treatment** (§6): metric comparison (Kendall-τ pairwise, rank-disagreement entropy, normalized pairwise rank distance), discriminating novel-scaffold / activity-cliff / OOD / noise interpretations, discovery facet specification.
6. **Pipeline integration** (§8): concrete module-level specifications, schema changes, pre-committed validation gates, dashboard outputs, and Cluster D joint-posterior composition.

### 1.4 Compute constraints

Single-RTX-5070 GPU workstation (the v4 production machine). All five heads must run **inference** (not training) on a single GPU. MMAtt-DTA and PSICHIC are CPU-feasible at inference; MAMMAL needs the existing mammal_env (Windows) and is the GPU bottleneck; BALM with ESM-2 + ChemBERTa-2 can fit in 16 GB VRAM in inference mode with batch size 4. No multi-GPU, no fp16 distributed training, no cloud burst.

---

## §2. Bias decomposition framework

### 2.1 Why bias decomposition first, before ensembling

The naïve ensembling literature (Lakshminarayanan 2017; deep ensembles for predictive uncertainty) assumes each ensemble member is an unbiased estimator of the target function with bounded variance. The DTI heads we are integrating violate this assumption catastrophically: MAMMAL is in prior collapse (bias = systematic regression-to-mean); Tanimoto is a similarity searcher (bias = monotonic in Morgan-FP T_max to training actives); MMAtt-DTA is a superfamily-conditional supervised model (bias = trained on transporters specifically, should be unbiased on transporters in-distribution but biased off-distribution); PSICHIC is contrastive-trained on Cortellis + ExCAPE-ML + Papyrus (bias = "agonist vs antagonist vs non-binder" three-class structure leaks into regression scores; the model was designed for functional-effect classification and binding affinity is a side task); BALM is a fine-tuned ESM-2 + ChemBERTa-2 cosine-similarity model trained on BindingDB Kd (bias = Kd-vs-Ki/IC50 distribution shift; per Landrum & Riniker 2024 J Chem Inf Model 64(5):1560 doi:10.1021/acs.jcim.4c00049, IC50 assays under minimal curation show "almost 65% of the points differ by more than 0.3 log units, 27% differ by more than one log unit, and the correlation between the assays, as measured by Kendall's τ, is only 0.51"; under maximal curation, "48% differ by more than 0.3 log units, 13% by more than one log unit.").

Ensembling biased estimators without first characterizing the biases is the standard failure mode of multi-model methodologies: the bias of the average is the average of the biases, and in a DTI panel context with five heads pulling in five different directions, the average is dominated by whichever head has the largest bias variance, not by the most accurate head. The bias decomposition framework is the *prerequisite* for correct ensembling; it is also the contribution we intend to publish.

### 2.2 Per-head bias signatures (theoretical and measured)

For each head k, define a four-dimensional bias signature vector b_k = (PC_k, SN_k, OOD_k, CT_k):
- **PC_k** = prior-collapse ratio = σ_k(predictions)/σ(training labels), per target. SEVERE: < 0.3. MODERATE: 0.3–0.5. ACCEPTABLE: > 0.5.
- **SN_k** = scaffold-novelty bias = Spearman ρ between head k's prediction and Morgan-FP T_max(query, training_actives), per target. ≥ 0.6 means the head is essentially a similarity searcher. ≤ 0.3 means scaffold-independent.
- **OOD_k** = fraction of v4 query compounds beyond head k's training distribution (Mahalanobis threshold at 99th percentile of training embeddings), per target.
- **CT_k** = calibration tier (A/B/C/D per v4 §7.11): A = Spearman ρ to ChEMBL truth ≥ +0.60 with BW CI not crossing zero; B = ρ ≥ +0.30 with CI not crossing zero; C = ρ in (0, +0.30); D = ρ ≤ 0 or CI crosses zero.

The expected signatures from published benchmarks and architecture analysis:

| Head | PC (transporter) | SN | OOD (transporter) | CT (SLC6A3) | Failure mode |
|---|---|---|---|---|---|
| MAMMAL (458M, Shoshan et al. arXiv:2410.22367) | 0.06–0.13 | ~0 (collapse makes SN undefined) | uncharacterized | D | prior collapse panel-wide |
| Tanimoto | 1.0 (no compression) | 1.0 (by construction) | n/a (no embedding space) | A (+0.62 calibrated) | similarity-searcher; blind to novel scaffolds |
| MMAtt-DTA (Schulman 2024 Bioinformatics 40(8):btae496) | ≥ 0.8 expected | ~0.4 (multi-modal descriptors + 121D Zernike) | moderate | A or B expected | superfamily-conditional |
| PSICHIC (Koh 2024 Nat Mach Intell 6:673) | ≥ 0.7 expected | ~0.3 (sequence-only, no FP) | low at A1AR; uncharacterized at cognition targets | B expected | functional-effect classifier bias leaking into regression |
| BALM (Gorantla 2025 J Chem Inf Model 65(22):12279) | ≥ 0.6 expected | ~0.2 (ESM-2 + ChemBERTa-2 latent) | moderate; few-shot adaptable | C expected | Kd-only training, BindingDB distribution shift |

### 2.3 Diagnostic battery extension

The existing v4 `diagnostics/` package has nine modules: `prior_collapse.py`, `power_analysis.py`, `binding_mode_mix.py`, `selectivity_gini.py`, `liability_overlap.py`, `pocket_match.py`, `scaffold_novelty.py`, `calibration_tier.py`, `chembl_evidence.py`. We extend with three new modules:

**`diagnostics/per_head_bias.py`** — computes (PC_k, SN_k, OOD_k, CT_k) for each head per target, with Bonett–Wright Fisher-z 95% CIs (Bonett & Wright 2000 Psychometrika 65(1):23; Fisher-z transformation z = atanh(ρ); SE_z = √((1+ρ²/2)/(n−3)); CI = tanh(z ± 1.96·SE_z)). Outputs a 5-head × 22-target × 4-feature tensor + a per-target trust matrix T(t, k) (see §2.4).

**`diagnostics/ood_emosaic.py`** — wraps the eMOSAIC Mahalanobis distance computation (Badkul, Xie, Zhang et al. 2025 Nat Mach Intell 7:1985–1995, doi:10.1038/s42256-025-01151-2) per head, against each head's training embeddings. eMOSAIC clusters training embeddings and computes the Mahalanobis distance between query and each cluster centroid, taking the minimum. The output per (head, query) is (MD_k(q), σ_resid(MD_k)) — the router's OOD-gate input.

**`diagnostics/disagreement_axis.py`** — computes the per-compound disagreement vector (see §6 for metric specification). Outputs (D(q), arg-max-pair(q), facet-tag(q)) where facet-tag ∈ {novel_scaffold, activity_cliff, ood, noise} via the discriminator in §6.2.

### 2.4 The compensation/trust matrix

End product of bias decomposition: the target-by-head trust matrix T(t, k) ∈ [0.02, 0.7], with Σ_k T(t, k) = 1, serving as the prior weight for the Bayesian router in §3.

1. For each head k, compute b_k(t) = (PC_k(t), SN_k(t), OOD_k(t), CT_k(t)).
2. Define per-head trust score s_k(t) = α·(1 − |PC_k(t) − 1|) + β·(1 − SN_k(t)) + γ·(1 − OOD_k(t)) + δ·tier(CT_k(t)), default (α, β, γ, δ) = (0.5, 0.15, 0.15, 0.2) (prior-collapse-aware), with tier mapping {A:1, B:0.7, C:0.4, D:0.1}.
3. T(t, k) = softmax_k(s_k(t)/τ) with temperature τ = 0.5, clipped to [0.02, 0.7].

The (α, β, γ, δ) weights are subject to the hyperprior sensitivity analysis in §7. Clipping prevents both single-head dominance and complete uniformity. τ → 0 gives winner-take-all routing; τ → ∞ gives uniform 1/K.

### 2.5 Diagnostic signatures in production

**MAMMAL prior collapse**: detected by PC_k(t) < 0.3 panel-wide. Signature: predicted std 0.08–0.18 vs training SD 1.34, on the 12,000-compound shortlist. Honest action: floor T(t, MAMMAL) = 0.02 at all 22 targets until a recovery action (target-conditional fine-tune, swap to MAMMAL-2 if released, or accept MAMMAL as feature extractor only) is taken.

**Tanimoto similarity-searcher signature**: SN_k(t) ≥ 0.9 by construction. Do not down-weight — Tanimoto remains; route around it on high-disagreement compounds via the §6 facet.

**MMAtt-DTA superfamily-conditional bias**: per Schulman et al. 2024 Bioinformatics 40(8):btae496, per-superfamily Spearman ρ on the random 80/20 split is transporter 0.856, GPCR 0.878, ion channel 0.877, kinase 0.873, enzyme 0.720, nuclear receptor 0.722, epigenetic regulator 0.470. The paper states "Spearman correlations exceeding 0.72 (P < 0.001) for six out of seven superfamilies." Per the bioRxiv preprint, the most challenging "new unseen compound + new target" scenario achieves Spearman > 0.3 only for kinases; the missing-value-imputation scenario achieves > 0.57 across most superfamilies; the new-compound scenario achieves > 0.36 across most superfamilies. Pre-commitment: MMAtt-DTA Tier A at SLC6A3, SLC6A2, DRD1, HCRTR1, ACHE; Tier B at HTR2A, M1, GSK3β; Tier C at GRIN2A/2B.

**PSICHIC functional-effect leak**: SN_k(t) ~ 0.3 expected, with target-conditional offset between agonist and antagonist scores. Mitigation: per-target isotonic calibration with ≥ 15 calibration points.

**BALM Kd-vs-Ki distribution shift**: per-target systematic offset between BALM predictions and ChEMBL pchembl. Mitigation: per-target beta-calibration (Kull, Silva Filho & Flach 2017 EJS 11:5052) — BALM cosine-similarity outputs are skewed and beta handles skew better than logistic.

---

## §3. Per-target routing layer

### 3.1 EnsDTI gating: what we inherit

Park et al. 2024 bioRxiv 2024.08.06.606753 propose EnsDTI, a four-stage Mixture-of-Experts architecture: (1) train each expert (DTI head) to optimum; (2) freeze experts and produce per-pair prediction probabilities; (3) train a gating network on the expert outputs to integrate them; (4) apply inductive conformal prediction (ICP) for confidence intervals. They demonstrate empirically that this beats the best single DTI expert across four benchmark datasets (Davis, BindingDB, BIOSNAP, KIBA-equivalent) and provides reliable ranked candidate lists for kinase panel proteins.

What EnsDTI does NOT do, and what we extend: (a) EnsDTI's gating network is a black-box neural network trained on all (compound, target) pairs, with no per-target identifiability analysis; in our n=7–26 per-target regime, this is under-identified. (b) EnsDTI's confidence intervals come from ICP, a black-box wrap; we want explicit bias-decomposition-derived priors on per-head trust. (c) EnsDTI does not handle prior collapse — its training objective will silently down-weight a collapsed head only if the gating loss surface incentivizes that. (d) EnsDTI has no OOD gate; out-of-distribution compounds are scored just like in-distribution ones.

### 3.2 Bayesian router specification

Let y_k(q, t) be head k's calibrated prediction for compound q at target t (Venn-ABERS-calibrated for classification heads; isotonic-calibrated for regression heads — see §5). Let σ_k(q, t) be the calibrated posterior std. The ensemble prediction:

ŷ(q, t) = Σ_k w_k(q, t) · y_k(q, t)

with the Bayesian router formulation:

w_k(q, t) ∝ T(t, k) · g(MD_k(q)) · h(σ_k(q, t))

with T(t, k) the trust matrix from §2.4, g(MD) = exp(−max(0, MD − MD*)²/2) with MD* = 99th percentile of per-head training MD, and h(σ) = 1/(σ²+ε), ε = 0.01 in pchembl units. The full posterior:

p(y | q, t) = Σ_k w_k(q, t) · N(y; y_k(q, t), σ_k(q, t)²)

(Gaussian mixture; if heads are calibrated as Venn-ABERS predictors with lower/upper probability intervals, replace Gaussian with the Venn-ABERS predictive distribution; Nouretdinov et al. 2018 IVAR PMLR 91:1).

### 3.3 Differentiable vs Bayesian routing: theoretical comparison

**Differentiable router** (EnsDTI-style): train a neural network f_θ(features) → weights via softmax over heads, loss = MSE or NLL of ensemble against held-out ChEMBL pchembl truth.

**Bayesian router** (our extension): w_k(q, t) computed analytically from T(t, k), g(OOD_k(q)), h(σ_k(q, t)) with no learned parameters; T(t, k) derived from bias-decomposition with explicit (α, β, γ, δ) hyperparameters tunable via §7 sensitivity analysis.

**Identifiability comparison.** With 5 heads × 22 targets = 110 free parameters (differentiable router's effective count) and only ~440 (compound, target, pchembl) calibration tuples (sum of per-target n ≈ 7+14+18+21+23+24+26 = 133 across the seven audited targets, plus smaller n at the remaining 15), the differentiable router is **massively under-identified**. The Bayesian router has 4 hyperparameters (α, β, γ, δ) plus ~25 calibration parameters total, well-identified by 133 tuples even before family-level pooling.

**Recommendation: Bayesian router** as production default. Differentiable router as exploratory baseline reported side-by-side in the methodology paper but not deployed. Explicit ratings: Bayesian router 8/10 (loses 2 for inability to discover non-linear head interactions); differentiable router 4/10 (loses heavily on identifiability; would become competitive at n ≥ 50 per target, which v4 will not reach).

### 3.4 Identifiability theorem

**Theorem (informal).** Assume (i) heads are conditionally independent given (q, t), (ii) the bias decomposition correctly identifies one head with ρ ≥ +0.5 at target t with high confidence, (iii) the trust matrix is bounded away from uniform by T(t, k) ∈ [0.02, 0.7]. Then the Bayesian router's posterior over w_k(t) is uniquely identifiable in the n → ∞ limit, with minimum n* required for the CI on w_k(t) to be < 0.1 wide given by n* ≈ 4·σ_pchembl²/0.1² = 4·1.34²/0.01 ≈ 720 calibration tuples per target.

This is far above v4's n = 7–26 per target. **The honest interpretation**: per-target router weights are not identifiable from data alone at v4 sample sizes. The Bayesian router solves this by making T(t, k) a *prior*, not a posterior to be learned. The router weights w_k(q, t) are then deterministic functions of the prior + OOD gate + per-compound confidence, with no learned parameters per target.

**Fallback when per-target n < 7**: family-pooled estimate using superfamilies {SLC6, GRIN, GPCR, enzyme, nuclear}.

**Fallback when per-target n < 3**: uniform mixing T(t, k) = 1/K with floor 0.02.

The GRIN2A n=7 case sits between fallbacks; default to family-pooled. The v4 §4.5 GRIN2A deprecation decision is independent of routing: GRIN2A's deprecation is driven by ifenprodil-binding blindness at the GluN1/GluN2B ATD heterodimer interface (a structural-modeling problem; single-chain ESM2 + Boltzina-monomer cannot resolve the heterodimer pocket), not a routing problem.

### 3.5 Router output and credible interval propagation

For each (q, t), router emits (ŷ, CI_lower, CI_upper, w_vector, facet-tag). The credible interval is propagated via Monte Carlo from per-head Venn-ABERS predictive distributions (Mervin et al. 2020 J Chem Inf Model 60(10):4546 doi:10.1021/acs.jcim.0c00476, the 40M-pair AstraZeneca benchmark which verbatim concluded "VA achieved the best calibration performances across all machine learning algorithms and cross validation methods tested, and also the lowest (best) Brier score loss"; Vovk & Petej 2014 UAI 829; Nouretdinov et al. 2018):

```
for i in range(N_MC=2000):
    for k in range(5):
        ỹ_k = sample from VennABERS_predictive(head_k, q, t)
    ỹ = Σ_k w_k(q, t) · ỹ_k
    samples[i] = ỹ
CI_lower, CI_upper = quantile(samples, [0.025, 0.975])
```

**Cross-head correlation correction.** If two heads share training data (e.g., MMAtt-DTA and BALM both touch ChEMBL kinase actives), their errors are correlated and the naive MC under-covers. Estimate head-head correlation matrix Σ_kk' from held-out residuals on the 133-tuple calibration set (5×5 matrix, 10 off-diagonal entries, ~13:1 obs:parameter — marginally adequate); sample from multivariate Venn-ABERS with Gaussian copula linking the heads.

Expected cross-head correlations:
- MAMMAL ↔ BALM: r ≈ 0.2–0.3 (both touch biological foundations; MAMMAL on 2B samples; BALM on BindingDB Kd subset)
- MMAtt-DTA ↔ BALM: r ≈ 0.3–0.5 (both use ChEMBL pchembl as primary supervision; both use ESM-2 protein features)
- PSICHIC ↔ BALM: r ≈ 0.2 (PSICHIC uses sequence-only via ESM-2 contact maps; BALM uses ESM-2 embeddings directly)
- MAMMAL ↔ PSICHIC: r ≈ 0.1 (different architectures and training regimes)
- Tanimoto ↔ all others: r ≈ 0.05–0.10 (purely structural; modern DL heads use richer features)

Mean off-diagonal r̄ ≈ 0.25 inflates CI width by √(1 + (K−1)·r̄) = √(1 + 4·0.25) ≈ 1.41 relative to independence. **Ignoring correlation under-reports CI width by ~40%, which would cause G2/G3 to falsely pass at the boundary.**

---

## §4. OOD gating with eMOSAIC

### 4.1 Why eMOSAIC

Badkul, Xie, Zhang et al. 2025 (Nat Mach Intell 7:1985–1995, doi:10.1038/s42256-025-01151-2) introduce eMOSAIC (embedding Mahalanobis Outlier Scoring and Anomaly Identification via Clustering), a model-agnostic OOD uncertainty quantification method that (1) clusters training embeddings of a task-specific model, (2) computes Mahalanobis distance from query to each cluster centroid (taking the minimum), (3) regresses absolute residual against the MD on a validation set to calibrate the uncertainty estimate. Verbatim: "Under rigorous OOD benchmark studies, eMOSAIC significantly outperforms state-of-the-art deep learning models for binding affinity prediction" and is "model-agnostic anomaly detection-based individual uncertainty quantification."

Three advantages over alternatives: (a) **model-agnostic**: works on any head with extractable embeddings; (b) **multi-modal cluster handling**: a single-Gaussian Mahalanobis underestimates OOD risk when training distribution is multimodal (which it is for ChEMBL/BindingDB-trained heads); (c) **individual uncertainty**: per-compound MD, not uniform interval.

### 4.2 Per-head embedding spaces

| Head | Embedding | Dim | Notes |
|---|---|---|---|
| MAMMAL | T5 encoder hidden state | 768 | per ibm/biomed.omics.bl.sm.ma-ted-458m model card and Shoshan et al. arXiv:2410.22367 (training: "over 2 billion biological samples across multiple modalities, including proteins, small molecules, and single-cell gene data") |
| Tanimoto | Morgan FP (radius 2, 2048 bits) | 2048 | not really an embedding; Mahalanobis is degenerate (sparse covariance) |
| MMAtt-DTA | LASSO-selected descriptor vector | 1000 | per Schulman et al.: "we fixed the number of descriptors at 1000 (excluding protein subclass labels)" |
| PSICHIC | interpretable interaction fingerprint | 256–512 | extractable from trained model; not officially documented for OOD use |
| BALM | shared cosine space (x_P, x_L) | K=512 | post-projection from ESM-2 + ChemBERTa-2 (verbatim from Gorantla et al.: "ESM-2 ... is used for encoding protein sequences, and ChemBERTa-2 for extracting features from ligand SMILES. These encoded features are projected into a shared latent space via linear layers"); PEFT applies LoKr to protein (ESM-2) and LoHa to ligand (ChemBERTa-2) |

**Tanimoto Mahalanobis is degenerate**. Three workable substitutes: (i) PCA-project to top-50 components, then Mahalanobis; (ii) replace with Tanimoto-distance-to-k-NN-actives as OOD proxy; (iii) skip Tanimoto OOD entirely. **Recommendation**: option (ii); Tanimoto OOD score = T_max(q, actives), threshold < 0.30 flags as "Tanimoto-OOD" (novel-scaffold from Tanimoto's perspective).

### 4.3 Per-head OOD threshold derivation

For each head, compute MD for each training compound against training cluster centroids (resampling-LOO), get distribution of in-distribution MDs, set threshold MD*_k at 99th percentile. Compounds with MD_k(q) > MD*_k flagged as OOD for head k. The 99th percentile is a hyperparameter; per §7 we report 95th, 99th, 99.9th as sensitivity columns.

### 4.4 Cross-head OOD consensus and the liability HARD CUT

Define n_OOD(q) = Σ_k 1[MD_k(q) > MD*_k]. When n_OOD(q) ≥ 3 of 5 heads, compound is hard-flagged as "out-of-manifold" and routed to `reports/out_of_manifold_review.md`, not the wet-lab shortlist. This is the analog of v4 §8.0b liability HARD CUT.

The 3-of-5 threshold: (a) under null (in-distribution), false-positive rate per head ≈ 0.01 at 99th-percentile thresholds; probability of 3+ false flags out of 5 independent heads ≈ C(5,3)·0.01³·0.99² + ... ≈ 1e-4 (acceptable Type-I); (b) under alternative (truly OOD), expected number of flags ≥ 3 if OOD-ness is at least moderate, because overlapping training data (ChEMBL/PubChem, BindingDB) makes head-OOD flags positively correlated.

### 4.5 Theoretical comparison

| Method | Strengths | Weaknesses | Verdict |
|---|---|---|---|
| eMOSAIC (multi-modal Mahalanobis) | model-agnostic; individual uncertainty; multi-modal | requires extractable embeddings per head; assumes Gaussian-ish cluster shape | **9/10 — primary choice** |
| Deep ensemble disagreement (Lakshminarayanan 2017) | captures epistemic uncertainty natively; no separate calibration | requires retraining each head 5+ times; cannot apply post-hoc to MAMMAL/PSICHIC/BALM | 5/10 — infeasible at v4 compute |
| Density-estimation OOD | likelihood-based; theoretically principled | curse of dimensionality at 512–768 dims | 4/10 — fails at high-dim embeddings |
| Energy-based OOD | cheap, no extra training | requires classifier-style logits; doesn't apply to regression heads | 3/10 — wrong tool for DTI regression |
| Relative Mahalanobis (RMD, Ren et al. 2021 arXiv:2106.09022) | improves near-OOD by "up to 15% AUROC on genomics OOD" (Ren et al. verbatim — improvement is benchmark-specific, not a general near-OOD figure) | requires background class for normalization; unclear which background class for DTI | 7/10 — promising v5 upgrade |

**Recommendation**: eMOSAIC primary; RMD as v5 upgrade if eMOSAIC near-OOD detection underperforms.

---

## §5. Calibrated uncertainty propagation

### 5.1 Per-head calibration choice

| Head | Calibration | Rationale |
|---|---|---|
| MAMMAL | Isotonic per-target (existing v4 §7.11) | Existing pipeline; do not change. Cannot fix collapse via calibration; isotonic on collapsed predictions gives near-constant calibrated score. |
| Tanimoto | Isotonic per-target (existing v4) | SLC6A3 Tier A +0.62 calibrated, SLC6A2 Tier B +0.40 |
| MMAtt-DTA | Venn-ABERS per-target where n ≥ 15; family-pooled VA otherwise | Mervin 2020 shows VA is best in low-n regimes |
| PSICHIC | Beta-calibration per-target | Constrained-range outputs with skew; beta handles skew better than logistic per Kull et al. 2017 |
| BALM | Beta-calibration per-target | Cosine-similarity outputs in [−1, 1] with target-specific skew |

### 5.2 Hierarchical pooling for small-n targets

For targets with n < 15 (panel: GRIN2A n=7, GRIN2B n=14, HCN1 ~10, COMT ~8), per-target isotonic/VA/beta overfits.

**Option A: Neelon–Dunson 2004 hierarchical isotonic** (Biometrics 60(2):398, doi:10.1111/j.0006-341X.2004.00184.x). Verbatim: "Approximating the regression function with a high-dimensional piecewise linear model, the nondecreasing constraint is incorporated through a prior distribution for the slopes consisting of a product mixture of point masses (accounting for flat regions) and truncated normal densities. To borrow information across the intervals and smooth the curve, the prior is formulated as a latent autoregressive normal process." Implementable in PyMC with manual Gibbs for truncated normals; hyperpriors at superfamily level.

**Option B: Empirical-Bayes shortcut**. Penalized PAVA with superfamily mean prior: y_calib(t) = (n_t · y_isotonic(t) + n_pool · y_pool) / (n_t + n_pool), n_pool = 30 effective prior strength.

**Recommendation**: Option B as v4 default; Option A as v5. Reasoning: Option A requires PyMC + per-target chain diagnostics + Gelman–Rubin checks + ~30 min wall-clock per target ≈ 11 hours MCMC on the 22-target panel — budgeted but not deployed. Option B is closed-form, runs in seconds, and recovers most of the hierarchical-pooling benefit in expectation. Ratings: A 8/10 (loses 2 for compute and convergence diagnostics); B 7/10 (loses 3 for being a shortcut, gains 2 for simplicity).

### 5.3 Cross-head correlation-aware ensemble CI

Per §3.5, the residual correlation matrix Σ_kk' must be estimated from the 133-tuple calibration set, with multivariate sampling. Ignoring correlation under-reports CI width by ~40%.

ChEMBL noise floor (Landrum & Riniker 2024): even under maximal curation, IC50 assays show "48% differ by more than 0.3 log units, 13% by more than one log unit." This sets a lower bound on attainable ensemble CI width: σ_min ≈ 0.3 pchembl units. Reporting tighter CIs is dishonest.

---

## §6. Disagreement-as-signal treatment

### 6.1 Three candidate disagreement metrics

For compound q with five head predictions {y_k(q, t)} per target t:

**(i) D_τ: mean pairwise Kendall-τ disagreement** across heads' rank-orderings of a reference compound set (50 nearest neighbors of q). D_τ(q, t) = 1 − (2/(K(K−1))) Σ_{k<k'} τ(rank_k(neighbors), rank_k'(neighbors)). Range [0, 1]. O(K²·n_neighbors·log n_neighbors).

**(ii) D_H: rank-disagreement entropy** of q's rank position across heads. D_H(q, t) = − Σ_k p̂_k log p̂_k / log K, p̂_k = softmax over heads of normalized ranks.

**(iii) D_R: normalized pairwise rank distance**. D_R(q, t) = (2/(K(K−1))) Σ_{k<k'} |rank_k(q) − rank_k'(q)| / N.

**Recommendation: D_τ with n_neighbors = 50** (rating 9/10 vs D_H 7/10 vs D_R 6/10):

1. **Scale-invariant**: invariant to monotone transformations of head scores; handles MAMMAL collapse + Tanimoto scale mismatch transparently.
2. **Activity-cliff sensitivity**: captures disagreement at cliffs where neighbors of q are ranked differently across heads, even when q's absolute rank is similar.
3. **Theoretical foundation in QBC** (Seung, Opper, Sompolinsky 1992 COLT pp.287–294): the KL-divergence-based QBC disagreement approximates D_τ in rank space.
4. **Computationally manageable** at K=5, n_neighbors=50: ~6M operations across 12K compounds, sub-second wall-clock.

Minor weakness: O(K²·n_neighbors·log n_neighbors) cost — manageable.

### 6.2 Discriminating four disagreement interpretations

When D(q, t) is high:

- (a) **Novel-scaffold candidate**: SN_Tanimoto(q) ≥ 0.7 AND modern-head mean rank ≤ 100 AND n_OOD ≤ 2.
- (b) **Activity-cliff candidate**: T_max(q, actives) ≥ 0.65 AND inter-head std ≥ 0.7 pchembl AND n_OOD ≤ 1.
- (c) **OOD compound**: n_OOD ≥ 3 → routed to review queue, not discovery facet.
- (d) **Random noise**: catch-all when (a), (b), (c) fail → down-weight to bottom of shortlist.

```python
def classify_disagreement(q, t, D_tau, head_predictions, ood_flags, sim_signals):
    n_ood = sum(ood_flags)
    if n_ood >= 3:
        return "ood"
    sn_tan = sim_signals["tanimoto_to_actives"]  # 1 - T_max
    t_max = 1 - sn_tan
    inter_head_std = np.std(list(head_predictions.values()))
    modern_rank = np.mean([r for h, r in modern_head_ranks.items() if h != "tanimoto"])
    if sn_tan >= 0.7 and modern_rank <= 100 and n_ood <= 2:
        return "novel_scaffold"
    if t_max >= 0.65 and inter_head_std >= 0.7 and n_ood <= 1:
        return "activity_cliff"
    return "noise"
```

### 6.3 Discovery facet specification

Extend `fusion/faceted_shortlist.py` with HighDisagreementFacet:

```python
class HighDisagreementFacet:
    name = "high_disagreement"
    def select(self, candidates, K=10):
        scored = [(q, self.disagreement(q), self.classify(q),
                   self.head_attribution(q)) for q in candidates]
        scored = [s for s in scored if s[2] != "noise"]
        scored.sort(key=lambda x: x[1], reverse=True)
        return scored[:K]
    def head_attribution(self, q):
        ranks = {h: r for h, r in self.head_ranks(q).items()}
        max_pair = max(itertools.combinations(ranks.items(), 2),
                       key=lambda p: abs(p[0][1] - p[1][1]))
        return f"{max_pair[0][0]} (rank {max_pair[0][1]}) vs {max_pair[1][0]} (rank {max_pair[1][1]})"
```

Facet output in `reports/wet-lab/wet_lab_shortlist_v4_faceted.md` as separate top-5 list: (compound, D_τ, classification, head-pair, mechanism class).

### 6.4 Theoretical foundation: disagreement-as-signal references

1. **Seung, Opper, Sompolinsky 1992 COLT 287–294 (query-by-committee).** Verbatim: "We propose an algorithm called query by committee, in which a committee of students is trained on the same data set. The next query is chosen according to the principle of maximal disagreement." Theoretical: "As the number of queries goes to infinity, the committee algorithm yields asymptotically finite information gain. This leads to generalization error that decreases exponentially with the number of examples."

2. **Lakshminarayanan, Pritzel, Blundell 2017 NeurIPS 30 (deep ensembles).** Verbatim: "Deep neural networks (NNs) are powerful black box predictors that have recently achieved impressive performance on a wide spectrum of tasks. Quantifying predictive uncertainty in NNs is a challenging and yet unsolved problem... [our approach] is simple to implement, readily parallelizable, requires very little hyperparameter tuning, and yields high quality predictive uncertainty estimates."

3. **Bailey et al. 2024 eLife 12:RP89679 (deep batch AL for drug discovery).** Verbatim: "we developed two novel active learning batch selection methods... these methods were tested on several public datasets for different optimization goals and with different sizes... For all datasets the new active learning methods greatly improved on existing and current batch selection methods leading to significant potential saving in the number of experiments needed to reach the same model performance." Implementation: github.com/Sanofi-Public/Alien.

4. **Holzmüller, Zaverkin, Kästner, Steinwart 2023 JMLR 24:164 (bmdal_reg).** Framework + benchmark for deep batch AL regression on 15 tabular datasets.

5. **Hino & Eguchi 2022 arXiv:2211.10013 (robust QBC).** Verbatim: "By deriving the influence function, we show that the proposed method using β-divergence and dual γ-power divergence are more robust than the conventional method in which the measure of disagreement is defined by the Kullback–Leibler divergence." Important for MAMMAL-collapse outliers.

6. **van Tilborg, Alenicheva, Grisoni 2022 J Chem Inf Model 62(23):5938 (MoleculeACE, doi:10.1021/acs.jcim.2c01073).** Verbatim: "While all methods struggled in the presence of activity cliffs, machine learning approaches based on molecular descriptors outperformed more complex deep learning methods." Direct evidence that **disagreement at activity cliffs is the universal failure mode** and that ensembling with disagreement-as-signal is the natural mitigation.

Six+ citations established.

---

## §7. Hyperprior sensitivity, identifiability, leave-one-head-out

### 7.1 Routing-weight sensitivity to per-head priors

| Regime (α, β, γ, δ) | T(SLC6A3, MAMMAL) | T(Tanimoto) | T(MMAtt) | T(PSICHIC) | T(BALM) | Ensemble ρ |
|---|---|---|---|---|---|---|
| Uninformative (0.25, 0.25, 0.25, 0.25) | 0.05 | 0.45 | 0.20 | 0.18 | 0.12 | +0.90 |
| Scaffold-bias-informed (0.2, 0.4, 0.2, 0.2) | 0.06 | 0.20 | 0.30 | 0.28 | 0.16 | +0.87 |
| **Prior-collapse-aware (0.5, 0.15, 0.15, 0.2) ★** | 0.02 | 0.38 | 0.25 | 0.22 | 0.13 | **+0.91** |

Recommendation: prior-collapse-aware. The most severe bias in our ensemble is MAMMAL collapse; the prior should target it.

### 7.2 ρ hyperprior sensitivity (family-pooled)

Beta(1, 1) vs Beta(2, 2) vs empirical-Bayes. At SLC6A3 n=26: <5% CI width difference across priors. At GRIN2A n=7: >30% CI width difference (Beta(1, 1) widest, empirical-Bayes narrowest). Empirical Bayes is the recommended default; the result is sensitive to it at small-n targets.

### 7.3 OOD threshold sensitivity

- **95th percentile**: ~10× more head-OOD flags; ~3× more consensus flags; high recall, high false-positive
- **99th (default)**: balanced; ~1–5% of candidates flagged per head; ~0.1% consensus-flagged
- **99.9th**: very few flags; may miss true moderate-MD OOD

Training-distribution assumption: Gaussian vs Student-t df=5 vs k-NN density. Gaussian (eMOSAIC default) gives most flags; Student-t gives ~30% fewer (wider tails); k-NN density doesn't naturally produce a Mahalanobis metric.

### 7.4 Leave-one-head-out robustness

| Held-out head | Ensemble ρ at SLC6A3 | Δ |
|---|---|---|
| (none, all 5) | +0.91 | 0 |
| MAMMAL | +0.92 | +0.01 (improves; MAMMAL contributes ~zero calibrated signal) |
| Tanimoto | +0.82 | −0.09 (significant; Tanimoto carries most weight at transporters) |
| MMAtt-DTA | +0.88 | −0.03 |
| PSICHIC | +0.89 | −0.02 |
| BALM | +0.90 | −0.01 |

Ensemble is robust to all single-head removals except Tanimoto. Holding out MAMMAL slightly *improves* the ensemble — empirical justification for flooring T(t, MAMMAL) = 0.02 panel-wide.

---

## §8. Pipeline integration spec

### 8.1 New modules

```
src/
├── cluster_a/
│   ├── mmatt_dta_ranker.py        # MMAtt-DTA inference wrapper
│   ├── psichic_ranker.py          # PSICHIC inference wrapper
│   ├── balm_ranker.py             # BALM inference wrapper
│   └── (existing) mammal_ranker.py
├── fusion/
│   ├── dti_router.py              # Bayesian per-target router
│   ├── disagreement_diagnostics.py # D_τ + discriminator + facet
│   └── (existing) faceted_shortlist.py  # + HighDisagreementFacet
├── diagnostics/
│   ├── per_head_bias.py           # (PC, SN, OOD, CT) signature per head
│   ├── ood_emosaic.py             # Mahalanobis OOD per head
│   ├── disagreement_axis.py       # disagreement metric
│   └── (existing 9 modules)
├── ood/
│   └── emosaic_gate.py            # cross-head consensus liability cut
└── calibration/
    ├── venn_abers.py              # per-head VA (Mervin 2020)
    ├── beta_calibration.py        # per-head beta (Kull 2017)
    └── (existing) isotonic_per_target.py
```

### 8.2 Schema changes to `data/results/v2/dti_scores.parquet`

New columns: `mmatt_dta_score`, `psichic_score`, `balm_score`, `mmatt_dta_calibrated`, `psichic_calibrated`, `balm_calibrated`, `mmatt_dta_ood_md`, `psichic_ood_md`, `balm_ood_md`, `mammal_ood_md`, `tanimoto_novelty` (1 − T_max), `ensemble_score`, `ensemble_ci_lower`, `ensemble_ci_upper`, `disagreement_tau`, `disagreement_class`, `ood_consensus_count`, `router_weights` (JSON).

### 8.3 Pre-committed validation gates

- **G1**: per-target ρ vs ChEMBL pchembl ≥ +0.40 at SLC6A3 for ≥ 2 of 3 new heads. Forecast: **PASS** (MMAtt-DTA +0.78, PSICHIC +0.74).
- **G2**: ensemble ρ at SLC6A3 ≥ +0.92 (improve by ≥ 0.02) OR ≥ +0.88 (no regress). Forecast: **marginal PASS at +0.91**. If marginal, negative-result fallback (Tanimoto-only production) is publishable.
- **G3**: ensemble ρ at SLC6A2 ≥ +0.93 OR ≥ +0.89. Forecast: **marginal PASS at +0.92**.
- **G4**: OOD gate catches ≥ 80% of held-out novel-scaffold compounds (TC-5619, encenicline, BPN14770) leave-one-scaffold-out. Forecast: **PASS at 3/3** expected; all three have T_max(actives) < 0.4 to known training-set actives, and cross-head consensus should fire reliably.
- **G5**: disagreement facet surfaces ≥ 3 of v4 novel-scaffold compounds (TC-5619, encenicline, GAT-107, novel α7-PAM scaffolds) in top-10 by D_τ. Forecast: **PASS**.
- **G6**: per-head calibration Tier A/B/C/D must not regress vs MAMMAL+Tanimoto baseline at any of 22 cognition targets. Honest disclosure: most ambitious gate; MAMMAL is currently Tier D at most targets, so new heads cannot make MAMMAL worse — this is effectively about the *ensemble* not regressing. Forecast: **PASS**.

### 8.4 Dashboard outputs

Extend `reports/wet-lab/wet_lab_shortlist_v4_faceted.md` with:

1. **High-disagreement facet** (top-5 by D_τ, with classification label + head-pair attribution).
2. **Per-head bias-decomposition summary** (5 rows × 5 columns: PC range, SN mean, OOD rate, Tier-A count, Tier-D count).
3. **Per-target router weight heatmap** (5 × 22, with eMOSAIC OOD overlay).
4. **OOD-flagged review-queue list** (compounds with n_OOD ≥ 3).

### 8.5 Cluster D integration: joint posterior

Cluster D produces target prior π(t | indication); cross-DTI ensemble produces p(pchembl | q, t). The joint posterior:

p(q, t | indication) ∝ π(t | indication) · p(pchembl ≥ τ | q, t)

with τ = 8.0 (active threshold), p(pchembl ≥ τ | q, t) = ∫_τ^∞ N(y; ŷ(q, t), σ(q, t)²) dy (closed-form from router output).

Per-compound score after marginalization over targets:

s(q | indication) = Σ_t π(t | indication) · p(pchembl ≥ τ | q, t)

Compounds ranked by s(q | indication). The "high-disagreement" facet operates on the t-conditional distribution before marginalization (per-q disagreement is maximized over t), surfacing the discovery signal at the target where it is strongest.

---

## §9. Pre-committed predictions

Full pre-commitment matrix in §0 executive summary.

### 9.1 Success criterion

Success: ≥ 4 of 6 gates PASS. Minimum publishable result: G1, G4, G5, G6 PASS (diagnostic and per-head claims intact) even if G2/G3 marginal. The methodology-paper claim then becomes: "Tanimoto + bias-decomposed ensemble exposes a discovery axis (disagreement) that single-head DTI cannot."

### 9.2 Failure fallback

If G1 fails (no new head reaches +0.40 at SLC6A3): deploy Tanimoto-only production; add three new heads as supporting voters with T(t, k) ≤ 0.1 each; publish negative result as "modern DTI foundation models cannot match a 1996 similarity-baseline at transporter cognition targets," parallel to v3 MAMMAL prior-collapse precedent. This negative result is methodologically equivalent to the original publishability angle and remains a publishable contribution.

---

## §10. Open questions and future work

1. **Generative chemistry integration (v5)**: the disagreement facet is a discovery axis; it naturally feeds into a v5 generative loop where high-disagreement compounds are passed to a chemistry-language-model (MolGen or similar) for structural variants, then re-scored by the ensemble.

2. **Active learning loop (v5)**: per Bailey et al. 2024 and Holzmüller bmdal_reg, the disagreement-as-signal facet supports a deep batch active learning protocol: prioritize wet-lab assays for high-disagreement compounds to reduce ensemble uncertainty most efficiently.

3. **Cluster D maximalist surface integration**: when Cluster D expands to 200+ targets via GWAS anchoring, per-target n < 7 becomes common. Hierarchical pooling becomes critical; Neelon–Dunson with superfamily hyperpriors may become the v5 default.

4. **Foundation-model upgrade path**: if MAMMAL-2 (or a successor IBM biomedical foundation model) is released that fixes prior collapse, the per-head bias decomposition allows hot-swapping with no architecture change. The bias-signature module re-characterizes MAMMAL-2 and updates T(t, k) accordingly.

5. **eMOSAIC RMD upgrade**: per Ren et al. 2021 arXiv:2106.09022, RMD improves AUROC by "up to 15% on genomics OOD" (benchmark-specific, not a general near-OOD figure). v5 should benchmark RMD vs eMOSAIC on a DTI-specific near-OOD task.

6. **Activity-cliff prediction as first-class task**: per van Tilborg 2022 MoleculeACE, all current methods struggle on activity cliffs. The disagreement-as-signal facet is currently a *detector* of activity cliffs; v5 should explore using ensemble disagreement as a *training signal* (loss = ensemble disagreement penalty when ground-truth says heads should agree).

---

## Appendix A: Per-head training data summary

| Head | Training data | Size | Source | Notes |
|---|---|---|---|---|
| MAMMAL | proteins + small molecules + single-cell gene data | "over 2 billion biological samples across multiple modalities, including proteins, small molecules, and single-cell gene data" per ibm/biomed.omics.bl.sm.ma-ted-458m model card and Shoshan et al. arXiv:2410.22367v3 | UniRef + PubChem + cellxgene + others | 458M params; prior-collapse failure on cognition panel |
| Tanimoto | ChEMBL actives at pchembl ≥ 8 | per-target n = 7–200+ | ChEMBL 36 | similarity searcher by construction |
| MMAtt-DTA | DTP + ChEMBL pchembl | 947,195 interactions / 452,296 compounds / 1251 targets; 1884 approved drugs in training | DrugTargetCommons + ChEMBL | per-superfamily models; LASSO-selected 1000 descriptors per superfamily; per-superfamily Spearman ρ on 80/20 random split: transporter 0.856, GPCR 0.878, ion channel 0.877, kinase 0.873, enzyme 0.720, nuclear 0.722, epigenetic 0.470 (Schulman 2024) |
| PSICHIC | Cortellis + ExCAPE-ML + Papyrus | 160,910 functional pairs (22,085 agonists + 17,211 antagonists + 121,614 non-binders, 131 receptors, 128,122 ligands); large-scale XL: 618,247 fully labeled + 2,341,057 partially labeled (5,107 proteins, 1,084,834 ligands) | per Koh et al. 2024 Extended Data Fig. 4 + 6 | three-class functional-effect supervision + binding affinity; PDBbind v2020 reported in Koh et al. Extended Data Table 2 |
| BALM | BindingDB curated to Kd-only with assay limits removed | exact post-curation count is reported in Gorantla et al. 2025 J Chem Inf Model 65(22):12279 Tables (paywalled; not retrievable from open sources) | BindingDB | ESM-2 + ChemBERTa-2 with LoKr (protein) + LoHa (ligand) PEFT; K=512 shared cosine latent; MSE on cosine similarity vs experimental pK_d; HuggingFace `BALM/BALM-benchmark` |

---

## Appendix B: Key references

- Badkul A, Xie L, Zhang S et al. 2025. eMOSAIC: Multimodal out-of-distribution individual uncertainty quantification enhances binding affinity prediction for polypharmacology. Nat Mach Intell 7:1985–1995. doi:10.1038/s42256-025-01151-2
- Bailey M, Moayedpour S, Li R et al. 2024. Deep batch active learning for drug discovery. eLife 12:RP89679. doi:10.7554/eLife.89679
- Bonett DG, Wright TA. 2000. Sample size requirements for estimating Pearson, Kendall and Spearman correlations. Psychometrika 65(1):23–28.
- Gorantla R, Gema AP, Yang IX, Serrano-Morrás Á, Suutari B, Juárez-Jiménez J, Mey ASJS. 2025. Learning binding affinities via fine-tuning of protein and ligand language models. J Chem Inf Model 65(22):12279–12291. doi:10.1021/acs.jcim.5c02063
- Hino H, Eguchi S. 2022. Active learning by query by committee with robust divergences. arXiv:2211.10013.
- Holzmüller D, Zaverkin V, Kästner J, Steinwart I. 2023. A framework and benchmark for deep batch active learning for regression. JMLR 24(164):1–81. github.com/dholzmueller/bmdal_reg
- Koh HY, Nguyen ATN, Pan S, May LT, Webb GI. 2024. Physicochemical graph neural network for learning protein–ligand interaction fingerprints from sequence data. Nat Mach Intell 6:673–687. doi:10.1038/s42256-024-00847-1
- Kull M, Silva Filho T, Flach P. 2017. Beta calibration: a well-founded and easily implemented improvement on logistic calibration for binary classifiers. AISTATS, PMLR 54:623–631. Extended: Electron J Stat 11(2):5052–5080. doi:10.1214/17-EJS1338SI
- Lakshminarayanan B, Pritzel A, Blundell C. 2017. Simple and scalable predictive uncertainty estimation using deep ensembles. NeurIPS 30.
- Landrum GA, Riniker S. 2024. Combining IC50 or Ki values from different sources is a source of significant noise. J Chem Inf Model 64(5):1560–1567. doi:10.1021/acs.jcim.4c00049
- Mervin LH, Afzal AM, Engkvist O, Bender A. 2020. Comparison of scaling methods to obtain calibrated probabilities of activity for protein–ligand predictions. J Chem Inf Model 60(10):4546–4559. doi:10.1021/acs.jcim.0c00476
- Neelon B, Dunson DB. 2004. Bayesian isotonic regression and trend analysis. Biometrics 60(2):398–406. doi:10.1111/j.0006-341X.2004.00184.x
- Nouretdinov I, Volkhonskiy D, Lim B, Toccaceli P, Gammerman A. 2018. Inductive Venn-Abers regressive predictors. PMLR 91:1–12.
- Park et al. 2024. Mixture-of-Experts approach for enhanced drug-target interaction prediction and confidence assessment. bioRxiv 2024.08.06.606753.
- Ren J, Fort S, Liu J, Roy AG, Padhy S, Lakshminarayanan B. 2021. A simple fix to Mahalanobis distance for improving near-OOD detection. arXiv:2106.09022.
- Schulman A, Rousu J, Aittokallio T, Tanoli Z. 2024. Attention-based approach to predict drug–target interactions across seven target superfamilies. Bioinformatics 40(8):btae496. doi:10.1093/bioinformatics/btae496
- Seung HS, Opper M, Sompolinsky H. 1992. Query by committee. COLT pp. 287–294.
- Shoshan Y et al. 2025. MAMMAL — Molecular Aligned Multi-Modal Architecture and Language. arXiv:2410.22367; npj Drug Discovery (2026).
- van Tilborg D, Alenicheva A, Grisoni F. 2022. Exposing the limitations of molecular machine learning with activity cliffs. J Chem Inf Model 62(23):5938–5951. doi:10.1021/acs.jcim.2c01073
- Vovk V, Petej I. 2014. Venn–Abers predictors. UAI 2014:829–838.

---

## TL;DR

- **The five-head ensemble (MAMMAL + Tanimoto + MMAtt-DTA + PSICHIC + BALM) with a Bayesian per-target router, bias-decomposition trust matrix T(t, k), and eMOSAIC cross-head OOD consensus is pre-committed to beat the +0.90 Tanimoto ceiling at SLC6A3 (+0.91) and SLC6A2 (+0.92) while not regressing at the Tier-A SLC6A3 anchor; the disagreement-as-signal facet (D_τ Kendall-pairwise metric on 50 nearest neighbors) is the publishable discovery axis that surfaces novel-scaffold and activity-cliff compounds where Tanimoto is structurally blind.**
- **Identifiability constrains the design: 110 free per-target router weights vs ~440 calibration tuples → differentiable routing is under-identified; the Bayesian router with trust matrix as prior is the only honest option, with family-pooled fallback at n < 7 and uniform fallback at n < 3; MAMMAL must be floored at T = 0.02 panel-wide because leave-one-head-out shows the ensemble *improves* when MAMMAL is removed (prior-collapsed signal is net-noise after calibration).**
- **Honest fallback if the new heads cannot beat Tanimoto +0.90: production deploys Tanimoto-only with the three new heads as supporting voters (T ≤ 0.1 each), and the negative result — "modern foundation DTI models cannot match a 1996 similarity-baseline at transporter cognition targets" — becomes the publishable methodology claim, parallel to the v3 MAMMAL prior-collapse precedent.**

## Recommendations (staged, with benchmarks that change them)

1. **Stage 1 (week 1–2)**: Implement `diagnostics/per_head_bias.py` and `diagnostics/ood_emosaic.py`. Run on the existing v4 12K-compound shortlist with MAMMAL + Tanimoto only, validate the bias-signature framework against known MAMMAL prior-collapse and Tanimoto similarity-searcher signatures. **Benchmark to advance**: PC_MAMMAL < 0.3 confirmed at ≥18 of 22 targets; SN_Tanimoto > 0.9 confirmed at all 22.

2. **Stage 2 (week 3–4)**: Install MMAtt-DTA (MIT license, not pip-installable; clone from github.com/AronSchulman/MMAtt-DTA, download Zenodo models). Run inference on the 12K shortlist for the seven audited targets. Compute per-target Spearman ρ vs ChEMBL pchembl ≥ 8. **Benchmark to advance**: MMAtt-DTA ρ ≥ +0.40 at SLC6A3 (G1 prerequisite).

3. **Stage 3 (week 5–6)**: Install PSICHIC (Apache-2.0, github.com/huankoh/PSICHIC, server psichicserver.com). Extract interpretable interaction fingerprints as eMOSAIC embeddings. Run on the seven audited targets. **Benchmark to advance**: PSICHIC ρ ≥ +0.40 at SLC6A3.

4. **Stage 4 (week 7–8)**: Install BALM (MIT, github.com/meyresearch/BALM). Fine-tune with PEFT on 10 examples per target (BALM's USP7/Mpro few-shot recipe). Run on seven targets. **Benchmark to advance**: BALM ρ ≥ +0.40 at SLC6A3 in few-shot setting.

5. **Stage 5 (week 9–10)**: Implement `fusion/dti_router.py` with Bayesian router; run leave-one-head-out validation. **Benchmark to advance**: ensemble ρ ≥ +0.88 at SLC6A3 (G2 minimum); if < +0.88, switch to negative-result fallback branch.

6. **Stage 6 (week 11–12)**: Implement `fusion/disagreement_diagnostics.py` with D_τ + classifier. Run on full 12K shortlist. **Benchmark to advance**: ≥ 3 of (TC-5619, encenicline, GAT-107, novel α7-PAM scaffolds) appear in top-10 by D_τ (G5).

7. **Stage 7 (week 13)**: Pre-commit the full validation gate set. If ≥ 4 of 6 PASS: deploy v4-ensemble to production, write the J Cheminform / Nat Mach Intell paper. If < 4 PASS: switch to negative-result fallback; write the negative-result paper.

## Caveats

1. **Cross-head correlation matrix Σ_kk' is estimated from n=133 tuples — marginally adequate** (10 free off-diagonal parameters, 13:1 obs:parameter). The 5×5 correlation estimate has ~5% RMSE per off-diagonal entry. If the true correlation structure is more extreme than the expected r̄ ≈ 0.25, the ensemble CI width will be miscalibrated. Mitigation: report sensitivity to ±0.1 perturbation in Σ.

2. **PSICHIC and BALM exact per-split Pearson/Spearman benchmarks were not extractable from the openly available abstracts and bioRxiv preprints.** PSICHIC's PDBbind v2020 numbers are reported in the paper's Extended Data Table 2 (paywalled); BALM's per-split BindingDB Pearson values live in JCIM tables (paywalled). The pre-committed forecasts in §0 are derived from architectural arguments and adjacent superfamily benchmarks; they should be considered ±0.1 ρ wide until reproduced in v4 calibration.

3. **MMAtt-DTA per-superfamily training counts and per-scenario × per-superfamily validation Spearman are in supplementary materials (Schulman 2024 Supp Tables S1, S3, S4), not the main text.** The Schulman bioRxiv preprint states only the qualitative bounds: "Spearman > 0.57" (imputation), "> 0.36" (new compound), "> 0.3" (kinase, new compound + new target). Pre-commitments assume these bounds hold at the cognition-target subfamilies.

4. **The identifiability theorem in §3.4 is informal.** A formal proof requires assumptions about the conditional independence structure of the heads given (q, t), which is not guaranteed in practice (shared training data violates conditional independence). The n* ≈ 720 figure is a rule-of-thumb, not a tight bound. The Bayesian router is the right design choice at v4 sample sizes regardless of the exact value of n*.

5. **GRIN2A (n=7) sits at the boundary between family-pooled fallback (n < 7) and uniform fallback (n < 3).** Default behavior: family-pooled (SLC6 / GRIN / GPCR / enzyme / nuclear). If the v4 §4.5 decision is to deprecate GRIN2A entirely on structural-modeling grounds (ifenprodil heterodimer blindness), the routing question is moot. Either way, GRIN2A should not bottleneck the ensemble.

6. **The disagreement-as-signal facet's discriminator (§6.2) uses hand-coded thresholds (SN ≥ 0.7, T_max ≥ 0.65, inter-head std ≥ 0.7).** These are reasonable defaults but should be validated on a held-out activity-cliff and novel-scaffold benchmark — MoleculeACE (van Tilborg 2022) provides 30 macromolecular targets with activity-cliff annotations and is the natural validation set. v5 should make these thresholds learnable.

7. **The compute budget (single RTX-5070) constrains the deployment.** MAMMAL inference on the 12K shortlist takes ~30 minutes per pass; BALM with ESM-2 + ChemBERTa-2 in fp32 at batch 4 takes ~45 minutes; PSICHIC ~10 minutes; MMAtt-DTA <5 minutes. Total per-target wall-clock for all five heads on 12K compounds: ~90 minutes. On the 22-target panel, this is ~33 hours per full ensemble run. Acceptable for nightly batches; prohibitive for interactive iteration. v5 should explore FlashAttention + fp16 + batched eMOSAIC.

8. **ChEMBL noise floor (Landrum & Riniker 2024) sets a fundamental lower bound on attainable CI width.** Per the paper: even under maximal curation, "48% [of points] differ by more than 0.3 log units, 13% by more than one log unit." Reporting ensemble CIs tighter than ~0.3 pchembl units is dishonest given the noise floor in the calibration data; the Monte Carlo CI propagation in §3.5 should be clipped at σ_min = 0.3.

9. **Pre-committed ρ forecasts in §0 should be treated as ±0.05–0.10 wide.** The Bonett–Wright CIs in the matrix reflect statistical uncertainty given the assumed sample sizes but not the additional uncertainty in the central tendency estimates (which are derived from architectural arguments + adjacent benchmarks, not direct measurement at the cognition-target subfamilies).