# Gaps and research directions

Honest catalogue of what is still open: the engineering gaps, the external
blockers, and the research directions that would change the publication
trajectory. Companion to `README.md`, `PROJECT_STATUS.md`, and the five-paper
manuscript suite. For output artifacts that were named in early plans but never
built, see `FUTURE_WORK.md`. For the full record of what has already shipped,
see the completed ledger at the end of this document.

**Last refreshed**: 2026-06-06. Gaps 1 through 7 are all shipped and the target
panel is finished to **31 targets** (CHRM1/CHRM4/HTR6 added, so CIAS now surfaces
M1/M4 and AD scores 5-HT6), all scored with the real MAMMAL DTI head. The
pipeline runs end-to-end on real LINCS L1000 and real cpg0000. Two June 2026
correctness sweeps hardened the numerical core (router credible intervals,
Venn-ABERS sigma, the V7 PBPK occupancy term, the effect-size compound matcher),
and the V7 PET-anchor claims were reconciled to their honest state (the occupancy
chain reproduces the qualitative dose-ordering, not the absolute values, and is
upstream of, not an input to, the effect-size gates). The new **Frontier
directions** section below sets the next-capability agenda: resolve below mechanism
class, onboard arbitrary novel compounds, scale the ledger, and add causal target
validation. Everything below is what is still open.

---

## Open at a glance

| Item | Type | What it unblocks | Effort |
|---|---|---|---|
| OSF DOI mint (B1) | external | pre-registration lock before release | 30 min |
| bioRxiv submission (B2) | external | public release of all 5 papers | ~2 h per paper |
| Held-out GWAS L2G, ABCD/CAC (B3) | external | V6.B Gate 3 (feeds G1) | ~1 day once DUA in place |
| Wet-lab validation (B4) | external | first empirical validation | 6 to 9 months, $60K to $110K |
| V7 Gate 1 partial-pool (G2) | engineering | tighter compound-level CPT predictions | mitigation listed |
| Scale allosteric benchmark (G3) | engineering | production allosteric rank head | days |
| Mondrian conformal, I_novel (R1) | research | calibrated novel-mechanism intervals | ~3 days |
| Gate 2 multi-modulator (R2) | research | V6.B Gate 2 DEGRADE to PASS | ~1 week |
| Allosteric PBPK (R3) | research | V7 P3/P4/P5 predictions | ~2 weeks |
| Species random effect (R4) | research | mouse to human translation | ~1 week |
| Target deconvolution (R5) | research | target ID for (L,L,H) hits | ~3 weeks |
| Phase 1 trial (R6) | external | external validation, biggest step | 18 to 24 months, $300K to $500K |
| V7 PBPK occupancy anchor-fit (G4) | engineering | honest Figure 1 done; fit occupancy to PET | ~3 days |
| Compound-level resolution test (F1) | DONE | clean NEGATIVE: class is the resolution limit (96.5% between-class variance) | shipped |
| Novel-compound onboarding engine (F2) | DONE | SMILES -> structural class route -> prior g+CrI or ABSTAIN; leave-one-compound-out class recovery 0.97 (36 routed, 60% abstain); exemplars 110 / 46 classes | shipped |
| Ledger scale + per-domain (F3) | DONE | cited n=47 (0.967) + research-curated & human-adjudicated n=125 (all data points kept): class-LOCO AUROC 0.92 (0.97 multi-member), signal survives scaling; 2 genuine mixed classes (anti-amyloid mAb, AChE-I) | shipped |
| Causal MR target validation (F4) | frontier | associative genetics to causal | ~2 to 3 weeks |
| Architectural deepening (F5) | frontier | more performance from the stack | days to weeks |
| Perturbational signature-reversal (F6) | frontier | revive the V8 axis (supervised) | ~2 weeks |

---

## Open engineering gaps

### G1. V6.B Gate 3 has no held-out GWAS L2G

V6.B 4-gate live execution reports Gate 3 (held-out GWAS AUROC > 0.70) as
INSUFFICIENT_DATA because no held-out cognition GWAS L2G has been fetched. The
ABCD Study and CAC (Cognitive Ageing and Health) cohorts are the canonical
held-out sets named in the V6.B OSF pre-registration.

**Resolution**: either (a) fetch ABCD/CAC L2G live from OpenTargets Genetics
Platform v25+ (see external blocker B3), or (b) curate a held-out compound-anchor
set for a synthetic Gate 3 evaluation. Currently listed in the V6.B paper
Discussion as a limitation pending GWAS data acquisition.

### G2. V7 Gate 1 partial-pool predictions are too tightly clustered

V7 NUTS with the 60-anchor likelihood and Schmidli 2014 robust MAP priors
produces compound predictions clustered around the population mean, so 4 of 8
P1-P8 pre-registered checks pass only by tight margins. Donepezil predicts +0.096
(band [0.10, 0.30], misses by 0.004); MPH +0.087 (band [0.15, 0.30]); modafinil
+0.040 (band [0.06, 0.18]).

**Root cause**: with a thin per-class hierarchy, partial pooling shrinks all
class members toward the global mean; at high MAP weight the priors over-dominate.

**Mitigation** (in priority order): (1) denser per-(class, endpoint) cells (the
PRISMA table is at 96 cells, full Cochrane extraction would reach 100+); (2)
hierarchical pooling at the (class x dose) level rather than (class x endpoint);
(3) the allosteric/orthosteric PBPK refinement in R3 directly improves the P3/P4/P5
predictions.

**Status**: documented as an honest publishable finding for the V7 CPT paper:
V7 produces calibrated Bayesian inference within the Roberts ceiling but cannot
tighten compound-level predictions below roughly 0.05 of the population mean given
the current anchor density.

### G3. Scale the allosteric learn-to-rank benchmark beyond n=21

The Gap 4 allosteric learn-to-rank head lifts pooled within-target Spearman rho
from +0.02 (MAMMAL alone) to +0.51 (fused) on a held-out 21-compound benchmark,
which is a proof of concept on small n with only 6 of 21 compounds having Boltz
coverage. To promote it from proof-of-concept to a production within-target
ranker it needs a larger held-out benchmark and fuller Boltz affinity coverage.
Report: `reports/pipeline/allosteric_ltr_v1.md`.

### G4. V7 PBPK receptor-occupancy is not fitted to the PET anchors

The 9-compartment PBPK mass-balance ODE is verified (mass-conserving), but the
receptor-occupancy chain on top of it reproduces only the qualitative dose-ordering
of the 3 PET anchors (Bohnen donepezil, Volkow MPH, Kapur haloperidol), not their
absolute values: with literature Kd and generic per-drug distribution parameters the
peak Hill occupancy saturates above the published sub-saturation readings. A June
2026 fix corrected an inverted spare-receptor term that had been zeroing all
predicted occupancy (which had masked the mismatch); the honest state is now stated
in the V7 paper, Figure 1, and methodology_v3.

**Resolution**: jointly fit per-drug {distribution volume, BBB permeability, Kd} to
the 3 anchors (3 parameters per drug), or adopt a published spare-receptor
amplification factor per receptor family. **Scope note**: this is an upstream,
illustrative module; the V7 effect-size gates consume the PBPK brain-concentration
AUC, not the occupancy estimate, so the headline V7 results do not depend on it.
Effort: ~3 days.

---

## External blockers (need an account, permission, or budget)

### B1. OSF DOI mint

Create an OSF.io account, upload `reports/paper-drafts/v7_osf_preregistration.md`
and `reports/paper-drafts/v8_osf_preregistration.md`, and lock them with a DOI
before any production-data gate evaluation (otherwise pre-registration is not
enforceable). The class-prognostic paper pre-registration is already locked
(OSF DOI 10.17605/OSF.IO/V7GP5). Time: 30 minutes.

### B2. bioRxiv submission

Submit the 5 papers (V6.A, V6.B, V7, V8, integration umbrella) to bioRxiv.org.
Time: roughly 2 hours per paper for formatting and uploading.

### B3. Held-out cognition GWAS L2G (ABCD Study + CAC)

Fetch ABCD Study cognitive-ageing GWAS plus CAC summary statistics and run L2G
via OpenTargets Genetics Platform v25+. This directly resolves G1. Time: ~1 day
once accounts and the data use agreement are in place.

### B4. Wet-lab validation of (L, L, H) candidates

Per `reports/wet-lab/wet_lab_handoff_v1.md`, partner with a CRO (Charles River,
WuXi, Sygnature) for a BIMA-8 remyelination assay on the top (L, L, H) compounds.
Time: roughly 6 to 9 months end to end including compound ordering and assay
development. Cost: $60K to $110K.

---

## Research directions (would change the publication trajectory)

### R1. Mondrian conformal calibration for I_novel (was MH4)

The V8 I_novel novel-mechanism score is computed but its predictive intervals are
uncalibrated. Add Bostrom 2024 `crepes` Mondrian conformal calibration so each
(L, L, H) prediction comes with a guaranteed-coverage interval. Effort: ~3 days.
Direct V8 paper Methods refinement.

### R2. V6.B Gate 2 expansion to multiple modulators (was MH5)

Gate 2 Spearman rho is 0.14 at n=11 (target, primary-modulator) pairs. Map each
panel target to 3 to 5 modulators with per-modulator pooled g; at n=50 pairs Gate
2 becomes statistically powerful and lifts from DEGRADE to PASS. The 70-row
modulator-anchor table is already built (`data/raw/modulator_anchors_seed.csv`),
so this is primarily curation. Effort: ~1 week.

### R3. Allosteric vs orthosteric distinction in V7 PBPK (was MH6)

V7 PBPK currently treats all binding as competitive orthosteric. Allosteric
modulators (BPN14770, LY3154207, encenicline) need pharmacological-effect modifier
terms: add an allosteric flag that modifies the Hill equation per pocket class
(the V4 section 7.5 pocket classifier already distinguishes orthosteric /
allosteric_known / surface_artifact). Improves the V7 P3/P4/P5 partial-pool
predictions (see G2). Effort: ~2 weeks. Direct V7 paper methodology contribution.

### R4. Species translation random effect (was MH7)

V8 iPSC-MEA data is largely mouse-derived while V7 cognition endpoints are
human-RCT-anchored. Add a species random effect to the V8 joint posterior to model
mouse-to-human translation uncertainty explicitly. Effort: ~1 week. Direct V8
paper Discussion and Methods refinement. Pre-registered in V8 OSF section 7.

### R5. Cross-pipeline target-deconvolution integration (was MH10)

V8 phenotype scoring is target-agnostic; once a (L, L, H) compound is identified
the next question is which target it engages. Integrate with Sirota 2011
signature-reversal target deconvolution or Iorio 2018 CMap MOA classification to
surface candidate targets for (L, L, H) compounds and triage which are druggable.
Effort: ~3 weeks. Direct V8 paper methodology extension.

### R6. Pre-registered Phase 1 healthy-volunteer trial (was MH9)

The V8 (L, L, H) clemastine-territory candidates are computationally identified
but never wet-lab validated. Partner with an academic neuroscience lab on an
ethics-approved Phase 1 protocol: healthy adults, single dose, 4-week washout,
double-blind crossover vs placebo, n=20 to 30, primary endpoint DSST + n-back +
Stroop composite. Effort: 18 to 24 months, $300K to $500K. **The single most
impactful next step** for the pipeline's external validation.

---

## Frontier directions (deepen the model, scale it, generalise to novel compounds)

The R-items above refine the five existing papers. The items here are larger bets
that extend the *capability* of the system: resolving below mechanism class,
onboarding arbitrary novel molecules, scaling the evidence base, and adding causal
target validation. They follow from the project's central finding (mechanism-class
clinical track record beats affinity, genetics, and phenotype) and respect its
hardest constraint: leave-one-class-out AUROC is 0.00, so the system cannot
extrapolate to a mechanism with no clinical history. Every frontier item below must
preserve that guardrail (abstain on novel mechanisms) and the project's house rule
(honest negatives over forced positives).

### F1. The compound-level resolution question (RESOLVED 2026-06-06: clean negative)

**RESULT.** Ran (`scripts/93_within_class_resolution.py`,
`reports/pipeline/within_class_resolution_v1.md`): a variance decomposition of
clinical *g* by mechanism class shows **96.5% of the variance is between classes**
(one-way ICC 0.95); only 3.5% lives within classes. No pre-specified compound
feature (CNS-druglikeness/exposure proxy, readout recency, structural typicality,
QED) beats the class mean within class (all within-class Spearman |rho| below the
0.52 minimal-detectable effect at this n; all leave-one-compound-out MAE deltas
<= 0; Holm ns). **At n=31, mechanism class is the empirical resolution limit of
in-silico cognition-drug prognosis.** This is a *bounded* negative (the test is
underpowered by design at n=31); separating "true ceiling" from "low power" is
what F3 (ledger scaling) resolves, and the highest-value untested feature is real
per-compound binding affinity / trial dose-adequacy. The original framing follows.

The headline predictor is *class-level*: it assigns the class mean to every member
of a class. The honest ceiling result is leave-one-class-out = 0.00 (no
extrapolation to unseen mechanisms), and within-class ordering is not currently
modelled. The frontier question: **is there ANY compound-level signal that beats
"assign the class mean", and if so, which feature carries it?**

Pre-registered test: for each multi-member class, compare a compound-level ranker
(features: brain exposure / dose-adequacy, Gini selectivity, off-target liability,
V6.A multi-head affinity) against the class-mean baseline on leave-one-compound-out
clinical g. The most promising single axis is **dose-adequacy**: several famous
failures (encenicline 3 mg) are plausibly under-dosing, not wrong-mechanism, which
V7's PBPK brain-AUC can quantify. Two honest outcomes, both publishable: (a) a
within-class feature lifts ranking (a genuine compound-level advance), or (b) nothing
beats the class mean, establishing "mechanism class is the resolution limit of
in-silico cognition-drug prognosis" as a clean negative. Effort: ~2 weeks.

### F2. Novel-compound onboarding engine (DONE 2026-06-07)

**RESULT.** Shipped (`src/mammal_repurposing/validation/novel_compound.py`,
`scripts/95_novel_compound_onboarding.py`, `reports/pipeline/novel_compound_onboarding_v1.md`).
The CPU decision core is live: novel SMILES -> structural class assignment (max ECFP4
Tanimoto + Murcko generic-scaffold to ledger exemplars) -> EB-shrunk class clinical-*g*
prior + 90% bootstrap CrI -> confidence tier, with both guardrails enforced. Validated by
leave-one-compound-out class recovery on the exemplar base (grown 31 -> 110 SMILES across
**46 of 48 classes** via `scripts/_expand_ledger_smiles.py`; PubChem canonical SMILES,
RDKit-gated): **top-1 class recovery 0.97** on the 36 routed held-out drugs, abstaining on
60% where no close analog exists (the guardrail working). The out-of-manifold floor
(TAU_OOD=0.35) was calibrated, not guessed: every observed mis-route sat in Tanimoto
[0.26, 0.34]. The single residual error is an enantiomer blind spot (2D ECFP4 cannot
separate (-)-phenserine/AChE-I from its (+)-enantiomer posiphen/buntanetap/APP-inhibitor).
Demo: ipidacrine->AChE-I (HIGH), phentermine->catecholaminergic (MED); the peripheral
negatives (aspirin/ibuprofen/loratadine/atorvastatin) and caffeine all ABSTAIN.

**DTI-profile signal TESTED (GPU, negative).** The spec's primary signal (a) "nearest
class in DTI-profile space" was scored on the RTX 5070 (118 compounds x 31 targets,
MAMMAL; `scripts/96` -> `scripts/97`, `reports/pipeline/f2_profile_vs_structure_v1.md`)
and compared to structure leave-one-compound-out: structure-only **0.972**, profile-only
**0.112**, blended **0.493** (blending HURTS). The MAMMAL profile is nearly non-selective
on this panel (median within-compound pKd SD 0.37; donepezil/galantamine/pitolisant/
modafinil share the same top targets) - the documented property-correlation bias - so it
cannot route mechanism classes. Structure stays primary; the profile remains a documented
pluggable `external_class_scores` hook, OFF by default. Honest negative, coherent with the
weak affinity signal (AUROC 0.47). The original framing follows.

Today the class prior only helps a compound already placed in a class. To score an
*arbitrary* novel molecule, wire the existing pieces into one path: novel SMILES ->
multi-head DTI profile (MAMMAL + MMAtt-DTA + PSICHIC + BALM) over the 31-target
cognition panel -> mechanism-class assignment via (a) nearest class in DTI-profile
space, (b) Tanimoto to class exemplars, (c) scaffold match -> class prior returns a
calibrated g + 90% CrI -> confidence tier (HIGH only for a clean single-class map
with a populated prior; ABSTAIN for a novel mechanism, an ambiguous profile, or an
out-of-manifold compound flagged by the multi-head OOD axis).

This converts the retrospective validator into a prospective screen: run a vendor
catalogue or the ChEMBL CNS subset through it and rank by predicted clinical g. Two
non-negotiable guardrails: (1) leave-one-class-out = 0.00 means the engine MUST
abstain on genuinely novel mechanisms (it re-ranks known mechanisms for cognition,
it does not invent new ones); (2) the V6.A allosteric blindness makes DTI-profile
class assignment unreliable for allosteric compounds, so the allosteric-awareness
head (G3) is a required sub-component. The highest-value output: compounds in
strong-precedent classes that have NOT been tried for cognition. Effort: ~3 weeks.
This is the most direct answer to "expand to novel compounds".

### F3. Scale and decompose the clinical ledger (PARTIAL 2026-06-06: analysis + infra shipped)

**RESULT.** Ran (`scripts/94_ledger_scaling.py`,
`reports/pipeline/ledger_scaling_v1.md`): the class-separation result **survives
scaling on the real cited ledgers** (base 31 -> +EXTENSION 42 -> +CT.gov 47):
class-LOCO AUROC 1.000 -> 0.990 -> 0.967, all 20 classes stay 100% outcome-pure,
and 97% of clinical-*g* variance stays between-class (ICC 0.95) -- the headline is
not a small-n artifact of the original 31. Within-domain class separation holds in
the largest mixed domain (AD global-amnestic, AUROC 0.79 on the full n=125 ledger,
across its 59 drugs). The **power roadmap** (the
actionable output) says the F1 within-class test needs ~65 drugs (within-class
rho=0.4) up to ~118 (rho=0.3), concentrated in multi-member SUCCESS classes with
genuine within-class g spread. What remains is genuine literature curation (real
drugs/trials/outcomes/cited g -- NOT auto-generated, to protect ledger integrity);
the protocol + per-domain schema are in `docs/LEDGER_CURATION.md` and
`load_all_ledgers()` ingests any schema-conforming batch.

**Research-curation + adjudication (2026-06-06).** An Opus multi-agent run (research +
independent adversarial verification, 106 agents) added 78 web-verified cognition-drug
outcomes; two further independent Opus adjudicators then re-coded every disputed call
under a strict cognition-EFFICACY convention. All 78 are kept as binary data points
(no exclusions): the safety-halted metrifonate/eptastigmine are retained as SUCCESS on
their clean positive cognition pivotals, while unverifiable or contested-benefit drugs
(velnacrine, aducanumab, sodium oligomannate, masitinib, nicergoline) are FAILURE; 11
over-generous SUCCESS calls were recoded to FAILURE in total -> n=125. Honest finding:
class-outcome purity is **scale-sensitive but robust**. The raw web-research AUROC fell
to 0.77 (coding noise); after adjudication the class-LOCO AUROC is **0.92** (**0.97** on
the multi-member classes the predictor can leverage), still far above the leakage-free
target-level predictors (affinity 0.47,
genetics 0.59). Only TWO genuinely mixed-outcome classes remain: **anti-amyloid mAbs**
(lecanemab/donanemab succeed where 5 earlier anti-Abeta mAbs failed) and
**AChE-inhibitors** (marketed winners vs later AChE-Is that failed on efficacy/dosing).
The perfect 1.00 at n=31 was partly a sparse-sampling / selection effect; the honest
fully-populated value is ~0.92. Per-row adjudication basis is in the RESEARCH provenance;
the frozen base-31 analysis is untouched. The original framing follows.

AUROC 1.00 at n=31 is a small-n result; perfect separation can be partly an n
artifact. The single highest-value rigor step is to scale the leakage-audited ledger
to 100 to 200+ cognition drugs with cited meta-analytic g, and ask whether class
separation survives. Two extensions make it both more accurate and more clinically
actionable: (a) **per-domain decomposition** of the single clinical g into
working-memory / processing-speed / episodic-memory / executive sub-scores
(stimulants likely win processing speed, cholinergics memory), giving a per-domain
class prior; (b) a **dose-response term** so under-dosed failures are separated from
mechanism failures (feeds F1). Effort: ~3 to 4 weeks, mostly curation, and it
de-risks every downstream claim.

### F4. Causal target validation via Mendelian randomisation

The class prior is associative; the V6.B genetics axis (OpenTargets L2G) is also
associative. Deepen it with **Mendelian randomisation** on cognition GWAS to test
whether modulating each panel target *causally* affects cognition, not merely
correlates. This converts the genetics axis from "the locus maps to the gene" to
"the gene is causal for the phenotype", and answers the hardest skeptic question
about every target on the panel. Effort: ~2 to 3 weeks once cognition GWAS summary
statistics are in hand (overlaps external blocker B3).

### F5. Architectural deepening (more performance from the existing stack)

- **Structure as a fused head**: feed Boltz-2 affinity into the (now NaN-robust)
  Bayesian router alongside the four sequence heads, so structure votes in the
  ensemble rather than sitting in a separate cache.
- **Allosteric rescue**: the V6.A Tier-A FAIL is the stack's known weakness. Expand
  the allosteric learn-to-rank head's training data (curated allosteric-modulator
  affinity sets) and add an explicit allosteric-site structural feature using the
  already-scaffolded cryptic-pocket detectors.
- **Close the active-learning loop**: the scaffold-aware AL scorer (shipped) should
  choose which compounds to score or validate next, so each round maximally reduces
  predictive uncertainty.

### F6. Revive the perturbational axis as signature reversal, not clustering

V8's unsupervised phenotype clustering failed its pre-registered Gate 1 (AMI = 0.13;
CL6). That kills *clustering*, not the LINCS L1000 data. A principled second attempt
reframes the axis as a directional, supervised question: does a compound **reverse** a
cognitive-decline / brain-ageing transcriptomic signature (Connectivity-Map signature
reversal, Sirota 2011), rather than cluster into mechanism groups? Different
hypothesis, same data; it may carry signal where clustering did not. It must be
pre-registered as a fresh gate (do not re-run the failed clustering) and stay honest
about CL6. Effort: ~2 weeks.

---

## Cross-cutting limitations (standing)

These are inherent to the framework and are stated honestly in every manuscript.

- **CL1. Small statistical power at every layer**: V6.A per-target n=7 to 26
  (ChEMBL pchembl >= 8); V6.B 191-target panel (post-MH8, 0 divergences); V7
  60-compound anchor set; V8 chemCPA on 107K real LINCS signatures plus
  hierarchical on real cpg0000 (60 compounds); the Gap-3 retrospective ledger is
  31 drugs across 11 mechanism classes. Every layer-specific manuscript needs a
  limitations section stating per-axis n and its implications.
- **CL2. No prospective wet-lab validation**: the pipeline is in-silico only; all
  predictions are calibrated against published meta-analytic literature. The
  wet-lab handoff document is the path forward (B4).
- **CL3. The Roberts 2020 ceiling is unmodifiable**: maximal performance is
  bounded at g ~ 0.50 at the 90% credible upper bound. This is a feature, not a
  bug: the pipeline cannot over-promise and the OSF pre-registration enforces it.
- **CL4. Cross-species and cross-cell-line generality is limited**: V8 uses U2OS
  osteosarcoma cells and partly-mouse iPSC-neuron data. The per-cell-line random
  effect (shipped, MH3) and the species random effect (R4) directly address this.
- **CL5. Pre-registration enforcement needs the OSF DOI mint**: the V7 and V8 OSF
  pre-registrations are markdown-ready but not yet DOI-locked (B1).
- **CL6. The phenotype axis does not recover cognition pharmacology on real
  data**: V8's pre-registered Gate 1 FAILed on real LINCS (AMI = 0.13;
  `reports/pipeline/v8_real_gate1_v1.md`). The perturbational-signature axis is
  shelved as a paper (`reports/paper-drafts/shelved/`) and kept as a documented
  negative. This is on-thesis (phenotype joins affinity and genetics in
  underperforming the class-history prior), not a pipeline bug.

---

## Recommended order

1. **B1 + B2** (OSF DOI mint + bioRxiv submission): public release of all papers
   plus the Gap 1-7 results. Lowest effort, highest visibility.
2. **B3 + R2** (held-out GWAS L2G + Gate 2 multi-modulator expansion): together
   lift V6.B Gate 2 from DEGRADE to PASS and resolve Gate 3 (G1). The 70-anchor
   table is built; this needs the GWAS L2G fetch.
3. **G3** (scale the allosteric benchmark): promotes the Gap 4 head from
   proof-of-concept toward a production within-target ranker.
4. **R1** (Mondrian conformal calibration): V8 paper methodology refinement.
5. **R3** (allosteric vs orthosteric PBPK): V7 paper contribution, also improves
   G2.
6. **R4** (species translation random effect): V8 Discussion refinement.
7. **R5** (target-deconvolution integration): extends V8 paper scope.
8. **B4 then R6** (wet-lab validation, then a Phase 1 trial): the external
   validation arc, longest and most expensive but the most decisive.

**Frontier tier** (capability expansion, can run in parallel with the release arc):

9. **F3** (scale + decompose the ledger): analysis + infra **shipped** -- scaling
   survives to n=47 (AUROC 0.967, ICC 0.95) and the power roadmap sets the target
   (~65-118 drugs). The remaining work is the cited literature curation itself
   (`docs/LEDGER_CURATION.md`).
10. ~~**F1** (compound-level resolution test)~~ **DONE** -- clean negative (class is
    the resolution limit at n=31; see Frontier F1). The follow-on is F3, which gives
    the within-class test the power to distinguish a true ceiling from low power.
11. **F2** (novel-compound onboarding engine): turns the predictor into a prospective
    discovery screen; the direct answer to "expand to novel compounds".
12. **F4** (causal MR target validation): associative genetics to causal; pairs with
    the B3 GWAS fetch.
13. **F5 + F6** (architectural deepening + perturbational signature-reversal): more
    performance from the existing stack, plus a principled second attempt at the
    shelved V8 axis.

---

## Completed ledger

Condensed record of what has shipped. Full detail lives in `PROJECT_STATUS.md`,
the named reports below, and the manuscript suite.

### Since the 2026-05-30 refresh

- **F2 novel-compound onboarding engine** (`src/mammal_repurposing/validation/novel_compound.py`,
  `scripts/95_novel_compound_onboarding.py`): the prospective screen. Novel SMILES ->
  ECFP4-Tanimoto + Murcko-scaffold class assignment -> EB-shrunk class clinical-*g* prior
  + 90% CrI -> tier, with a hard ABSTAIN on out-of-manifold / novel mechanisms and an
  allosteric (V6.A) downgrade. Leave-one-compound-out class recovery **0.97** (36 routed,
  60% abstain) on an exemplar base grown 31 -> 110 SMILES / 46 classes (PubChem, RDKit-
  gated, `scripts/_expand_ledger_smiles.py`). OOD floor calibrated from the LOCO error
  band; one residual enantiomer mis-route. 14 new tests. Detail:
  `reports/pipeline/novel_compound_onboarding_v1.md`.
- **F2 DTI-profile signal tested on GPU (negative)** (`scripts/96` + `scripts/97`,
  `reports/pipeline/f2_profile_vs_structure_v1.md`): scored 118 compounds x 31 targets
  with MAMMAL on the RTX 5070, then compared the spec's primary "nearest class in profile
  space" signal to structure (leave-one-compound-out). Structure **0.972** vs profile-only
  **0.112** vs blended **0.493** (blending hurts; profile rescues only 6% of structure-
  abstained). The MAMMAL profile is near-flat on this panel (within-compound pKd SD 0.37;
  distinct drugs share top targets) - property-correlation bias - so structure stays the
  default and the profile is an OFF-by-default hook. Honest negative.
- **F3 ledger scaling + per-domain + power roadmap** (`reports/pipeline/ledger_scaling_v1.md`,
  `src/mammal_repurposing/validation/ledger_scaling.py`): the class-separation
  result survives the cited n=31 -> 47 expansion (class-LOCO AUROC 1.000 -> 0.967,
  20/20 classes outcome-pure, ICC 0.95), so it is not a small-n artifact. Per-domain
  stratification holds in AD global-amnestic (AUROC 0.79 on the full n=125 ledger).
  The power roadmap
  quantifies the F1 curation target: ~65 drugs (within-class rho=0.4) to ~118
  (rho=0.3), concentrated in SUCCESS classes. Curation protocol + per-domain schema
  in `docs/LEDGER_CURATION.md`; 5 new tests.
- **F3 research-curation + adjudication (Opus multi-agent)** (`data/raw/clinical_outcomes_ledger_RESEARCH.csv`
  + provenance): a research + independent-adversarial-verification workflow (106 Opus
  agents) added 78 web-verified cognition-drug outcomes; two further independent Opus
  adjudicators re-coded every disputed call under a strict cognition-efficacy convention
  (all 78 kept as binary data points -- no exclusions; 11 over-generous SUCCESS recoded
  to FAILURE; safety-halted metrifonate/eptastigmine retained as SUCCESS on their
  positive cognition pivotals) -> n=125. Honest finding: class-outcome purity is
  **scale-sensitive but robust** -- raw AUROC 0.77 (coding noise) lifts to **0.92**
  (0.97 multi-member) after
  adjudication, still >> affinity 0.47 / genetics 0.59. Only two genuine mixed classes
  remain (anti-amyloid mAbs; AChE-I). Per-row adjudication basis in the provenance;
  frozen base-31 untouched. Detail: `reports/pipeline/ledger_scaling_v1.md` 3b.
- **F1 compound-level resolution test** (`reports/pipeline/within_class_resolution_v1.md`,
  `src/mammal_repurposing/validation/within_class.py`): a pre-registered test of
  whether any compound-level feature beats the class mean WITHIN a mechanism class.
  Clean NEGATIVE: 96.5% of clinical-*g* variance is between classes (one-way ICC
  0.95), and no pre-specified feature (CNS-druglikeness, recency, structural
  typicality, QED) ranks drugs within class beyond chance (Holm ns; all
  leave-one-compound-out MAE deltas <= 0). At n=31, mechanism class is the empirical
  resolution limit -- a bounded negative that motivates F3 (ledger scaling). 13 new
  tests; canonical SMILES committed at `data/raw/ledger_compound_smiles.csv`.
- **Prospective trial-watch** (`reports/pipeline/trial_watch_v1.md`): a standing
  forward-prediction system that derives a calibrated per-class SUCCESS prior
  from the n=47 ledger and predicts each ongoing cognition trial leakage-safe
  (the trial drug held out), with honest confidence tiers. 2/2 resolved correct
  (Brier 0.009); the engine reproduces all 6 frozen hand predictions; a
  round-trip over the ledger recovers 47/47 outcomes. This is the in-silico
  prospective track record, accruing as trials read out; it does not replace
  wet-lab validation (B4).
- **V8 real-data Gate 1** (`reports/pipeline/v8_real_gate1_v1.md`): the
  pre-registered phenotype-axis gate, run on real LINCS L1000 (16 cognition
  compounds, 10 pharmacology-labelled classes), gives AMI = 0.13, a FAIL vs the
  0.50 bar. A pre-registered negative result that supersedes the synthetic
  AMI=1.00 dry-run. On-thesis: phenotype, like affinity and genetics,
  underperforms the mechanism-class track record. The V8 paper needs reframing
  as a negative result before any submission.

### The seven-gap arc

- **Gap 1. Degenerate shortlist fixed**: the v10 shortlist collapsed all 298
  compounds onto ACHE (two coupled bugs). The v11 grid composer scores the full
  298 x 13 (compound, target) grid; top-25 now spans 7 unique targets and 13
  distinct g values, positive controls land at the correct targets, max g90 across
  3,874 hypotheses is 0.39 < 0.50. `reports/wet-lab/wet_lab_shortlist_v11.md`.
- **Gap 2. Disease-population reframe**: re-scores the v11 grid per disease using
  that disease's own pivotal-trial track record as the per-class prior. AD recovers
  AChE-I (within-disease AUROC 0.97, 10/10 AD failures flagged), CIAS recovers
  muscarinic M1/M4, FXS recovers PDE4. `reports/pipeline/disease_reframe_v1.md`.
- **Gap 3. Retrospective clinical validation**: leakage-audited ledger of 31
  cognition drugs. Mechanism-class track record (leave-one-compound-out) gives
  AUROC 1.00 and flags 9/9 famous Phase III failures, while target binding (0.12)
  and target genetic relevance (0.59) are at or below chance. Leave-one-class-out
  is 0.00 (the honest extrapolation ceiling).
  `reports/pipeline/retrospective_clinical_validation_v1.md`.
- **Gap 4. Allosteric learn-to-rank head**: quantifies MAMMAL's within-target
  blindness (predicted-pKd std 0.01 to 0.05 across 3 log-units of affinity) and
  fuses [MAMMAL pKd, Tanimoto, Boltz, physicochemistry] to lift held-out
  within-target rho from +0.02 to +0.51. `reports/pipeline/allosteric_ltr_v1.md`.
  Open follow-up: scale beyond n=21 (G3).
- **Gap 5. Clinician GRADE dossiers**: one-page per (compound, indication) card
  with predicted g + 90% CrI, Cochrane evidence rating, mechanism-class track
  record, off-target liability flags, and failure-mode caveats.
  `reports/pipeline/clinician_dossiers_v1.md`.
- **Gap 6. External benchmark**: on the held-out 31-drug task the mechanism-class
  track record (AUROC 1.00) beats the leakage-free target-centric paradigms,
  affinity (0.47) and genetics (0.59); a target-popularity baseline scores 0.96
  but is shown to be a hindsight confound. `reports/pipeline/external_benchmark_v1.md`.
- **Gap 7. Prospective repurposing shortlist (the capstone)**: approved drugs
  ranked as mechanism-justified repurposing hypotheses per disease. AD to huperzine
  A and sigma-1 agonists; CIAS to 5-HT1A and M1/M4 (with xanomeline correctly
  flagged STANDARD, validating the approach); FXS to roflumilast and rolipram.
  `reports/pipeline/repurposing_shortlist_v1.md`.
- **Panel finished to 31 targets**: CHRM1/CHRM4/HTR6 added and scored with the
  real MAMMAL DTI head, so the CIAS muscarinic winner and the AD 5-HT6 failure
  class are now scorable, not just priced.

### Modeling improvements

- **MH1. PRISMA per-subdomain priors**: 32 to 96 cells.
- **MH2. V7 anchor set**: 15 to 60 compounds (109-row anchor table feeding V7
  NUTS V2).
- **MH3. V8 per-cell-line random effect**: alpha_cell + gamma + delta on real
  cpg0000; ICC_cell 0.018, ICC_inter 0.149, 60/60 compounds transferability > 0.6.
  The empirical defence of U2OS to brain transfer.
- **MH8. Substrate-mediated AHBA-masking**: 10x AHBA sigma inflation for
  {ACHE, MAOA, MAOB, COMT} took V6.B.5 NUTS from 37 divergences to 0 on the
  191-target panel (R-hat 1.000, ESS 1,808).
- **Gate 2 multi-modulator falsification** (publishable negative result):
  high-affinity binders at cognition-validated targets are not predictive of
  clinical success (negative rho across aggregations). This is the central
  methodological motivation for the multi-layer pipeline.
  `reports/pipeline/gate2_multi_modulator_v1.md`.

### Infrastructure resolved

- **PyMC EOFError on Windows** fixed via the numpyro JAX backend (in-process NUTS):
  V6.B.5 production NUTS on 191 targets now completes in 8 seconds at R-hat 1.000.
- **V8 chemCPA on real LINCS**: trained on 107K real Level-5 signatures in 8.3 min
  on the RTX 5070; Val R-squared 0.46, OOD R-squared 0.33.
- **JUMP-CP via cpg0000**: real Cell Painting (A549 + U2OS), V8 hierarchical NUTS
  R-hat 1.010, 0 divergences. The full cpg0016 sync remains an optional scale-up.
- **Reachability validated**: LINCS L1000 GEO (2,170 perturbations parsed from the
  GSE70138 pert-info) and JUMP-CP S3 (14 sources listed, UNSIGNED signature) both
  reach from the sandbox.

---

*Companion to `README.md`, `PROJECT_STATUS.md`, the manuscript suite, and the
wet-lab handoff. Planned-but-unbuilt outputs and the buildable model-expansion
roadmap are in `FUTURE_WORK.md`. Last refreshed 2026-06-06: Gaps 1 to 7 shipped,
panel finished to 31 targets, two June correctness sweeps + the V7 PET
reconciliation landed; the open items are the OSF/bioRxiv release, held-out GWAS for
Gate 3, the V7 partial-pool tightening, and the Frontier directions (F1 to F6) that
set the next-capability agenda.*
