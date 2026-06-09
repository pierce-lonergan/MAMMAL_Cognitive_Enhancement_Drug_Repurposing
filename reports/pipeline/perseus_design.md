# PERSEUS - Persistence-aware Pro-cognition Engine (design + v1 status)

Source: an adversarially-verified Opus multi-agent research synthesis (2026-06-07; 12 web
research lanes each independently verified, + a repo data audit + chief-architect
synthesis; 26 agents). The verification layer refuted 6 fabricated/mis-cited numbers
before they entered the design (appendix). This document is the engineering spec and the
status of the buildable v1.

## v2 revisions (post expert review, 2026-06-07)

An expert review judged the framework strong (~9/10) but the v1 artifact weaker (~5/10):
the L3 substrate was inert on the real data (constant `transient_signaling` because it
keyed off the misrouted structural class), the class-prior artifact had reappeared in
`evidence_design`, and several integrity bugs survived. All were fixed:

1. **L3 fires from MECHANISM, not the structural class** - substrate now comes from the
   curated persistence axis (so L3 and L5 are one coherent call and can no longer
   contradict), and the ordinal is collapsed to **3 honest tiers** -
   `transient < plasticity_window < ablative` - with a separate `self_maintaining` flag.
   Reversible epigenetic/NRF2 chemotypes are CAPABILITY FLAGS only (not promoted), because
   transient engagement does not self-maintain (appendix #6: BDNF->TrkB is self-limiting).
   Only ablative (senolytic / aggregate-clearing) is durable by construction.
2. **evidence_design de-broadcast** - a new weakest tier `class_extrapolation` carries a
   class fact borrowed onto a member with no per-compound study (the anorectics no longer
   inherit methylphenidate's `randomized_discontinuation`).
3. **NULL split** into NULL_SYMPTOMATIC / TESTED_NEGATIVE / EXCLUDE_NOT_COGNITION; the
   NSAIDs are re-bucketed out of "symptomatic" into not-cognition.
4. **L0-mismatch guard** - a `known_mechanism_class` that contradicts the structural route
   (fluoxetine=SSRI, duloxetine=SNRI, rasagiline/selegiline=MAO-B) forces the SYMPTOMATIC
   head to ABSTAIN instead of emitting a wrong-class prior.
5. **Prodrugs** scored as the active species (serdexmethylphenidate -> dexmethylphenidate);
   novel prodrugs not in the curated map are a documented limitation.
6. **Governor fixes** - anti-amyloid mAbs downgraded to `longitudinal_followup` (parallel-
   group cannot reach delayed-start under the engine's own rule); selegiline corrected off
   `delayed_start_rct` (DATATOP/Sano were not delayed-start - only rasagiline/ADAGIO is);
   the fenoprofen-calcium templating leak fixed by neutralised-parent salt dedup.
7. **Salts collapsed** to a neutralised parent before scoring.

**Honest novelty (demoted from the v1 overclaim).** PERSEUS is NOT yet a "validated
computable persistence predictor." It is a **CNS-gated, abstain-by-default, evidence-design-
governed wrapper that reports symptomatic and persistence as separate outputs and refuses to
call anything demonstrated** - a guardrail on MAMMAL/F2 that correctly stops a symptomatic
+0.40 from masquerading as durability. The control panel (13/13) is a set of CONSISTENCY
unit tests, not out-of-sample validation. The genuine evaluation results so far are both
SPECIFICITY measures, in two directions:
- `scripts/101_perseus_eval.py` -> **0 / 15 persistence false positives** (reversible-
  enhancer + persistence-illusion negative panels);
- `scripts/102_persistence_groundtruth_eval.py` -> **0 / 14 over-claims** against a cited
  persistence-DESIGN ledger (`data/raw/persistence_ground_truth.csv`: delayed-start /
  randomized-discontinuation / washout / parallel-group readouts, labelled by what the trial
  actually showed). An over-claim = the verdict asserting MORE durability than the trial-
  design label supports - the directional error that matters for an abstain-by-default
  system. Report: `reports/pipeline/perseus_eval_groundtruth_v1.md`.

The bidirectional eval earned its keep on the first run: it caught **selegiline** over-
claiming CONTESTED. Root cause - DATATOP / Sano 1997 are NOT delayed-start designs; the
time-to-levodopa endpoint is confounded by selegiline's own symptomatic effect (the exact
confound delayed-start was invented to remove). Fix: selegiline demoted CONTESTED ->
tested_negative, leaving rasagiline/ADAGIO and the fluoxetine plasticity window as the only
live threads. Coverage-accuracy over the evidence-design rank holds 1.00 non-over-claim
accuracy at every threshold. The **label budget is ~381 confirmed delayed-start positives**
(1% prior, +/-0.1) before recall is estimable - reported as a first-class deliverable.
Sensitivity / PPV remain unidentifiable without a NON-EMPTY positive persistence ledger + PU
/ leave-one-mechanism-out estimator (the next deliverable). The "first structure-computable
persistence prior" claim is withdrawn until the persistence-target DTI module + that
estimator land.

## Thesis

Treat "helps cognition" and "persists after cessation" as **two orthogonal outputs, never
one score**. The symptomatic question is well-posed (route a SMILES to a mechanism class,
emit the class's validated clinical-g prior, class-LOCO AUROC ~0.92). The persistence
question has a **near-empty positive class** (durable post-cessation cognitive gain in
healthy people is essentially undemonstrated), so we do NOT train a high-capacity
persistence classifier. Instead persistence is a calibrated, falsifiable, **abstain-by-
default multi-gate AND**: a non-null call requires (1) free-brain exposure, (2) a STATE-
changing (not tone-changing) mechanism, and (3) where trial evidence exists, an evidence-
design tier strong enough to support it. The headline output is honest abstention; the
novel capability is flagging a state-changing MECHANISM as a persistence HYPOTHESIS even
with no trial, while refusing to call it demonstrated.

## Layered architecture (status)

| layer | purpose | v1 status |
|---|---|---|
| L0 structure router + applicability gate | route SMILES -> known cognition class or ABSTAIN | SHIPPED (`validation/novel_compound.py`) |
| **L1 free-brain CNS-exposure gate** | PASS/FAIL/ABSTAIN: CNS-MPO-like + permanent-charge/peptide vetoes | **SHIPPED** (`engine/cns_exposure.py`) |
| L2 symptomatic pro-cognition head | class clinical-g prior + CrI + tier (the validated half) | SHIPPED (reuses `novel_compound` / `retrospective`) |
| **L3 mechanism reversibility** | 5-level persistence-substrate ordinal (tone vs state) + structural alerts | **SHIPPED** (`engine/reversibility.py` + 2 cited CSVs) |
| **L5 evidence axis + composer** | curated status + evidence-design tier; AND with abstain-by-default | **SHIPPED** (`engine/perseus.py` + `validation/persistence.py`) |
| L4 permissive-window firewall | cap plasticity-enablers (TrkB/5-HT2A/PNN) to "durable IFF paired training" | v1 via the persistence axis `plasticity_gated` status; full TrkB-TMD / 5-HT2A-access head is FUTURE |

The **persistence-substrate ordinal** (persistence-capability is monotone): `transient_signaling`
(0, reuptake/AChE/orthosteric = TONE) < `durable_transcriptional` (1, TrkB/NRF2/ISR) <
`structural_ecm` (2, PNN/ECM) < `self_propagating_epigenetic` (3, HDAC/DNMT) <
`ablative_cell_population` (4, senolytic / aggregate clearance). State-changing is
necessary-not-sufficient.

## v1 result (control panel + F2 shortlist)

`scripts/100_perseus.py` -> `reports/pipeline/perseus_v1.md` + `perseus_scored.csv`.

- **Control panel 12/12 as expected**: reversible enhancers (methylphenidate/modafinil/
  donepezil) -> symptomatic tier but NULL persistence; misroutes (neostigmine/demecarium/
  distigmine/difelikefalin) -> EXCLUDE_NO_CNS at the L1 gate; HDACi/NRF2 exemplars
  (vorinostat/dimethyl-fumarate/sulforaphane) -> CANDIDATE_MECHANISTIC; out-of-manifold /
  CNS-borderline (caffeine, entinostat) -> ABSTAIN. Nothing is DEMONSTRATED_HEALTHY.
- **F2 shortlist re-scored** (31 compounds): the uniform +0.40 symptomatic prior splits into
  4 CNS-excluded, 13 not-cognition-excluded, 11 null/symptomatic, 1 TESTED_NEGATIVE
  (selegiline, demoted off delayed-start - see eval below), 1 CONTESTED (rasagiline;
  delayed-start), 1 WINDOW_CONDITIONAL (fluoxetine). The persistence head is a model output,
  not a buried caveat.

## Evaluation framework (for a near-empty positive class)

Ordinary held-out AUROC on the persistence label is meaningless (recall is un-estimable as
the positive prior -> 0). The plan: (1) LABELS only from genuine persistence designs
(delayed-start RCT, randomized-discontinuation, washout follow-up); everything else is
UNLABELED, not negative. (2) PU / one-class with an externally-supplied class-prior
interval, sanity-bounded by the empirical ~0 base rate. (3) leave-one-MECHANISM-out CV
(grouped, not random) to avoid chemotype memorization and leave-one-out distributional
leakage. (4) prior-corrected calibration (Saerens) + conformal/Venn-ABERS. (5) **abstention
is the headline metric** - a coverage-accuracy curve vs two negative-control panels:
reversible enhancers (must score non-persistent) and a "persistence-illusion" panel
(single washout signal that LATER failed a definitive trial, e.g. exenatide-PD). (6)
falsifiable per-compound hooks (BDNF Val66Met genotype-dependent efficacy; a scheduling
prediction for plasticity-enablers). Deliverable framing: "at X% abstention, Y PPV at the
realistic ~1% prior" - which no current repurposing predictor reports for persistence.

## Novelty vs SOTA

- Two orthogonal heads reported separately (SOTA - BRDKRM, Cheng/CMap, DRIAD, TxGNN,
  aging-clocks - collapse symptomatic and disease-modifying into one score).
- First structure-/target-COMPUTABLE persistence prior (BRDKRM is corpus-bound to KG nodes;
  PERSEUS scores an arbitrary SMILES).
- Abstention as the headline output, base-rate-aware, with a coverage-accuracy curve.
- Permissive-vs-instructive firewall (L4): plasticity-window-opening is directionally
  neutral and durability is conditional on paired experience.
- Persistence-substrate taxonomy (state vs tone) as the backbone - re-points the question
  from "does it reverse a signature" to "does the pathway change cell STATE".
- Evidence-design-tiered confidence governor: parallel-group on-drug benefit contributes
  ZERO to persistence confidence; only delayed-start/discontinuation can raise it.

## v2.8 four-frontier batch (the roadmap items below, mostly shipped)

- **F1 HDACi pulsed self-maintaining tier** (`reversibility.py`, Nat Genet 2025): a curated
  pulsed-HDACi now earns a conditional epigenetic-memory WINDOW hypothesis instead of a flat
  abstain; default behaviour unchanged, never auto-durable.
- **F2 AMPA-PAM + PNN/ECM channels** (`scripts/114`): pre-registered negatives CONFIRMED -
  GRIA1 orthosteric 0.09, AMPA-PAM 0.26 (allosteric blindness), MMP9/PNN 0.42, all fail the
  size-matched gate. With TrkB-TMD (0.23), this COMPLETES the durability-node negative map:
  only senolytic-ablative (BCL2) is sequence-DTI-recoverable; TrkB-TMD, AMPA-PAM and PNN/ECM
  are all off-axis - which is exactly why durability is routed through the L4 permeability
  window and the engine stays abstain-by-default on the rest.
- **F3 pluggable efflux** (`free_exposure.py`): the Stage-3 efflux feature accepts an external
  ADMET-AI Pgp probability (guarded hook, Didziapetris fallback); both ADMET-AI-Pgp and a true
  unbound-Kp,uu spine are documented as data/dep-blocked upgrade paths, not faked.
- **F4 ledger expansion** (16->19): +tabernanthalog/MDMA/7,8-DHF (cited); recall 8/16 with the
  Jeffreys CI tightened to 0.27-0.73, neuroplasticity 6/6, off-channel mechanisms honestly 0.

## v2.9 dependency-resolution batch (D1-D4) - external deps resolved, claims re-tested

The standing "data/dep-blocked upgrade paths" from v2.8 were actually installed and exercised
(no new fabrication; every honest-negative and scope caveat survives). Env discipline held: the
heavy ML deps went into the system Python where torch is unused, and Boltz went into an isolated
3.12 venv, so MAMMAL's nightly cu128 / sm_120 torch was never touched.

- **D1 ADMET-AI Pgp Stage-3 retrain - SHIPPED** (`scripts/115` cache + `scripts/111` head-to-head,
  `free_exposure.py`, `cns_exposure.py`; report `cns_exposure_kpuu_v1.md`). Resolves the F3 efflux
  path end to end. ADMET-AI Pgp_Broccatelli (cached over all 1058 B3DB cpds, mean 0.337) vs the
  Didziapetris rule on ONE Bemis-Murcko scaffold split: **R2 0.276 vs 0.214 (+0.063)**, RMSE 0.59
  vs 0.61, conformal coverage 0.85 vs 0.84. Measurably better on logBB, but the gain sits on the
  PROXY and is dwarfed by the logBB-vs-Kp,uu residual the conformal band already carries, and
  Stage-3 only fires on confident P-gp-substrate downgrades. DECISION: the CI-safe, dependency-free
  rule model stays the DEFAULT; the ADMET variant ships as an OPT-IN (`PERSEUS_STAGE3_ADMET=1`,
  loaded only when admet_ai is importable; the model carries a `use_admet_ai` featurization-contract
  flag so train/inference features always match). Bonus rigor: the from-scratch numpy Mondrian
  split-conformal was cross-validated against `crepes` - max half-width disagreement **0.003 logBB**
  (0.000 on the per-category strata), validating the hand-rolled conformal math.
- **D2 Kp,uu data dependency - RESOLVED AS EXTERNALLY BLOCKED (documented, not faked).** A
  literature + web search found no clean, openly-licensed, machine-readable unbound-Kp,uu-with-SMILES
  set: Friden 2009 and Loryan/Morales 2024 live in reuse-restricted paper SI (PDF/XLSX needing manual
  extraction); the one larger collection (CMD-FGKpuu) is unlicensed and PARP-specific. So PERSEUS
  ships NO Kp,uu model - logBB stays the honest CC0 proxy and the abstain-wide band carries the
  residual. Swapping in a licensed Kp,uu spine (same featurizer + conformal pipeline; only the
  training table changes) remains the single highest-value future upgrade once reuse permission
  is in hand. Documented in `cns_exposure_kpuu_v1.md`.
- **D3 Boltz isolated install + honest TrkB-TMD test - RESOLVED, limitation stands.** **Boltz 2.2.1
  installed and import-verified** in `.venv-boltz312` (isolated 3.12 env). This makes the v2.7
  "even a Boltz second opinion is insufficient" claim empirical, not hypothetical: Boltz predicts an
  apo, single-copy, solvated fold and has no representation of the lipid bilayer, the cholesterol
  co-factor, or the crossed-dimer quaternary state that the TrkB-TMD site requires (Casarotto 2021;
  Cordeiro 2024) - so the off-axis verdict is robust to the "just run a structure model" rebuttal
  (informational gap, not a tooling gap). The legitimate future Boltz use is the **BCL2 senolytic
  channel** (soluble, well-defined groove - the one durability node where structural assumptions
  hold). Documented in `trkb_tmd_sitesplit_v1.md`.
- **D4 triptan-precision audit fix - SHIPPED** (`psychoplastogen.py`, commit 7a73750). Adversarial
  self-audit found a REAL L4 false positive: zolmitriptan (a 5-HT1B/1D triptan, NOT a 5-HT2A
  psychoplastogen) slipped the TPSA<=60 gate via its oxazolidinone. Added a triptan-pharmacophore
  veto (sulfonamide / cyclic carbamate) so triptans never reach the permeability gate; sumatriptan +
  zolmitriptan now window-negative, all psychedelics intact, decoys (dopamine/melatonin/
  methamphetamine/venlafaxine) still negative.

Dependency-resolution status (for reproduction): system Python 3.13 now has admet-ai 2.0.1 +
crepes 0.9.0 + mapie (torch untouched, used only for the ADMET variant + conformal cross-check);
`.venv-boltz312` has boltz 2.2.1. Full suite 619 passed / 2 skipped.

## Remaining roadmap (FUTURE_WORK)

1. **L1 Stage-3 free exposure (SHIPPED v2.6, efflux-aware conformal logBB).**
   `engine/free_exposure.py`, `scripts/110-111`, wired into `cns_exposure.py`. Replaces the
   Stage-3 ABSTAIN stub with a LightGBM logBB regressor on B3DB (CC0, 1058 cpds) + a cited
   Didziapetris P-gp efflux feature (the model-within-a-model lever) + Mondrian split-conformal
   bands (numpy) + a kNN applicability-domain abstain. Bemis-Murcko scaffold split (no analog
   leakage): R2 0.23, RMSE 0.60, conformal coverage 0.85 - honestly at/below the field ceiling
   (public CNS-penetration models plateau ~0.3-0.6; the deliverable is the calibrated interval
   + abstain rule, not a high R2). The Stage-3 gate refines a physchem PASS efflux-awarely:
   only a confidently predicted P-gp SUBSTRATE with sub-threshold logBB downgrades PASS->ABSTAIN
   (CNS-penetrant small molecules untouched, so the symptomatic + L4 plasticity heads are
   preserved - psychoplastogen recall held at 7/13). **HONEST SCOPE:** the trained target is
   logBB (total brain:plasma), a passive-penetration proxy; the unbound Kp,uu that finally
   governs free exposure is efflux-dominated and its public data are tiny/license-encumbered,
   so true Kp,uu remains the documented residual gap (swap in ADMET-AI Pgp + a licensed Friden/
   Morales Kp,uu spine to close it). Report: `reports/pipeline/cns_exposure_kpuu_v1.md`.
2. **Persistence-target DTI module (SHIPPED v1, calibration-gated; honest-negative).**
   `engine/persistence_dti.py`, `scripts/103-105`, `reports/pipeline/persistence_dti_*.md`.
   A 9-target substrate panel (BCL2/BCL-xL = ablative; HDAC1/2/6, DNMT1, EHMT2, KEAP1 =
   capability; NTRK2 = plasticity_window) is MAMMAL-scored so L3 can read substrate from
   PREDICTED target engagement for any chemical. Crucially it is CALIBRATION-GATED with two
   gates: a target may contribute only if (a) MAMMAL ranks its known engagers above matched
   non-engagers (AUROC + permutation-p) AND (b) the engagers clear a SPECIFICITY-first
   threshold (so the channel separates per-compound, not just in rank). **The decisive
   methodological move is SIZE-MATCHED negative controls.** MAMMAL's pKd is heavily
   molecular-weight-driven (corr(MW, pKd) = 0.61 over the 23 non-engagers; 0.89-0.91 within
   BCL2/BCL-xL/HDAC2), so a naive small-molecule-only negative pool spuriously passed 5/9
   targets (held-out 2/3 BH3-mimetics "recovered"). With a SIZE-MATCHED negative pool
   (129-671 Da) only **2/9 pass: BCL2 (AUROC 1.00) and BCL-xL (0.90)** - both the ablative
   (senolytic) tier; every capability/plasticity channel collapses toward chance
   (HDAC2 0.60, EHMT2 0.73, NTRK2 0.72, perm-p>0.05). Held-out demo: 0 non-persistence leaks
   and 0 flavonoid-senolytic false-durables (the channel correctly stays silent on
   fisetin/quercetin, which are caught instead by the L3 Tanimoto detector), but 0/3 held-out
   moderately-sized BH3-mimetics clear the size-matched threshold. **Honest conclusion:**
   MAMMAL gives a narrow, size-confounded population-RANKING signal for direct BCL2-family
   BH3-mimicry, NOT yet a general per-compound structure-computable persistence-substrate
   detector - so the persistence head stays abstain-by-default and the DTI substrate channel
   is wired but gated (only BCL2/BCL-xL usable, conservative threshold). The size-matched
   control is the deliverable: it prevents a molecular-weight artifact from masquerading as a
   persistence predictor. **MW-residualized re-calibration (`scripts/106`, CPU) answers "no
   signal, or signal masked by size?":** residualizing each score against a size-line fit on
   the non-engagers and re-calibrating RESCUES **NTRK2/TrkB** (raw AUROC 0.72 fail ->
   residualized 0.82 PASS, engagers fully in-MW-range so valid) - a genuine size-INDEPENDENT
   plasticity-substrate signal (MAMMAL recognises TRK inhibitors beyond their size); confirms
   the capability channels (HDAC/DNMT/EHMT/KEAP1, residualized AUROC <=0.62) as real
   size-artifacts. **BCL2 DE-ENTANGLED (v2.4, BH3-mimetic-sized negatives).** The first
   residualization run could not de-confound BCL2/BCL-xL because their BH3-mimetic engagers
   were larger than every non-engager (extrapolation). Adding 10 large non-BCL2 negatives
   into the engager size range (HIV/HCV antivirals, rifamycins, paclitaxel, digoxin; MW
   705-889, pool now 129-889 Da, confound r=0.73) raises BCL2's engager-in-MW-range to 0.67,
   and **BCL2 now SURVIVES residualization (raw 0.98 + residualized 0.83, both PASS)** - a
   genuine size-INDEPENDENT confirmation that MAMMAL recognises BH3-mimetic BCL2 engagement
   beyond molecular weight. BCL-xL stays size-entangled (raw 0.77 / residualized 0.67 fail,
   in-range 0.75); NTRK2 stays residualized-confirmed (0.81, in-range 1.00). **Net:** the
   durable axis now has ONE confirmed size-independent member (BCL2, ablative) plus a
   confirmed size-independent plasticity channel (NTRK2); BCL-xL and the reversible capability
   channels carry no size-free signal. Figure: `scripts/108` ->
   `reports/figures/persistence_dti/size_confound_deentanglement.png` (BCL2 engagers above the
   size line; raw-vs-residualized AUROC per channel). Next: gate operative inference on residualized-pass
   (keeps BCL2 + NTRK2, drops the size-entangled BCL-xL); the verified positive-compound
   ledger (psychoplastogens) to finally measure SENSITIVITY; a persistence-target head
   fine-tune.
3. **L4 psychoplastogen window (SHIPPED v2.5, intracellular-access leg).**
   `engine/psychoplastogen.py`, wired into `score()`; `scripts/107`. Encodes the Vargas 2023
   mechanism structurally: a serotonergic/monoaminergic-agonist SCAFFOLD (tryptamine with a
   small amine; psychedelic phenethylamine with >=2 aromatic OMe/halo; ergoline) AND
   intracellular ACCESS (TPSA<=60, HBD<=2 - the permeability determinants), gated on CNS PASS,
   emits a permissive WINDOW_CONDITIONAL (never auto-durable). Decisive validation: serotonin
   is window-NEGATIVE while DMT is window-POSITIVE despite near-identical 5-HT2A affinity - the
   discriminator is membrane permeability, OFF the DTI axis. On the verified positive ledger
   this lifts PERSEUS recall 0/13 -> 7/13 (neuroplasticity domain 5/5: psilocin/LSD/DMT/
   5-MeO-DMT/DOI) at FPR 0/15. The TrkB-transmembrane-domain leg (ketamine/fluoxetine site) is
   now CONFIRMED off-axis by a pre-registered site-split calibration (`scripts/112`, v2.7):
   NTRK2 tmd_wedge AUROC 0.23 (below chance), ecd 0.44, atp_pocket 0.59 - MAMMAL is blind to all
   three NTRK2 sites and worst at the durability TMD site, so durability is routed through the
   off-axis L4 window not a TrkB DTI score; even a Boltz-2 second opinion is insufficient
   (cholesterol-dependent crossed-dimer-in-lipid, Casarotto 2021).
4. **Persistence ground-truth ledger + empty-positive evaluator (SHIPPED v2.5).** Verified
   positive ledger recovered from the deep-research workflow (`data/raw/persistence_positive_
   ledger.csv`, 16 cited compounds). `validation/persistence_pu_eval.py` + `scripts/109`:
   sensitivity 0.54 (Jeffreys 95% CI 0.28-0.78) with FPR 0/15 (Jeffreys upper 0.15) and a
   PPV-vs-prior curve (0.03 at a 1% prior under the conservative upper FPR) - the bidirectional
   metric the empty-positive class previously made unidentifiable. SAR caveat surfaced (the
   positives are trial-availability-biased). v2.7 adds grouped leave-one-MECHANISM-out
   (`scripts/113`): serotonergic-psychedelic recall 0.86 (6/7), iboga 1/1, 0.00 on every
   non-serotonergic mechanism (correctly off-channel) - the window is an unfitted rule so LOMO
   confirms scaffold-generalization not chemotype memorization; plus a label-shift deployment
   transport (Saerens 2002 / Lipton 2018) making the rare-event PPV trap explicit (proper
   conformal needs a continuous score the categorical head lacks - it lives in Stage-3).
5. **Population x regime + bias-provenance covariates**: force ABSTAIN on healthy-young-
   durable-gain; report all-evidence vs bias-hardened scores.

## Appendix - claims the adversarial layer refuted (integrity record)

The verification agents caught and corrected 6 mis-cited / fabricated figures before they
entered the design: (1) ketamine single-dose relapse mis-attributed an 18-day median that
is actually the repeated-infusion figure; (2) "XSum is the best CMap connectivity method"
inverted (the cited paper finds Zhang best; method choice is contested); (3) an exenatide-PD
"cognitive 5.0 / p=0.006" that does not exist (the real endpoint was MOTOR -3.5; cognition
was non-significant) - the qualitative persistence-illusion narrative survives; (4) a
tolerance/beta-arrestin mechanism mis-cited to a paper that invokes lipophilicity instead;
(5) anti-amyloid "ADAS-Cog SMD 1.06" wrong by an order of magnitude (actual ~ -0.1 - the
biomarker/clinical disconnect is stronger with the correct tiny number); (6) a BDNF->TrkB
"bistable switch = formal definition of persistence" overstated (the cited model is self-
limiting, not bistable). All corrected; none load-bearing claims were left fabricated.
