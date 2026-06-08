# PERSEUS - Persistence-aware Pro-cognition Engine (design + v1 status)

Source: an adversarially-verified Opus multi-agent research synthesis (2026-06-07; 12 web
research lanes each independently verified, + a repo data audit + chief-architect
synthesis; 26 agents). The verification layer refuted 6 fabricated/mis-cited numbers
before they entered the design (appendix). This document is the engineering spec and the
status of the buildable v1.

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
- **F2 shortlist re-scored**: the uniform +0.40 symptomatic prior splits into 5 CNS-
  excluded, 11 not-cognition-excluded, 12 null/symptomatic, 2 CONTESTED (selegiline,
  rasagiline; delayed-start), 1 WINDOW_CONDITIONAL (fluoxetine). The persistence head is a
  model output, not a buried caveat.

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

## Remaining roadmap (FUTURE_WORK)

1. **L1 Stage-3 free exposure (Kp,uu)**: fit an efflux-aware unbound brain/plasma regressor
   (small public rat Kp,uu / B3DB set) + conformal applicability band; until then Stage 3
   returns ABSTAIN.
2. **Persistence-target DTI module**: add HDAC1/2/3/6, DNMT1/3A/3B, EHMT1/2, KEAP1, MTOR,
   eIF2B, BCL2/BCL-xL, SRC, NTRK2 to the panel + MAMMAL-score them, so L3 can read substrate
   from predicted target engagement (not only class + structural alert) for any chemical.
3. **L4 full**: TrkB transmembrane-domain binding head (psychoplastogen anchors) +
   intracellular-5-HT2A access (passive permeability x 5-HT2A affinity).
4. **Persistence ground-truth ledger + PU/LOMO evaluator**: curated delayed-start /
   discontinuation / washout cognition outcomes + the empty-positive evaluator and the
   "persistence-illusion" negative-control suite.
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
