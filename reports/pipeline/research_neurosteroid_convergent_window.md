# PERSEUS L4 expansion research: GABA-A neurosteroids + the convergent rapid-plasticity pathway

Research lane for the PERSEUS persistence engine. Question: can the L4 "psychoplastogen window"
(currently serotonergic-only, encoding Vargas 2023 membrane-permeability-to-intracellular-5-HT2A)
be extended to the non-serotonergic rapid-acting antidepressants (GABA-A neurosteroids, NMDA/
ketamine, muscarinic/scopolamine), and is there a single computable cross-class durability
signature or must each upstream class get its own structural window?

Author: deep-research synthesis, web-grounded, 2026-06-08. Every empirical claim carries an
author/year/venue citation. Honest negatives are reported as such. Spot-check any single row
before relying on it.

This lane sits ON TOP of, and is consistent with, the existing engine analysis in
`persistence_engineering_gaps.md` (which already established that the BDNF/mTORC1 cascade is a
DOWNSTREAM EFFECTOR, not a drug target, and that scoring affinity for MTOR/RPTOR is a category
error). It extends that finding to the two classes the gaps doc did not work through in
per-compound trial detail: the GABA-A neurosteroids and scopolamine.

---

## TL;DR (the crux answers)

1. **Neurosteroid durability is REAL but SHALLOW and CONTESTED, and it is NOT clean
   post-cessation disease-modification.** Brexanolone (IV, 60 h) shows within-patient symptom
   improvement sustained to the Day 30 follow-up (Meltzer-Brody 2018, Lancet). Zuranolone (oral,
   fixed 14-day course) shows high WITHIN-RESPONDER retention (86.1% of the Day-15 HAM-D
   improvement retained at Day 42 in WATERFALL; Clayton 2023, Am J Psychiatry) BUT the
   drug-versus-placebo SEPARATION decayed fast: in WATERFALL the between-group difference was only
   nominally significant through Day 12, and the companion MOUNTAIN trial MISSED its Day-15
   primary endpoint outright (Clayton 2023, J Clin Psychiatry). Clinically the drug is used as an
   EPISODIC, repeat-as-needed therapy, not a one-course cure (SHORELINE: median time to first
   repeat course 135-249 days; Cutler/Bonthapally 2023, J Clin Psychiatry). So the honest label is
   "rapid, with a multi-week tail and a high relapse/re-treatment rate," NOT "durable after
   cessation." This is closer to a long-offset symptomatic effect than to a delayed-start /
   randomized-discontinuation-grade durability claim.

2. **The convergent pathway is REAL, well-cited, and exactly as the prompt frames it
   (Zanos/Duman 2018; Kavalali & Monteggia 2012/2020): chemically diverse rapid-acting
   antidepressants converge downstream on a glutamate burst -> AMPA/GluA1 -> BDNF -> TrkB ->
   mTORC1 -> rapid synaptogenesis cascade.** BUT scopolamine and ketamine are explicitly placed on
   that cascade, whereas for the NEUROSTEROIDS the same authors write that "the intersection of
   these agents with the mechanisms underlying the rapid response to glutamatergic agents remains
   to be identified" (Duman, Shinohara, Fogaca, Hare 2019, Mol Psychiatry). Independent work even
   suggests the neurosteroid antidepressant effect is partly BDNF-INDEPENDENT (GABA-A PAM +
   membrane progesterone receptor + PXR; Frye 2014, Front Cell Neurosci; Zhang 2021, PMC7231971).
   So the convergence is neither universal nor settled even at the biology level.

3. **VERDICT: there is NO single structure-computable cross-class durability signature.** The
   convergence is purely DOWNSTREAM and INTRACELLULAR (it lives in BDNF/TrkB/mTORC1/GluA1, none of
   which is the drug's own target), so it is invisible from the chemical structure of the upstream
   ligand. The repo's own grouped leave-one-mechanism-out audit already proves this empirically:
   the serotonergic L4 window scores 0.00 recall on gaba_neurosteroid, muscarinic_antagonist, and
   nmda_dissociative (`perseus_lomo_transport_v1.md`, `scripts/113`). PERSEUS must therefore build
   PER-MECHANISM structural windows where one is justified, and ABSTAIN (route to the
   evidence/trial layer) where one is not. The cross-class object that DOES exist is not a
   structural feature; it is the L1 free-brain-exposure gate plus the evidence-design governor that
   are already shared.

4. **Per-class handle:** scopolamine -> a narrow, defensible structural window IS justifiable (a
   CNS-penetrant tertiary-amine muscarinic-antagonist tropane scaffold), but it is a
   mechanism-CLASS detector, not a durability detector, and it should fire only a
   plasticity-window hypothesis, never durable-by-construction. Neurosteroids -> NO honest
   structural durability window is justified; the pregnane/GABA-A-PAM scaffold is structurally
   recognizable but maps to a SHALLOW, contested, episodic-relapsing effect, so the correct output
   is ABSTAIN with an evidence-layer annotation. NMDA/ketamine -> already covered by the gaps doc
   as off-axis (TrkB-TMD); same abstain logic.

---

## PART A: neurosteroid post-cessation durability (honest, per-compound)

### A.0 Framing the question correctly

The PERSEUS persistence axis (`persistence_axis_v1.md`) judges DURABILITY-AFTER-CESSATION, scored
by trial DESIGN: delayed-start RCT (tier 7) > randomized discontinuation (6) > longitudinal
follow-up (5) > washout observation (4) > ... A "the drug works while you take it and the benefit
has a tail after you stop" result is a WASHOUT-OBSERVATION at best, and a fixed-course drug that is
then re-dosed on relapse is, by construction, telling you the effect did NOT durably hold. The
neurosteroids must be read against that bar, not against the (much lower) bar of "is the on-drug
effect real" (it is).

### A.1 Brexanolone / allopregnanolone (IV, ~60 h infusion; PPD)

- **Mechanism**: brexanolone is an IV formulation of the endogenous neurosteroid allopregnanolone,
  a positive allosteric modulator (PAM) of synaptic AND extrasynaptic GABA-A receptors
  (Meltzer-Brody 2018, Lancet; Frontiers review, Cornett 2021, Front Psychiatry 12:699740).
- **Pivotal evidence**: Meltzer-Brody S, Colquhoun H, Riesenberg R, et al. "Brexanolone injection
  in post-partum depression: two multicentre, double-blind, randomised, placebo-controlled, phase
  3 trials." Lancet 2018;392(10152):1058-1070 (PMID 30177236). Two trials (Study 1 enrolled
  HAM-D >=26; Study 2 HAM-D 20-25). Primary endpoint = change in HAM-D total at 60 hours.
  - Study 1: BRX60 -19.5 vs placebo -14.0 (difference -5.5, p=0.0013); BRX90 -17.7 vs -14.0
    (difference -3.7, p=0.0252).
  - Study 2: BRX90 -14.6 vs placebo -12.1 (difference -2.5, p=0.0160).
- **Durability after the infusion ends**: patients were followed to Day 30. The abstract states
  the effect had "rapid onset of action and durable treatment response during the study period"
  (Meltzer-Brody 2018). Pooled HUMMINGBIRD-program analyses report improvement maintained through
  Day 30 (Deligiannidis 2023, J Affect Disord, pooled analysis). A 2026 translational substudy
  found BDNF elevated at 6 h, 7 days, AND 30 days post-infusion, tracking sustained symptom
  improvement (Translational Psychiatry 2026, s41398-026-03834-9).
- **HONEST CAVEAT (load-bearing)**: I could NOT find a clean published statement that the
  brexanolone-versus-placebo HAM-D SEPARATION remained STATISTICALLY SIGNIFICANT at Day 30 in the
  pivotal trials (the trials were powered and reported at the 60-hour primary endpoint; the Day-30
  data are described as "sustained" within-group but the controlled between-group Day-30
  significance is not the headline result and I did not verify it as significant). So the strongest
  HONEST claim is: within-patient benefit is sustained to 30 days post-infusion (a
  WASHOUT-OBSERVATION tier signal), with BDNF as a candidate biomarker; whether the CONTROLLED
  difference persists to 30 days is not something I verified. Do not over-state this as
  delayed-start-grade durability.

### A.2 Zuranolone / SAGE-217 (oral, fixed 14-day course; MDD and PPD)

This is the central case the prompt flags, and the honest answer is the most nuanced.

- **Mechanism**: zuranolone is an oral neuroactive steroid GABA-A receptor PAM (synaptic +
  extrasynaptic), dosed as a FIXED 14-day once-daily course (Clayton 2023, Am J Psychiatry).
- **The durability question = does the effect persist after the 14-day course ends, or decay?**
  The pivotal data give a SPLIT answer depending on whether you measure within-responder retention
  or drug-versus-placebo separation:

  **WATERFALL (positive trial)** -- Clayton AH, Lasser R, Parikh SV, et al. "Zuranolone for the
  Treatment of Adults With Major Depressive Disorder: A Randomized, Placebo-Controlled Phase 3
  Trial." Am J Psychiatry 2023;180(9):676-684 (DOI 10.1176/appi.ajp.20220459). N=543 (534
  analyzed), zuranolone 50 mg vs placebo for 14 days, followed to Day 42.
  - Primary endpoint (Day 15 HAM-D change): -14.1 vs -12.3, statistically significant.
  - Onset by Day 3: -9.8 vs -6.8.
  - VERBATIM on durability: improvements "were sustained at all visits throughout the treatment
    and follow-up periods of the study (through day 42), with the difference remaining nominally
    significant through day 12." -> i.e. after dosing stopped, the placebo group kept improving and
    the BETWEEN-GROUP separation lost even nominal significance after Day 12.
  - Within-responder framing (the company-favoured number): Day-15 responders retained on average
    86.1% of their HAM-D improvement at Day 42, and 87.6% on MADRS (Biogen/Sage Phase 3 release,
    2021; reported in the WATERFALL/SKYLARK communications).

  **MOUNTAIN (negative trial)** -- Clayton AH, et al. "Zuranolone in Major Depressive Disorder:
  Results From MOUNTAIN -- A Phase 3, Multicenter, Double-Blind, Randomized, Placebo-Controlled
  Trial." J Clin Psychiatry 2023 (DOI 10.4088/JCP.22m14445). Zuranolone 20/30 mg vs placebo,
  14 days, observation to Day 42, extended follow-up to Day 182.
  - Day-15 PRIMARY ENDPOINT MISSED: HDRS-17 -12.5 (30 mg) vs -11.1 (placebo), p=0.116. Improvement
    vs placebo was significant only at Days 3, 8, 12 (i.e. ON-DRUG), not at the Day-15 primary or
    later.
  - Both MOUNTAIN and WATERFALL showed NO statistically significant zuranolone-vs-placebo
    difference at Day 42 (meta-analytic summary, Frontiers Pharmacol 2023, 10.3389/fphar.2023.1334694;
    Cao 2024 meta-analysis, PMC11079210).

  **SKYLARK (PPD, positive)** -- Phase 3 in postpartum depression, zuranolone 50 mg x14 days,
  followed to Day 45; met primary (Day-15 HAMD-17) and key secondaries, sustained through Day 45
  within-group (Sage/Biogen 2022 release; Deligiannidis 2023, Am J Psychiatry). PPD-specific.

  **CORAL** -- Parikh SV, et al. Neuropsychopharmacology 2024 (s41386-023-01751-9): zuranolone 50 mg
  CO-INITIATED with a standard antidepressant vs placebo+ADT for 14 days, then open-label ADT to
  Day 42. Tests rapid onset as an ADT adjunct, not single-agent post-cessation durability.

  **SHORELINE (open-label, the durability reality check)** -- Cutler AJ, et al. "Long-Term Safety
  and Efficacy of Initial and Repeat Treatment Courses With Zuranolone..." J Clin Psychiatry 2023
  (PMID 38153320). Patients who respond can receive REPEAT 14-day courses "as needed." Median time
  to first repeat course: 135 days (30 mg cohort) to 249 days (50 mg cohort); ~76.8% needed only
  1-2 total courses in up to a year. The US PPD label (approved 2023) is a 14-day course with
  repeat courses as needed.

- **HONEST VERDICT on zuranolone durability**: the effect is RAPID and has a multi-week tail
  WITHIN responders, but it is NOT post-cessation disease-modification. (a) The controlled
  drug-vs-placebo separation collapses within ~1-2 weeks of stopping (WATERFALL nominally
  significant only through Day 12; MOUNTAIN never separated at the Day-15 primary). (b) The drug is
  designed and labeled as an EPISODIC, repeat-as-needed therapy, which is the operational
  definition of an effect that does NOT durably hold (you re-dose because relapse is expected). On
  the PERSEUS evidence-design ladder this is a WASHOUT-OBSERVATION-tier within-group signal with a
  POSITIVE-then-decaying controlled contrast -- explicitly NOT delayed-start or
  randomized-discontinuation. The 86.1% retention figure is a within-responder conditional number
  and must not be reported as if it were a controlled durability effect size.

### A.3 Why the neurosteroid biology argues AGAINST a clean durability window

Even mechanistically, the neurosteroids do not sit cleanly on the durable-synaptogenesis cascade:
- Allopregnanolone's antidepressant action is attributed to GABA-A PAM activity PLUS membrane
  progesterone receptors PLUS the pregnane X receptor (PXR) in the VTA (Frye 2014, Front Cell
  Neurosci 8:106, PMC3988369), not primarily to a BDNF/TrkB/mTORC1 growth program.
- Hippocampal BDNF changes have been DISSOCIATED from the antidepressant-like behavioural effect
  for allopregnanolone (rapid BDNF regulation observed "even in the absence of an associated
  antidepressant-like effect"; Zhang 2021 review, PMC7231971). So the one node that would couple
  neurosteroids to the durable cascade is, for this class, not reliably load-bearing.
- This is consistent with the clinical picture (shallow, episodic, relapsing) rather than the
  psychedelic picture (single-dose spine increases persisting >=1 month; Shao 2021, Neuron).

---

## PART B: the convergent rapid-plasticity pathway and its structure-computability

### B.1 The convergence hypothesis is real and well-cited

The hypothesis the prompt names (Duman; Kavalali & Monteggia) is solidly established in review:

- **Zanos P, Thompson SM, Duman RS, Zarate CA Jr, Gould TD. "Convergent Mechanisms Underlying
  Rapid Antidepressant Action." CNS Drugs 2018 Mar;32(3):197-227** (PMID 29516301,
  DOI 10.1007/s40263-018-0492-x). THE canonical convergence reference. Names ketamine,
  scopolamine, GLYX-13/rapastinel, GluN2B-NAMs, (2R,6R)-HNK, NMDA glycine-site modulators, mGluR2/3
  antagonists, and GABA-A modulators, and states they converge on "potentiation of excitatory
  synapses" via mTOR activation, enhanced BDNF/TrkB signalling, AMPAR activation, enhanced protein
  synthesis, and strengthened excitatory synapses.

- **Duman RS, Shinohara R, Fogaca MV, Hare B. "Neurobiology of rapid-acting antidepressants:
  convergent effects on GluA1-synaptic function." Mol Psychiatry 2019;24:1816-1832**
  (PMC6754322; already cited in `persistence_engineering_gaps.md`). VERBATIM: "although rapid
  antidepressant actions are produced by multiple classes of agents... there are some convergent,
  downstream mechanisms," including "increased BDNF release and/or expression, activation of
  protein synthesis pathways (i.e., mTORC1 and eEF2 kinase), increased expression of synaptic
  proteins (GluA1, PSD95, and synapsin), and increased synaptic number and function in the mPFC."

- **Kavalali ET, Monteggia LM. "Synaptic Mechanisms Underlying Rapid Antidepressant Action of
  Ketamine." Am J Psychiatry 2012;169:1150-1156** (PMID 23534055) and Kavalali & Monteggia 2020
  (Neuron, "Targeting Homeostatic Synaptic Plasticity for Treatment of Mood Disorders"): the
  synaptic-plasticity / BDNF-translation framework underpinning the convergence.

- **Anacker C. "New Insight Into the Mechanisms of Fast-Acting Antidepressants: What We Learn From
  Scopolamine." Biol Psychiatry 2018;83(1):e5-e7** (PMID 29173709): states scopolamine and ketamine
  CONVERGE -- both reduce Ca2+ in GABAergic interneurons -> decreased GABA -> increased glutamate ->
  BDNF/TrkB -> ERK1/2 + Akt -> mTOR -> synaptogenesis.

So the cascade glutamate burst -> AMPA/GluA1 -> BDNF -> TrkB -> mTORC1 -> rapid spine formation is
the well-supported common downstream effector for ketamine and scopolamine, and is the proposed
(though, per Duman 2019, not-yet-confirmed) integration point for the neurosteroids.

### B.2 The convergence is DOWNSTREAM and INTRACELLULAR, hence NOT structure-computable from the ligand

This is the crux. The shared cascade is, node by node, a property of the NEURON, not of the drug:
- The upstream RECEPTORS are chemically unrelated and mutually exclusive: NMDA channel pore
  (ketamine), M1/M-family muscarinic orthosteric site (scopolamine), GABA-A neurosteroid allosteric
  site (allopregnanolone/zuranolone), intracellular 5-HT2A (psychedelics). A structure encodes
  WHICH receptor it can hit, not the downstream program that receptor triggers.
- The convergent nodes -- BDNF, TrkB, mTORC1, GluA1, PSD95, synapsin -- are DOWNSTREAM EFFECTORS,
  i.e. NOT drug targets. The repo's `persistence_engineering_gaps.md` already states this precisely:
  scoring a compound's affinity for MTOR/RPTOR/RICTOR is a "category error"; the cascade "belongs in
  the transcriptomic/pathway axis (Cluster D / LINCS), not the DTI panel." That finding generalizes
  one-for-one to the cross-class case: there is no ligand-structure feature that reads out "this
  molecule will, via whatever receptor it hits, drive a BDNF->mTORC1->GluA1 growth program."
- Each upstream class's durability is gated by a DIFFERENT, class-specific physical bottleneck that
  is NOT shared:
  - psychedelics: membrane permeability to the intracellular 5-HT2A pool (Vargas 2023, Science
    379:700-706) -- an ADMET/physchem property, the ONE that is cleanly structure-computable and is
    exactly what the existing L4 window exploits.
  - ketamine/fluoxetine (TrkB-TMD leg): a cholesterol-dependent crossed-dimer transmembrane wedge
    (Casarotto 2021, Cell 184:1299; Cordeiro 2024, Nat Commun 15) -- requires the lipid bilayer +
    dimer + cholesterol; information-theoretically absent from a 1D sequence and not a clean
    small-molecule structural alert either.
  - scopolamine: M-receptor antagonism + CNS penetration -- structure-recognizable as a CLASS, but
    the durability is downstream, not in the ligand.
  - neurosteroids: GABA-A PAM + mPR + PXR -- structure-recognizable as a pregnane CLASS, but maps to
    a shallow/episodic effect, so the scaffold does NOT predict durability.
- Therefore the SAME downstream cascade is reached by FOUR different, non-interchangeable structural
  routes, each with its own gate. A single feature that fires on all four would have to encode the
  downstream convergence point, which is not in the molecule.

### B.3 The repo already proves this empirically (do not re-derive)

`perseus_lomo_transport_v1.md` (`scripts/113`, grouped leave-one-mechanism-out): the serotonergic
L4 window scores

| mechanism class | recall |
|---|---|
| serotonergic_psychedelic | 0.88 (7/8) |
| iboga_atypical | 1.00 (1/1) |
| **gaba_neurosteroid** | **0.00 (0/1)** |
| **muscarinic_antagonist** | **0.00 (0/1)** |
| **nmda_dissociative** | **0.00 (0/2)** |
| **trkb_agonist** | **0.00 (0/1)** |

The doc already (correctly) reads these zeros as "off-channel, not false negatives of this rule."
That IS the empirical demonstration that the convergence does not transfer a serotonergic
structural feature to the other classes. The persistence-DTI side concurs: only BCL2 (ablative) is
size-independently recoverable and NTRK2/TrkB is rescued only after MW-residualization, while the
TrkB-TMD durability site is below chance (AUROC 0.23; `trkb_tmd_sitesplit_v1.md`, `scripts/112`).

---

## DELIVERABLES

### 1. Neurosteroid durability evidence (honest, per-compound, with how-long-after-dosing)

| compound | route / course | best durability evidence | how long after dosing | honest design tier | verdict |
|---|---|---|---|---|---|
| Brexanolone (allopregnanolone) | IV, ~60 h | Meltzer-Brody 2018 Lancet 392:1058 (PMID 30177236); within-group HAM-D sustained to Day 30; BDNF elevated to Day 30 (Transl Psychiatry 2026) | 30 days post-infusion (within-group) | washout_observation (controlled Day-30 separation NOT verified by me) | shallow post-cessation tail; NOT disease-modification |
| Zuranolone / SAGE-217 (MDD) | oral, 14-day course | WATERFALL: Clayton 2023 AJP 180:676; Day-15 primary met; within-responder 86.1% HAM-D retained at Day 42 | controlled separation only nominally sig through Day 12; within-responder to Day 42 | washout_observation, POSITIVE-then-decaying | rapid + multi-week tail; effect does NOT durably hold (episodic re-dosing) |
| Zuranolone (MDD, second pivotal) | oral, 14-day | MOUNTAIN: Clayton 2023 JCP; Day-15 primary MISSED (p=0.116); no Day-42 separation | on-drug only (Days 3-12) | parallel_rct, negative at primary | confirms shallow/inconsistent durability |
| Zuranolone (PPD) | oral, 14-day | SKYLARK: Deligiannidis 2023 AJP; met Day-15 primary, within-group sustained to Day 45 | 45 days within-group (PPD) | washout_observation | PPD-specific; same episodic caveat |
| Zuranolone (real-world durability) | oral, repeat as needed | SHORELINE: Cutler 2023 JCP (PMID 38153320); median 135-249 days to first repeat course | n/a (relapse/re-treat design) | n/a | operational proof the effect is episodic, not durable |

### 2. Convergent-pathway evidence and key references

- Zanos, Thompson, Duman, Zarate, Gould 2018, CNS Drugs 32(3):197-227 (PMID 29516301) -- canonical
  convergence statement; names ketamine, scopolamine, GABA-A modulators, rapastinel.
- Duman, Shinohara, Fogaca, Hare 2019, Mol Psychiatry 24:1816-1832 (PMC6754322) -- convergent
  GluA1/mTORC1/BDNF; explicitly says neurosteroid intersection with the glutamatergic cascade
  "remains to be identified."
- Kavalali & Monteggia 2012, Am J Psychiatry 169:1150-1156 (PMID 23534055); Kavalali & Monteggia
  2020, Neuron -- synaptic-plasticity / BDNF-translation framework.
- Anacker 2018, Biol Psychiatry 83(1):e5-e7 (PMID 29173709) -- scopolamine-ketamine convergence,
  spelled out node by node.
- Casarotto 2021, Cell 184:1299; Cordeiro 2024, Nat Commun 15 -- TrkB-TMD leg (off the sequence
  axis).
- Vargas 2023, Science 379:700-706 -- intracellular 5-HT2A / membrane permeability (the serotonergic
  leg PERSEUS already exploits).
- Shao 2021, Neuron (psilocybin spines >=1 month) -- the durable structural read-out.

### 3. VERDICT: single computable cross-class signature, or per-mechanism windows?

**No single structure-computable cross-class durability signature exists. PERSEUS must build
per-mechanism structural windows (one per upstream class) where one is justified, and ABSTAIN
elsewhere.**

Justification:
- The convergence is purely DOWNSTREAM and INTRACELLULAR (BDNF/TrkB/mTORC1/GluA1). Those nodes are
  effectors, not drug targets; nothing in the ligand's structure reads them out. (Already the
  repo's stated position for the within-class case; it generalizes verbatim to cross-class.)
- The four upstream classes hit chemically unrelated, mutually exclusive receptors, and each
  class's durability is gated by a DIFFERENT physical bottleneck (membrane permeability for
  psychedelics; cholesterol-dependent TMD wedge for the TrkB leg; M-antagonism+CNS for scopolamine;
  GABA-A PAM pregnane for neurosteroids). There is no common molecular feature to compute.
- Empirically the repo's own grouped LOMO already shows 0.00 transfer of the serotonergic window to
  gaba_neurosteroid / muscarinic / nmda classes, and the persistence-DTI calibration shows the
  durability nodes are off MAMMAL's axis (TrkB-TMD AUROC 0.23).
- The genuine CROSS-CLASS objects in PERSEUS are NOT structural features: they are (a) the L1
  free-brain-exposure gate (every rapid-acting antidepressant must reach the brain) and (b) the
  evidence-design governor (only delayed-start / randomized-discontinuation can raise persistence
  confidence). Those already span classes and should remain the shared backbone. The structural
  layer must stay class-typed.

This is the consistent, abstain-by-default answer and it preserves the engine's central guardrail:
a window OPENS plasticity; durability is conditional and must be demonstrated by trial design, not
inferred from a downstream cascade the molecule does not encode.

### 4. Scopolamine and neurosteroids: structure-computable handle, or abstain?

**Scopolamine (muscarinic): a NARROW structural CLASS window is defensible, but it is NOT a
durability detector.**
- Durability evidence is real and off-drug: Furey & Drevets 2006 (Arch Gen Psychiatry 63:1121,
  PMID 17015814) crossover -- subjects who received scopolamine FIRST retained a 38% MADRS
  reduction from baseline THROUGH the subsequent placebo block (a genuine carryover-after-cessation
  signal); replicated Drevets & Furey 2010 (Biol Psychiatry 67:432, PMID 20074703). Note the
  later parallel-group dose-response study (Jaffe-style, Trials 2020, PMC7011244) is more equivocal,
  so the durability is real but not uncontested.
- Structure-computable handle: scopolamine is a CNS-penetrant TERTIARY-amine TROPANE muscarinic
  ANTAGONIST. A tropane-scaffold + tertiary-amine + CNS-PASS detector is implementable in the same
  RDKit/SMARTS style as `psychoplastogen.py`. The decisive precision veto mirrors the existing
  module's quaternary-charge logic: the QUATERNARY tropane congeners (e.g. methscopolamine,
  ipratropium, tiotropium) are BBB-impermeant and antidepressant-inactive, and the L1 gate's
  permanent-charge veto already kills them -- so the window would correctly separate scopolamine
  (tertiary, permeant) from methscopolamine (quaternary, impermeant), exactly the serotonin-vs-DMT
  style discriminator the project favours.
- BUT: this window detects the MECHANISM CLASS, not durability. The durability is downstream
  (Anacker 2018). So it must fire only a `plasticity_window` hypothesis (permissive, never
  auto-durable) AND it would be precision-fragile (most tertiary-amine muscarinic antagonists --
  atropine, benztropine, oxybutynin, many tricyclics -- are NOT rapid antidepressants; antagonism
  at M-subtypes and regional selectivity are not structure-derivable). Recommendation: implement it
  ONLY as a low-precision class flag that ROUTES TO THE EVIDENCE LAYER, or abstain. Do not present
  it as a durability predictor.

**Neurosteroids (GABA-A): ABSTAIN. No honest structural durability window is justified.**
- The pregnane / neuroactive-steroid scaffold (allopregnanolone, zuranolone, ganaxolone) IS
  structurally recognizable (a 3-alpha-hydroxy-5-alpha-pregnan-20-one-like steroid, a GABA-A PAM
  pharmacophore). So a scaffold detector is buildable.
- But it would map to an effect that is SHALLOW, CONTESTED, and EPISODIC-RELAPSING (Part A), i.e.
  the scaffold does NOT predict post-cessation durability -- it predicts a rapid effect that decays.
  Firing a durability/plasticity-window flag on a pregnane scaffold would be the precise kind of
  false positive PERSEUS exists to prevent. And the biology is BDNF-dissociable for this class
  (Zhang 2021; Frye 2014), so it does not even cleanly join the convergent cascade.
- Correct output: ABSTAIN at the structural layer; annotate at the EVIDENCE layer as
  `washout_observation` (positive-then-decaying for zuranolone; within-group-to-Day-30 for
  brexanolone), explicitly NOT delayed-start/discontinuation. This is exactly the
  `class_extrapolation` / evidence-design machinery the engine already has.

**Ketamine / NMDA: ABSTAIN at structure (per existing gaps doc).** TrkB-TMD durability site is
off-axis (Casarotto 2021; AUROC 0.23 in `scripts/112`). Annotate at the evidence layer
(R-ketamine, nitrous oxide, lanicemine are in the ledger with their real, mostly washout-tier
designs).

### 5. Concrete proposed design + pre-registered falsifier

**Design (consistent with the abstain-by-default, per-mechanism, evidence-governed architecture):**

A. **Do NOT build a unified cross-class structural durability head.** Record this as a documented
   NEGATIVE (the convergence is downstream/intracellular, not ligand-computable), citing Zanos 2018,
   Duman 2019 (the "remains to be identified" line), and the in-repo LOMO 0.00 transfer. This is a
   publishable methods point that extends the existing allosteric-blindness / off-axis story.

B. **Add a small set of MECHANISM-CLASS structural ROUTERS (not durability detectors),** each
   emitting only a class tag + `route_to_evidence_layer`, reusing the `psychoplastogen.py` pattern
   (RDKit SMARTS + CNS gate + precision vetoes):
   - `muscarinic_tropane` flag: tropane + tertiary amine + CNS-PASS; quaternary-tropane veto
     (already covered by L1 permanent-charge veto). Output: `plasticity_window` HYPOTHESIS, low
     confidence, with the caveat that durability is downstream and unproven from structure.
   - `neuroactive_steroid` flag: 3-alpha-hydroxy-pregnan-20-one-like GABA-A-PAM scaffold. Output:
     ABSTAIN on durability + evidence-layer annotation `washout_observation, decaying`. NOT a
     window-positive.
   These extend the persistence ground-truth ledger schema (`data/raw/persistence_ground_truth.csv`:
   compound, mechanism_class, persistence_design, persistence_label, basis, source,
   structure_scoreable) with the rows above, and slot mechanism_class into the existing grouped-LOMO
   loop so recall stays mechanism-resolved.

C. **Keep the cross-class burden on the SHARED non-structural layers** (L1 free-exposure gate + the
   evidence-design governor), which already span all classes correctly.

D. **Ledger hygiene:** the existing positive ledger (`persistence_positive_ledger.csv`) lists
   zuranolone and scopolamine as "verified." Re-tag both honestly: scopolamine =
   washout_observation (Furey 2006 carryover, real but with a more equivocal later parallel-group
   study); zuranolone = washout_observation, POSITIVE-then-DECAYING (controlled separation lost
   after ~Day 12; episodic re-dosing). Neither is delayed-start/discontinuation. This prevents the
   ledger from inflating the persistence positive class with shallow/episodic effects, which would
   poison the PU sensitivity estimate.

**Pre-registered falsifiers (state the expected result BEFORE running, per perseus_design.md):**

- FALSIFIER 1 (structural-window justification): build the `neuroactive_steroid` scaffold detector
  and the `muscarinic_tropane` detector, then run them against the persistence GROUND-TRUTH ledger.
  PRE-REGISTERED EXPECTATION: neither flag correlates with the `persistence_label`
  (delayed-start/discontinuation positive) better than chance -- i.e. the scaffolds detect the
  mechanism class but carry NO durability information. If a scaffold DID predict
  delayed-start-grade durability above chance, that would FALSIFY the "no structural durability
  window" verdict for that class and justify promoting it from a class router to a real window.

- FALSIFIER 2 (cross-class non-transfer, confirmatory): extend the grouped-LOMO
  (`scripts/113`) so that ANY single structural feature trained/defined on one mechanism class is
  applied to held-out classes. PRE-REGISTERED EXPECTATION: cross-class recall stays at/near 0.00
  for gaba_neurosteroid, muscarinic, nmda (as the serotonergic window already does). A
  cross-class recall meaningfully above 0 from a single structural feature would FALSIFY the
  "downstream-only convergence, not structure-computable" verdict.

- FALSIFIER 3 (neurosteroid durability honesty): if a future randomized-discontinuation or
  delayed-start trial of zuranolone/brexanolone shows the CONTROLLED benefit persists after
  cessation (not just within-group), that would upgrade the neurosteroid row from
  washout_observation to a higher tier and would justify revisiting the abstain decision. Current
  evidence (controlled separation lost by ~Day 12; episodic re-dosing) predicts this will NOT
  happen for the single-course design.

---

## Honesty ledger for this lane

- Verified primary citations (author/year/venue/PMID where available): Meltzer-Brody 2018 Lancet
  (PMID 30177236); Clayton 2023 AJP WATERFALL (DOI 10.1176/appi.ajp.20220459); Clayton 2023 JCP
  MOUNTAIN (DOI 10.4088/JCP.22m14445); Cutler 2023 JCP SHORELINE (PMID 38153320); Furey & Drevets
  2006 Arch Gen Psychiatry (PMID 17015814); Zanos 2018 CNS Drugs (PMID 29516301); Duman 2019 Mol
  Psychiatry (PMC6754322); Kavalali & Monteggia 2012 AJP (PMID 23534055); Anacker 2018 Biol
  Psychiatry (PMID 29173709); Vargas 2023 Science 379:700; Casarotto 2021 Cell 184:1299; Shao 2021
  Neuron; Frye 2014 Front Cell Neurosci (PMC3988369).
- FLAGGED UNCERTAINTIES (not fabricated, explicitly unresolved): (a) I did NOT verify that the
  brexanolone-vs-placebo HAM-D difference was STATISTICALLY SIGNIFICANT at Day 30 in the pivotal
  trials -- only that within-group benefit was "sustained" to Day 30; the controlled Day-30
  contrast is unverified. (b) Exact page numbers / volume for MOUNTAIN (JCP) and the SKYLARK
  primary paper were taken from secondary summaries; confirm against the primary PDFs before
  citing in a manuscript. (c) The "86.1% / 87.6% retention" figures are from Sage/Biogen Phase 3
  communications (within-responder conditional numbers), not a peer-reviewed controlled effect
  size; treat accordingly. (d) Deligiannidis 2023 AJP (SKYLARK) and Deligiannidis 2023 J Affect
  Disord (brexanolone pooled) author/venue attributions should be re-checked against the primary
  records.
- Nothing in this lane should be read as a positive durable-cognition repurposing hit. The
  high-value output is a rigorous NEGATIVE (no cross-class structural signature) plus an honest
  re-tagging of two ledger rows that were over-credited as durable.
