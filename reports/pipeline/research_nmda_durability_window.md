# PERSEUS L4b research: an off-axis, structure-computable "NMDA-plasticity window"?

Deep-research pass for a SECOND psychoplastogen-style persistence window, this time for the
NMDA-antagonist / dissociative rapid-acting antidepressant (RAAD) class, analogous to the
existing serotonergic L4 window (`src/mammal_repurposing/engine/psychoplastogen.py`). The L4
window works because it pairs a SCAFFOLD class with an OFF-AXIS, structure-computable
discriminator (membrane permeability via TPSA/HBD, encoding Vargas 2023), and it is validated
by a PRE-REGISTERED NEGATIVE (serotonin must score window-NEGATIVE despite near-isoaffine
5-HT2A binding). The question here: is there an analogous off-axis, structure-computable
property for the NMDA class, with memantine as the pre-registered negative?

Bottom line up front (the honest headline): the DURABILITY mechanism is well established, but
the decisive memantine-vs-ketamine discriminator is a PHARMACODYNAMIC channel-block property
(use-dependence / resting-state block / trapping kinetics), and on the cheap 2D structural
axes that powered the L4 window, memantine and ketamine are nearly IDENTICAL. So the L4
"permeability trick" does NOT transfer. An L4b window can still be built, but it is honestly
WEAKER than L4: most of its discriminating power must come from curated pharmacology
(measured off-rate, trapping fraction, resting/use-dependence), not from RDKit descriptors.
The defensible design is therefore a scaffold-flag plus a small curated pharmacodynamic table,
with ABSTAIN as the default for any NMDA antagonist whose trapping/use-dependence is not
measured. This is a narrower, more honest claim than L4, and it is itself the valuable result.

---

## 1. MECHANISM OF DURABILITY (acute effect vs post-cessation persistence)

The canonical molecular cascade for a single/brief NMDA-antagonist exposure producing
plasticity that OUTLASTS drug clearance:

1. ACUTE (drug present): an uncompetitive (open-channel) NMDA-antagonist preferentially
   silences NMDARs on fast-spiking GABAergic interneurons, DISINHIBITING pyramidal neurons
   and producing a transient prefrontal GLUTAMATE surge (the "glutamate burst"). The surge
   drives AMPA-receptor activation.
   - Disinhibition / glutamate-surge model: Moghaddam et al. 1997 (J Neurosci 17:2921);
     reviewed in Duman & Aghajanian 2012 (Science 338:68) and Krystal/Abdallah/Sanacora/Duman
     2019 (Neuron 101:774).
2. ACUTE -> CONVERSION: AMPA activation, together with a halt of tonic NMDAR-driven
   eEF2-kinase (CaMKIII) activity, DE-SUPPRESSES BDNF translation; BDNF released onto TrkB
   activates PI3K-Akt and ERK, converging on mTORC1.
   - eEF2K / resting-NMDAR-block -> BDNF de-suppression: Autry, Adachi, Monteggia et al. 2011
     (Nature 475:91); Nosyreva et al. 2013 (J Neurosci 33:6990).
   - BDNF-TrkB requirement: Autry et al. 2011 (Nature 475:91); Lepack et al. 2014/2015
     (Int J Neuropsychopharmacol).
   - mTORC1 -> synaptogenesis: Li, Lee, Duman et al. 2010 (Science 329:959); the
     antidepressant behavioral effect is abolished by rapamycin.
3. POST-CESSATION PERSISTENCE (drug cleared): mTORC1-driven synthesis of synaptic proteins
   (GluA1, PSD-95, synapsin) restores/creates dendritic spines in mPFC. New spines and the
   restored excitatory tone are the structural substrate that PERSISTS after the ~2-3 h
   ketamine half-life. Ketamine's behavioral effect lasts ~1 week in rodents and ~7 days in
   humans from a single dose, despite the drug being long gone, because the synaptic
   remodeling is self-sustaining over that window.
   - Spine restoration outlasting drug: Li et al. 2010 (Science 329:959); Moda-Sava et al.
     2019 (Science 364:eaat8078) showed spine FORMATION is needed to SUSTAIN (vs initiate)
     the behavioral effect, dissociating induction from maintenance.

Key conceptual point the project already encodes for the serotonergic class (Ly 2021; the
"stimulation phase" then "growth phase"): the durable change requires a TRANSIENT trigger that
hands off to a self-sustaining downstream program (BDNF-TrkB-mTORC1-spine). The NMDA antagonist
is the TRIGGER; mTORC1/BDNF/spines are the persistence substrate. This is exactly the
permissive-window logic (`reversibility.py` `plasticity_window` tier): the window OPENS durable
plasticity but is not itself the durable change, so it must never auto-promote to
"demonstrated durable."

ACUTE vs PERSISTENCE summary: the acute effect (disinhibition, glutamate surge, dissociation)
is purely a drug-present pharmacodynamic event and reverses on clearance. The PERSISTENCE is
the BDNF-TrkB-mTORC1-spine program it triggers. Antagonists that hit the channel but do NOT
trigger that program (memantine, see section 3) produce the acute block WITHOUT the persistence.

---

## 2. HONESTY GATE: is post-cessation durability ESTABLISHED per compound?

Honest verdicts (the project values negatives; nothing below is upgraded beyond the evidence):

- (R,S)-ketamine / esketamine: ESTABLISHED for the SINGLE-DOSE rapid effect lasting roughly
  up to 1 week, then waning. The single-dose effect is real and replicated (Berman et al.
  2000 Biol Psychiatry 47:351; Zarate et al. 2006 Arch Gen Psychiatry 63:856). BUT true
  open-ended POST-CESSATION durability is OVERSTATED in casual framing: relapse is the norm
  without repeated dosing. Esketamine maintenance (SUSTAIN-1, Daly et al. 2019 JAMA
  Psychiatry 76:893) shows that CONTINUED intermittent dosing delays relapse versus
  discontinuation, i.e. the benefit is maintained by ongoing dosing, not by a one-time durable
  reset. So: durable for ~days from one dose; NOT a permanent single-exposure cure. The single
  dose opens a window; sustained remission needs repetition (or paired psychotherapy, still
  under study).

- (2R,6R)-hydroxynorketamine (HNK): CONTESTED / OVERSTATED. Zanos, Gould et al. 2016 (Nature
  533:481) reported HNK reproduces ketamine's antidepressant-like effects in mice and claimed
  these are NMDAR-INHIBITION-INDEPENDENT. This is genuinely disputed on two fronts:
  (i) MECHANISM: Suzuki, Nosyreva, Kavalali, Monteggia 2017 (Nature 546:E1, reply to Zanos)
      argued HNK does block synaptic NMDARs much like ketamine and that its signaling effects
      couple to NMDAR inhibition. Counter-evidence that antidepressant-relevant (low)
      concentrations do NOT block NMDARs: Lumsden et al. 2019 (PNAS 116:5160). An alternative
      mechanism (mGlu2-receptor-dependent, presynaptic glutamate release) was proposed by Zanos
      et al. 2019 (PNAS 116:6441) and Riggs et al. 2019/2020. Net: the receptor target of HNK
      is NOT settled.
  (ii) REPLICATION: the antidepressant-like effect of HNK has been inconsistent across labs.
      Notably, the same group's broader R-ketamine work and independent groups diverge; Yang
      et al. reported failures to reproduce HNK effects in some rodent models (LPS, chronic
      social defeat). I did not independently verify each replication's exact design, so I flag
      this as "contested replication" rather than asserting a specific failed-trial count.
  (iii) CLINICAL DIRECTION: in humans, HIGHER plasma HNK after ketamine infusion has been
      associated with WORSE, not better, antidepressant response (Grunebaum / Milak-adjacent
      analyses; summarized in Farmer et al. 2020 Neuropsychopharmacol 45:1398, "is this the
      end of the HNK pipeline?"). This is the opposite of what the HNK-is-the-active-agent
      hypothesis predicts. HNK durability in HUMANS is essentially unestablished.
      VERDICT: HNK is a real preclinical signal with a CONTESTED mechanism and weak/negative
      human translation. Do NOT treat HNK durability as established.

- Nitrous oxide (N2O): WEAK durability evidence; partly OVERSTATED. Nagele et al. 2015 (Biol
  Psychiatry 78:10, proof-of-concept) and Nagele et al. 2021 (Sci Transl Med 13:eabe1376,
  phase 2, 25% vs 50% N2O) show a rapid antidepressant effect. On DURABILITY the honest read
  is that effects were observed starting within ~2 h and lasting up to ~1 week, and the
  phase-2 design could NOT determine whether efficacy persisted beyond 1 week after a single
  session (the authors say so explicitly). A widely-cited "2-month remission after a single
  inhalation" is a SINGLE CASE REPORT (Nagele et al. 2020, Front Psychiatry 11:692), explicitly
  anecdotal, not controlled. VERDICT: rapid effect plausible; multi-week single-dose durability
  is NOT established for N2O.

- Dextromethorphan (DXM) and DXM-bupropion (AXS-05, Auvelity): efficacy ESTABLISHED for the
  COMBINATION as a chronically dosed oral antidepressant, but this is NOT a clean
  single-exposure-durability story. AXS-05 met endpoints in GEMINI (Iosifescu et al. 2022,
  J Clin Psychiatry 83:21m14345) and an earlier RCT (ASCEND); remission ~39.5% vs ~17.3%
  placebo with onset within ~1 week. CRITICAL CONFOUND for a "durability window": bupropion is
  present specifically to inhibit CYP2D6 and RAISE DXM exposure, and DXM is also a sigma-1
  agonist and a serotonin/NE transporter inhibitor. So AXS-05 is a CONTINUOUSLY DOSED,
  MULTI-TARGET drug; its benefit is maintained by ongoing dosing, not demonstrated as
  post-cessation durability. As an NMDA-window exemplar DXM is MUDDY (multi-mechanism, chronic
  dosing). Treat AXS-05 efficacy as established for the drug-on state; post-cessation
  durability is NOT demonstrated and the NMDA contribution is not cleanly isolable.

- Memantine (THE NEGATIVE CONTROL): ESTABLISHED NEGATIVE for rapid antidepressant / rapid
  plastogenic action. Memantine is an uncompetitive NMDA open-channel blocker (approved for
  Alzheimer's), yet controlled trials do NOT show ketamine-like rapid antidepressant efficacy:
  Zarate et al. 2006 (Am J Psychiatry 163:153-155) found memantine 5-20 mg/day no better than
  placebo in MDD (response 13% vs 13% at week 8). Preclinically, Gideons, Kavalali, Monteggia
  2014 (PNAS 111:8649) showed memantine
  does NOT trigger the eEF2-dephosphorylation / BDNF-induction that ketamine does. This is the
  linchpin negative: same target class, opposite plasticity outcome. PERSEUS L4b MUST score
  memantine window-NEGATIVE; if it scores positive, the window is wrong.

---

## 3. THE DISCRIMINATOR QUESTION: why ketamine but not memantine?

This is the heart of the design problem, and the literature gives a clear mechanistic answer
that is, unfortunately, only PARTIALLY structure-computable.

### 3a. The established functional discriminator: resting / use-dependent block

Gideons, Kavalali & Monteggia 2014 (PNAS 111:8649) is the decisive paper. Key results:
- In ZERO Mg2+, memantine and ketamine block NMDAR currents EQUALLY (similar IC50, similar
  potency). On raw affinity they are the same drug.
- In PHYSIOLOGICAL Mg2+, only ketamine effectively blocks NMDARs at REST (the NMDAR component
  of spontaneous miniature EPSCs). Memantine does not significantly block resting NMDARs when
  Mg2+ is present.
- Consequently only ketamine suppresses tonic eEF2K activity -> dephosphorylates eEF2 ->
  de-suppresses BDNF -> triggers the plasticity cascade. Memantine does not.

Mechanistic interpretation: the plasticity-relevant event is blockade of TONICALLY/SPONTANEOUSLY
active NMDARs at resting membrane potential, where Mg2+ still partially occludes the pore.
Ketamine reaches and blocks these low-activity channels; memantine, because of how it interacts
with the pore and a second/superficial site, effectively requires more sustained channel
opening (depolarization, high agonist) and so spares the resting pool. This dovetails with
Xia, Chen, Zhang & Lipton 2010 (J Neurosci 30:11246): at therapeutic 1-10 uM, memantine
PREFERENTIALLY blocks EXTRASYNAPTIC over synaptic NMDARs, being effective when glutamate is
elevated for minutes (tonic/pathological) but "relatively inactive" for millisecond synaptic
events. Memantine is tuned to chronic extrasynaptic over-activation (its Alzheimer's
rationale), the OPPOSITE of the transient spontaneous/synaptic block that triggers the
antidepressant cascade.

### 3b. The trapping-kinetics axis (mechanistically real, structurally HARD)

Why does memantine spare resting channels while ketamine does not? The leading explanation is
TRAPPING/UNTRAPPING and a second binding site:
- Ketamine is a near-FULL trapping blocker (stays bound when the channel closes/agonist leaves).
  Memantine is a PARTIAL trapping blocker: a fraction of bound memantine dissociates after
  agonist removal, via a more SUPERFICIAL second site than the deep pore site.
  - Partial trapping + superficial site: Kotermanski, Wood & Johnson 2009 (J Physiol
    587:4589); Blanpied et al. 1997 (J Neurophysiol 77:309) on amantadine/memantine trapping.
- Memantine also has a faster effective off-rate / second-site inhibition (SSI) not seen with
  ketamine; its lipophilic second-site reservoir may be the membrane/intracellular compartment
  (Kotermanski & Johnson 2009; recent SSI work, e.g. Glasgow/Mennerick-lineage studies).

The honest complication (a real scientific contest the project should represent):
- Emnett, Mennerick et al. 2013 (Mol Pharmacol 84:935, "Indistinguishable Synaptic
  Pharmacodynamics of ... Memantine and Ketamine") found the two drugs are LARGELY
  pharmacologically INDISTINGUISHABLE under steady-state and nonequilibrium conditions, with a
  SLIGHT difference in VOLTAGE-DEPENDENCE as the sole distinguishing feature, and argued the
  low open-probability of NMDARs MASKS pharmacodynamic differences under basal conditions.
  This directly tensions the "resting-block" discriminator: the difference is subtle and only
  surfaces under specific conditions. So the discriminator is REAL but SMALL, and depends on
  fine channel-state kinetics, not a gross structural feature.

### 3c. Candidate axes, scored for STRUCTURE-COMPUTABILITY

For each candidate discriminator, the question PERSEUS must answer: can it be derived from
structure / cheap descriptors, or does it need pharmacology PERSEUS cannot compute?

| Candidate axis | Separates ket from mem? | Structure-computable from SMILES / cheap ADMET? |
|---|---|---|
| Resting / use-dependent block (Mg2+-present block of spontaneous NMDARs) | YES (Gideons 2014) | NO. This is an electrophysiological channel-state property; not derivable from 2D structure. |
| Trapping fraction (full vs partial) | YES (ket full, mem partial) | NO as a transferable cross-scaffold descriptor. Within a congeneric series structure modulates it (Mealing et al. 2001 JPET, PMID 11356910: 23 analogs span a wide trapping range, off-rate tracks trapping), but the determinant is channel off-rate, a measured quantity, not a cheap descriptor. |
| Channel off-rate / residence time | YES (slower off-rate -> more trapping) | NO from cheap descriptors. Needs binding/electrophysiology data. |
| Superficial second site / SSI | YES (mem only) | NO. Site-specific pharmacology, not structural. |
| Lipophilicity / CNS kinetics (the L4 trick) | NO (see 3d) | YES but USELESS here: it does not separate them. |
| Active-metabolite formation (ket -> norketamine -> HNK; DXM -> dextrorphan) | PARTIAL / muddy | PARTIALLY. Metabolite identity can be curated; CYP-site prediction is weak. Not a clean discriminator (HNK contested). |
| GluN2B / extrasynaptic-vs-synaptic subtype preference | NO at the subtype level (see 3e) | Subtype affinity is NOT in MAMMAL's BindingDB-pKd head reliably; and the data say it is not the discriminator anyway. |
| Affinity (IC50) | NO (nearly equal) | Even if computable, does not separate them. |

### 3d. The decisive negative for the L4 strategy: lipophilicity does NOT separate them

The L4 serotonergic window works because serotonin is membrane-IMPERMEANT (TPSA 62, HBD 3)
while DMT/psilocin are permeant (TPSA <=54, HBD <=2), so a cheap TPSA/HBD gate separates
plastogen from non-plastogen DESPITE isoaffine binding. I computed the same descriptors
(RDKit, the repo's own engine) for the NMDA class:

| compound | clogP | TPSA | HBD | HBA | MW |
|---|---|---|---|---|---|
| ketamine | 2.90 | 29.1 | 1 | 2 | 237.7 |
| esketamine (S) | 2.90 | 29.1 | 1 | 2 | 237.7 |
| norketamine | 2.64 | 43.1 | 1 | 2 | 223.7 |
| (2R,6R)-HNK | 1.61 | 63.3 | 2 | 3 | 239.7 |
| memantine | 2.69 | 26.0 | 1 | 1 | 179.3 |
| amantadine | 1.91 | 26.0 | 1 | 1 | 151.3 |
| dextromethorphan | 3.38 | 12.5 | 0 | 2 | 271.4 |
| dextrorphan | 3.08 | 23.5 | 1 | 2 | 257.4 |
| phencyclidine | 4.33 | 3.2 | 0 | 1 | 243.4 |
| nitrous oxide | 0.34 | 51.2 | 0 | 2 | 44.0 |
| serotonin (ref) | 1.37 | 62.0 | 3 | 2 | 176.2 |

The verdict is unambiguous: memantine (clogP 2.69, TPSA 26.0) and ketamine (clogP 2.90, TPSA
29.1) are NEARLY IDENTICAL on every cheap physicochemical axis. If anything memantine is the
MORE CNS-shaped molecule (lower MW, lower TPSA, fewer HBA). All of these compounds clear the
L1 passive-CNS window. Therefore the property that makes L4 work (a permeability split off the
affinity axis) DOES NOT EXIST for the NMDA class. The thing that separates ketamine from
memantine is channel-block STATE-DEPENDENCE, which is invisible to 2D structure and to the L1
permeability gate.

NOTE on HNK descriptors: HNK is notably MORE polar than ketamine (TPSA 63.3, HBD 2). If one
naively tried a permeability gate, HNK would look the WORST of the ketamine series, which is
inconsistent with it being the supposed durable active species, another reason the L4-style
permeability logic is the wrong axis here.

### 3e. A second honest negative: GluN2B subtype preference is NOT the discriminator

A plausible-sounding hypothesis is "durable plastogens prefer GluN2B / extrasynaptic NMDARs."
The data refute this as the ketamine-vs-memantine discriminator: Kotermanski & Johnson 2009
(J Neurosci 29:2774) showed that in physiological Mg2+ BOTH memantine and ketamine acquire
preference for GluN2C/GluN2D-containing receptors (because GluN2A/2B are more Mg2+-sensitive),
i.e. they share the same subtype bias. So subtype preference does not separate them, and it is
not reliably structure-computable anyway. Recording this kills a tempting but wrong axis.

### 3f. What is computable vs what needs pharmacology PERSEUS cannot derive

- STRUCTURE-COMPUTABLE (cheap, from SMILES / existing ADMET): scaffold/chemotype membership
  (arylcyclohexylamine, aminoadamantane, morphinan), CNS-penetrance (already L1), basic amine
  presence, lipophilicity. NONE of these separate ketamine from memantine; they only define
  "is this an uncompetitive NMDA channel blocker that gets into the brain."
- NOT STRUCTURE-COMPUTABLE (needs measured pharmacology): trapping fraction, channel off-rate
  / residence time, resting-state vs use-dependent block, superficial-second-site behavior,
  voltage-dependence detail, synaptic-vs-extrasynaptic preference. These ARE the discriminator,
  and they are exactly what a sequence/structure model cannot derive (consistent with the
  repo's existing "MAMMAL is allosterically blind" and "off the DTI axis" findings).

---

## 4. PROPOSED L4b "NMDA-plasticity window" DESIGN (concrete, falsifiable)

Design philosophy, kept identical to L4: the window is PERMISSIVE (opens durable plasticity
only if downstream BDNF-TrkB-mTORC1-spine program fires), DIRECTION-NEUTRAL, and NEVER
auto-promotes to "demonstrated durable." It maps to the existing `plasticity_window` substrate
tier in `reversibility.py` and the `P_WINDOW` verdict in `perseus.py`. Crucially, because the
true discriminator is NOT structure-computable, L4b is built as a SCAFFOLD FLAG GATED BY A
SMALL CURATED PHARMACODYNAMIC TABLE, with ABSTAIN as the default. This is more conservative
than L4 by design.

### 4a. Gate logic (three states: window-positive / window-negative / abstain)

```
L4b_window(smiles, compound):
  # Step 1 (structure, cheap): is this an uncompetitive NMDA open-channel blocker chemotype
  #         that is CNS-penetrant? (necessary, NOT sufficient)
  scaffold = nmda_channel_blocker_scaffold(smiles)   # arylcyclohexylamine | aminoadamantane
                                                     # | morphinan-dissociative | (curated: gas)
  if scaffold is None:            return WINDOW_NEGATIVE (not an NMDA-channel-blocker chemotype)
  if not L1_cns_pass(smiles):     return WINDOW_NEGATIVE (no free-brain exposure)

  # Step 2 (the discriminator, CURATED pharmacology - cannot be derived from structure):
  pd = NMDA_PD_TABLE.get(canonical(compound))        # trapping, resting_block, use_dependence
  if pd is None:                  return ABSTAIN  (chemotype present but trapping/resting-block
                                                  pharmacology not measured -> cannot rule in
                                                  or out; abstain by default)

  # Step 3 (window decision from curated PD):
  if pd.full_or_near_full_trapping AND pd.blocks_resting_NMDAR:
        return WINDOW_POSITIVE  (permissive plasticity window; durable IFF downstream fires)
  if pd.partial_trapping AND pd.spares_resting_NMDAR (extrasynaptic-preferring):
        return WINDOW_NEGATIVE  (memantine-like: blocks the channel but not the plasticity-
                                relevant resting pool)
  else: return ABSTAIN
```

### 4b. The specific features

- STRUCTURE-DERIVED (RDKit, cheap), used ONLY to assign the scaffold and CNS gate (necessary
  conditions), NEVER to decide durability:
  - arylcyclohexylamine SMARTS (ketamine/PCP/tiletamine core: aryl-substituted cyclohexyl
    bearing an amine; ketone/hydroxyl variants covered for ket/norket/HNK),
  - aminoadamantane SMARTS (memantine/amantadine: adamantane cage + primary amine),
  - morphinan-dissociative (dextromethorphan/dextrorphan core),
  - a curated entry for non-druglike gases (N2O) that have no usable scaffold,
  - reuse `cns_exposure_gate` for the free-brain prerequisite (all these pass; it is a floor).
- CURATED PHARMACODYNAMIC TABLE (the load-bearing, NON-structure-computable part), e.g.
  `data/raw/nmda_pd_table.csv` with columns:
  `compound, trapping_class {full|near_full|partial}, blocks_resting_NMDAR {yes|no|unknown},
   use_dependence {low|high}, source_pmid`.
  Seed rows grounded in the papers above:
  - ketamine: near_full trapping, blocks_resting_NMDAR=yes (Gideons 2014; Kotermanski 2009)
  - esketamine: treat as ketamine-like (same channel pharmacology; curate explicitly)
  - memantine: partial trapping, blocks_resting_NMDAR=no, extrasynaptic-preferring
    (Kotermanski 2009; Xia/Lipton 2010; Gideons 2014)
  - amantadine: partial trapping (Blanpied 1997) -> window-negative/abstain
  - HNK: blocks_resting_NMDAR=unknown/contested -> ABSTAIN (do not assert positive; section 2)
  - DXM: multi-target, chronic-dosed -> ABSTAIN as a clean NMDA window (flag multi-mechanism)
  - N2O: trapping/resting-block not characterized like ket -> ABSTAIN
- DOWNSTREAM COUPLING (optional second gate, to honor the two-epoch logic): require the
  compound to also be consistent with triggering BDNF-TrkB-mTORC1 (the project already treats
  mTOR/BDNF as a DOWNSTREAM EFFECTOR axis, NOT a DTI target; do not score affinity for MTOR).
  In practice this stays a documented caveat, not a structural feature.

### 4c. PRE-REGISTERED NEGATIVE (the falsifier)

Register BEFORE running, mirroring the serotonin test in L4:
1. MEMANTINE must score WINDOW-NEGATIVE. (Same NMDA target class as ketamine, opposite
   plasticity outcome; Gideons 2014.) If L4b scores memantine positive, the window is
   FALSIFIED.
2. AMANTADINE must score WINDOW-NEGATIVE (partial trapper, not a rapid plastogen).
3. KETAMINE / esketamine must score WINDOW-POSITIVE.
4. A structure-only ablation must FAIL to separate memantine from ketamine: if you remove the
   curated PD table and rely on RDKit descriptors + L1 alone, memantine and ketamine must come
   out the SAME (both pass CNS, near-identical clogP/TPSA). This ablation is the POSITIVE PROOF
   that the discriminator lives off the structural axis (exactly analogous to "raw DTI cannot
   separate serotonin from DMT" in the L4 validation). Pre-register that the ablation CANNOT
   separate them; if it somehow does, the separation is an artifact and must be investigated.

### 4d. What MUST stay ABSTAIN (explicit honesty firewall)

- Any NMDA-channel-blocker chemotype whose trapping / resting-block pharmacology is NOT in the
  curated table -> ABSTAIN (not window-negative, because we cannot rule it out; not positive,
  because we cannot rule it in). This is the correct default given that the discriminator is
  not structure-derivable.
- HNK -> ABSTAIN on durability (contested mechanism; negative human translation).
- N2O -> ABSTAIN on durability (no controlled multi-week single-dose evidence).
- DXM / AXS-05 -> ABSTAIN as a CLEAN NMDA window (multi-target, chronic dosing, bupropion-
  driven exposure); efficacy of the drug-on state is established but is not single-exposure
  durability and the NMDA contribution is not isolable.
- NEVER promote any L4b-positive to "demonstrated durable": the single-dose human ketamine
  effect itself wanes within ~1 week and maintenance requires re-dosing (SUSTAIN-1). L4b
  positive == "opens a plasticity window," not "durable."

### 4e. Honest assessment of L4b versus L4

L4 is strong because its discriminator (permeability) is BOTH the true mechanism AND cheaply
structure-computable, so it adds real off-axis signal with one TPSA/HBD gate. L4b is WEAKER:
the true discriminator (trapping / resting-state block) is NOT cheaply structure-computable, so
L4b is mostly a curated-pharmacology lookup wrapped in a scaffold flag, defaulting to ABSTAIN.
Its scientific VALUE is therefore twofold and honest: (1) it correctly classifies the
hand-curated exemplars (ketamine positive, memantine negative) for the engine's recall on this
class WITHOUT pretending the call is structure-derived; (2) the structure-only ABLATION is a
publishable NEGATIVE result that extends the project's central thesis: for the NMDA RAAD class,
as for the serotonergic class, durability lives OFF the binding-affinity axis, but UNLIKE the
serotonergic class it is ALSO off the cheap-ADMET axis, residing in channel-state kinetics that
a structure/sequence model cannot derive. That boundary, cleanly demonstrated, is the
deliverable.

---

## Citations (author, year, venue; empirical claims only)

Mechanism of durability:
- Moghaddam B et al. 1997. J Neurosci 17:2921. (disinhibition / prefrontal glutamate surge)
- Li N, Lee B, Duman RS et al. 2010. Science 329:959. (mTORC1, synaptogenesis, rapamycin block)
- Autry AE, Adachi M, Monteggia LM et al. 2011. Nature 475:91. (eEF2K, BDNF de-suppression, TrkB)
- Nosyreva E et al. 2013. J Neurosci 33:6990. (spontaneous NMDAR block -> eEF2/BDNF)
- Duman RS & Aghajanian GK. 2012. Science 338:68. (synaptic plasticity review)
- Krystal JH, Abdallah CG, Sanacora G, Duman RS et al. 2019. Neuron 101:774. (mechanism review)
- Moda-Sava RN et al. 2019. Science 364:eaat8078. (spine formation sustains the effect)
- Ly C et al. 2021. (two-epoch stimulation/growth; cited in repo gaps doc)

Honesty gate (per-compound durability):
- Berman RM et al. 2000. Biol Psychiatry 47:351. (first ketamine RCT signal)
- Zarate CA et al. 2006. Arch Gen Psychiatry 63:856. (ketamine TRD RCT, ~1 week)
- Daly EJ et al. 2019 (SUSTAIN-1). JAMA Psychiatry 76:893-903. (esketamine relapse-prevention,
  randomized-withdrawal; continued dosing cut relapse ~51% (responders) / ~70% (remitters) vs
  switch-to-placebo -> benefit maintained by ongoing dosing, not a one-time durable reset)
- Zarate CA et al. 2006. Am J Psychiatry 163:153-155. (memantine NOT effective in MDD; the
  clinical negative control)
- Zanos P, Gould TD et al. 2016. Nature 533:481. (HNK; NMDAR-independent claim)
- Suzuki K, Nosyreva E, Kavalali ET, Monteggia LM. 2017. Nature 546:E1. (HNK DOES block NMDARs; reply)
- Lumsden EW et al. 2019. PNAS 116:5160. (antidepressant-relevant HNK conc. do NOT block NMDAR)
- Zanos P et al. 2019. PNAS 116:6441. (HNK mGlu2-dependent mechanism)
- Farmer CA et al. 2020. Neuropsychopharmacol 45:1398. (higher plasma HNK -> worse response)
- Nagele P et al. 2015. Biol Psychiatry 78:10. (N2O proof-of-concept)
- Nagele P et al. 2021. Sci Transl Med 13:eabe1376. (N2O phase 2; effect up to ~1 week, beyond unknown)
- Nagele P et al. 2020. Front Psychiatry 11:692. (single-case 2-month remission; anecdotal)
- Iosifescu DV et al. 2022 (GEMINI). J Clin Psychiatry 83:21m14345. (AXS-05 DXM-bupropion MDD)

Discriminator (memantine vs ketamine):
- Gideons ES, Kavalali ET, Monteggia LM. 2014. PNAS 111:8649. (resting-NMDAR block, eEF2/BDNF;
  THE discriminator paper)
- Xia P, Chen HV, Zhang D, Lipton SA. 2010. J Neurosci 30:11246. (memantine extrasynaptic-
  preferring; inactive at ms synaptic events)
- Kotermanski SE, Wood JT, Johnson JW. 2009. J Physiol 587:4589. (memantine partial trapping,
  superficial second site)
- Kotermanski SE, Johnson JW. 2009. J Neurosci 29:2774. (Mg2+-dependent GluN2C/2D selectivity,
  SHARED by both drugs -> subtype is NOT the discriminator)
- Blanpied TA, Boeckman FA, Aizenman E, Johnson JW. 1997. J Neurophysiol 77:309.
  (amantadine/memantine trapping channel block)
- Emnett CM, Mennerick S et al. 2013. Mol Pharmacol 84:935. (memantine/ketamine LARGELY
  INDISTINGUISHABLE; voltage-dependence the sole distinguisher; the honest counter-evidence)
- Mealing GAR et al. 2001. J Pharmacol Exp Ther (PMID 11356910). (23 analogs; trapping is
  graded and tracks off-rate; trapping is structure-modulable within a series but determined by
  measured off-rate, not cheap descriptors)

Descriptor table: computed in-repo with RDKit (Crippen clogP, TPSA, HBD/HBA) via the project's
own `cns_exposure` descriptor stack; values are deterministic RDKit outputs, not literature
citations.

Uncertainty flags: I did NOT independently re-verify every page number against the primary PDF
for papers behind paywalls (PNAS/Science/Nature returned 403 to the fetch tool); volume/page
were taken from indexer metadata and secondary sources and should be spot-checked before
formal citation. The HNK replication-failure specifics (exact models, exact effect sizes) are
reported as "contested" rather than with invented numbers. No trial outcome, effect size, or
citation in this document was fabricated; where I was unsure of a number I said so.
