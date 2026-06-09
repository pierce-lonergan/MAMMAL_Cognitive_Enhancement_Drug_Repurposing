# L4b NMDA trapping-kinetics curated lookup table (citation-verified)

Purpose: a VERIFIED, per-compound curated pharmacodynamic table for the PERSEUS L4b
"NMDA-plasticity / durability window." It encodes the one discriminator that the lane-A research
report (`research_nmda_durability_window.md`) established as decisive and NON-structure-computable:
channel TRAPPING / resting-state, use-dependent block. The scientific thesis (Gideons 2014): the
property that separates a durable rapid-antidepressant NMDA antagonist (ketamine) from a
non-durable one (memantine) is whether the blocker reaches and blocks resting/spontaneous NMDARs
in physiological Mg2+ (the BDNF-triggering property), driven by full vs partial trapping, NOT by
2D structure (ketamine and memantine are near-identical on clogP/TPSA/HBD).

This document is a citation-verification pass. Every PMID/DOI below was looked up directly on
PubMed, the publisher, or Europe PMC during curation (NOT taken from indexer metadata). Where a
claim could not be verified it is marked UNVERIFIED with the reason. No PMID, DOI, or numeric
value here was fabricated.

## Pre-registered window rule (applied verbatim)

```
WINDOW   iff (blocks_resting_NMDAR == yes) AND (durable_rapid_antidepressant in {established, preclinical_only})
NEGATIVE iff (durable_rapid_antidepressant == established_negative) OR (blocks_resting_NMDAR == no)
ABSTAIN  otherwise   (e.g. a trapper whose post-cessation durability is not established: PCP, MK-801, HNK, N2O)
```

Precedence note: NEGATIVE is evaluated as a hard floor. If a compound is both
established_negative AND would otherwise abstain, it is NEGATIVE. If blocks_resting_NMDAR == no it
is NEGATIVE regardless of the antidepressant column. WINDOW requires BOTH the resting-block AND a
positive durability tier. Everything else ABSTAINS. This makes ABSTAIN the honest default for
"trapper, but durability after cessation not demonstrated."

---

## VERIFIED CITATION LEDGER (looked up this pass)

Mechanism / trapping / resting-block primary sources:

| short | PMID | DOI | venue (verified) | what it establishes |
|---|---|---|---|---|
| Gideons 2014 | 24912158 | 10.1073/pnas.1323920111 | PNAS 111(23):8649-8654 | THE discriminator: in physiological Mg2+ ketamine blocks resting NMDAR currents and drives eEF2 dephos / BDNF; memantine is a poor resting-NMDAR blocker and does neither |
| Kotermanski & Johnson 2009 | 19261873 | 10.1523/JNEUROSCI.3703-08.2009 | J Neurosci 29(9):2774-2779 | physiological 1 mM Mg2+ cuts memantine block ~20-fold near rest; in Mg2+ both memantine AND ketamine shift toward GluN2C/2D preference (subtype is SHARED, not the discriminator) |
| Kotermanski, Wood & Johnson 2009 | 19687120 | 10.1113/jphysiol.2009.176297 | J Physiol 587(Pt 19):4589-4603 | memantine is a PARTIAL trapper via a superficial second (non-trapping) site; ketamine is essentially FULLY trapped (little superficial-site binding) |
| Blanpied 1997 | 9120573 | 10.1152/jn.1997.77.1.309 | J Neurophysiol 77(1):309-323 | amantadine AND memantine are trapping channel blockers (trapped after channel closure/agonist unbinding); IC50 amantadine 39 uM, memantine 1.4 uM at -67 mV |
| Mealing 1999 | 9862772 | (no DOI on record) | J Pharmacol Exp Ther 288(1):204-210 | three blockers with SIMILAR block kinetics differ in DEGREE of trapping (AR-R15896AR/lanicemine-class low-trapping vs ketamine/memantine); trapping is graded, not all-or-none |
| Mealing 2001 | 11356910 | (no DOI on record) | J Pharmacol Exp Ther 297(3):906-914 | 23 structural analogs span a wide trapping range; off-rate positively correlates with trapping; trapping is structure-modulable WITHIN a series but determined by measured off-rate |
| Huettner & Bean 1988 | 2448800 | (no DOI on record) | PNAS 85(4):1307-1311 | MK-801 selectively binds OPEN NMDA channels and is TRAPPED when the channel closes; block persists long after washout (canonical full-trapping demonstration; also the PCP-site/foundational reference) |
| Sanacora 2014 (lanicemine) | 24126931 | 10.1038/mp.2013.130 | Mol Psychiatry 19(9):978-985 | lanicemine (AZD6765) is a LOW-trapping NMDA channel blocker: at steady state ketamine ~86% trapped vs lanicemine ~54%; low trapping linked to minimal psychotomimetic effect |
| Jevtovic-Todorovic 1998 (N2O) | 9546794 | 10.1038/nm0498-460 | Nat Med 4(4):460-463 | nitrous oxide is a noncompetitive NMDA antagonist; its block is much FASTER and more easily REVERSIBLE than ketamine's (i.e. fast off-rate, low-trapping behavior) |
| Werling 2007 (DXM) | 17848867 | 10.1016/j.expneurol.2007.06.027 | Exp Neurol 207(2):248-257 | dextromethorphan is a LOW-affinity uncompetitive NMDA channel blocker (plus sigma-1 agonist, multi-target) |

Per-compound durability / clinical primary sources:

| short | PMID | DOI | venue (verified) | what it establishes |
|---|---|---|---|---|
| Daly 2019 (SUSTAIN-1) | 31166571 | 10.1001/jamapsychiatry.2019.1189 | JAMA Psychiatry 76(9):893-903 | esketamine nasal spray + oral AD cut relapse risk ~51% vs placebo (randomized-withdrawal); benefit MAINTAINED by ongoing dosing |
| Zarate 2006 (memantine MDD) | 16390905 | 10.1176/appi.ajp.163.1.153 | Am J Psychiatry 163(1):153-155 | memantine 5-20 mg/day NOT effective in MDD vs placebo (no treatment effect on MADRS) -- the clinical NEGATIVE control |
| Zanos 2016 (HNK) | 27144355 | 10.1038/nature17998 | Nature 533:481-486 | (2R,6R)-HNK reproduces ketamine-like antidepressant actions in mice and CLAIMS they are NMDAR-inhibition-INDEPENDENT (AMPAR-dependent) |
| Suzuki 2017 (HNK reply) | 28640258 | 10.1038/nature22084 | Nature 546(7659):E1-E3 | rebuttal: HNK DOES inhibit synaptic NMDARs and couples signaling to NMDAR inhibition (mechanism CONTESTED) |
| Lumsden 2019 (HNK) | 30796190 | 10.1073/pnas.1816071116 | PNAS 116(11):5160-5169 | antidepressant-RELEVANT (low) concentrations of (2R,6R)-HNK do NOT block NMDAR function (resting-block at AD-relevant dose: NO/unknown) |
| Zanos 2019 (HNK mGlu2) | (not separately re-verified this pass) | 10.1073/pnas.1819540116 | PNAS 116:6441-... | proposes HNK acts via an mGlu2-dependent presynaptic mechanism (alternative, non-NMDA target) |
| Nagele 2021 (N2O) | 34108247 | 10.1126/scitranslmed.abe1376 | Sci Transl Med 13(597):eabe1376 | N2O phase 2 crossover (n=24): 25% and 50% N2O improved TRD vs placebo; effect observed over "several weeks"; durability BEYOND ~1 week from a single session not established by design |
| Sanacora 2017 (lanicemine ph2b) | 27681442 | 10.1016/j.jad.2016.08.058 | J Affect Disord (2017) | large adjunctive phase 2b: repeated-dose lanicemine showed NO superiority to placebo at any timepoint; clinical development terminated for lack of efficacy |

DOIs marked "(no DOI on record)" are pre-DOI-era articles (1988-2001); they are cited by PMID,
which is the stable identifier. The Werling 2007 DOI and the Sanacora 2014/2017 volume-page
details were taken from the publisher/secondary listings during search and are flagged below as
needing a one-line spot-check before any formal manuscript use; the PMIDs themselves are verified.

---

## PER-COMPOUND REASONING

### 1. Ketamine (racemic, R,S)
- trapping_class: full_trapping. Ketamine is essentially fully trapped (Kotermanski/Wood/Johnson
  2009, J Physiol, PMID 19687120: little superficial-site binding, so bound ketamine stays when
  the channel closes; Mealing 1999 PMID 9862772: high trapping relative to low-trappers).
- blocks_resting_NMDAR: yes. Gideons 2014 (PMID 24912158): in physiological Mg2+ ketamine blocks
  the NMDAR component of spontaneous/resting transmission and drives eEF2 dephosphorylation ->
  BDNF de-suppression.
- use_dependence: moderate. Uncompetitive open-channel blocker (needs channel opening) but its
  resting-block at physiological Mg2+ means it is NOT purely high-use-dependent like a pure
  open-channel-only blocker; it reaches the low-activity/resting pool.
- durable_rapid_antidepressant: established. Single-dose rapid effect replicated; esketamine
  randomized-withdrawal maintenance (Daly 2019, PMID 31166571). (Honest caveat from lane A:
  durable for ~days-to-1-week from one dose; sustained remission needs re-dosing. The pre-reg
  tier "established" is about the rapid-antidepressant-with-plasticity claim, which is solid.)
- window_verdict: WINDOW (blocks_resting_NMDAR yes AND durable established).
- best trapping/resting primary: PMID 24912158 (Gideons 2014).

### 2. Esketamine (S-ketamine)
- trapping_class: full_trapping (same channel pharmacology as racemate; the S-enantiomer is the
  higher-NMDAR-affinity enantiomer). Curated as ketamine-like; no enantiomer-specific trapping
  paper re-verified this pass (flag: trapping classed by read-across to ketamine, which is
  standard and well supported, but is read-across not a separate measurement).
- blocks_resting_NMDAR: yes (read-across to ketamine; Gideons 2014 used racemic/S-active block).
- use_dependence: moderate (as ketamine).
- durable_rapid_antidepressant: established. FDA-approved (2019) for TRD; SUSTAIN-1 relapse
  prevention (Daly 2019, PMID 31166571).
- window_verdict: WINDOW.
- best trapping/resting primary: PMID 24912158 (Gideons 2014).

### 3. Arketamine (R-ketamine)
- trapping_class: full_trapping (arylcyclohexylamine channel blocker, same pore site as ketamine;
  classed by read-across). NOTE the genuine nuance: arketamine has WEAKER NMDAR affinity than
  esketamine yet equal/greater preclinical antidepressant effect, which is part of the argument
  that its antidepressant action may be partly NMDAR-INDEPENDENT (Hashimoto-lineage work). So the
  trapping class is read-across at the channel; the antidepressant mechanism is debated.
- blocks_resting_NMDAR: yes (read-across to ketamine channel block; lower affinity than S). Flag:
  not independently measured at resting state in the verified sources; assigned by congener
  read-across with a contested-mechanism caveat.
- use_dependence: moderate.
- durable_rapid_antidepressant: preclinical_only. Greater/longer-lasting antidepressant than
  esketamine in rodent models with fewer side effects; human RCT evidence still emerging, not
  established to the esketamine standard. (Per pre-reg rule, preclinical_only is a WINDOW-eligible
  tier.)
- window_verdict: WINDOW (blocks_resting_NMDAR yes AND durable preclinical_only). Flagged as the
  weakest WINDOW call because both the resting-block and the trapping are read-across, and the
  mechanism is argued by some to be NMDAR-independent. If PERSEUS prefers maximal conservatism,
  arketamine is the one row a reviewer could legitimately move to ABSTAIN; under the rule as
  written it is WINDOW.
- best trapping/resting primary: PMID 24912158 (Gideons 2014, parent-compound resting block).

### 4. (2R,6R)-hydroxynorketamine (HNK)
- trapping_class: low_trapping. As a much lower-affinity, weak channel blocker its trapping is
  low; but more importantly its channel action at antidepressant-relevant concentrations is
  itself in dispute (see below). Classed low_trapping with a contested flag.
- blocks_resting_NMDAR: unknown (genuinely CONTESTED). Zanos 2016 (PMID 27144355) says HNK's
  antidepressant action is NMDAR-INHIBITION-INDEPENDENT; Suzuki 2017 (PMID 28640258) says HNK
  DOES inhibit synaptic NMDARs; Lumsden 2019 (PMID 30796190) shows AD-relevant (low) HNK
  concentrations do NOT block NMDAR. Net: resting-block at therapeutic exposure is not
  established -> unknown.
- use_dependence: low (weak/uncertain channel block).
- durable_rapid_antidepressant: contested. Preclinical signal disputed across labs; higher human
  plasma HNK after ketamine has been associated with WORSE response; human durability essentially
  unestablished.
- window_verdict: ABSTAIN. The rule: not established_negative, and blocks_resting_NMDAR is unknown
  (not "no"), and durability is "contested" (not in {established, preclinical_only}) -> ABSTAIN.
  This is the correct honesty firewall outcome.
- best trapping/resting primary: PMID 30796190 (Lumsden 2019 -- AD-relevant conc do not block
  NMDAR), with PMID 27144355 / PMID 28640258 as the two sides of the contest.

### 5. Memantine (THE pre-registered NEGATIVE)
- trapping_class: partial_trapping. Kotermanski/Wood/Johnson 2009 (PMID 19687120): memantine
  binds a superficial second (non-trapping) site so a fraction dissociates after agonist removal
  = partial trapping. Blanpied 1997 (PMID 9120573) established memantine as a trapping blocker;
  the 2009 J Physiol work refined it to PARTIAL.
- blocks_resting_NMDAR: no. Gideons 2014 (PMID 24912158): memantine is a POOR blocker of resting
  NMDAR currents in physiological Mg2+ and does NOT trigger eEF2 dephos / BDNF. Kotermanski &
  Johnson 2009 (PMID 19261873): physiological Mg2+ cuts memantine block ~20-fold near rest.
- use_dependence: high. Effectively requires sustained channel opening (depolarization / elevated
  glutamate); tuned to chronic extrasynaptic over-activation, spares the resting/spontaneous pool.
- durable_rapid_antidepressant: established_negative. Zarate 2006 (PMID 16390905): memantine not
  effective in MDD vs placebo.
- window_verdict: NEGATIVE (established_negative AND blocks_resting_NMDAR no -- doubly negative;
  the linchpin falsifier of the whole window).
- best trapping/resting primary: PMID 24912158 (Gideons 2014); PMID 19687120 (partial trapping).

### 6. Amantadine
- trapping_class: partial_trapping. Blanpied 1997 (PMID 9120573) showed amantadine is a trapping
  channel blocker (low affinity, IC50 39 uM); later work (Blanpied/Johnson lineage) characterizes
  it as accelerating channel closure / partial trapping akin to but weaker than memantine. Classed
  partial_trapping by read-across to the aminoadamantane mechanism with the Blanpied trapping
  measurement.
- blocks_resting_NMDAR: no. Read-across to memantine (same aminoadamantane low-affinity,
  Mg2+-sensitive, use-dependent class; even weaker affinity than memantine). Flag: not directly
  measured at resting state in the verified sources; assigned no by class read-across + much lower
  affinity than ketamine.
- use_dependence: high (low-affinity uncompetitive blocker; even more use-dependent than
  memantine given ~28x lower potency).
- durable_rapid_antidepressant: not_studied (no controlled rapid-antidepressant evidence; used
  for Parkinson/influenza, not a validated rapid antidepressant). Pre-reg falsifier #2 expects
  amantadine WINDOW-NEGATIVE.
- window_verdict: NEGATIVE (blocks_resting_NMDAR no -> NEGATIVE by the rule). Satisfies the
  pre-registered "amantadine must be window-negative" falsifier.
- best trapping/resting primary: PMID 9120573 (Blanpied 1997).

### 7. Dextromethorphan (DXM)
- trapping_class: low_trapping. Low-affinity uncompetitive NMDA channel blocker (Werling 2007,
  PMID 17848867); low-affinity uncompetitive blockers in this class have fast off-rates and low
  trapping (Mealing 1999/2001). Multi-target (sigma-1 agonist, SERT inhibitor).
- blocks_resting_NMDAR: unknown. Not characterized for resting-state block the way ketamine is in
  the verified sources; its low affinity makes strong resting-pool block unlikely but this was not
  directly measured -> unknown.
- use_dependence: moderate-to-high (low-affinity uncompetitive); recorded as moderate with a note
  that block is open-channel/use-dependent and confounded by non-NMDA targets.
- durable_rapid_antidepressant: contested. Efficacy is established for the chronically-dosed
  DXM-bupropion COMBINATION (AXS-05), but that is multi-target, continuously dosed, and
  bupropion-boosted; clean single-exposure post-cessation NMDA durability is not demonstrated and
  the NMDA contribution is not isolable. Recorded as contested (muddy), not established.
- window_verdict: ABSTAIN. blocks_resting_NMDAR unknown (not "no") and durability contested (not
  established/preclinical_only) -> ABSTAIN. (Multi-mechanism, chronic dosing flag.)
- best trapping/resting primary: PMID 17848867 (Werling 2007).

### 8. Nitrous oxide (N2O)
- trapping_class: low_trapping. Jevtovic-Todorovic 1998 (PMID 9546794): N2O is a noncompetitive
  NMDA antagonist whose block is much FASTER and more easily REVERSIBLE than ketamine's -- i.e.
  fast off-rate / low-trapping behavior, the opposite of a full trapper. (It is a gas with no
  arylcyclohexylamine/aminoadamantane scaffold; channel-block kinetics, not structure, define it.)
- blocks_resting_NMDAR: unknown. Its rapid, readily-reversible noncompetitive block has not been
  characterized as a sustained resting-NMDAR block in the Gideons sense in the verified sources ->
  unknown.
- use_dependence: moderate (noncompetitive; fast on/off).
- durable_rapid_antidepressant: not_studied (for durability). Rapid effect is real (Nagele 2021,
  PMID 34108247) but multi-week single-session durability is explicitly NOT established by the
  trial design; the "2-month remission" claim is a single anecdotal case report. Recorded as
  not_studied for POST-CESSATION durability (rapid effect established; durability unestablished).
- window_verdict: ABSTAIN. blocks_resting_NMDAR unknown (not "no"); durability not established ->
  ABSTAIN. (If one took the strict view that a fast, easily-reversible blocker plainly does NOT
  block the resting pool, N2O would flip to NEGATIVE via blocks_resting_NMDAR=no; the verified
  evidence does not directly assert resting-block either way, so the honest call is unknown ->
  ABSTAIN. Flagged for reviewer.)
- best trapping/resting primary: PMID 9546794 (Jevtovic-Todorovic 1998).

### 9. Phencyclidine (PCP)
- trapping_class: full_trapping. Canonical high-affinity open-channel blocker trapped on channel
  closure (Huettner & Bean 1988, PMID 2448800, the foundational MK-801/PCP-site trapping work;
  MacDonald et al. 1991 J Physiol corroborates pronounced trapping for PCP).
- blocks_resting_NMDAR: unknown. Not characterized for resting/spontaneous block in the Gideons
  sense in the verified sources (PCP is studied as a psychotomimetic/neurotox tool, not an
  antidepressant) -> unknown.
- use_dependence: high (classic high-affinity open-channel-only blocker; needs channel opening).
- durable_rapid_antidepressant: not_studied (psychotomimetic drug of abuse; never developed as a
  therapeutic antidepressant).
- window_verdict: ABSTAIN. A full trapper but durability not established and resting-block
  unknown -> ABSTAIN (exactly the "trapper, durability not established" abstain case named in the
  pre-reg rule).
- best trapping/resting primary: PMID 2448800 (Huettner & Bean 1988).

### 10. Dizocilpine (MK-801)
- trapping_class: full_trapping. THE canonical full-trapping blocker: Huettner & Bean 1988 (PMID
  2448800) showed MK-801 selectively binds open channels and is trapped on closure, block
  persisting long after washout.
- blocks_resting_NMDAR: unknown. Not characterized as a resting/spontaneous-pool blocker in the
  Gideons antidepressant sense in the verified sources (used as a pharmacological tool) -> unknown.
- use_dependence: high (high-affinity open-channel blocker).
- durable_rapid_antidepressant: not_studied (research tool / not a therapeutic).
- window_verdict: ABSTAIN. Full trapper, durability not established -> ABSTAIN.
- best trapping/resting primary: PMID 2448800 (Huettner & Bean 1988).

### 11. Lanicemine (AZD6765)
- trapping_class: low_trapping. Sanacora 2014 (PMID 24126931): explicitly a LOW-trapping NMDA
  channel blocker; at steady state ketamine ~86% trapped vs lanicemine ~54% (low trapping was the
  design rationale for fewer psychotomimetic effects). Mealing 1999 (PMID 9862772) characterized
  the AR-R15896AR/lanicemine low-trapping class.
- blocks_resting_NMDAR: no/low. Low trapping is posited to PRESERVE use-dependent block and bias
  block toward high-tonic-activity elements rather than the resting pool; combined with its
  clinical failure this is recorded as no (does not effectively block the resting/plasticity pool
  the way ketamine does). Flag: "no" here is inferred from the low-trapping mechanism +
  clinical-negative, not from a direct Gideons-style resting-mEPSC measurement of lanicemine.
- use_dependence: high (low-trapping = block preferentially retained under ongoing/tonic
  activity; spares low-activity/resting channels).
- durable_rapid_antidepressant: established_negative. The large adjunctive phase 2b (Sanacora
  2017, PMID 27681442) showed NO separation from placebo and development was terminated. (An
  earlier small study suggested benefit; the definitive larger trial was negative.) This makes
  lanicemine a SECOND clinical-negative control alongside memantine, and a clean test that
  "NMDA channel block alone, without ketamine-like trapping/resting-block, is not enough."
- window_verdict: NEGATIVE (established_negative AND blocks_resting_NMDAR no). Reinforces the
  thesis: low-trapping blocker that spares the resting pool and clinically failed.
- best trapping/resting primary: PMID 24126931 (Sanacora 2014).

---

## SUMMARY TABLE

| compound | trapping_class | blocks_resting_NMDAR | use_dependence | durable_rapid_antidepressant | window_verdict | best PMID/DOI |
|---|---|---|---|---|---|---|
| ketamine (racemic) | full_trapping | yes | moderate | established | WINDOW | PMID 24912158 |
| esketamine (S) | full_trapping | yes | moderate | established | WINDOW | PMID 24912158 |
| arketamine (R) | full_trapping | yes (read-across; contested) | moderate | preclinical_only | WINDOW | PMID 24912158 |
| (2R,6R)-HNK | low_trapping | unknown (contested) | low | contested | ABSTAIN | PMID 30796190 |
| memantine | partial_trapping | no | high | established_negative | NEGATIVE | PMID 24912158 |
| amantadine | partial_trapping | no (read-across) | high | not_studied | NEGATIVE | PMID 9120573 |
| dextromethorphan | low_trapping | unknown | moderate | contested | ABSTAIN | PMID 17848867 |
| nitrous oxide | low_trapping | unknown | moderate | not_studied | ABSTAIN | PMID 9546794 |
| phencyclidine (PCP) | full_trapping | unknown | high | not_studied | ABSTAIN | PMID 2448800 |
| dizocilpine (MK-801) | full_trapping | unknown | high | not_studied | ABSTAIN | PMID 2448800 |
| lanicemine (AZD6765) | low_trapping | no | high | established_negative | NEGATIVE | PMID 24126931 |

Pre-registered falsifier check (from the lane-A design): memantine NEGATIVE (pass), amantadine
NEGATIVE (pass), ketamine/esketamine WINDOW (pass). Lanicemine adds a SECOND clinical-negative
that the rule also scores NEGATIVE for the right mechanistic reason (low-trapping, spares resting
pool, failed in trials). PCP/MK-801/HNK/N2O all ABSTAIN as intended (trappers or weak blockers
whose post-cessation durability is not established).

---

## INTEGRITY / UNCERTAINTY FLAGS (read before trusting any single cell)

VERIFIED THIS PASS (PMID confirmed on PubMed / Europe PMC / publisher):
24912158, 19261873, 19687120, 9120573, 9862772, 11356910, 2448800, 24126931, 27681442,
27144355, 28640258, 30796190, 9546794, 34108247, 17848867, 31166571, 16390905.

READ-ACROSS (not an independent per-compound measurement; assigned by congener class, standard
but flagged):
- esketamine trapping + resting-block: read-across to racemic/ketamine (S is the higher-affinity
  enantiomer; well supported but not a separate esketamine electrophysiology paper here).
- arketamine trapping + resting-block: read-across to ketamine channel block; its antidepressant
  MECHANISM is contested (argued partly NMDAR-independent despite weaker affinity). Weakest WINDOW
  call; a conservative reviewer could move it to ABSTAIN.
- amantadine blocks_resting_NMDAR=no: read-across to the aminoadamantane class + its ~28x lower
  potency than memantine; not a direct resting-mEPSC measurement in the verified sources.
- lanicemine blocks_resting_NMDAR=no: inferred from low-trapping mechanism (PMID 24126931) +
  clinical-negative (PMID 27681442), not a direct Gideons-style resting-block recording of
  lanicemine.
- DXM, N2O, PCP, MK-801 blocks_resting_NMDAR=unknown: none has a verified Gideons-style
  resting/spontaneous-NMDAR antidepressant-context measurement; "unknown" is the honest call and
  it routes them to ABSTAIN (except where established_negative would override, which is not the
  case for these four).

NUMERIC CLAIMS used (all from verified sources, quoted not invented): Mg2+ cuts memantine block
~20-fold near rest (PMID 19261873); IC50 amantadine 39 uM / memantine 1.4 uM at -67 mV (PMID
9120573); ketamine ~86% vs lanicemine ~54% trapping (PMID 24126931); esketamine relapse risk
reduction ~51% (PMID 31166571). The lane-A report's "memantine MDD response 13% vs 13%" figure was
NOT reproduced here; the verified Zarate 2006 abstract states "no treatment effect / not
effective" (PMID 16390905) without my confirming that specific percentage, so only the qualitative
negative is asserted.

SPOT-CHECK BEFORE MANUSCRIPT (PMID verified; volume/page or DOI taken from publisher/secondary
listing during search, low risk but not double-keyed): Sanacora 2014 Mol Psychiatry volume/pages
(19:978-985) and DOI 10.1038/mp.2013.130; Sanacora 2017 J Affect Disord volume/pages and DOI
10.1016/j.jad.2016.08.058; Werling 2007 Exp Neurol volume/pages (207:248-257) and DOI; Zanos 2019
PNAS (10.1073/pnas.1819540116) was carried from the lane-A report and not re-verified this pass.

NOT FABRICATED: no PMID, DOI, journal, or number above was invented. Where a value or mechanism
was not directly verifiable it is marked unknown / read-across / spot-check, and the window_verdict
was derived only from the pre-registered rule applied to those honest values.

CONCEPTUAL CAVEAT (carried from lane A, still true): this table is the LOAD-BEARING,
NON-structure-computable part of L4b. RDKit descriptors cannot separate ketamine (clogP 2.90 /
TPSA 29.1) from memantine (clogP 2.69 / TPSA 26.0); the discriminator encoded here lives in
channel-state kinetics (trapping / resting-block), which a structure or sequence model cannot
derive. A WINDOW verdict means "opens a plasticity window," never "demonstrated durable" -- even
ketamine's single-dose human effect wanes within ~1 week and maintenance requires re-dosing
(PMID 31166571).
