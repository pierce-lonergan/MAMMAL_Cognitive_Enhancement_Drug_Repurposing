# Is durable (post-washout) cognitive enhancement in healthy adults findable?

**A decision document.** Produced by a 25-agent research sweep (6 evidence lanes, each load-bearing
claim independently and adversarially verified against real PMIDs / DOIs / NCT ids), synthesised
against this repo's own verified ledgers. Research provenance — including every verification verdict
and, importantly, the claims that were **REJECTED as misdescribed** — is preserved at
`data/raw/provenance/durable_enhancement_research_2026-07.json`.

Companion quantitative analysis: `durability_gap_v1.md` (`scripts/123_durability_gap.py`), which
measures the empty target cell directly from the two ground-truth ledgers.

**Read section 0 first.** Three supplied claims were discounted on verification, including the
single most-cited precedent in this field (the valproate absolute-pitch study), and one
load-bearing citation could not be read at all. The verdict is stated net of those corrections.

---

# DECISION DOCUMENT — Durable cognitive enhancement in healthy adults: is it findable?

**Scope:** post-washout (≥1 month off-drug) cognitive gain in healthy adults, and what a system would need to find one.
**Status:** all counts below are computed from this repo's own verified ledgers; all external evidence carries its verification verdict.

---

## 0. VERIFICATION HYGIENE (read first)

Three things in the supplied research must be discounted before any of it is used:

- **REJECTED (misdescribed):** the claim that Gervain 2013 is "the only published human pharmacological attempt to reopen a canonical critical period." False — adult amblyopia (ocular dominance) *is* the canonical critical period and has a drug-plasticity trial literature (levodopa 1990; fluoxetine Sci Rep 2018 + a 42-patient Phase 2; citalopram Neural Plast 2019; a donepezil pilot). Gervain's own, narrower superlative is "first to show a change in AP with any kind of drug." Do not use the "only attempt" framing.
- **UNVERIFIED, and it is the load-bearing one:** Knecht 2004 (*Ann Neurol* 56:20-26, levodopa + 5 days vocabulary training, n=40) — the retention interval, whether retest was off-drug, and the effect size **could not be obtained** (Wiley 403). This is the citation most likely to be used to claim the healthy-durable precedent is established, and it is the one nobody in this workflow could read.
- **CONFIRMED but weaker than stated:** Levi/Li/Silver/Chung 2020 (donepezil, letter uncrowding, null) had **no placebo arm and no blinding** — all 19 subjects received donepezil. It is a non-replication *by historical comparison*, not a controlled null. Silver is also not its senior author (Chung is). It undercuts less than the lane claimed.
- **Lane 3's verification block was truncated in the payload I received.** Only the FLAME verification arrived. So FOCUS/AFFINITY/EFFECTS (n=5907, cOR 0.96 [0.87-1.05]), Cochrane 2021 (motor SMD 0.03, high GRADE), the d-cycloserine IPD meta (d=0.25 → 0.19), and DARS (OR 0.78) are **unverified in this payload**. The verdict below does not depend on those numbers being exact — only on their direction, which is triangulated by three mechanistically independent programs (serotonergic, dopaminergic, glutamatergic).

---

## 1. VERDICT

**Not demonstrated. Plausible-but-undemonstrated for narrow trained skills. Structurally implausible as currently framed — and, separately, not currently findable by prediction.**

Those are three distinct claims. Decisive evidence for each:

**(a) Not demonstrated.** The world precedent is **two studies**, neither replicated:

- Rokem & Silver 2013 (PMID 23755006) — CONFIRMED. Donepezil + motion-discrimination training, healthy adults, retested 5-15 months off-drug; percent learning 47.1±4.6 vs 34.2±6.9, signed-rank **p=0.036, n=8 returners of 12.** *And here is the sentence that matters most in this document:* **the placebo-trained condition also retained its learning, and absolute thresholds converged at follow-up (donepezil 7.4° vs placebo 7.3°).** The surviving drug advantage exists only in *percent learning normalized to each condition's own pre-training baseline* — and the paper reports a near-significant donepezil-vs-placebo pre-training threshold imbalance (rank p=0.05, driven by one subject at z=2.42; p=0.1 excluding him), which mechanically inflates a baseline-normalized gain. So the one existence proof does not show that the drug produced durability. It shows **training produced durability and the drug scaled its magnitude**, in n=8, on a baseline-normalized measure with a baseline imbalance.
- Shellshear 2015 (PMID 25900350) — levodopa + word learning, 1-month off-drug recognition advantage **for semantically-described items only**, recall null, effect size unobtainable (403), authors call it preliminary. That is the textbook multiplicity signature that dies on replication.
- Chamoun 2017 is not a third: n=9 at 4-14 month follow-up, donepezil +86% (p=0.043) vs control +45% (p=0.068) — a **difference of significances, not a significant difference**, with no between-group contrast reported.

**(b) Structurally implausible as framed ("a molecule that durably raises cognition").** Every **drug-alone** post-washout test in healthy people is null, including the largest and best-designed one: Rucker 2022 (PMID 35090363), n=89, psilocybin 25/10 mg vs placebo, CANTAB global composite and every domain **null at day 29**, with cognition as a *designed* endpoint. Meanwhile the same drug **reverses sign across plasticity assays in the same species and cortical area**: donepezil 5 mg *reduced* the monocular-deprivation ocular-dominance shift in healthy adults, t(11)=-4.9, p<0.001 (Sheynin 2019, PMID 30766471, CONFIRMED) — opposite to its effect on perceptual learning. A property that flips sign by readout is not a molecular property.

**(c) Not currently findable by prediction — this repo has already measured that.**

- `data\raw\persistence_positive_ledger.csv`: **19 verified** post-washout entries (mood 10, neuroplasticity 6, cognition 3). Of the 3 cognition entries, **0 are healthy** — Cerebrolysin (vascular dementia), Ibogaine (TBI), 7,8-DHF (AD model). The population classifier is deliberately biased *toward* calling entries healthy, and the cell is still empty.
- `data\raw\healthy_adult_cognition_ledger.csv`: 42 rows, and **no washout / retention / durability column exists** — the ledger that knows who enhances carries zero information about persistence.
- The two ledgers are **disjoint on the only axis that matters**. No join produces a durable-healthy-cognition label.
- Even the *acute* label is unpredictable: the mechanism-class prior that hit AUROC 1.00 on disease pivotal trials **collapses to 0.52 (perm p=0.45)** against healthy-adult ground truth. Methylphenidate and d-amphetamine are the same class with the same SMD (0.21) and opposite outcomes. The only separator is a coarse "is it a stimulant" gate (AUROC 0.83) that a pure statistical-power proxy (`n_studies`) beats at 0.88.
- Label budget already computed here: **~381 confirmed positive delayed-start readouts** before PERSEUS recall is estimable to ±0.1 at a 1% prior. And at that prior, with the engine's Jeffreys-upper FPR, **PPV = 0.03**.

**The scale evidence closes the escape route.** Where the drug×training hypothesis *has* been tested at scale, it failed: FLAME's motor effect (n=118, measured **on drug** at day 90 — it never contained post-washout evidence) went to cOR 0.96 across n=5907, Cochrane motor SMD 0.03 at high GRADE; AFFINITY's 12-month timepoint, 6 months after cessation, cOR 0.93. Dopaminergic: DARS n=593, OR 0.78 *favouring placebo*. Glutamatergic, the purest test — d-cycloserine does nothing alone and acts only when paired with the learning session — 21 RCTs, n=1047: peak **d=0.25 at post-treatment, decaying to d=0.19 with the CI crossing zero at follow-up.**

That last number is the whole verdict in one line: **the best-powered test of drug×training in humans peaks at d≈0.25, which is *inside* this project's already-established acute ceiling of 0.1-0.35 (max clean healthy SMD 0.34, nicotine).** Pairing a plasticity drug with training has not bought a single point of effect-size headroom, and it has not bought durability.

---

## 2. THE REFRAME

"Which molecule durably enhances cognition?" is the wrong question — not because it is uninteresting, but because **it names a variable that does not exist.** Durability is not a property of the molecule. It is a property of a four-way tuple, and the molecule is the *weakest* term in it.

**Old target:** `f(molecule) → durable cognitive SMD in healthy adults`. Label count in the target cell: **0**. Not scarce — empty.

**New target, stated precisely:**

```
f(compound, plasticity_assay, experience_protocol, retention_interval) → r
where r = (post-washout retained gain | drug + training)
        / (post-washout retained gain | placebo + IDENTICAL training)
```

Four things change, each forced by a specific piece of evidence:

1. **The unit of prediction becomes (compound × assay), not compound.** Forced by Sheynin 2019: same drug, same species, same cortex, opposite sign. Any compound-level "window score" — including this repo's L4 window — currently has no well-defined truth value.
2. **The outcome becomes a retention *ratio*, not a level.** Forced by Gilleen 2014: modafinil + 10 training sessions bought learning **rate** (AUC 24.7 vs 12.3, ES 0.84) with placebo catching up by the last day, no group difference in retention decay (F=0.2, p=0.65), and **zero transfer**. A pipeline that scores the learning curve will manufacture a rate effect and read it as durability.
3. **Both arms must contain training; drug-alone is a known-null control, not an arm of interest.** Forced by Rucker 2022 (n=89 null) plus every washout/discontinuation row in `data\raw\persistence_ground_truth.csv`.
4. **`experience_protocol` is a first-class predicted variable, with timing.** Forced by the rodent canon: chondroitinase required reverse lid-suture (Pizzorusso 2006), valproate/butyrate required reverse lid-suture (Silingardi 2010), and dark exposure had to **precede** occlusion removal (He/Quinlan 2007). Three independent reopener classes, same architectural requirement, and one of them is *order-dependent*. Timing is mechanism, not logistics.

**The screenable sub-target the codebase can actually chase today** is one level down and has a non-empty label set: `g(compound, assay) → P(opens a plasticity window in that assay)`. Durability is then delegated to the experience protocol, where the rodent evidence says it lives. This is already the engine's architecture — `psychoplastogen.py` marks a window as "direction-NEUTRAL and durable ONLY if paired with experience," and `perseus.py` reserves `DEMONSTRATED_HEALTHY` as an explicitly empty verdict. **The reframe is not a pivot; it is admitting that the engine's own firewall was right and the ranking layer was asking the wrong question.**

---

## 3. THE FIVE REQUIREMENTS

| # | Requirement | Exists today | Missing | Gap type |
|---|---|---|---|---|
| **R1** | **A non-empty, population-resolved durability label**: post-washout cognitive outcome in *healthy* people | 19 verified persistence rows, **0 in the target cell**; 2 world studies (n=8; one conditional subgroup) | Any usable number of positives. The label budget for estimable recall is ~381; even 20 would change the situation | **DATA — and it is a gap in the world, not in curation.** No amount of compute or curation fills it. Only human trials do, and the only registered one (NCT07226141, valproate + patching) is in **children 8-17**, primary completion **Sept 2027** |
| **R2** | **An assay-indexed window label**: does compound X open a plasticity window *in assay A* | Rodent OD canon (chABC, VPA, butyrate, dark exposure); human patching OD-shift assay (Sheynin, Min 2023 with a real Glx biomarker); 6 preclinical spine-density rows | No assay-keyed table anywhere in the repo. L4 is an *unfitted structural rule* on one channel: serotonergic recall 0.88, and **0.00 on NMDA, GABA-neurosteroid, muscarinic, TrkB, neurogenic, entactogen** | **DATA (curation, cheap) + BIOLOGY (the assay-dependence is real, per Sheynin)** |
| **R3** | **A representation that can see the durability lever** | MAMMAL DTI (BindingDB pKd), Kp,uu / free-exposure model, L4 permeability window (Vargas 2023 intracellular-pool logic; FPR 1/31, recall 6/6 on serotonergic decoys) | Allosteric/PAM sensitivity, PNN/ECM remodeling, intracellular-pool quantitation. **Measured, not speculated:** AMPA-PAM AUROC 0.26, AMPA orthosteric 0.09, MMP9 0.42 — all perm-p > 0.7 | **METHOD.** The BindingDB-pKd head is allosterically blind, and the plasticity levers are allosteric or extracellular-matrix |
| **R4** | **The experience protocol as a data field** | Nothing. Zero occurrences of a training/experience/protocol variable outside three docstrings | A schema: task, training dose, timing vs Tmax, ordering, sleep/consolidation, retention interval — the axis on which the 2 positives differ from every null | **DATA/SPEC (cheap) + BIOLOGY (ordering is load-bearing)** |
| **R5** | **A decision rule that knows when a positive is real** | **This one is satisfied.** Abstain-by-default PERSEUS (0/14 over-claims), Jeffreys not Wald CIs, grouped LOMO, label-shift transport, `ledger_guard.py` enforcing the stated inclusion rule (it caught the single l-theanine row the headline hinged on), `power_analysis.py` | Nothing structural | **Neither — and this is why the verdict is trustworthy.** The project's best-built component's honest output on this question is *abstain* |

---

## 4. THE BINDING CONSTRAINT

**R1.** The durable-healthy-cognition label gates every other requirement, and it is **not addressable by this codebase.**

R2, R3 and R4 are all instrumental — they exist to produce a ranked candidate whose value is realized only through R1. R5 is already built, and its verdict is abstention. So the chain terminates on a label that only new human trials can create, and the field has run **zero** such trials since 2013.

The operational corollary, stated bluntly: **the project cannot become the thing that finds the agent. It can become the thing that makes the finding-attempt cheap, well-aimed, and falsifiable.** The highest-ranked requirement this codebase can actually move is **R2**, because (i) it is fillable by curation, (ii) until the window label is assay-indexed, the engine's only live durability channel has no truth value at all, and (iii) an assay-indexed window screen is the artifact you hand a trialist.

Do not confuse these. R1 gates the **claim**. R2 gates the **build**.

---

## 5. FALSIFIABLE BUILD PLAN

Ordered cheapest-and-most-decisive first. Next free script number is 124. All are pandas/RDKit, CPU, hours not weeks. Every step has a kill criterion I expect to fire on at least one of them.

### B1 — The paired-vs-unpaired contrast. Test the reframe itself before building on it.
`C:\Users\Pierce Lonergan\Documents\GitHub\MAMMAL_Cognitive_Enhancement_Drug_Repurposing\scripts\124_paired_experience_contrast.py` → `data\raw\paired_experience_ledger.csv`

Assemble every human study with a post-washout cognitive/functional endpoint and record `(paired_experience, effect_at_washout, n, design)`. Rows available *now* from the verified research — **paired:** Rokem&Silver (+, n=8), Shellshear (+conditional), Chamoun (+preliminary, n=9), Walker-Batson (+, n=21, fails multiplicity), SSRI-amblyopia meta (+0.09 logMAR, subclinical), Gilleen (0 level), DCS meta (0.19, CI crosses 0, n=1047), FOCUS/AFFINITY/EFFECTS (0, n=5907), DARS (0/negative, n=593), levodopa+patching (0, n=139), vortioxetine+26wk training (0 at endpoint, n=100). **Unpaired:** Rucker psilocybin (0, n=89), esketamine (0), donepezil/galantamine washout (0), MPH/guanfacine discontinuation (0).

- **Success (pre-registered):** the paired-vs-unpaired difference in post-washout effect is significant by n-weighted permutation **and** survives dropping every study with n<25.
- **Kill:** not significant, or carried entirely by n<25 studies → the drug×experience hypothesis has no support at the level of evidence the project can assemble; the deliverable becomes the negative + the assay screen, and no further predictor work is licensed.
- **Honest prediction: the kill criterion probably fires.** Every paired study with n>100 is null-to-subclinical; every paired positive is n≤21 (Shellshear's n is unknown). Run it anyway — that is the point of a kill criterion, and this is the cheapest way to find out whether the reframe survives its own evidence.

### B2 — Make the window label assay-indexed, or delete it as a ranking signal.
`scripts\125_window_assay_index.py` → `data\raw\plasticity_window_assays.csv`; adds a required `assay` key to `src\mammal_repurposing\engine\psychoplastogen.py`; re-runs the decoy scan (`scripts\117_window_decoy_scan.py`) per assay family.

- **Success:** ≥40 (compound, assay) rows across ≥3 assay families, ≥8 compounds appearing in ≥2 families, and L4 agreement with the label beats the size-matched permutation gate (the same gate `scripts\114_ampa_pnn_channels.py` uses, perm-p < 0.05) in at least one family.
- **Kill:** among the multi-assay compounds, the window label's sign is inconsistent across families at ≥30% → "does compound X open a window" is not a well-defined object; **demote L4 from a ranking signal to a per-assay annotation** and stop scoring compounds on it. Also kill if L4 fails the permutation gate in *every* family.
- Decisive because Sheynin 2019 is either an outlier or the rule, and the entire compound-level screen hinges on which.

### B3 — Enter the two precedents honestly, then check the engine does not break.
Extend `EVIDENCE_RANK` in `src\mammal_repurposing\validation\persistence.py` with `paired_training_washout_retest` (rank 5, below `randomized_discontinuation`=6, with the within-subject/carryover risk documented). Add Rokem&Silver 2013 and Shellshear 2015 to `data\raw\persistence_positive_ledger.csv` flagged `UNREPLICATED`, Chamoun 2017 as `preliminary`. Add the decisive negatives to `data\raw\persistence_ground_truth.csv`: Rucker 2022 (drug-alone, null day 29), Gilleen 2014 (rate-not-level), Sheynin 2019 (sign flip), Levi/Chung 2020 **labelled as a non-replication by historical comparison, not a placebo-controlled null**. Extend `ledger_guard.py` to require, on any healthy+cognition+durable row: retention interval, off-drug flag, paired-experience flag, replication count.

- **Success:** PERSEUS still emits `DEMONSTRATED_HEALTHY` for **zero** compounds, over-claim rate stays 0, and recall moves no further than its existing Jeffreys CI [0.27, 0.73].
- **Kill:** adding 2 unreplicated positives flips any compound to `DEMONSTRATED_HEALTHY`, or moves recall outside [0.27, 0.73] → the persistence head is label-fragile at n≈2; freeze it and ship the ledger alone.
- Tests already exist: `tests\test_persistence.py`, `tests\test_perseus.py`, `tests\test_durability_gap.py`, `tests\test_ledger_guard.py`.

### B4 — The N* gate: compute whether the replication is even reachable.
Extend `src\mammal_repurposing\diagnostics\power_analysis.py`; report via `scripts\126_replication_power.py`.

From Rokem & Silver's reported statistics (SEMs 4.6/6.9 at n=8 → SDs ≈13.0/19.5, difference 12.9, between-group d ≈ 0.78), the required n per arm for 80% power at α=.05 two-sided is approximately: **d=0.78 → 26/arm; 50%-shrunk d=0.39 → 103/arm; at the d-cycloserine meta-analytic ceiling d=0.25 → 251/arm.** The design must be **between-group**, not the cheaper within-subject crossover — the one drug×training crossover in the literature (Gervain arm 2) died of carryover, and its blind was penetrated 17/18.

- **Success:** N* ≤ 60/arm under the shrunk effect → write a trial-design memo; the engine's contribution is candidate selection, and that is a real deliverable.
- **Kill:** N* > 250/arm → the human replication is out of reach for any single-site academic design; drop "find the agent" and ship "find the assay" plus the negative result.
- Note where those numbers land: **if the true effect is the meta-analytic augmentation ceiling, the study is ~500 people to durably improve performance on one motion-discrimination task, with no transfer.** That is the honest price, and it should be on the record before anyone commits.

### B5 — Census the world's pipeline, and lock forward predictions on whatever exists.
`scripts\127_healthy_durable_census.py`, reusing `src\mammal_repurposing\fetchers\clinicaltrials.py` (CT.gov v2) and `src\mammal_repurposing\reporting\trial_watch.py`.

Query for healthy-volunteer trials with a cognitive endpoint **after end of dosing**; classify each by paired-training arm and post-washout retest; emit a locked, time-stamped PERSEUS prediction (`WINDOW_CONDITIONAL` vs `NULL_SYMPTOMATIC`) per trial before readout.

- **Success:** ≥10 eligible trials → the field is about to generate the missing label, and the project's job is to have falsifiable predictions on record first. That is the forward test no retrospective AUROC can be.
- **Kill:** <3 eligible trials worldwide → the label will not arrive on any timescale worth waiting for. Publish the emptiness (`reports\pipeline\durability_gap_v1.md` is already 80% of that paper) and stop building predictors for an empty cell.
- **Expected:** near-zero. Nobody has run this design since 2013; NCT07226141 is pediatric and reads out Sept 2027; the human 10-day dark-exposure protocol (NCT02685423, n=8) has been active-not-recruiting **since 2016** with no results.

---

## 6. WHAT WOULD MAKE ME WRONG

Specific, checkable observations — in descending order of how much they would move the verdict.

1. **Transfer.** Any drug×training result in healthy adults where the retained off-drug gain **generalizes off the trained task.** Rokem & Silver is direction- and location-specific; Gilleen found zero transfer to MCCB/Cogstate. Without transfer, "durable cognitive enhancement" means "one trained perceptual skill retained," which is a far weaker claim than the program is built around. A demonstrated transfer result is the single largest possible update, and nothing in the current literature comes close.
2. **A properly powered between-group replication.** Pre-registered, placebo-controlled, **between-group**, n≥30/arm (per B4), healthy adults, drug + training vs placebo + **identical** training, retested ≥1 month off drug, with the retained **level** as the primary endpoint — not the learning curve — and the drug arm winning. That converts the verdict from "structurally thin" to "demonstrated for that task."
3. **An assay-invariant window.** If B2's multi-assay compounds show consistent sign across assay families, then the window *is* a molecular property, Sheynin 2019 is an outlier, and a compound-level screen is well-defined. That reopens "which molecule" as a legitimate prediction target and makes R2 tractable by curation alone.
4. **Knecht 2004's full text.** If it contains a **≥1-month off-drug** retention test with a real effect size in n=40 randomized parallel groups, the healthy precedent goes from two fragile studies to one properly designed one, and the verdict softens materially. **This is the cheapest possible falsification of my verdict — it is a paywall, not an experiment.** Get the PDF.
5. **A second mechanism class entering the cell.** A human HDAC or PNN/ECM window-opener with a post-washout cognitive readout in healthy adults. One class is an anecdote; two mechanistically independent classes is a mechanism. (Note the current floor: PNN degradation has **never** been done in a human brain — the one approved chondroitinase, condoliase, is intradiscal only, and chABC is a thermally unstable bacterial enzyme requiring local sustained delivery. The approved-drug-exists framing is a false comfort.)

**What would *not* make me wrong**, and this matters because the repo is full of it: more acute effects (the ceiling is measured at 0.34); more durable results in **patient** populations (restoration of a degraded system is a different mechanism — there is a deficit to reverse); more preclinical spine-density persistence (6 of the 19 verified rows already are); more durable **mood** effects (10 of 19). None of those touch the cell. The cell needs healthy **and** cognitive **and** off-drug, and it needs the paired-experience arm.

---

### Bottom line for the record

Durable cognitive enhancement in healthy adults is **not demonstrated**, is **plausible only for narrow trained skills**, and is **not currently findable by prediction** — because the label cell is empty from both directions, because the mechanism-class predictor that worked on disease collapses to chance here (0.52), and because the best-powered human test of the exact hypothesis peaks at d≈0.25, inside the acute ceiling it was supposed to exceed. The single existence proof (n=8) shows training producing durability with a drug scaling its baseline-normalized magnitude — not a drug producing durability.

The right move is not a better model. It is to (i) test the reframe against the project's own assembled evidence (B1), (ii) make the one live durability channel well-defined or delete it (B2), (iii) put the honest price of the replication on the record (B4), and (iv) if the world is not generating the label, **publish the emptiness as the finding** (B5). `reports\pipeline\durability_gap_v1.md` already says this. The build plan above is how to make it falsifiable rather than merely correct.