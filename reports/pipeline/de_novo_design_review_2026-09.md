# Predicting novel compounds: what would the prediction mean?

The pipeline ranks known compounds. The natural next ask is to predict novel ones. This reviews what
that would take, and concludes that the useful move is not the obvious one.

Six lanes, 48 items, 14 adversarial verifications completed (13 confirmed, 1 misdescribed). Unlike
the previous survey, the verification phase ran: the workflow defect that reported empty shells has
been fixed, and the counts here are real.

## The structural fact about every generative method

**A generative model is an amplifier of a scoring function. None of them contain biology.**

REINVENT4 maximises a user-supplied score. Pocket-conditioned diffusion models (Pocket2Mol,
TargetDiff, DiffSBDD, ResGen) maximise pocket-conditioned atom likelihood and are then filtered by
docking. GFlowNets sample proportional to a reward. Active-learning loops select against an oracle.

So the generative question is never "which architecture". It is **"what is the objective, and is it
measured?"**

This project's objectives are measured:

| target | DTI head AUROC | permutation p |
| --- | ---: | ---: |
| AMPA PAM | **0.26** | 0.945 |
| AMPA orthosteric | **0.09** | 0.992 |
| MMP9 | 0.42 | 0.717 |

Source: `reports/pipeline/ampa_pnn_channels_v1.md`, from `scripts/114_ampa_pnn_channels.py`.

CORRECTED, same day, by `reports/pipeline/dti_scale_lohi_v1.md`. This section originally read that
two of those numbers are below 0.5, that this is not "weak" but **anti-correlated**, and that a
generator pointed at an anti-correlated scorer would produce confidently anti-active molecules. Two
things were wrong with that.

First the n. Those AUROCs were measured on three and four anchor compounds, which cannot resolve
anything. Second, and more importantly, the anti-correlation claim does not survive measurement.
Re-running the head at ChEMBL scale across twelve cognition targets, 120 actives against 120 hard
negatives each, gives a pooled mean AUROC of **0.468** and **zero of twelve** targets clearing the
project's gate. Five targets do sit significantly below chance, including GRIA1 at 0.387 with an
interval of [0.319, 0.459]. But a paired within-compound test, comparing each molecule's score at a
target it binds against its score at a target it does not, shows the apparent reversal is an
artifact of target-level score offsets: it is significant raw and vanishes under within-target
centering (48 of 100 higher at the real target, sign test p = 0.76).

So the scorer is **uninformative, not inverted**. The conclusion for generation is unchanged and the
reasoning is stronger, because an uninformative objective bears no relationship to activity in
either direction. There is no sign to flip and nothing to exploit: optimising it hard produces
molecules whose activity is simply unconstrained.

## A generated molecule can never acquire this project's evidence type

The inclusion rule is a clean healthy-adult meta-analysis whose interval lies entirely above zero. A
newly designed molecule has no trials, so it can never satisfy that rule, and no generated prediction
can ever be scored right or wrong.

That converts a falsifiable repurposing project into an unfalsifiable design project. The concern is
not aesthetic: the project's whole claim to rigour is that its statements have kill criteria.

## The famous success is the pathology this project already measured

GENTRL/DDR1 (DOI 10.1038/s41587-019-0224-x) is the canonical AI-design result: 6 compounds
synthesised, 2 potent, 46 days. Walters and Murcko (DOI 10.1038/s41587-020-0418-2) showed that the
best compound's behaviour followed from its structural similarity to ponatinib, an already-known DDR1
inhibitor. The authors' reply conceded that in-depth validation was not the goal.

That is the same failure this project measured in its own engine this week: 0.972 class recovery
achieved on near-neighbours, and on the novelty-ceiling test modafinil "discovering" armodafinil at a
Tanimoto of exactly 1.000, its own enantiomer.

Two further cautions from the benchmark literature, both verified:

- **PoseCheck** (arXiv:2308.07413) redocked five leading structure-based generators and found their
  poses physically worse than the reference ligands: median strain energies of 194.9, 592.2, 1241.7,
  1243.1 and 18,693.8 kcal/mol against a 102.5 reference, with median redocking RMSD above 2.0 A for
  every method.
- A 16-model benchmark (arXiv:2406.03403) found **AutoGrow4**, a 2D genetic algorithm that never sees
  the protein in 3D, dominates the 3D methods on optimisation ability.
- The naming convention is itself misleading: TacoGFN reports a "52.63% Novel Hit Rate" where the
  "hit" is a docking threshold, not an assay.

## What survives, and it inverts the direction

Three experiments answer the generation question by measurement, before a line of generator code is
written, and two are likely to answer it NO. All are zero cost.

**1. The Renz control-model exploitation test.** Point a generator at this project's own ranker as
the reward, then re-score the optimised molecules with independently trained control models on
disjoint splits (DOI 10.1016/j.ddtec.2020.09.003, PMID 33386095; code at
github.com/ml-jku/mgenerators-failure-modes). If the score collapses under the control, the ranker is
exploitable, and every downstream "novel enhancer" is an artifact of the scorer rather than a
property of chemistry. A documented instance of exactly this: a REINVENT4 run with a docking-only
reward grew long flexible aliphatic linkers to manufacture protein contacts and inflate the score
while destroying ligand efficiency (DOI 10.1021/acs.jcim.5c01203, PMID 40893044).

**2. AddCarbon as the null model.** Insert one carbon at random into a training molecule. By
construction it emits only valid, "novel", unique molecules and scores extremely well on standard
novelty metrics (Renz et al.). Any generator this project builds must beat AddCarbon by a
pre-registered margin on the identical metric suite. Without that control, every novelty number is
uninterpretable. This is the doc-290 lesson (rank-calibrate the baseline) applied to generation.

**3. A Lo-Hi split of the existing DTI head.** DONE, see `reports/pipeline/dti_scale_lohi_v1.md`.
It answered no, and more cheaply than expected. A single Tanimoto boundary turned out not to
separate familiar from unfamiliar in ChEMBL at all, so the split was run as a similarity gradient
instead. Mean AUROC across quartiles of neighbourhood density rises monotonically, 0.40 to 0.46 to
0.50 to 0.51, which confirms a near-neighbour effect and simultaneously makes it academic: the
effect tops out AT chance. Distance from a target's earliest ligands does nothing (0.47, 0.46, 0.48,
0.47).

Note the relationship between two numbers. The onboarding engine's abstention threshold is 0.35. The
Lo-Hi splitter marks **0.4** as the boundary of the hard regime. The engine therefore operates
entirely inside the easy regime by construction, which is another way of saying what the novelty
ceiling measured directly.

## The alternative that actually fixes the binding constraint

The adversarial lane's recommendation, and it is the best idea in this survey:

> Freeze the current ranker and **pre-register dated, prospective ranked predictions** against
> healthy-adult cognition trials and meta-analyses that will read out anyway.

The project cannot generate new human data. But the world is generating it continuously, and a dated
prediction registered before a readout is a scoreable claim that needs no wet lab, no lab budget, and
no new evidence type. It is the one route that produces a falsifiable statement about compounds this
project has never assessed.

That is the honest answer to "predict novel compounds": this project cannot defensibly predict novel
**molecules**, but it can make novel **predictions**, and those are what a reader would actually
credit.

## Sequencing

1. ~~Run the Lo-Hi re-split of the DTI head.~~ DONE 2026-09-20. Decisive: 0 of 12 targets rank,
   pooled mean AUROC 0.468.
2. Run the Renz exploitation test against the project's own ranker. A collapse is publishable.
3. Implement AddCarbon as the mandatory null before any generative metric is reported.
4. Pre-register dated prospective predictions against readouts already scheduled.
5. Do not build a generator to find enhancers.

## One item failed verification

The time-split and scaffold-split recommendation came back REJECT_MISDESCRIBED.

CORRECTION, same day. An earlier version of this section stated that all five of that item's
identifiers resolved correctly, that its characterisation of Sheridan 2013 was exact, and that the
rejection fell on a sub-claim rather than the core. Those specifics were written from a workflow
summary that is no longer retrievable, and I could not reproduce the verifier's verdict text from
the session record when I went back to check. They are withdrawn rather than left standing on a
source I cannot produce.

What IS on the record, recovered verbatim from the session log, is the item's own text, and it is
accurate on every point I can check independently: Sheridan 2013 (J Chem Inf Model 53(4):783-790,
DOI 10.1021/ci400084k, PMID 23521722) is correctly described as a retrospective comparison against
Merck's actual prospective predictions rather than a wet-lab test; MoleculeNet (DOI
10.1039/C7SC02664A, PMID 29629118), Bemis-Murcko (DOI 10.1021/jm9602928), SIMPD (DOI
10.1186/s13321-023-00787-9) and LIT-PCBA (DOI 10.1021/acs.jcim.0c00155) are correctly attributed.

So the split methodology stands on its own sources and remains correct practice. What I cannot tell
you is which part of that item the verifier objected to. Treat the item as unverified rather than as
either confirmed or discredited.

Written by hand; this report has no generator. It summarises a research sweep whose structured
output is not committed.
