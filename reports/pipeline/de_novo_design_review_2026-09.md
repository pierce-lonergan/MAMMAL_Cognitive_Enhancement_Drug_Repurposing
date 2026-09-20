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

Two of those are below 0.5, which is not "weak". It is **anti-correlated**. A generator pointed at an
anti-correlated scorer produces confidently anti-active molecules, faster and in greater diversity
than any human could. That is the whole argument in one line, and it is measured rather than feared.

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

**3. A Lo-Hi split of the existing DTI head.** Re-split its ChEMBL data at ECFP4 Tanimoto below 0.4
and measure whether it predicts anything at all outside near-neighbour range. Given it already scores
0.26 and 0.09 at those sites, the expected answer is no, and that is a cheap, decisive end to
the generation question.

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

1. Run the Lo-Hi re-split of the DTI head. Cheapest, and likely decisive.
2. Run the Renz exploitation test against the project's own ranker. A collapse is publishable.
3. Implement AddCarbon as the mandatory null before any generative metric is reported.
4. Pre-register dated prospective predictions against readouts already scheduled.
5. Do not build a generator to find enhancers.

## One item failed verification

The time-split and scaffold-split recommendation came back REJECT_MISDESCRIBED. All five of its
identifiers resolved correctly and its characterisation of Sheridan 2013 was exact; the rejection was
on a sub-claim rather than the core, and the split methodology itself stands. Recorded here so the
rejection is not read as discrediting split-based validation, which remains correct practice.

Written by hand; this report has no generator. It summarises a research sweep whose structured
output is not committed.
