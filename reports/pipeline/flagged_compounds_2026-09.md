# What this system actually flags, and how much of each flag to believe

Written by hand on 2026-09-20 against the repository's own data. Every number below is reproducible
from `data/raw/healthy_adult_cognition_ledger.csv`, the archive predicates, or the named report.

The short answer: **the curated evidence flags four independent molecules with verified ACUTE
effects, and zero compounds with verified durable effects.** The structural engine flags thirty-one
more, and its output should not be acted on; the last section shows why, using its own rows.

## Tier 1: verified acute enhancers

Clean meta-analysis, genuinely healthy population, pooled interval entirely above zero, and a point
estimate reaching the project's meaningful-effect floor of g = 0.20.

| compound | g | 95% CI | k | population | source |
|---|---:|---:|---:|---|---|
| l-theanine | 0.35 | [0.10, 0.61] | 4 | healthy | Payne 2025 Nutrition Reviews, PMID 40314930 |
| nicotine | 0.34 | [0.18, 0.50] | 9 | non-smokers | Heishman 2010 Psychopharmacology, PMID 20414766 |
| caffeine | 0.28 | [0.21, 0.36] | 31 | rested healthy | Klove & Petersen 2025, PMID 40335666 |
| methylphenidate | 0.21 | [0.09, 0.32] | 24 | healthy | Roberts 2020 Eur Neuropsychopharmacol, PMID 32709551 |

Two more pass the same filter but are NOT independent, because both contain caffeine and caffeine is
already labelled: caffeine plus taurine (0.53 [0.03, 1.07]) and theanine plus caffeine (0.33 [0.13,
0.54]). The ledger flags both in `red_flags` rather than quietly counting them. The taurine row has
a second defect recorded against it: its interval is a 95% credible interval, not a confidence
interval, so it is a type mismatch against the stated rule.

Caveats the ledger carries on the four, in its own words:

- **caffeine** is the most replicated at k = 31, and its exact CI is paywalled; the point estimate
  is verified, the interval is not. It is also the only compound whose interval sits entirely at or
  above the 0.20 floor.
- **methylphenidate** has the best single domain at recall 0.43 [0.20, 0.65] with I-squared 0, and
  Egger is non-significant overall.
- **nicotine** holds in non-smokers specifically, with all six domains showing intervals above zero
  (0.16 to 0.44). This is not withdrawal relief.
- **l-theanine** is borderline and marked so: only choice reaction time is significant, and simple
  reaction time is null. It rests on k = 4.

## Tier 2: passes on magnitude, fails on evidence quality

Interval above zero and g at or above 0.20, but the pooled estimate is not a clean healthy-adult
meta-analysis. Listed because they are the honest next candidates, not because they are established.

| compound | g | 95% CI | k | why it does not reach tier 1 |
|---|---:|---:|---:|---|
| mixed polyphenols | 1.47 | [0.62, 2.32] | 4 | contested tier, k = 4, implausibly large |
| ashwagandha | 0.52 | [0.27, 0.78] | 6 | contested tier; population is mixed/stressed, not healthy |
| anthocyanins | 0.46 | [0.30, 0.63] | 11 | mixed population |
| resveratrol | 0.39 | [0.08, 0.70] | 3 | mixed population, k = 3 |
| exogenous ketones | 0.36 | [0.18, 0.54] | n/a | mixed population |
| oxytocin | 0.29 | [0.15, 0.43] | 7 | contested tier |
| cocoa flavanols | 0.22 | [0.01, 0.43] | 11 | mixed population; lower bound touches zero |

## Tier 3: real but trivial

Interval excludes zero, so the effect is probably real, and the point estimate is below the 0.20
floor, so it is not of practical size. Reporting these as enhancements would be technically true and
substantively misleading.

flavonoids total 0.148 [0.098, 0.198] (k = 80; quoted to 3 decimals because rounding its upper bound to 0.20 would make it look as though it reaches the floor) - intranasal oxytocin 0.13 [0.02, 0.24] (k = 23) -
**modafinil 0.12 [0.02, 0.21]** (k = 14) - multivitamin/mineral 0.07 [0.03, 0.11] (k = 3)

Modafinil belongs in this tier and that is worth stating plainly, because it is one of the two
best-known "smart drugs". Its ledger note records the pooled effect as TOST-equivalent to zero,
trivially small.

## Tier 0: durable enhancement, which is the actual goal

**Nothing.** The machine-checkable predicate returns:

> `durable_healthy_rows(1)` -> False. "0 row(s) meeting the full durability standard vs required 1.
> 5 row(s) are off-drug and healthy at all; 1 carry an unqualified positive direction; rejected on
> the standard: donepezil(unreplicated)."

Every compound in every tier above is an ACUTE, on-drug effect. Take the drug, perform better while
it is in you, return to baseline. Across the whole corpus, five rows are even off-drug and healthy,
one of those points positive, and it is a single unreplicated donepezil result. That is G1, and it
is why this project exists.

## The impairers, which are also real results

The best-evidenced findings in the entire ledger are negative, and several are larger than any
enhancement.

psilocybin -1.13 [-1.70, -0.57] - scopolamine -0.86 [-1.08, -0.64] - melatonin -0.74 [-1.03, -0.45]
- psilocybin/LSD microdosing -0.34 [-0.62, -0.06] - acute alcohol -0.30 [-0.39, -0.21] -
dehydration -0.21 [-0.31, -0.11]

All six are clean meta-analyses. Note that microdosing has its own row and its own estimate: it is
not a small version of the full-dose effect and must not be conflated with it.

## What the structural engine flags, and why not to act on it

`reports/pipeline/f2_catalogue_shortlist_v1.md` screens 2,235 approved drugs and returns **31
repurposing hypotheses**. The report is honest that "the predicted g is the assigned class's prior,
a model output, not a measured outcome". Read the rows and it is worse than that caveat implies.

Every compound assigned to a class receives the IDENTICAL prediction. All seventeen
catecholaminergic hits get +0.40 [+0.35, +0.45] with P(success) 0.90. All twelve cholinesterase
hits get +0.20 [+0.12, +0.28] with P(success) 0.56. The compound contributes nothing to the number;
only its class assignment does, and that assignment is made by ECFP4 structural similarity.

That similarity step is the one measured this month, and it does not work:

- `dti_scale_lohi_v1.md`: 0 of 12 targets rank, pooled AUROC 0.468, five significantly below chance
  at n = 120 per side.
- `novelty_ceiling_v1.md`: against the healthy-adult evidence base, 9 of 9 compounds abstain and the
  maximum pairwise Tanimoto is 0.250.

The shortlist shows what that produces. Its `AChE_inhibitor` class contains **codeine, oxycodone,
hydrocodone, dihydrocodeine and fenoprofen** - four opioids and an NSAID, none of them
cholinesterase inhibitors - each carrying a predicted +0.20 and P(success) 0.56. Its
`wake_promoting` class contains **hydroxyzine**, a first-generation sedating antihistamine, at
predicted +0.36 and P(success) 0.79.

The hydroxyzine row is the sharpest case, because the project already holds the contradicting
evidence. Its own targeted positive hunt confirmed both antihistamine generations as cognitive
IMPAIRERS; they are absent from the ledger for a metric reason (weighted mean difference rather than
standardised mean difference), not because the impairment was in doubt. So the engine predicts a
79% chance of cognitive enhancement for a drug class this project has separately verified impairs
cognition.

**Nothing on that shortlist should be acted on.** It is a class-prior lookup wearing a
compound-level interface, and the routing step underneath it has been measured at chance.

## What has a readout coming

From `data/raw/prospective_healthy_adult.csv`, dated predictions frozen 2026-09-20 under a rule
committed beforehand:

- **caffeine**, five pending readouts: four trials plus PROSPERO review CRD420251163859. Predicted
  POSITIVE on all five.
- **l-theanine**, PROSPERO review CRD420251160876. Predicted POSITIVE.
- **nicotine**, NCT07408180, primary completion 2026-03-16. Predicted POSITIVE.
- **caffeine plus taurine**, PROSPERO review CRD420261499408, registered thirteen days ago.
- **methylphenidate**: no live prediction. The two pending reviews that cover it pool it with four
  or five other stimulants, so neither yields a single-compound estimate.
- **modafinil**: no live prediction anywhere. One candidate was overturned as an early-psychosis
  study, the other has a BOLD primary rather than a cognitive one.

## Candidates the system had never considered, surfaced 2026-09-20

The eight-lane sweep that produced the readout registry searched by INTERVENTION NAME, from a
hand-written compound list. A critic pointed out that such a design can only find what someone
already suspected. A second pass searched ClinicalTrials.gov from the OUTCOME side with no
intervention filter at all, and opened registries the first pass never reached.

It confirmed **77 further pending readouts**, more than the first pass found, with **zero**
identifiers failing to resolve. Across 81 aggregated rows not one pass-1 CONFIRMED identifier
reappeared, so the two passes produced near-disjoint sets. On ClinicalTrials.gov alone, where both
passes searched, the outcome-side query independently re-found 19 of 40 pass-1 NCTs, which is a
recall check on both. The honest conclusion is that the first pass captured at most about **45%**
of the reachable corpus, and that the failure was of two separate kinds: the compound-list framing,
and registries nobody opened.

None of the following is evidence of efficacy. Every one is an ONGOING study with no result. They
are listed because a compound-name search could not have generated them, which is the point.

**Mechanistically informative, in rough order of what they would settle:**

- **Paraxanthine**, the primary human metabolite of caffeine, on chronic 6-week dosing. Caffeine is
  this project's best-evidenced enhancer, and this trial separates "caffeine works" from "the
  adenosine-antagonist mechanism downstream of caffeine works". It abstains in the registry
  precisely because mapping it onto the caffeine row would destroy the only question it asks.
- **Equol**, the gut-microbial metabolite of daidzein, administered DIRECTLY. Only 10 to 20% of
  Western populations produce it, and that producer/non-producer split is what has confounded soy
  isoflavone trials for decades. Giving the metabolite bypasses the confound.
- **Reboxetine head-to-head against methylphenidate** under standardised mental fatigue: a
  selective noradrenaline reuptake inhibitor against the ledger's own labelled enhancer.
- **Pseudoephedrine head-to-head against modafinil**, funded by the Japan Anti-Doping Agency.
- **Levetiracetam**, an anti-seizure drug, in cognitively healthy mid-life APOE e4 carriers.
- **Daridorexant and suvorexant**, dual orexin receptor antagonists, on overnight memory
  consolidation, with suvorexant run against zolpidem to dissociate orexin from GABA-A.
- **GT-002**, a first-in-class GABA-A PARTIAL positive allosteric modulator, against oxazepam.
- **Creatine under physiological stress rather than at rest**: single 0.2 g/kg dose across 21 hours
  of sleep deprivation, and loading under acute sleep restriction. The ledger's creatine row is
  NULL at rest, so a stressor-conditioned effect would be new information rather than a
  contradiction.
- **Magnesium L-threonate** with a PRE-SPECIFIED sex-divergence hypothesis under a mental-fatigue
  challenge.
- **Caffeine mouth rinse**, tasted and expectorated without ingestion, alongside an espresso
  mouthwash arm. A route manipulation that isolates oral and trigeminal signalling from systemic
  pharmacokinetics, and the reason both abstain rather than mapping to the ingested caffeine row.

**Longevity and frontier agents now in healthy adults:** fisetin (senolytic, DSST endpoint),
alpha-ketoglutarate over 24 weeks, L-serine 6 g/day in 65 to 85 year olds, nicotinamide riboside,
rapamycin, tocotrienol-rich fraction, and CX-50, which is orally ingested bovine adipose
extracellular vesicles at 80 g/day scored on the NIH Toolbox. Also **AFA-281**, an undisclosed
orally available small molecule in Phase 1 with choice reaction time as a named endpoint, which no
compound list can contain because its structure is not public.

**And one that is not a compound at all:** oral faecal microbiota transplant from young, physically
active donors, in healthy adults. Alongside a fully-provided hazelnut-oil versus palm-oil diet and
a quadrivalent influenza vaccine used as an acute cognitive decrement model, these are the clearest
demonstrations of what outcome-side searching buys: the generator for the first pass was "things
already believed to be nootropics", which is precisely the set with no information value.

**What this does NOT change.** Not one of these has a result. None of them moves Tier 0, which is
still empty. They are hypotheses with dates attached, and the registry exists so that when they
read out, what this project predicted beforehand is already on the record.

## The honest bottom line

If the question is which compounds have potential to make you sharper this afternoon, the answer is
caffeine, and then nicotine, methylphenidate and l-theanine with smaller and less certain effects.

If the question is which compounds have potential to raise cognitive ability durably, past washout,
the answer this system returns is **none, and it has looked hard**. That is a measured result rather
than an absence of effort: the durability predicate is machine-checked on every run, the search for
new positives across seven mechanism classes returned zero enhancers and nine confirmed impairers,
and the structural route to novel candidates has been refuted and closed.

Written by hand; this report has no generator. Sources are named inline.
