# Could this pipeline ever have discovered what it already knows?

Before building anything that predicts NOVEL compounds, two prior questions bound what such a thing could do here. Both are answerable from data already in the repository.

## 1. How much of the evidence base is even a designable molecule?

The primary analysis set is **n = 21**. Of those, **12** are not a single small molecule at all and are therefore permanently outside the reach of de novo design, whatever the model:

| compound | why it cannot be designed |
| --- | --- |
| guarana | botanical extract; many constituents, caffeine among them |
| ginkgo_biloba | botanical extract (EGb 761 and relatives); flavonoid and terpene mixture |
| bacopa_monnieri | botanical extract; bacoside mixture |
| multivitamin_mineral | mixture of many vitamins and minerals; no single structure exists |
| b_vitamins | a class of several distinct molecules dosed together |
| folic_acid | single molecule, but a vitamin repletion exposure rather than a design target |
| menopausal_hormone_therapy | several regimens of different hormones |
| psilocybin_lsd_microdosing | two distinct molecules pooled as one exposure |
| dietary_nitrate | an ion delivered in food, not a designed molecule |
| insulin_intranasal | 51-residue protein; outside small-molecule generative scope |
| oxytocin_intranasal | nine-residue cyclic peptide; outside small-molecule generative scope |
| fruit_derived_polyphenols | a polyphenol CLASS, not a compound |

**The enhancer class is the part that matters**, because a generator would be conditioned on the positives. There are **6** labelled enhancers, of which **4** are single molecules with a resolved structure: methylphenidate, modafinil, caffeine, nicotine.

Conditioning a generative model on 4 molecules is not discovery. Whatever it produced would be analogues of those few, and they are already known, already scheduled or already consumed at scale.

## 2. The novelty ceiling

`scripts/95_novel_compound_onboarding.py` routes a novel SMILES to a known mechanism class only when its maximum Tanimoto to a known chemotype clears **0.35**, and abstains otherwise, on the stated grounds that it cannot invent a mechanism. On its own demo, six of eight novel compounds abstained.

So: hold each known compound out and present it as if it were novel. Does the engine clear its own threshold against the remaining actives, or does it refuse to score a compound we already know works?

### Against the FULL exemplar base (110 compounds, what the engine actually does)

| compound | label | max Tanimoto | nearest exemplar | verdict |
| --- | :---: | ---: | --- | --- |
| modafinil | 1 | 1.000 | armodafinil | routed |
| dextroamphetamine | 0 | 0.400 | lisdexamfetamine | routed |
| l_theanine | 0 | 0.379 | D-serine | routed |
| testosterone | 0 | 0.288 | ganaxolone | **ABSTAIN** |
| cannabidiol | 0 | 0.274 | nabilone | **ABSTAIN** |
| nicotine | 1 | 0.264 | Org 26576 | **ABSTAIN** |
| creatine | 0 | 0.259 | sarcosine | **ABSTAIN** |
| methylphenidate | 1 | 0.240 | modafinil | **ABSTAIN** |
| caffeine | 1 | 0.150 | ABT-288 | **ABSTAIN** |

**6 of 9** would be abstained on against the full base. Of the **3** that route, **1** match at or above 0.95, which is the compound finding ITSELF under another name: modafinil to armodafinil at 1.000.

Read the routed rows carefully, because two of them are not discoveries. Modafinil matches armodafinil, which is its own single enantiomer, at a Tanimoto of exactly 1.000 since the fingerprint ignores chirality. Dextroamphetamine matches lisdexamfetamine, which is its own lysine prodrug. Subtracting those leaves ONE non-trivial route in the whole set, and it barely clears the threshold.

### Against the healthy-adult evidence base only

The stricter and more relevant comparison, because it is the healthy-adult evidence that defines this project's outcome. A compound routed on its resemblance to a CNS drug with no healthy-adult cognitive evidence has been matched to a chemotype, not to a demonstrated effect.

| compound | label | max Tanimoto to any other | nearest | verdict |
| --- | :---: | ---: | --- | --- |
| dextroamphetamine | 0 | 0.250 | modafinil | **ABSTAIN** |
| modafinil | 1 | 0.250 | dextroamphetamine | **ABSTAIN** |
| methylphenidate | 1 | 0.240 | modafinil | **ABSTAIN** |
| creatine | 0 | 0.184 | l_theanine | **ABSTAIN** |
| l_theanine | 0 | 0.184 | creatine | **ABSTAIN** |
| nicotine | 1 | 0.179 | methylphenidate | **ABSTAIN** |
| testosterone | 0 | 0.154 | methylphenidate | **ABSTAIN** |
| cannabidiol | 0 | 0.127 | testosterone | **ABSTAIN** |
| caffeine | 1 | 0.122 | nicotine | **ABSTAIN** |

**9 of 9** structurally resolved compounds would be abstained on, including **4 of 4** known ENHANCERS.

Median maximum Tanimoto among the enhancers is **0.209**.

## What this bounds

The similarity gate is not a tuning parameter that could be relaxed. It is what keeps the engine honest: below it, the engine has no evidence that the compound belongs to any precedented mechanism, and routing it anyway would be inventing a mechanism. Lowering the threshold does not create knowledge, it only stops recording ignorance.

So a similarity-gated screen built on this evidence base can only find compounds that resemble the handful it already has. That is interpolation. It is a legitimate and useful thing to do, and it should be described as analogue search rather than discovery.

The deeper limit is the one in section 1, and no model fixes it: most of what this project has established as working in healthy adults is not a molecule anyone could design.

Generated by `scripts/132_novelty_ceiling.py`.