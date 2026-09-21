# When can the forward registry detect anything?

140 verified pending readouts, **47 dated calls**. Before anything is
built on this, the prior question: at this size and balance, what could it distinguish
from a constant predictor at all?

## The sample size is not the number of calls

The rule and the constant predictor are scored on the SAME rows, so the comparison is
paired and McNemar is the right test. That has a sharp consequence. On every row where the
rule itself says NULL_EFFECT, the two predictors agree, the pair is CONCORDANT, and it
contributes nothing to separating them no matter how the trial turns out.

- dated calls: **47**
- concordant with the constant predictor (NULL_EFFECT): **16** - these can never contribute
- **discordant, and therefore the real sample size: 31** (25 POSITIVE, 6 NEGATIVE)

This is the same arithmetic that emptied the patient arm's track record, seen in advance
rather than after the fact.

## What the registry can detect at full resolution

With all 31 discordant pairs resolved, a one-sided McNemar test at alpha 0.05 needs **21 wins out of 31** (68%) to reject the null that the rule is no better than the constant.

| true per-pair win rate | power at n = 31 | n needed for 80% power |
|---|---:|---:|
| 0.55 | 11% | 620 |
| 0.60 | 25% | 158 |
| 0.65 | 46% | 69 |
| 0.70 | 69% | 37 |
| 0.75 | 87% | 23 |
| 0.80 | 97% | 18 |
| 0.90 | 100% | 8 |

## The timeline, which is the surprise

The registry is not a 2028 problem. Cumulative discordant pairs whose stated readout date
falls at or before each year:

| by end of | discordant pairs due | power at true 0.70 | power at true 0.80 |
|---|---:|---:|---:|
| 2024 | 2 | 0% | 0% |
| 2025 | 10 | 15% | 38% |
| 2026 **(today)** | 26 | 63% | 94% |
| 2027 | 27 | 58% | 93% |
| 2028 | 27 | 58% | 93% |
| 2029 | 27 | 58% | 93% |
| 2030 | 27 | 58% | 93% |

**26 of the 31 discordant pairs already have a stated readout date at or before the end of 2026.** The rate-limiting step is therefore NOT waiting for trials to finish. It is finding out whether the results have been published, which is work that can be done now.

## Reading

**The registry is one fifth the size it looks.** 140 readouts, 47 calls, and 31 rows that
can actually separate the rule from a constant predictor. Any future headline must be
quoted against 31, not 140.

**It is already decisive IF the rule is good.** 26 of the 31 discordant pairs have a stated
readout date at or before the end of this year. At n = 26 the power to detect a true
per-pair win rate of 0.80 is 94%. So if the healthy-adult ledger genuinely predicts new
readouts well, that is measurable NOW and does not need a single new trial. The bottleneck
is not the calendar, it is finding out which of those 26 have published.

**It is hopeless if the rule is only slightly good.** At a true win rate of 0.60, detecting
it at 80% power needs 158 discordant pairs, five times what two search passes produced. A
modest edge is not reachable by searching harder; it is out of range of this design.

**So the targeting criterion changes, and this is the actionable part.** Scaling the
registry does NOT mean adding readouts. A readout on which the rule says NULL_EFFECT is
concordant with the constant predictor and contributes exactly nothing, and 16 of
the current 47 calls are in that position. Growth has to come from readouts on
compounds whose ledger interval lies clearly ABOVE or BELOW zero, because only those
generate discordant pairs. Searching for more trials on compounds whose ledger row spans
zero is wasted effort, however many it finds.

## How much more searching, and is it worth it

Measured end-to-end yield over both search passes: 140 confirmed readouts produced 47 calls and 31 discordant pairs, so 22% of a readout's worth survives to the test. Of readouts that map to a ledger compound at all (51 of 140, 36%), 61% are discordant.

| to reach | needed for | more discordant pairs | untargeted readouts | targeted |
|---|---|---:|---:|---:|
| n = 18 | 80% power at true 0.80 | 0 | none, already there | none |
| n = 37 | 80% power at true 0.70 | 6 | ~27 | ~16 |
| n = 158 | 80% power at true 0.60 | 127 | ~574 | ~349 |

So the decision is tiered, and which tier applies is not yet known:

- **If the rule is strong (true win rate around 0.80), the search is already finished.** 31 discordant pairs exceed the 18 needed, and 26 of them are already due. Nothing is gained by searching more; everything is gained by resolving what is held.
- **If it is moderate (0.70), one targeted pass closes the gap.** Six more discordant pairs, which is roughly 16 to 27 more confirmed readouts.
- **If it is weak (0.60), the design cannot reach it.** 127 more discordant pairs means roughly 350 to 575 more confirmed readouts, two and a half to four times the total output of two full search passes, for an edge small enough that it would not change any advice.

**Which means resolution must come before more searching.** Until the held rows are scored nobody knows which tier this is, and three of the four possible answers make further searching either unnecessary or pointless.

### If a targeted pass is run, these are the only compounds worth searching

Exactly 19 ledger compounds generate a discordant call.

- expect a benefit (13): anthocyanins, ashwagandha, caffeine, caffeine_plus_taurine, cocoa_flavanols, exogenous_ketones, l_theanine, methylphenidate, nicotine, oxytocin, polyphenol_mixed, resveratrol, theanine_plus_caffeine
- expect a decrement (6): alcohol_acute, dehydration, melatonin, psilocybin, psilocybin_lsd_microdosing, scopolamine

Every other compound in the ledger yields NULL_EFFECT and is worth nothing to the test. Note the trap in that list: flavonoids_total, modafinil, multivitamin_mineral, oxytocin_intranasal all have an interval entirely above zero and STILL call NULL_EFFECT, because their point estimate sits below the 0.20 floor. A targeting rule written as "interval excludes zero" rather than "the rule departs from NULL_EFFECT" would waste effort on all four.

Generated by `scripts/136_registry_power.py`.
