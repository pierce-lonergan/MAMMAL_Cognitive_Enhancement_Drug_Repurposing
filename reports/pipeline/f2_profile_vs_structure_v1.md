# F2 - DTI-profile vs structure for class routing

**Question.** The shipped F2 engine routes by 2D structure (Tanimoto + Murcko scaffold) and abstains when no close analog exists (60% of held-out drugs). Does MAMMAL's learned DTI profile over the cognition panel - the spec's primary signal (a) - route better, and does it RESCUE the structure-abstained compounds? Reproduced by `scripts/96` (GPU scoring) + `scripts/97`.

Profiled ledger compounds: **110** over **31** panel targets (`data/results/f2_dti_profiles.parquet`, MAMMAL on RTX 5070). Evaluable held-out compounds (class keeps a sibling): **89**.

## Leave-one-compound-out class recovery

| mode | top-1 recovery | routed | abstains? |
|---|---|---|---|
| structure-only (shipped) | 0.972 | 36/89 | yes (Tanimoto < 0.35) |
| profile-only (nearest centroid) | 0.112 | 89/89 | no (argmax) |
| blended 50/50 | 0.493 | 71/89 | yes |

Profile-only recovery on the *same* compounds structure routed: **0.194**.

**Rescue test.** Of the 53 compounds structure ABSTAINED on, the profile signal alone recovers the true class for **6%**. The profile does NOT reliably rescue structure-abstained compounds; the MAMMAL affinity profile is too noisy to route these (consistent with the weak affinity-vs-outcome signal, AUROC 0.47). Structure stays primary.

## Why - the profiles are nearly non-selective

Median within-compound spread of predicted pKd across the 31 panel targets is only **0.37** log units: MAMMAL assigns almost the same affinity to every (compound, target) pair on this panel, and structurally distinct drugs (donepezil, galantamine, rivastigmine, pitolisant, modafinil) all share the same top panel targets (CHRM4 / HRH3 / CHRNA7). After within-target z-scoring the residual is mostly noise, so the cross-target profile cannot separate mechanism classes. This is the documented MAMMAL property-correlation bias (predictions track bulk molecular properties more than compound-specific binding); 2D structure, which encodes the actual chemotype, routes far better. A further limit: many ledger classes act on targets absent from the 31-panel (BACE1, GSK3, gamma-secretase, ...), so their profiles cannot be informative even in principle.

## Verdict

**Structure stays primary.** The DTI-profile signal does not beat 2D structure for class routing here; it remains a documented optional hook, not a default. This is itself an honest finding: MAMMAL's binding profile is a weaker class discriminator than chemical structure.
