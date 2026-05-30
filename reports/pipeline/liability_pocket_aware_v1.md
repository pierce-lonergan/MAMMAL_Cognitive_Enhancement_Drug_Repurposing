# Pocket-Conditioned Liability Gate v1 (§8.13)

Demo run with **synthetic pose scenario = `allosteric_known`**. Real-grid operation awaits §7.17 pose-saving Boltz wrapper.

Per the V4 §8.13 design + research/4-tier/archived/Pocket-Conditioned-Boltz2.md §3.3, the absolute-mode §8.0b CUT is too aggressive when the predicted pose binds OUTSIDE the orthosteric pocket. This gate applies literature-grounded demotion rules:

- **5-HT2B**: Roth 2007 valvulopathy class warning applies to orthosteric agonists (fen-phen / pergolide / cabergoline pattern); allosteric NAMs don't trigger.
- **hERG**: Dumotier & Urban 2024 — central-pore Y652/F656/T623 binding is the classical block; allosteric or vestibule binding is materially lower risk.
- **HRH1**: Gray 2015 anticholinergic dementia risk is orthosteric-antagonist mediated; allosteric H1 ligands have no documented cognition-impairing precedent.
- **CB1**: Topol 2010 CRESCENDO rimonabant neuropsych AEs are orthosteric-only; allosteric NAMs have different clinical risk profile.
- **CHRM1 / OPRM1 / MAOA**: demoted by the same logic.

## Demotable pocket classes per gene

| Gene | Demotable pocket classes |
|---|---|
| HTR2B | allosteric_known, allosteric_putative |
| KCNH2 | allosteric_known, allosteric_putative, no_pocket_match |
| HRH1 | allosteric_known, allosteric_putative |
| CNR1 | allosteric_known, allosteric_putative |
| CHRM1 | allosteric_known, allosteric_putative |
| OPRM1 | allosteric_known, allosteric_putative |
| MAOA | allosteric_known, allosteric_putative |

## Verdict shift (synthetic demo)

| Status | Before | After (`pocket_aware`) | Δ |
|---|---|---|---|
| CUT | 14 | 0 | -14 |
| FLAG | 21 | 35 | +14 |
| PASS | 80 | 80 | +0 |

_Compounds touched_: 14 with non-zero demotions.

## Per-compound demotion detail

| Compound | Original | Pocket-aware | Demotions |
|---|---|---|---|
| 2bact | CUT | FLAG | KCNH2:T1_CUT→T2_FLAG_pocket_aware |
| aripiprazole | CUT | FLAG | CHRM1:T1_CUT→T2_FLAG_pocket_aware; MAOA:T1_CUT→T2_FLAG_pocket_aware |
| bpn14770 | CUT | FLAG | OPRM1:T1_CUT→T2_FLAG_pocket_aware; CNR1:T1_CUT→T2_FLAG_pocket_aware; CHRM1:T1_CUT→T2_FLAG_pocket_aware; MAOA:T1_CUT→T2_FLAG_pocket_aware |
| hydroxyzine | CUT | FLAG | KCNH2:T1_CUT→T2_FLAG_pocket_aware |
| lemborexant | CUT | FLAG | HTR2B:T1_CUT→T2_FLAG_pocket_aware; KCNH2:T1_CUT→T2_FLAG_pocket_aware; HRH1:T1_CUT→T2_FLAG_pocket_aware; CNR1:T1_CUT→T2_FLAG_pocket_aware; CHRM1:T1_CUT→T2_FLAG_pocket_aware; MAOA:T1_CUT→T2_FLAG_pocket_aware |
| lm22a-4 | CUT | FLAG | HTR2B:T1_CUT→T2_FLAG_pocket_aware; OPRM1:T1_CUT→T2_FLAG_pocket_aware; HRH1:T1_CUT→T2_FLAG_pocket_aware; CNR1:T1_CUT→T2_FLAG_pocket_aware; CHRM1:T1_CUT→T2_FLAG_pocket_aware |
| lurasidone | CUT | FLAG | HTR2B:T1_CUT→T2_FLAG_pocket_aware; OPRM1:T1_CUT→T2_FLAG_pocket_aware; HRH1:T1_CUT→T2_FLAG_pocket_aware; CNR1:T1_CUT→T2_FLAG_pocket_aware; CHRM1:T1_CUT→T2_FLAG_pocket_aware; MAOA:T1_CUT→T2_FLAG_pocket_aware |
| methylene blue | CUT | FLAG | HRH1:T1_CUT→T2_FLAG_pocket_aware |
| paroxetine | CUT | FLAG | MAOA:T1_CUT→T2_FLAG_pocket_aware |
| risperidone | CUT | FLAG | HTR2B:T1_CUT→T2_FLAG_pocket_aware; HRH1:T1_CUT→T2_FLAG_pocket_aware; CHRM1:T1_CUT→T2_FLAG_pocket_aware; MAOA:T1_CUT→T2_FLAG_pocket_aware |
| suvorexant | CUT | FLAG | OPRM1:T1_CUT→T2_FLAG_pocket_aware; CNR1:T1_CUT→T2_FLAG_pocket_aware; CHRM1:T1_CUT→T2_FLAG_pocket_aware; MAOA:T1_CUT→T2_FLAG_pocket_aware |
| tc-5619 | CUT | FLAG | HTR2B:T1_CUT→T2_FLAG_pocket_aware; CNR1:T1_CUT→T2_FLAG_pocket_aware; CHRM1:T1_CUT→T2_FLAG_pocket_aware; MAOA:T1_CUT→T2_FLAG_pocket_aware |
| tulrampator | CUT | FLAG | HTR2B:T1_CUT→T2_FLAG_pocket_aware; OPRM1:T1_CUT→T2_FLAG_pocket_aware; HRH1:T1_CUT→T2_FLAG_pocket_aware; CNR1:T1_CUT→T2_FLAG_pocket_aware; CHRM1:T1_CUT→T2_FLAG_pocket_aware; MAOA:T1_CUT→T2_FLAG_pocket_aware |
| xen-1101 | CUT | FLAG | KCNH2:T1_CUT→T2_FLAG_pocket_aware |

## How to operationalise

This script's `pose_df` is **synthetic** — every CUT compound is assigned the scenario pocket-class. To run the gate against real Boltz-2 poses:

```python
# 1. Save mmCIF poses during the Boltz sweep (§7.17 unblocks this).
# 2. Extract heavy-atom centroid xyz per pose.
# 3. For each (compound, target_gene) classify via §7.5:
from mammal_repurposing.pockets.pocket_classifier import classify_pose
from mammal_repurposing.pockets.pocket_database import load_pocket_database
db = load_pocket_database('data/pockets/centroids/')
pose_df = (pose_centroids_df.apply(lambda r: 
    classify_pose(r[['x','y','z']].values, r['target_gene'], db),
    axis=1))
# 4. Then call pocket_aware_liability_gate(liability_gates_df, pose_df).
```

Currently only 7 of 22 cognition targets have curated centroids (CHRNA7, ACHE, HRH3, DRD1, PDE4D, SIGMAR1, GRIN2B). To extend to the full 44-target liability panel, either: (a) reuse a single canonical orthosteric pocket per liability target via UniProt → reference PDB lookup, OR (b) extend the curated DB to 44 targets (~2 weeks of centroid curation).

---

Generated by `scripts/39_v5_pocket_conditional_liability.py` (scenario=allosteric_known).