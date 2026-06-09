# Stage-3 logBB proxy vs measured Kp,uu,brain - anchor validation

Quantifies how well the shipped logBB Stage-3 model (the honest proxy) tracks REAL measured unbound brain exposure, using 10 named marketed drugs with measured Kp,uu,brain parsed from the Heliyon 2024 SI (the only OA experimental source; CC-BY-NC-ND, structures otherwise withheld - see kpuu_data_acquisition_guide.md). Reproduced by `scripts/119_kpuu_anchor_validation.py`. These 10 are a held-out YARDSTICK only - nothing is trained on them, and n=10 so the CI is wide.

## Result (n=10)

- Spearman(predicted logBB, measured Kp,uu) = **0.503**
- Spearman(predicted logBB, log10 Kp,uu) = **0.503**
- logBB gate (PASS) vs true Kp,uu >= 0.3 agreement: **7/10**

Honest read: a positive Spearman means the logBB proxy rank-tracks true unbound exposure to a useful degree even though it is not Kp,uu; a weak/!=1 value is the measured COST of using the proxy and the quantified case for obtaining a licensed Kp,uu spine. n=10 is an anchor, not a benchmark.

## Anchors (measured Kp,uu facts; SMILES from PubChem)

| drug | measured Kp,uu | pred logBB | P-gp | gate | Kp,uu>=0.3 |
|---|---|---|---|---|---|
| methylphenidate | 3.43 | -0.107 | nonsubstrate | PASS | yes |
| hydroxyzine | 1.51 | 0.167 | nonsubstrate | PASS | yes |
| sertraline | 1.44 | 1.59 | nonsubstrate | PASS | yes |
| haloperidol | 1.06 | 1.266 | nonsubstrate | PASS | yes |
| propoxyphene | 0.85 | 0.729 | nonsubstrate | PASS | yes |
| phenacetin | 0.55 | -0.573 | nonsubstrate | ABSTAIN | yes |
| meprobamate | 0.42 | -0.537 | uncertain | ABSTAIN | yes |
| risperidone | 0.26 | -0.039 | uncertain | ABSTAIN | no |
| zolpidem | 0.24 | 0.141 | nonsubstrate | PASS | no |
| sulpiride | 0.07 | -1.05 | uncertain | ABSTAIN | no |

Source of measured Kp,uu: Wu et al. 2024, Heliyon e24304 (PMC10828645), SI Tables S4/S5; values are facts of public marketed drugs, cited, not a redistribution of the CC-BY-NC-ND compiled table. SMILES independently from PubChem.
