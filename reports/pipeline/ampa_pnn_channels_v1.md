# AMPA-PAM + PNN/ECM persistence channels - pre-registered calibration

Two candidate structure-derivable durability channels, calibrated against the size-matched negative pool. Reproduced by `scripts/114_ampa_pnn_channels.py`.

**PRE-REGISTERED:** AMPA-PAM FAILS (the durability lever is allosteric potentiation; MAMMAL's BindingDB-pKd head is allosterically blind - v1 CHRNA7 PAM std 0.029). AMPA orthosteric may rank; MMP9/PNN is exploratory and must clear the same gate.

| target | channel/site | engagers | AUROC | perm-p | passes |
|---|---|---|---|---|---|
| GRIA1 | AMPA/orthosteric | glutamate, quisqualate, kainic acid | 0.09 | 0.992 | **fail** |
| GRIA1 | AMPA/pam | cyclothiazide, aniracetam, CX-516, CX-717 | 0.26 | 0.945 | **fail** |
| MMP9 | PNN_ECM/mmp_inhibitor | marimastat, batimastat, prinomastat, ilomastat, doxycycline | 0.42 | 0.717 | **fail** |

## Reading

CONFIRMED: MAMMAL FAILS the AMPA-PAM (allosteric-potentiation) site - the durability lever is invisible to a BindingDB-pKd head, extending the v1 allosteric-blindness audit to the AMPA channel. The orthosteric/MMP9 rows are reported for completeness; any channel that does not clear the size-matched + permutation gate is NOT wired into the persistence head (abstain-by-default). Net: AMPA durability (PAM/phospho/mTORC1-coupled) and PNN/ECM remodeling are confirmed OFF the sequence-DTI axis, consistent with routing plasticity through the L4 permeability window and reserving these tiers for a structure/allosteric-aware second opinion.
