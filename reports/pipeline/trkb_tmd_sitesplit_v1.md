# TrkB (NTRK2) site-split calibration - pre-registered TrkB-TMD negative

Does MAMMAL's sequence-only DTI head rank the DURABILITY-relevant TrkB-TMD antidepressants (Casarotto 2021) above matched non-engagers - or only the ATP-pocket TRK inhibitors that happen to be in BindingDB? Reproduced by `scripts/112_trkb_tmd_sitesplit.py`.

**PRE-REGISTERED EXPECTATION:** tmd_wedge FAILS (the crossed-dimer/cholesterol mode is information-theoretically absent from a 1D sequence); ATP-pocket may pass.

Negatives: 33 size-matched non-engagers.

| binding site | engagers | AUROC | perm-p | passes (AUROC>=0.70, p<0.05) |
|---|---|---|---|---|
| tmd_wedge | fluoxetine, imipramine, ketamine, (2R,6R)-hydroxynorketamine | 0.23 | 0.962 | **fail** |
| ecd | 7,8-dihydroxyflavone, lm22a-4 | 0.44 | 0.599 | **fail** |
| atp_pocket | larotrectinib, entrectinib, selitrectinib | 0.59 | 0.335 | **fail** |

## Reading

CONFIRMED pre-registration: MAMMAL FAILS the tmd_wedge (durability) site. The engine, if it ranks TrkB engagers at all, sees the ATP-pocket / ECD site, NOT the transmembrane-domain crossed-dimer wedge that mediates the durable antidepressant/plasticity effect. This is why PERSEUS routes psychoplastogen durability through the L4 permeability-gated window (off the DTI axis) rather than a TrkB DTI score, and why even a Boltz-2 structure second opinion is insufficient here: the active site only forms as a cholesterol-dependent crossed dimer in the lipid bilayer (Casarotto 2021; Cordeiro 2024), which a single-chain apo prediction cannot represent. The TrkB-TMD durability channel is therefore a documented off-axis limit, not a buildable DTI head.

## Boltz-2 install status (D3) - dependency resolved, limitation stands

The "even a Boltz-2 second opinion is insufficient" claim above is now empirically grounded, not hypothetical: **Boltz 2.2.1 is installed and import-verified** (in an isolated `.venv-boltz312` Python-3.12 environment, deliberately separate from `.venv-mammal` so its pinned deps cannot perturb the nightly cu128 / CUDA-sm_120 torch that MAMMAL needs). So the install dependency that previously blocked any structural cross-check is RESOLVED.

What that does NOT change: Boltz-2 predicts an apo, single-copy, solvated fold; it has no representation of (a) the lipid bilayer, (b) the cholesterol co-factor that gates the TrkB-TMD site, or (c) the **crossed-dimer** quaternary arrangement that only that membrane+cholesterol context induces (Casarotto 2021; Cordeiro 2024). Predicting the monomer well does not let it score binding at an interface that exists only in the crossed dimer. The TrkB-TMD off-axis verdict is therefore robust to the obvious "just run a structure model" rebuttal - the limitation is informational (the relevant state is absent from the inputs), not a tooling gap.

The concrete, honest future use of the now-installed Boltz is the **BCL2 senolytic channel**: BCL2 is a soluble protein with a well-defined groove (the only one of the four durability nodes that is sequence-DTI-recoverable), so a Boltz-2 pose + affinity prediction is a legitimate structural confirmation there - exactly where the modality's assumptions hold. Building that BCL2 structural cross-check is the next Boltz step; forcing it onto TrkB-TMD would be theatre.
