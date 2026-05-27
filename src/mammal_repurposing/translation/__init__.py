"""V7 — Clinical Effect-Size Translation Function.

Three subsystems consume V6.A pchembl posteriors + V6.B Cluster D θ̄
posteriors + PBPK exposure to produce predicted healthy-adult cognition
Hedges' g with credible intervals:

  - pbpk.py — 9-compartment JAX/diffrax ODE (gut/plasma/peripheral/cortex/
    striatum/hippocampus/basal-forebrain/brainstem/CSF); Watson 1989
    receptor-occupancy-with-reserve; PET-anchored to Bohnen 2005 donepezil
    19.1% cortical AChE, Volkow 1998 MPH DAT 12-74% across 5-60 mg,
    Kapur 2000 haloperidol D2 striatal EC50.
  - prisma_priors.py — Schmidli 2014 robust meta-analytic-predictive (MAP)
    priors for 12 mechanism classes extracted from Roberts 2020 + MetaPsy +
    Cochrane.
  - effect_size_model.py — PyMC 3-level hierarchical model with Cluster D
    θ̄ multiplicative gate β_target[t_c] = θ̄_{t_c} · β_raw_target[t_c] +
    5 failure-mode moderators + sigmoid translation.

Target venue: Clinical Pharmacology & Therapeutics (Wiley, IF 7.3) or
CPT:PSP for negative-result fallback per V4 §13.Y.
"""

from __future__ import annotations
