"""Diagnostic protocol for the MAMMAL_ONLY_INVERTED finding (Phase A.7 calibration).

Implements the 5 diagnostics + lateral checks from
research/4-tier/Diagnosing MAMMAL DTI Anti-Correlation.md:

  prior_collapse        — does MAMMAL return its training prior at this target?
  power_analysis        — Bonett-Wright Fisher-z CI for the observed ρ
  scaffold_saturation   — Diagnostic A: library overlap with BindingDB top cluster
  distribution_shift    — Diagnostic B: K-S + Wasserstein on pchembl distributions
  tanimoto_correlation  — Diagnostic D (HIGHEST VALUE): ρ(MAMMAL_pKd, max_Tanimoto)
  temporal_strat        — Lateral 6.2: pre-2015 vs post-2015 ChEMBL split
  binding_mode_mix      — Lateral 6.1: parse mechanism_of_action / action_type

ESM2 attention (Diagnostic C) and Boltz positive controls (Diagnostic E)
are GPU-bound and live elsewhere — esm2_attention.py and boltz_positive_control.py
will be added when the overnight Boltz sweep frees the GPU.

Reference: Carroll/Newman RTI series literature for tropane scaffold context;
Karakas et al. 2011 Nature for GluN2B ATD dimer interface; Bonett & Wright
2000 Psychometrika for Fisher-z CI; Zhang et al. 2025 ICLR Oral for low-
Tanimoto degradation; Sundar & Colwell 2020 J Chem Inf Model for scaffold
saturation; Dablander 2023 J Cheminform for activity-cliff failure modes.
"""
