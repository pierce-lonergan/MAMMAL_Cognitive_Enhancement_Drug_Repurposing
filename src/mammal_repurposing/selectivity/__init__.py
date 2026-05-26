"""Selectivity scoring for the cognition panel.

Implements §7.4 of the V3 forward plan, with concrete spec from
research/4-tier/Graczyk-Style Gini + S(10x) Selectivity Scoring and
Multi-Class Faceted Wet-Lab Shortlist Generator ....md.

Modules:
  gini_scorecard.py  — Graczyk Gini (2007, J Med Chem 50:5773) + S(10x)
                       Karaman-style + selectivity vector with rank-percentile
                       substitution at MAMMAL_ONLY_INVERTED targets
  categorize.py      — mono/dual/poly/flat assignment + provenance flags
  cross_cluster.py   — MAMMAL-Gini vs Boltz-Gini diagnostic
"""

from .gini_scorecard import (
    gini,
    s_10x,
    selectivity_vector,
    bootstrap_gini_ci,
    score_compound,
    score_panel,
    PANEL_22,
    INVERTED,
    PriorByTarget,
)
from .categorize import categorize, SelectivityCategory

__all__ = [
    "gini",
    "s_10x",
    "selectivity_vector",
    "bootstrap_gini_ci",
    "score_compound",
    "score_panel",
    "categorize",
    "SelectivityCategory",
    "PANEL_22",
    "INVERTED",
    "PriorByTarget",
]
