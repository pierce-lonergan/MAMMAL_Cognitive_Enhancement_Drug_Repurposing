"""Per-candidate provenance: cluster contributions, disagreement diagnosis, narrative.

Creative addition to the v2 hybrid (not in the research doc). For every compound
that survives to the final ranking, we record WHICH cluster placed it where +
WHY it survived the gates. Drives three outputs:
    - provenance.parquet     (machine-readable, per-row)
    - disagreement_report.md (model-disagreement archetypes per the doc §7)
    - funnel_narrative.md    (per-candidate prose explanation for top-N)
"""
