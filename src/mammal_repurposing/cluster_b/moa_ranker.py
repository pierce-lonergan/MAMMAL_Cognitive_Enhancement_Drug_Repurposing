"""§8.7 — Mechanism-of-Action (MoA) preference ranker.

Promotes the v4 `diagnostics/binding_mode_mix.py` diagnostic to a fusion-input
ranker. For each (compound, target) it looks up the ChEMBL `action_type` and
`mechanism_of_action` columns, then scores how well the compound's annotated
MoA matches the preferred-MoA-for-cognition for that target.

Per-target preferences are hand-curated from the V3/V4 cognition-target panel
literature. The score is in [0, 1]:
    1.0  — annotated MoA matches the target's preferred MoA for cognition
    0.5  — MoA is mechanism-relevant but not the preferred direction
    0.0  — wrong-direction (antagonist where agonist is wanted) or unknown

Used as a fifth ranker in the 5-cluster RRF when --add-moa-ranker is passed.

Reference: ChEMBL action_type taxonomy; Bowes 2012 Nat Rev Drug Discov 11:909
on Safety-44 MoA classes; cognition-target MoA preferences from V3 doc §2.

NOTE: this is a *prior* — the compound's MoA at OTHER targets (in §8.0b
liability panel) is independent. If a CHRNA7 PAM is also a 5-HT2B agonist
that's captured in liability, not here.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd

logger = logging.getLogger(__name__)


# Per-target preferred MoA for cognitive enhancement.
# Keys are gene symbols (matching PANEL_22). Values are dict of:
#   {action_type_string: score_in_[0,1]}.
# action_type values are ChEMBL controlled vocabulary (UPPERCASE).
COGNITION_PREFERRED_MOA: dict[str, dict[str, float]] = {
    # Cholinergic — AChE inhibition raises ACh tone (donepezil/galantamine/rivastigmine)
    "ACHE":    {"INHIBITOR": 1.0, "NEGATIVE ALLOSTERIC MODULATOR": 0.7},
    # α7 nAChR — PAMs (encenicline, TC-5619 family) preferred over agonists
    "CHRNA7":  {"POSITIVE ALLOSTERIC MODULATOR": 1.0,
                "AGONIST": 0.7, "PARTIAL AGONIST": 0.7,
                "ANTAGONIST": 0.0},
    # AMPA receptors — PAMs (ampakines) preferred
    "GRIA1":   {"POSITIVE ALLOSTERIC MODULATOR": 1.0,
                "AGONIST": 0.6, "PARTIAL AGONIST": 0.6,
                "ANTAGONIST": 0.0},
    "GRIA2":   {"POSITIVE ALLOSTERIC MODULATOR": 1.0,
                "AGONIST": 0.6, "ANTAGONIST": 0.0},
    "GRIA3":   {"POSITIVE ALLOSTERIC MODULATOR": 1.0,
                "AGONIST": 0.6, "ANTAGONIST": 0.0},
    "GRIA4":   {"POSITIVE ALLOSTERIC MODULATOR": 1.0,
                "AGONIST": 0.6, "ANTAGONIST": 0.0},
    # NMDA — GluN2B-selective antagonists / NAMs (memantine / ifenprodil class)
    "GRIN2A":  {"NEGATIVE ALLOSTERIC MODULATOR": 1.0, "ANTAGONIST": 0.7,
                "AGONIST": 0.0},
    "GRIN2B":  {"NEGATIVE ALLOSTERIC MODULATOR": 1.0, "ANTAGONIST": 0.7,
                "AGONIST": 0.0},
    # Dopaminergic — D1 PAMs / D1 partial agonists for working memory
    "DRD1":    {"POSITIVE ALLOSTERIC MODULATOR": 1.0,
                "PARTIAL AGONIST": 0.8, "AGONIST": 0.6},
    # DAT — inhibitors (methylphenidate / modafinil)
    "SLC6A3":  {"INHIBITOR": 1.0, "RELEASER": 0.7,
                "BLOCKER": 1.0},
    # NET — inhibitors (atomoxetine / reboxetine)
    "SLC6A2":  {"INHIBITOR": 1.0, "BLOCKER": 1.0},
    # α2A — agonists (guanfacine, clonidine)
    "ADRA2A":  {"AGONIST": 1.0, "PARTIAL AGONIST": 0.8,
                "ANTAGONIST": 0.2},
    # H3 — inverse agonists / antagonists (pitolisant)
    "HRH3":    {"INVERSE AGONIST": 1.0, "ANTAGONIST": 1.0,
                "AGONIST": 0.0},
    # Orexin receptors — direction depends on indication:
    #   narcolepsy: AGONIST; insomnia: ANTAGONIST.
    # For COGNITION (wake-promotion path), AGONIST is preferred.
    "HCRTR1":  {"AGONIST": 1.0, "PARTIAL AGONIST": 0.8,
                "ANTAGONIST": 0.3},
    "HCRTR2":  {"AGONIST": 1.0, "PARTIAL AGONIST": 0.8,
                "ANTAGONIST": 0.3},
    # PDE inhibitors raise cAMP/cGMP (rolipram, BPN14770, PF-04447943)
    "PDE4D":   {"INHIBITOR": 1.0},
    "PDE9A":   {"INHIBITOR": 1.0},
    # TrkB — agonists / PAMs (7,8-DHF, LM22A-4)
    "NTRK2":   {"AGONIST": 1.0, "POSITIVE ALLOSTERIC MODULATOR": 1.0,
                "PARTIAL AGONIST": 0.7, "ANTAGONIST": 0.0},
    # Sigma-1 — agonists (pridopidine, blarcamesine, fluvoxamine)
    "SIGMAR1": {"AGONIST": 1.0, "PARTIAL AGONIST": 0.8,
                "ANTAGONIST": 0.2},
    # Kv7 (KCNQ2/3) — openers / activators (retigabine, XEN-1101)
    "KCNQ2":   {"OPENER": 1.0, "ACTIVATOR": 1.0,
                "POSITIVE ALLOSTERIC MODULATOR": 0.8,
                "INHIBITOR": 0.0, "BLOCKER": 0.0},
    "KCNQ3":   {"OPENER": 1.0, "ACTIVATOR": 1.0,
                "POSITIVE ALLOSTERIC MODULATOR": 0.8,
                "INHIBITOR": 0.0, "BLOCKER": 0.0},
    # HCN1 — blockers (ivabradine direction, but for cognition this is unusual)
    # Conservative default: don't downweight; mild preference for INHIBITOR
    "HCN1":    {"INHIBITOR": 0.5, "BLOCKER": 0.5},
}

# Soft default for any (compound, target) where MoA annotation exists but the
# annotated action_type isn't in the preferred map.
DEFAULT_NEUTRAL_SCORE = 0.5


@dataclass
class MoaRankerConfig:
    """Configuration for the MoA ranker.

    `default_score`: score for compounds with no ChEMBL MoA annotation at the
        target — neutral (0.5) rather than penalising lack of data.
    """
    default_score: float = 0.5
    unknown_score: float = 0.4


def _normalize_action_type(at: str | None) -> str:
    if not isinstance(at, str):
        return ""
    return at.strip().upper()


def score_moa_at_target(
    target_gene: str,
    action_type: str | None,
    mechanism_of_action: str | None = None,
    config: MoaRankerConfig | None = None,
) -> float:
    """Return MoA score in [0, 1] for one (compound, target) given its
    annotated `action_type` / `mechanism_of_action`.

    Resolution order:
      1. action_type matches preferred → score from table
      2. action_type present but not in table → DEFAULT_NEUTRAL_SCORE
      3. action_type missing → unknown_score (slightly below neutral)
    """
    cfg = config or MoaRankerConfig()
    preferred = COGNITION_PREFERRED_MOA.get(target_gene)
    if preferred is None:
        return cfg.default_score

    at_norm = _normalize_action_type(action_type)
    if at_norm and at_norm in preferred:
        return float(preferred[at_norm])

    # Fall-back: scan mechanism_of_action free text for keywords
    if isinstance(mechanism_of_action, str):
        moa_lower = mechanism_of_action.lower()
        for at_key, sc in preferred.items():
            if at_key.lower() in moa_lower:
                return float(sc)

    if at_norm:
        return float(DEFAULT_NEUTRAL_SCORE)
    return float(cfg.unknown_score)


def build_moa_ranker_long(
    library_df: pd.DataFrame,                 # cols: compound_name, [inchikey]
    target_uniprots: list[str],
    gene_by_uniprot: dict[str, str],
    chembl_moa_loader,                         # uniprot -> DataFrame with columns inchikey, action_type, mechanism_of_action
    config: MoaRankerConfig | None = None,
    ranker_name: str = "cluster_b_moa",
) -> pd.DataFrame:
    """Score every (target, compound) in the grid via the MoA preference table.

    The compound→inchikey lookup uses library_df['inchikey'] if present, else
    the compound name is matched against the ChEMBL annotation table by
    fuzzy join on lowercased name (last resort; usually inchikey).

    Returns long-format DataFrame compatible with the RRF fusion input shape:
        target_uniprot, compound_name, predicted_pkd, ranker_name
    """
    cfg = config or MoaRankerConfig()
    rows: list[pd.DataFrame] = []

    # Build a per-compound inchikey lookup if available
    have_inchi = "inchikey" in library_df.columns
    if have_inchi:
        inchi_map = dict(zip(
            library_df["compound_name"].str.lower().str.strip(),
            library_df["inchikey"].astype(str),
        ))
    else:
        inchi_map = {}

    for u in target_uniprots:
        gene = gene_by_uniprot.get(u, "")
        moa_df = chembl_moa_loader(u)
        if moa_df is None or moa_df.empty:
            # No annotated MoA for this target; fill all compounds with the default
            scores = [cfg.default_score] * len(library_df)
            rows.append(pd.DataFrame({
                "target_uniprot": u,
                "compound_name": library_df["compound_name"].tolist(),
                "predicted_pkd": scores,
                "ranker_name": ranker_name,
            }))
            continue

        # Build per-inchikey best-MoA-score map
        moa_df_norm = moa_df.copy()
        # Plain dict + (-inf) sentinel, NOT defaultdict: indexing a defaultdict
        # materializes the key at unknown_score (0.4) BEFORE the comparison, so a
        # wrong-direction annotation (sc=0.0, e.g. an antagonist at a target that wants
        # agonism) never beats 0.4 and gets floored to 0.4 == "no annotation" instead of
        # the correct 0.0. The (-inf) sentinel lets the first real score (incl. 0.0) win,
        # and unannotated compounds still fall through to unknown_score below (the key is
        # only present when a real annotation set it).
        per_compound_score: dict[str, float] = {}
        for _, r in moa_df_norm.iterrows():
            ik = str(r.get("inchikey") or "").strip()
            if not ik:
                continue
            sc = score_moa_at_target(
                gene,
                r.get("action_type"),
                r.get("mechanism_of_action"),
                cfg,
            )
            # Keep the most-favourable annotation per compound
            if sc > per_compound_score.get(ik, float("-inf")):
                per_compound_score[ik] = sc

        # Score library compounds
        target_scores: list[float] = []
        for _, lib_r in library_df.iterrows():
            name_lc = str(lib_r["compound_name"]).lower().strip()
            ik = inchi_map.get(name_lc, "")
            if ik and ik in per_compound_score:
                target_scores.append(per_compound_score[ik])
            else:
                target_scores.append(cfg.unknown_score)

        rows.append(pd.DataFrame({
            "target_uniprot": u,
            "compound_name": library_df["compound_name"].tolist(),
            "predicted_pkd": target_scores,
            "ranker_name": ranker_name,
        }))
        logger.info("  %s (%s): %d compounds scored; %d distinct MoA annotations",
                    gene, u, len(target_scores), len(per_compound_score))
    return pd.concat(rows, ignore_index=True)
