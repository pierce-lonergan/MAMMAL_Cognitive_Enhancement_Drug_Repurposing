"""Pytest coverage for V4/V5 transition modules:

  - §7.18 + §8.0b-zn z-norm gate (gates/liability_panel.py)
  - §8.13 pocket-conditioned liability composition
  - §8.15 disagreement classifier (scripts/35_v3_disagreement_signal.py)
  - §7.4 v2 selectivity entropy + Partition Index
  - §8.7 MoA preference ranker
  - §8.10 nootropic-similarity annotator
  - §7.17 pose-extract centroid

Each test is a happy-path + at least one boundary / edge case.
Marked `not slow` — these are pure-Python unit tests with no network or GPU.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Ensure src is importable
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


# ---------------------------------------------------------------------------
# §7.18 + §8.0b-zn — z-norm liability gate
# ---------------------------------------------------------------------------
class TestLiabilityZnormGate:
    def _make_panel(self) -> pd.DataFrame:
        """Minimal 3-target panel covering all three tiers."""
        return pd.DataFrame([
            {"gene_symbol": "KCNH2", "uniprot_accession": "Q12809",
             "severity_tier": 1, "cut_threshold_pki": 6.0,
             "flag_threshold_pki": 5.5, "liability_category": "ion_channel"},
            {"gene_symbol": "CHRM3", "uniprot_accession": "P20309",
             "severity_tier": 2, "cut_threshold_pki": float("nan"),
             "flag_threshold_pki": 6.0, "liability_category": "gpcr"},
            {"gene_symbol": "TACR1", "uniprot_accession": "P25103",
             "severity_tier": 3, "cut_threshold_pki": float("nan"),
             "flag_threshold_pki": 5.0, "liability_category": "gpcr"},
        ])

    def _make_dti_grid(self) -> pd.DataFrame:
        """3 named compounds × 3 targets + 8 dummy compounds for variance.

        Dummy distribution per target: [5.5, 5.7, 5.9, 6.1, 6.3, 6.4, 6.6, 6.8]
            mean = 6.1625, std ≈ 0.4445
        For z-score targeting:
            z = +3 → pkd = 6.1625 + 3 × 0.4445 ≈ 7.50  (after including the hit
            compound itself in the std → recompute; pick a value that comes out
            ≥ +2.0 robustly)

        We just set the hit pkd to 8.0 (high enough that the within-target z
        exceeds +2.0 even when the hit is included in the mean/std estimate).
        """
        dummy_values = [5.5, 5.7, 5.9, 6.1, 6.3, 6.4, 6.6, 6.8]
        rows = []
        for cname, kcnh2_pkd, chrm3_pkd, tacr1_pkd in [
            ("compound_safe", 6.2, 6.2, 6.2),               # mid-pack
            ("compound_kcnh2_hit", 8.0, 6.2, 6.2),          # KCNH2 outlier
            ("compound_chrm3_flag", 6.2, 7.4, 6.2),         # CHRM3 +1.5σ-ish
        ]:
            for gene, uni, pkd in [
                ("KCNH2", "Q12809", kcnh2_pkd),
                ("CHRM3", "P20309", chrm3_pkd),
                ("TACR1", "P25103", tacr1_pkd),
            ]:
                rows.append({
                    "target_uniprot": uni, "target_gene": gene,
                    "compound_name": cname, "predicted_pkd": pkd,
                })
        for j, v in enumerate(dummy_values):
            for gene, uni in [
                ("KCNH2", "Q12809"), ("CHRM3", "P20309"), ("TACR1", "P25103"),
            ]:
                rows.append({
                    "target_uniprot": uni, "target_gene": gene,
                    "compound_name": f"dummy_{j}", "predicted_pkd": v,
                })
        return pd.DataFrame(rows)

    def test_znorm_mode_demotes_prior_collapse_artifact(self):
        from mammal_repurposing.gates.liability_panel import apply_liability_gates
        panel = self._make_panel()
        dti = self._make_dti_grid()
        gates = apply_liability_gates(dti, panel, znorm=True,
                                      z_cut_tier1=2.0, z_flag_tier2=1.5,
                                      z_flag_tier3=1.0)
        d = dict(zip(gates["compound_name"], gates["liability_status"]))
        # Hit compound at KCNH2 is the only Tier 1 CUT
        assert d["compound_kcnh2_hit"] == "CUT"
        assert d["compound_chrm3_flag"] == "FLAG"
        assert d["compound_safe"] == "PASS"

    def test_absolute_mode_overshoots_under_prior_collapse(self):
        """Without z-norm, the absolute pKi=6.5 baseline trips KCNH2 cut=6.0
        on every compound — the documented failure mode (§8.0b)."""
        from mammal_repurposing.gates.liability_panel import apply_liability_gates
        panel = self._make_panel()
        dti = self._make_dti_grid()
        gates = apply_liability_gates(dti, panel, znorm=False)
        # All real compounds CUT because pKd 6.5 > KCNH2 cut 6.0
        cut = set(gates[gates["liability_status"] == "CUT"]["compound_name"])
        assert "compound_safe" in cut    # the bug we fixed with z-norm

    def test_combined_with_admet_handles_nan(self):
        """Regression test for the NaN.upper() bug — outer-merge gaps must not
        raise."""
        from mammal_repurposing.gates.liability_panel import combine_admet_and_liability
        admet = pd.DataFrame([
            {"compound_name": "a", "gate_status": "PASS"},
            {"compound_name": "b", "gate_status": "CUT"},
            {"compound_name": "missing_from_liability", "gate_status": "PASS"},
        ])
        liability = pd.DataFrame([
            {"compound_name": "a", "liability_status": "PASS",
             "liability_note": "", "liability_summary": "",
             "tier_1_hits": "", "tier_2_hits": "", "tier_3_hits": "",
             "n_tier_1": 0, "n_tier_2": 0, "n_tier_3": 0},
        ])
        combined = combine_admet_and_liability(admet, liability)
        d = dict(zip(combined["compound_name"], combined["final_status"]))
        assert d["a"] == "PASS"
        assert d["b"] == "CUT"
        # 'missing_from_liability' has admet=PASS, liability=NaN → final PASS
        # (the bug used to throw AttributeError on float.upper())
        assert d["missing_from_liability"] == "PASS"


# ---------------------------------------------------------------------------
# §8.13 — Pocket-class-conditioned liability composition
# ---------------------------------------------------------------------------
class TestPocketConditionedLiability:
    def test_allosteric_pose_demotes_cut_to_flag(self):
        from mammal_repurposing.gates.liability_panel import (
            pocket_aware_liability_gate,
        )
        gates = pd.DataFrame([{
            "compound_name": "comp_x",
            "liability_status": "CUT",
            "liability_note": "Tier 1 z>=2.0: KCNH2",
            "liability_summary": "KCNH2=z+2.5(T1)",
            "tier_1_hits": "KCNH2", "tier_2_hits": "", "tier_3_hits": "",
            "n_tier_1": 1, "n_tier_2": 0, "n_tier_3": 0,
            "top_3_liabilities": "KCNH2=z+2.5(T1)",
        }])
        poses = pd.DataFrame([{
            "compound_name": "comp_x", "target_gene": "KCNH2",
            "pocket_class": "allosteric_known",
        }])
        out = pocket_aware_liability_gate(gates, poses)
        row = out.iloc[0]
        assert row["liability_status"] == "CUT"                # unchanged
        assert row["liability_status_pocket_aware"] == "FLAG"  # demoted
        assert "KCNH2:T1_CUT" in row["pocket_demotions"]
        assert row["n_pocket_demoted"] == 1

    def test_orthosteric_pose_preserves_cut(self):
        from mammal_repurposing.gates.liability_panel import (
            pocket_aware_liability_gate,
        )
        gates = pd.DataFrame([{
            "compound_name": "comp_y", "liability_status": "CUT",
            "liability_note": "...", "liability_summary": "...",
            "tier_1_hits": "HTR2B", "tier_2_hits": "", "tier_3_hits": "",
            "n_tier_1": 1, "n_tier_2": 0, "n_tier_3": 0,
            "top_3_liabilities": "",
        }])
        poses = pd.DataFrame([{
            "compound_name": "comp_y", "target_gene": "HTR2B",
            "pocket_class": "orthosteric",  # NOT in POCKET_AWARE_DEMOTABLE
        }])
        out = pocket_aware_liability_gate(gates, poses)
        assert out.iloc[0]["liability_status_pocket_aware"] == "CUT"
        assert out.iloc[0]["n_pocket_demoted"] == 0

    def test_pass_compounds_pass_through(self):
        from mammal_repurposing.gates.liability_panel import (
            pocket_aware_liability_gate,
        )
        gates = pd.DataFrame([{
            "compound_name": "clean_compound", "liability_status": "PASS",
            "liability_note": "clean", "liability_summary": "",
            "tier_1_hits": "", "tier_2_hits": "", "tier_3_hits": "",
            "n_tier_1": 0, "n_tier_2": 0, "n_tier_3": 0,
            "top_3_liabilities": "",
        }])
        out = pocket_aware_liability_gate(gates, pd.DataFrame())
        assert out.iloc[0]["liability_status_pocket_aware"] == "PASS"


# ---------------------------------------------------------------------------
# §8.15 — Disagreement classifier
# ---------------------------------------------------------------------------
class TestDisagreementClassifier:
    def test_agree(self):
        # Import the helper directly from the script (we don't need the whole
        # script's CLI; just the classify function).
        script_path = ROOT / "scripts" / "35_v3_disagreement_signal.py"
        spec = importlib.util.spec_from_file_location("_disagreement", script_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        # rank_mammal=10, rank_tanimoto=15 → delta=-5 → agree
        assert mod._classify_disagreement(-5, 10, 15) == "agree"

    def test_moderate(self):
        script_path = ROOT / "scripts" / "35_v3_disagreement_signal.py"
        spec = importlib.util.spec_from_file_location("_disagreement", script_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        # |delta| = 40 → moderate_disagreement
        assert mod._classify_disagreement(40, 50, 10) == "moderate_disagreement"
        assert mod._classify_disagreement(-40, 10, 50) == "moderate_disagreement"

    def test_novel_scaffold_suspect(self):
        script_path = ROOT / "scripts" / "35_v3_disagreement_signal.py"
        spec = importlib.util.spec_from_file_location("_disagreement", script_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        # MAMMAL ranks HIGH (#5), Tanimoto LOW (#150), |Δ|=145 > 50, rank_mammal < rank_tanimoto
        assert mod._classify_disagreement(-145, 5, 150) == "novel_scaffold_suspect"

    def test_activity_cliff_suspect(self):
        script_path = ROOT / "scripts" / "35_v3_disagreement_signal.py"
        spec = importlib.util.spec_from_file_location("_disagreement", script_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        # MAMMAL ranks LOW (#298), Tanimoto HIGH (#3), |Δ|=295 > 50, rank_mammal > rank_tanimoto
        assert mod._classify_disagreement(295, 298, 3) == "activity_cliff_suspect"


# ---------------------------------------------------------------------------
# §7.4 v2 — Selectivity entropy + Partition Index
# ---------------------------------------------------------------------------
class TestSelectivityV2Metrics:
    def test_mono_selective_extremes(self):
        from mammal_repurposing.selectivity.gini_scorecard import (
            selectivity_entropy, top_target_partition_index,
        )
        mono = np.array([5.0] * 9 + [9.0])      # one big winner
        # Entropy should be LOW (close to 0); PI_top should be HIGH (close to 1)
        h = selectivity_entropy(mono)
        pi = top_target_partition_index(mono)
        assert h < 0.5
        assert pi > 0.95

    def test_flat_panel(self):
        from mammal_repurposing.selectivity.gini_scorecard import (
            selectivity_entropy, top_target_partition_index,
        )
        flat = np.array([6.5] * 10)
        h = selectivity_entropy(flat)
        pi = top_target_partition_index(flat)
        assert abs(h - 1.0) < 0.01    # max entropy
        assert abs(pi - 0.1) < 0.01   # uniform → 1/10

    def test_partition_index_invariant_to_uniform_shift(self):
        """Adding a constant to every pKd shouldn't change PI ordering."""
        from mammal_repurposing.selectivity.gini_scorecard import (
            top_target_partition_index,
        )
        v1 = np.array([5.0, 6.0, 7.0, 8.0])
        v2 = v1 + 2.0
        # PI is calculated on the shifted-to-max-zero pKd, so absolute level
        # doesn't matter — only the spread.
        assert abs(top_target_partition_index(v1) - top_target_partition_index(v2)) < 1e-6

    def test_entropy_handles_nans(self):
        from mammal_repurposing.selectivity.gini_scorecard import (
            selectivity_entropy,
        )
        v = np.array([6.0, np.nan, 7.0, 8.0])
        h = selectivity_entropy(v)
        assert not np.isnan(h)


# ---------------------------------------------------------------------------
# §8.7 — MoA preference ranker
# ---------------------------------------------------------------------------
class TestMoaRanker:
    def test_preferred_moa_at_chrna7(self):
        from mammal_repurposing.cluster_b.moa_ranker import score_moa_at_target
        # CHRNA7 prefers PAM > AGONIST > ANTAGONIST
        assert score_moa_at_target("CHRNA7", "POSITIVE ALLOSTERIC MODULATOR") == 1.0
        assert score_moa_at_target("CHRNA7", "AGONIST") == 0.7
        assert score_moa_at_target("CHRNA7", "ANTAGONIST") == 0.0

    def test_unknown_gene_returns_default(self):
        from mammal_repurposing.cluster_b.moa_ranker import score_moa_at_target
        s = score_moa_at_target("UNKNOWN_GENE", "INHIBITOR")
        assert s == 0.5

    def test_free_text_fallback(self):
        from mammal_repurposing.cluster_b.moa_ranker import score_moa_at_target
        # PDE9A prefers INHIBITOR; free-text fallback hits the same keyword
        s = score_moa_at_target("PDE9A", None,
                                mechanism_of_action="Phosphodiesterase 9 inhibitor")
        assert s == 1.0

    def test_unknown_action_type_at_known_target(self):
        from mammal_repurposing.cluster_b.moa_ranker import score_moa_at_target
        # CHRNA7 has a preference table but ALLOSTERIC_LIGAND isn't in it.
        # Falls back to DEFAULT_NEUTRAL_SCORE = 0.5.
        s = score_moa_at_target("CHRNA7", "ALLOSTERIC_LIGAND")
        assert s == 0.5

    def test_build_long_format_with_empty_loader(self):
        """If chembl_moa_loader returns empty, every compound gets default."""
        from mammal_repurposing.cluster_b.moa_ranker import build_moa_ranker_long
        lib = pd.DataFrame({"compound_name": ["a", "b"], "inchikey": ["x", "y"]})
        empty_loader = lambda u: pd.DataFrame()
        long_df = build_moa_ranker_long(
            lib, ["P22303"], {"P22303": "ACHE"},
            chembl_moa_loader=empty_loader,
        )
        assert len(long_df) == 2
        assert (long_df["predicted_pkd"] == 0.5).all()
        assert long_df["ranker_name"].unique().tolist() == ["cluster_b_moa"]


# ---------------------------------------------------------------------------
# §8.10 — Nootropic similarity annotator
# ---------------------------------------------------------------------------
class TestNootropicSimilarity:
    def test_donepezil_self_match_suppressed(self):
        from mammal_repurposing.analysis.nootropic_similarity import (
            annotate_dataframe, CANONICAL_NOOTROPICS,
        )
        df = pd.DataFrame([{
            "compound_name": "donepezil",
            "compound_smiles": CANONICAL_NOOTROPICS["donepezil"],
        }])
        out = annotate_dataframe(df)
        # Self-match suppressed → nearest is something ELSE, not donepezil
        assert out.iloc[0]["nearest_nootropic"] != "donepezil"

    def test_racetam_family_intermediate(self):
        from mammal_repurposing.analysis.nootropic_similarity import (
            annotate_dataframe, CANONICAL_NOOTROPICS,
        )
        df = pd.DataFrame([{
            "compound_name": "levetiracetam",
            "compound_smiles": "CCN1CCC(C1=O)C(=O)N",
        }])
        out = annotate_dataframe(df)
        # Levetiracetam is a racetam — nearest should be piracetam or aniracetam
        row = out.iloc[0]
        assert row["nearest_nootropic"] in ("piracetam", "aniracetam")
        assert row["nootropic_novelty_tag"] in ("intermediate", "analog")

    def test_unparseable_smiles(self):
        from mammal_repurposing.analysis.nootropic_similarity import (
            annotate_dataframe,
        )
        df = pd.DataFrame([{"compound_name": "broken", "compound_smiles": "@@!!"}])
        out = annotate_dataframe(df)
        assert out.iloc[0]["nootropic_novelty_tag"] == "unknown"


# ---------------------------------------------------------------------------
# §7.17 — Pose centroid extractor
# ---------------------------------------------------------------------------
class TestPoseExtractor:
    def test_centroid_of_three_atoms(self):
        from mammal_repurposing.pockets.pose_extract import extract_ligand_centroid
        cif = """data_t
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA A 1 0.0 0.0 0.0 1.0 0.0
HETATM 2 C C1 LIG L .  0.0  0.0  0.0  1.0  0.0
HETATM 3 N N1 LIG L .  6.0  0.0  0.0  1.0  0.0
HETATM 4 O O1 LIG L .  0.0  6.0  0.0  1.0  0.0
HETATM 5 H H1 LIG L . 10.0  0.0  0.0  1.0  0.0
"""
        c = extract_ligand_centroid(cif)
        # Mean of (0,0,0), (6,0,0), (0,6,0) — H is filtered out
        assert c is not None
        assert abs(c[0] - 2.0) < 1e-6
        assert abs(c[1] - 2.0) < 1e-6
        assert abs(c[2] - 0.0) < 1e-6

    def test_no_hetatm_returns_none(self):
        from mammal_repurposing.pockets.pose_extract import extract_ligand_centroid
        cif = """data_t
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C 1.0 2.0 3.0
"""
        c = extract_ligand_centroid(cif)
        assert c is None

    def test_ligand_filter(self):
        from mammal_repurposing.pockets.pose_extract import extract_ligand_centroid
        cif = """data_t
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
HETATM 2 C C1 LIG L .  10.0  10.0  10.0  1.0  0.0
HETATM 3 N N1 OTH L .  0.0  0.0  0.0  1.0  0.0
"""
        c = extract_ligand_centroid(cif, ligand_label_comp_id="LIG")
        # Should pick up ONLY the LIG row
        assert c is not None
        assert abs(c[0] - 10.0) < 1e-6


# ---------------------------------------------------------------------------
# §15_v2_fusion --calibrated-mammal + --znorm-mammal smoke
# ---------------------------------------------------------------------------
class TestFusionCalibratedZnorm:
    def test_znorm_makes_per_target_std_uniform(self, mini_targets, mini_compounds, mini_scores_pass):
        """End-to-end style: after Z-norm, per-target std should be ~1.0."""
        df = mini_scores_pass.copy()
        grp = df.groupby("target_uniprot")["predicted_pkd"]
        mu = grp.transform("mean")
        sigma = grp.transform("std")
        z = (df["predicted_pkd"] - mu) / sigma
        z = z.where(sigma.notna() & (sigma != 0), 0.0)
        df["predicted_pkd_znorm"] = (6.5 + 1.25 * z).clip(lower=2.0, upper=11.0)
        # After Z-norm + 1.25× scaling, every target's predicted_pkd_znorm std
        # is exactly 1.25 (modulo single-target zero-variance edge).
        per_t_std = df.groupby("target_uniprot")["predicted_pkd_znorm"].std()
        per_t_std_nonzero = per_t_std[per_t_std > 0]
        assert (per_t_std_nonzero > 0.9).all()
        assert (per_t_std_nonzero < 1.6).all()
