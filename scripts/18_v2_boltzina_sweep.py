"""V2 Phase 0.4 — Boltz-2 affinity for top-N surviving compounds × 22 targets.

Reads ADMET gates (PASS + selected FLAG), takes the top-N MAMMAL hits per
compound to limit the call count, then runs Boltz-2 affinity for every
surviving (target, compound) pair.

Wall-clock budget: 5-30 s per pair on RTX 5070 → ~1500 pairs × 15 s ≈ 6 hours.
Incremental flush every 25 pairs; `--resume` mandatory.

Output: data/results/v2/boltzina_affinity.parquet
"""

from __future__ import annotations

import argparse
import datetime as dt
import logging
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.cluster_a.boltzina import free_gpu, score_affinity  # noqa: E402
from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    DTI_SCORES_PARQUET,
    RESULTS_DIR,
    TARGETS_PARQUET,
    ensure_dirs,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("boltzina_sweep")

V2_DIR = RESULTS_DIR / "v2"
DEFAULT_OUT = V2_DIR / "boltzina_affinity.parquet"
ADMET_GATES = V2_DIR / "admet_gates.parquet"


def _build_pair_list(
    *,
    targets: pd.DataFrame,
    compounds: pd.DataFrame,
    gates: pd.DataFrame,
    mammal_scores: pd.DataFrame | None,
    include_flag: bool,
    top_n_per_compound: int | None,
) -> list[tuple[str, str, str, str]]:
    """Build the (target_uniprot, sequence, compound_name, smiles) call list.

    Surviving compounds = gate_status in {PASS} ∪ (FLAG if include_flag).
    Optionally restrict to top-N MAMMAL targets per surviving compound (most
    of the structural signal is concentrated where MAMMAL already predicts
    high affinity).
    """
    allowed = {"PASS"} | ({"FLAG"} if include_flag else set())
    surviving = gates[gates["gate_status"].isin(allowed)]
    surviving_names = set(surviving["compound_name"].str.lower().str.strip())
    logger.info("Surviving ADMET-PASS+FLAG: %d compounds.", len(surviving_names))

    cmp_lookup = compounds.set_index(
        compounds["name"].str.lower().str.strip()
    )[["smiles"]].to_dict()["smiles"]
    tgt_lookup = targets.set_index("uniprot")[["sequence"]].to_dict()["sequence"]

    pairs: list[tuple[str, str, str, str]] = []

    if mammal_scores is not None and top_n_per_compound:
        # For each surviving compound, pick its top-N targets by MAMMAL pKd
        df = mammal_scores[
            mammal_scores["compound_name"].str.lower().str.strip().isin(surviving_names)
        ].copy()
        df["lc_name"] = df["compound_name"].str.lower().str.strip()
        topn = (
            df.sort_values(["lc_name", "predicted_pkd"], ascending=[True, False])
              .groupby("lc_name", group_keys=False)
              .head(top_n_per_compound)
        )
        for _, row in topn.iterrows():
            tgt = row["target_uniprot"]
            seq = tgt_lookup.get(tgt)
            name = row["compound_name"]
            smi = cmp_lookup.get(name.lower().strip())
            if seq and smi:
                pairs.append((tgt, seq, name, smi))
    else:
        # Full cross-product
        for tgt, seq in tgt_lookup.items():
            for lc_name, smi in cmp_lookup.items():
                if lc_name not in surviving_names:
                    continue
                # Recover original casing of name from compounds df
                orig = compounds[compounds["name"].str.lower().str.strip() == lc_name]
                if orig.empty:
                    continue
                pairs.append((tgt, seq, orig.iloc[0]["name"], smi))

    return pairs


def _filter_resume(
    pairs: list[tuple[str, str, str, str]],
    out_path: Path,
) -> tuple[list[tuple[str, str, str, str]], pd.DataFrame]:
    if not out_path.exists():
        return pairs, pd.DataFrame()
    existing = pd.read_parquet(out_path)
    done = set(zip(
        existing["target_uniprot"].astype(str),
        existing["compound_name"].astype(str).str.lower().str.strip(),
    ))
    remaining = [p for p in pairs if (p[0], p[2].lower().strip()) not in done]
    logger.info("Resume: skipping %d already-scored pairs; %d to go.",
                len(pairs) - len(remaining), len(remaining))
    return remaining, existing


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--gates", type=Path, default=ADMET_GATES)
    parser.add_argument("--mammal-scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--top-n-per-compound", type=int, default=10,
                        help="Per surviving compound, score against its top-N MAMMAL targets")
    parser.add_argument("--include-flag", action="store_true", default=True,
                        help="Include ADMET-FLAG compounds in the sweep (default true)")
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--flush-every", type=int, default=25)
    parser.add_argument("--limit", type=int, default=None,
                        help="Smoke test: process only first N pairs")
    parser.add_argument("--device", choices=["cuda", "cpu"], default="cuda")
    parser.add_argument("--no-msa", action="store_true",
                        help="Disable MSA server for faster but lower-quality predictions")
    parser.add_argument("--targets-filter", nargs="+", default=None,
                        help="Only score these UniProt IDs (e.g. P36544 Q08499 Q9Y5N1)")
    parser.add_argument("--include-allosteric-panel", action="store_true",
                        help="Force-include allosteric audit ligands at their cognate target "
                             "even if not in top-N MAMMAL hits")
    args = parser.parse_args()

    ensure_dirs()
    V2_DIR.mkdir(parents=True, exist_ok=True)

    targets = pd.read_parquet(args.targets)
    compounds = pd.read_parquet(args.compounds)
    gates = pd.read_parquet(args.gates)
    mammal_scores = pd.read_parquet(args.mammal_scores) if args.mammal_scores.exists() else None
    logger.info("Loaded: %d targets, %d compounds, %d gated compounds, mammal=%s",
                len(targets), len(compounds), len(gates),
                "yes" if mammal_scores is not None else "no")

    # --targets-filter restricts the panel to a subset (e.g. for focused gate tests)
    if args.targets_filter:
        targets = targets[targets["uniprot"].isin(args.targets_filter)].copy()
        logger.info("Targets filter applied: %d remaining (%s)",
                    len(targets), list(targets["uniprot"]))

    pairs = _build_pair_list(
        targets=targets, compounds=compounds, gates=gates,
        mammal_scores=mammal_scores,
        include_flag=args.include_flag,
        top_n_per_compound=args.top_n_per_compound,
    )

    # Force-include allosteric audit ligands at their cognate target.
    if args.include_allosteric_panel:
        allosteric_extra = {
            "P36544": ["galantamine", "encenicline", "tc-5619"],            # CHRNA7
            "Q08499": ["bpn14770", "zatolmilast", "rolipram"],              # PDE4D
            "P42261": ["tulrampator", "cx-516", "cx-717", "aniracetam"],    # GRIA1
            "P42262": ["tulrampator", "cx-516", "cx-717", "aniracetam"],
            "P42263": ["tulrampator", "cx-516", "cx-717", "aniracetam"],
            "P48058": ["tulrampator", "cx-516", "cx-717", "aniracetam"],
            "Q9Y5N1": ["pitolisant", "cep-26401"],                          # HRH3
            "P22303": ["donepezil", "rivastigmine", "galantamine"],         # ACHE
            "Q01959": ["methylphenidate", "modafinil", "d-amphetamine"],    # DAT
            "P23975": ["atomoxetine", "methylphenidate"],                   # NET
            "P21728": ["aripiprazole"],                                     # DRD1
        }
        cmp_lookup = compounds.set_index(
            compounds["name"].str.lower().str.strip()
        )[["smiles", "name"]].apply(tuple, axis=1).to_dict()
        tgt_lookup = targets.set_index("uniprot")[["sequence"]].to_dict()["sequence"]
        existing_keys = {(tgt, name.lower().strip()) for tgt, _, name, _ in pairs}
        extra_added = 0
        for tgt, comp_names in allosteric_extra.items():
            seq = tgt_lookup.get(tgt)
            if seq is None:
                continue
            for name in comp_names:
                lc = name.lower().strip()
                if (tgt, lc) in existing_keys:
                    continue
                entry = cmp_lookup.get(lc)
                if entry is None:
                    continue
                smi, orig_name = entry
                pairs.append((tgt, seq, orig_name, smi))
                extra_added += 1
        logger.info("Allosteric panel: added %d extra (target,compound) pairs.", extra_added)

    logger.info("Built call list: %d (target, compound) pairs.", len(pairs))

    if args.resume:
        pairs, existing = _filter_resume(pairs, args.out)
    else:
        existing = pd.DataFrame()

    if args.limit:
        pairs = pairs[:args.limit]
        logger.info("Limit applied: only %d pairs", len(pairs))

    if not pairs:
        logger.info("Nothing to do; output is up to date.")
        return 0

    buffer: list[dict] = []
    completed: list[pd.DataFrame] = [existing] if not existing.empty else []
    written = 0
    last_target: str | None = None
    t0 = time.perf_counter()

    try:
        for i, (tgt, seq, name, smi) in enumerate(pairs):
            if last_target and tgt != last_target:
                free_gpu()
            pair_start = time.perf_counter()
            try:
                r = score_affinity(
                    target_uniprot=tgt, sequence=seq,
                    compound_name=name, smiles=smi,
                    device=args.device, use_cache=True,
                    use_msa_server=not args.no_msa,
                )
                buffer.append(r.as_dict())
            except Exception as e:
                logger.warning("Boltz failed (%s, %s): %s; recording NaN.", tgt, name, e)
                buffer.append({
                    "target_uniprot": tgt,
                    "compound_name": name,
                    "smiles": smi,
                    "affinity_pred_value": float("nan"),
                    "affinity_probability_binary": float("nan"),
                    "pose_plddt": None,
                    "mode": "full",
                    "scored_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                })
            pair_elapsed = time.perf_counter() - pair_start
            done = i + 1
            avg = (time.perf_counter() - t0) / done
            remaining_secs = avg * (len(pairs) - done)
            logger.info("[%d/%d] %s + %s -> %.1fs (avg %.1fs, ETA %.0f min)",
                        done, len(pairs), tgt, name, pair_elapsed, avg, remaining_secs / 60)
            last_target = tgt

            if (done % args.flush_every) == 0:
                _flush(buffer, completed, args.out)
                written += len(buffer)
                buffer.clear()
    finally:
        if buffer:
            _flush(buffer, completed, args.out)
            written += len(buffer)

    logger.info("Wrote %d new Boltz affinity rows; total %s rows in %s",
                written, "?", args.out)
    return 0


def _flush(buffer: list[dict], completed: list[pd.DataFrame], out_path: Path) -> None:
    if not buffer:
        return
    new_df = pd.DataFrame(buffer)
    completed.append(new_df)
    combined = pd.concat(completed, ignore_index=True, sort=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    combined.to_parquet(tmp, index=False)
    tmp.replace(out_path)
    completed.clear()
    completed.append(combined)
    logger.info("Flushed: %d new rows; total in file %d.", len(new_df), len(combined))


if __name__ == "__main__":
    raise SystemExit(main())
