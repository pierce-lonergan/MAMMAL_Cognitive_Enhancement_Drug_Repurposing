"""Download the real LINCS L1000 Phase-2 (GSE70138) Level-5 data.

Reproducibility helper: fetches the three GSE70138 supplementary files that the
V8 pipeline reads (the real-data Gate 1 in scripts/92 and the chemCPA training in
scripts/73) into data/cache/lincs/. The data itself is NOT committed (it is
~5 GB); this script is the committed, idempotent way to obtain it from GEO so a
fresh clone can reproduce the V8 results.

Files (from GEO series GSE70138 supplementary):
  GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz  (~5 GB)
  GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz                     (~2 MB)
  GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz                    (~80 KB)
The pert-info file is saved locally as GSE70138_compoundinfo.txt.gz, the name the
pipeline reads (it carries pert_id / pert_iname / canonical_smiles / inchi_key).

The script is idempotent: it skips any file already present, and decompresses the
Level-5 GCTX only if the .gctx is missing. GSE92742 (Phase 1, ~7-8 GB) is a
larger optional add-on and is not fetched here.

Usage:
  python scripts/download_lincs.py            # download anything missing
  python scripts/download_lincs.py --check    # report presence only, no download
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LINCS = ROOT / "data" / "cache" / "lincs"
GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl"
GEO_LANDING = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138"

# (GEO supplementary filename, local filename). Local name may alias the GEO name
# where the pipeline expects a shorter name.
FILES = [
    ("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz",
     "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"),
    ("GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
     "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"),
    ("GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
     "GSE70138_compoundinfo.txt.gz"),
]
GCTX_GZ = "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
GCTX = GCTX_GZ[:-3]  # strip .gz


def _human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}GB"


def _download(url: str, dest: Path) -> bool:
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url) as r:  # nosec - fixed NCBI host
            total = int(r.headers.get("Content-Length", 0))
            done = 0
            with open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
                    done += len(chunk)
                    if total:
                        pct = 100 * done / total
                        print(f"\r  {dest.name}: {_human(done)}/{_human(total)} "
                              f"({pct:.0f}%)", end="", flush=True)
        print()
        tmp.rename(dest)
        return True
    except Exception as e:  # noqa: BLE001
        print(f"\n  ERROR downloading {url}: {e}")
        if tmp.exists():
            tmp.unlink()
        return False


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="report presence only; do not download")
    args = ap.parse_args(argv)
    LINCS.mkdir(parents=True, exist_ok=True)

    print(f"LINCS cache: {LINCS}")
    gctx = LINCS / GCTX
    have_gctx = gctx.exists()
    missing = []
    for geo_name, local_name in FILES:
        dest = LINCS / local_name
        if dest.exists():
            print(f"  present: {local_name} ({_human(dest.stat().st_size)})")
        elif local_name == GCTX_GZ and have_gctx:
            # the .gz archive is only needed to produce the .gctx; once the
            # decompressed .gctx exists the .gz is redundant (don't re-download 5 GB)
            print(f"  present (decompressed): {GCTX}; .gz not needed")
        else:
            print(f"  MISSING: {local_name}")
            missing.append((geo_name, dest))
    print(f"  decompressed GCTX: {'present' if have_gctx else 'MISSING'}")

    if args.check:
        ok = not missing and have_gctx
        print("\nall present" if ok else
              f"\n{len(missing)} archive(s) missing"
              f"{'' if have_gctx else ' + GCTX not decompressed'}")
        return 0 if ok else 1

    for geo_name, dest in missing:
        url = f"{GEO_SUPPL}/{geo_name}"
        print(f"\ndownloading {geo_name} ...")
        if not _download(url, dest):
            print(f"  could not fetch automatically. Get it manually from the GEO "
                  f"supplementary listing:\n    {GEO_LANDING}")
            return 2

    # decompress the Level-5 GCTX if needed (the pipeline reads the .gctx)
    gz = LINCS / GCTX_GZ
    if not have_gctx and gz.exists():
        print(f"\ndecompressing {GCTX_GZ} (~5 GB, takes a minute) ...")
        with gzip.open(gz, "rb") as fi, open(gctx, "wb") as fo:
            shutil.copyfileobj(fi, fo, length=1 << 24)
        print(f"  wrote {gctx.name} ({_human(gctx.stat().st_size)})")

    print("\nLINCS ready. Reproduce V8 with:\n"
          "  python scripts/92_v8_real_gate1.py   # real-data Gate 1\n"
          "  python scripts/73_chemcpa_real_lincs_training.py --scale cognition")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
