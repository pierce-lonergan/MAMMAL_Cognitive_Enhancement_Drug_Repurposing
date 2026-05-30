"""One-off: replace em-dashes with context-appropriate punctuation in the
public-facing docs (manuscript, cover letter, OSF pre-reg, bioRxiv note).

Rules (then a human review fixes any remaining awkward spot):
  0. A table-cell em-dash  ` | — | `  ->  ` | n/a | `  (it meant "not computed").
  1. A PAIRED aside  ` — aside — ` :
       - short aside (<= 35 chars, no internal comma)  ->  `, aside, `  (commas)
       - longer / comma-containing aside               ->  ` (aside) `  (parens)
  2. A remaining SINGLE ` — `  ->  `: `  (the "statement: elaboration" pattern that
     dominates this writing).
  3. Light fixes: a colon immediately before a conjunction/relativiser that joins a
     clause ( so / but / which / whereas ...) reads better as a comma.
En-dashes (numeric ranges 0.2-0.5, 2014-2022, M1-M4, 150(2-3)) are left untouched:
correct typography, not an AI tell.

Usage: python scripts/util_deemdash.py
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = [
    "reports/manuscript_class_prognostic_biorxiv.md",
    "reports/cover_letter_journal.md",
    "reports/osf_preregistration_class_prognostic.md",
    "reports/biorxiv_submission_note.md",
]
EM = "—"


def _pair_sub(m: re.Match) -> str:
    aside = m.group(1).strip()
    if len(aside) <= 35 and "," not in aside:
        return f", {aside}, "
    return f" ({aside}) "


def convert(text: str) -> tuple[str, int]:
    before = text.count(EM)
    # 0. table cell
    text = re.sub(rf"\|\s*{EM}\s*\|", "| n/a |", text)
    # 1. paired aside (non-greedy, up to 180 chars, no embedded em-dash)
    text = re.sub(rf"{EM}\s*([^{EM}\n]{{1,180}}?)\s*{EM}", _pair_sub, text)
    # 2. remaining single em-dash -> colon
    text = re.sub(rf"\s*{EM}\s*", ": ", text)
    # 3. colon before a clause-joining word reads better as a comma
    text = re.sub(r":\s+(so |but |yet |whereas |which |confirming |reflecting )",
                  lambda m: ", " + m.group(1), text)
    # tidy parentheses spacing and stray doubles
    text = text.replace("( ", "(").replace(" )", ")")
    text = re.sub(r",\s*,", ", ", text)
    text = re.sub(r"[ \t]{2,}", " ", text)
    return text, before


def main() -> int:
    total = 0
    for rel in FILES:
        p = ROOT / rel
        new, n = convert(p.read_text(encoding="utf-8"))
        p.write_text(new, encoding="utf-8")
        total += n
        print(f"{rel}: {n} converted, {new.count(EM)} remaining")
    print(f"total converted: {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
