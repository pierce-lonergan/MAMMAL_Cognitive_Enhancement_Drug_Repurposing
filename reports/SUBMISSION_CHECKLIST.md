# Submission checklist — OSF pre-registration + bioRxiv preprint

**Author:** Pierce Lonergan · ORCID [0009-0008-4235-396X](https://orcid.org/0009-0008-4235-396X)
**Both accounts already created and ORCID-linked.**

The two artifacts are written and committed:
- **Manuscript:** `reports/manuscript_class_prognostic_biorxiv.md`
- **Pre-registration:** `reports/osf_preregistration_class_prognostic.md`
- **Flagship figure (Fig. 1):** `figures/flagship/thesis_synthesis.png` (300-ready at 200 dpi; regenerate at higher dpi with `scripts/83_flagship_figure.py` if needed)

---

## Step 1 — OSF pre-registration (do this FIRST, before the preprint goes public)

1. osf.io → **Create new Preregistration** → template **"OSF Preregistration"** (matches the doc's section structure).
2. Paste each section from `osf_preregistration_class_prognostic.md` into the matching field. The doc is already organised as Study information / Design / Sampling / Variables & analysis / etc.
3. Attach the frozen repository link and the manuscript as supplementary files.
4. **Register** (this timestamps and locks it — required so §6's prospective predictions count as pre-registered). Mint the DOI.
5. Confirm ORCID 0009-0008-4235-396X is listed as contributor.
6. Copy the registration DOI back into the manuscript's "Data and code availability" and the repo `PROJECT_STATUS.md`.

## Step 2 — bioRxiv preprint

1. **The PDF is already generated**: `reports/manuscript_class_prognostic.pdf` (typeset via `pandoc manuscript_class_prognostic_biorxiv.md -o manuscript_class_prognostic.pdf --pdf-engine=xelatex --resource-path=.:.. -V geometry:margin=1in -V mainfont="DejaVu Serif" -V colorlinks=true -V header-includes='\usepackage{xurl}'` (the font flag is required for the Greek ρ/θ glyphs; `colorlinks` + `xurl` make the ORCID/repo links clickable and wrap long URLs so the header no longer overflows) — the figure is embedded and Unicode renders correctly). Regenerate with that command after any manuscript edit. Current PDF is 9 pages (methods doubled in review round 2).
2. biorxiv.org → **Submit** → category **Neuroscience** (secondary: **Pharmacology and Toxicology**).
3. Title, author (Pierce Lonergan), ORCID (auto-linked), abstract (paste the Abstract section).
4. Upload manuscript PDF + Fig. 1 (`thesis_synthesis.png`) as a separate figure file.
5. License: **CC-BY-4.0** (matches the repo). Declare no conflicts, no funding.
6. Link the OSF pre-registration DOI in the manuscript and in the bioRxiv "Supplementary/related" field.

## Step 3 — after posting

- Add the bioRxiv DOI to `README.md`, `PROJECT_STATUS.md`, and `CITATIONS.bib`.
- Update the repo citation block author/ORCID (currently lists a model co-author placeholder — set to Pierce Lonergan + ORCID).

---

## Before submission — verify the bibliography

The manuscript reference list is author/journal/year only (deliberately, to avoid asserting volume/page/DOI from memory). **Reconcile each of the 8 references against `CITATIONS.bib` and the source DOI before posting** — bioRxiv does not enforce reference accuracy, but reviewers will. The robustness supplement and figure scripts have no such caveat (every number recomputes from data).

## Pre-flight quality gate (all currently ✅)

- [x] Every headline number recomputes from committed data (`scripts/83` regenerates Fig. 1 from source).
- [x] Leakage audit stated per predictor; permutation p + bootstrap CI reported.
- [x] Honest-scope / limitations section present (n=31, Roberts ceiling, no wet-lab, LOCO ceiling = 0).
- [x] The popularity confound is disclosed, not hidden.
- [x] Pre-registration cleanly separates completed (exploratory) from prospective (committed) analyses.
- [x] 503 automated tests pass.

## Notes / honest caveats for the cover letter

- This is a single-author, compute-modest (one consumer GPU) study; its strength is the *leakage discipline* and the *negative-and-useful* central result, not scale.
- The retrospective AUROC of 1.00 should be presented as a *contrast* (class ≫ affinity/genetics) driven by mechanism-class outcome-homogeneity, **not** as a predictive oracle — the manuscript and figure already frame it this way. Reviewers will probe this; the framing is defensible because the leave-one-class-out ceiling (0.00) is reported alongside.
