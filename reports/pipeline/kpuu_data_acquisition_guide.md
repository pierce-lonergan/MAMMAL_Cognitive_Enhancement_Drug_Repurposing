# How to obtain Kp,uu,brain data + the licenses - explicit directions

The Stage-3 free-exposure regressor (`engine/free_exposure.py`) currently trains on logBB (B3DB,
CC0) as an HONEST proxy. The principled target is unbound brain-to-plasma Kp,uu,brain. This is the
single highest-value upgrade. Below are the concrete, ranked routes to legitimately get the data,
the licenses each needs, and how it drops into the existing pipeline. Verify the specifics (prices,
exact license terms, current URLs) before relying on them - terms change.

## The one legal idea that unblocks most of this

Raw numeric facts (a compound's measured Kp,uu value) are generally NOT copyrightable in the US;
what is protected is the specific compiled TABLE / its creative selection-and-arrangement, and in
the EU a sui-generis database right can apply. Practical consequence:
- TRAINING a model on extracted numeric values from a published paper for internal research is
  low-risk and routine.
- REDISTRIBUTING the dataset in this open-source repo needs an actual license (CC-BY / CC0) OR
  written permission. So the cheap routes below split into "openly licensed -> can vendor in the
  repo" vs "extract-for-training-only / request-permission-to-redistribute".
- When in doubt, confirm with your institution's library/legal; this guide is not legal advice.

## Route 1 - openly licensed data you can vendor directly (fastest, free)

1a. **Open-access (CC-BY) Kp,uu papers - extract the SI table, redistribute with attribution.**
   - Loryan et al. 2022, "Unbound Brain-to-Plasma Partition Coefficient, Kp,uu,brain - a Game
     Changing Parameter for CNS Drug Discovery and Development," Pharm Res (OA, PMC9246790,
     DOI 10.1007/s11095-022-03246-6). Compiles/reviews Kp,uu values; check its tables + SI license.
   - Heliyon 2024, "Accurate prediction of Kp,uu,brain based on experimental measurement of
     Kp,brain ..." (OA, PMC10828645, ScienceDirect S2405844024003359). Heliyon is CC-BY; its
     experimental Kp,brain / Kp,uu table is reusable WITH attribution.
   - Frontiers Drug Discov 2024, "Application of machine learning to predict unbound drug
     bioavailability in the brain" (CC-BY). Check for a deposited dataset.
   STEP: open the article -> Supporting Information / Data Availability -> download the XLSX/CSV ->
   confirm the article license is CC-BY or CC0 (stated on the page) -> if yes, vendor it into
   `data/raw/` with a provenance header (DOI + license) and cite it.

1b. **Polaris / Biogen ADME-Fang (open) for the EFFLUX + binding COMPONENTS of Kp,uu.**
   - `polarishub.io/datasets/biogen/adme-fang-v1` (Fang et al. 2023, J Med Chem). Openly licensed;
     includes MDR1-MDCK efflux ratio + human/rat plasma protein binding - the two dominant drivers
     of Kp,uu. Not Kp,uu itself, but lets you build a stronger efflux/binding-aware proxy now,
     dependency-free. `pip install polaris-lib`; `po.load_dataset("biogen/adme-fang-v1")`.

1c. **CarbonAI / public BBB sets** already partly used (B3DB CC0). These are logBB/BBB+/-, not
   Kp,uu - keep as the proxy floor, not the target.

## Route 2 - permission-gated published data (free, needs a request)

2a. **Friden et al. 2009, J Med Chem 52:6233-6243** - the seminal ~41-drug Kp,uu set (ACS,
   paywalled SI). Two legitimate paths:
   - Extract the numeric values for internal model training (facts; low risk) and cite; do NOT
     re-host the ACS table verbatim in the public repo without permission.
   - To REDISTRIBUTE: request reuse via ACS RightsLink (article page -> "Rights & Permissions" ->
     select reuse type) - academic/non-commercial reuse of a data table is routinely granted, often
     free. See pubs.acs.org/page/copyright/permissions.html.
2b. **Email the corresponding authors for the compiled dataset.** The Uppsala group
   (Hammarlund-Udenaes / Irena Loryan / Markus Friden) authored most of the public Kp,uu work and
   often share curated tables for academic research. Template:

   > Subject: Request to reuse Kp,uu,brain dataset for an open-source CNS-penetration model
   > Dear Prof. ___, I am building an open-source (non-commercial) machine-learning model of
   > unbound brain exposure for cognitive-drug repurposing. Would you permit reuse of the
   > Kp,uu,brain values from [paper, DOI], and if possible share a machine-readable table (SMILES +
   > Kp,uu)? I will cite the source and, with your permission, redistribute under CC-BY. Happy to
   > sign a data-use agreement. Thank you - [name, affiliation].

2c. **A signed Data Use Agreement (DUA) / Material(Data) Transfer Agreement** is the formal
   instrument if an academic group or company shares non-public data - route it through your
   institution's research-contracts office.

## Route 3 - generate / buy the data (definitive; you then OWN the license)

The cleanest way to get exactly the compounds you care about, with full rights to use and
redistribute: COMMISSION the measurement at a CRO. You own the resulting data outright.
- The assay (what to order): Kp,uu = Kp x (fu,brain / fu,plasma), where
  - Kp = total brain:plasma AUC ratio from a rodent PK study (or steady-state), and
  - fu,brain = unbound fraction by brain-homogenate EQUILIBRIUM DIALYSIS, fu,plasma by plasma
    equilibrium dialysis. (For efflux confirmation, add an MDR1-MDCK or Pgp-KO mouse Kp,uu.)
- Vendors that offer the brain-homogenate fu / Kp,uu service (request a quote; confirm current
  scope/price): **Pharmaron** (DMPK "Binding / Partitioning" e-store lists rat/mouse brain
  homogenate fu + Kp,uu), **Cyprotex / Evotec**, **WuXi AppTec**, **Sai Life Sciences**,
  **Charles River**. Indicative scale: brain-homogenate fu is cheap per compound; a full in-vivo
  Kp,uu (animals + bioanalysis) is the costly part - budget per compound, so prioritise a focused
  panel (the psychoplastogen / candidate set) rather than thousands.
- A no-animal first cut: order MDR1-MDCK efflux ratio + brain/plasma fu in vitro and combine with a
  literature Kp - cheaper, gives an efflux-aware estimate to expand the training set.

## Route 4 - commercial / consortium databases (paid license)

- **DrugBank** - free academic license by application (drugbank.com), paid commercial license via
  OMx; has some PK but Kp,uu coverage is thin - not a primary Kp,uu source.
- **Certara / Simcyp**, **Schrodinger** - PBPK + brain models and internally-curated Kp,uu, but the
  underlying data are proprietary (the Schrodinger JCIM 2023 physics-based Kp,uu model,
  DOI 10.1021/acs.jcim.3c00150, used internal data not shared). License = vendor contract.

## Recommended sequence (cheapest first)

1. NOW (free, this week): pull the OA tables (Loryan 2022, Heliyon 2024) + the Polaris ADME-Fang
   efflux/PPB set; confirm each license; vendor the CC-BY/CC0 ones into `data/raw/` with provenance.
   This alone may give a few-hundred-row efflux-aware set to retrain Stage-3 against a real
   Kp,uu/efflux target rather than logBB.
2. PARALLEL (free, ~weeks): RightsLink request + author email for Friden 2009 and any other
   reuse-restricted Kp,uu tables; on grant, add them under CC-BY.
3. IF the model still under-performs on YOUR chemical space: commission a focused CRO Kp,uu panel
   (Route 3) on the psychoplastogen / candidate compounds - data you fully own and can redistribute.

## How it drops into the existing pipeline (no architecture change)

`engine/free_exposure.py` is already built target-agnostic: same RDKit featurizer + LightGBM +
Mondrian split-conformal + kNN applicability domain. To switch from logBB to Kp,uu you only swap the
TRAINING TABLE:
1. write the obtained data to `data/raw/kpuu_train.csv` (columns: `smiles, kpuu` [or `log_kpuu`],
   `source`, `license`),
2. point `scripts/111_train_logbb.py` (or a `scripts/119_train_kpuu.py` copy) at it and set the
   PASS threshold to the Kp,uu CNS cutoff (Kp,uu >= 0.3 is the field convention for adequate free
   exposure),
3. retrain; the conformal band + abstain rule carry the (larger) Kp,uu uncertainty honestly.
The efflux feature is already pluggable (rule / ADMET-AI Pgp), so an efflux-measured Kp,uu set
slots straight in.

INTEGRITY: never paste a Kp,uu value into the repo without its real source + license recorded; do
not fabricate values to fill gaps - abstain (the conformal AD already does this) on compounds with
no measured/again-licensed data.

Sources: Loryan 2022 (DOI 10.1007/s11095-022-03246-6, PMC9246790); Heliyon 2024 (PMC10828645,
S2405844024003359); Friden 2009 (J Med Chem 52:6233-6243); Biogen ADME-Fang (polarishub.io/
datasets/biogen/adme-fang-v1; Fang 2023); ACS permissions (pubs.acs.org/page/copyright/
permissions.html); Pharmaron DMPK binding/partitioning service; Schrodinger Kp,uu (DOI 10.1021/
acs.jcim.3c00150).
