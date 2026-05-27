# Nootropic-Similarity Annotation v1 (§8.10)

Per-compound max Tanimoto similarity to each of the canonical nootropic chemotypes. Uses ECFP4 / Morgan-2 / 2048 bits — same fingerprint as the §A.4 Tanimoto-to-actives ranker, so scores are directly comparable across the pipeline.

Canonical set (14 compounds): aniracetam, atomoxetine, bupropion, d-amphetamine, donepezil, fluoxetine, galantamine, memantine, methylphenidate, modafinil, piracetam, pitolisant, rivastigmine, rolipram

## Tag thresholds

- **novel_scaffold**: T_max < 0.30 to every canonical nootropic — structurally distinct from the existing field; IP-novel candidate
- **intermediate**: 0.30 ≤ T_max ≤ 0.85 — scaffold-related but not a direct analog
- **analog**: T_max > 0.85 — essentially a structural analog and likely patent-encumbered
- **unknown**: SMILES failed to parse

## Pipeline-wide tag distribution

- **novel_scaffold**: 279 (93.6%)
- **intermediate**: 19 (6.4%)

## Top 50 candidates by RRF, with nootropic-similarity annotation

| # | Compound | Tier | RRF | Nearest nootropic | T | Novelty tag |
|---|---|---|---|---|---|---|
| 1 | d-amphetamine | positive_control | 0.830 | modafinil | 0.25 | `novel_scaffold` |
| 2 | methylphenidate | positive_control | 0.805 | modafinil | 0.24 | `novel_scaffold` |
| 3 | bupropion | extended_cns | 0.802 | rivastigmine | 0.24 | `novel_scaffold` |
| 4 | aniracetam | named_in_research | 0.766 | piracetam | 0.30 | `intermediate` |
| 5 | rasagiline | extended_cns | 0.737 | methylphenidate | 0.18 | `novel_scaffold` |
| 6 | pridopidine | named_in_research | 0.730 | donepezil | 0.22 | `novel_scaffold` |
| 7 | levetiracetam | extended_cns | 0.729 | piracetam | 0.44 | `intermediate` |
| 8 | lisdexamfetamine | extended_cns | 0.727 | d-amphetamine | 0.40 | `intermediate` |
| 9 | lanicemine | named_in_research | 0.714 | d-amphetamine | 0.32 | `intermediate` |
| 10 | pramiracetam | extended_cns | 0.712 | piracetam | 0.45 | `intermediate` |
| 11 | cx-717 | named_in_research | 0.711 | aniracetam | 0.23 | `novel_scaffold` |
| 12 | rivastigmine | positive_control | 0.709 | bupropion | 0.24 | `novel_scaffold` |
| 13 | selegiline | extended_cns | 0.709 | d-amphetamine | 0.42 | `intermediate` |
| 14 | cx-516 | named_in_research | 0.702 | aniracetam | 0.26 | `novel_scaffold` |
| 15 | guanfacine | extended_cns | 0.695 | modafinil | 0.23 | `novel_scaffold` |
| 16 | piracetam | extended_cns | 0.685 | aniracetam | 0.30 | `intermediate` |
| 17 | modafinil | positive_control | 0.684 | d-amphetamine | 0.25 | `novel_scaffold` |
| 18 | atomoxetine | positive_control | 0.680 | fluoxetine | 0.55 | `intermediate` |
| 19 | donepezil | positive_control | 0.678 | rolipram | 0.24 | `novel_scaffold` |
| 20 | chembl1255723 | chembl_expanded | 0.668 | methylphenidate | 0.17 | `novel_scaffold` |
| 21 | pramipexole | extended_cns | 0.667 | donepezil | 0.09 | `novel_scaffold` |
| 22 | chembl1256414 | chembl_expanded | 0.665 | methylphenidate | 0.15 | `novel_scaffold` |
| 23 | cep-26401 | named_in_research | 0.665 | pitolisant | 0.21 | `novel_scaffold` |
| 24 | chembl302231 | chembl_expanded | 0.663 | d-amphetamine | 0.10 | `novel_scaffold` |
| 25 | chembl1256378 | chembl_expanded | 0.663 | methylphenidate | 0.16 | `novel_scaffold` |
| 26 | atenolol | negative_control | 0.658 | modafinil | 0.21 | `novel_scaffold` |
| 27 | isrib | named_in_research | 0.657 | pitolisant | 0.23 | `novel_scaffold` |
| 28 | levodopa | extended_cns | 0.656 | d-amphetamine | 0.28 | `novel_scaffold` |
| 29 | enalapril | negative_control | 0.656 | methylphenidate | 0.28 | `novel_scaffold` |
| 30 | rolipram | positive_control | 0.655 | donepezil | 0.23 | `novel_scaffold` |
| 31 | (r,s)-ampa | chembl_expanded | 0.654 | d-amphetamine | 0.20 | `novel_scaffold` |
| 32 | chembl4228464 | chembl_expanded | 0.652 | donepezil | 0.19 | `novel_scaffold` |
| 33 | chembl608151 | chembl_expanded | 0.652 | rolipram | 0.14 | `novel_scaffold` |
| 34 | pitolisant | positive_control | 0.652 | donepezil | 0.18 | `novel_scaffold` |
| 35 | ibuprofen | negative_control | 0.652 | d-amphetamine | 0.33 | `intermediate` |
| 36 | chembl91184 | chembl_expanded | 0.651 | methylphenidate | 0.19 | `novel_scaffold` |
| 37 | clomipramine | extended_cns | 0.650 | pitolisant | 0.20 | `novel_scaffold` |
| 38 | risperidone | extended_cns | 0.650 | pitolisant | 0.19 | `novel_scaffold` |
| 39 | chembl292558 | chembl_expanded | 0.648 | d-amphetamine | 0.10 | `novel_scaffold` |
| 40 | chembl292924 | chembl_expanded | 0.648 | d-amphetamine | 0.10 | `novel_scaffold` |
| 41 | (s)-ampa | chembl_expanded | 0.648 | d-amphetamine | 0.20 | `novel_scaffold` |
| 42 | clonidine | extended_cns | 0.647 | bupropion | 0.15 | `novel_scaffold` |
| 43 | dizocilpine | chembl_expanded | 0.645 | methylphenidate | 0.18 | `novel_scaffold` |
| 44 | oxiracetam | extended_cns | 0.644 | piracetam | 0.47 | `intermediate` |
| 45 | propranolol | extended_cns | 0.644 | atomoxetine | 0.24 | `novel_scaffold` |
| 46 | riluzole | named_in_research | 0.641 | fluoxetine | 0.16 | `novel_scaffold` |
| 47 | memantine | named_in_research | 0.640 | methylphenidate | 0.08 | `novel_scaffold` |
| 48 | clemastine | named_in_research | 0.635 | pitolisant | 0.30 | `novel_scaffold` |
| 49 | ropinirole | extended_cns | 0.634 | rivastigmine | 0.16 | `novel_scaffold` |
| 50 | danavorexton | named_in_research | 0.633 | methylphenidate | 0.26 | `novel_scaffold` |

## Novel-scaffold candidates (T_max < 0.30, n=279)

These compounds have no close structural neighbour in the canonical nootropic set — IP-novel chemotypes that the pipeline surfaced via mechanism-driven evidence alone.

| Compound | RRF | Nearest nootropic | T | Tier |
|---|---|---|---|---|
| d-amphetamine | 0.830 | modafinil | 0.25 | positive_control |
| methylphenidate | 0.805 | modafinil | 0.24 | positive_control |
| bupropion | 0.802 | rivastigmine | 0.24 | extended_cns |
| rasagiline | 0.737 | methylphenidate | 0.18 | extended_cns |
| pridopidine | 0.730 | donepezil | 0.22 | named_in_research |
| cx-717 | 0.711 | aniracetam | 0.23 | named_in_research |
| rivastigmine | 0.709 | bupropion | 0.24 | positive_control |
| cx-516 | 0.702 | aniracetam | 0.26 | named_in_research |
| guanfacine | 0.695 | modafinil | 0.23 | extended_cns |
| modafinil | 0.684 | d-amphetamine | 0.25 | positive_control |
| donepezil | 0.678 | rolipram | 0.24 | positive_control |
| chembl1255723 | 0.668 | methylphenidate | 0.17 | chembl_expanded |
| pramipexole | 0.667 | donepezil | 0.09 | extended_cns |
| chembl1256414 | 0.665 | methylphenidate | 0.15 | chembl_expanded |
| cep-26401 | 0.665 | pitolisant | 0.21 | named_in_research |
| chembl302231 | 0.663 | d-amphetamine | 0.10 | chembl_expanded |
| chembl1256378 | 0.663 | methylphenidate | 0.16 | chembl_expanded |
| atenolol | 0.658 | modafinil | 0.21 | negative_control |
| isrib | 0.657 | pitolisant | 0.23 | named_in_research |
| levodopa | 0.656 | d-amphetamine | 0.28 | extended_cns |
| enalapril | 0.656 | methylphenidate | 0.28 | negative_control |
| rolipram | 0.655 | donepezil | 0.23 | positive_control |
| (r,s)-ampa | 0.654 | d-amphetamine | 0.20 | chembl_expanded |
| chembl4228464 | 0.652 | donepezil | 0.19 | chembl_expanded |
| chembl608151 | 0.652 | rolipram | 0.14 | chembl_expanded |
| pitolisant | 0.652 | donepezil | 0.18 | positive_control |
| chembl91184 | 0.651 | methylphenidate | 0.19 | chembl_expanded |
| clomipramine | 0.650 | pitolisant | 0.20 | extended_cns |
| risperidone | 0.650 | pitolisant | 0.19 | extended_cns |
| chembl292558 | 0.648 | d-amphetamine | 0.10 | chembl_expanded |
| chembl292924 | 0.648 | d-amphetamine | 0.10 | chembl_expanded |
| (s)-ampa | 0.648 | d-amphetamine | 0.20 | chembl_expanded |
| clonidine | 0.647 | bupropion | 0.15 | extended_cns |
| dizocilpine | 0.645 | methylphenidate | 0.18 | chembl_expanded |
| propranolol | 0.644 | atomoxetine | 0.24 | extended_cns |
| riluzole | 0.641 | fluoxetine | 0.16 | named_in_research |
| memantine | 0.640 | methylphenidate | 0.08 | named_in_research |
| clemastine | 0.635 | pitolisant | 0.30 | named_in_research |
| ropinirole | 0.634 | rivastigmine | 0.16 | extended_cns |
| danavorexton | 0.633 | methylphenidate | 0.26 | named_in_research |

---

Generated by `scripts/37_v5_nootropic_similarity.py`.