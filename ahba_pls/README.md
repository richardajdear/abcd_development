# ahba_pls — imaging transcriptomics of adolescent cortical thinning in ABCD

**Question.** Can the "transcriptomic signature of adolescent thinning" —
NSPN-PLS2 in Whitaker, Vértes et al. 2016 (*PNAS*, doi:10.1073/pnas.1601745113)
and AHBA-C3 in Dear et al. 2024 (*Nat Neurosci*, doi:10.1038/s41593-024-01624-4) —
be re-derived *from scratch* by imaging transcriptomics on the much larger ABCD
longitudinal sample, and does the resulting gene set carry SCZ and MDD GWAS
signal?

This directory is self-contained: code, inputs, results and figures for the
analysis live here, and this README is the working log. Upstream inputs come
from the settled 7.0 mixed-model runs in `../out/` and from the cleaned AHBA
matrices in `~/Git/AHBA_updated`.

## Hypotheses

| # | Hypothesis | Status |
|:-:|:---|:---|
| H1 | An ABCD-derived PLS component is correlated with NSPN-PLS2 and AHBA-C3 in **both** DK spatial scores and gene weights | planned |
| H2 | Its gene weights are enriched for SCZ and MDD GWAS genes (PNAS-style permutation on prioritised genes; MAGMA gene-property) | planned |
| H3 | Its gene weights are associated with the genes driving variation in global dCT in the ABCD GWAS (gene-property on ABCD `global_slope`) | **deferred** — awaits a clean ABCD GWAS; specified in `FOLLOWUP_GENETICS.md` |
| H4 | Projecting subject-level slopes onto the component gives a GWAS phenotype more heritable than `global_slope` | **deferred** — cluster work; specified in `FOLLOWUP_GENETICS.md`. Prior: `slope_projC3` h² = 0.121 ± 0.057 vs `global_slope` 0.115 ± 0.056 (hpc_v2 v3 run), so expect a modest gain at best |

## Design

**X (genes).** AHBA left-hemisphere Desikan–Killiany expression from
`AHBA_updated/outputs/expression_levels/ahba_dk_lh_native_ds{0,25,50}.csv`
(33 regions × 16,009 / 12,007 / 8,005 genes; donor-native parcellation,
abagen defaults as in Dear et al. 2024). Genes z-scored across regions.
The ds75 (4,003-gene) level is **not** used — too few genes for enrichment.

**Y (imaging).** Bilateral (LH/RH-averaged) DK maps from the settled runs:

| map | source run | term |
|:---|:---|:---|
| dCT — total thinning rate (mm/yr) | `out/thickness_dsk_70_139406217085` | age slope, no global covariate |
| CT — baseline thickness | same | intercept |
| dT1T2 — T1w/T2w rate | `out/t1t2_ratio_dsk_70_9a62dde44370` | age slope |
| T1T2 — T1w/T2w intercept | same | intercept |
| slopePC1–3 — inter-individual slope covariance maps | `docs/developmental_maps_noglobal.csv` | PCA loadings |

**Y-matrix options** (each fit at ds0 / ds25 / ds50):

1. `dCT` only — single-Y PLS; component 1 ≡ gene–map Pearson correlations (the baseline)
2. `dCT + CT` — closest analogue of the PNAS CT/dCT pair
3. `dCT + dT1T2` — two rate-of-change measures
4. `CT + dCT + T1T2 + dT1T2` — full four-variable PNAS design with T1w/T2w (a myelin proxy) replacing MT
5. `slopePC1–3` — exploratory: covariance maps rather than group means

**Method.** PLS-SVD of X'Y (PLSC). Component significance by spin permutation
of Y (DK spin geometry from `src/abcd/spatial.py`); gene weights as
region-bootstrap Z-scores (PNAS 2016 convention); components sign-aligned so
that the dCT loading is positive. **Limit to keep in view:** 33 regions, not
NSPN's 308 — power is bounded and DS-filter sensitivity matters.

**References for H1.** NSPN-PLS2 DK scores (`AHBA/data/whitakervertes2016_complete.csv`,
68 DK regions, also carries NSPN CT/dCT/MT/dMT maps) and bootstrapped gene weights
(`AHBA/data/whitakervertes2016_genes.csv`); AHBA C1–C3 gene weights
(`../data/weights.csv`, 7,973 genes) and DK scores (`../data/ahba_dme_dsk_scores.csv`).

**GWAS gene sets for H2.** SCZ: Trubetskoy 2022 prioritised / fine-mapped /
locus pool (`AHBA/data/gwas/`, and the sets used in `hpc/` MAGMA). MDD: Adams 2025.
MAGMA gene-level Z for SCZ and MDD: `../hpc/work/results/magma/{SCZ,MDD}.genes.raw` — present on this
machine but **gitignored**, so on a fresh clone they must be re-fetched from the cluster before
`08_magma_gene_property.py` will run.

## Layout

```
ahba_pls/
  README.md                 this file — plan, log, results index
  FOLLOWUP_GENETICS.md      H3/H4 specification for the cluster agent (Phase 4)
  code/                     pls.py and the numbered analysis scripts (python);
                            fig1_lead_signature.R / fig2_enrichment.R (figures, R)
  data/                     aligned X/Y inputs (checkpoints) and reference signatures
  results/                  tables (tsv/csv)
  figures/                  publication figures
  hpc/                      HPC-ready job specs (not run here)
  imaging_transcriptomics.qmd   explanatory notebook (Phase 5)
```

## Work log

- [x] **Phase 1** — inputs: Y-matrix, aligned X at 3 DS levels, NSPN-PLS2 and C3 references, GWAS gene sets (`code/01_*`–`03_*`; notes in `data/INPUTS.md`, `data/Y_MAPS.md`, `data/reference/*.md`, `data/reference/gene_sets/GENE_SETS.md`)
- [x] **Phase 2** — PLS fits (5 options × 3 DS), spin/bootstrap inference, concordance with PLS2/C3, lead signature (`code/pls.py`, `04_fit_pls.py`, `05_concordance.py`, `06_lead_signature.py`)
- [x] **Phase 3** — SCZ/MDD enrichment: permutation + MAGMA gene-property, with C3/PLS2 benchmarks (`07_permutation_enrichment.py`, `08_magma_gene_property.py`, `09_summarise_enrichment.py`; MAGMA macOS binary in `../tools/bin/magma_mac/`)
- [x] **Phase 4** — [`FOLLOWUP_GENETICS.md`](FOLLOWUP_GENETICS.md) (H3, H4 specification for the cluster agent) and exports in `hpc/` (`code/10_export_for_hpc.py`)
- [x] **Phase 5** — notebook [`imaging_transcriptomics.qmd`](imaging_transcriptomics.qmd), this results index, and a pointer row in the top-level README (`tests/test_readme.py` passes)

## Results index

### Phase 1 — inputs (facts worth knowing)

- Uncovered AHBA region: `lh_frontalpole` → 33 regions enter PLS. LH/RH agreement of the seven Y maps r = 0.96–0.995 (`results/y_map_lr_agreement.csv`).
- dCT and CT reproduce `docs/developmental_maps_noglobal.csv` to <1e-5. The T1w/T2w `rh_temporalpole` fit is singular/non-converged (kept; flagged in `data/Y_MAPS.md`).
- NSPN gene symbols are 2016-vintage: 1,310 renamed via mygene (rule: unique current symbol present in ds0, no collision) → overlap with ds0/25/50 = 14,865 / 11,175 / 7,454 of 20,710 (`data/reference/NSPN_REFERENCE.md`). 26 Excel-date-corrupted symbols dropped.
- C1–C3 scores recomputed from `weights.csv` on the current DK matrices agree with the shipped DK scores at ρ = 0.99 / 0.94 / 0.91 (C1/C2/C3); the recomputed versions are used for like-for-like comparison.
- SCZ (Trubetskoy 2022) and MDD (Adams 2025 Table S21) sets replicate the `hpc/` MAGMA definitions in symbol space; MAGMA gene-level Z for SCZ/MDD mapped to symbols at 99.3 %. Brain-expressed backgrounds: 13,763 / 10,322 / 6,931 genes.

### Phase 2 — the transcriptomic signature (H1: **supported**)

`results/pls_components.tsv`, `concordance_scores.tsv`, `concordance_weights.tsv`, `weight_stability_dCT_components.csv`; figure `figures/fig_lead_signature.png`.

| option (Y) | dCT-carrying component | cov. explained | spin p (ds0/25/50) | scores ρ vs NSPN-PLS2 | scores ρ vs C3 (recomputed) | weights ρ vs PLS2_z | weights ρ vs C3 |
|:--|:--|--:|:--|--:|--:|--:|--:|
| 1. dCT | PLS1 | 1.00 | 0.06 / 0.08 / 0.10 | 0.45 (p 0.08) | 0.24 (n.s.) | 0.60 | 0.58 (but also **0.59 with C1**) |
| 2. dCT + CT | **PLS2** (dCT salience 0.94) | 0.13 | **0.002 / 0.004 / 0.006** | **0.76 (p 0.001)** | **0.85 (p 0.001)** | **0.61** | **0.77** (C1: 0.21) |
| 3. dCT + dT1T2 | PLS2 (0.90) | 0.16 | 0.017 / 0.024 / 0.028 | 0.57 (p 0.01) | 0.87 (p 0.001) | 0.47 | 0.78 |
| 4. CT + dCT + T1T2 + dT1T2 | PLS2 (0.87) | 0.09 | 0.08 / 0.11 / 0.13 | 0.73 (p 0.001) | 0.78 (p 0.001) | 0.56 | 0.64 |
| 5. slopePC1–3 | PLS2 (= slopePC3) | 0.27 | 0.001 | 0.74 (p 0.001) | 0.49 (p 0.01) | 0.68 | 0.46 |

Signs above are reported in the *thinning* orientation (positive = higher expression where thinning is faster), which is the orientation in which NSPN-PLS2 and C3 are defined; `pls.py` itself fixes the dCT salience positive, so raw outputs have the opposite sign.

Reading:
- **dCT alone does not isolate the signature.** Its gene vector loads equally on C1 (the dominant static axis) and C3, and it is only marginal under the spin null. The thinning map carries the static transcriptional gradient because thinning rate covaries with baseline thickness / myelination.
- **Adding CT as a second Y column absorbs the static axis into PLS1** (opt2 PLS1 weights ρ = −0.92 with C1) and leaves a spin-significant PLS2 that is the C3-like thinning signature at both levels — the same structure as NSPN, where the thinning signal also appeared as the second component.
- The dCT-carrying components of options 2, 3 and 4 are one signature (pairwise weight ρ 0.83–0.94) and are stable across DS filters (weight ρ = 0.995–0.999 between adjacent filters, 0.981–0.998 between the ds0/ds50 extremes), so gene filtering is not the limiting factor here as it was for deriving C3 by PCA.
- The ABCD thinning map itself matches the NSPN thinning map (ρ = 0.64, p_spin = 0.002), so the replication is not an artefact of a shared atlas.

**Lead signature = option 2, PLS2, ds25** (`results/lead_signature_weights.tsv`, `lead_signature_scores.csv`). Top-decile overlap with C3 top decile 5.3×, with NSPN-PLS2 4.1× (hypergeometric p ≈ 0). Top genes: *AGBL4, CHRM3, CPNE8, PRSS12, LRRTM4, RAPGEF2, CCK*; bottom: *CAPN2, GIT2, PDE9A, HDAC1, MCM3, HLA-E*. Cell-class markers (Seidlitz 2020 compilation): neuronal (Ex, In) strongly positive, all glial/vascular classes negative (`lead_celltype_enrichment.tsv`; z values are descriptive — genes are not independent).

### Phase 3 — SCZ / MDD enrichment of the gene weights (H2: **supported by the continuous test, nominal in strength; the signal is a subset of C3's**)

`results/enrichment_summary.tsv` (one row per vector × disorder), `permutation_enrichment.tsv`, `magma_gene_property.tsv`, `magma_conditional.tsv`, `magma_runs/`; figure `figures/fig_enrichment.png`.

Two tests, shared brain-expressed universe (AHBA genes with a MAGMA Z in both disorders). Sign convention throughout: positive = genes weighted towards *faster* thinning carry more disorder GWAS signal.

**A. PNAS-style set permutation** (mean weight of set genes vs 10,000 size- and gene-length-matched random sets):

| weight vector | SCZ prioritised (120) | SCZ locus pool (685) | MDD high-confidence (308) | MDD pool (4,600) |
|:--|--:|--:|--:|--:|
| ABCD PLS2 ds0 / ds25 / ds50 | z 1.0 / 0.0 / −0.1 | 1.8 / 1.0 / 1.0 | **2.8 / 2.4 / 1.5** (p 0.005 / 0.015 / 0.12) | 1.4 / 1.8 / 1.6 |
| static PLS1 (control) | 1.5 | 2.3 | −2.4 | 0.9 |
| AHBA C3 | 1.1 | 2.4 | **3.4** | 1.9 |
| NSPN PLS2 | 0.2 | 1.4 | 0.1 | 0.3 |

No SCZ set is enriched in the thinning signature by this test; the MDD high-confidence set (Adams 2025) is nominally enriched, in the same direction as for C3. The static/C1-like axis carries the SCZ locus-pool signal (z ≈ 2.1–2.3) as strongly as anything, and is *depleted* for MDD high-confidence genes, so the SCZ pool result is not specific to thinning.

**B. MAGMA gene-property** (disorder gene Z regressed on the continuous weight, gene-gene LD correlations from the `.genes.raw`, MAGMA's internal gene-size/density covariates):

| weight vector | SCZ β_std (p) | MDD β_std (p) |
|:--|--:|--:|
| ABCD PLS2 ds0 | 0.037 (4e-4) | 0.026 (0.011) |
| ABCD PLS2 ds25 | 0.028 (0.022) | 0.028 (0.021) |
| ABCD PLS2 ds50 | 0.034 (0.034) | 0.041 (0.006) |
| dCT alone (opt 1) | 0.041 (0.001) | 0.031 (0.010) |
| dCT+dT1T2 PLS2 (opt 3) | 0.033 (0.008) | 0.039 (0.001) |
| static PLS1 (control) | 0.026 (0.037) | 0.015 (0.20) |
| AHBA C1 | 0.026 (0.09) | 0.008 (0.58) |
| **AHBA C3** | **0.064 (7e-5)** | **0.075 (1e-6)** |
| NSPN PLS2 | 0.033 (4e-4) | 0.015 (0.11) |

The thinning signature is positively associated with both SCZ and MDD gene-level association at every DS level and for every Y-option that isolates it — effect sizes ≈ half of C3's. Top-decile indicators are weaker (opt2 ds0 SCZ p = 0.03; others n.s.), i.e. the effect is spread along the continuum rather than concentrated in the extreme genes.

**C. Conditional models** (shared universe n = 6,672): ABCD PLS2 keeps its association when conditioned on C1 (SCZ p = 0.004, MDD 0.002) or on the static PLS1 (0.008, 0.002), so it is not the static gradient. But conditioned on **C3 it drops to zero** (SCZ β = −0.01 p = 0.62; MDD −0.03 p = 0.17), whereas C3 conditioned on ABCD PLS2 retains its signal (SCZ p = 0.006, MDD 1e-5). ABCD PLS2 and NSPN PLS2 are mutually attenuating rather than nested: in the joint model neither SCZ coefficient survives (ABCD PLS2 p = 0.36, NSPN PLS2 p = 0.085, from p = 0.014 and 0.004 alone), while for MDD only ABCD PLS2 retains its association (p = 0.015 vs 0.87).

Reading: the data-driven ABCD signature recovers the disorder-relevant part of C3, but does not add to it — C3 (derived from AHBA alone, with careful gene filtering) remains the sharper gene ranking for SCZ/MDD. This is the expected order: the ABCD component is a 33-region projection and can only capture the part of C3 that is expressed in the DK-resolution thinning map. The gain from ABCD is therefore in *validation* (an independent, 20× larger imaging sample reproduces the signature and its disorder enrichment) rather than in a sharper gene list.

### Phase 4 — exports for the cluster (H3, H4 deferred)

`FOLLOWUP_GENETICS.md` specifies both arms with the conventions of `hpc_v2`/`hpc_v3`
and states the priors to beat (`global_slope` h² = 0.115 ± 0.056, `slope_projC3` 0.121 ± 0.057;
every MAGMA test in `hpc_v3` §10 was null). Exports, all in the thinning orientation:

| file | what |
|:--|:--|
| `hpc/lead_pls2_dk_scores_68.csv` | lead-component regional scores, 68 DK labels (both frontal poles NA) — input to the H4 projection phenotype |
| `hpc/lead_pls2_gene_covar_entrez.txt` | MAGMA `--gene-covar` file, Entrez-keyed: `thinning_Z_ds0/25/50` + `AHBA_C3` for conditioning |
| `hpc/lead_pls2_gene_weights_symbol.tsv` | the same weights by gene symbol, with rank and decile |

### Parcellation control — is the weaker enrichment a DK-resolution artefact? (**no**)

`code/11_parcellation_test.py` → `results/parcellation_test.tsv`, `parcellation_test_weights.tsv`.

C3's weights were fitted on the HCP-MMP matrix (`hcp_3d.csv`, **137** left parcels ×
7,973 genes); this project's PLS ran on DK (**33** regions) — 4.2× fewer observations per
gene weight, so granularity was a live explanation for C3's ~2× larger gene-property
effect. Controlled test: the same 7,973 genes, the same PCA, parcellation the only
difference; plus an *oracle* PLS whose Y is the C3 score map itself at DK resolution
(what a 33-region PLS can recover with a perfect phenotype). One shared universe of
6,839 genes, so the βs are directly comparable.

| gene weights | regions | SCZ β_std (se) | MDD β_std (se) | ρ with shipped C3 |
|:--|--:|--:|--:|--:|
| C3 shipped (DME, HCP) | 137 | 0.061 (0.016) | 0.076 (0.016) | — |
| C3 by PCA, HCP | 137 | 0.059 (0.016) | 0.066 (0.015) | 0.92 |
| **C3 by PCA, DK** | **33** | **0.047 (0.016)** | **0.080 (0.016)** | 0.80 |
| oracle PLS (Y = C3 map), DK | 33 | 0.057 (0.016) | 0.063 (0.016) | 0.91 |
| ABCD PLS2 (Y = thinning), DK | 33 | 0.040 (0.016) | 0.041 (0.015) | 0.77 |

**Coarsening the parcellation costs almost nothing.** Re-deriving C3 at 33 regions keeps
the full effect (HCP-vs-DK difference |z| ≤ 0.6 for both disorders), and an oracle PLS at
33 regions keeps it too — so neither the parcellation nor the PLS step is the bottleneck.
What costs ~40 % of the effect is that **thinning is an imperfect proxy for the axis**:
the ABCD weights correlate ρ = 0.77 with C3, and simple regression dilution predicts
β ≈ ρ²·β_C3 = 0.036 (SCZ) / 0.044 (MDD) against 0.040 / 0.041 observed. The attenuation
is fully accounted for by that fidelity, with nothing left over for resolution.

Implication: **a finer imaging parcellation should not be expected to close the gap.** It
would be worth doing if finer parcels made the *thinning map itself* a better proxy for
C3 (DK parcels straddle the gradient), and that is cheap to check first — HCP-MMP
thickness for ABCD 5.1 already exists at `ABCD/abcd-data-release-5.1/processed/` (360
ROIs; 10,778 / 7,092 / 2,800 subjects at v0 / v2 / v4) and `src/abcd/io.py` already loads
it (`parcellation: hcp`). If the dCT–C3 score correlation rises appreciably above the
DK value (ρ = −0.55 for the rate map, 0.85 for the PLS scores), a 7.0 re-parcellation on
the cluster is justified; if it does not, the limit is the phenotype, not the atlas. Note
7.0 has no `processed/` directory, so 7.0 HCP thickness needs the surface parcellation
re-run.

### HCP-MMP arm — fitting the PLS in C3's own parcellation (**adds MDD signal, not SCZ**)

`code/12_hcp_pls.py` → `results/hcp_pls_{components,weights,scores}`, `hcp_concordance.tsv`,
`hcp_vs_dk_enrichment.tsv`, `hcp_vs_dk_conditional.tsv`; figure `figures/fig_hcp_vs_dk.png`.

The parcellation control above said a coarse atlas costs nothing *when the phenotype is
already the axis*, and left open the one case that would matter: whether finer parcels make
the thinning map a better proxy. HCP-MMP thickness for 7.0 now exists (commit 3961821,
`configs/ct_70_hcp_noglobal_mv2.yaml`), so this tests it directly.

**The run.** `thickness_hcp_70_fec93121f0dd` — the settled specification on HCP-MMP.
Assemble keeps **5,947 subjects** with ≥2 QC-passing visits (DK: 8,192); 360 regions fit in
2 min, 4/360 singular, 4 non-converged, 0 errors. The PLS uses the 137 left parcels with
AHBA donor coverage, bilateral (lh/rh mean), with the medial-wall parcel `H` dropped.

**The component is the same axis, now measured in its native space.** PLS2 carries dCT
(salience 0.99), explains 14% of the cross-covariance, spin p = 0.012 (5,000 rotations of
the complete 180-parcel map), bootstrap reproducibility 0.90. Scores vs C3
ρ = 0.69 (p_spin < 0.001, 137 parcels), gene weights ρ = 0.67; weights vs the DK signature
ρ = 0.73. The DK score correlation of 0.85 was over 33 regions — 0.69 over 137 is the more
honest number, not a loss of signal.

**Enrichment** (one shared MAGMA universe of 6,063 genes; the DK comparator is
refitted on `dk_3d.csv`, the same abagen build and the same 7,973 genes as the HCP matrix, so
parcellation is not confounded with pipeline):

| gene weights | parcels | SCZ β_std (p) | MDD β_std (p) |
|:--|--:|--:|--:|
| **HCP-MMP PLS2** | **137** | 0.044 (0.012) | **0.070** (2e-05) |
| DK PLS2, matched genes | 33 | 0.046 (0.007) | 0.055 (0.0007) |
| HCP-MMP, dCT alone | 137 | 0.053 (0.0022) | 0.079 (2e-06) |
| AHBA C3 | 137 | 0.068 (9e-05) | 0.082 (9e-07) |

Because the two weight vectors are correlated (ρ = 0.73), the comparison is made inside
MAGMA rather than by differencing βs:

| model | SCZ β (p) | MDD β (p) |
|:--|--:|--:|
| HCP PLS2 \| DK PLS2 | 0.020 (0.45) | **0.065 (0.0092)** |
| DK PLS2 \| HCP PLS2 | 0.031 (0.23) | 0.007 (0.78) |
| HCP PLS2 \| AHBA C3 | -0.003 (0.89) | 0.029 (0.19) |
| AHBA C3 \| HCP PLS2 | 0.071 (0.0026) | 0.062 (0.0055) |

**Reading.** For **MDD** the finer parcellation adds real signal: β rises 0.055 → 0.070, HCP
survives conditioning on DK (p = 0.009) and DK does not survive conditioning on HCP
(p = 0.78) — the DK version is a degraded copy of the HCP one. For **SCZ** there is no gain;
the two attenuate each other and neither dominates. Both remain absorbed by C3, which keeps
its association conditioned on either. So the HCP arm improves the *ABCD-derived* ranking for
MDD without changing the headline conclusion that C3 is the sharper ranking.

Also notable: in HCP space **dCT alone** (option 1) is the strongest ABCD vector
(MDD 0.079, SCZ 0.053). In DK its enrichment was already competitive (SCZ 0.041, MDD 0.031,
vs 0.037/0.026 for the lead at ds0); what disqualified it there was not the enrichment but
that it fails to *isolate* the axis — its weights load equally on C1 (ρ = −0.59) and C3
(−0.58) — and that its component is marginal under the spin null. The same caution applies
here: spin p = 0.27, so the map-level covariance with expression is not spatially specific.
Treat the single-Y option as a gene-level result without a spatial claim.

**Caveats.** (i) 5,947 subjects, not 8,192 — 8,001 sessions are missing from the parcellated
table (3,435 never reached by the array job, 5,439 `mri_surf2surf` stubs), so this is not a
like-for-like sample comparison with the DK run. (ii) Release vintage, checked against the
tables this repo reads: the DK thickness table has **4,086** six-year sessions against 7,612
six-year FreeSurfer sessions, i.e. it is 6.0-sized; the covariate table `ab_g_dyn` has 5,056
six-year rows. Of the **2,646** sessions the HCP table has that DK lacks (2,639 of them
six-year), only **12** have an age row and 10 a QC row — so the extra scans HCP uniquely
offers are unusable until the 7.0 tabulated release lands. A re-run with the missing
parcellations plus 7.0 covariates is the version to trust.

## Reproducing

Analysis (python, env `ahba-pls`): `code/01_*` → `code/12_*` in order (`11_` is the
parcellation control, `12_` the HCP-MMP arm; `12_` needs the
`thickness_hcp_70_*` run — `ABCD_CONFIG=ct_70_hcp_noglobal_mv2 python -m abcd.assemble`
then `Rscript R/fit_lmm.R --run-dir out/thickness_hcp_70_* --cores 8`). Figures (R, env
`ahba-pls-r`): `Rscript code/fig1_lead_signature.R`, `fig2_enrichment.R`, `fig3_hcp_vs_dk.R`.
The full command list with timings is at the end of `imaging_transcriptomics.qmd`.
`.gitignore` here excludes the regenerable intermediates (aligned X matrices, full
per-option weight tables, spin nulls, MAGMA run directories).

## Decisions and caveats

- 2026-09-12: ds75 dropped (user). H3 deferred until a clean ABCD GWAS exists; H4 deferred to a cluster run (user).
- All PLS runs on the 33 AHBA-covered bilateral DK regions; the missing region is recorded in `data/`.
- 2026-09-12: figures are built in R (`ggplot2` + `patchwork` + `ggseg` 2.2.1 from CRAN, env `ahba-pls-r`), not matplotlib. ggseg's R atlas renders frontal and temporal pole correctly, unlike the python polygon set. The R scripts fit nothing and hardcode nothing: every statistic printed on a figure is read from the `results/` table it plots. DK polygons are cached to `data/dk_polygons.csv`.
- Figure conventions: one-line interpretation per panel subtitle, overall interpretation as the figure subtitle, methods as short bullets in the caption.
