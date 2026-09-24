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
locus pool (`AHBA/data/gwas/`, and the sets used in `legacy/hpc/` MAGMA). MDD: Adams 2025.
MAGMA gene-level Z for SCZ and MDD: `../genetic_analysis/inputs/magma/{SCZ,MDD}.genes.raw` — present on this
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
- [x] **Phase 4** — [`FOLLOWUP_GENETICS.md`](FOLLOWUP_GENETICS.md) and exports (`code/10_export_for_hpc.py`); **H3 run on the cluster and null**, H4 still open
- [x] **Phase 5** — notebook [`imaging_transcriptomics.qmd`](imaging_transcriptomics.qmd), this results index, and a pointer row in the top-level README (`tests/test_readme.py` passes)

## Results index

### Phase 1 — inputs (facts worth knowing)

- Uncovered AHBA region: `lh_frontalpole` → 33 regions enter PLS. LH/RH agreement of the seven Y maps r = 0.962–0.995 (`results/y_map_lr_agreement.csv`).
- The Y maps come from the 7.0-tabulated runs: **8,716 children / 26,949 sessions** (DK thickness, all 68 hemisphere-region models converged, none singular) and the T1w/T2w run on the same sample, where `rh_temporalpole` is singular/non-converged (kept; flagged in `data/Y_MAPS.md`).
- dCT and CT reproduce `docs/developmental_maps_noglobal.csv` to 4.9e-08 and 4.8e-06 (ρ = 1.000 for both).
- NSPN gene symbols are 2016-vintage: 1,310 renamed via mygene (rule: unique current symbol present in ds0, no collision) → overlap with ds0/25/50 = 14,865 / 11,175 / 7,454 of 20,710 (`data/reference/NSPN_REFERENCE.md`). 26 Excel-date-corrupted symbols dropped.
- C1–C3 scores recomputed from `weights.csv` on the current DK matrices agree with the shipped DK scores at ρ = 0.99 / 0.94 / 0.91 (C1/C2/C3); the recomputed versions are used for like-for-like comparison.
- SCZ (Trubetskoy 2022) and MDD (Adams 2025 Table S21) sets replicate the `legacy/hpc/` MAGMA definitions in symbol space; MAGMA gene-level Z for SCZ/MDD mapped to symbols at 99.3 %. Brain-expressed backgrounds: 13,763 / 10,322 / 6,931 genes.

### Phase 2 — the transcriptomic signature (H1: **supported**)

`results/pls_components.tsv`, `concordance_scores.tsv`, `concordance_weights.tsv`, `weight_stability_dCT_components.csv`; figure `figures/fig_lead_signature.png`.

| option (Y) | dCT-carrying component | cov. explained | spin p (ds0/25/50) | scores ρ vs NSPN-PLS2 | scores ρ vs C3 (recomputed) | weights ρ vs PLS2_z | weights ρ vs C3 |
|:--|:--|--:|:--|--:|--:|--:|--:|
| 1. dCT | PLS1 | 1.00 | 0.065 / 0.082 / 0.108 | 0.49 (p 0.048) | 0.28 (p 0.185) | 0.61 | 0.59 (C1: 0.56) |
| 2. dCT + CT | **PLS2** (dCT salience 0.94) | 0.13 | **0.002 / 0.002 / 0.004** | **0.75 (p 0.001)** | **0.83 (p 0.001)** | **0.62** | **0.76** (C1: -0.20) |
| 3. dCT + dT1T2 | PLS2 (dCT salience 0.90) | 0.17 | 0.015 / 0.021 / 0.026 | 0.57 (p 0.010) | 0.87 (p 0.001) | 0.48 | 0.78 (C1: -0.34) |
| 4. CT + dCT + T1T2 + dT1T2 | PLS2 (dCT salience 0.86) | 0.09 | 0.078 / 0.098 / 0.119 | 0.76 (p 0.001) | 0.78 (p 0.001) | 0.57 | 0.62 (C1: -0.07) |
| 5. slopePC1–3 | PLS2 | 0.27 | 0.001 / 0.000 / 0.000 | 0.74 (p 0.001) | 0.52 (p 0.007) | 0.68 | 0.48 (C1: -0.09) |

All ρ and saliences at ds25 unless the column says otherwise. Signs are reported in the *thinning*
orientation (positive = higher expression where thinning is faster), which is the orientation in which
NSPN-PLS2 and C3 are defined; `pls.py` itself fixes the dCT salience positive, so raw outputs carry the
opposite sign.

Reading:
- **dCT alone does not isolate the signature.** Its gene vector loads about equally on C3 (0.59) and on
  C1, the dominant static axis (0.56), and it is only marginal under the spin null
  (p = 0.065 / 0.082 / 0.108). The thinning map carries the static transcriptional gradient because
  thinning rate covaries with baseline thickness and myelination.
- **Adding CT as a second Y column absorbs the static axis into PLS1** (opt2 PLS1 weights ρ = 0.92 with C1)
  and leaves a spin-significant PLS2 that is the C3-like thinning signature at every DS level — the same
  structure as NSPN, where the thinning signal also appeared as the second component.
- The dCT-carrying components of options 2, 3 and 4 are one signature (pairwise weight ρ = 0.83–0.94
  at matched DS) and each is near-invariant to the DS filter (ρ = 0.980–0.999 across ds0/ds25/ds50), so gene
  filtering is not the limiting factor here as it was for deriving C3 by PCA.
- The ABCD thinning map itself matches the NSPN thinning map (ρ = 0.66, p_spin = 0.002), so the
  replication is not an artefact of a shared atlas.

**Lead signature = option 2, PLS2, ds25** (`results/lead_signature_weights.tsv`, `lead_signature_scores.csv`).
Top-decile overlap with the C3 top decile 5.2×, with NSPN-PLS2 4.2× (hypergeometric p ≈ 0). Top genes:
*AGBL4*, *PRSS12*, *LKAAEAR1*, *CHRM3*, *SEC24D*, *WIPF2*, *CPNE8*; bottom: *CAPN2*, *PHACTR2*, *ITGAV*, *HLA-B*, *BRD2*, *HLA-C*. Cell-class markers: neuronal classes strongly positive, glial and
vascular classes negative (`lead_celltype_enrichment.tsv`, and the eight-ranking comparison below).

### Phase 3 — SCZ / MDD enrichment of the gene weights (H2: **supported by the continuous test, nominal in strength; the signal is a subset of C3's**)

`results/enrichment_summary.tsv` (one row per vector × disorder), `permutation_enrichment.tsv`, `magma_gene_property.tsv`, `magma_conditional.tsv`, `magma_runs/`; figure `figures/fig_enrichment.png`.

Two tests, shared brain-expressed universe (AHBA genes with a MAGMA Z in both disorders). Sign convention throughout: positive = genes weighted towards *faster* thinning carry more disorder GWAS signal.

**A. PNAS-style set permutation** (mean weight of set genes vs 10,000 size- and gene-length-matched random sets):

| weight vector | SCZ prioritised (92) | SCZ locus pool (394) | MDD high-confidence (202) | MDD pool (1626) |
|:--|--:|--:|--:|--:|
| ABCD PLS2 ds0 / ds25 / ds50 | z 1.0 / -0.0 / -0.1 | 1.8 / 1.0 / 1.0 | **2.8 / 2.4 / 1.5** (p 0.005 / 0.016 / 0.132) | 1.4 / 1.7 / 1.5 |
| static PLS1 (control) | 1.5 | 2.3 | -2.4 | 0.9 |
| AHBA C3 | 1.1 | 2.4 | **3.4** | 1.9 |
| NSPN PLS2 | 0.2 | 1.4 | 0.1 | 0.3 |

No SCZ set is enriched in the thinning signature by this test; the MDD high-confidence set (Adams 2025)
is nominally enriched at ds0 and ds25, in the same direction as for C3. The static/C1-like axis carries
the SCZ locus-pool signal (z = 2.3) as strongly as anything, and is *depleted* for MDD
high-confidence genes (-2.4), so the SCZ pool result is not specific to thinning.

**B. MAGMA gene-property** (disorder gene Z regressed on the continuous weight, gene-gene LD correlations from the `.genes.raw`, MAGMA's internal gene-size/density covariates):

| weight vector | SCZ β_std (p) | MDD β_std (p) |
|:--|--:|--:|
| ABCD PLS2 ds0 | 0.037 (5e-4) | 0.025 (0.014) |
| ABCD PLS2 ds25 | 0.027 (0.026) | 0.026 (0.027) |
| ABCD PLS2 ds50 | 0.034 (0.036) | 0.039 (0.01) |
| dCT alone (opt 1) | 0.040 (0.001) | 0.030 (0.012) |
| dCT+dT1T2 PLS2 (opt 3) | 0.032 (0.011) | 0.038 (0.002) |
| static PLS1 (control) | 0.026 (0.036) | 0.015 (0.208) |
| AHBA C1 | 0.026 (0.093) | 0.008 (0.585) |
| **AHBA C3** | **0.064 (7e-5)** | **0.075 (1e-6)** |
| NSPN PLS2 | 0.033 (4e-4) | 0.015 (0.105) |

The thinning signature is positively associated with both SCZ and MDD gene-level association at every DS
level and for every Y-option that isolates it — effect sizes ≈ half of C3's. Top-decile indicators are
weaker (opt2 ds0 SCZ p = 0.027; others n.s.), i.e. the effect is spread along the continuum rather
than concentrated in the extreme genes.

**C. Conditional models** (shared universe n = 6,672 SCZ / 6,662 MDD): ABCD PLS2 keeps its association
when conditioned on C1 (SCZ p = 0.004, MDD 0.003) or on the static PLS1 (0.009, 0.003), so it is
not the static gradient. But conditioned on **C3 it drops to zero** (SCZ β = -0.012 p = 0.64;
MDD -0.035 p = 0.13), whereas C3 conditioned on ABCD PLS2 retains its signal
(SCZ β = 0.069 p = 0.0056; MDD 0.106 p = 7e-6). ABCD PLS2 and NSPN PLS2 are mutually
attenuating rather than nested: in the joint model neither SCZ coefficient survives (ABCD PLS2
p = 0.376, NSPN PLS2 p = 0.086, from p = 0.014 and 0.0042 alone), while for MDD only
ABCD PLS2 retains its association (p = 0.024 vs 0.93).

Reading: the data-driven ABCD signature recovers the disorder-relevant part of C3, but does not add to it
— C3 (derived from AHBA alone, with careful gene filtering) remains the sharper gene ranking for SCZ/MDD.
This is the expected order: the ABCD component is a 33-region projection and can only capture the part of
C3 that is expressed in the DK-resolution thinning map. The gain from ABCD is therefore in *validation*
(an independent imaging sample roughly 30× the size of NSPN's — 8,716 children against ~300 —
reproduces the signature and its disorder enrichment) rather
than in a sharper gene list.

### Phase 4 — exports for the cluster, and H3's answer (**null**)

`FOLLOWUP_GENETICS.md` specifies both arms with the conventions of `legacy/hpc_v2`/`hpc_v3`. Exports,
all in the thinning orientation:

| file | what |
|:--|:--|
| `hpc/lead_pls2_dk_scores_68.csv` | lead-component regional scores, 68 DK labels (both frontal poles NA) — input to the H4 projection phenotype |
| `hpc/lead_pls2_gene_covar_entrez.txt` | MAGMA `--gene-covar` file, Entrez-keyed: `thinning_Z_ds0/25/50` + `AHBA_C3` for conditioning |
| `hpc/lead_pls2_gene_weights_symbol.tsv` | the same weights by gene symbol, with rank and decile |

**H3 — do the signature's genes drive the ABCD's own thinning GWAS? No.** The cluster agent ran it as
step 8 of the 7.0-tabulated genetics (`genetic_analysis/step8_ahba_pls_h3.sh`,
`work/results_70tab/magma_ahba_pls_h3/table_h3.tsv`): each phenotype's EUR gene-level Z regressed on
`thinning_Z_ds0/25/50` and on C3, marginally and conditioned on C3, across five phenotypes — 35 tests,
**none below p = 0.05**. For `global_slope` the marginal βs are -0.002 / +0.008 / +0.013
(p 0.84 / 0.35 / 0.23); conditioned on C3 they rise to +0.029–+0.033 (p 0.061–0.072) — the same
near-threshold signal three times over, since the DS columns are near-duplicates. That is the expected
result rather than evidence against the signature: `global_slope` is heritable in GREML
(h² = 0.166 ± 0.045, n = 6,011) but has no GWAS hit and an LDSC h² z of only 0.57–0.60, so there is
no detectable common-variant signal for the gene weights to line up with. **H4** (a projection phenotype with
higher heritability) remains open and is the informative follow-up.

The cluster agent then repeated the whole pipeline on the HCP-MMP export (commit `3e80410`, `work/results_70tab_hcp/`). H3 is null there too — 50 rows, smallest p = 0.066; for `global_slope` the marginal βs are +0.009–+0.016 (p 0.081–0.215), and conditioning on C3 does not help. Worth noting for H4: on HCP the thinning GWAS's LDSC h² z rises to 1.78 from 0.57 on DK, so a parcellation-aware phenotype is the more promising route than a better gene weighting.

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
the full effect (HCP-vs-DK difference |z| ≤ 0.63 for both disorders), and an oracle PLS at
33 regions keeps it too — so neither the parcellation nor the PLS step is the bottleneck.
What costs ~40 % of the effect is that **thinning is an imperfect proxy for the axis**:
the ABCD weights correlate ρ = 0.77 with C3, and simple regression dilution predicts
β ≈ ρ²·β_C3 = 0.036 (SCZ) / 0.044 (MDD) against 0.040 / 0.041 observed. The attenuation
is fully accounted for by that fidelity, with nothing left over for resolution.

Implication: **a finer imaging parcellation should not be expected to close the gap** by giving the
weights more observations. It would only help if finer parcels made the *thinning map itself* a better
proxy for the axis (DK parcels straddle the gradient) — which is exactly what the next section tests,
and for MDD that is what happens.

### HCP-MMP arm — fitting the PLS in C3's own parcellation (**adds MDD signal, not SCZ**)

`code/12_hcp_pls.py` → `results/hcp_pls_{components,weights,scores}`, `hcp_concordance.tsv`,
`hcp_vs_dk_enrichment.tsv`, `hcp_vs_dk_conditional.tsv`; figure `figures/fig_hcp_vs_dk.png`.

The parcellation control above said a coarse atlas costs nothing *when the phenotype is
already the axis*, and left open the one case that would matter: whether finer parcels make
the thinning map a better proxy. This tests it directly.

**The run.** `thickness_hcp_70_aa6e91efba82` — the settled specification on HCP-MMP, refitted
after the 2026-09-14 parcellation backfill. It keeps **8,716 children** with ≥2 QC-passing
visits and 26,946 sessions, i.e. **exactly the DK sample** (8,716); all 358 hemisphere-region
models converged with none singular (the earlier 6,537-child run had 4 singular and 4
non-converged of 360). Parcel `H` (FreeSurfer medial wall) is excluded by the run config.
The PLS uses the 137 left parcels with AHBA donor coverage, bilateral (lh/rh mean).

**The component is the same axis, now measured in its native space.** PLS2 carries dCT
(salience 0.99), explains 14% of the cross-covariance, spin p = 0.011 (5,000 rotations of
the complete 180-parcel map), bootstrap reproducibility 0.90. Scores vs C3
ρ = 0.70 (p_spin = 0.0004, 137 parcels), gene weights ρ = 0.68; weights vs the DK
signature ρ = 0.73. The DK score correlation of 0.83 was over 33 regions — 0.70 over
137 is the more honest number, not a loss of signal.

**Enrichment** (one shared MAGMA universe of 6,063 SCZ / 6,051 MDD genes; the DK comparator is
refitted on `dk_3d.csv`, the same abagen build and the same 7,973 genes as the HCP matrix, so
parcellation is not confounded with pipeline):

| gene weights | parcels | SCZ β_std (p) | MDD β_std (p) |
|:--|--:|--:|--:|
| **HCP-MMP PLS2 (dCT+CT)** | **137** | 0.045 (0.01) | **0.069** (3e-5) |
| DK PLS2 (dCT+CT) | 33 | 0.046 (0.008) | 0.053 (0.001) |
| HCP-MMP, dCT alone | 137 | 0.054 (0.002) | 0.078 (3e-6) |
| DK, dCT alone | 33 | 0.065 (2e-4) | 0.059 (4e-4) |
| AHBA C3 | 137 | 0.068 (9e-5) | 0.082 (9e-7) |
| AHBA C1 (static) | 137 | 0.037 (0.029) | 0.010 (0.551) |

All six rows are on the **matched** gene basis (see *Gene bases* above): the same 7,973 genes taken
from `hcp_3d.csv` or `dk_3d.csv`, so an HCP row differs from its DK counterpart by parcellation
alone. The cross-build AHBA_updated DK vectors are in the table under `basis = "AHBA_updated
native-DK"` and are not plotted.

Because the two weight vectors are correlated (ρ = 0.73), the comparison is made inside
MAGMA rather than by differencing βs:

| model | SCZ β (p) | MDD β (p) |
|:--|--:|--:|
| HCP PLS2 \| DK PLS2 | 0.023 (0.39) | **0.067 (0.008)** |
| DK PLS2 \| HCP PLS2 | 0.029 (0.27) | 0.003 (0.91) |
| HCP PLS2 \| AHBA C3 | -0.002 (0.94) | 0.027 (0.23) |
| AHBA C3 \| HCP PLS2 | 0.070 (0.0032) | 0.063 (0.0047) |

**Reading.** For **MDD** the finer parcellation adds real signal: β rises 0.053 → 0.069, HCP
survives conditioning on DK (p = 0.008) and DK does not survive conditioning on HCP
(p = 0.91) — the DK version is a degraded copy of the HCP one. For **SCZ** there is no gain;
the two attenuate each other and neither dominates (0.39 and 0.27 jointly). Both remain absorbed
by C3, which keeps its association conditioned on either. So the HCP arm improves the
*ABCD-derived* ranking for MDD without changing the headline conclusion that C3 is the sharper
ranking.

Also notable: in HCP space **dCT alone** (option 1) is the strongest ABCD vector
(MDD 0.078, SCZ 0.054). What disqualifies it is not the enrichment but that it fails to
*isolate* the axis — in DK its weights load about equally on C1 and C3 — and that its component
is marginal under the spin null here too (p = 0.27), so the map-level covariance with expression
is not spatially specific. Treat the single-Y option as a gene-level result without a spatial claim.

**Sample caveat retired.** The earlier version of this section ran on 6,537 children because the
parcellated table covered only 24,921 of 33,825 FreeSurfer sessions and the covariate tables were
6.0-vintage. Both are fixed: coverage is 33,795/33,825 after `tools/hcp_backfill.sbatch`, the 7.0
tabulation supplies age and QC for the extra sessions, and the HCP and DK arms now run on the same
8,716 children. **Every conclusion above survived the change unchanged** — MDD β 0.069 against
0.070 on the smaller sample, the joint-model verdicts identical — which is the strongest available
evidence that the MDD gain is about the atlas and not about who is in the sample.

### Two summary figures

Built from single-universe re-runs so every number on them is comparable:
`code/14_magma_all_options.py` (one MAGMA covar file with all nine rankings →
`results/magma_all_marginal.tsv`, `magma_all_joint.tsv`) and
`code/15_celltype_all_options.py` (→ `results/celltype_all_options.tsv`).

| figure | script | content |
|:--|:--|:--|
| `figures/fig_signature_both.png` | `code/fig1_signature_both.R` | the signature derived in **both** parcellations: brain maps (DK row, HCP-MMP row), six small concordance panels, and the cell-class profile of three rankings |
| `figures/fig_enrichment_combined.png` | `code/fig2_enrichment_combined.R` | MAGMA only: every ranking's SCZ/MDD β alone, then the seven head-to-head joint models |
| `figures/fig_celltypes.png` | `code/fig4_celltypes.R` | cell-class marker enrichment of all nine rankings, plus the astrocyte/oligodendrocyte plane |
| **`figures/fig_signature_designs.png`** | `code/fig1_signature_designs.R` | **the version to read**: brain maps in both parcellations, then two square pair matrices (regional maps, gene vectors) with DK below the diagonal and HCP-MMP above |
| `figures/fig_signature_dct_only.png` | `code/fig1_signature_dct_only.R` | the dCT-alone design on its own, with its cell-class panel |

HCP-MMP brain rendering needs `ggsegGlasser`, which is not on CRAN for this R
version; it is installed from GitHub into `ahba_pls/.Rlib` (gitignored) —
`remotes::install_github("ggseg/ggsegGlasser", lib = "ahba_pls/.Rlib")`. The
polygons are cached to `data/hcp_polygons.csv`, so the figures render without it.

**Single-universe enrichment** (n = 6,672 SCZ / 6,662 MDD genes) reproduces the separate runs:
all ABCD-derived rankings are enriched for both disorders, NSPN-PLS2 reaches SCZ only
(MDD p = 0.10) and C1 neither. The atlas pairing uses the gene-matched DK vector (see
*Gene bases*): jointly, in this run the HCP ranking keeps its MDD association (0.060, p = 0.012) while
the matched DK ranking does not (0.007, p = 0.76); for SCZ neither survives (0.33 and 0.22). The
HCP ranking also beats NSPN-PLS2 for MDD, survives C1 for both disorders, and loses to C3.

### All three Y designs in one figure — and the four-feature design

`code/17_design_grid.py` → `results/design_grid_{scores,weights,components,concordance}`;
figure `figures/fig_signature_designs.png` (`code/fig1_signature_designs.R`). This is the version of
Figure 1 to read: the cell-class panel has moved to its own figure, and in its place are two square
pair matrices covering **every** pairwise comparison among the regional maps and among the gene
vectors (see *NSPN PLS2 in HCP-MMP space* below). The four-feature design is computed but not
plotted — its component is not spin-significant and it exists in DK only — and its numbers stay in
the table below and in `design_grid_*.tsv`.

The four-feature design (CT + dCT + T1w/T2w + dT1w/T2w) is the closest available analogue of the
NSPN PNAS design. Its C3-aligned component is **selected in the script**, not assumed:
|rho| with C3 weights: PLS2 0.62 > PLS4 0.54. It is DK-only — ABCD tabulates T1w/T2w in Desikan space and the
locally-derived HCP-MMP parcellation covers thickness only, so a 137-parcel T1w/T2w map needs a new
surface run.

| design | component | spin p | cov. expl. | scores ρ vs C3 (p_spin) | scores ρ vs NSPN PLS2 | weights ρ vs C3 | weights ρ vs **C1** | weights ρ vs NSPN_z |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|
| dCT + CT, DK | PLS2 | 0.002 | 0.13 | 0.83 (0.001) | 0.75 | 0.76 | -0.20 | 0.62 |
| dCT alone, DK | PLS1 | 0.082 | 1.00 | 0.28 (0.185) | 0.49 | 0.59 | **+0.56** | 0.61 |
| **four features, DK** | PLS2 | 0.098 | 0.09 | **0.78** (0.001) | **0.76** | 0.62 | **-0.07** | 0.57 |
| dCT + CT, HCP | PLS2 | 0.011 | 0.14 | 0.70 (0.0004) | — | 0.68 | -0.11 | 0.35 |
| dCT alone, HCP | PLS1 | 0.266 | 1.00 | 0.55 (0.0004) | — | 0.61 | +0.28 | 0.37 |

Reading, in the order the figure supports it:

- **What the design buys is separation from C1, not agreement with C3.** All three DK designs reach
  weights ρ 0.59–0.76 with C3, but their C1 loadings run -0.07 (four features) and
  -0.20 (dCT + CT) against +0.56 for dCT alone. Any static map in Y does the job; the
  four-feature design, with two static and two rate maps, separates best.
- **The four-feature design agrees with NSPN PLS2 best of all** (0.76 on scores, 0.57 on weights),
  which is what you would expect since it is the closest analogue of the PNAS Y matrix — two
  thickness-like and two myelin-proxy columns, with T1w/T2w standing in for MT.
- **But its own component is not spin-significant** (p = 0.098), like the single-Y fit (0.082) and
  unlike the two-Y fits (0.002 DK, 0.011 HCP). It explains only 9 % of the cross-covariance, so the
  spin null on the singular value is unforgiving. The lead signature stays the two-Y fit; the
  four-feature component is the better *NSPN analogue* and a usable gene ranking, not a replacement.

### NSPN PLS2 in HCP-MMP space, and every pair in one figure

`code/18_nspn_to_hcp.py` → `data/reference/nspn_hcp_maps.csv`, `nspn_hcp_coverage.csv`,
`NSPN_HCP.md`; pairs from `code/17_design_grid.py` → `design_grid_pairs_{scores,weights}.tsv`;
figure `figures/fig_signature_designs.png` (`code/fig1_signature_designs.R`).

**The conversion.** The NSPN maps are published on the 308-region subdivision of DK, and both
that parcellation and HCP-MMP1 ship fsaverage annot files, so the transfer is parcels → vertices →
parcels — the route in `~/Git/AHBA/notebooks/MT_whitakervertes.ipynb`. Reimplemented here with
nibabel so the 308 values are matched to annot parcels **by name**
(`<region>_part<N>` against hemi / region / n_sub_regions) instead of by row position, asserted to
be a bijection. 357 of 360 HCP parcels keep ≥50 % vertex coverage
(179 regions after bilateral averaging); `L_H`, `R_H` and `R_Pir` fall below it.

Two things worth knowing about the source table and the result:

- **The published 308-region table uses −99 as a missing-value sentinel** (two PLS2 cells:
  left lateraloccipital part 8, right parahippocampal part 1). Left in place it survives the vertex
  averaging and destroys any Pearson correlation while leaving Spearman almost intact.
- **`lh.aparc.annot` prefixes its parcel names (`lh_bankssts`) and `rh.aparc.annot` does not**
  (`bankssts`). The first version of the DK validation grouped on the raw annot name, so the two
  hemispheres became separate keys and the number it reported compared *right-hemisphere* values
  against the published *bilateral* map — and `lh_unknown` / `lh_corpuscallosum` escaped the
  background filter. Fixed by normalising the hemisphere prefix before any matching; the validation
  is now bilateral-vs-bilateral and PLS2 agreement rises from r = 0.97 to **r = 0.985** (CT, MT and
  their deltas 0.998–0.9997).

### Is the resampling lossy? Measured: no — `code/19_resample_check.py`

An earlier version of this section claimed the resampled map's "effective resolution is the 308
parcellation, not 180, because neighbouring HCP parcels inherit one 308-value". **That reasoning was
wrong twice over** and the check below replaces it: 308 is the *bilateral* count, so the target grid
is the finer one, and an HCP parcel is a vertex-weighted average of however many source parcels it
overlaps rather than a copy of one.

| measurement | value |
|:--|:--|
| parcels per hemisphere | 152 in the 308 scheme, 180 in HCP-MMP, 34 in DK |
| 308-parcels contributing per HCP parcel (median) | 4 (>=5% of vertices: 3) |
| dominant source share per HCP parcel | 0.52 (0.40-0.67) |
| HCP parcels taking >90 % from one source | 12 of 180 |
| round trip HCP → 308 → HCP, ABCD dCT+CT | r = 0.909, rho = 0.877 (n = 137) |
| round trip HCP → 308 → HCP, AHBA C3 | r = 0.899, rho = 0.893 (n = 137) |
| 308 → DK direct vs 308 → HCP → DK | r = 0.961 (n = 34) |
| label alignment with the HCP arm | 137 shared of 179 and 137 |

A round trip through the 308 grid returns an HCP-native map at r ≈ 0.90, so resampling alone could
attenuate ρ = 0.75 to about 0.68 — not to 0.21. And the attenuation is specific to NSPN pairs
(|ρ|_HCP / |ρ|_DK: median 0.47 (0.28-0.69, n = 4) for pairs involving NSPN PLS2 against
median 1.20 (0.58-4.92, n = 6) for the rest), so it is not the general parcellation effect either.

**What it actually is: the ABCD component differs between parcellations.** Pushing each HCP-space
map down to DK regions and comparing with the native DK version of the same quantity:

| quantity | ρ, HCP version pushed to DK vs native DK |
|:--|--:|
| ABCD dCT+CT scores | +0.74 |
| AHBA C3 | +0.91 |
| dCT rate | +0.85 |

The reference and the imaging map survive the change of parcellation well; the *component* does not
— despite its gene weights agreeing at ρ = 0.90 between the two fits, so it is the score map that
moves rather than the ranking. And NSPN PLS2 tracks the DK version: the HCP-space component pushed
down to DK agrees with the NSPN DK map at only ρ = +0.29, against +0.75 for the native DK component.
So the weak HCP–NSPN cells say that the HCP fit is a somewhat different spatial axis — a fact about
the two ABCD fits rather than about the NSPN reference or the conversion.

### Basis or parcellation? The 0.74 decomposed — `code/20_basis_vs_parcellation.py`

The section above leaves one question: the dCT+CT component's score map correlates only ρ = 0.74
between its DK and HCP-MMP fits while its gene weights agree at 0.90. Those two fits differ in *two*
ways — the expression basis (AHBA_updated ds25, 12,007 genes vs abagen-data `hcp_3d.csv` on the
shipped 7,973) and the parcellation (33 vs 137). A third fit separates them: the same option-2 PLS on
the same 33 DK regions but with the HCP fit's expression basis (`dk_3d.csv` restricted to those
7,973 genes — the matched X of `12_hcp_pls.py`, refitted here so its *scores* exist too).

| contrast | what varies | scores ρ | gene weights ρ |
|:--|:--|--:|--:|
| native DK vs matched DK | expression basis only | **0.986** | 0.975 |
| matched DK vs HCP→DK | parcellation only | **0.763** | 0.749 |
| native DK vs HCP→DK | both | 0.737 | 0.733 |

**My hypothesis was wrong: the expression basis accounts for essentially none of it.** Swapping
abagen build and gene list while holding the parcellation fixed leaves the score map intact
(ρ = 0.986) and the gene ranking nearly so; swapping the parcellation while holding the basis
fixed reproduces the whole gap (0.763). The same holds for the NSPN comparison — the matched-basis
component agrees with NSPN PLS2 at ρ = 0.76, indistinguishable from the native DK fit's 0.75,
while the HCP fit pushed down to DK gives 0.47.

So fitting the PLS at 137 parcels genuinely finds a somewhat different spatial axis, not a
differently-processed version of the same one. Two caveats on the reading:

- **Which NSPN map is used matters at n = 33.** The published DK table and the 308-map resampled to
  DK agree at r = 0.985, but against the same component they give ρ = 0.75 and
  0.61 respectively — a 0.15 spread from the reference version alone. The figure uses the
  published table for DK cells and the resampled map for HCP cells, so part of its 0.75 → 0.21
  asymmetry is that version difference; both columns are in `basis_vs_parcellation.tsv`.
- Pushing the HCP map down to DK recovers some but not all of the agreement (0.47 at 33 regions
  against 0.21 at 137), so the disagreement is partly, but not only, at fine spatial scale.

### snRNAseq PC1 in the gene matrix

The snRNA-seq maturation axis from
[`transcriptional_maturation`](https://github.com/richardajdear/transcriptional_maturation) is
vendored here as `../data/velmeshev_PC1_gene_loadings.csv` (17,631 genes, two independent datasets
agreeing at ρ = 0.54), and `src/abcd/genemaps.snrnaseq_pc1` already loaded it. It is now a column
of the gene-weight matrix — the Herring version plotted, the U01 version in
`design_grid_pairs_weights.tsv`:

| vs snRNAseq PC1 (gene weights) | ρ |
|:--|--:|
| AHBA C3 | 0.47 |
| ABCD dCT+CT, DK | 0.31 |
| ABCD dCT+CT, HCP-MMP | 0.29 |
| ABCD dCT alone, DK | 0.25 |
| NSPN PLS2 | 0.10 |
| AHBA C1 (control) | -0.05 |

The maturation axis sits closest to C3 (0.47), then the ABCD signatures
(0.25–0.31), and is essentially unrelated to NSPN PLS2 (0.10) and to C1
(-0.05). Since C3 and NSPN PLS2 themselves agree at 0.55, the maturation content of this
family of axes sits in the C3 direction rather than in what NSPN PLS2 adds.

**The figure is now two square pair matrices**, regional maps (panel a, five variables) and gene
vectors (panel b, five, including snRNAseq PC1), each with DK below the diagonal and HCP-MMP above it, ρ in bold inside
every cell and the spin p under it. This is possible only because NSPN PLS2 exists in both
parcellations. Each matrix shades the row and column of its odd one out — dCT rate in panel a (an
imaging map, not a gene-derived axis) and snRNAseq PC1 in panel b (the only axis not derived from
the AHBA) — and DK is drawn in light red against HCP-MMP's blue. AHBA C1 is not a column of the gene
matrix — it is the static-gradient control
rather than one of the axes being compared — but its loadings are still computed in
`design_grid_pairs_weights.tsv` and quoted where they matter (the dCT-alone C1 loading is the
reason that design fails at 33 regions). The NSPN-vs-C3 cell of panel b is
parcellation-independent and marked "=".

### Gene bases — which comparisons are clean, and why

Three different gene sets are in play, and mixing them silently confounds parcellation with
pipeline. Stated once here because every cross-parcellation claim depends on it:

| basis | matrix | genes | used for |
|:--|:--|--:|:--|
| **matched** | `AHBA/data/abagen-data/expression/{hcp,dk}_3d.csv` restricted to `data/weights.csv` | 7,973 | the HCP-vs-DK arm — both atlases, same genes, same abagen build |
| AHBA_updated native-DK | `AHBA_updated/outputs/.../ahba_dk_lh_native_ds{0,25,50}.csv` | 16,009 / 12,005 / 8,005 | the main DK analysis (phases 1–3) |
| shipped C1–C3 | `data/weights.csv` | 7,973 | the published components |

The 7,973 genes of `weights.csv` are **not** an unfiltered set: they are exactly the columns of
`hcp_3d_ds5.csv`, i.e. the top half of 15,946 genes by differential stability computed **in HCP
space** (verified by set comparison). So the HCP arm carries a DS-50 filter of its own. DS filters
are parcellation-specific, which is why that list overlaps `dk_3d_ds5` by only 88 % and the
AHBA_updated `ds50` matrix by 87 %. An AHBA_updated `ds50` row therefore differs from the HCP fit in
**three** ways at once — parcellation, abagen build, gene identity — and is not a parcellation
comparator. (An earlier version of `fig_hcp_vs_dk.png` carried such a row, chosen because its gene
*count* was closest; it has been removed.)

Standardisation applied:

- `hcp_vs_dk_enrichment.tsv` now carries a **`basis`** column, and `fig_hcp_vs_dk.png` plots only
  the matched-basis rows — so on that figure parcellation is the only difference between an HCP row
  and its DK counterpart. Every HCP design now has a matched DK counterpart:
  `12_hcp_pls.py` fits the single-Y design on `dk_3d.csv` too (`DK_PLS1_dCTonly_matchedX`) and
  exports both matched vectors to `results/dk_matched_weights.tsv` for reuse.
- `14_magma_all_options.py` adds `ABCD_PLS2_DKmatched` and the **atlas joint model now uses it**:
  in `magma_all_joint.tsv` (nine covariates, n = 6,662 MDD / 6,672 SCZ) HCP vs DK on matched genes gives
  MDD β = 0.060 (p = 0.012) for HCP against 0.007 (p = 0.76) for DK, and for SCZ neither survives
  (0.33 and 0.22) — the same verdict as the cross-build version, now without the confound.
  The arm's own three-covariate run (`hcp_vs_dk_conditional.tsv`, n = 6,051) is a separate MAGMA
  model on a smaller universe and gives the same reading with slightly different coefficients
  (MDD HCP 0.067, p = 0.008; DK 0.003, p = 0.91; SCZ 0.39 and 0.27) — quote one or
  the other, never a mix.
- `15_celltype_all_options.py` adds the matched DK vector as a ninth ranking, which **settles the
  astrocyte question**: on the HCP gene basis the DK fit still gives Astro z = -5.8 and
  Oligo -5.8, essentially identical to the AHBA_updated DK fit (-5.2, -5.4) and opposite
  to HCP (+5.0, -12.1). The flip is parcellation, not the gene set or the abagen build.

One thing the matched single-Y row makes visible: **for SCZ the dCT-alone fits lead every PLS2 fit
in either atlas** — DK dCT-alone β = 0.065 (p = 0.0002) and HCP dCT-alone 0.054, against
0.046 and 0.045 for the two-Y fits, with only C3 (0.068) above them. For MDD the single-Y fits still lead,
but by much less — DK 0.059 vs 0.053 (a gap of 0.006 against 0.019 for SCZ) and HCP 0.078 vs 0.069.
Read with the dCT-only section below — those components are not spin-significant, so this is a
statement about gene rankings, not about spatial maps.

### The dCT-only variant — why the single-Y design is kept as a control, not a result

`code/16_dct_only.py` → `results/dct_only_*`; figure `figures/fig_signature_dct_only.png`
(`code/fig1_signature_dct_only.R`). Same panel structure as `fig_signature_both.png` so the two
read side by side; Y is the thinning rate alone, so the PLS has a single component and PLS1 *is*
the vector of gene–map correlations.

| | DK, 33 regions | HCP-MMP, 137 parcels |
|:--|--:|--:|
| component spin p (5,000 rotations) | 0.082 | 0.266 |
| bootstrap reproducibility | 0.83 | 0.88 |
| scores vs C3 (ρ, p_spin) | 0.28, 0.185 | 0.55, 0.0004 |
| weights vs C3 | 0.59 | 0.61 |
| weights vs **C1** (the static axis) | **0.56** | **0.28** |
| weights vs the two-Y signature | 0.70 | 0.90 |

Two things are clearer here than in the DK-only analysis:

- **The C1 confound is a coarse-atlas problem.** At 33 regions the dCT-only vector loads on the
  static gradient as heavily as on C3 (0.56 vs 0.59) — the reason the two-Y design was adopted.
  At 137 parcels the C1 loading collapses to 0.28 while the C3 loading holds at 0.61, and the
  vector is ρ = 0.90 with the HCP two-Y signature. With enough parcels, thinning rate alone
  nearly recovers the axis; the thinning map and the static map are less collinear when parcels
  are small enough not to straddle the gradient.
- **Neither single-Y component is spin-significant** (0.082 DK, 0.266 HCP), even though the HCP
  *scores* track C3 at p_spin = 0.0004. The component's total covariance with expression is within
  what spatial autocorrelation produces; what is spatially specific is the part aligned with C3.
  So the dCT-only vector is usable as a gene ranking (it is the strongest ABCD vector for both
  disorders in HCP space — see the arm above) but carries no independent spatial claim.

Cell classes agree with the main figure's reading: neuronal up, glial down in both parcellations,
with astrocytes again splitting by atlas (+1.4 in HCP-MMP, p = 0.16, against
-12.1 in DK) — i.e. the flip does not depend on the Y-matrix design.

### Cell-class profiles of all nine rankings — the astrocyte flip is a parcellation effect

`code/15_celltype_all_options.py` → `results/celltype_all_options.tsv`;
figure `figures/fig_celltypes.png` (`code/fig4_celltypes.R`). Supersedes
`13_celltype_compare.py` (three rankings) for figure purposes. Every ranking is
tested on **one shared universe** of 7,338 genes, so a difference between two
columns cannot come from a difference in universe; running each on its own
universe instead gives z agreeing at Spearman 0.994 over all 81 cells (both are
in the table).

| class | PLS2 HCP | dCT HCP | **PLS2 DK, HCP basis** | PLS2 DK | dCT DK | dCT+dT1T2 DK | C3 | NSPN PLS2 | C1 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| Neuro-Ex | +12.6 | +12.2 | +15.4 | +13.6 | +13.8 | +12.8 | +19.3 | +9.7 | +2.5 |
| **Astro** | **+5.0** | **+1.4** | **-5.8** | -5.2 | -12.1 | -2.6 | -8.9 | -11.3 | -11.9 |
| **Oligo** | **-12.1** | **-9.2** | **-5.8** | -5.4 | -1.6 | -7.0 | -13.2 | +1.8 | +5.5 |
| Micro | -11.3 | -11.9 | -11.3 | -10.0 | -11.7 | -8.3 | -8.8 | -1.6 | -4.1 |

(full 9 classes × 9 rankings with p_perm in the table and figure; Endo, OPC and Per omitted here)

- **Neuronal up, microglia/endothelia down in every thinning-derived ranking** —
  the Y-matrix option matters far less than the atlas (profiles correlate ρ 0.65–0.98
  among the 6 ABCD vectors; ρ 0.98 between the two HCP-MMP fits).
- **Astrocytes split by atlas, not by option — and not by gene set.** Both HCP-MMP rankings are
  the only ones in the astrocyte-positive half (+5.0, p < 0.001; +1.4, p = 0.16); every DK
  ranking and every published component is negative (-2.6 to -12.1). The decisive control is
  the DK fit on the **HCP gene basis** (-5.8): identical to the other DK fits, so neither the
  abagen build nor the DS gene set produces the flip.
- **Oligodendrocytes move the same way with the atlas**: -12.1 in HCP-MMP against
  -5.4 in DK, i.e. the HCP version sits next to C3 (-13.2) while the DK
  version does not. On this axis the finer parcellation moves the ABCD signature
  *towards* C3 even though its SCZ/MDD β does not overtake it.
- **C1 is the control that behaves differently**, as it should: weakly neuronal
  (+2.5) and oligodendrocyte-**positive** (+5.5).
- Caveat: marker genes are co-expressed, so the independent-gene null is
  anti-conservative. The signs and the column-to-column pattern are the result;
  the magnitudes are not calibrated.

**Open question.** Whether the astrocyte flip reflects real biology resolved by
137 parcels or a parcel-size artefact is not settled here. The cheap test is
whether it survives restricting the HCP fit to parcels inside DK regions that
are astrocyte-marker-rich, or a leave-one-region-out refit. It is not a
sample-size effect: it is unchanged (+5.1 → +5.0) between the 6,537-child and
8,716-child HCP runs.

### Disorder panel: the signature against every disorder gene analysis — `code/21_magma_disorder_panel.py`

**What has been tested so far.** Yes — the dCT+CT PLS2 gene weights have been tested against SCZ and
MDD by MAGMA gene-property analysis since Phase 3 (scripts 08 and 14), but only against the two gene
analyses shipped in `genetic_analysis/inputs/magma/`: `SCZ.genes.raw` = **PGC3 primary** and
`MDD.genes.raw` = **MDD2025 div**. Both are multi-ancestry GWAS whose gene analysis was run against
1000G EUR alone — the LD mismatch `genetic_analysis/README_HPC.md` step 10 documents. Script 21
replaces them with the LD-matched gene analyses CSD3 already produced, and runs whichever are present.

**Files needed.** All four exist on CSD3 from steps 10–11, and were copied into
`genetic_analysis/inputs/magma/` on 2026-09-24 (they are gitignored; on a fresh clone repeat this):

```sh
R=rajd2@login-q-1.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/magma_scz2025/genes
for g in SCZ25_META SCZ25_EUR PGC3_EUR MDD_EUR; do
  scp "$R/$g.genes.raw" "$R/$g.genes.out" genetic_analysis/inputs/magma/
done
```

| gene analysis | GWAS | LD reference | role |
|:--|:--|:--|:--|
| `SCZ25_META` | 2025 SCZ, AFR + EUR + EAS | per-ancestry 1000G panels, then `magma --meta` | **primary SCZ** |
| `SCZ25_EUR` | 2025 SCZ, European | `g1000_eur` (matched) | EUR sensitivity |
| `PGC3_EUR` | PGC3 european | `g1000_eur` (matched) | old-GWAS comparator |
| `MDD_EUR` | MDD2025 eur | `g1000_eur` (matched) | **primary MDD** |

**Does the covariate analysis need an LD reference?** Not at this step. The gene-property model
corrects for correlations between genes, but those correlations are read from the `.genes.raw` file
itself — they were computed once, at the *gene analysis* step, from that analysis's reference panel.
So the LD choice is made upstream, once per disorder GWAS, and the rule there is that the panel must
match the **GWAS sample's** ancestry, not ours. (The test itself is two-sided: MAGMA's default for gene
properties, unlike gene sets.)

**Should the ABCD be the LD reference?** No, for the disorder side. The reference has to match the
people whose genotypes produced the summary statistics; the ABCD's ancestry mix is not the 2025 SCZ
meta's (which has no AMR cohort and is ~15 % AFR by Neff), and a single ABCD-wide panel would be a
mismatch of the same kind as 1000G EUR on a multi-ancestry file. It would also carry relatedness
(siblings and twins) and imputation error into every gene's LD. The principled route for a
multi-ancestry GWAS is the one step 10 took and the MAGMA manual prescribes: gene analysis per
ancestry against the matching 1000G panel, then meta-analyse the gene Z (`SCZ25_META`). The panel's
size is not the constraint for gene-level tests (1000G EUR is 503 people), so there is nothing a
larger in-sample panel would buy. ABCD LD *is* the right choice in one place — the gene analysis of
**our own** ABCD GWAS — but that is the phenotype side, which this gene-property test does not use.

**Results, all six gene analyses** (β_std, two-sided p; each vector on its own universe, max
n = 16,723–17,367 genes; MAGMA's gene-property default is two-sided, so a negative coefficient can be
significant).

| gene weights | SCZ25_META | SCZ25_EUR | PGC3_EUR | PGC3_primary | MDD_EUR | MDD_div |
|:--|--:|--:|--:|--:|--:|--:|
| ABCD PLS2, DK (lead) | +0.030 (0.0084) | +0.018 (0.12) | +0.023 (0.051) | +0.023 (0.058) | +0.020 (0.082) | +0.028 (0.015) |
| ABCD PLS2, HCP-MMP | +0.036 (0.018) | +0.030 (0.053) | +0.044 (0.0047) | +0.038 (0.016) | +0.051 (0.00062) | +0.055 (0.0002) |
| ABCD PLS2, DK on HCP genes | +0.040 (0.0085) | +0.025 (0.11) | +0.036 (0.018) | +0.040 (0.012) | +0.036 (0.015) | +0.043 (0.0036) |
| AHBA C3 | +0.060 (8.6e-05) | +0.052 (0.001) | +0.056 (0.00037) | +0.058 (0.0003) | +0.071 (1.7e-06) | +0.072 (1.5e-06) |
| NSPN PLS2 | +0.027 (0.00089) | +0.024 (0.0053) | +0.024 (0.0069) | +0.031 (0.00074) | +0.014 (0.1) | +0.017 (0.055) |
| AHBA C1 (control) | +0.047 (0.0016) | +0.047 (0.0024) | +0.028 (0.066) | +0.032 (0.039) | +0.013 (0.37) | +0.013 (0.39) |

| joint model | coefficient | SCZ25_META | SCZ25_EUR | PGC3_EUR | PGC3_primary | MDD_EUR | MDD_div |
|:--|:--|--:|--:|--:|--:|--:|--:|
| ABCD DK + C1 | ABCD PLS2, DK (lead) | +0.038 (0.014) | +0.025 (0.12) | +0.032 (0.043) | +0.035 (0.03) | +0.034 (0.026) | +0.041 (0.0063) |
| ABCD DK + C1 | AHBA C1 (control) | +0.056 (0.00026) | +0.053 (0.00088) | +0.035 (0.027) | +0.041 (0.01) | +0.020 (0.18) | +0.021 (0.16) |
| ABCD DK + C3 | ABCD PLS2, DK (lead) | -0.035 (0.12) | -0.056 (0.019) | -0.034 (0.14) | -0.030 (0.21) | -0.054 (0.015) | -0.038 (0.082) |
| ABCD DK + C3 | AHBA C3 | +0.084 (0.00025) | +0.093 (9.2e-05) | +0.080 (0.00066) | +0.077 (0.0016) | +0.114 (4.6e-07) | +0.103 (4.6e-06) |
| ABCD HCP + C3 | ABCD PLS2, HCP-MMP | -0.009 (0.68) | -0.009 (0.66) | +0.012 (0.57) | -0.001 (0.96) | +0.005 (0.8) | +0.013 (0.52) |
| ABCD HCP + C3 | AHBA C3 | +0.066 (0.0015) | +0.058 (0.0072) | +0.048 (0.025) | +0.058 (0.0071) | +0.068 (0.00079) | +0.063 (0.0018) |
| ABCD HCP + ABCD DKmatched | ABCD PLS2, HCP-MMP | +0.013 (0.56) | +0.026 (0.27) | +0.038 (0.11) | +0.019 (0.44) | +0.055 (0.015) | +0.052 (0.02) |
| ABCD HCP + ABCD DKmatched | ABCD PLS2, DK on HCP genes | +0.030 (0.19) | +0.005 (0.82) | +0.008 (0.73) | +0.026 (0.28) | -0.005 (0.83) | +0.004 (0.87) |

Reading, primary files first (`SCZ25_META`, `MDD_EUR`):

1. **The 2025 SCZ GWAS strengthens the lead signature's SCZ association**: β = 0.030, p = 0.008, against
   0.023, p = 0.058 on PGC3 primary. The gain comes from the multi-ancestry meta, not from LD matching: the
   2025 European GWAS alone gives 0.018 (p = 0.12), and PGC3 european vs PGC3 primary barely differ
   (0.023 vs 0.023). Fixing the LD mismatch changes little for these tests.
2. **MDD: the HCP-MMP version carries it.** On `MDD_EUR` the HCP signature gives 0.051 (p = 0.0006) and
   keeps its association with the gene-matched DK fit in the model (0.055, p = 0.015), while the DK
   lead is marginal (0.020, p = 0.08). For SCZ neither parcellation survives the other.
3. **The static gradient C1 is itself SCZ-associated in the 2025 GWAS** (0.047, p = 0.0016; it was 0.032
   on PGC3), and not MDD-associated. The lead signature survives it (joint 0.038, p = 0.014), so it is
   not a proxy for C1.
4. **C3 still dominates.** C3 is the strongest ranking for every disorder file, and with C3 in the
   model neither ABCD version adds anything. The DK lead's joint coefficient turns *negative*, and
   nominally significantly so for `SCZ25_EUR` (-0.056, p = 0.019) and `MDD_EUR`
   (-0.054, p = 0.015). The two vectors correlate at ρ ≈ 0.76, so this is the usual
   suppression pattern in a collinear pair: the part of the ABCD ranking orthogonal to C3 leans
   slightly against disorder risk. It is not evidence of an opposing biological axis, and would need
   replication before being read as one.
5. **NSPN PLS2** is SCZ-associated in every SCZ file but not MDD-associated in the matched file
   (0.014, p = 0.10), matching the original PNAS claim, which was SCZ-specific.

**Tooling note.** The x86_64 macOS MAGMA stopped running here (no Rosetta), so `tools/bin/magma_src/`
now holds a native arm64 build from the v1.10 source; it reproduces `magma_all_marginal.tsv` exactly.
Build recipe in `tools/bin/README.md`.

### HCP-MMP summary slide — `code/22_hcp_summary.py` → `code/fig5_hcp_summary.R`

`figures/fig_hcp_summary.png` (16:9) summarises the HCP-MMP arm alone: the option-2 fit (Y = dCT + CT,
137 parcels) with **both** components against AHBA C1 and C3. Script 22 repeats that fit to recover the
PLS1 region scores 12_hcp_pls.py did not save, and asserts the refit's gene weights and PLS2 scores
match the saved ones. Tables: `results/hcp_summary_{maps.csv,map_pairs,gene_pairs,gene_weights,sets}.tsv`.

**Orientation.** PLS1 is flipped to align with AHBA C1, so positive = *thinner* baseline cortex (this
applies to `hcp_summary_*` and to `ABCD_PLS1_HCP` in `magma_disorder_panel.tsv`; `hcp_pls_weights.tsv`
keeps the raw SVD sign). PLS2 is in the thinning orientation used everywhere else.

- **PLS1 is the static axis and is C1**: ρ = -0.59 with CT, -0.11 with dCT; 0.99 with C1 across
  regions and 0.95 across genes. SCZ 0.030, p = 0.044 (C1: 0.047, p = 0.0016); MDD n.s.
- **PLS2 is the thinning axis and is C3**: ρ = -0.58 with dCT (dCT is negative mm/yr, so this is faster thinning), -0.01 with CT; 0.70 / 0.68
  with C3. SCZ 0.036, p = 0.018, MDD 0.055, p = 0.0002; C3 is stronger for both
  (0.060, p = 8.6e-05 / 0.072, p = 1.5e-06).
- **Panel e uses `MDD_div`**, so both disorders are multi-ancestry GWAS. The two files differ in LD
  handling: SCZ25_META is per-ancestry against matched panels, while MDD div is the multi-ancestry file
  against 1000G EUR. That makes little difference here (PLS2: 0.051, p = 0.00062 on MDD_EUR).
- **Panel layout** (letters as on the slide): a maps (CT / PLS1 / C1 over dCT / PLS2 / C3), b CT vs
  dCT, c components vs CT and dCT, d regions and genes vs C1/C3, e MAGMA, f cell classes and layers;
  e and f list PLS1, C1, PLS2, C3 top to bottom. The methods live here, not on the slide.
- **CT vs dCT (b)**: ρ = 0.09, p_spin = 0.63 over all 179 parcels. Baseline thickness and thinning rate
  are spatially unrelated, which is why Y = [dCT, CT] separates cleanly into a static and a thinning
  component.
- **Colour scales (a)**: CT and dCT are white-anchored — CT white (thinnest) to blue (thickest), dCT white
  (0 mm/yr, no thinning) to red (fastest thinning). Component and AHBA maps use Spectral, so they are not
  confused with panel f's red–blue enrichment scale. No panel carries a subtitle; points are uncoloured.
- **Parcel systems**: the Glasser 2016 cortex lookup is kept in `hcp_summary*_maps.csv` but no longer drawn,
  because colouring the points made the scatters harder to read. It lives in
  `data/reference/hcp_cortices/hcp_parcel_systems.csv`, from `HCP-MMP1_UniqueRegionList.csv` `Cortex_ID` with
  the groupings of `HCP-MMP1_cortices.txt`: cortices 1–12 sensorimotor (51 of 137 parcels), 13–22
  association (86). Downloaded 2026-09-24 from bitbucket.org/dpat/tools, linked from neuroimaging-core-docs.
- **Fading (alpha 0.3)**: b–d region panels at spin p ≥ 0.05; e at p ≥ 0.05; f at BH q ≥ 0.05. Gene-level ρ
  has no valid p (7,973 co-expressed genes make every ρ nominally significant), so each gene panel is
  faded with its matching region pair.
- **Panel f**: square size = |z| (capped at 20), fill = direction. Enrichment z uses an
  independent-gene null, which marker co-expression makes anti-conservative. Layers = Maynard 2021,
  FDR < 0.05 and t > 0, the `which='maynard'` definition of `AHBA/code/enrichments_data.get_layer_genes`.
- **Global slope**: tested (`slope_{1lmm,perregion}_{HCP,DK}` in `magma_disorder_panel.tsv`) and left off
  the slide — none of PLS1, PLS2, C1, C3 is associated with any of the four gene analyses (smallest p = 0.28).

- **Other details that were on the slide**: ABCD 7.0, 8,716 children and 26,946 sessions (`hcp_run_provenance.tsv`);
  CT and dCT are the intercept and age slope of one mixed model per parcel; 179 bilateral HCP-MMP
  parcels, 137 with AHBA coverage; PLS on the 7,973-gene C1–C3 gene set, spin p of the singular values
  PLS1 < 0.001 and PLS2 = 0.011; map ρ are Spearman with 5,000-rotation spin tests; MAGMA gene-property
  tests are two-sided.

#### Three versions by AHBA matrix

`22_hcp_summary.py <variant>` and `fig5_hcp_summary.R <variant>` build the same slide from three AHBA matrices
(`~/Git/AHBA/data/abagen-data/expression/`), fitting the components from scratch each time. The 3d_ds5
fit reproduces 12_hcp_pls.py's saved weights, bootstrap Z and PLS2 scores (asserted).

- `3d_ds5` → `figures/fig_hcp_summary.png`: `hcp_3d_ds5.csv`, the ≥3-donor region filter plus the DS5 gene
  filter. This is the matrix C1–C3 were fitted on, and it is identical to `hcp_3d.csv` restricted to those genes.
- `ds5` → `fig_hcp_summary_ds5.png`: `hcp_ds5.csv`, no region filter but its own DS5 gene set. Its 7,973
  genes share 7,862 with C1–C3.
- `base` → `fig_hcp_summary_base.png`: `hcp_base.csv`, no filter.

In every version parcel ids 1–180 are the left-hemisphere `regionID`s of `HCP-MMP1_UniqueRegionList.csv`,
checked against the 137-parcel id→label table. Parcels with no samples and genes with any missing
value are dropped (one gene in each unfiltered matrix). AHBA C1/C3 are always the published maps and
weights, so the region comparisons with C1/C3 use their 137 parcels even when PLS covers 177.

| | 3d_ds5 | ds5 | base |
|:--|--:|--:|--:|
| matrix | hcp_3d_ds5.csv | hcp_ds5.csv | hcp_base.csv |
| parcels × genes | 137 × 7,973 | 177 × 7,972 | 177 × 15,636 |
| PLS1 / PLS2 cov | 86% / 14% | 86% / 14% | 81% / 19% |
| PLS2 spin p | 0.0106 | 0.0006 | 0.0028 |
| PLS1–C1 regions / genes | 0.99 / 0.95 | 0.99 / 0.94 | 0.97 / 0.90 |
| PLS2–C3 regions / genes | 0.70 / 0.68 | 0.72 / 0.71 | 0.73 / 0.69 |
| PLS2–dCT | -0.58 | -0.55 | -0.52 |
| PLS2 SCZ β (p) | 0.036 (0.018) | 0.036 (0.017) | 0.033 (0.0006) |
| PLS2 MDD β (p) | 0.055 (0.0002) | 0.054 (0.00027) | 0.040 (7.1e-05) |
| PLS1 SCZ β (p) | 0.030 (0.044) | 0.030 (0.045) | 0.019 (0.04) |
| PLS2 astro z | +5.3 | +3.9 | +2.9 |
| min p, global slope | 0.28 | 0.29 | 0.29 |

- **The signature is robust to both filters.** Without the region filter the components become more
  C3-like and more spin-significant (PLS2 p 0.011 → 0.0006). Removing the gene filter as well keeps
  the same structure: PLS2–C3 ρ = 0.73 on regions and 0.69 on genes.
- **MAGMA β is smaller without the gene filter, but p is smaller too.** The `base` vectors are tested on
  roughly twice the genes (own universe), so β_std is diluted by the lower-stability genes while the
  larger n tightens the estimate. It is the same pattern as the DK signature on its wider gene list.
- **The astrocyte loading shrinks as filters are removed** (+5.3 → +2.9), moving toward C3's negative value.
- **Global slope stays null** in all three versions.

## Reproducing

Analysis (python, env `ahba-pls`): `code/01_*` → `code/22_*` in order. `11_` is the parcellation
control, `12_` the HCP-MMP arm, `14_`/`15_` the single-universe enrichment and cell-class tables that
the summary figures read, `16_` the dCT-only variant, `17_` the design grid and all pairwise statistics, `18_` the NSPN→HCP-MMP
conversion (needs `nibabel`, and the fsaverage annot files in `~/Git/AHBA/data/parcellations/`; run it
before `17_`), `19_` the resampling audit behind that conversion and `20_` the basis-vs-parcellation decomposition (both seconds; `code/surface.py` holds the shared annot helpers). `12_` needs the backfilled HCP run
(`ABCD_CONFIG=ct_70_hcp_noglobal_mv2 python -m abcd.assemble` then
`Rscript R/fit_lmm.R --cores 8`, which lands in `out/thickness_hcp_70_aa6e91efba82`); everything else
reads the DK and T1w/T2w runs listed in `tools/rerun_local.sh`.

Figures (R, env `ahba-pls-r`, with `ahba_pls/.Rlib` on `.libPaths()` for `ggsegGlasser`):

```
Rscript code/fig1_lead_signature.R      # DK-only lead signature
Rscript code/fig1_signature_both.R      # both parcellations + cell classes
Rscript code/fig2_enrichment.R          # permutation + MAGMA, DK
Rscript code/fig2_enrichment_combined.R # MAGMA only, all rankings
Rscript code/fig3_hcp_vs_dk.R           # the HCP vs DK arm
Rscript code/fig4_celltypes.R           # cell classes, all nine rankings
Rscript code/fig1_signature_designs.R   # brain maps + the two pair matrices
Rscript code/fig1_signature_dct_only.R  # the dCT-only design on its own
```

MAGMA v1.10 binaries and their provenance are in `../tools/bin/README.md`; the disorder
`.genes.raw` files come from `../genetic_analysis/inputs/magma/` and are gitignored, so a fresh
clone must copy them from the cluster before `08_`, `11_`, `12_` or `14_` will run.
`.gitignore` here excludes the regenerable intermediates (aligned X matrices, full per-option
weight tables, spin nulls, MAGMA run directories).

## Decisions and caveats

- 2026-09-12: ds75 dropped (user). H3 deferred until a clean ABCD GWAS exists; H4 deferred to a cluster run (user).
- All PLS runs on the 33 AHBA-covered bilateral DK regions; the missing region is recorded in `data/`.
- 2026-09-12: figures are built in R (`ggplot2` + `patchwork` + `ggseg` 2.2.1 from CRAN, env `ahba-pls-r`), not matplotlib. ggseg's R atlas renders frontal and temporal pole correctly, unlike the python polygon set. The R scripts fit nothing and hardcode nothing: every statistic printed on a figure is read from the `results/` table it plots. DK polygons are cached to `data/dk_polygons.csv`.
- Figure conventions: one-line interpretation per panel subtitle, overall interpretation as the figure subtitle, methods as short bullets in the caption.
- 2026-09-15: the HCP-MMP arm moved to the backfilled run `thickness_hcp_70_aa6e91efba82` (8,716 children, the DK sample); the 6,537-child run `thickness_hcp_70_fec93121f0dd` is superseded and kept only for the stability comparison quoted above.

## Data vintage — what this directory now runs on

Two corrections landed after the first version of this analysis, and **neither changed a conclusion**.

### 1. 2026-09-14 — the tables were 6.0, not 7.0

Every input map and every result here was regenerated from runs on the true 7.0 tabulation
(`../docs/RERUN_7.0_TABULATED.md`). The DK runs kept their run ids (the config hash does not encode
data vintage) but now contain **8,716 children / 26,949 sessions**, with the six-year wave complete.
Scripts `01`–`15` and the R figures were re-run unchanged apart from path updates (MAGMA disorder files
now at `../genetic_analysis/inputs/magma/`, gene-set sources under `../legacy/hpc/work/`).

**Nothing changed in kind; every number moved by a rounding step.** The group maps that feed the PLS
correlate ρ = 0.996 with the 6.0-vintage maps, which is both the expected outcome and a useful
stability check on the whole chain:

| quantity (lead signature = option 2, PLS2, ds0 unless stated) | 6.0-vintage | 7.0 tables |
|:--|--:|--:|
| PLS2 spin p (singular value) | 0.0016 | 0.0022 |
| PLS2 scores vs NSPN-PLS2 (ρ, p_spin) | −0.69, 0.002 | −0.73, 0.001 |
| gene weights vs NSPN-PLS2 / C3 (ρ over shared genes) | 0.61 / 0.76 | 0.62 / 0.76 |
| top-decile overlap with C3, fold | 5.3 | 5.2 |
| MAGMA gene-property β_std, SCZ: ds0 / ds25 / ds50 | 0.037 / 0.028 / 0.034 | 0.037 / 0.027 / 0.034 |
| MAGMA gene-property β_std, MDD: ds0 / ds25 / ds50 | 0.026 / 0.028 / 0.041 | 0.025 / 0.026 / 0.039 |
| MDD high-confidence set, permutation z (ds0) | 2.82 | 2.80 |
| HCP-MMP PLS2 scores vs C3 (ρ, 137 parcels) | 0.69 | 0.70 |
| HCP-MMP PLS2 MDD β_std / conditional on DK, p | 0.070 / 0.009 | 0.069 / 0.008 |
| DK PLS2 (matched X) MDD β_std / conditional on HCP, p | 0.055 / 0.78 | 0.053 / 0.91 |
| SCZ β_std HCP vs DK (matched X) | 0.044 vs 0.046 | 0.045 vs 0.046 |
| astrocyte z, HCP-MMP PLS2 / DK PLS2 / C3 | +5.1 / −5.3 / −8.9 | +5.1 / −5.2 / −8.9 |
| oligodendrocyte z, HCP-MMP PLS2 / DK PLS2 / C3 | −12.3 / −5.7 / −13.2 | −12.0 / −5.4 / −13.2 |

Sources: `results/pls_components.tsv`, `concordance_scores.tsv`,
`lead_overlap_with_references.tsv`, `enrichment_summary.tsv`, `hcp_concordance.tsv`,
`hcp_vs_dk_conditional.tsv`, `celltype_all_options.tsv`.

### 2. 2026-09-15 — the HCP-MMP arm refitted on the backfilled parcellation

The HCP-MMP arm had been running on a partial parcellation (24,921 of 33,825 FreeSurfer sessions →
6,537 children) with 6.0-vintage covariates for the extra scans. After `tools/hcp_backfill.sbatch`
(33,795/33,825 sessions) and the 7.0 covariates, the arm was refitted on
`out/thickness_hcp_70_aa6e91efba82` — **8,716 children, the same sample as DK**, 179 regions with
parcel `H` excluded at fit time, all 358 hemisphere-region models converged and none singular.

| quantity (HCP-MMP arm) | 6,537 children | 8,716 children |
|:--|--:|--:|
| PLS2 spin p / bootstrap reproducibility | 0.010 / 0.898 | 0.011 / 0.899 |
| scores vs C3 (ρ, 137 parcels) | 0.702 | 0.701 |
| gene weights vs C3 / vs DK PLS2 | 0.677 / 0.733 | 0.678 / 0.733 |
| MDD β_std, HCP PLS2 / DK matched | 0.069 / 0.053 | 0.069 / 0.053 |
| MDD joint: HCP \| DK / DK \| HCP (p) | 0.008 / 0.91 | 0.008 / 0.91 |
| SCZ β_std, HCP PLS2 / DK matched | 0.045 / 0.046 | 0.045 / 0.046 |
| astrocyte z, HCP PLS2 | +5.1 | +5.0 |
| oligodendrocyte z, HCP PLS2 | -12.0 | -12.1 |

This is the more useful of the two checks, because the sample grew by a third and the HCP arm's whole
point — the MDD gain over DK — could have been a sample artefact. It was not: adding 2,179 children
moves the MDD β by 0.0001 and leaves every joint-model verdict intact. The prose and tables above
are the 8,716-child versions throughout; the figures were re-rendered from them.

Still deferred: **H4** (`FOLLOWUP_GENETICS.md`), a projection phenotype for the ABCD GWAS. **H3** was
run on the cluster and is null (Phase 4 above).
