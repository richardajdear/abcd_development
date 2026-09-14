# Re-run on the 7.0 tabulated tables — what changed (2026-09-14)

**Read this before `REPORT_7.0.md`.** Every table and figure in `docs/` was
regenerated on 2026-09-14 from runs on the true ABCD 7.0 tabulated tables. The
prose of `REPORT_7.0.md` was written between 2026-07-30 and 2026-09-12 against
runs that, it turned out, used the **6.0** tabulation under a directory
labelled 7.0. The design decisions in that report (no global covariate, no
family random effect, ≥2 visits, the QC stack) were re-checked on the new
tables and stand; the sample sizes and most numbers in its prose are the old
ones. This document gives the new headline numbers and the old-vs-new
comparison; `docs/vintage_comparison.csv` is the source
(`tools/compare_vintage.py`).

## 1. What the 7.0 tabulation adds

The FreeSurfer surfaces on CSD3 and the genotypes were always 7.0. Only the
*tabulated* tables — thickness, QC flags, ages, site, scanner, family ids,
ancestry PCs — were 6.0 (data freeze February 2024). The 7.0 tabulation (freeze
1 August 2025, released 13 May 2026) covers the rest of the six-year imaging
wave; imaging still ends at the six-year visit, so the gain is entirely
six-year scans:

| | 6.0 tables | 7.0 tables |
|:--|--:|--:|
| six-year thickness rows in `mr_y_smri__thk__dsk` | 4,086 | 7,607 |
| scans with imaging (all waves) | 30,276 | 33,794 |
| dropped by the release T1 flag / defect count / Philips | 682 / 905 / 3,346 | 338 / 1,238 / 3,802 |
| scans entering the settled model | 23,409 | 26,949 |
| six-year scans in the model | 3,539 | 6,510 |
| children with ≥2 QC-passing scans | 8,192 | **8,716** |
| children with 2 / 3 / 4 scans | 2,997 / 3,365 / 1,830 | 2,204 / 3,507 / **3,005** |
| mean scans per child | 2.86 | 3.09 |
| families | 6,871 | 7,285 |
| phenotyped and genotyped (in `gn_y_genrel`) | 8,082 | **8,596** |
| age centre (years) / oldest scan | 12.44 / 17.7 | 12.80 / 18.0 |
| median slope reliability (BLUP) | 0.157 | **0.211** |
| effective N for the slope (Σ reliability) | 1,361 | **1,752** |

531 children are new to the sample and 7 dropped out of it (their scans were
re-flagged by the revised release QC). 2,478 of the 8,185 shared children
gained a scan.

**The group map is unchanged; the subject-level phenotype is not.** The 68
regional thinning rates correlate ρ = 0.996 (Spearman) between the two fits,
and the intercept map 0.9999. `global_slope` on the shared children correlates
r = 0.92 — the extra scans re-estimate each child's slope. That is the point of
the exercise: heritability and polygenic-score analyses spend their power on
subject-level precision, and that is what improved.

## 2. Other differences between the two tabulations that touch this pipeline

- `ab_g_dyn__design_mr__manufact` is an integer code in 7.0 (1 GE, 2 Philips,
  3 Siemens) where 6.0 carried vendor strings. `Release70Adapter.scanner`
  decodes it; without that the Philips exclusion matched nothing. The mapping
  was established on the 30,367 imaging visits present in both tables.
- The release T1 inclusion flag and the manual FreeSurfer QC scores were
  revised (~900 sessions); the defect counts changed for 11 sessions; regional
  thickness values changed for at most 25 sessions per region.
- Ancestry PCs (`ab_g_stc__gen_pc__*`) were recomputed for all 11,670 genotyped
  children, and `gn_y_genrel` family/birth ids re-coded. `design_site` is
  unchanged. `ab_g_stc__design_id__group` was removed (the adapter tolerates it).
- The 7.0 tables also carry seven-year non-imaging sessions (`ses-07A`), which
  the visit map drops as before.

## 3. Code and layout changes made for the re-run

- `abcd-data-release-6.0/` holds the old tables; `abcd-data-release-7.0/` the
  new ones plus `processed/hcp/` (moved from the former `abcd-7.0/`).
- `io.Release60Adapter` (`release: "6.0"`) and
  `configs/ct_60_noglobal_mv2_genetic.yaml` reproduce the old sample
  exactly (model table identical to the archived run except the `release`
  column), so code changes and data changes are separable.
- Every run manifest now records `data_vintage` (rows per session of the
  imaging and covariate tables) and assembly refuses a 7.0 run from tables
  with fewer than 7,000 six-year rows (`Release70Adapter.assert_vintage`).
- The 6.0-vintage run directories are archived under
  `out/legacy_6.0_tabulated/` (gitignored, local). Their run ids collide with
  the new ones because the config hash does not encode data vintage; the
  manifest's `data_vintage` key tells them apart.
- `tools/rerun_local.sh` runs every specification the report needs; `R/fit_lmm.R`
  and the Makefile now honour `PY` for the interpreter.
- `hpc/`, `hpc_v2/`, `hpc_v3/` moved to `legacy/`; the genetics re-run is
  specified in `genetic_analysis/README_HPC.md`.

## 4. Headline numbers on the 7.0 tables

Settled specification (`ct_70_noglobal_mv2_genetic`): 8,716 children, 26,949
scans, 68 regions, 19 sites; every regional fit converged, none singular
(`fit_summary.csv`).

**Reliability and the visit filter** (`reliability_grid.csv`): median slope
reliability 0.230 at ≥2 visits (effective N 1752),
0.242 at ≥3 (n = 6,512, effective N 1504),
0.262 at 4 (n = 3,005, effective N 793).
The ≥2-visit filter remains the right one. Against release 5.1 at the same
filter the effective N is 2.1× (`handoff_release_comparison.csv`).

**Map-to-transcriptome** (`ahba_vs_maps_noglobal.csv`, 34 bilateral regions,
spin nulls): absolute thinning rate vs C3 ρ = -0.546, p_spin = 0.002; vs C2
ρ = -0.330, p_spin = 0.016; slope PC3 vs C2 ρ = +0.854, p_spin = 0.001; slope PC2 vs C1
ρ = -0.812, p_spin = 0.001. Identical to the 6.0-vintage values to three decimals —
the group maps did not move.

**Falconer heritability** (`h2_candidate_phenotypes.csv`; 273 MZ and 448 DZ
pairs from `gn_y_genrel`, family effect off, 95 % bootstrap CI):

| phenotype | h² |
|:--|:--|
| baseline thickness (control) | 0.73 [0.57, 0.90] |
| global mean slope | 0.48 [0.24, 0.71] |
| slope PC1 | 0.50 [0.19, 0.79] |
| slope PC3 | 0.35 [0.06, 0.64] |
| slope PC2 | 0.23 [-0.08, 0.51] |

The global slope's Falconer h² rises with reliability (0.44 on the 6.0 tables);
with the family random effect on, the control returns h² = 1.46 — the
artefact §4.2 of the report describes, reproduced (`h2_family_effect_contrast.csv`).

**Site and scanner** (`site_scanner_icc.csv`): ICC of the global slope 0.029
by site (19 levels), 0.001 by scanner manufacturer.

**Phenotype priority for the GWAS export** (`gwas_phenotype_priority.csv`) is
unchanged in order: global mean slope (primary), slope PC3, slope PC2, slope
PC1, with baseline thickness as the positive control. The export at
`out/thickness_dsk_70_139406217085/gcta_inputs/` has 8,716 rows, FID = family
id, and the ten 7.0 ancestry PCs; 8,596 of the children are genotyped.

**HCP-MMP arm** (`ct_70_hcp_noglobal_mv2`, derived 360-parcel table, same
partial session coverage as before): the 7.0 covariates give the six-year HCP
sessions an age and a QC row, so the run grows from 5,947 to **6,537**
children and from 15,969 to 18,633 scans (six-year scans 2,475 → 4,517);
3 of 360 parcels singular. It stays secondary to DK until the re-parcellation
is complete.

**Imaging transcriptomics** (`../ahba_pls/`, all fifteen scripts and six figures
re-run): nothing changes in kind. The lead PLS component's concordance with
NSPN-PLS2 (ρ = −0.73 vs −0.69) and C3 (weights ρ = 0.76), its SCZ/MDD
gene-property betas (0.037 / 0.025 at ds0), the HCP-vs-DK MDD contrast and the
astrocyte sign flip by atlas all reproduce to a rounding step; the dated
section at the end of `../ahba_pls/README.md` has the side-by-side table.
