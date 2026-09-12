# abcd_development

**Which genes drive adolescent cortical development?**

Imaging genetics of *longitudinal change* in the ABCD study. The phenotype is a
per-subject **rate** of cortical thinning — a random slope of thickness on age
estimated from repeated scans — not thickness at one timepoint.

The hypothesis: genes driving this rate are (i) enriched for GWAS signal from
disorders with adolescent onset, schizophrenia and major depression, and
(ii) connected to this group's prior transcriptional results (AHBA components
C1–C3, and the leading component of an independent snRNA-seq analysis).

Why the design matters: the genetic architecture of a developmental *rate* need
not resemble that of a static measure. Almost all published brain-imaging GWAS
use cross-sectional phenotypes, which average over exactly the variation of
interest.

## Status — 2026-09-12

Genetics now runs on the **full 7.0 sample** (11,670 genotyped, all ancestries;
cross-ancestry GRM built), and the PRS work has settled into a **four-method
grid** — C+T, PRS-CS, SBayesR, SBayesRC — with the discovery GWAS matched to
the target arm and a control battery (ASD, two ALZ releases ± APOE, polygenic
education). Canonical results: `hpc_v2/work/results_v2/prs_final/table_main.tsv`
(documented in [`hpc_v2/README_HPC.md`](hpc_v2/README_HPC.md) §14; one-slide
summary [`hpc_v3/slide_prs_methods.png`](hpc_v3/slide_prs_methods.png)).

Alongside it, `ahba_pls/` adds the **imaging-transcriptomics arm**: the
NSPN-PLS2 / AHBA-C3 transcriptomic signature of adolescent thinning, re-derived
from the ABCD maps and tested against SCZ and MDD genetics
([`ahba_pls/README.md`](ahba_pls/README.md)). It replicates, but it does not
improve on C3 as a gene ranking — see *Imaging transcriptomics* below.

- **SCZ PRS → faster thinning is the robust result**: significant in 3 of 4
  methods in *both* the EUR arm (n = 4,116; β −0.035 to −0.047 SD/SD) and the
  within-ancestry-standardised pooled arm (n = 8,082), and it survives
  family-level FDR under 3 of 4 method choices
  (`hpc_v3/prs_tables/prs_fdr_sensitivity.tsv`).
- **MDD is suggestive, not established**: same (negative) direction under every
  method, arm and design, but significant only in the pooled arm (2/4 methods
  after within-ancestry standardisation); the EUR arm has ~46 % power for the
  observed effect size, so its null is uninformative.
- **The controls behave**: ASD null in all four methods (its one earlier hit
  was a mismatched-stratum artifact); AD-without-APOE is a C+T-only signal that
  vanishes under joint SNP modelling and in the proxy-free Kunkle release;
  polygenic education runs the *opposite* direction to the disorders.
- **Do not quote pre-grid PRS numbers** — everything before
  `prs_final/` predates the age-covariate fix, the allele-frequency fix and
  the discovery-to-arm matching (`hpc_v2/README_HPC.md` §14.5–14.8).

## Where to look

| you want | go to |
|:---|:---|
| **the cluster genetics: state, results, and the current task** | [`hpc_v2/README_HPC.md`](hpc_v2/README_HPC.md) (v1 pipeline: `hpc/README_HPC.md`) |
| the findings, their caveats and the corrections | [`docs/REPORT_7.0.md`](docs/REPORT_7.0.md) |
| how the mixed model works and why | [`notebooks/01_longitudinal_model.qmd`](notebooks/01_longitudinal_model.qmd) |
| spatial nulls and the map-to-gene tests | [`notebooks/02_maps_and_genes.qmd`](notebooks/02_maps_and_genes.qmd) |
| heritability and phenotype choice | [`notebooks/04_heritability.qmd`](notebooks/04_heritability.qmd) |
| the imaging-transcriptomics PLS study (NSPN-PLS2 / AHBA-C3 re-derivation, SCZ & MDD enrichment) | [`ahba_pls/README.md`](ahba_pls/README.md), notebook [`ahba_pls/imaging_transcriptomics.qmd`](ahba_pls/imaging_transcriptomics.qmd) |
| the 5.1 draft this supersedes | [`docs/REPORT_5.1_legacy.md`](docs/REPORT_5.1_legacy.md) |

## Key findings so far

**Phenotype and modelling**

- **Slope reliability, not sample size, binds.** Median regional slope
  reliability is 0.149 at ≥2 visits, 0.229 at ≥3. Effective N rises from 833 on
  release 5.1 to 1,361 on 7.0 at the same filter.
- **Use `min_visits: 2`, not 3.** The ≥3 filter buys per-subject precision but
  discards 37 % of subjects; effective N falls and the group map is unchanged
  (ρ = 0.998).
- **The family random effect must be omitted for genetic phenotypes.** With it
  enabled the baseline-thickness control returns *h²* = 1.46 — impossible. It
  leaves the group map exactly unchanged (ρ = 1.000) while destroying the
  subject-level genetic signal, which is why it went undetected in the 5.1 draft.
- **The developmental map is transcriptionally patterned.** Absolute thinning
  rate vs AHBA C3: ρ = −0.546, p_spin = 0.0016, reproducing the 5.1 result across
  a release change and a pipeline rewrite. The strongest couplings are on the
  slope components (PC3–C2 ρ = +0.854, PC2–C1 ρ = −0.812).

**Imaging transcriptomics** (`ahba_pls/`, 2026-09-12)

- **The NSPN-PLS2 / AHBA-C3 "signature of adolescent thinning" re-derives from
  the ABCD.** PLS of AHBA expression on the ABCD thinning map recovers a
  component that matches both prior signatures in regional scores (ρ = 0.76 with
  NSPN-PLS2, 0.85 with C3; p_spin = 0.001) and in gene weights (ρ = 0.61 and
  0.77 over 11.2k / 7.9k genes), with the same neuronal-up / glial-down cell-class
  profile. The ABCD and NSPN thinning maps themselves agree (ρ = 0.64,
  p_spin = 0.002), so this is replication in a 20× larger cohort, not atlas
  circularity.
- **Thinning rate alone does not isolate it.** With dCT as the only Y variable
  the component is marginal (p_spin = 0.06–0.10) and loads as heavily on C1, the
  static expression gradient, as on C3. Adding baseline thickness as a second Y
  column pushes the static axis into PLS1 and leaves the thinning signature as
  PLS2 — the same structure as the PNAS 2016 design. It is stable across
  differential-stability filters (ρ ≥ 0.99), so gene filtering is not the
  constraint it was for deriving C3 by PCA.
- **It carries SCZ and MDD risk, but adds nothing to C3.** MAGMA gene-property
  regression is positive for both disorders at every filter (SCZ β 0.028–0.037,
  MDD 0.026–0.041) at about half C3's effect size; conditioning on C3 removes it
  entirely, while C3 survives conditioning on it. The prioritised-gene
  permutation test detects only the MDD high-confidence set (z = 2.4) — and the
  SCZ locus-pool signal sits on the static axis, not on thinning. **The
  contribution is a developmental warrant for C3, not a better gene list.**
- Still open (specified for the cluster in
  [`ahba_pls/FOLLOWUP_GENETICS.md`](ahba_pls/FOLLOWUP_GENETICS.md)): whether the
  signature weights relate to ABCD's own thinning GWAS, and whether projecting
  subject slopes onto the component beats `global_slope` for heritability.

**Genetics (v1 run: release 4.0 genotypes, EUR only, N = 4,119 — the PRS rows
are superseded by the 7.0 grid in Status above; h²/rg/GWAS conclusions stand)**

- **h² of `baseline_thickness` = 0.575 ± 0.147**, independently corroborated by
  LDSC at 0.584 ± 0.131 — a joint validation of the GRM, ID alignment and
  covariates. `global_slope` h² = 0.228 ± 0.141.
- **No genome-wide-significant loci**, as expected at this N. λ_GC 1.00–1.03.
- **Genetic correlation is uninformative for the slope**, not merely null: h² z =
  1.45 against LDSC's z > 4 guidance. An interpretable rg needs N ≈ 11,300 —
  above ABCD's phenotyped ceiling of 8,192, so rg is out of reach in ABCD alone.
- **SCZ PRS → faster thinning** is the lead worth pursuing, and it is
  disorder-specific: the regional SCZ and MDD PRS maps are uncorrelated.
- **T1w/T2w** was tested as an alternative metric: more heritable (0.595 vs
  0.444) but 3.4× more site-confounded, with no PRS association. A different
  phenotype, not a replication.

**Where the sample goes.** 11,868 enrolled → 8,192 with ≥2 usable visits →
5,678 genotyped (EUR only) → 4,119 analysed. The genotype shortfall is an
**ancestry restriction and release vintage, not QC attrition**: ABCD's 7.0
curated genotypes cover 11,670 of 11,868, with the ~198 difference accounted for
by subjects never genotyped. Details in
[`hpc/README_HPC.md`](hpc/README_HPC.md) §4.

## Running it

```bash
export ABCD_CONFIG=ct_70_noglobal_mv2_genetic     # the settled specification
make all                                          # assemble -> fit -> phenotype -> gcta
```

`make help` prints the active config and the run directory it resolves to. Every
analysis choice is a `RunConfig` field that hashes into the `run_id`, so runs
cannot silently overwrite each other and `out/<run_id>/config.yaml` records what
produced the numbers.

The same four steps run directly (needs `PYTHONPATH=src` or `pip install -e .`):

```bash
python -m abcd.assemble              # tidy long table
Rscript R/fit_lmm.R --cores 8        # per-region lme4 fits
python -m abcd.phenotype             # BLUPs + reliability
python -m abcd.gcta_export           # GCTA/MAGMA inputs
```

`ABCD_ROOT` is **not** required: 7.0 is vendored at `abcd-data-release-7.0/`
(gitignored — access-controlled) and 5.1 is found under `~/Git/ABCD`.

## Layout

```
src/abcd/        # Python: assembly, QC, phenotypes, spatial stats, gene work
R/               # model fitting (lme4)
configs/         # one YAML per specification; the run_id is a hash of it
docs/            # REPORT_7.0.md + every table and figure it cites
tools/           # regenerators for every table and figure in the report
hpc/             # SLURM genetics pipeline  -- see hpc/README_HPC.md
  legacy/          superseded HPC docs; read only if the repo contradicts the current one
ahba_pls/        # imaging transcriptomics: PLS of AHBA expression on the ABCD
                 # thinning maps, and its SCZ/MDD enrichment -- self-contained,
                 # see ahba_pls/README.md
notebooks/       # explanatory documents, not analysis scripts
tests/           # 152 tests, incl. provenance and README checks
```

The Python/R seam is Parquet in `out/<run_id>/`. Computation is kept separate
from plotting throughout.

| config | what it is for |
|:---|:---|
| `ct_70_noglobal_mv2_genetic.yaml` | **the settled specification** |
| `ct_70_noglobal_mv{2,3,4}.yaml` | the visit-filter comparison |
| `ct_70_global_mv3_genetic.yaml` | global-covariate contrast |
| `ct_51_noglobal_mv2_matched.yaml` | 5.1 on the same specification |
| `t1t2_70_noglobal_mv2_genetic.yaml` | T1w/T2w ratio, matched to the settled spec |

## Reproducing the report

Every table and figure in `docs/` is produced by script from committed data —
never from a run directory:

```bash
python tools/regen_report_tables.py     # spatial/covariance/site tables
python tools/regen_h2_tables.py         # heritability tables (bootstrap; slow)
python tools/regen_report_figures.py    # table-based figures
python tools/regen_brain_maps.py        # DK surface maps
python -m pytest tests/ -q              # 152 tests
```

This is enforced, not trusted: `tests/test_docs_provenance.py` checks that every
cited figure and table exists and that every figure on disk has a generator. It
exists because 11 figures were once drawn in ad-hoc cells and silently survived
corrections to the numbers underneath them.

**If you work on this repo through Claude Science**, note that editing a file on
disk does not update its artifact, and neither does committing — so a corrected
report can sit on disk while a stale copy is what gets read. Run
`make audit-artifacts` at the end of any session that edits deliverables.

## Related prior work

The transcriptional targets tested here are inputs to this project, not part of
it: **AHBA C1–C3** from this group's previous work, and the **snRNA-seq
maturation axis** in
[`transcriptional_maturation`](https://github.com/richardajdear/transcriptional_maturation).
