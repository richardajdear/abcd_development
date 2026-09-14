# genetic_analysis/ — the ABCD thinning-rate genetics, re-run on the 7.0 tabulation

**For an agent starting fresh on CSD3.** Written 2026-09-14. This directory
replaces `hpc/`, `hpc_v2/` and `hpc_v3/`, which are now under `legacy/` and
which you do not need to read: everything they established that still matters
is in §4 and §5 below. Where a legacy section number is cited it is for audit,
not because you have to go there.

**Status: NOT STARTED on the cluster.** Do not submit anything until the local
re-run is confirmed and pushed to `main` (README.md "Status"). The phenotype
export it produces is the input to every step here.

---

## 0. The one-paragraph version

The target phenotype is `global_slope`: each child's rate of cortical thinning
across adolescence, a random slope of thickness on age from a mixed model
fitted to 2–4 MRI scans per child (`src/abcd/`, `R/fit_lmm.R`). The hypothesis
is that genetic risk for schizophrenia (SCZ) and depression (MDD) predicts a
faster rate. Three cluster passes (v1 GCTA/fastGWA, v2 GENESIS + a four-method
PRS grid, v3 region-subset phenotypes) established one robust result — **SCZ
polygenic score → faster thinning, β ≈ −0.035 to −0.047 SD/SD, significant in
3 of 4 PRS methods in both the European and the pooled arm** — and one
suggestive one (MDD, same direction, significant only pooled). Everything was
bounded by the phenotype's low SNP heritability (LDSC h² z ≈ 1.2), so no GWAS
hit, no interpretable genetic correlation, and no gain from region selection.
**All of that was computed on the 6.0 tabulated tables**, mislabelled as 7.0.
The true 7.0 tabulation adds 3,520 six-year scans; the longitudinal sample
grows from 8,192 to 8,716 children, four-visit children from 1,830 to 3,005,
and the slope's effective N by 29 %. The genotype side is unchanged. This
directory re-runs every *phenotype-dependent* step on the new export and
reuses every *genotype-only* product from v2.

## 1. Why we are doing it again — the vintage problem

| | 6.0 tables (what every legacy result used) | 7.0 tables (now at `abcd-data-release-7.0/`) |
|:--|--:|--:|
| data freeze | February 2024 | 1 August 2025 |
| six-year thickness rows | 4,086 | 7,607 |
| scans entering the settled model | 23,409 | 26,949 |
| children with ≥2 QC-passing scans | 8,192 | 8,716 |
| … of whom genotyped (in `gn_y_genrel`) | 8,082 | **8,596** |
| children with 4 scans | 1,830 | 3,005 |
| mean scans per child | 2.86 | 3.09 |
| median slope reliability | 0.157 | 0.211 |
| effective N for the slope (Σ reliability) | 1,361 | 1,752 |
| age centre of the model (years) | 12.44 | 12.80 |
| group thinning map, Spearman old vs new | — | 0.996 |
| `global_slope`, Pearson on the 8,185 shared children | — | 0.921 |

Source: `docs/vintage_comparison.csv`, written by `tools/compare_vintage.py`
from the two run directories. The map is the same map; the *subject-level*
phenotype moved (r = 0.92), because 2,478 of the shared children gained a
scan. That is the point: the PRS and heritability analyses spend their power
on subject-level precision, which is what improved.

Two consequences you must carry into the analysis:

1. **The 6.0-vintage phenotype files on CSD3 are now wrong inputs.**
   `$ABCD_HPC_ROOT/pheno/`, `legacy/hpc/work/pheno_allanc*/`,
   `legacy/hpc_v2/work/pheno_v3/` — do not point anything at them except to
   reproduce a legacy number.
2. **The age centre moved** (12.44 → 12.80 years), so `baseline_thickness`
   (the intercept phenotype) is thickness at a slightly older age than before.
   It remains the positive control; its h² and the SCZ locus-pool result
   should reproduce in kind, not to the third decimal.

Also changed in the 7.0 tables and relevant here: the ancestry PCs in
`ab_g_stc__gen_pc__01..32` were **recomputed** for all 11,670 genotyped
children (every value differs from 6.0), the genetic family/birth ids in
`gn_y_genrel` were re-coded, and ~3,600 pi-hat values changed. The export's
`covar_quant.txt` carries the 7.0 PCs; the v2 pipeline substitutes PC-AiR PCs
in GENESIS regardless, but `prs_assoc.R` and GCTA read the export's.

## 2. Where things are

### 2.1 In the repo (after `git pull` on CSD3)

```
genetic_analysis/
  README_HPC.md            this file -- append your run log to it (§8)
  config.sh                shared paths; sourced by every script.  NEW outputs go to
                           $ABCD_HPC_ROOT/results_70tab, phenotypes come from
                           $ABCD_HPC_ROOT/pheno_70tab; GDS + kinship are READ from results_v2
  config.local.sh.example  copy to config.local.sh (gitignored) and fill the CSD3 block
  00_check_inputs.sh       preflight -- run it before every submission
  run_all.sh               submits 03 -> 04 -> collect with afterok dependencies
  03_null_model.sbatch     GENESIS fitNullModel, one task per phenotype
  04_assoc.sbatch          GENESIS assocTestSingle, phenotypes x 22 chromosomes
  06_prs_family.sbatch     Fulker between/within-family PRS decomposition
  R/                       the GENESIS R code (01_make_gds, 02_kinship are NOT to be re-run)
  envs/genesis_env.yml     the R env spec (already built on CSD3, §2.2)
  setup/                   the PRS grid, LDSC/MAGMA drivers and helpers (table below)
  inputs/magma/            SCZ.genes.raw, MDD.genes.raw (+ .out): disorder-side MAGMA gene
                           results, staged for ahba_pls/.  Gitignored; copy from
                           $ABCD_HPC_ROOT/results/magma/ if missing
  work/                    gitignored scratch for this run; only *_summary.tsv and
                           table_*.tsv under work/results/ are tracked
```

Every script here is a **verbatim copy** of its legacy original with the
directory references rewritten (`hpc_v2/` → `genetic_analysis/`,
`hpc_v3/…` → `genetic_analysis/setup/…`, read-only legacy inputs → `legacy/…`).
Nothing has been executed since the copy. Treat the first run of each as a
smoke test: `DRY_RUN=1`, then one array task, then the array.

| copied file | origin | what it does | must-check before use |
|:--|:--|:--|:--|
| `config.sh`, `config.local.sh.example`, `00_check_inputs.sh`, `run_all.sh`, `03_null_model.sbatch`, `04_assoc.sbatch`, `06_prs_family.sbatch`, `R/*.R`, `envs/genesis_env.yml` | `legacy/hpc_v2/` | the GENESIS pipeline | `PHENO_DIR`, `OUT_V2`, `KIN_DIR` defaults (§2.1 header); `run_all.sh` step list still includes 01/02/05 — do not pass them |
| `setup/paths.sh` | `hpc_v2/work/setup/` | absolute path table for the PRS grid | `V2ROOT` now points at `…/genetic_analysis`; its `work/inputs/` symlink tree must be **recreated** (`ln -s` each entry of `legacy/hpc_v2/work/inputs/`), and `PHENO`/`COVQ`/`COVC`/`MANIF` must point at the new export |
| `setup/prs_final.sbatch`, `setup/collect_final.py`, `setup/collect_prs_tables.py`, `setup/standardise_within_ancestry.py`, `setup/within_ancestry_all.sbatch` | `hpc_v2/work/setup/` | four-method × trait-arm PRS grid and its tables | scoring is phenotype-independent (§3); only the `prs_assoc.R` layer must re-run — read the script and skip or short-circuit the weight/score steps |
| `setup/normalise_gwas.py`, `setup/make_grm_famid.sh` | `hpc_v2/work/setup/` | sumstat normaliser (not needed unless a GWAS is added); FID-as-family-id GRM view for GCTA | `make_grm_famid.sh` must be re-run on the new export: it rewrites `.grm.id` column 1 in the GRM's row order |
| `setup/eur_magma_ldsc.sbatch`, `setup/eur_prio_gsa.sbatch` | `hpc_v2/work/setup/` | EUR-arm MAGMA + LDSC; prioritised gene-set panel | they call `legacy/hpc/04_magma.sbatch`, `05_ldsc_rg.sbatch`, `work/prioritised_gsa.sbatch` with `MAGMA_DIR`/`LDSC_DIR`/`GWAS_DIR` overrides — point those at `results_70tab/` (fresh directories, §4 rule 9) |
| `setup/01_prs_ct.sbatch`, `02_ldsc_magma.sbatch`, `03_prs_paired.sbatch`, `04_collect.sbatch`, `05_prio_gsa.sbatch`, `prs_paired_delta.py` | `legacy/hpc_v3/` | C+T PRS over an existing score dir; LDSC+MAGMA into fresh dirs; paired Δβ / min-p permutation / random-region null; the collector with `TMPDIR` on rds | written for 9 phenotypes incl. the v3 four; drop those names |
| `setup/align_export.py` | `hpc_v3/align_export_v3.py` | fixes the three silent export mismatches (§5 item 1) | it **copies v2's covariate files** — that was right when PCs were identical; now the 7.0 export has its own PCs, so keep the export's covariates and only fix IID spelling, subset and scale |

`tools/prs_assoc.R` (population PRS association, `(1 | family_id)` mixed model)
is unchanged and stays in `tools/`.

### 2.2 On CSD3

Repo: `~/rds/hpc-work/abcd_development` (`/home/rajd2/rds/hpc-work/abcd_development`).
Working root: `ABCD_HPC_ROOT=~/rds/hpc-work/ABCD`. `/home` is nearly full — never
write there. Account `VERTES-SL3-CPU`, partition `icelake` (every header says so;
the cluster default `cclake` once cost 14 h).

**After pulling this commit, the untracked work trees are stranded at the old
paths** — git does not move ignored files. Run once:

```bash
cd ~/rds/hpc-work/abcd_development
for d in hpc hpc_v2 hpc_v3; do
  [ -d "$d/work" ] && mkdir -p "legacy/$d" && mv "$d/work" "legacy/$d/work"
  [ -d "$d/slurm" ] && mv "$d/slurm" "legacy/$d/slurm"
  [ -f "$d/config.local.sh" ] && mv "$d/config.local.sh" "legacy/$d/config.local.sh"
  rmdir "$d" 2>/dev/null || true
done
```

Then check that `legacy/hpc/work/envs/abcd` (micromamba, Python 3.11 + pandas/
pyarrow) and `~/rds/hpc-work/envs/genesis` (R 4.4.3, GENESIS 2.36.0, built from
`envs/genesis_env.yml`) still activate; `tools/hcp_extract.sbatch` reads
`ABCD_CONDA_ENV` if the first has moved.

| asset | path | notes |
|:--|:--|:--|
| 7.0 array genotypes | `/rds/project/rds-CeXlNYOYMxw/Data_Genetics/genotype_microarray/smokescreen/merged_chroms` | 11,670 × 515,228, hg19, FID = IID. Used for kinship only |
| imputed (v1 conversion) | `$ABCD_HPC_ROOT/genotype_imputed/imp_union_chr{CHR}` (22 filesets, rsID-keyed, 25.95 M variants); PRS fileset `abcd_imp_prs` (7.07 M) | association scan and scoring |
| GDS | `$ABCD_HPC_ROOT/results_v2/gds/` (73 GB) | reuse |
| kinship / PC-AiR / PC-Relate | `$ABCD_HPC_ROOT/results_v2/kinship/` (`pcair_pcs.tsv`, `pcrelate.rds`, sparse matrix), EUR arm `kinship_eur/` | reuse; density 0.004, largest block 726 |
| dense imputed GRM | `$ABCD_HPC_ROOT/results/grm_imp_pooled/abcd_imp` (+ family-id view via `make_grm_famid.sh`) | for GCTA GREML |
| PC-AiR unrelated set | `Data_Genetics/genotype_microarray/genesis/unrelateds_individuals.txt` (8,181) | the REML keep list |
| EUR arm | `$ABCD_HPC_ROOT/keep/eur_anchor.keep` (5,656 genotyped) | 4,116 of the old 8,082 were EUR; recount on the new export |
| strata | `$ABCD_HPC_ROOT/keep/strata_k4.tsv` | optional `group.var` |
| PRS score profiles | C+T: `$ABCD_HPC_ROOT/results/prs_imp/score_<DIS>_<thr>.profile`; PRS-CS / SBayesR / SBayesRC scores under `$ABCD_HPC_ROOT/results_v2/` per `setup/paths.sh` | **reuse** — a score never touches our phenotype |
| discovery sumstats (normalised) | `legacy/hpc_v2/work/inputs/gwas/` (symlinks to `$ABCD_HPC_ROOT/sumstats/`) | SCZ primary/european, MDD div/eur, ASD, ALZ (Wightman, ± APOE), ALZ_IGAP (Kunkle), EA |
| references | `magma/reference_data/g1000_eur`, `ldsc/eur_w_ld_chr` + `w_hm3.snplist`, `inputs/ref/{ldblk_ukbb_eur,sbayesr_ldm,sbayesrc_eigen,annot_baseline2.2.txt}`, `maf_union.tsv` | |
| disorder MAGMA gene results | `$ABCD_HPC_ROOT/results/magma/{SCZ,MDD}.genes.raw` | reuse; copies in `genetic_analysis/inputs/magma/` |
| binaries | `legacy/hpc/work/bin/{gcta64,magma,plink,plink2}`; GCTB via `inputs/bin/gctb`; PRScs via `inputs/bin/PRScs` | PLINK **1.9** is required for `--clump`/`--score` |
| LDSC env | `~/.conda/envs/ldsc` (Python 2.7 Broad LDSC) | pandas ≥ 3 breaks ldsc 2.0.1 |

## 3. What is reused and what must be re-run

**Reused unchanged (genotype-only):** GDS, KING/PC-AiR/PC-Relate kinship and
PCs, both GRMs, ancestry keep lists and strata, every PRS score profile and
every posterior-weight file, normalised sumstats, LD references, MAGMA gene
annotation and the disorder-side `.genes.raw`, LDSC LD scores.

**Must re-run (phenotype-dependent), in this order:**

| step | what | script | cost last time |
|:--|:--|:--|:--|
| 1 | phenotype export → aligned `$PHENO_DIR` (two FID conventions) | local `python -m abcd.gcta_export`; `setup/align_export.py` (adapted); `setup/make_grm_famid.sh` | minutes |
| 2 | preflight | `00_check_inputs.sh` | seconds |
| 3 | GENESIS null models, pooled and EUR | `03_null_model.sbatch` (array = phenotypes) | ~1 min/task |
| 4 | GENESIS association, pooled and EUR; collect with `TMPDIR` on rds | `04_assoc.sbatch` (5 × 22), `setup/04_collect.sbatch` | 3–31 min/task, ~2 h wall; collect 20 min |
| 5 | GCTA GREML h² on the dense imputed GRM, PC-AiR unrelated keep | `legacy/hpc/02_reml.sbatch` pattern with `PHENO`/`OUT` overrides | minutes/phenotype |
| 6 | PRS association layer: four methods × trait-arms, matched strata; within-ancestry standardisation; within-family; min-p permutation; controls | `setup/prs_final.sbatch` (association only), `setup/collect_final.py`, `06_prs_family.sbatch`, `setup/03_prs_paired.sbatch` | ~10 min each once scores exist |
| 7 | LDSC h² + rg (EUR arm), MAGMA gene + gene-property both directions, prioritised gene-set panel | `setup/02_ldsc_magma.sbatch`, `setup/05_prio_gsa.sbatch` | 3.5 h; 3 min |
| 8 | ahba_pls H3: PLS-signature gene-property on the new phenotype `.genes.raw` | `ahba_pls/FOLLOWUP_GENETICS.md` | minutes |

Do **not** re-run: 01 (GDS), 02 (kinship), any PRS *scoring*, the Zaitlen
two-GRM REML (§5 item 11), LAVA (§5 item 12), Experiment A's phenotypes.

## 4. Rules that cost weeks to learn — do not re-derive

1. **IDs.** ABCD spells a child three ways (`NDAR_INVxxxxxxxx`,
   `sub-NDARINVxxxxxxxx`, `sub-xxxxxxxx`). Every `.fam`, `.grm.id` and
   `.profile` on CSD3 uses `sub-xxxxxxxx`. GENESIS normalises to the 8-char
   token; **GCTA and PLINK match FID+IID literally and return zero subjects
   without an error.** Join on the token; write the `.fam` spelling.
2. **Two FID conventions, two files.** The genotype `.fam` has FID = IID, so
   the real family id must come from the phenotype export. GCTA wants
   FID matching the GRM (`make_grm_famid.sh` makes a family-id view of the
   `.grm.id` **in the GRM's row order** — never reorder it); `prs_assoc.R` and
   the Fulker step want FID = family_id. Keep both exports side by side.
3. **Scale.** Every β in this project is per phenotype SD. `global_slope` has
   native variance ~1e-6; z-score every column over the
   phenotyped-and-genotyped analysis set (8,596 now), *after* subsetting.
4. **Ancestry matching decides the stratum.** Multi-ancestry discovery GWAS
   (SCZ primary, MDD div) → pooled target; European discovery (SCZ european,
   MDD eur, ASD, both ALZ, EA) → EUR target. European GWAS → pooled target is
   *confounded* (SE inflates with threshold density); it manufactured a false
   ASD "hit" once. `collect_final.py` emits a `matched` flag — read only
   matched cells.
5. **Allele frequencies come from the target sample** (`maf_union.tsv`), never
   borrowed from a reference panel; borrowing silently dropped 12 % of ASD SNPs
   and created a spurious specificity gradient.
6. **Multiplicity over nested C+T thresholds is family-block permutation of
   min-p**, not Bonferroni (8 thresholds ≈ 4.1 effective tests). Or use a
   continuous-shrinkage score and avoid it.
7. **Controls are mandatory and do not discriminate on their own.** Run ASD,
   ALZ (± APOE, and Kunkle), and EA alongside SCZ/MDD from the start. ALZ is
   itself associated with `global_slope` in the population arm (~80 % of SCZ's
   magnitude) and attenuates within families; SCZ and MDD do not. EA runs the
   opposite direction. Report the panel, not the target alone.
8. **Report λ_GC with every scan.** v1's 20 "hits" at λ = 1.107 were
   stratification; v2 at λ ≈ 1.01–1.04 has zero. A hit count without λ is
   uninterpretable.
9. **Fresh output directories for LDSC/MAGMA.** Their collectors rewrite
   `*_summary.tsv` to cover exactly the manifest they were handed; pointing
   them at a legacy directory overwrites the published table without error.
   `results_70tab/` is the new root; symlink disorder-side intermediates in.
10. **`sbatch --export` splits on commas**; pass cell lists positionally.
    `DRY_RUN=1` once *submitted* a whole pipeline of no-ops. Verify outputs
    (row counts, `n_chr == 22`), never `sacct` state — 13 jobs once reported
    COMPLETED having produced nothing. Collectors that `fread(cmd=)` need
    `TMPDIR` on rds, not node-local `/tmp`.
11. **The Zaitlen two-GRM "pedigree h²" is not obtainable from this design**:
    the bK GRM's diagonal is PC-AiR PC1 (r = −0.99), and deleting 97 % of
    relative pairs moved `h2_ped` by 0.0000. The canonical h² is v1's GCTA
    GREML on the dense imputed GRM over the PC-AiR unrelated set. Do not
    report `reml_zaitlen*`.
12. **LAVA local rg is anti-conservative** (11 % of tests at p < 0.05 genome-
    wide; the ZEB2 claim was retracted). Skip it unless you permutation-
    calibrate it.
13. **rg needs h² z ≳ 4 on both sides.** `global_slope` sits at ~1.2, so rg is
    *uninformative, not null*; say so and report the z. `baseline_thickness`
    (z ≈ 4.5) is the one phenotype where rg means something — it is the
    positive control for the LDSC step. PRS is the primary disorder test for
    exactly this reason: the discovery GWAS carries the statistical burden.
14. **The `age_c` covariate was silently absent from every legacy PRS model**
    (the export names it `baseline_age`). It is now derived in
    `tools/prs_assoc.R`, `R/06_prs_family.R`, `R/07_prs_conditional.R`; the
    effect on β was ≤ 0.006. Keep it in.
15. **Region selection is not a lever.** Top-ΔCT and top-C3 subset means share
    80–90 % of variance with `global_slope`, are *less* heritable, and the C3
    set is weaker than 80 % of random 8-region sets. Projections are worse.
    Do not build new phenotypes this way; the whole-cortex mean stays primary.
16. **Never commit per-subject genotype-derived files** (`.profile`, `.sscore`,
    `pcair_pcs.tsv`, per-individual scores). Summary tables only. The repo is
    public.
17. **`scratch/hpc_test*` are synthetic fixtures with planted effects.** A
    phenotype file there loads cleanly and gives plausible wrong answers.

## 5. Legacy numbers to compare against (6.0-vintage phenotype, n = 8,082 / EUR 4,116)

| readout | value | file |
|:--|:--|:--|
| SCZ PRS → `global_slope`, EUR: C+T / PRS-CS / SBayesR / SBayesRC | −0.047\* / −0.027 / −0.035\* / −0.036\* SD/SD | `legacy/hpc_v2/work/results_v2/prs_final/table_main.tsv` |
| same, pooled (within-ancestry standardised) | −0.033\* / −0.015 / −0.029\* / −0.025\* | same |
| MDD PRS → `global_slope`, pooled: C+T / SBayesR | −0.044 / −0.079\* | same |
| ASD (EUR, 4 methods) | +0.030 / +0.016 / −0.005 / −0.004 (all null) | same |
| ALZ, EUR, C+T | −0.049\*; noAPOE −0.044\*; Kunkle −0.040 | same |
| min-p permutation, SCZ → `global_slope`, EUR | p_perm 0.007 (Bonferroni 0.016) | `legacy/hpc_v3/prs_tables/prs_minp_permutation.tsv` |
| within-family SCZ, EUR (686 pairs) | β_W −0.057 (p 0.30) vs β_B −0.046 (p 0.004), p_diff 0.85; ~18 % power | `legacy/hpc_v2/work/results_v2/prs_family/` |
| GCTA GREML h² (n 5,649 unrelated): baseline / global_slope / PC2 / PC1 / PC3 | 0.246±0.049 / **0.137±0.046** / 0.161 / 0.075 / 0.041 | `legacy/hpc/work/results/reml_imp_pooled/reml_summary.tsv` |
| LDSC EUR: h² z baseline / global_slope; rg(SCZ, global_slope) | 4.49 / 1.16; −0.149 ± 0.109 (−0.168 ± 0.113 with the European SCZ file) | `legacy/hpc_v2/work/results_v2/ldsc_eur*/` |
| GWAS: hits / λ_GC, pooled v2 | 0 / 1.014 (`global_slope`), 1.036 (`baseline_thickness`) | `legacy/hpc_v2/work/results_v2/assoc/gwas_summary.tsv` |
| MAGMA `SCZ_locus_pool` → baseline / global_slope, EUR | β 0.187, p 1.1e-4 / β 0.124, p 7.4e-3 | `legacy/hpc_v2/work/results_v2/magma_prio_eur*/` |
| AHBA C1–C3 gene-property on the ABCD GWAS | null in the matched arm | same |

\* threshold-adjusted p < 0.05. The direction to expect: same signs, smaller
SEs. The SCZ result surviving with a tighter interval is the headline; MDD
crossing into significance under more than one method would be new. Anything
that flips sign is a bug until proven otherwise.

## 6. Step-by-step

### Step 0 — bring the data and the export over (from the laptop)

The tabulated tables are not on rds (the `derivatives/tabulated/` copy is 6.0).
They are gitignored in the repo, so they travel by rsync into the repo's
release directories on CSD3, where `src/abcd/paths.py` finds them:

```bash
cd ~/Git/abcd_development
rsync -avh --progress --exclude '.DS_Store' \
  abcd-data-release-7.0/ login-q-1.hpc.cam.ac.uk:rds/hpc-work/abcd_development/abcd-data-release-7.0/
rsync -avh --progress --exclude '.DS_Store' \
  abcd-data-release-6.0/ login-q-1.hpc.cam.ac.uk:rds/hpc-work/abcd_development/abcd-data-release-6.0/
```

The first command also carries `processed/hcp/` (64 MB, the HCP-MMP table) so
the cluster copy matches the laptop. If `abcd-7.0/` still exists on CSD3 from
the HCP extraction, remove it after checking it holds nothing else — two
release directories for one release would let `paths.release_dir` pick either.

Then the export (written locally by `tools/rerun_local.sh`; 8,716 rows,
FID = family id, IID = `NDAR_INV…`, 10 ancestry PCs from the 7.0 `ab_g_stc`):

```bash
rsync -avh out/thickness_dsk_70_139406217085/gcta_inputs/ \
  login-q-1.hpc.cam.ac.uk:rds/hpc-work/ABCD/pheno_70tab/raw/
```

### Step 1 — align the export (on CSD3)

Produce, under `$ABCD_HPC_ROOT/pheno_70tab/`:

- `phenotypes_gcta.txt` with IID in the `.fam` spelling, **subset to the
  genotyped** (join on the 8-char token against `merged_chroms.fam`; expect
  8,596), **then** z-scored per column; FID = family id. Plus `covar_quant.txt`
  (baseline_age, n_visits, PC1–10 from the export — *not* copied from v2),
  `covar_categorical.txt` (sex, site), `phenotype_manifest.tsv`.
- a sibling `pheno_70tab_fidiid/` with FID = IID for GCTA, and the family-id
  GRM view from `setup/make_grm_famid.sh`.
- Record: n, native-unit mean/SD per phenotype (so a β can be put back in
  mm/yr), number of multi-member families, EUR-arm n.

Adapt `setup/align_export.py`: keep its IID/subset/scale logic, drop the
covariate copy. Sanity: the five columns must correlate > 0.9 with the legacy
export on shared children (`global_slope` r ≈ 0.92 locally; `baseline_thickness`
higher), not reproduce it.

### Step 2 — preflight

`cp config.local.sh.example config.local.sh`, fill the CSD3 block (the commented
values there are 2026-08 and must be re-verified against the filesystem), then
`bash 00_check_inputs.sh`. It checks `.bed` integrity, the 22 imputed filesets,
the ID join, the manifest, and that the R has GENESIS.

### Step 3–4 — GENESIS scan, both arms

`bash run_all.sh 03 04` (pooled). For the EUR arm export
`KIN_DIR=$V2_LEGACY/kinship_eur`, `NULL_DIR`/`ASSOC_DIR` suffixed `_eur`, and the
EUR keep list, as `legacy/hpc_v3/README_HPC.md` §Status did. Collect with
`setup/04_collect.sbatch` (`TMPDIR` on rds, 48 G). Check `gwas_summary.tsv`:
five rows per arm, `n_chr == 22`, λ_GC in 0.99–1.05, n ≈ 8,596 / EUR n.

### Step 5 — SNP heritability

GCTA GREML, dense imputed GRM, PC-AiR unrelated keep, covariates sex + site +
baseline_age + n_visits + 10 PCs — the v1 canonical design
(`legacy/hpc/02_reml.sbatch` with `PHENO`/`COVAR_*`/`OUT` overrides). Expect the
unrelated n to rise from 5,649 with the sample. Report h² ± SE for the five
phenotypes next to the legacy row.

### Step 6 — PRS

Scores exist. Run only the association layer: `tools/prs_assoc.R` per method ×
trait-arm (via `setup/prs_final.sbatch` if it can be told to skip scoring;
otherwise call `prs_assoc.R` directly over the existing score directories),
`setup/standardise_within_ancestry.py` for the pooled arm, `setup/collect_final.py`
→ `work/results/prs_final/table_{main,all,family}.tsv`. Then `06_prs_family.sbatch`
(Fulker), `setup/03_prs_paired.sbatch` (min-p permutation; the paired-Δβ and
random-region parts are for the retired v3 phenotypes and can be left out), and
the control panel (ASD, ALZ, ALZ_noAPOE, ALZ_IGAP, EA) in every table.

### Step 7 — LDSC and MAGMA

`setup/02_ldsc_magma.sbatch` with the manifest reduced to the five phenotypes and
output under `results_70tab/{ldsc_eur,magma_eur}/`; then `setup/05_prio_gsa.sbatch`
for the SCZ locus pool / MDD pool panel into `results_70tab/magma_prio_eur/`.
Report h² z first; rg only where z supports it.

### Step 8 — the imaging-transcriptomics tie-in

`ahba_pls/FOLLOWUP_GENETICS.md` H3: MAGMA gene-property of the new phenotype
`.genes.raw` on `ahba_pls/hpc/lead_pls2_gene_covar_entrez.txt` (three DS columns
marginal, and conditioned on `AHBA_C3`). Expect null; report it as such.
(The PLS weights themselves are being refreshed locally on the 7.0 maps — use
the committed file at the time you run.)

### Later — HCP-MMP

The 360-parcel thickness table exists (`abcd-data-release-7.0/processed/hcp/`)
but covers 24,921 of 33,825 sessions until the re-parcellation
(`processed/hcp/sessions_to_reparcellate.txt`, 8,874 sessions) is done by
R. Romero-Garcia. When a DK-vs-HCP genetics comparison is wanted, the export
comes from `configs/ct_70_hcp_noglobal_mv2.yaml` and follows the same steps.

## 7. Definition of done

1. A results table in `work/results/summary_70tab.tsv` (tracked) with one row
   per (readout, phenotype, arm, method): estimate, SE, p, n, and the legacy
   value alongside — every row in §5 has a partner.
2. λ_GC and n per scan; h² ± SE with the unrelated n; LDSC h² z per phenotype
   and the explicit rg verdict.
3. The PRS panel with all five controls, matched strata only, permutation-
   corrected, plus the within-family decomposition with its power caveat.
4. A dated run log appended to this file (§8): job ids, scripts, n, wall time,
   failures with diagnosis — the same contract the legacy READMEs kept.
5. A one-paragraph statement of what changed against §5 and what it licenses.
6. Nothing per-subject committed; `git status` clean of `work/` except the
   summary tables.

## 8. Run log

*(empty — nothing has been run on the cluster on the 7.0 tabulation)*
