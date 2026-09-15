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

### 2026-09-14 — bring-up, and step 1 (phenotype export) submitted

Pulled `f9a8806`. Steps 0 and 1 are done or submitted; **nothing in steps 2–8
has been run yet.** Read the two decisions marked **DEVIATION** before
continuing — both change what a later step should point at.

**§2.2 migration — done, but not with the snippet as written.** The 7.0 commit
already tracks summary tables under `legacy/hpc/work/` and
`legacy/hpc_v2/work/`, so `mv "$d/work" "legacy/$d/work"` would have nested the
untracked tree at `legacy/hpc/work/work`. It was merged entry by entry instead
(move-only; a file identical to a tracked copy is parked, never deleted):
**1,938 entries moved, 0 conflicts, 0 duplicates**, then the empty skeletons
removed. `hpc/`, `hpc_v2/` and `hpc_v3/` are gone and `git status` is clean.

Three things the migration broke, all now fixed — check these first if a path
error appears:

1. The repo-root `out` symlink pointed at `hpc/work/out`; repointed to
   `legacy/hpc/work/out`.
2. **152 symlinks broke** — every one an absolute path into the pre-move
   `hpc/`, `hpc_v2/` or `hpc_v3/`. All repaired, verified 0 broken. This is the
   migration's most dangerous side effect and is worth re-checking after any
   further move, because most of the casualties are **assets §3 says to
   reuse**, and a missing score file is the kind of thing that surfaces as a
   quietly smaller N rather than an error:
   - 16 of the 33 entries in `legacy/hpc_v2/work/inputs/` (the tree every
     `setup/` script resolves through);
   - 134 under `legacy/hpc_v2/work/results_v2/`, including the `prs_ct_v3`
     C+T score profiles and weights, every `prs_final/*/_zanc/` within-ancestry
     weight file across all four PRS methods, the `ldsc_eur_v3` munged
     sumstats, and `legacy/hpc_v2/work/grm/abcd_imp_famid.grm.*`;
   - 2 needed a second pass, because their target was itself a link repaired
     later in the first.

   A full scan of `legacy/` finds 203 further broken links, all inside the
   conda package cache at `legacy/hpc/work/mamba/pkgs/`. Those are ordinary
   relative links between unextracted package files, were broken before the
   migration, and are not on any analysis path.
3. The relocated `legacy/hpc/work/envs/abcdR` R 4.5.3 install (this is the
   **lme4** env, distinct from the GENESIS env) has its prefix compiled in.
   `lib/R/bin/R` and `etc/Renviron` were rewritten in place (originals kept as
   `*.pre70bak`), and `Rscript` — an ELF that ignores `R_HOME` — was replaced
   by a shell shim forwarding to `R`, with the binary kept as
   `Rscript.elf.pre70bak`. The micromamba `envs/abcd` Python 3.11 relocated
   cleanly; `~/rds/hpc-work/envs/genesis` (R 4.4.3, GENESIS 2.36.0,
   SeqArray 1.46.0) is outside the repo and was unaffected.

**Roots — §2.2's `$ABCD_HPC_ROOT/...` paths do not exist as written.** There is
no single working root on this filesystem: v1's was `hpc/work` and v2's was
`hpc_v2/work`, both inside the repo (now under `legacy/`). `~/rds/hpc-work/ABCD`
is an unrelated 2024–25 tree with none of these assets and is **not** used.
`genetic_analysis/config.local.sh` gives this run its own root
(`genetic_analysis/work`) and names every reused asset against its real legacy
root. Three values in `config.local.sh.example`'s CSD3 block were wrong and are
corrected there: the imputed template is `imp_union_chr{CHR}` (not
`abcd_imp_chr{CHR}`), the dense GRM is `results/grm_imp_pooled` (not
`results/grm_imp`), and the keep lists are under `results/ancestry/` (not
`keep/`). Every reused asset in §2.2 was confirmed present: GDS (23 files),
kinship and `kinship_eur`, 22 `imp_union` filesets, the dense GRM, 89 PRS score
files, sumstats, references, and `gcta64` / `plink` 1.9.0-b.7.7 / `plink2` /
`magma` 1.10.

**A `make_grm_famid.sh` trap, fixed in config.** The script does
`DST="$GRM_FULL"`, so pointing `GRM_FULL` at v1's GRM would have rewritten v1's
own `.grm.id` in place, against the read-only policy. `GRM_FULL` is now the new
view `work/grm/abcd_imp_famid` and `GRM_SRC` is v1's. `FAMMAP` points at the
**new export's** `family_map.tsv` (written by `align_export.py`), not v1's —
`gn_y_genrel`'s family ids were re-coded in 7.0, so v1's map is a wrong input.
`V1_ROOT` is also defined, which that script expects and `config.sh` does not
set.

**Step 0 — the release tables are already on CSD3 and are the true 7.0.**
Verified by the row count that separates the vintages: `mr_y_smri__thk__dsk.tsv`
has **7,607** six-year rows (6.0 has 4,086), 33,794 rows total. `abcd-7.0/`
still exists beside `abcd-data-release-7.0/`, but `_RELEASE_DIR_PATTERNS` tries
`abcd-data-release-{r}` first so the correct directory always wins. It was left
in place because it also holds the HCP area/volume tables and the long
anatomical-stats parquet, which the release directory does not; the HCP arm is
"Later" in this plan.

**DEVIATION 1 — the export was regenerated on CSD3, not rsynced.** The export
sitting at `out/thickness_dsk_70_139406217085` was the **6.0-vintage** one
(8,192 subjects, assembled 2026-08-10, still carrying the retired v3
phenotypes). The laptop's 7.0 run has the *same run id* — the config hash is
unchanged because only the underlying tables moved — so it would have been
overwritten silently. It is archived at
`out/legacy_6.0_tabulated/thickness_dsk_70_139406217085`, mirroring the
laptop's own convention in `docs/vintage_comparison.csv`.

The laptop is not reachable from here, and the 7.0 tables are, so
`step1_export_pheno.sbatch` (**job 35547946**, queued at the time of writing)
re-runs `assemble → fit_lmm → phenotype → gcta_export` for
`ct_70_noglobal_mv2_genetic` on CSD3. This is only legitimate if it reproduces
the local run, so **`step1_verify_export.py` gates it**: sample counts, family
count, per-visit scan counts, the age centre and the 5-phenotype export shape
must match `docs/vintage_comparison.csv` exactly, and `global_slope` must
correlate with the archived 6.0 export at r ≈ 0.9212 over the 8,185 shared
children. lme4 BLUPs may differ in the last digits between machines, which is
why the sample construction is checked exactly and the phenotype by
correlation. **If that gate fails, discard the export and rsync the laptop's
(Step 0) — do not proceed to steps 2–8.**

**DEVIATION 2 — `setup/paths.sh` now writes to `results_70tab`.** It had
`RES="$V2ROOT/work/results_v2"`, which under the new `V2ROOT` is a fresh empty
directory but is misleadingly named; it now follows `$OUT_V2`, per §4 rule 9.
`work/inputs/` was rebuilt as §2.1 requires: 33 reused entries relinked to
their real targets plus the 4 reference files, with `pheno/*` pointing at
`work/pheno_70tab/` (dangling until job 35547946 lands — expected).

`.gitignore` gained a `work/results_70tab/` block mirroring the `work/results/`
contract, so summary and `table_*` files are trackable and nothing per-subject
can be. Note that §7 item 1's `summary_70tab.tsv` sits at the root of the
results directory, which the `*/*_summary.tsv` rule does not reach; it is named
explicitly.

New scripts, named `step1_*` so they cannot be confused with `run_all.sh`'s
`01` (GDS, which must not be re-run): `step1_export_pheno.sbatch` and
`step1_verify_export.py`. `setup/align_export.py` was rewritten per §2.1: it
keeps the IID-spelling, genotyped-subset and z-score logic, **drops the v2
covariate copy** (the 7.0 PCs differ from 6.0's, so the export's own covariates
are the right ones), writes both FID conventions side by side, and records
native mean/SD per phenotype in `align_report.tsv` so a standardised β can be
put back into mm/yr.

**A `make_grm_famid.sh` path bug, fixed.** It sourced
`$(dirname $0)/../../config.sh`, which was right at `hpc_v2/work/setup/` but
resolves to the repo root from `genetic_analysis/setup/`; it now sources
`../config.sh`. No other copied script has the same pattern.

**`setup/align_export.py` was validated before the new export exists**, by
running it on the archived 6.0 export: it reproduces v2's published
`pheno_allanc_prs/phenotypes_gcta.txt` **byte-exactly** (max |diff| = 0 across
all five settled z-scored columns), the same FIDs, and the same 8,082-row
family map, and independently recovers three published numbers — 8,192
phenotyped → 8,082 genotyped, EUR arm 4,116, and native mean `global_slope`
−0.01892 mm/yr. The 6.0 export carried **no ancestry PCs** (4 covariate
columns), which is what made v3 copy v2's. That is fixed for 7.0: the release
adapter reads 10 PCs for 11,860 children from `abcd-data-release-7.0/g/stc/`,
so the new export carries its own. Checked directly, **100 % of PC cells differ
between 6.0 and 7.0**, confirming §1 and settling that the v2 covariate copy
had to go.

**Chained, not left hanging.** `step1_align_and_preflight.sbatch`
(**job 35548496**) is queued `--dependency=afterok:35547946` and runs
`step1_verify_export.py` → `align_export.py` → `make_grm_famid.sh` →
`00_check_inputs.sh`, with `set -e` so the verification gate stops the chain.
It deliberately stops at the preflight: steps 3–8 are hours of compute and
should be launched after reading its log. The icelake queue was 365 pending /
651 running at submission, so expect a wait.

**Job 35547946 FAILED at 00:01:54; fixed and resubmitted as 35558401**
(with 35558402 chained behind it; 35548496 was cancelled by its own `afterok`,
which is the dependency working).

`assemble` succeeded, and its manifest was checked against
`docs/vintage_comparison.csv` on the spot: **all 12 assemble-stage quantities
match the local run exactly** — 8,716 subjects, 26,949 scans in the model,
7,285 families, the 2/3/4-visit split 2,204 / 3,507 / 3,005, all four per-visit
scan counts, and the age centre to four decimals (12.7973). Sample construction
on CSD3 is therefore identical to the laptop's, which is the part of
DEVIATION 1 that could have been wrong. What the full gate still has to settle
is the lme4 side: the BLUPs, via the r ≈ 0.9212 check.

The `PY` fix was also verified directly rather than by resubmitting and
hoping: with `PY` exported, R's `resolve_run_dir` returns the right run
directory with status 0.

`fit_lmm.R` then failed. It resolves the run directory by shelling out to
`python -m abcd.run_dir` via `Sys.getenv("PY", unset = "python")`, and the
sbatch set `PY` **without exporting it**, so R fell back to bare `python` and
got CSD3's system Python 3.7, which has no `typing.Literal` and cannot
`import abcd`. `R/fit_lmm.R` line 66 warns about exactly this. The sbatch now
exports `PY` and `RSCRIPT` as absolute paths and passes `--run-dir` explicitly
so the resolution cannot happen at all. Worth remembering: an interpreter that
is right for a direct call can still be wrong for a *nested* one, and
`set -u` does not catch an unexported variable.

**Job 35558401 COMPLETED. Step 1 is done and DEVIATION 1 is settled.**

`step1_verify_export.py`: **all 18 checks pass**, including the two that
decide it — `subjects_shared` 8,185 and `global_slope` r vs the 6.0 export
**0.9212**, matching `docs/vintage_comparison.csv` to four decimals. lme4 on
CSD3 therefore reproduces the laptop's BLUPs, not just its sample. The export
carries its own **10 ancestry PCs** (120 subjects NA, which GCTA drops), so the
v2 covariate copy is gone for good.

`setup/align_export.py` on the verified export:

| | |
|:--|--:|
| phenotyped | 8,716 |
| phenotyped **and genotyped** (the analysis set) | **8,596** |
| families / multi-member families | 7,198 / 1,350 |
| EUR arm (recounted, was 4,116 on 6.0) | **4,308** |
| `global_slope` native mean, SD (mm/yr) | −0.018696, 0.0013222 |

The 8,596 is exactly what §2.2 predicted. Native mean/SD for all five
phenotypes are in `work/pheno_70tab/align_report.tsv`: multiply a standardised
β by the SD there for mm/yr. Both FID conventions were written
(`pheno_70tab`, `pheno_70tab_fidiid`), plus `family_map.tsv`.

`setup/make_grm_famid.sh`: 11,670 rows, 8,596 FIDs replaced with family ids,
IID column and order preserved, and v1's own `.grm.id` verified untouched. The
4 `work/inputs/pheno/` links now resolve; 0 broken in that tree.

**Step 2 PASSED (job 35558402).** `preflight PASSED`, with every gate green:
the array fileset and **all 22 imputed filesets** bed/bim/fam-consistent
(25.95 M variants, n = 11,670 throughout), **ID token join 8,596 of 8,596**,
1,350 multi-member families, the manifest covering 5 phenotypes in 5 exported
columns, and the GENESIS R carrying GENESIS/SNPRelate/GWASTools plus the
lme4 stack step 06 needs.

**Step 3 started as the §2.1 smoke-test ladder, not the full array.**
`DRY_RUN=1 bash run_all.sh 03 04` was inspected first and is correct: array
1–5 for the null models, 1–110 for the association scan (5 phenotypes × 22
chromosomes), a `v2_collect` pass chained `afterok` behind the whole array, and
every path resolving into `work/results_70tab/` — nothing points at a legacy
output directory (rule 9). Note the DRY_RUN trap in rule 10 is already fixed in
this copy: `run_all.sh`'s own `submit()` honours it and pins `DRY_RUN=0` into
real submissions.

**Job 35558883** is `03_null_model.sbatch` with `--array=1` only —
`baseline_thickness`, the positive control — as the single-task rung of the
ladder before the array.

**The smoke test earned its keep — it caught a silent design change.**
35558883 completed and wrote a converged null model with the right shape
(n = 8,596, 14 covariates = sex + site + baseline_age + n_visits + 10 PCs).
Checked against v2's published summary rather than just eyeballed:

| | 6.0-vintage v2 | first 7.0 attempt | after the fix |
|:--|--:|--:|--:|
| n | 8,082 | 8,596 | 8,596 |
| `group_var` | none | **stratum** | none |
| `varcomp_kin` | 0.8416 | 0.8419 | — |
| `prop_kin` | 0.879 | **0.635** | — |

`varcomp_kin` agreeing to three decimals across a 6 % larger sample is a good
sign. But `group_var` was **stratum**, because `config.local.sh` had set
`STRATA_FILE`. §2.2 lists strata as an "optional `group.var`" and v2 ran
without it; switching it on makes GENESIS fit a separate residual variance per
ancestry stratum, which is why `prop_kin` moved from 0.879 to 0.635. That is a
different model, and §5 is only a benchmark if the design matches. `STRATA_FILE`
is now left unset, with the reasoning recorded in `config.local.sh` and the
line kept commented for a stratified sensitivity run. **Job 35559026** is the
re-run; it should come back with `group_var = none` and `prop_kin` near 0.879.

The general lesson, since this class of thing has cost this project weeks
before: a config value that is merely *available* is not thereby *wanted*. The
optional knobs in §2.2 default to off because v1/v2 ran with them off.

35559026 came back `group_var = none`, `varcomp_kin` 0.8460 against v2's
0.8416 and `prop_kin` 0.8852 against 0.8793, with `varcomp_resid` **lower**
(0.1098 vs 0.1155) — the design now matches and the larger sample shows up as
less residual variance, which is the direction §5 predicts.

**Step 3–4 pooled arm SUBMITTED**: `bash run_all.sh 03 04` gave
**35559124** (null, array 1–5) → **35559125** (assoc, array 1–110,
`afterok`) → **35559126** (`v2_collect`, `afterok` behind the whole array).
Roughly 2 h wall.

**The EUR arm recipe, worked out and checked against the filesystem** (do not
re-derive it):

```bash
KIN_DIR=$V2_LEGACY/kinship_eur \
NULL_DIR=$OUT_V2/nullmodel_eur ASSOC_DIR=$OUT_V2/assoc_eur \
  bash genetic_analysis/run_all.sh 03 04
```

`03_null_model.sbatch` takes both the sparse kinship and the PC-AiR PCs from
`$KIN_DIR` (lines 30–31), and neither 03 nor 04 has a keep-list argument: the
arm **is** the kinship's sample set, intersected with the phenotype file. So
`kinship_eur` (n = 5,656, 4,573 PC-AiR unrelated, density 0.00025, largest
block 5) is what makes it the EUR arm. **`GDS_DIR` must stay pooled** —
`results_v2/gds_eur/` holds only `array.gds`, no per-chromosome imputed GDS, and
GENESIS subsets to the null model's samples anyway.

**Benchmarks for reading the output** (from
`legacy/hpc_v2/work/results_v2/assoc*/gwas_summary.tsv`, 6.0-vintage, and note
`n_mean` is below the analysis n because of per-SNP MAF/missingness):

| | pooled | EUR |
|:--|--:|--:|
| `n_snps` | 9,413,281 | 7,469,278 |
| `n_mean` | 7,893 | 4,017 |
| λ_GC `global_slope` / `baseline_thickness` | 1.0139 / 1.0363 | 1.0001 / 1.0270 |
| genome-wide hits | 0 | 0 |

Expect the same shape with larger `n_mean`. **The EUR arm n is 4,308 phenotyped
and genotyped, not §5's 4,116** — that row is 6.0-vintage. Check `n_chr == 22`
and λ_GC in 0.99–1.05 before reading anything else (rules 8 and 10).

**All five pooled null models COMPLETED and verified** (35559124). Design
matches v2 exactly on every row — n = 8,596, 14 covariates, `group_var` none,
`inv_norm` FALSE — and the variance components carry the first real signal
that the 7.0 phenotype is better, not merely bigger:

| phenotype | `varcomp_kin` 6.0 → 7.0 | `varcomp_resid` 6.0 → 7.0 | `prop_kin` 6.0 → 7.0 |
|:--|:--|:--|:--|
| baseline_thickness | 0.8416 → 0.8460 | 0.1155 → **0.1098** | 0.879 → 0.885 |
| global_slope | 0.3400 → **0.3760** | 0.5415 → **0.4785** | 0.386 → **0.440** |
| slope_PC1 | 0.3131 → 0.3181 | 0.6165 → 0.5903 | 0.337 → 0.350 |
| slope_PC2 | 0.3446 → 0.3069 | 0.6244 → 0.6547 | 0.356 → 0.319 |
| slope_PC3 | 0.3598 → 0.3371 | 0.5688 → 0.5798 | 0.387 → 0.368 |

`global_slope` is the one that matters and it moves the right way on all three:
more variance carried by kinship, **12 % less residual**, and `prop_kin` up
from 0.386 to 0.440. That is what "the slope's effective N rose 29 %" should
look like at the null-model stage. The two PCs that go the other way are
second-order and not the primary; do not over-read them before the scan.

This is not yet h² — it is the sparse-kinship decomposition of the GENESIS null
model, which includes shared-family variance. The canonical h² is still step 5's
GCTA GREML on the dense imputed GRM over the PC-AiR unrelated set (rule 11).

**EUR null models submitted** as **job 35559228** (`03` only, not `04`): the
null step is cheap and independent, so it validates the EUR wiring — above all
the claim that the kinship *is* the arm — before committing another 110
association tasks. Its n should come back near 4,308 rather than 8,596.

**EUR null models COMPLETED and the wiring is confirmed** (35559228): **n =
4,308** on every phenotype, exactly the figure `align_export.py` predicted and
not the pooled 8,596 — so the arm really is the kinship's sample set, with no
keep-list flag involved. Design matches (14 covariates, `group_var` none).

| phenotype | `varcomp_kin` 6.0 → 7.0 | `varcomp_resid` 6.0 → 7.0 | `prop_kin` 6.0 → 7.0 |
|:--|:--|:--|:--|
| baseline_thickness | 0.8315 → 0.8341 | 0.1086 → 0.1098 | 0.884 → 0.884 |
| global_slope | 0.3099 → **0.3654** | 0.5141 → **0.4166** | 0.376 → **0.467** |
| slope_PC1 | 0.2555 → 0.3083 | 0.6387 → 0.5484 | 0.286 → 0.360 |
| slope_PC2 | 0.3403 → 0.3259 | 0.6110 → 0.5900 | 0.358 → 0.356 |
| slope_PC3 | 0.4131 → 0.4121 | 0.5174 → 0.4809 | 0.444 → 0.461 |

The pattern is the one to want: `baseline_thickness`, the positive control,
barely moves (0.884 → 0.884) because it was already well measured, while
`global_slope` gains sharply — **19 % less residual variance** and `prop_kin`
0.376 → 0.467. The gain is larger in EUR than pooled. The 7.0 tables bought
precision on exactly the phenotype the hypothesis is about.

**EUR scan SUBMITTED**: **35559981** (assoc, array 1–110) → **35559982**
(`v2_collect`). It was submitted once its null models were verified rather than
alongside them, so a wiring error could not have cost 110 tasks twice.

**A partial `gwas_summary.tsv` appears mid-run and must not be read.** The
inline per-phenotype collect fires when a phenotype's chr22 task lands, so with
one chromosome done it wrote `n_chr = 1` rows (λ_GC 1.0419 / 1.0133 — close to
the legacy full-genome 1.0363 / 1.0139, but on 121,046 variants, not 9.4 M).
`run_all.sh` chains the final collect exactly for this reason. **Check
`n_chr == 22` before reading any row** (rule 10).

**POOLED SCAN COMPLETE AND VALID** (35559125 all 110 tasks COMPLETED, zero
failures; collected by 35559126). **`n_chr == 22` on all five rows**, so this
one is readable:

| phenotype | λ_GC 6.0 → 7.0 | hits 6.0 → 7.0 | n_mean 6.0 → 7.0 |
|:--|:--|:--|:--|
| baseline_thickness | 1.0363 → 1.0388 | 0 → 0 | 7,893 → 8,395 |
| global_slope | 1.0139 → **1.0265** | 0 → **0** | 7,893 → 8,395 |
| slope_PC1 | 1.0142 → 1.0196 | 0 → 0 | 7,893 → 8,395 |
| slope_PC2 | 1.0358 → 1.0470 | 1 → 1 | 7,893 → 8,395 |
| slope_PC3 | 1.0048 → 1.0149 | 2 → 0 | 7,893 → 8,395 |

9,426,204 variants against v2's 9,413,281, and n_mean up 6.4 %. **Every λ_GC is
inside the 0.99–1.05 band** §6 Step 3–4 requires, so rule 8 is satisfied and
the scan is interpretable. λ_GC rose a little across the board (global_slope
1.0139 → 1.0265) — expected, since λ_GC grows with power at fixed
stratification, and it stays in band.

**No genome-wide hit on `global_slope`, again.** That is the §0/§5 prediction,
not a disappointment: the phenotype's h² z is ~1.2, so a 6 % bigger sample was
never going to produce one. `slope_PC3`'s two legacy "hits" did **not**
reproduce (2 → 0), which is what unreplicated borderline signals do and is
worth stating plainly. `slope_PC2` keeps a single hit (p 4.6e-9 vs the legacy
2.9e-8); it is a secondary target, and one hit at λ_GC 1.047 is not a finding
to lead with.

**EUR SCAN COMPLETE AND VALID** (35559981 + 35559982). `n_chr == 22` on all
five, 7,470,569 variants, n_mean 4,017 → **4,205**:

| phenotype | λ_GC 6.0 → 7.0 | hits 6.0 → 7.0 |
|:--|:--|:--|
| baseline_thickness | 1.0270 → 1.0225 | 0 → 0 |
| global_slope | 1.0001 → **0.9951** | 0 → 0 |
| slope_PC1 | 1.0051 → 0.9996 | 0 → 0 |
| slope_PC2 | 1.0186 → 1.0232 | 0 → 0 |
| slope_PC3 | 0.9953 → 0.9970 | 0 → 0 |

The EUR arm is **better controlled than pooled** (λ_GC 0.995–1.023 against
1.015–1.047), which is what a single-ancestry target should look like, and no
hits anywhere in either vintage.

**Steps 3–4 are DONE for both arms: 696 jobs, zero failures.** Nothing in the
GENESIS track needs re-running.

The honest summary of the scan: **the extra 3,520 six-year scans bought
precision, not discovery.** λ_GC in band everywhere, n up 6.4 % pooled and
4.7 % EUR, no genome-wide hit on `global_slope` in either arm, and the two
borderline `slope_PC3` hits from 6.0 gone. §0 predicted exactly this — the
phenotype's h² z ≈ 1.2 bounds what a GWAS of this size can do. **The
disorder-side test is the PRS (§4 rule 13), which is Step 6, and that is where
the precision gain should actually show up.**

**STEP 5 COMPLETE — GCTA GREML, and this is the substantive result so far**
(job 35562491, all 5 tasks, `work/results_70tab/reml_imp_pooled/`). n = 6,011
unrelated, up from 5,649, exactly the overlap predicted before submitting.

| phenotype | h² ± SE 6.0 | h² ± SE 7.0 | z 6.0 → 7.0 | p 7.0 |
|:--|:--|:--|:--|:--|
| baseline_thickness | 0.2465 ± 0.0489 | 0.2234 ± 0.0464 | 5.04 → 4.82 | 4.2e-07 |
| global_slope | 0.1373 ± 0.0464 | **0.1656 ± 0.0446** | **2.96 → 3.71** | 5.6e-05 |
| slope_PC2 | 0.1605 ± 0.0464 | 0.2022 ± 0.0447 | 3.46 → 4.53 | 8.6e-07 |
| slope_PC3 | 0.0410 ± 0.0448 | 0.1135 ± 0.0437 | 0.92 → 2.60 | 3.4e-03 |
| slope_PC1 | 0.0753 ± 0.0456 | 0.1110 ± 0.0438 | 1.65 → 2.53 | 4.3e-03 |

**Every slope phenotype's h² rose, and the SEs barely moved.** That is the
signature of *less measurement error*, not of a bigger sample: a 6 % larger n
cannot move an SE from 0.0464 to 0.0446, but de-attenuating a noisy phenotype
does raise the point estimate. `global_slope` goes 0.137 → 0.166 with p
improving twenty-fold (1.2e-3 → 5.6e-5). `slope_PC3` was indistinguishable from
zero at 6.0 (z 0.92) and is now nominally non-zero (z 2.60) — it should no
longer be described as unheritable.

**`baseline_thickness` fell slightly (0.2465 → 0.2234) and that is expected,
not a regression.** §1 says the age centre moved 12.44 → 12.80, so the
intercept phenotype is thickness at an older age and should "reproduce in kind,
not to the third decimal". It does: z 4.82, p 4.2e-7, still the strongest of
the five. Do not read the drop as a problem.

**What it licenses, and what it does not.** This raises the ceiling on every
downstream test — h² is what bounds them. But rule 13's threshold is about the
**LDSC** h² z (1.16 at 6.0), not this GREML z, and LDSC on ~4.2 k EUR will stay
well below GREML on 6 k unrelated. Step 7 must report its own z before any rg
is quoted; do not carry 3.71 across.

**STEP 6 (PRS association layer) COMPLETE — this is the headline.** Jobs
35562596 (one-cell smoke test) and 35562748 (40 cells), 160 output files, zero
failures, collected by `setup/collect_final.py` into
`work/results_70tab/prs_final/table_{main,all,family}.tsv`.

Scores were **reused, not recomputed** (§3: "a score never touches our
phenotype"). All 40 cells and their `_zanc` copies were verified present under
`legacy/hpc_v2/work/results_v2/prs_final/` first. `step6_prs_assoc.sbatch` runs
only `prs_assoc.R` and `R/06_prs_family.R` over those directories with the new
$PHENO — the `--prs-dir` is the legacy tree, every `--out` is `results_70tab/`.
Running `setup/prs_final.sbatch` instead would have re-clumped and re-run
PRS-CS and GCTB for days to rebuild files we already had.

**Every single SE got tighter — all 23 matched SCZ/MDD cells, without
exception.** §5's prediction was "same signs, smaller SEs". No sign flipped
(rule 5's bug-until-proven-otherwise test passes) and no SE widened.

`global_slope`, ancestry-matched cells only, β (SD/SD), threshold-adjusted p:

| trait arm | method | β 6.0 → 7.0 | SE 6.0 → 7.0 | p_adj 7.0 |
|:--|:--|:--|:--|--:|
| SCZ_pooled (zanc) | C+T | −0.0328 → −0.0337 | 0.0108 → 0.0104 | 0.0093 * |
| | PRS-CS | −0.0149 → **−0.0242** | 0.0109 → 0.0104 | 0.020 * |
| | SBayesR | −0.0287 → −0.0303 | 0.0109 → 0.0105 | 0.0039 * |
| | SBayesRC | −0.0244 → −0.0289 | 0.0108 → 0.0104 | 0.0054 * |
| SCZ_eur | C+T | −0.0467 → −0.0382 | 0.0152 → 0.0147 | 0.073 |
| | PRS-CS | −0.0265 → −0.0261 | 0.0158 → 0.0146 | 0.074 |
| | SBayesR | −0.0347 → −0.0266 | 0.0153 → 0.0147 | 0.071 |
| | SBayesRC | −0.0353 → −0.0270 | 0.0152 → 0.0146 | 0.065 |
| MDD_pooled (zanc) | C+T | −0.0204 → −0.0216 | 0.0109 → 0.0105 | 0.31 |
| | PRS-CS | −0.0188 → **−0.0206** | 0.0109 → 0.0104 | 0.048 * |
| | SBayesR | −0.0284 → −0.0290 | 0.0110 → 0.0105 | 0.0060 * |
| | SBayesRC | −0.0215 → **−0.0223** | 0.0109 → 0.0104 | 0.033 * |

**Three things changed against §5, and one of them is the answer to the
question this re-run was built to ask.**

1. **MDD crossed into significance under three methods in the pooled arm**
   (PRS-CS, SBayesR, SBayesRC), where 6.0 had only SBayesR. §5 said in advance:
   "MDD crossing into significance under more than one method would be new."
   It did. MDD is no longer only "suggestive".
2. **SCZ pooled is now significant in 4 of 4 methods**, up from 3 of 4 —
   PRS-CS was the null one at 6.0 (−0.0149) and nearly doubled (−0.0242).
3. **SCZ EUR dropped from 3 of 4 significant to 0 of 4** — and this must not be
   buried. The betas moved *toward zero* (C+T −0.0467 → −0.0382, SBayesR
   −0.0347 → −0.0266) while the SEs tightened, so all four now sit at
   p_adj 0.065–0.074: just above threshold, not reversed, same sign. With EUR n
   rising only 4,116 → 4,308 and the subject-level phenotype only r = 0.92
   correlated with its 6.0 self, a drift of this size is within what changing
   the phenotype can do. **Do not report "SCZ replicates in both arms".** The
   honest statement is that SCZ is robust in the pooled arm across all four
   methods and borderline in EUR.

**The control panel behaves exactly as rule 7 documents**, which is what makes
the above readable at all: ASD null in all four methods (+0.005 to +0.034),
ALZ still associated at ~80 % of SCZ's magnitude (C+T −0.0401*, PRS-CS
−0.0363*), and EA running the **opposite** direction (+0.030 to +0.034, three
of four nominally significant). A panel where the controls misbehaved would
invalidate the SCZ/MDD reading; this one does not.

**Still outstanding in Step 6:** the min-p permutation
(`setup/03_prs_paired.sbatch`, §5 benchmark p_perm 0.007 for SCZ EUR) and the
Fulker within-family decomposition tables are written
(`table_family.tsv`, 2,200 rows) but not yet read against §5's 686-pair,
~18 %-power caveat.

### 2026-09-15 — Step 6 continued: within-family, min-p, and a challenge to the SCZ-EUR drop

**Within-family (Fulker) decomposition read against §5.** `table_family.tsv`
has no `matched` column, so select with rule 4 by hand (SCZ/MDD `_pooled` →
`full`, `_eur` → `EUR`). Pairs rose 686 → **726** (EUR) and 1,339 → **1,449**
(pooled). `global_slope`, raw score, C+T at its best threshold:

| arm | method | β_W 6.0 → 7.0 (SE) | β_B 6.0 → 7.0 (p_B 7.0) | p_diff 7.0 |
|:--|:--|:--|:--|--:|
| SCZ_eur | C+T (thr 1) | −0.084 → −0.065 (0.051) | −0.044 → −0.036 (0.019) | 0.58 |
| SCZ_eur | SBayesR | −0.015 → −0.005 (0.048) | −0.037 → −0.029 (0.062) | 0.64 |
| SCZ_pooled | C+T (0p5) | +0.002 → +0.022 (0.047) | −0.045 → −0.048 (0.0010) | 0.15 |
| SCZ_pooled | SBayesR | +0.018 → +0.033 (0.067) | −0.080 → −0.080 (0.0006) | 0.11 |
| MDD_pooled | SBayesR | −0.156 → −0.095 (0.077) | −0.073 → −0.073 (0.0034) | 0.78 |

The verdict is unchanged from §5 and for the same reason: **the within-family
SEs are ~0.05 (EUR) to ~0.07 (pooled), three to four times the between-family
SEs, so no β_W is distinguishable from its β_B** (every p_diff ≥ 0.11) and none
is distinguishable from zero. In the pooled arm the SCZ within-family estimate
sits at or slightly above zero while β_B is strongly negative — the same
pattern as 6.0 (+0.002 / +0.018), which is what stratification *or* an
underpowered within test both look like; §5's "~18 % power" caveat stands with
726 pairs as it did with 686. Report β_W alongside β_B with that caveat, not as
a replication or a refutation.

**Min-p permutation is C+T-only, and the §5 benchmark was not
ancestry-matched.** The permutation (rule 6) corrects for choosing the best of
8 nested clumping thresholds; PRS-CS / SBayesR / SBayesRC each yield one score
and need no such correction. Tracing the §5 row (p_perm 0.007, SCZ EUR) back
through `legacy/hpc_v3/prs_tables/` → `prs_ct_v3` → v1's `prs_imp/SCZ_*.weights`
→ `sumstats/SCZ.tsv` shows it was built from the **PGC3 *primary*
(multi-ancestry) file and read on the EUR target** — the rule-4 mismatch that
loses power but is not confounded. `step6_minp_permutation.py` (job 35578468,
2,000 family-block permutations, seed as the module) therefore runs two
designs: *matched* (`prs_final/CT/<arm>` in the rule-4 stratum) and
*legacy-like* (`prs_ct_v3`, both strata) so the benchmark can be compared
like-for-like. Output: `prs_final/table_minp_permutation.tsv`.

**The SCZ-EUR drop (3/4 → 0/4 significant) is being decomposed, not
accepted.** Checked what changed between the two runs' inputs on the 8,081
shared children: family partition identical (6,788 families both vintages);
`site` differs for **one** child; the 7.0 PCs correlate 0.9988 with 6.0's;
`n_visits` changed for 2,432 (the children who gained a scan); and
`baseline_age` changed for *all* 8,081 — which looked alarming until it turned
out to be `age_c` = age − model age-centre, and the centre moved 12.4406 →
12.7973: a constant shift of 0.357 (sd 0.006 among children with unchanged
visits) that `prs_assoc.R` re-centres away. So the only material differences
are the **phenotype values** (r = 0.92) and the **480 added children** (192
EUR). Job 35578218 runs `prs_assoc.R` on the four SCZ_eur cells across a
ladder that changes one thing at a time:

| variant | phenotype | covariates + family ids | subjects |
|:--|:--|:--|--:|
| A | 6.0 | 6.0 | 8,082 |
| C | **7.0** | 6.0 | 8,081 shared |
| B | 7.0 | **7.0** | 8,081 shared |
| D (reported) | 7.0 | 7.0 | **8,596** |

A must reproduce the 6.0 table (pipeline check; the only expected deviation is
rule 14's `age_c`, ≤ 0.006 in β). A → C isolates the phenotype change, C → B
the covariate/PC change, B → D the sample growth. Note the user's suggestion —
7.0 phenotype on the 6.0 subjects — is variant B, and it should **not** come
out identical to 6.0, because the phenotype values themselves moved for the
2,478 children who gained a six-year scan; identity is expected only for A.

**Rule 10, met in the wild.** The first ladder submission (35578218) came back
COMPLETED in 20 s with no outputs and no log: its inputs, outputs and log path
were all under the session scratchpad in `/tmp`, which is **node-local** — the
compute node saw an empty directory, every `prs_assoc.R` call failed behind the
loop's `|| true`, and the log went to that node's own `/tmp`. Anything a job
reads or writes must be on rds. Rebuilt under `genetic_analysis/work/scz_eur_check/`
and resubmitted as 35578568.

**The ladder (35578568) settles the SCZ-EUR question: the drop is real, it is
the phenotype, and the pipeline is exact.** `global_slope`, SCZ_eur → EUR,
C+T at threshold 1 (the 6.0 best) and the three single-score methods:

| variant | n | C+T β (SE) | PRS-CS β | SBayesR β | SBayesRC β |
|:--|--:|:--|:--|:--|:--|
| 6.0 table | 4,116 | −0.0467 (0.0152) | −0.0265 | −0.0347 | −0.0353 |
| A  6.0 pheno, 6.0 cov | 4,116 | **−0.0467 (0.0152)** | **−0.0265** | **−0.0347** | **−0.0353** |
| C  7.0 pheno, 6.0 cov, shared | 4,115 | −0.0406 (0.0150) | −0.0221 | −0.0259 | −0.0252 |
| B  7.0 pheno, 7.0 cov, shared | 4,115 | −0.0415 (0.0150) | −0.0241 | −0.0250 | −0.0249 |
| D  7.0 full (reported) | 4,308 | −0.0382 (0.0147) | −0.0261 | −0.0266 | −0.0270 |

1. **A reproduces the 6.0 table to every printed digit**, all four methods:
   the new scripts, the reused scores and the EUR keep list are doing exactly
   what v2 did. (It also shows the 6.0 table already carried rule 14's `age_c`
   fix — there is not even a ≤ 0.006 deviation.) There is no bug to find.
2. **A → C is the whole drop.** Same 4,115 children, same covariates, same
   family ids; only the phenotype values changed — and β moves −0.0467 →
   −0.0406 (C+T), −0.0347 → −0.0259 (SBayesR), −0.0353 → −0.0252 (SBayesRC).
   Every subsequent rung is noise around that.
3. **C → B (7.0 PCs, 7.0 `age_c`, 7.0 family ids) changes nothing** (≤ 0.002),
   as the input diagnostics predicted.
4. **B → D (adding 192 EUR children) does not rescue it**: C+T slips a little
   further (−0.0415 → −0.0382), PRS-CS firms slightly (−0.0241 → −0.0261), the
   SBayes pair are flat; the SEs tighten as they should.

Is a shift of 0.006–0.010 in β surprising for a phenotype that correlates
r = 0.92 with its former self? Roughly not: for the same predictor, two
outcomes with correlation ρ differ in standardised β by about
SE·√(2(1−ρ)) ≈ 0.015·√0.16 ≈ 0.006 as a one-SD movement, so these are 1–1.7 SD
shifts in the direction of attenuation, in the arm that at 6.0 carried the
headline p-value — i.e. the arm most exposed to winner's curse. The pooled arm,
which has 2× the n and was less selected, held or strengthened under the same
phenotype change. **Answer to the challenge: yes, it is real; it is not an
artefact; it is what the extra 3,520 six-year scans did to the EUR slope
estimates.** The honest statement stands: SCZ robust pooled (4/4), borderline
EUR (p_adj 0.065–0.074), same sign everywhere.

**Min-p permutation COMPLETE** (35578468; 2,000 family-block permutations,
`prs_final/table_minp_permutation.tsv`, 100 rows). `global_slope`, C+T only:

| design | arm → stratum | n | min p | **p_perm** | p_Bonferroni | 6.0 p_perm |
|:--|:--|--:|--:|--:|--:|--:|
| matched | SCZ_pooled → full | 8,596 | 0.0011 | **0.0065** | 0.0086 | — |
| matched | SCZ_eur → EUR | 4,308 | 0.0155 | **0.053** | 0.124 | — |
| matched | MDD_pooled → full | 8,596 | 0.0276 | 0.099 | 0.221 | — |
| matched | ALZ → EUR | 4,308 | 0.0050 | 0.026 | 0.040 | — |
| matched | ASD → EUR | 4,308 | 0.0169 | 0.077 | 0.135 | — |
| legacy-like | SCZ (primary) → EUR | 4,308 | 0.0177 | 0.061 | 0.141 | **0.007** |
| legacy-like | SCZ (primary) → full | 8,596 | 0.0017 | 0.011 | 0.014 | 0.013 |
| legacy-like | MDD (div) → full | 8,596 | 0.0134 | 0.047 | 0.108 | 0.068 |
| legacy-like | ASD → full | 8,596 | 0.0003 | 0.003 | 0.002 | 0.013 |

Four readings:

1. **SCZ pooled survives the correction that matters** (p_perm 0.0065), so the
   C+T headline is not a best-of-eight artefact. Permutation and Bonferroni
   agree here because the signal is well clear of the threshold.
2. **SCZ EUR is borderline under permutation too** (0.053 matched, 0.061
   legacy-like), exactly where the mixed-model p_adj put it (0.065–0.074). The
   §5 benchmark of 0.007 was the *legacy-like* design, and the same design on
   the 7.0 phenotype gives 0.061 — the ladder's attenuation, seen a third way.
   Every route to the EUR arm now says "borderline"; none says "gone".
3. **MDD's new significance does not come from C+T.** MDD pooled C+T fails
   permutation (0.099), as it failed p_adj in the association table (0.31).
   The 3-of-4 result rests on PRS-CS, SBayesR and SBayesRC, each a single score
   with **no threshold multiplicity to correct**. State it that way.
4. **Rule 4's false ASD "hit" reappears on cue.** ASD is a European-only GWAS;
   read on the *pooled* target (the confounded direction) it gives p_perm
   0.003 — "stronger" than the 6.0 artefact (0.013). In the matched EUR
   stratum it is 0.077, null. This row exists in the table only because the
   legacy-like design was run for comparison; it is a demonstration of the
   confound, not a result.

**Step 6 is complete** in every part §6 lists: association layer (40 cells,
4 methods × 10 trait-arms, both strata, raw and within-ancestry-standardised),
Fulker within-family decomposition (2,200 rows), min-p permutation (two
designs), and the five controls in every table.

### 2026-09-15 — Step 7 submitted; HCP-MMP arm opened in parallel

**Step 7 (LDSC + MAGMA + prioritised panel) — job 35581017.**
`step7_ldsc_magma.sbatch` is `setup/eur_magma_ldsc.sbatch` + `eur_prio_gsa.sbatch`
in one job: stage the five EUR `.sumstats.tsv.gz` into the `.fastGWA` layout,
then v1's own `04_magma.sbatch` / `05_ldsc_rg.sbatch` / `work/prioritised_gsa.sbatch`
verbatim with `GWAS_DIR`/`MAGMA_DIR`/`LDSC_DIR`/`MANIFEST` redirected into
`results_70tab/{magma_eur,ldsc_eur,magma_prio_eur}` (rule 9). Three migration
facts it depends on, each verified before submitting: (i) `legacy/hpc/config.sh`
derives `REPO_ROOT` from its own location, which is now `legacy/`, so with
`REPO_ROOT` **unset** v1's config resolves every asset under `legacy/hpc/work`
— all 12 checked present; genetic_analysis/config.sh does not export it, and the
job unsets it anyway; (ii) `legacy/hpc/work/bin/{ldsc.py,munge_sumstats.py,Rscript}`
were bash wrappers with the pre-move path compiled in — fixed in place,
originals `*.pre70bak` (six other v1 sbatch files still carry the old absolute
path in comments/log paths; none is on this run's path); (iii)
`prioritised_gsa.sbatch` defaults `MAGMA` to `$REPO/hpc/work/bin/magma`, so it
is passed explicitly. **Read `h2_z` first** (rule 13): the 6.0 EUR values were
1.16 (`global_slope`) and 4.49 (`baseline_thickness`).

**HCP-MMP arm — the backfill merged, and the pipeline made parcellation-agnostic.**

*Merge.* `origin/hcp-backfill` (252c05d, three commits on top of f9a8806)
fast-forwarded onto `main`: `src/abcd/{dk_stats,hcp_stats,io}.py`, the backfill
tools, validators and `docs/hcp_census/`. No overlap with this run's edits. The
backfilled tables were already in place in this checkout's
`abcd-data-release-7.0/processed/hcp/` (**33,795 sessions**, written
2026-09-14 16:27, 365 columns; 33,792 of the 33,794 DK sessions present).
`abcd-7.0/processed/hcp/` (24,921 sessions, 2026-09-12) is now a stale
duplicate; harmless because `_RELEASE_DIR_PATTERNS` prefers the release
directory, but it should go.

*A data trap, caught before fitting.* The HCP adapter keeps all 360 parcels and
the hippocampal parcel **`H` has thickness exactly 0 in 7,792 (lh) and 1,945
(rh) sessions** (`mris_anatomical_stats` writes 0 for a parcel with no
vertices; only 5 zero cells exist outside H). A zero is not a missing value,
so `complete_regions` QC would pass it and it would sit inside the 360-parcel
whole-cortex mean and the slope PCs. There was no way to exclude a region. Two
code changes, both merged-code-adjacent and covered by the existing tests:
`RunConfig.exclude_regions` (applied in `assemble` after load, before QC) and,
in `Release70Adapter._imaging_hcp`, exact zeros → NaN so QC sees the stray
cells. `configs/ct_70_hcp_noglobal_mv2.yaml` gets `exclude_regions: [H]` →
**358 parcels**. The field participates in the hash — it changes the phenotype
— **but an empty value is hash-neutral**, because adding a field otherwise
re-hashed every existing run: the settled DK id `thickness_dsk_70_139406217085`
was verified to still resolve. HCP run id: `thickness_hcp_70_aa6e91efba82`.
`tests/test_config.py` + `test_assemble.py`: 2 failed / 7 errors **before and
after**, all `DataRootError: Release 5.1 not found` — environmental, not mine.

*Layout: one switch, two roots.* `config.local.sh` now has `PARC` (`dsk`
default, `hcp`), which sets `ABCD_CONFIG_PARC`, `PHENO_DIR` and `OUT_V2`:

| | DK | HCP-MMP |
|:--|:--|:--|
| local config | `ct_70_noglobal_mv2_genetic` | `ct_70_hcp_noglobal_mv2` |
| phenotype export | `work/pheno_70tab/` | `work/pheno_70tab_hcp/` |
| results root | `work/results_70tab/` | `work/results_70tab_hcp/` |
| genotype-only assets (GRM view, kinship, scores) | shared | shared |

Every step reads `$PHENO_DIR` and writes `$OUT_V2`, so
`PARC=hcp bash run_all.sh 03 04`, `PARC=hcp sbatch step5_reml.sbatch`, etc.,
run the identical code on the other atlas. The DK directories keep their
unsuffixed names because they pre-date the switch and are cited throughout
this log; the layout is otherwise symmetric. `.gitignore` uses
`results_70tab*` so both roots follow the same summary-tables-only contract.
`step1_export_pheno.sbatch` and `step1_align_and_preflight.sbatch` take the
config and run directory from `PARC`; `step1_verify_export.py --parc hcp` gates
the HCP export against the **DK run** (same children, scans, QC and model —
only the atlas differs — so n within 1 %, 358 regions, `exclude_regions == [H]`,
and `global_slope` / `baseline_thickness` r > 0.95 with DK on shared children)
rather than against `docs/vintage_comparison.csv`, which is DK-specific.

**HCP jobs: 35581499/35581500 were CANCELLED by me and resubmitted as
35581806 (export) → 35581807 (`afterok`: verify → align → GRM view →
preflight).** The first export's log showed `R version 4.4.3` — the GENESIS R,
not the abcdR 4.5.3 lme4 environment the DK export was fitted with. Cause:
making the script parcellation-aware meant sourcing `config.sh` first, which
sets `RSCRIPT` (to GENESIS) from `config.local.sh`, so the script's own
`${RSCRIPT:-abcdR}` fallback silently kept the wrong interpreter. It would have
crashed anyway — the GENESIS R lacks `arrow` and `yaml`, both required by
`R/fit_lmm.R` — but a silent success on a different lme4 build would have been
worse: the DK-vs-HCP comparison would then confound atlas with software. The
lme4 interpreter is now pinned (`LMER_RSCRIPT`, default abcdR) independently of
`RSCRIPT`. Same lesson as the `PY` incident: an interpreter that is right for
one script is wrong for another, and a `:-` fallback cannot tell.

**HCP export COMPLETED (35581806) and gated.** Sample construction is the
DK run's to the child: 8,716 subjects, 7,285 families, 10,183 after QC, age
centre 12.797; 26,946 scans against DK's 26,949 (the three lost are the stray
zero cells the adapter now treats as missing). `excluded regions ['H']: 2
labels dropped` → **358 regions**. Export: 5 phenotypes, 10 PCs, lme4 R 4.5.3
— the same build as DK.

The gate's first run **failed on one row**: `global_slope` r vs DK = 0.9484
against a 0.95 line I had set by judgement; `baseline_thickness` was 0.9777.
The chained job 35581807 stopped there, as designed. Diagnosed before touching
the threshold:

| | DK ↔ HCP |
|:--|--:|
| `baseline_thickness` r (Pearson / Spearman / 1 % trimmed) | 0.978 / 0.975 / 0.979 |
| `global_slope` r | 0.948 / 0.940 / 0.952 |
| `global_slope` r by visits: 2 / 3 / 4 | **0.942 / 0.948 / 0.953** |
| slope reliability, median, 2/3/4 visits: DK | 0.119 / 0.210 / 0.242 |
| same, HCP | 0.150 / 0.213 / 0.246 |
| `global_slope` native mean (mm/yr): DK / HCP | −0.01870 / −0.01993 |
| `slope_PC1` / `PC2` / `PC3` \|r\| | 0.877 / 0.788 / 0.509 |

Agreement **rises with visits**, i.e. with slope reliability — the fingerprint
of measurement noise in a noisy quantity, not of a systematic difference —
and trimming 1 % barely moves it. The two whole-cortex means are different
weightings of the same cortex (68 large DK regions vs 358 small parcels, both
unweighted), so identity was never the expectation; the native means differ
accordingly. The slope **PCs are atlas-specific decompositions and are not
the same phenotype across atlases** (PC1 sign flips; PC3 shares only a
quarter of its variance) — a point the comparison must respect: compare
`global_slope` and `baseline_thickness` across atlases; treat the PCs as
within-atlas readouts. Gate now: `global_slope` r > 0.90 (below that, a bug —
wrong sessions, wrong QC, a zero parcel leaking in — becomes the likelier
reading; for scale, the 6.0 → 7.0 change of the *same* DK phenotype was
0.92), `baseline_thickness` > 0.95, PCs reported but not gated. **15/15 pass.**
Chain resubmitted as **35582792** (verify → align → GRM view → preflight).

**HCP steps 1–2 DONE (35582792): gate 15/15, align, GRM view, `preflight
PASSED`.** The HCP analysis set is the DK one exactly — 8,716 → **8,596**
phenotyped-and-genotyped, 7,198 families (1,350 multi-member), **EUR 4,308**,
ID token join 8,596 of 8,596 — so every DK-vs-HCP difference downstream is the
atlas, not the sample. Native units (`pheno_70tab_hcp/align_report.tsv`):
`global_slope` mean −0.019930 mm/yr, SD 0.0013090 (DK: −0.018696, 0.0013222);
`baseline_thickness` 2.7440 ± 0.0657 (DK 2.7121 ± 0.0630) — the HCP whole-cortex
mean weights 358 small parcels equally, so it is a slightly different average of
the same cortex.

**HCP steps 3–6 SUBMITTED**, the identical scripts under `PARC=hcp`, all into
`work/results_70tab_hcp/`:

| step | jobs |
|:--|:--|
| 3–4 pooled: null → assoc (110) → collect | 35583560 → 35583561 → 35583562 |
| 3–4 EUR (`KIN_DIR=kinship_eur`, `_eur` dirs) | 35583571 → 35583572 → 35583573 |
| 5 GCTA GREML | 35583574 (array 1–5) |
| 6 PRS association layer, 40 cells | 35583575 |
| 6 min-p permutation (`--pheno-dir pheno_70tab_hcp`) | 35583576 |

Both scan arms were submitted together rather than pooled-then-EUR as for DK,
because the arm wiring (kinship = arm, pooled GDS) was proven there.
`compare_parcellations.py` (new) builds
`results_70tab_hcp/compare/table_dk_vs_hcp.tsv` from whatever summaries exist in
both roots — one row per (readout, phenotype, arm, method) with DK and HCP
columns — and flags the slope PCs `comparable=no`; smoke-tested on the DK-only
state (165 rows). HCP step 7 follows once `assoc_eur` collects.

**HCP first readouts (null models 10/10, GREML 4/5) — the whole-cortex
phenotypes are atlas-independent, the PCs are not.**

GREML on the same 6,011 unrelated children:

| phenotype | DK h² (z) | HCP h² (z) |
|:--|:--|:--|
| baseline_thickness | 0.2234 ± 0.0464 (4.81) | 0.2247 ± 0.0465 (4.84) |
| global_slope | 0.1656 ± 0.0446 (3.71) | 0.1560 ± 0.0446 (3.50) |
| slope_PC2 | 0.2022 (4.52) | 0.0882 (2.03) |
| slope_PC3 | 0.1135 (2.60) | 0.2290 (5.07) |

The two whole-cortex estimates agree to within a quarter of an SE; SNP h² of
cortical thickness and of its thinning rate does not depend on how the cortex
is carved. The PCs reorder between atlases (DK's heritable component is PC2,
HCP's is PC3), exactly as their cross-atlas |r| of 0.79 / 0.51 said they would
— they are within-atlas readouts and their rows are marked `comparable=no`.

Null models (GENESIS sparse-kinship decomposition, pooled / EUR `prop_kin`):
baseline 0.885 / 0.884 (DK) vs 0.900 / 0.904 (HCP); **`global_slope` 0.440 /
0.467 (DK) vs 0.561 / 0.592 (HCP)**, with `varcomp_resid` 0.478 → 0.381.
The HCP whole-cortex slope carries less residual noise — consistent with its
higher 2-visit reliability (0.150 vs 0.119: an unweighted mean over 358
parcels averages out more per-parcel measurement error than one over 68). But
note what did **not** move: GREML h² (0.166 → 0.156). The extra variance the
sparse kinship absorbs is the shared-family part, not additive SNP variance
— the same rule-11 distinction as before, now visible as a DK–HCP contrast.
If the PRS betas come out tighter in HCP, this is the reason to expect it.

**DK STEP 7 COMPLETE** (35581017, 1 h 32 min; `ldsc_eur/`, `magma_eur/`,
`magma_prio_eur/` under `results_70tab/`; 10 / 60 / 133+112 rows).

*LDSC, EUR arm — h² z first (rule 13):*

| phenotype | h²_obs ± SE 6.0 → 7.0 | **h² z 6.0 → 7.0** | rg(SCZ) 6.0 → 7.0 |
|:--|:--|:--|:--|
| baseline_thickness | 0.467 ± 0.104 → 0.421 ± 0.097 | 4.49 → **4.34** | 0.026 ± 0.052 → 0.024 ± 0.053 |
| global_slope | 0.124 ± 0.107 → 0.055 ± 0.097 | 1.16 → **0.57** | −0.149 ± 0.109 → −0.164 ± 0.200 |
| slope_PC2 | 0.377 → 0.271 | 3.46 → 2.60 | −0.052 → −0.109 |

`global_slope`'s LDSC h² z **fell** (1.16 → 0.57) while its GREML h² **rose**
(0.137 → 0.166, z 3.71). These are not in tension: LDSC on ~4.2 k EUR
children has an SE of ~0.10 on h², so both vintages are within one SE of each
other and of zero; GREML on 6 k unrelated children with a dense GRM is the
canonical estimate (rule 11). The consequence is rule 13's, stated plainly:
**rg between `global_slope` and SCZ/MDD is uninformative, not null** — the
rg(SCZ) SE nearly doubled (0.109 → 0.200) because it scales with 1/√h²_z. Do
not quote −0.16 as a genetic correlation. `baseline_thickness` remains the
positive control (z 4.34) and its rg with SCZ and MDD is a stable, tight null
(0.024 ± 0.053; −0.011 ± 0.055).

*MAGMA prioritised gene sets (EUR, marginal):*

| set → phenotype | β (p) 6.0 | β (p) 7.0 |
|:--|:--|:--|
| SCZ_locus_pool (455) → baseline_thickness | 0.187 (1.1e-4) | 0.172 (4.2e-4) |
| SCZ_locus_pool → global_slope | 0.124 (7.4e-3) | **0.143 (3.0e-3)** |
| SCZ_prioritised (101) → global_slope | −0.008 (0.54) | 0.019 (0.42) |
| MDD_pool (1,881) → global_slope | 0.029 (0.13) | 0.009 (0.37) |

The §5 positive control reproduces in kind (`SCZ_locus_pool` → baseline), and
its `global_slope` counterpart **strengthened** (p 0.0074 → 0.0030) — the one
place in step 7 where the extra precision shows. The 101 fine-mapped
prioritised genes and the MDD pools stay null, as at 6.0.

*MAGMA gene-property, AHBA C1–C3 (±) and snRNA PC1 on the matched EUR
scans:* **null for `global_slope` in both directions**, as §5 records
(all p > 0.19). Across the whole 60-row table only four rows sit below 0.05,
none survives any correction, and two are the reverse-direction
disorder-GWAS → `baseline_thickness` rows already present at 6.0
(SCZ 0.012 → 0.0075; MDD 2.1e-4 → 5.1e-4) — i.e. genes whose association
with baseline thickness is higher also carry more SCZ/MDD signal, the
thickness-not-thinning pattern seen throughout.

`build_summary_70tab.py` (new) now writes §7 item 1:
`results_70tab/summary_70tab.tsv`, **659 rows, 559 with a 6.0 partner** (the
100 without are the matched-design min-p rows, the ALZ_IGAP / EA controls
and the within-family rows v2 never tabulated). Fixed on first run: LDSC h²
rows were emitted once per disorder file and cross-multiplied on the legacy
merge; the disorder is now in the readout name.

**DK STEP 8 COMPLETE (ahba_pls H3) — null, as the prior said.**
`step8_ahba_pls_h3.sh` (new; seconds, login node) regresses each phenotype's
EUR gene-level Z from step 7's `.genes.raw` on the NSPN-PLS2 / thinning
signature weights (`thinning_Z_ds0/25/50`) and on AHBA C3 — 4 marginal models
+ the 3 DS columns each conditioned on `AHBA_C3` — for the five phenotypes: 35
tests, `results_70tab/magma_ahba_pls_h3/table_h3.tsv` (50 rows incl. the
conditioning variable's own coefficient).

Two MAGMA traps on the way, both fixed in the script: (i) `--gene-covar` on the
4-column file **silently dropped every column containing NA** ("variable
removed during preprocessing") — ds25/ds50/C3 have 3,447 / 6,846 / 6,862 NA of
13,789 genes — so the first run tested ds0 only and the conditional model had
nothing to condition on; the file is now split per variable with NA genes
removed, which is also literally what H3 specifies. (ii) The conditional
`.gsa.out` carries an extra `MODEL` column; the collector is header-driven.

| model | `global_slope`: ds0 / ds25 / ds50 β_std (p) | AHBA_C3 β_std (p) |
|:--|:--|:--|
| marginal | −0.002 (0.84) / 0.008 (0.35) / 0.013 (0.23) | −0.004 (0.70) |
| conditioned on C3 | 0.029 (0.072) / 0.031 (0.061) / 0.033 (0.067) | −0.026 (0.11) |

Every other phenotype: all p > 0.21. **No cell of the 35 is below 0.05.**
The only thing near threshold is the conditional `global_slope` row —
the same signal three times (the DS columns are near-duplicates of one
weighting), positive, p 0.06–0.07, uncorrected, on a GWAS whose LDSC h² z is
0.57 — which is the definition of not a result. Per FOLLOWUP_GENETICS.md:
reported as null; not hunted.

**The DK arm is now complete through every step of §3 (1–8).** Its §7 table
is `results_70tab/summary_70tab.tsv`.

**Next:** read the HCP scans (λ_GC, `n_chr`); `PARC=hcp` step 7 then step 8;
then the DK-vs-HCP table and the HCP `summary_70tab.tsv`.
