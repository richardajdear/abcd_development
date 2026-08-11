# README_HPC — state of the cluster run

Living document for the CSD3 execution of `hpc/`. `README.md` states what the
pipeline *is*; this file states where it *is up to*, what was found broken, and
what a person picking it up next should type.

**Last updated:** 2026-08-11
**Status:** **complete end to end.** Steps 00–06 have all run genome-wide on the
real genotypes and produced results — §12. Seven defects were found and fixed
along the way (§2, §2b, §6c, §6d, §9b, §12b, §13), five of which produced a
*silently wrong or silently empty* result rather than an error.

The headline numbers: SNP-h² for `baseline_thickness` is **0.575 (GCTA) / 0.584
(LDSC)**, no genome-wide-significant loci for any phenotype at N = 4,119, and no
genetic correlation with SCZ or MDD. The sample is EUR-only and about half the
size `README.md`'s power section assumes (§5b, §11.1), which is the single
biggest constraint on everything above.

Two follow-up analyses were added after the main run: MAGMA competitive tests of
every phenotype against the **SCZ prioritised genes** (Trubetskoy 2022, §13) and
the **MDD high-confidence genes** (Cell 2024, §14). The like-for-like comparison
— SCZ's 120 prioritised vs MDD's 308 high-confidence — is null on both sides for
every phenotype. The one signal is that `baseline_thickness` is enriched in the
broad 462-gene SCZ *locus pool* (p = 1.9e-04), carried by the genes
prioritisation did **not** nominate.

**Bottom line across every method tried** — LDSC rg, PRS, and two prioritised
gene sets — the developmental *slope* phenotypes show no relationship to SCZ or
MDD genetics. Everything non-null sits on `baseline_thickness`, and is an
unsigned-overlap test rather than a directional one.

---

## 1. What changed relative to `README.md`'s assumptions

`README.md` says nothing here has touched real genotypes and that the account is
`VERTES-SL2-CPU`, the partition `cclake`, the tree
`/home/rajd2/rds/hpc-work/ABCD/`. Six of its assumptions do not hold on this
cluster.

| assumption in `README.md` | reality | consequence |
|---|---|---|
| account `VERTES-SL2-CPU` | **exhausted** — 109,592 h used of 109,592 h; `sbatch` returns `AssocGrpCPUMinutesLimit` | use `VERTES-SL3-CPU` (155,042 h left). Set in `config.local.sh` |
| partition `cclake` | **effectively decommissioned** — 641 of 665 nodes sit in permanent MAINT reservations named `Retired_Cclake`, `cclake-to-remove` and `Arcus_Hypervisors`; 2 nodes idle | jobs queue forever. Use `icelake` (552 nodes, live). Set in `config.local.sh` |
| genotypes at `$ABCD_HPC_ROOT/genotype/ABCD_release_7.0_QCed` | no such file. The usable set is `ABCD/genetics_copy/abcd_eur` — 5,678 EUR subjects × 13,697,177 SNPs, hg19 | `GENO` repointed |
| `01_grm` must build all three GRMs (~2 h) | the full GRM **and** the 0.05 sparse GRM already exist for exactly these 5,678 subjects | only `--grm-cutoff` is left to run |
| `gcta_export` writes IDs the `.fam` matches | it writes `FID=10210 IID=NDAR_INV005V6D2C`; the `.fam` says `FID=IID=sub-NDARINV005V6D2C` | **both** columns mismatch → GCTA would analyse zero subjects, silently. See §4 |
| the target sample has N ≈ 8,000 | 8,192 subjects are phenotyped but only **5,678 are genotyped**, and those are EUR-only | every power statement in `README.md` is optimistic by roughly the ratio; see §11 |

## 2. Defect found and fixed: `require_paths GRM_SPARSE`

`hpc/config.sh` validated `GRM_SPARSE` by looking for `$GRM_SPARSE.grm.bin` or
`.grm.gz`. A sparse GRM has neither: `gcta --make-bK-sparse` writes `.grm.sp`
(a three-column list of retained pairs) plus `.grm.id`. So

```
ERROR: GRM_SPARSE=.../abcd_eur_sp -- neither .../abcd_eur_sp.grm.bin nor .../abcd_eur_sp.grm.gz found
```

on any real sparse GRM, and `03_gwas.sbatch` aborted in `require_paths` before
running a single command. The local test harness did not catch this because it
never reached step 03 with a GCTA-built sparse GRM.

Fixed in `hpc/config.sh` by giving `GRM_SPARSE` its own branch that checks
`.grm.sp`. This is the one change made to a tracked pipeline file.

## 2b. Defect found and worked around: LDSC vs pandas

`munge_sumstats.py` crashed on every input with

```
KeyError: '[-1 -1 -1 ... -1 -1 -1] not in index'
  at  dat.loc[~jj, [i for i in dat.columns if i != 'SNP']] = float('nan')
```

Cause: `jj` is created as a bool Series and then assigned into with
`jj[ii] = match`, which upcasts it to **object** dtype (NaN outside `ii`).
pandas ≥ 0.23 does not recognise an object array as a boolean mask, falls back
to label lookup, and every label misses. This is the same class of problem
`README.md` anticipates for the PyPI `ldsc` wheel, but it bites the Broad
repo copy on this cluster's py2.7 env (pandas 0.23.4) too.

`work/ldsc_patched/` is a copy of the Broad tree with one line changed —
`_mask = ~np.asarray(jj.fillna(False).astype(bool))` — and `LDSC`/`LDSC_MUNGE`
point at it. After the patch, munging a 200k-row slice of the SCZ file gives
mean χ² = 2.147, λ_GC = 1.838, 69 genome-wide-significant SNPs: sensible for
SCZ, so the columns and the signed-statistic spec are being read correctly.

`ldsc.py --rg` was also exercised and reads the LD scores and writes the summary
block that `05_ldsc_rg.sbatch`'s awk parses. The probe's rg is `NA`, which is an
artefact of the probe rather than a defect: **`SCZ.tsv` is not sorted by
chromosome** (it begins at chr8), so a head-truncated SCZ slice and a
head-truncated MDD slice share zero SNPs. Do not use `head` to make LDSC test
inputs from these files.

Running `05_ldsc_rg.sbatch`'s collector awk over that probe log confirms it
parses the block correctly, and exposes one small logic flaw worth fixing: when
LDSC fails and writes an all-`NA` row, `h2se+0` is 0, so `z` becomes `NA` and

```awk
under = ((z!="NA" && z<4) ? "yes" : "no")
```

labels the row `underpowered=no` — i.e. a *failed* rg is reported as adequately
powered. It should be `unknown` (or `yes`), since the whole point of the column
is to stop a null rg being read as evidence of no overlap.

## 3. Environment built (all under `hpc/work/`, gitignored)

Nothing usable was on `PATH`: no `gcta64`, `magma`, `plink`, no Python ≥ 3.10
(system default is 3.7.4), and system R 4.5.3 has neither `lme4` nor `arrow`.
`/home` is at 48.6 GB of a 52.4 GB quota, so everything went to
`/rds/user/rajd2` instead.

| what | where | note |
|---|---|---|
| `work/bin/gcta64` | → `ABCD/gcta-1.94.1-.../gcta64` | v1.94.1 |
| `work/bin/magma` | → `hpc-work/magma/magma` | v1.10 |
| `work/bin/plink2` | → `hpc-work/plink2` | v2.00a6LM |
| `work/bin/plink` | downloaded | v1.9.0-b.7.7 — **needed**: `06_prs` uses 1.9's `--clump`/`.clumped` and `--score f 1 2 3 sum`; PLINK 2 has `--clump` but a different output format |
| `work/bin/{ldsc,munge_sumstats}.py` | wrappers | exec `work/ldsc_patched/` (Broad LDSC, Python 2) under the `~/.conda/envs/ldsc` py2.7 env — see §2b |
| `work/bin/Rscript` | wrapper | conda R 4.4.3 with `lme4`, `lmerTest`, `arrow`, `dplyr`, `tidyr`, `purrr`, `optparse`, `yaml`, `data.table` |
| `work/envs/abcd` | micromamba | Python 3.11 + pandas/pyarrow/scipy/nibabel/yaml — runs the `abcd` package |
| `work/envs/abcdR` | micromamba | the R above |

`work/mamba/` is the micromamba package cache and can be deleted once the envs
are built.

## 4. The ID-alignment problem, and `work/align_ids.py`

`abcd.gcta_export` documents that "the genotype `.fam` follows the genetics
convention" and normalises IIDs to `NDAR_INV…`, writing the family ID into FID.
The `.fam` in use here follows neither half of that:

```
gcta_export:  FID=10210               IID=NDAR_INV005V6D2C
abcd_eur.fam: FID=sub-NDARINV005V6D2C IID=sub-NDARINV005V6D2C
```

GCTA identifies an individual by the **FID+IID pair**, so this is a
zero-subject analysis that reports no error. Worse, `00_check_inputs.sh` would
not catch it: its ID-intersection check compares column 2 only, so an aligned
IID with a mismatched FID reads as a clean intersection.

`work/align_ids.py` rewrites both columns of the three GCTA input files to the
`.fam`'s own spelling, joining on the 8-character NDAR token — the one part
invariant across all three spellings ABCD uses (`NDAR_INVxxxxxxxx`,
`sub-NDARINVxxxxxxxx`, and 7.0's bare `sub-xxxxxxxx`). It runs as the last step
of `work/fit_pheno.sbatch`.

Two things worth deciding upstream rather than patching here:

- **`00_check_inputs.sh` should compare FID+IID, not IID.** Its comment names
  exactly this failure mode and then checks the wrong column.
- **`_to_genetics_id` is release- and site-specific.** Which convention the
  `.fam` uses is a property of the genotype file, not of ABCD; deriving the
  target form *from the `.fam`* (as `align_ids.py` does) removes the assumption.

## 5. Inputs located on CSD3

All validated by `require_paths`; every one of these exists and is readable.

| config variable | resolved path | contents |
|---|---|---|
| `GENO` | `ABCD/genetics_copy/abcd_eur` | 5,678 subj × 13,697,177 SNPs, hg19, EUR only, `--maf 0.001` |
| `GRM` | `ABCD/genetics/abcd_eur` | full GRM, 5,678 subj, GCTA 1.94.1 |
| `GRM_SPARSE` | `ABCD/genetics/abcd_eur_sp` | 6,934 pairs above 0.05 |
| `MAGMA_REF` | `magma/reference_data/g1000_eur` | 1000G EUR, 22,665,064 SNPs |
| `MAGMA_GENE_LOC` | `magma/gene_locations/NCBI37.3.gene.loc` | |
| `MAGMA_ANNOT` | `magma/snp_gene_annotation/ncbi37.window35-10.genes.annot` | prebuilt against that same `.bim`, so `04_magma`'s build branch is a no-op |
| `LD_REF` / `LD_WEIGHTS` | `ldsc/eur_w_ld_chr/` | 22 chromosomes |
| `HM3_SNPLIST` | `ldsc/eur_w_ld_chr/w_hm3.snplist` | |
| `GENESET_DIR` | `hpc/work/genesets` | written; see §6 |

The genotypes exist in two places with **byte-identical `.bed`** and different
`.fam` conventions: `rds-abcd/Data_Genetics/ABCD_hg19_allchrs_europeanonly.*`
(`FID=AB0000055 IID=NDAR_INVD5FWJDCY`) and `ABCD/genetics/abcd_eur.*`
(`FID=IID=sub-NDARINV…`). The latter is used, because every phenotype and
covariate file on this cluster is keyed the `sub-NDARINV` way, and because the
prebuilt GRMs are keyed to it. `genetics/abcd_eur.fam` itself is mode 0660 to a
group we are not in — hence `genetics_copy/`, which symlinks the `.bed`/`.bim`
and holds a readable `.fam`.

### Caveat on reusing the prebuilt GRM

The GRM was built by another user in Nov 2024 and its build log has since been
overwritten by a `--pca` run, so **the MAF and missingness filters behind it are
not recoverable**. It is the right N and the right subjects, and it is what the
sparse GRM was derived from, so it is fine for a dry run. If the published h²
needs documented provenance, run `01_grm` in full (~2 h) with `GRM` pointed at
`$ABCD_HPC_ROOT/results/grm/abcd_full`; the config makes that a one-line change.

Also present and unused: `ABCD/genetics/abcd_eur.eigenvec`, 50 PCs computed by
GCTA **within this EUR genotyped sample**. Those are arguably better ancestry
covariates than the 32 release PCs `gcta_export` pulls from `ab_g_stc`, which
were computed in the full multi-ancestry cohort. Not switched — flagging it.

## 5b. Provenance of the genotypes, and why only 5,678 subjects

The genotypes are **not** a subset of release 7.0. They are a 2020 imputation
with an ancestry filter on top, and the two losses compound.

### Dates

| file | date | what it is |
|---|---|---|
| `ABCD_chr{1..22,X}_hg19.*` | **2020-05-20 → 2020-07-30** | TOPMed-imputed, hg19, **10,072 subjects**, all ancestries. 60 GB, 23 chromosomes |
| `ABCD_hg19_allchrs_europeanonly.*` | log dated **2020-06-19** | the above merged, `--keep QC4_european_grm.grm.id`, `--maf 0.001` → **5,678 subjects × 13,697,177 SNPs** |
| `abcd_eur.*` | 2024-11-17 | byte-identical `.bed` to the previous row; only the `.fam` was rewritten to `sub-NDARINV` IDs (2024-11-22) |
| `Data_Genetics/ABCD_hg19_allchrs_europeanonly.*` | 2025-12-29 | another byte-identical copy of the same 2020 file |

The imputation log records `Working directory:
/mnt/b2/home4/arc/vw260/ABCD/ABCDgenotype/Genotype_postimputation/TOPMED/`,
i.e. it was produced elsewhere and copied in. So the 2024 and 2025 timestamps
are re-copies, not re-processing: there is one genotype dataset here, from 2020.

Content date versus arrival date, from `stat`: the per-chromosome files have
**mtime 2020-05-20/21** (`.bed`, `.fam`) and **2020-06-20 / 2020-07-30**
(`.bim`, i.e. the variant IDs were re-annotated two months later), but
**ctime 2025-02-06** — the day the inode was written on this filesystem. mtime
is the data's age; ctime is when it was copied here. Nothing was reprocessed.

### Which ABCD release these genotypes are

Cross-referencing against release 4.0's own genotyping manifest
(`Data_Phenotype/release4/genomics_sample03.txt`) identifies them precisely:

```
release-4 genomics manifest subjects   10,217
local 2020 imputed fileset             10,072      <- essentially the same cohort
7.0 genotyped (has genetic PCs)        11,663

local \ release-4 manifest                109
release-4 manifest \ local                254
```

So the local files are the **release-4-era genotyping**. And of the 1,627
subjects genotyped in 7.0 but absent locally:

```
  in the release-4 manifest               253   existed then; we lack them
  NOT in the release-4 manifest         1,374   genotyped after release 4
```

i.e. five-sixths of the gap is subjects who were simply **not yet genotyped in
2020**. This confirms the suspicion that prompted the check: these are old data
from before ABCD's genotyping was complete. In 7.0 it effectively is complete —
11,663 of 11,868 subjects (98.3 %), and 8,077 of our 8,192 phenotyped (98.6 %).

### Loss 1 — the genotype release is old

Release 7.0 tabulates genetic PCs and a genetic family ID for **11,663** of its
11,868 subjects; those are the genotyped ones. The local 2020 set has 10,072.

```
7.0 genotyped                    11,663
local 2020 imputation            10,072
  of which also in 7.0           10,036   (36 presumably withdrawn since)
7.0 genotyped, absent locally     1,627   <-- lost to the old release
```

### Loss 2 — the European-ancestry filter

`10,072 → 5,678` is 56.4 %. Checked against the release's own genetic PCs and
self-reported race, the filter is genuinely **genetic-ancestry-based, and
strict**:

| NIH race | in EUR set | not in EUR set | kept |
|---|---|---|---|
| 2 White | 5,490 | 1,324 | 80.6 % |
| 3 Black | 1 | 1,539 | 0.0 % |
| 8 Multiple | 147 | 947 | 13.4 % |
| 13 Other | 18 | 375 | 4.6 % |
| 4, 5, 6 (Asian, NHPI, AIAN) | 16 | 206 | 7.2 % |

and on `ab_g_stc__gen_pc__01` the kept group has mean +0.006 with SD **0.001**,
against −0.008 with SD 0.012 for the rest — a very tight cluster, an order of
magnitude less dispersed. That is a distance-to-EUR-centroid cut, not a QC
accident: it removes admixed individuals as well as non-European ones, which is
why a fifth of self-reported White subjects are also gone.

Caveat on certainty: the keep-list is named `QC4_european_grm.grm.id`, so it
sits downstream of an upstream QC chain whose logs are **not** on this cluster.
The numbers are consistent with ancestry doing essentially all of the work, but
a call-rate / heterozygosity / sex-check contribution cannot be excluded from
the files present.

### What it costs this analysis

Of the 8,192 phenotyped subjects:

| genotype source | usable N | note |
|---|---|---|
| 7.0 genotype release | **8,077** | the ceiling — not on this cluster |
| local 2020, all ancestries | **7,111** | `ABCD_chr*_hg19.*`, present and readable |
| local 2020, EUR only | **4,126** | what `config.local.sh` currently points at |

For precedent: a 2024 lab imaging GWAS
(`genetics/ABCD/GWAS/Imaging_genetics/2024_SubCorVol/*.log`) used exactly this
EUR fileset and finished at N = 4,986 after covariate matching, so ~4–5k is the
established working sample for this genotype set — not a mistake in our setup.

### The lever

**The all-ancestry per-chromosome files are right there and readable**, so
merging them raises the analysable sample from 4,126 to 7,111 — a 72 % increase,
larger than anything else available in this pipeline. What that costs:

- a fresh full GRM and sparse GRM over 10,072 subjects (hours, not minutes);
- `02_reml`'s SNP-h² becomes questionable in a multi-ancestry GRM, where
  ancestry differences enter the additive term. The conventional answers are to
  run REML within ancestry, or to keep REML on the EUR subset and use the full
  sample only for `03_gwas` (fastGWA + PCs handles structure) and `06_prs`;
- conversely it **repairs** `06_prs`: the "EUR primary, full multi-ancestry
  sensitivity" design currently has no full sample, so its `full` stratum is the
  EUR stratum relabelled.

No newer genotype release exists anywhere under `rds-abcd-CeXlNYOYMxw` or
`rds-genetics_hpc-Nl99R8pHODQ` — searched. Obtaining 7.0's genotype package from
NDA is the only route to the 8,077 ceiling.

## 6. Disorder GWAS needed normalising: `work/prepare_sumstats.sh`

Neither file has the shape `config.sh`'s defaults describe.

- **SCZ** `PGC3_SCZ_wave3.primary…vcf.fixed.tsv` is a PGC *sumstats VCF*: 73
  `##` metadata lines precede the header. Every consumer in `hpc/` reads the
  header from line 1 — MAGMA `--pval`, LDSC `munge_sumstats`, and `06_prs`'s
  awk — so all three would have taken `##fileFormat=PGCsumstatsVCFv1.0` as the
  header. It also has **no N column** (it carries `NCAS`, `NCON`, `NEFFDIV2`),
  against `04_magma`'s hard-coded `ncol=N` and config's `SCZ_N_SPEC="--N-col NEFF"`.
- **MDD** `pgc-mdd2025_…formatted.tsv` has the right columns but spells the
  sample size `n`, lowercase, against `ncol=N` / `--N-col N`.

`work/prepare_sumstats.sh` writes normalised copies to `work/sumstats/{SCZ,MDD}.tsv`
with a real header on line 1 and an `N` column that is the **effective** sample
size (SCZ: `2 × NEFFDIV2`), and `config.local.sh` points at those. 7,585,069 SCZ
rows and 5,918,521 MDD rows.

Gene-set covariates were exported with
`python -m abcd.magma_export --out-dir hpc/work/genesets`, which needs
`ABCD_AHBA_MAGMA_DIR=/rds/user/rajd2/hpc-work/magma`:

- `ahba_components_posneg.txt` — 7,643 genes × 6 covariates (C1–C3, ± halves), 0 dropped in Entrez mapping
- `snrnaseq_pc1.txt` — 17,618 genes × 4 covariates, 1 dropped

## 6b. Checked and found *not* to be a problem

Recorded because each looked like a defect and would otherwise be re-litigated.

- **GCTA and the header row.** `gcta_export` writes `FID IID <names>` headers,
  and GCTA's `--pheno`/`--qcovar`/`--covar` are conventionally headerless. GCTA
  1.94.1 reads the header as one extra individual (`5679`) and then drops it on
  the GRM intersection (`5678 individuals are in common`), for all three file
  types. Verified directly against the real GRM. No stripping needed.
- **`--mpheno` offset.** The manifest's `mpheno = column_index - 1` matches
  GCTA's convention of counting phenotype columns after FID/IID, confirmed by
  the same run.
- **SNP identifiers.** 11,826,979 of 13,697,177 `abcd_eur` variants (86.3 %)
  carry rsIDs, so the rsID-keyed reference data works for the bulk of the scan.
  See §11.4 for what the remaining 13.7 % costs.

## 6c. Defect found: a bare `python` on CSD3 is 3.7.4

`R/fit_lmm.R` resolves the run directory by shelling out to
`python -m abcd.run_dir`, deliberately, so the config hash has one
implementation rather than a second one in R. On CSD3 a bare `python` is the
system 3.7.4, which cannot even import the package —
`from typing import Literal` needs 3.8 — so the R step dies with

```
Error: could not resolve a run directory. Either pass --run-dir, or set
  export ABCD_CONFIG=ct_70_genetic
```

which names neither the real cause nor the interpreter. The same root cause
makes `make help` fail, and with it the two repo tests noted in §8. Fixed by
exporting `PATH="$W/envs/abcd/bin:$PATH"` before the R call in
`work/fit_pheno.sbatch`.

Worth hardening upstream: `fit_lmm.R` could report the python it invoked and
its stderr's first line, and the `Makefile` could honour a `PYTHON` variable.

## 6d. Defect found and fixed: `run_all.sh`'s SLURM mode sourced nothing

The most dangerous defect found in this exercise, because it **succeeds**.

Every step begins `source "$(dirname "$0")/config.sh"`. That is correct when
run as `bash hpc/02_reml.sbatch`, which is how the local harness and the smoke
test invoke them. Under `sbatch` it is not: SLURM runs a **copy** of the script
at `/var/spool/slurm/slurmd/jobNNN/slurm_script`, so `dirname "$0"` is the spool
directory and `config.sh` is not there.

The failure is silent in the worst way, because `set -euo pipefail` lives
*inside* `config.sh` and so is never enabled either. Each job then runs to the
end with every function undefined:

```
/var/spool/slurm/slurmd/job33429747/slurm_script: line 19: .../config.sh: No such file or directory
line 20: require_paths: command not found
line 21: ensure_dirs: command not found
line 28: phenotype_names: command not found
array index 1 exceeds 0 phenotypes; nothing to do
```

and exits **0**. The first full-pipeline submission returned 13 jobs all marked
`COMPLETED` by `sacct`, in 1–4 seconds each, having produced nothing. Anyone
trusting the job state rather than the outputs would have concluded the pipeline
had run.

Fixed in all seven scripts plus `run_all.sh`: `set -e` first, then search
`$(dirname "$0")`, `$PWD` (which `run_all.sh` already sets correctly via
`--chdir`) and `$SLURM_SUBMIT_DIR`, and `exit 2` with a named error if none has
`config.sh`. Verified from the repo root, from `hpc/`, by absolute path from
`/tmp`, and against a copy placed in a spool-like directory.

## 7. SLURM: resource requests that do not fit the partition

`icelake` has `MaxMemPerCPU=3370` MB (`cclake` had 3410, `icelake-himem` has
6760). A `--mem` above `cpus × 3370` is not granted as asked — SLURM widens the
CPU request or rejects the job. Every script except `04_magma` asks for more
memory per CPU than icelake allows:

| script | requests | icelake cap | verdict |
|---|---|---|---|
| `01_grm` | 16 cpu, 64 G | 52.7 G | **over** |
| `02_reml` | 8 cpu, 32 G | 26.3 G | **over** |
| `03_gwas` | 16 cpu, 64 G | 52.7 G | **over** |
| `04_magma` | 16 cpu, 32 G | 52.7 G | ok |
| `05_ldsc` | 4 cpu, 16 G | 13.2 G | **over** |
| `06_prs` | 8 cpu, 32 G | 26.3 G | **over** |

Measured rather than assumed: a probe job asking `--cpus-per-task=4 --mem=16G`
on icelake came back with `ReqTRES=cpu=4,mem=16G` but
`AllocTRES=cpu=5,mem=16G,billing=5`. So SLURM **silently widens the CPU
allocation** (4 → 5, because 5 × 3370 ≥ 16 G) rather than rejecting. The
consequence is not a failure, it is a **billing surprise**: every step above
quietly costs more CPU-hours than its header asks for, and `05_ldsc`'s 4 cpu /
16 G is charged as 5.

Either drop `--mem` and let `cpus × 3370` supply it, or submit to
`icelake-himem`. Not yet changed in the tracked scripts — it is a site fact, and
`run_all.sh` already overrides `--account`/`--partition` at submit time, so
`--partition icelake-himem` is available without editing them. The `work/*.sbatch`
helpers written here omit `--mem` for exactly this reason.

### What the full run will actually cost

fastGWA is far cheaper than the script headers imply. The lab's 2024 ABCD GWAS
over this same fileset — 13,697,177 SNPs, 4,986 subjects — finished in
**2 min 9 s on a single thread**. So `03_gwas` genome-wide over 5 phenotypes is
minutes, not the 8 h requested, and the full run's cost is dominated by MAGMA's
genome-wide SCZ/MDD gene analyses and by PLINK clumping.

**SL3 caps wall time at 12 h** (QOS `cpu2`, `MaxWall=12:00:00`, 448 CPUs per
user). Every tracked script is within that — `01_grm`'s 12 h request is exactly
at the limit — but it is a smaller ceiling than SL2 had.

Separately, `04_magma.sbatch` is the only step that still hard-codes
`#SBATCH -A VERTES-SL2-CPU`, `-p cclake` and `-o slurm/…` in its header — a
leftover from the pre-rewrite version. Under `run_all.sh` the command-line
`--account` wins, so it works; submitted directly with `sbatch 04_magma.sbatch`
it fails on the exhausted account and on the missing relative `slurm/` directory.

## 8. Local phenotype pipeline — reproduces the report

The 7.0 release is the tabulated tree at
`rds-abcd-CeXlNYOYMxw/derivatives/tabulated`, symlinked to
`work/abcd_root/abcd-7.0` so `paths.release_dir` finds it with
`ABCD_ROOT=hpc/work/abcd_root`. Every table the 7.0 adapter needs is present
(`mr_y_smri__thk__dsk`, `ab_g_stc`, `ab_g_dyn`, `mr_y_qc__incl`,
`mr_y_qc__post__aut`, `gn_y_genrel`).

`python -m abcd.assemble` under `ABCD_CONFIG=ct_70_noglobal_mv2_genetic`:

```
imaging loaded      rows=2,089,044  subjects=11,818
after QC            rows=1,748,667  subjects=10,126
sample rule applied rows=1,591,812  subjects=8,192
final: 1,591,812 rows, 8,192 subjects, 68 regions, 19 sites, 6,871 families
```

8,192 subjects is exactly the figure `REPORT_7.0` quotes, so this tree is the
release the report was built from.

`out/` is a symlink to `hpc/work/out` so the pipeline's own outputs stay inside
`hpc/` as required, without patching `paths.py`.

`tests/test_docs_provenance.py`, `test_config.py` and `test_readme.py` were run
after the `config.sh` edit: 28 pass, 2 fail. Both failures are environmental,
not caused by the edit — `test_make_help_resolves_a_run_dir` and
`test_readme_test_count_is_current` both shell out to `make help`, which invokes
a bare `python`, and on CSD3 that is the system 3.7.4 without the `abcd`
package. Prefix with `PATH=hpc/work/envs/abcd/bin:$PATH` to run them here.

## 9. Current state of the run

| job | what | state |
|---|---|---|
| `abcd_chr22` | chr22 subset of `abcd_eur` for the smoke test | **done** — 187,614 SNPs × 5,678 subjects, 4 s |
| `abcd_fit` | `fit_lmm` → `phenotype` → `gcta_export` → `align_ids` | **done** in 78 s (2nd attempt; 1st died on the bare-`python` bug, §6c) |
| `abcd_smoke` | `work/smoke_test.sbatch`, steps 01–06 | **done** — job 33414727, results below |
| full run, 1st attempt | jobs 33429667–33429671 | **void** — all 13 exited 0 having done nothing (§6d) |
| full run, 2nd attempt | jobs 33430067–33430071 | **done** for 02–05; `06_prs` failed (§12b) |
| `abcd_prs` rerun | job 33447215, after the §12b fix | **done** — 160 models fitted |

**Genome-wide results are in §12. The pipeline is now complete end to end.**

### `abcd_fit` output

68 regions fitted with `value ~ sex + age_c + (1 + age_c | subject) + (1 | site)`
over 1,591,812 rows: **0 singular, 0 non-converged, 0 errors**, 0.5 min on 16
cores. Reliabilities by exact visit count reproduce the report's pattern —
intercept 0.892 / 0.925 / 0.941 and slope 0.109 / 0.175 / 0.243 at 2 / 3 / 4
visits.

`gcta_export` wrote 8,192 rows × 5 phenotypes, with 10 ancestry PCs and 115
subjects missing them. `align_ids.py` then kept **4,126** and dropped 4,066 as
not genotyped — the EUR-subset cost quantified in §5b.

### Preflight

`00_check_inputs.sh` **passes**, including the manifest-consistency check
(`baseline_thickness → --mpheno 1` … `slope_PC1 → --mpheno 5`) and an ID
intersection of 4,126 of 4,126.

### Smoke test (job 33414727) — results

| step | result |
|---|---|
| `00_check_inputs` | pass |
| `01_grm` DRY_RUN | pass — command construction valid |
| `01_grm --grm-cutoff` | pass — **4,608 unrelated of 5,678** (1,070 removed) in 0.8 s |
| `02_reml` | pass — real h² for all five phenotypes, see below |
| `03_gwas` | **FAIL** — see §9b |
| `04_magma` | pass — full path incl. genome-wide SCZ and MDD gene analyses |
| `05_ldsc` | munging passed; `--rg` **hung** — see §9c |
| `06_prs` | not reached |

#### `02_reml` — the first real SNP-h² numbers in this project

These are genuine: the GRM is genome-wide (only the *association scan* was
restricted to chr22). N = 3,329 — the 4,126 phenotyped-and-genotyped subjects
intersected with the 4,608 unrelated and with complete covariates.

| phenotype | h² | SE | p | n |
|---|---|---|---|---|
| `baseline_thickness` | **0.575** | 0.147 | 4.3e-05 | 3,329 |
| `global_slope` | 0.228 | 0.141 | 0.051 | 3,329 |
| `slope_PC2` | 0.113 | 0.141 | 0.214 | 3,329 |
| `slope_PC3` | 0.000 | 0.141 | 0.500 | 3,329 |
| `slope_PC1` | 0.000 | 0.139 | 0.500 | 3,329 |

**The positive control works.** `baseline_thickness` returns SNP-h² = 0.58
(p = 4e-5), comfortably below the report's Falconer twin estimate of ~0.74–0.8
and above zero — which is exactly the expected relationship, since SNP-h²
captures only common-variant additive variance. Had this come back near zero or
above 1, the GRM, the covariates or the ID alignment would be wrong.

The developmental phenotypes behave as `README.md`'s power section predicts:
`global_slope` at 0.23 ± 0.14 is a z of 1.6 — suggestive, not significant — and
the slope PCs are indistinguishable from zero. With N = 3,329 rather than the
~8,000 assumed, and slope reliability ~0.11–0.24, this is the expected outcome
rather than a failure. It does mean §11's power caveat is now measured, not
predicted: **LDSC rg on the slope phenotypes cannot be interpretable**, since
their h² z-scores are ≈ 0–1.6 against LDSC's guidance of z > 4.

## 9b. Defect found and fixed: fastGWA's variance floor

`03_gwas` completed `baseline_thickness` and then died on `global_slope`:

```
Error: the Vp is below 1e-5. Please check: 1. Is there a scaling issue with the
phenotype? ...
```

fastGWA refuses any phenotype with phenotypic variance under 1e-5.
`global_slope` is a per-subject thinning rate in mm/year, and its variance is
**1.37e-06** — not degenerate, just measured in units that make it tiny. The
five exported phenotypes span eight orders of magnitude:

| phenotype | SD | variance | fastGWA |
|---|---|---|---|
| `baseline_thickness` | 0.0614 | 3.77e-03 | ok |
| `global_slope` | 0.00117 | **1.37e-06** | **refuses** |
| `slope_PC3` | 15.34 | 2.35e+02 | ok |
| `slope_PC2` | 16.39 | 2.69e+02 | ok |
| `slope_PC1` | 16.13 | 2.60e+02 | ok |

GCTA-REML has no such floor, which is why `02_reml` produced an estimate for
`global_slope` that `03_gwas` could not even start. Note that this would have
silently killed **the project's primary phenotype** — `README_7.0` ranks the
global mean slope first for GWAS.

`work/standardise_pheno.py` z-scores every phenotype column, keeping a
`.prescale` backup and printing the SDs needed to convert a BETA back to native
units. This is safe rather than a fudge: a linear model's p-values are invariant
under an affine transform of the outcome (BETA and SE rescale together), and h²
is a variance ratio, so nothing changes except that BETA is now per phenotype-SD
— a more useful unit anyway across five phenotypes whose native scales differ by
10⁸.

**Worth fixing upstream:** `abcd.gcta_export` should standardise on export, or
at minimum warn when a phenotype's variance is near fastGWA's floor.

## 9c. LDSC `--rg` hung on a filesystem stall

`05_ldsc` munged all three inputs correctly — SCZ and MDD genome-wide gave
1,039,708 HM3 SNPs for MDD at mean χ² = 2.244, λ_GC = 1.892, 3,177
genome-wide-significant SNPs, which is right for a well-powered MDD GWAS — and
then `ldsc.py --rg` produced a **zero-byte log and made no progress for two
hours**.

Diagnosed on the compute node rather than guessed at: the process sat at **0 %
CPU in state `I`**, blocked in an `fstat` (syscall 5) on an open fd pointing at
`/rds/user/rajd2/hpc-work/ldsc/eur_w_ld_chr/9.l2.ldscore.gz`, having read only
35 MB in 3,588 read syscalls. That is a Lustre stall on the shared LD-score
directory, not a code fault — the identical call completed in 10 s from a login
node earlier in the session.

Mitigation: `work/ldsc_ref/eur_w_ld_chr` is our own copy (45 MB, verified
byte-identical), and `LD_REF`/`LD_WEIGHTS`/`HM3_SNPLIST` now point at it. Same
filesystem, but a fresh inode and no sharing. If it recurs, note that
`run_all.sh` makes 04, 05 and 06 siblings, so a hung 05 does not block the other
two.

SL3 is a lower-priority QOS (`cpu2`, priority factor 0 fairshare), so queue
waits are long — 30 min for a 4-second job during this session. Keep wall-time
requests tight so jobs can backfill.

The helper scripts written for this run all live in `work/` and are all
site-specific, not pipeline stages:

| script | does |
|---|---|
| `work/prepare_sumstats.sh` | normalise the SCZ/MDD sumstats (§6) |
| `work/make_chr22_subset.sbatch` | build the smoke-test genotypes |
| `work/fit_pheno.sbatch` | `fit_lmm` → `phenotype` → `gcta_export` → `align_ids` |
| `work/align_ids.py` | FID/IID alignment (§4) |
| `work/smoke_test.sbatch` | run steps 01–06 in one job, continuing past failures |

### Smoke-test design

Chromosome 22 only, **all 5,678 subjects**. A chromosome rather than a subject
subsample because the GWAS cost is linear in SNPs and its correctness does not
depend on which ones, whereas dropping subjects would change the relatedness
structure and leave the sparse-GRM path untested. REML runs against the real
genome-wide GRM, so its h² estimates are real numbers, not fixture numbers;
only the association scan is restricted.

## 10. Next commands

```bash
cd /home/rajd2/rds/hpc-work/abcd_development

# preflight
bash hpc/00_check_inputs.sh

# rebuild the PRS-side phenotypes if work/pheno/ is regenerated (§12b)
hpc/work/envs/abcd/bin/python hpc/work/make_prs_pheno.py \
    hpc/work/pheno hpc/work/pheno_prs

# rerun one step (02..06); SLURM mode, afterok-chained
SBATCH_TIMELIMIT=08:00:00 bash hpc/run_all.sh 06

# full chain
SBATCH_TIMELIMIT=08:00:00 bash hpc/run_all.sh
```

**Verify outputs, not `sacct`.** §6d's failure mode was 13 jobs reporting
`COMPLETED` while producing nothing. Every step writes a `*_summary.tsv`; check
that it exists and has the expected row count before believing a green state.

## 11. Open questions for the next person

1. **The sample is EUR-only, and smaller than the power section assumes.**
   Only 5,678 of the 8,192 phenotyped subjects are genotyped, and the sole
   genotype fileset on this cluster — in both the `rds-abcd` and
   `rds-genetics_hpc` copies — is `europeanonly`. Two consequences:
   - `README.md`'s "N≈8,000" runs through the whole power argument (the h² z≈2
     estimate, the "0.5 % of variance is detectable" PRS claim). At N≈5,000, and
     less after the unrelated cut for REML, both are optimistic.
   - `06_prs`'s ancestry design — EUR primary, full multi-ancestry sample as
     sensitivity — has **no full sample to contrast against**. `prs_assoc.R`
     will fit its "full" stratum, but it is the EUR stratum with a different
     label. Either source multi-ancestry genotypes or drop that arm and say so.
2. **Ancestry PCs**: ~~release PCs or the in-sample EUR PCs already sitting in
   `abcd_eur.eigenvec`?~~ **RESOLVED — see §15.** Two corrections to what this
   entry used to say: there is no `abcd_eur.eigenvec` (it never existed, so
   in-sample PCs had to be computed), and although the release PCs are indeed
   nearly constant within EUR as suspected, swapping them for in-sample PCs
   moves h² by 0.009. The concern was real and the effect is negligible.
3. **GRM provenance**: reuse the undocumented prebuilt GRM, or spend 2 h
   rebuilding it with recorded filters? §5.
4. **SNP identifiers**: 86.3 % of `abcd_eur.bim` IDs are rsIDs; the other
   1,870,198 are `chr1:66861:C:T`-style. Those cannot match the g1000 reference,
   `w_hm3.snplist` or PGC weights, so MAGMA, LDSC and PRS all silently work off
   the rsID subset. That is acceptable — but it should be stated in the results
   rather than discovered later.

## 12. Genome-wide results (jobs 33430067–33430071)

Genotypes `abcd_eur` (13,697,177 SNPs × 5,678), phenotypes as z-scores (§9b),
**N = 4,119** after intersecting phenotypes, covariates and genotypes. All five
phenotypes are the standardised BLUPs, so every β below is per SD.

### `02_reml` — SNP-h², genome-wide GRM (N = 3,329 unrelated)

Unchanged from the smoke test, and correctly so: the smoke test already used the
genome-wide GRM and only restricted the *association scan*. Reproduced exactly.

| phenotype | h² | SE | p |
|---|---|---|---|
| `baseline_thickness` | **0.575** | 0.147 | 4.3e-05 |
| `global_slope` | 0.228 | 0.141 | 0.051 |
| `slope_PC2` | 0.113 | 0.141 | 0.214 |
| `slope_PC3` | 0.000 | 0.141 | 0.500 |
| `slope_PC1` | 0.000 | 0.139 | 0.500 |

### `03_gwas` — fastGWA-MLM, 8,777,886 SNPs after MAF > 0.01 and missingness < 0.10

| phenotype | λ_GC | min p | SNPs p<1e-5 | SNPs p<5e-8 |
|---|---|---|---|---|
| `baseline_thickness` | 1.032 | 4.4e-07 | 216 | **0** |
| `global_slope` | 1.006 | 8.1e-07 | 155 | **0** |
| `slope_PC2` | 1.023 | 2.5e-07 | 82 | 0 |
| `slope_PC1` | 1.013 | 2.4e-07 | 75 | 0 |
| `slope_PC3` | 1.000 | 1.0e-07 | 87 | 0 |

**No genome-wide-significant loci for any phenotype**, which at N = 4,119 is the
expected result rather than a surprise — the largest published cortical-thickness
GWASs need tens of thousands of subjects for their first hits. λ_GC of 1.00–1.03
says the mixed model is properly calibrated: no residual stratification or
cryptic relatedness inflating the null.

### `05_ldsc` — h² and genetic correlation with SCZ and MDD

The LDSC h² is the independent check that matters most here:

| phenotype | LDSC h²_obs | SE | z | REML h² | interpretable? |
|---|---|---|---|---|---|
| `baseline_thickness` | 0.584 | 0.131 | **4.45** | 0.575 | **yes** |
| `slope_PC2` | 0.340 | 0.103 | 3.30 | 0.113 | underpowered |
| `global_slope` | 0.165 | 0.114 | 1.45 | 0.228 | underpowered |
| `slope_PC1` | 0.084 | 0.102 | 0.82 | 0.000 | underpowered |
| `slope_PC3` | −0.243 | 0.088 | −2.75 | 0.000 | no (h² < 0) |

`baseline_thickness` gives **0.584 from LDSC against 0.575 from GCTA-REML** —
two methods with different assumptions, different samples (3,329 unrelated vs
4,119 all) and different estimators, agreeing to within 0.01. That is a strong
joint validation of the GRM, the ID alignment, the covariates and the sumstats
munging all at once.

| disorder | phenotype | rg | SE | p |
|---|---|---|---|---|
| SCZ | `baseline_thickness` | 0.036 | 0.047 | 0.45 |
| MDD | `baseline_thickness` | −0.019 | 0.049 | 0.70 |

**No genetic correlation with either disorder.** Only the `baseline_thickness`
rows are interpretable; every other row is flagged `underpowered` in
`ldsc_rg_summary.tsv` because its h² z is below LDSC's z > 4 guidance, and the
two `slope_PC3` rows are `NA` because its h² estimate is negative. This is §9's
predicted outcome, now measured.

### `04_magma` — gene-set and gene-property tests

The AHBA and snRNA-seq gene-set tests are null across all five phenotypes after
accounting for the number of tests: the smallest are `C3+` on `baseline_thickness`
(β = 0.145, p = 0.044) and `C2+` on `slope_PC1` (β = 0.123, p = 0.026), neither
of which survives correction for 6 components × 5 phenotypes.

The reverse-direction tests are the one place a genome-wide signal appears:

| test | β | SE | p | genes |
|---|---|---|---|---|
| `MDD_vs_baseline_thickness` | 0.047 | 0.0113 | **2.8e-05** | 18,239 |
| `SCZ_vs_baseline_thickness` | 0.025 | 0.0113 | 0.027 | 18,355 |
| all 8 slope/global tests | — | — | 0.14–0.74 | — |

**Read this carefully — it does not contradict the null rg.** MAGMA gene-level
`ZSTAT` is *unsigned*: it measures how strongly a gene is associated, not in
which direction. So this says genes with stronger thickness signal tend to have
stronger MDD signal, which is compatible with an LDSC rg of −0.019, because a
signed correlation averages to zero when the sharing is mixed in direction.

It is not a gene-size artefact: the `.gsa.out` header records
`CONDITIONED_INTERNAL = gene size, gene density, sample size, inverse mac` plus
logs of each. But `BETA_STD = 0.048` is a very small effect, it is one test in a
family of ten, and it is unreplicated by the better-powered method on the same
data. Treat it as a lead, not a finding.

### `06_prs` — polygenic scores (job 33447215, after the §12b fix)

**160 models fitted**: 2 disorders × 8 thresholds × 5 phenotypes × 2 strata.
EUR n = 3,725 in 3,103 families; "full" n = 4,126 in 3,455 families.

**Nothing survives correction — zero of 160 rows have `p_adj` < 0.05.** The
largest partial R² anywhere in the table is 0.0015.

The one pattern worth recording is SCZ → `global_slope`, because it is coherent
across thresholds rather than a single-threshold blip:

| p-threshold | SNPs | β (EUR) | p | β (full) | p |
|---|---|---|---|---|---|
| 5e-8 | 416 | +0.005 | 0.76 | +0.000 | 0.97 |
| 1e-5 | 1,465 | −0.001 | 0.96 | −0.003 | 0.87 |
| 0.001 | 8,496 | −0.006 | 0.69 | −0.007 | 0.62 |
| 0.01 | 28,843 | −0.018 | 0.26 | −0.020 | 0.19 |
| 0.05 | 74,479 | −0.027 | 0.081 | −0.029 | 0.056 |
| 0.1 | 114,004 | −0.026 | 0.095 | −0.027 | 0.078 |
| **0.5** | 291,110 | **−0.038** | **0.015** | **−0.040** | **0.009** |
| 1 | 390,325 | −0.036 | 0.024 | −0.037 | 0.014 |

`global_slope` is a thinning rate with raw mean −0.0189 mm/year, so a *negative*
β means **higher SCZ polygenic score → faster cortical thinning** — the
direction the literature predicts. The effect grows monotonically as the
threshold admits more of the polygenic tail, which is the signature of a real
polygenic signal rather than a few noisy top SNPs.

Against that: `p_adj` = 0.119 after Bonferroni across the 8 thresholds, and the
threshold correction is the honest one to apply since the threshold is a
nuisance choice, not a hypothesis. β = −0.038 SD per SD is 0.15 % of variance.
**This is a lead to power up, not a result.** Note also that it is the same
direction as, but independent of, the null LDSC rg — PRS association and rg are
differently powered, and at h² z = 1.45 for `global_slope` the rg was never going
to be informative.

The "full" stratum is **not** the multi-ancestry sensitivity analysis it was
designed to be — both strata come from the EUR-only genotypes, so "full" is the
EUR sample plus the 401 subjects outside the PC-space cut. See §11.1.

## 12b. Defect found and fixed: the family random effect was destroyed

`06_prs` clumped and scored both disorders across all 8 thresholds correctly —
16 `.profile` files — then died at the last step:

```
EUR subset: 3725 of 4126 subjects [PC-space cut: within 3.0 SD of median ...]
Error: no models fitted -- check ID overlap between scores and phenotypes
```

**The diagnosis in the error message is wrong**, which is what makes this one
worth recording. The ID overlap was complete: 4,126 of 4,126.

`tools/prs_assoc.R` fits `phenotype ~ scale(PRS) + ... + (1 | family_id)` and
sets `family_id := FID`, which is correct for `abcd.gcta_export`'s output, where
FID *is* the family id. But `align_ids.py` (§4) must overwrite FID with the
`.fam`'s spelling — where FID = IID — or GCTA matches nobody. After that rewrite
every subject is a family of one, `lmer` refuses the model ("number of levels of
each grouping factor must be < number of observations"), `prs_assoc.R`'s `try()`
skips it, and the loop ends with zero rows and a misleading message.

This mattered rather than being cosmetic: **657 of the 3,455 families in the
analysis sample contribute more than one child** (645 pairs, 11 triples, 1
quintuple). Dropping the random effect would have understated the standard error
on every PRS coefficient.

Fixed without touching the shared `tools/prs_assoc.R`, since the two consumers
genuinely want different FID columns:

- `work/align_ids.py` now writes **`family_map.tsv`** before overwriting FID.
- **`work/make_prs_pheno.py`** (new) writes `work/pheno_prs/` — the same three
  files with FID restored to the family id and IID untouched, so the merge
  against PLINK's `.profile` (which is by IID) is unaffected.
- `06_prs.sbatch` uses `pheno_prs/` when present and warns loudly when it is
  not, instead of failing 13 minutes of clumping later with the wrong reason.

Rerun as job 33447215: **160 models fitted**, 3,103 families in the EUR stratum.

### Also noted, not fixed

`prs_assoc.R` asks for a covariate named `age_c`; `covar_quant.txt` supplies
`baseline_age`. The script's `terms[terms %in% names(sub)]` filter silently drops
the missing name, so **the PRS models are not age-adjusted**. The impact is
small — the phenotypes are BLUPs from a model that already conditions on
`age_c`, so the age variance is largely gone before this step — but it is a
silent drop of a requested covariate and should be reconciled upstream.

## 13. SCZ prioritised-gene enrichment (jobs 33456034, 33456616)

Added after the main run, at the request to test the ABCD phenotype GWAS against
the **prioritised gene set of Trubetskoy et al. 2022** (Nature 604:502) rather
than only the AHBA / snRNA-seq sets of `04_magma`.

Source: figshare 19426775 file 35775617, `scz2022-Extended-Data-Table1.xlsx`,
md5 `8fee2faee10ddf2c9f8a16a54f2c5b84` (verified on download), kept at
`work/genesets/`. Sheet `Extended.Data.Table.1` is the 120 prioritised genes;
sheet `ST12 all criteria` is the 685-gene candidate pool they were chosen from.

Built by `work/make_prioritised_genesets.py`, run by `work/prioritised_gsa.sbatch`.
Both reuse the `.genes.raw` files `04_magma` already produced, so no gene
analysis was repeated — only the set step, which takes seconds.

### Why this is a fair test when the genome-wide rg was null

LDSC rg averages over ~18,000 genes, most carrying no SCZ signal at all, and
asks about *signed* sharing. A competitive gene-set test asks whether ABCD
signal is *concentrated* in one set chosen by orthogonal evidence — one degree
of freedom, and far better powered at N = 4,119. A null rg does not predict a
null here.

### ID mapping

104 of 120 prioritised genes map to Entrez IDs in `NCBI37.3.gene.loc` — **98 % of
the protein-coding ones**. The 16 losses are lincRNAs and clone-named
transcripts (`RP11-…`, `LINC01088`); NCBI37.3 contains 3 `LINC*` entries in
total, so those genes were never annotated in any of these GWAS and could not
have entered a gene-set test however they were spelled. One recoverable alias
(`GPR98` → `ADGRV1`) is handled explicitly.

### Sets

| set | n mapped | what it is |
|---|---|---|
| `SCZ_prioritised` | 104 | the 120 prioritised genes |
| `SCZ_prio_finemap` / `_smr` / `_rare` | 64 / 44 / 7 | the three evidence streams separately |
| `SCZ_locus_pool` | 462 | the 685-gene candidate pool |
| `SCZ_pool_not_prio` | 358 | **pool minus prioritised** — disjoint, so the comparison means something |

MHC is not a confound: the Trubetskoy tables exclude it, verified — zero genes
from any set fall in chr6:25–34 Mb. MAGMA additionally conditions internally on
gene size, gene density, inverse MAC and their logs. `--model direction=greater`
throughout, since the hypothesis is concentration *in* these genes.

### Controls

| trait | set | β | SE | p | |
|---|---|---|---|---|---|
| SCZ | `SCZ_prioritised` | 2.565 | 0.131 | 4.7e-85 | plumbing OK |
| SCZ | prioritised \| pool | 0.613 | 0.142 | 8.5e-06 | |
| MDD | `SCZ_prioritised` | 1.000 | 0.133 | 2.5e-14 | **independent** |
| MDD | prioritised \| pool | 0.455 | 0.149 | 1.2e-03 | |

Our `SCZ.tsv` **is** `PGC3_SCZ_wave3.primary.autosome.public.v3` — the same GWAS
the 120 genes were prioritised from — so the SCZ row is circular by construction
and tests only that the Entrez mapping and `--set-annot` work. MDD
(`pgc-mdd2025`, 415,568 cases) is the real control, and it matters: conditional
on the whole locus pool the prioritised genes still add signal (p = 1.2e-03), so
**the conditional contrast demonstrably has power to detect prioritisation
adding information** when it is there.

### Result: the enrichment is in the locus pool, not the prioritised genes

30 marginal tests (6 distinct sets × 5 phenotypes), Bonferroni 1.67e-03.

| phenotype | set | β | SE | p | |
|---|---|---|---|---|---|
| `baseline_thickness` | `SCZ_pool_not_prio` | 0.207 | 0.058 | **1.9e-04** | survives |
| `baseline_thickness` | `SCZ_locus_pool` | 0.183 | 0.051 | **1.9e-04** | survives |
| `baseline_thickness` | `SCZ_prioritised` | 0.067 | 0.091 | 0.23 | — |
| `baseline_thickness` | prioritised \| pool | −0.114 | 0.103 | 0.87 | — |
| `global_slope` | `SCZ_pool_not_prio` | 0.142 | 0.058 | 0.0074 | nominal only |
| `global_slope` | `SCZ_locus_pool` | 0.109 | 0.051 | 0.017 | nominal only |
| **`global_slope`** | **`SCZ_prioritised`** | **−0.006** | **0.091** | **0.53** | — |
| `global_slope` | prioritised \| pool | −0.128 | 0.104 | 0.89 | — |
| `slope_PC1/2/3` | every set | — | — | 0.03–0.98 | nothing survives |

**The same pattern in both phenotypes**: the 356 *non-prioritised* genes at SCZ
loci carry whatever enrichment exists, and the 104 prioritised genes carry none —
their conditional estimate is slightly negative in both cases.

### What that does and does not license

It does **not** license "prioritisation fails for cortical phenotypes". The
within-trait contrast between prioritised and pool-remainder is not significant
(for `global_slope`, difference 0.148 against a pooled SE of ~0.108, two-sided
p ≈ 0.2). We cannot distinguish "prioritised genes are less enriched" from
"we lack the power to see their enrichment".

It also does not license comparing ABCD's `BETA_STD` to MDD's. `BETA_STD` is
standardised on the *observed* gene-z scale, and a trait with little genetic
signal has gene z-scores that are mostly noise, so its achievable `BETA_STD` is
bounded by its own power. The MDD benchmark is not a biological yardstick for
ABCD.

The internally valid statement is the one above: within a single phenotype, on
one scale, at one power level, the enrichment sits in the pool remainder.

### Power on `global_slope` specifically

Asked for explicitly, and it is a clean null: β = −0.006 ± 0.091, 95 % CI
[−0.18, 0.17], two-sided p ≈ 0.95. That CI contains 0.183 — the effect
`baseline_thickness` shows — so it is uninformative rather than evidence of
absence. From the LDSC logs, mean χ² per SNP:

| | mean χ² | excess | h² z |
|---|---|---|---|
| SCZ (discovery) | 2.124 | 1.124 | — |
| `baseline_thickness` | 1.038 | 0.038 | 4.45 |
| `global_slope` | 1.015 | **0.015** | 1.45 |

`global_slope` carries ~40 % of the per-SNP polygenic signal `baseline_thickness`
does. A competitive test can only concentrate signal that exists in the gene
statistics; there is very little to concentrate. Same N = 4,119 constraint,
third appearance.

### Limitation neither pass fixes

SCZ loci are brain-expressed genes and cortical thickness is a brain phenotype,
so some enrichment is expected on those grounds alone. Calling the locus-pool
result *SCZ-specific* needs a negative control set — genes at loci of a
well-powered non-brain trait. That has not been built.

### Defect found and fixed: the conditional collector shifted every column

The first collector read `.gsa.out` by fixed column position. The conditional
run inserts a `MODEL` column the marginal run does not have, so every
conditional row shifted one column left: `SE` was reported as `P` and the real
`P` was dropped entirely. The summary table looked complete and was silently
wrong for half its rows — the conditional p-values in it were standard errors.

Fixed by parsing the header line for column names (`for (i=1;i<=NF;i++) c[$i]=i`)
instead of assuming a layout, and by skipping all `#` lines rather than a fixed
`tail -n +4` — the conditional output has five comment lines, one more than the
marginal, so the extra one was also being emitted as a data row.

## 14. MDD high-confidence gene enrichment (job 33460031)

The MDD half of §13, using the 308 high-confidence genes of the 2024 Cell
depression GWAS (doi 10.1016/j.cell.2024.12.002, PMC11829167), Table S21 sheet
`High-confidence Gene List`.

**Obtaining the file.** cell.com serves a Cloudflare challenge to `curl` and
pmc.ncbi.nlm.nih.gov serves a JavaScript proof-of-work challenge, so neither
direct link works from the cluster. The **Europe PMC `supplementaryFiles` API**
(`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11829167/supplementaryFiles`)
returns the whole supplementary set as a zip with no challenge — the right route
for any future paper. Kept at `work/genesets/mdd2024_TableS21.xlsx`; built by
`work/make_mdd_genesets.py`; run by the now-parameterised
`work/prioritised_gsa.sbatch` (`SETS=`, `TAG=`, `COND_SET=`, `OUT=`).

### Two asymmetries with the SCZ analysis

**Mapping is worse.** 213 of 308 map to Entrez in NCBI37.3 — **69 %**, against
87 % for SCZ. The MDD table uses modern Ensembl-era symbols with many antisense
and lincRNA entries (`RERE-AS1`, `LINC01134`) that the NCBI37.3 build predates.

**The pools are not comparable.** SCZ's pool is 685 genes at 239 genome-wide
loci — "which genes sit at these loci". MDD's Table S21 is 4,600 genes (1,884
mapped, a tenth of everything MAGMA tests) because two of its seven mapping
methods (fastBAT, H-MAGMA) are themselves genome-wide gene-based association
tests, not locus annotations. A competitive test against a set that large is a
weak contrast. **The 308-gene high-confidence set is the like-for-like
comparator to SCZ's 120 prioritised genes; the pools are not**, and the MDD pool
rows below should not be read against §13's.

### Controls

| trait | set | β | SE | p | |
|---|---|---|---|---|---|
| MDD | `MDD_highconf` | 2.116 | 0.081 | 1.6e-146 | circular |
| MDD | highconf \| pool | 1.309 | 0.083 | 1.5e-55 | circular |
| SCZ | `MDD_highconf` | 0.802 | 0.086 | **5.4e-21** | **independent** |
| SCZ | highconf \| pool | 0.392 | 0.090 | **6.8e-06** | **independent** |

The MDD rows are circular (these genes were derived from this GWAS). The SCZ
rows are the meaningful control and they are strongly positive, including
conditional on the 1,884-gene pool — so the design has power to detect
high-confidence genes adding information beyond pool membership, on a trait
that did not define them.

### Result: null across every ABCD phenotype

30 tests, Bonferroni 1.67e-03. Nothing approaches it.

| phenotype | `MDD_highconf` (213) | conditional on pool | best any MDD set |
|---|---|---|---|
| `baseline_thickness` | β 0.062, p 0.16 | β 0.048, p 0.23 | 0.11 (`hc_finemap`) |
| **`global_slope`** | **β −0.068, p 0.86** | **β −0.095, p 0.93** | 0.070 (`pool_not_hc`) |
| `slope_PC1` | β −0.002, p 0.51 | β −0.034, p 0.70 | 0.055 |
| `slope_PC2` | β 0.021, p 0.37 | β −0.013, p 0.58 | 0.044 |
| `slope_PC3` | β −0.094, p 0.94 | β −0.138, p 0.98 | 0.0086 (`pool_not_hc`) |

**No enrichment of MDD high-confidence genes in any ABCD phenotype**, and the
point estimates for `global_slope` are negative in both the marginal and the
conditional model.

This null is better powered than §13's SCZ prioritised test: 213 genes against
104 gives SE 0.062 against 0.091. For `global_slope` the 95 % CI is
[−0.19, 0.05] — it excludes enrichment above about 0.05, which is below the
0.109 the SCZ locus pool showed. So for MDD specifically this is closer to an
informative null than the SCZ test was.

### Contrast with SCZ, stated carefully

SCZ's *locus pool* was enriched in `baseline_thickness` (p = 1.9e-04); MDD's
pool is not (p = 0.167). It is tempting to read that as disorder specificity.
It is not safe to: the two "pools" are different objects built by different
methods at different scales, exactly as set out above. The comparison that is
like-for-like — SCZ's 120 prioritised versus MDD's 308 high-confidence — is
**null on both sides**, for every phenotype.

## 15. Two audit checks: ancestry PCs, and the N reaching LDSC

Both raised as "these could be silently wrong"; both chased to a definite
answer rather than an argument. One was a real defect with a negligible effect,
the other was already correct.

### 15a. The ancestry PCs are the wrong ones — and it does not matter

`src/abcd/io.py:711` takes the covariate PCs from `ab_g_stc__gen_pc__01..10`,
which the release computes on the **full multi-ancestry cohort**. Restricted to
the 4,126 EUR analysis subjects they are close to constant:

| | mean | SD | mean / SD |
|---|---|---|---|
| PC1 | 0.006310 | 0.000523 | **12.1** |
| PC2 | −0.005263 | 0.001125 | 4.7 |
| PC3 | −0.001102 | 0.001230 | 0.9 |

A PC that separates ancestry groups looks exactly like this inside one group.
Recomputing PCs *within* EUR from the existing GRM (`work/eur_pca_check.sbatch`,
job 33485665) gives PCs with **25x the spread** — SD 0.0128–0.0147 against
0.000523 — so they genuinely describe within-EUR variation the release PCs do
not.

Re-running REML on `baseline_thickness` with each:

| covariates | h² | SE | n |
|---|---|---|---|
| release PCs (multi-ancestry) | 0.5754 | 0.1474 | 3,329 |
| **in-sample EUR PCs** | **0.5663** | **0.1483** | 3,336 |

**A difference of 0.009 — 0.06 of one SE.** The covariates were wrong in
principle and irrelevant in practice, and the reason is visible in the
eigenvalues: in-sample PC1 explains **0.13 %** of variance, with PC1/PC2 =
7.23/3.92. There is no meaningful structure left inside the EUR subset for a PC
to correct.

So **h² = 0.575 is not an artefact of the ancestry covariates.** Against the
published 0.20–0.40 for mean cortical thickness, three things reconcile it
without invoking stratification:

- the SE of 0.147 gives a 95 % CI of **[0.29, 0.86]**, which already overlaps
  the published range;
- the LDSC intercept is **0.990 ± 0.006** and λ_GC is 1.03 — stratification
  normally pushes the intercept above 1, not to 1;
- `baseline_thickness` is a **BLUP intercept pooling 2–4 visits**, reliability
  0.892–0.941 (§9). Published estimates use single-timepoint thickness, whose
  measurement error attenuates h². A latent-trait phenotype should return a
  higher h², and that is the expected direction.

Two corrections to earlier notes in this file: §11.2 claimed in-sample PCs were
"already sitting in `abcd_eur.eigenvec`" — **that file does not exist**;
`genetics_copy/` holds only `abcd_eur.{bed,bim,fam}` and the sparse GRM.

For the all-ancestry GRM (§16) none of this reassurance carries over: there the
structure is real and in-sample PCs are required, not optional.

### 15b. The `N = 2 x NEFFDIV2` convention did reach LDSC

Checked through the whole chain rather than in the config alone, since a wrong
N rescales h² invisibly:

| trait | munge log | N seen | expected |
|---|---|---|---|
| SCZ | `--N-col N` → `N: Sample size` | max 170,115 | N_eff ≈ 171,860 ✓ |
| MDD | `--N-col N` → `N: Sample size` | 1,183,605 | study N_eff ✓ |
| ABCD | `--N-col N` → `N: Sample size` | 4,119, constant | matches fastGWA ✓ |

All three correct.

One refinement on the risk, though: a wrong N would **not** have moved rg. LDSC's
h² scales as 1/N and gencov as 1/√(N₁N₂), so in rg = gencov / √(h²₁ h²₂) the N
terms cancel exactly. **rg is invariant to N misspecification; only h² is
exposed.** The §12 rg nulls were never at risk from this.

## 16. Cross-ancestry GRM (jobs 33485653, 33485654) — SUBMITTED

Building the all-ancestry GRM deferred at §11.1, because N is the binding
constraint on every result in §12–§14.

**Verified gain** (joined on the 8-character NDAR token, not assumed):

| | subjects |
|---|---|
| phenotyped | 8,192 |
| phenotyped ∩ EUR genotyped | **4,126** (current) |
| phenotyped ∩ all-ancestry genotyped | **7,111** |
| gain | **+2,985 (+72 %)** |

`work/grm_allanc.sbatch` builds per-chromosome GRMs (array 1–22 at `%12`, 192
CPUs, inside the 448-CPU limit) over `ABCD_chr{1..22}_hg19` — 10,072 subjects,
18,943,024 autosomal SNPs, `--maf 0.01` to match `03_gwas`. Per chromosome
rather than merging ~60 GB of filesets: `--mgrm` combines exactly, as the
SNP-count-weighted mean, and runs in parallel.

`work/grm_allanc_merge.sbatch` then produces everything the pipeline consumes —
`abcd_all`, `abcd_all_sp` (sparse, 0.05), `abcd_all.unrel` (0.05 cutoff) and
**20 in-sample ancestry PCs** — so switching over is a config edit. It refuses
to merge an incomplete set rather than quietly building from fewer chromosomes;
that quiet-wrong failure mode has already occurred twice here (§6d, §12b).

Checked before building, both favourable:

- **IDs**: the all-ancestry `.fam` is `FID=AB0000055 IID=NDAR_INVD5FWJDCY`, so
  the IID already matches `gcta_export`'s native spelling and only FID needs the
  token join `align_ids.py` already does.
- **rsIDs**: 84.7 %, against 85.4 % for `abcd_eur` — MAGMA, LDSC and PRS stay
  as feasible as they are now.

### The caveat to weigh before trusting its h²

A single GRM pooled across ancestries assumes a common allele frequency and LD
structure, which is precisely what does not hold across these groups. In-sample
PCs help but do not repair it. If the all-ancestry h² differs materially from
the EUR estimate, the defensible design is **ancestry-stratified REML
meta-analysed across groups**, not one pooled GRM.
