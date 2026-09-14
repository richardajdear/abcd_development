# Shared configuration for the hpc_v2/ pipeline.  Sourced by every script;
# never run directly.
#
# hpc_v2 is the sibling-aware successor to hpc/: GENESIS (PC-AiR + PC-Relate)
# for a single pooled multi-ancestry GWAS that KEEPS relatives, a Zaitlen
# two-GRM REML that keeps them too, and a within-family PRS test.  hpc/ is kept
# untouched so its published results stay reproducible; nothing here overwrites
# a v1 output.
#
# DESIGN RULE (inherited from hpc/, do not regress): no path is hard-coded.
# Every one is `${VAR:-default}`; hpc_v2/config.local.sh (gitignored) is
# sourced FIRST so a root set there propagates through the derived defaults.
# Exporting ABCD_HPC_ROOT alone relocates everything.
#
# Usage:  source "$(dirname "$0")/config.sh"
#         require_paths GENO_ARRAY PHENO        # validate only what you use

set -euo pipefail

# ---------------------------------------------------------------------------
# Roots
# ---------------------------------------------------------------------------
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# Local overrides FIRST (deliberate: defaults below derive from each other).
if [[ -f "$REPO_ROOT/hpc_v2/config.local.sh" ]]; then
  # shellcheck source=/dev/null
  source "$REPO_ROOT/hpc_v2/config.local.sh"
fi

# On CSD3: RDS scratch, never /home (small quota, and it is nearly full).
ABCD_HPC_ROOT="${ABCD_HPC_ROOT:-$HOME/rds/hpc-work/ABCD}"

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
# ARRAY genotypes: ONE merged PLINK fileset (the 7.0 curated Smokescreen set,
# 11,670 x 515,228, hg19).  Used for kinship/PCA only -- KING, PC-AiR and
# PC-Relate want clean directly-genotyped SNPs, not imputed dosages.
GENO_ARRAY="${GENO_ARRAY:-$ABCD_HPC_ROOT/genotype_array/smokescreen_merged}"

# IMPUTED genotypes: per-chromosome PLINK filesets from the v1 conversion
# (TOPMed r3, rsIDs restored, R2>=0.8 & MAF-filtered at conversion).  Used for
# the association scan -- density matters there, not for kinship.
# {CHR} is substituted by geno_imp_prefix().  NOTE the two-step default:
# ${VAR:-abc{CHR}def} does NOT work -- bash ends the expansion at the brace
# closing {CHR}.  Keep placeholders out of ${...:-...} defaults.
GENO_IMP_DIR="${GENO_IMP_DIR:-$ABCD_HPC_ROOT/genotype_imputed}"
if [[ -z "${GENO_IMP_TPL:-}" ]]; then
  GENO_IMP_TPL='abcd_imp_chr{CHR}'
fi
geno_imp_prefix() {
  local chr="${1:?geno_imp_prefix needs a chromosome number}"
  local pat='{CHR}'
  printf '%s/%s\n' "$GENO_IMP_DIR" "${GENO_IMP_TPL//"$pat"/$chr}"
}

# Phenotypes/covariates written by `python -m abcd.gcta_export`.
# CRITICAL for step 06: the export's FID column carries the REAL family IDs.
# The genotype .fam files on CSD3 have FID = IID (family structure destroyed),
# so family membership must always come from $PHENO, never from a .fam.
PHENO_DIR="${PHENO_DIR:-$ABCD_HPC_ROOT/pheno}"
PHENO="${PHENO:-$PHENO_DIR/phenotypes_gcta.txt}"
COVAR_QUANT="${COVAR_QUANT:-$PHENO_DIR/covar_quant.txt}"
COVAR_CAT="${COVAR_CAT:-$PHENO_DIR/covar_categorical.txt}"
MANIFEST="${MANIFEST:-$PHENO_DIR/phenotype_manifest.tsv}"

# Dense GRM for the Zaitlen REML (step 05): the v1 pooled IMPUTED GRM over the
# FULL sample INCLUDING relatives (not any .unrel).  Point this at the GRM
# behind hpc/README_HPC.md §8.15 in config.local.sh.
GRM_FULL="${GRM_FULL:-$ABCD_HPC_ROOT/results/grm_imp/abcd_imp}"

# PRS score profiles from v1's step 06 on the imputed genotypes
# (results/prs_imp/score_<DISORDER>_<threshold>.profile).  v2 does not
# re-score; the within-family test consumes these.
PRS_PROFILE_DIR="${PRS_PROFILE_DIR:-$ABCD_HPC_ROOT/results/prs_imp}"

# EUR subject list (one ID per line, or FID IID) for the EUR-primary PRS arm.
EUR_IDS="${EUR_IDS:-$ABCD_HPC_ROOT/keep/eur_anchor.keep}"

# Optional ancestry-stratum file for heteroscedastic residuals in the null
# model (GENESIS group.var): header FID IID stratum -- v1 wrote strata_k4.tsv
# in this shape.  Empty = homoscedastic model.
STRATA_FILE="${STRATA_FILE:-}"

# ---------------------------------------------------------------------------
# Outputs (all under OUT_V2; nothing writes into v1's $OUT)
# ---------------------------------------------------------------------------
OUT_V2="${OUT_V2:-$ABCD_HPC_ROOT/results_v2}"
GDS_DIR="${GDS_DIR:-$OUT_V2/gds}"
KIN_DIR="${KIN_DIR:-$OUT_V2/kinship}"
NULL_DIR="${NULL_DIR:-$OUT_V2/nullmodel}"
ASSOC_DIR="${ASSOC_DIR:-$OUT_V2/assoc}"
REML2_DIR="${REML2_DIR:-$OUT_V2/reml_zaitlen}"
PRSFAM_DIR="${PRSFAM_DIR:-$OUT_V2/prs_family}"
LOG_DIR="${LOG_DIR:-$OUT_V2/logs}"

# ---------------------------------------------------------------------------
# Method parameters
# ---------------------------------------------------------------------------
# LD pruning ahead of KING/PC-AiR/PC-Relate (SNPRelate): MAF floor and r2.
KIN_MAF="${KIN_MAF:-0.05}"
KIN_LD_R2="${KIN_LD_R2:-0.1}"
# PC-AiR relatedness/divergence threshold: 2^(-11/2) ~ 0.022, the conventional
# 3rd-degree boundary.  Used for pcair()'s unrelated/related partition.
KIN_THRESH="${KIN_THRESH:-0.02209709}"

# Kinship threshold for making the null model's covariance matrix SPARSE.
# Separate from KIN_THRESH on purpose, because it controls a different and
# more dangerous thing.
#
# READ THIS BEFORE CHANGING IT.  pcrelateToMatrix(thresh=) does not zero
# individual entries -- it CLUSTERS samples by transitive closure (any pair
# above the threshold joins a cluster, and all within-cluster pairs are then
# kept) and zeroes only BETWEEN clusters.  So a chain of weak, noise-level
# pairs can merge the whole cohort into one block and hand fitNullModel a
# DENSE matrix, which at 11,670 subjects is ~1.1 GB and defeats the entire
# sparse design -- silently, because nothing errors.
#
# Measured on the synthetic fixture (1,500 subjects, only 4,977 pruned SNPs):
# 783 true sibling pairs, but 2,545 pairs cleared raw kinship 0.0221, and their
# transitive closure pulled 1,472 of 1,500 subjects into a single block ->
# 96 % dense.  On real data with ~100k+ pruned SNPs the noise is far smaller
# and this should not happen -- but 02_kinship.R now REPORTS the density and
# largest block, and warns above SPARSE_KIN_MAX_DENSITY, so you find out from
# the summary rather than from a memory error.
#
# NOTE the factor of 2: pcrelateToMatrix applies `thresh` AFTER scaling
# kinship by scaleKin=2, so this value is compared against 2*kinship.  The
# default 0.0442 therefore keeps pairs with raw kinship above ~0.0221 (3rd
# degree).  Raise it to 0.0884 to cluster only 2nd-degree-and-closer pairs if
# the density warning fires.
SPARSE_KIN_THRESH="${SPARSE_KIN_THRESH:-0.04419417}"
SPARSE_KIN_MAX_DENSITY="${SPARSE_KIN_MAX_DENSITY:-0.10}"
# Number of PC-AiR PCs entering the null model as fixed effects.
N_PCS="${N_PCS:-10}"
# Association-scan filters, applied by the collector (GENESIS reports
# frequency and MAC per variant; filtering at collection keeps the per-chrom
# scan outputs complete for later re-thresholding).
ASSOC_MAF="${ASSOC_MAF:-0.01}"
ASSOC_MAC="${ASSOC_MAC:-20}"
# Zaitlen second-GRM threshold: off-diagonals below this are zeroed in the bK
# GRM.  0.05 matches the v1 unrelated cutoff so the two h2 lines are readable
# against each other.
BK_THRESH="${BK_THRESH:-0.05}"

# ---------------------------------------------------------------------------
# Binaries / environments
# ---------------------------------------------------------------------------
# RSCRIPT must resolve to an R with GENESIS/SNPRelate/GWASTools/lme4 installed.
# On CSD3 build it once from hpc_v2/envs/genesis_env.yml (see README) and set
# RSCRIPT in config.local.sh.  A bare system Rscript will not have GENESIS.
RSCRIPT="${RSCRIPT:-$(command -v Rscript || echo Rscript)}"
GCTA="${GCTA:-$(command -v gcta64 || command -v gcta || echo "$HOME/bin/gcta64")}"
PLINK="${PLINK:-$(command -v plink || echo "$HOME/bin/plink")}"

# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------
SLURM_ACCOUNT="${SLURM_ACCOUNT:-VERTES-SL3-CPU}"
# icelake, EXPLICITLY: v1 lost fourteen hours to jobs landing on the cluster
# default partition (cclake) while it was drained.  Every sbatch header below
# also carries --partition; run_all.sh passes this value, which overrides.
SLURM_PARTITION="${SLURM_PARTITION:-icelake}"
THREADS="${THREADS:-8}"

# Set DRY_RUN=1 to echo commands without executing them.
DRY_RUN="${DRY_RUN:-0}"

# ---------------------------------------------------------------------------
# Helpers (same contracts as hpc/config.sh)
# ---------------------------------------------------------------------------
run() {
  if [[ "$DRY_RUN" == "1" ]]; then
    printf 'DRY_RUN: '; printf '%q ' "$@"; printf '\n'
  else
    printf '+ '; printf '%q ' "$@"; printf '\n'
    "$@"
  fi
}

# Fail early, naming the VARIABLE so the error says what to export.
# Prefix-valued variables are checked on the member file that must exist.
require_paths() {
  local missing=0 name value
  for name in "$@"; do
    value="${!name:-}"
    if [[ -z "$value" ]]; then
      echo "ERROR: $name is not set" >&2; missing=1; continue
    fi
    case "$name" in
      GENO_ARRAY)
        [[ -f "$value.bed" ]] || { echo "ERROR: $name=$value -- $value.bed not found" >&2; missing=1; } ;;
      GRM_FULL)
        if [[ ! -f "$value.grm.bin" ]]; then
          echo "ERROR: $name=$value -- $value.grm.bin not found (dense GRM, full sample)" >&2; missing=1
        elif [[ ! -f "$value.grm.id" ]]; then
          echo "ERROR: $name=$value -- $value.grm.id missing (GRM incomplete)" >&2; missing=1
        fi ;;
      RSCRIPT|GCTA|PLINK)
        command -v "$value" >/dev/null 2>&1 || [[ -x "$value" ]] \
          || { echo "ERROR: $name=$value is not executable" >&2; missing=1; } ;;
      *)
        [[ -e "$value" ]] || { echo "ERROR: $name=$value does not exist" >&2; missing=1; } ;;
    esac
  done
  [[ "$missing" == "0" ]] || {
    echo >&2
    echo "Set the variables above, or write $REPO_ROOT/hpc_v2/config.local.sh." >&2
    exit 2
  }
}

# Phenotype names / --mpheno indices from the manifest, never hand-maintained.
phenotype_names() {
  require_paths MANIFEST
  awk 'NR>1 {print $1}' "$MANIFEST"
}
mpheno_for() {
  local want="$1"
  awk -v w="$want" 'NR>1 && $1==w {print $2; found=1} END{if(!found) exit 3}' "$MANIFEST"
}

# Directory creation is explicit, never a side effect of sourcing.
ensure_dirs() {
  local d
  for d in "$@"; do mkdir -p "$d"; done
}

# R scripts live next to this file; resolve from BASH_SOURCE (real even under
# sbatch, where $0 is a spool copy).
RSCRIPT_DIR="${RSCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/R}"
