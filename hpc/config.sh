# Shared configuration for the hpc/ pipeline.  Sourced by every script; never
# run directly.
#
# DESIGN RULE: no path is hard-coded.  Every one is `${VAR:-default}`, so the
# whole pipeline relocates by exporting variables -- which is what lets the same
# scripts run on CSD3 and against a synthetic fixture on a laptop.  To point it
# somewhere new, either export the variables or write hpc/config.local.sh
# (gitignored, sourced first so that the defaults below derive from it; see the
# note at that line).  Exporting ABCD_HPC_ROOT alone relocates everything.
#
# Usage:  source "$(dirname "$0")/config.sh"
#         require_paths GENO_BED PHENO        # validate only what you use

set -euo pipefail

# ---------------------------------------------------------------------------
# Roots
# ---------------------------------------------------------------------------
# Repository containing this script (two levels up from hpc/config.sh).
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# Local overrides are sourced FIRST, not last, and that ordering is deliberate.
# Everything below is `${VAR:-default}` and the defaults *derive* from each
# other (GRM_DIR from OUT, OUT from ABCD_HPC_ROOT).  Sourcing overrides first
# means setting ABCD_HPC_ROOT alone relocates the entire tree; sourcing them
# last would require overriding every derived path individually, and the ones
# you forgot would keep pointing at the CSD3 defaults.
if [[ -f "$REPO_ROOT/hpc/config.local.sh" ]]; then
  # shellcheck source=/dev/null
  source "$REPO_ROOT/hpc/config.local.sh"
fi

# Where genotypes, GRMs and results live.  On CSD3 this is RDS scratch, not
# /home: GRMs for ~8k subjects are tens of GB and /home has a small quota.
ABCD_HPC_ROOT="${ABCD_HPC_ROOT:-$HOME/rds/hpc-work/ABCD}"

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
# PLINK binary fileset prefix (expects $GENO.bed/.bim/.fam).
GENO="${GENO:-$ABCD_HPC_ROOT/genotype/ABCD_release_7.0_QCed}"

# ALL-ANCESTRY genotypes, supplied per-chromosome rather than as one fileset.
# $GENO above is the single-ancestry (EUR) merged fileset every step defaults
# to; these two build the cross-ancestry GRM that raises N.
#
# GENO_ALLANC_DIR holds the per-chromosome filesets and GENO_ALLANC_TPL is the
# basename with {CHR} standing in for the chromosome number, because the naming
# differs between releases -- 4.0 shipped ABCD_chr1_hg19, and a 7.0 tree will
# not necessarily match.  Keeping the template a variable means a release
# change is one line in config.local.sh rather than an edit to the sbatch.
GENO_ALLANC_DIR="${GENO_ALLANC_DIR:-$ABCD_HPC_ROOT/genotype_allanc}"
# NOTE the two-step default.  Writing ${GENO_ALLANC_TPL:-ABCD_chr{CHR}_hg19}
# does NOT work: bash ends the parameter expansion at the first unquoted '}',
# which is the one closing {CHR}, so '_hg19}' is appended as literal text and
# an overridden template silently gains a '_hg19}' suffix.  Caught by the
# expansion test below; keep the placeholder out of ${...:-...} defaults.
if [[ -z "${GENO_ALLANC_TPL:-}" ]]; then
  GENO_ALLANC_TPL='ABCD_chr{CHR}_hg19'
fi

# Expand the template for one chromosome: geno_allanc_prefix 7 -> /path/ABCD_chr7_hg19
geno_allanc_prefix() {
  local chr="${1:?geno_allanc_prefix needs a chromosome number}"
  # Brace-literal substitution: ${var//\{CHR\}/...} escapes the OPENING brace
  # only, so the closing one is left in the output as a literal '}'. Assigning
  # the pattern to a variable first sidesteps the quoting entirely.
  local pat='{CHR}'
  printf '%s/%s\n' "$GENO_ALLANC_DIR" "${GENO_ALLANC_TPL//"$pat"/$chr}"
}

# Directory written by `python -m abcd.gcta_export`.
PHENO_DIR="${PHENO_DIR:-$ABCD_HPC_ROOT/pheno}"
PHENO="${PHENO:-$PHENO_DIR/phenotypes_gcta.txt}"
COVAR_QUANT="${COVAR_QUANT:-$PHENO_DIR/covar_quant.txt}"
COVAR_CAT="${COVAR_CAT:-$PHENO_DIR/covar_categorical.txt}"
MANIFEST="${MANIFEST:-$PHENO_DIR/phenotype_manifest.tsv}"

# Published GWAS summary statistics downloaded from the PGC (not ABCD data).
SUMSTATS_DIR="${SUMSTATS_DIR:-$ABCD_HPC_ROOT/sumstats}"

# MAGMA reference data: 1000G European LD panel and SNP-gene annotation.
MAGMA_REF_DIR="${MAGMA_REF_DIR:-$ABCD_HPC_ROOT/magma_ref}"
MAGMA_REF="${MAGMA_REF:-$MAGMA_REF_DIR/g1000_eur}"
MAGMA_GENE_LOC="${MAGMA_GENE_LOC:-$MAGMA_REF_DIR/NCBI37.3.gene.loc}"
# SNP->gene annotation produced by `magma --annotate window=35,10`; 04_magma
# builds it if absent, so this is an output path as much as an input one.
MAGMA_ANNOT="${MAGMA_ANNOT:-$MAGMA_REF_DIR/ncbi37.window35-10.genes.annot}"

# External disorder GWAS summary statistics for the reverse-direction test.
# Variables rather than literals inside 04_magma: these are the files most
# likely to be renamed or re-downloaded, and a missing one should be reported
# by name at validation rather than as a MAGMA parse error 40 minutes in.
# Each needs its SNP-id and p-value column names, which differ between consortia.
SCZ_SUMSTATS="${SCZ_SUMSTATS:-$SUMSTATS_DIR/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.fixed.tsv}"
SCZ_SNP_COL="${SCZ_SNP_COL:-ID}"
SCZ_P_COL="${SCZ_P_COL:-PVAL}"
# Effect allele and effect size, for PRS weights (06) and LDSC signing (05).
# The signed-sumstats spec is "COLUMN,null-value": log-odds are null 0, odds
# ratios null 1 -- getting this wrong silently flips the sign of rg.
SCZ_A1_COL="${SCZ_A1_COL:-A1}"
SCZ_EFFECT_COL="${SCZ_EFFECT_COL:-BETA}"
SCZ_SIGNED="${SCZ_SIGNED:-BETA,0}"
SCZ_N_SPEC="${SCZ_N_SPEC:---N-col NEFF}"
MDD_SUMSTATS="${MDD_SUMSTATS:-$SUMSTATS_DIR/pgc-mdd2025_no23andMe_div_v3-49-46-01_formatted.tsv}"
MDD_SNP_COL="${MDD_SNP_COL:-rsid}"
MDD_P_COL="${MDD_P_COL:-p_value}"
MDD_A1_COL="${MDD_A1_COL:-effect_allele}"
MDD_EFFECT_COL="${MDD_EFFECT_COL:-beta}"
MDD_SIGNED="${MDD_SIGNED:-beta,0}"
MDD_N_SPEC="${MDD_N_SPEC:---N-col N}"

# ---------------------------------------------------------------------------
# LDSC (step 05) reference data and binaries.
# ---------------------------------------------------------------------------
# @ is LDSC's own per-chromosome placeholder, not a shell construct: LDSC
# expands eur_w_ld_chr/@ to .../1.l2.ldscore.gz .. 22.l2.ldscore.gz.  Quote it.
LDSC_REF_DIR="${LDSC_REF_DIR:-$ABCD_HPC_ROOT/ldsc_ref}"
LD_REF="${LD_REF:-$LDSC_REF_DIR/eur_w_ld_chr/}"
LD_WEIGHTS="${LD_WEIGHTS:-$LDSC_REF_DIR/eur_w_ld_chr/}"
HM3_SNPLIST="${HM3_SNPLIST:-$LDSC_REF_DIR/w_hm3.snplist}"
LDSC="${LDSC:-$(command -v ldsc.py || echo "$HOME/bin/ldsc.py")}"
LDSC_MUNGE="${LDSC_MUNGE:-$(command -v munge_sumstats.py || echo "$HOME/bin/munge_sumstats.py")}"

# ---------------------------------------------------------------------------
# PRS (step 06): clumping and p-value thresholds.
# ---------------------------------------------------------------------------
# Conventional P+T settings: r2<0.1 within 250 kb.  CLUMP_P1 is 1 so that
# clumping defines independent loci across the whole p-value range and the
# thresholds below do the selecting -- setting it lower would silently cap
# every score at that threshold.
CLUMP_P1="${CLUMP_P1:-1}"
CLUMP_P2="${CLUMP_P2:-1}"
CLUMP_R2="${CLUMP_R2:-0.1}"
CLUMP_KB="${CLUMP_KB:-250}"
PRS_THRESHOLDS="${PRS_THRESHOLDS:-5e-8 1e-5 0.001 0.01 0.05 0.1 0.5 1}"
PLINK="${PLINK:-$(command -v plink || echo "$HOME/bin/plink")}"
RSCRIPT="${RSCRIPT:-$(command -v Rscript || echo Rscript)}"

# Gene-set / gene-covariate files written by `python -m abcd.magma_export`.
GENESET_DIR="${GENESET_DIR:-$ABCD_HPC_ROOT/genesets}"

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
OUT="${OUT:-$ABCD_HPC_ROOT/results}"
GRM_DIR="${GRM_DIR:-$OUT/grm}"
GRM="${GRM:-$GRM_DIR/abcd_full}"
# Relatedness-pruned GRM used by 02_reml.  A variable rather than "$GRM.unrel"
# spelled out in each script, so require_paths can validate it by name and a
# missing pruned GRM is reported as such instead of as a GCTA read error.
GRM_UNREL="${GRM_UNREL:-$GRM.unrel}"
GRM_SPARSE="${GRM_SPARSE:-$GRM_DIR/abcd_sparse}"

# All-ancestry GRM tree, written by work/grm_allanc*.sbatch.  Separate from the
# EUR paths above so both can coexist: the switch-over is done by pointing GRM,
# GRM_UNREL and GRM_SPARSE at these in config.local.sh, which keeps the EUR
# results reproducible rather than overwriting them.
GRM_ALLANC_DIR="${GRM_ALLANC_DIR:-$OUT/grm_allanc}"
GRM_ALLANC="${GRM_ALLANC:-$GRM_ALLANC_DIR/abcd_all}"
GRM_ALLANC_UNREL="${GRM_ALLANC_UNREL:-$GRM_ALLANC.unrel}"
GRM_ALLANC_SPARSE="${GRM_ALLANC_SPARSE:-$GRM_ALLANC_DIR/abcd_all_sp}"
GRM_ALLANC_PCA="${GRM_ALLANC_PCA:-$GRM_ALLANC_DIR/abcd_all_pca}"

REML_DIR="${REML_DIR:-$OUT/reml}"
GWAS_DIR="${GWAS_DIR:-$OUT/gwas}"
MAGMA_DIR="${MAGMA_DIR:-$OUT/magma}"
LDSC_DIR="${LDSC_DIR:-$OUT/ldsc}"
PRS_DIR="${PRS_DIR:-$OUT/prs}"
LOG_DIR="${LOG_DIR:-$OUT/logs}"

# ---------------------------------------------------------------------------
# Binaries
# ---------------------------------------------------------------------------
# Paths, not `module load`.  GCTA and MAGMA are user-installed on CSD3, and a
# silently different GCTA version changes REML convergence behaviour -- so the
# version in use should be visible in the config, not implied by the module
# environment.  `command -v` fallback keeps it working where they are on PATH.
GCTA="${GCTA:-$(command -v gcta64 || command -v gcta || echo "$HOME/bin/gcta64")}"
PLINK2="${PLINK2:-$(command -v plink2 || echo "$HOME/bin/plink2")}"
MAGMA="${MAGMA:-$(command -v magma || echo "$HOME/bin/magma")}"

# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------
SLURM_ACCOUNT="${SLURM_ACCOUNT:-VERTES-SL2-CPU}"
SLURM_PARTITION="${SLURM_PARTITION:-cclake}"
THREADS="${THREADS:-16}"

# GRM relatedness cutoff for the REML unrelated subset.  0.05 is conventional
# for SNP-h2: it removes the twin/sibling pairs whose shared environment would
# otherwise be read as additive genetic variance.  ABCD is a family study, so
# this discards a large fraction of the sample -- the effective N for REML is
# far below the exported 8,192 and the resulting SE is correspondingly wide.
GRM_CUTOFF="${GRM_CUTOFF:-0.05}"

# fastGWA sparse-GRM threshold (GCTA's own default).
SPARSE_CUTOFF="${SPARSE_CUTOFF:-0.05}"

# Minimum allele frequency for the GRM and the association scan.
#
# This is the only SNP filter here, because it is the only one GCTA implements:
# --geno (missingness) and --hwe are PLINK options and GCTA rejects them
# outright.  Do that filtering upstream if GENO is not already QC'd; see the
# note in 01_grm.sbatch.
MAF="${MAF:-0.01}"

# Set DRY_RUN=1 to echo commands without executing them.
DRY_RUN="${DRY_RUN:-0}"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Echo and run, or just echo under DRY_RUN.
run() {
  if [[ "$DRY_RUN" == "1" ]]; then
    printf 'DRY_RUN: '; printf '%q ' "$@"; printf '\n'
  else
    printf '+ '; printf '%q ' "$@"; printf '\n'
    "$@"
  fi
}

# Fail early with a message naming the variable, not just the missing path.
# Called with variable NAMES so the error can tell you what to export.
require_paths() {
  local missing=0 name value
  for name in "$@"; do
    value="${!name:-}"
    if [[ -z "$value" ]]; then
      echo "ERROR: $name is not set" >&2; missing=1; continue
    fi
    # Several of these variables are file-set PREFIXES, not files.  Checking
    # -e on the prefix always fails for them, which made 02/03 abort even with a
    # valid GRM present.  Check the member file that must exist.
    if [[ "$name" == "GENO" || "$name" == "MAGMA_REF" ]]; then
      [[ -f "$value.bed" ]] || { echo "ERROR: $name=$value -- $value.bed not found" >&2; missing=1; }
    elif [[ "$name" == GRM_SPARSE ]]; then
      # A sparse GRM is NOT a dense one with fewer entries: --make-bK-sparse
      # writes .grm.sp (a three-column list of retained pairs) and .grm.id, and
      # no .grm.bin at all.  Checking for .grm.bin here rejected every valid
      # sparse GRM, so 03_gwas aborted in require_paths before running.
      if [[ ! -f "$value.grm.sp" ]]; then
        echo "ERROR: $name=$value -- $value.grm.sp not found (--make-bK-sparse output)" >&2
        missing=1
      elif [[ ! -f "$value.grm.id" ]]; then
        echo "ERROR: $name=$value -- $value.grm.id missing (GRM is incomplete)" >&2
        missing=1
      fi
    elif [[ "$name" == GRM || "$name" == GRM_UNREL ]]; then
      # GCTA writes .grm.bin/.grm.id (binary) or .grm.gz/.grm.id (list).
      if [[ ! -f "$value.grm.bin" && ! -f "$value.grm.gz" ]]; then
        echo "ERROR: $name=$value -- neither $value.grm.bin nor $value.grm.gz found" >&2
        missing=1
      elif [[ ! -f "$value.grm.id" ]]; then
        echo "ERROR: $name=$value -- $value.grm.id missing (GRM is incomplete)" >&2
        missing=1
      fi
    elif [[ "$name" == GCTA || "$name" == PLINK || "$name" == PLINK2 \
         || "$name" == MAGMA || "$name" == LDSC || "$name" == LDSC_MUNGE ]]; then
      [[ -x "$value" ]] || { echo "ERROR: $name=$value is not executable" >&2; missing=1; }
    else
      [[ -e "$value" ]] || { echo "ERROR: $name=$value does not exist" >&2; missing=1; }
    fi
  done
  [[ "$missing" == "0" ]] || {
    echo >&2
    echo "Set the variables above, or write $REPO_ROOT/hpc/config.local.sh." >&2
    exit 2
  }
}

# Phenotype names in priority order, read from the manifest the export wrote.
# Never hand-maintained here: a list that disagrees with the manifest runs a
# GWAS under the wrong phenotype's name.
phenotype_names() {
  require_paths MANIFEST
  awk 'NR>1 {print $1}' "$MANIFEST"
}

# --mpheno index for one phenotype, from the same manifest.
mpheno_for() {
  local want="$1"
  awk -v w="$want" 'NR>1 && $1==w {print $2; found=1} END{if(!found) exit 3}' "$MANIFEST"
}

# Directory creation is a function, not a side effect of sourcing.  Sourcing
# this file must be harmless: with no overrides the defaults point at CSD3's
# $HOME/rds, and an eager `mkdir -p` there fails on any other machine and (under
# `set -e`) aborts before the script can report what was actually misconfigured.
# Scripts that write call ensure_dirs with the directories they need.
ensure_dirs() {
  local d
  for d in "$@"; do mkdir -p "$d"; done
}
