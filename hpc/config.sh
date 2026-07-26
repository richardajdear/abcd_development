#!/bin/bash
# Shared configuration. Every script sources this; override by editing here
# or by exporting the variables before submitting.
set -e -o pipefail

: "${ABCD_HPC_ROOT:=/home/rajd2/rds/hpc-work/ABCD}"
: "${SLURM_ACCOUNT:=VERTES-SL2-CPU}"
: "${SLURM_PARTITION:=cclake}"

GENO="$ABCD_HPC_ROOT/genotypes/abcd_imputed"     # PLINK prefix (bed/bim/fam)
PHENO_DIR="$ABCD_HPC_ROOT/gcta_inputs"
OUT="$ABCD_HPC_ROOT/out"
GRM="$OUT/grm/abcd"
GWAS_DIR="$ABCD_HPC_ROOT/gwas"
MAGMA_REF="$ABCD_HPC_ROOT/magma"

# Binaries. On CSD3 these were user-installed rather than modules; adjust if
# that has changed. Deliberately not `module load`-ing blindly: a silently
# different GCTA version changes REML convergence behaviour.
GCTA="${GCTA:-$HOME/bin/gcta64}"
MAGMA="${MAGMA:-$HOME/bin/magma}"
PLINK="${PLINK:-$HOME/bin/plink2}"

# Phenotypes to analyse. These are the column names in the phenotype file
# produced by `python -m abcd.gcta_export` (see hpc/README.md).
PHENOTYPES=(intercept slope)

mkdir -p "$OUT"/{grm,reml,gwas,magma} slurm
