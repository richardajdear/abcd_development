#!/bin/bash
#SBATCH --job-name=v2_env
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
# Build the GENESIS R environment for hpc_v2 (README_HPC.md section 4a).
# Runs under SLURM, not on the login node: the solve+download is long and
# /home has no room, so both the env and the package cache live on /rds.
set -euo pipefail

REPO=/home/rajd2/rds/hpc-work/abcd_development
MM="$REPO/hpc/work/bin/micromamba"
export MAMBA_ROOT_PREFIX="$REPO/hpc/work/mamba"
PREFIX=/rds/user/rajd2/hpc-work/envs/genesis

echo "=== micromamba $($MM --version) ; prefix=$PREFIX ; $(date) ==="
"$MM" create -y -p "$PREFIX" -f "$REPO/hpc_v2/envs/genesis_env.yml" \
  --channel-priority flexible
echo "=== solve done $(date) ==="
"$PREFIX/bin/Rscript" -e 'suppressMessages({library(GENESIS);library(SNPRelate);library(GWASTools);library(gdsfmt);library(lme4);library(lmerTest);library(data.table);library(optparse)}); cat("ALL PACKAGES LOAD OK\n"); print(R.version.string); print(packageVersion("GENESIS"))'
echo "=== finished $(date) ==="
