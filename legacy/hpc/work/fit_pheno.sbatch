#!/bin/bash
#SBATCH --job-name=abcd_fit
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:30:00
#SBATCH --output=/home/rajd2/rds/hpc-work/abcd_development/hpc/work/logs/%x_%j.log
#
# The three local stages that produce the GCTA inputs, run on a compute node
# because the login node kills long jobs:
#
#   R/fit_lmm.R          68 per-region lme4 fits  (the expensive one)
#   abcd.phenotype       BLUPs, reliability, slope PCA
#   abcd.gcta_export     phenotypes_gcta.txt + covariates + manifest
#   work/align_ids.py    FID/IID -> the .fam's spelling (see that file)
#
# `python -m abcd.assemble` is NOT here: it runs in seconds on the login node
# and its output is shared across every config of the same release.
set -euo pipefail

REPO=/home/rajd2/rds/hpc-work/abcd_development
W=$REPO/hpc/work
export ABCD_ROOT=$W/abcd_root
export PYTHONPATH=$REPO/src
export ABCD_CONFIG=${ABCD_CONFIG:-ct_70_noglobal_mv2_genetic}
PY=$W/envs/abcd/bin/python

# R/fit_lmm.R resolves the run directory by shelling out to `python -m
# abcd.run_dir` with a BARE `python` -- deliberately, so the config hash has a
# single implementation.  On CSD3 a bare `python` is the system 3.7.4, which
# cannot even import the package (`from typing import Literal` needs 3.8), so
# the R step died with "could not resolve a run directory".  Putting our env
# first on PATH is the fix; `make help` fails identically for the same reason.
export PATH="$W/envs/abcd/bin:$PATH"

cd "$REPO"
RUN_DIR=$("$PY" -m abcd.run_dir)
echo "run dir: $RUN_DIR"

echo "=== fit_lmm ($(date)) ==="
"$W/bin/Rscript" R/fit_lmm.R --run-dir "$RUN_DIR" --cores "${SLURM_CPUS_PER_TASK:-16}"

echo "=== phenotype ($(date)) ==="
"$PY" -m abcd.phenotype

echo "=== gcta_export ($(date)) ==="
"$PY" -m abcd.gcta_export --out-dir "$W/pheno"

echo "=== align_ids ($(date)) ==="
"$PY" "$W/align_ids.py" "$W/pheno" /rds/user/rajd2/hpc-work/ABCD/genetics_copy/abcd_eur

echo "=== done ($(date)) ==="
