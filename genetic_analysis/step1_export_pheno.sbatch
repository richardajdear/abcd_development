#!/bin/bash
#SBATCH --job-name=pheno70_export   # PARC=hcp sbatch ... for HCP-MMP
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=slurm/%x_%j.log
#
# Step 1 of README_HPC.md 3: produce the 7.0-tabulation phenotype export on
# CSD3 rather than rsyncing it from the laptop.
#
# The README's Step 0 assumes the export travels by rsync from the laptop.  The
# true 7.0 release tables are already on CSD3 (verified: 7,607 six-year
# thickness rows vs 4,086 in 6.0), so the export is reproducible here, and
# reproducing it is checkable against the committed local numbers in
# docs/vintage_comparison.csv (8,716 subjects, 26,949 scans in the model, 3,005
# four-visit children, age centre 12.7973, median slope reliability 0.211).
# step1_verify_export.py checks exactly those.
#
# This is the settled specification only (ct_70_noglobal_mv2_genetic); the
# other configs in tools/rerun_local.sh are local-report inputs, not genetics
# inputs, and are deliberately not re-run here.
#
#   sbatch genetic_analysis/step1_export_pheno.sbatch
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p slurm
set -u

export PYTHONPATH=src:${PYTHONPATH:-}
# Parcellation switch (config.local.sh): PARC=dsk (default) or PARC=hcp selects
# the local-pipeline config and, downstream, PHENO_DIR/OUT_V2.
source genetic_analysis/config.sh
export ABCD_CONFIG="${ABCD_CONFIG_PARC:?config.local.sh must set ABCD_CONFIG_PARC}"
echo "=== PARC=$PARC  ABCD_CONFIG=$ABCD_CONFIG"

# PY must be EXPORTED and absolute: R/fit_lmm.R resolves the run directory by
# calling `python -m abcd.run_dir` through Sys.getenv("PY", unset = "python"),
# so an unexported PY sends it to CSD3's system python 3.7, which has no
# typing.Literal and fails on `import abcd` (job 35547946).
export PY="${PY:-$PWD/legacy/hpc/work/envs/abcd/bin/python}"
# The lme4 fit must use the abcdR env (R 4.5.3), NOT $RSCRIPT from config.sh,
# which is the GENESIS R (4.4.3): sourcing config.sh above already set RSCRIPT,
# so a ${RSCRIPT:-...} fallback silently kept the wrong interpreter (job
# 35581499).  Both parcellations must be fitted with the same lme4 or the
# DK-vs-HCP comparison confounds atlas with software.
export RSCRIPT="${LMER_RSCRIPT:-$PWD/legacy/hpc/work/envs/abcdR/bin/Rscript}"
CORES="${SLURM_CPUS_PER_TASK:-8}"

echo "=== config $ABCD_CONFIG  $(date -u +%FT%TZ)"
"$PY" -V; "$RSCRIPT" -e 'cat(R.version.string, "\n")'

echo "=== assemble  $(date -u +%FT%TZ)"
"$PY" -m abcd.assemble

# The run directory is <hash>-named; resolve it from the config now that
# assemble has created it (run_dir refuses to resolve a directory that does not
# exist).  Passing it explicitly means fit_lmm.R never has to re-derive it.
RUN_DIR="$("$PY" -m abcd.run_dir)"
echo "run dir: $RUN_DIR"
echo "=== fit_lmm (${CORES} cores)  $(date -u +%FT%TZ)"
"$RSCRIPT" R/fit_lmm.R --run-dir "$RUN_DIR" --cores "$CORES"

echo "=== phenotype  $(date -u +%FT%TZ)"
"$PY" -m abcd.phenotype

echo "=== gcta_export  $(date -u +%FT%TZ)"
"$PY" -m abcd.gcta_export

echo "=== done  $(date -u +%FT%TZ)"
