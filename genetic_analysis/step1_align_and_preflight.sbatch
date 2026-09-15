#!/bin/bash
#SBATCH --job-name=pheno70_align
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm/%x_%j.log
#
# The rest of README_HPC.md 3 step 1, plus step 2 (preflight).  Submitted with
# --dependency=afterok on the export job, so it runs only if the export wrote.
#
#   1. step1_verify_export.py   gate the CSD3-regenerated export against the
#                               committed local numbers.  set -e stops the
#                               chain here if it fails -- which is the point.
#   2. setup/align_export.py    .fam spelling, genotyped subset, z-score, both
#                               FID conventions, family_map.tsv, align_report.
#   3. setup/make_grm_famid.sh  FID = family_id view of v1's dense GRM.
#   4. 00_check_inputs.sh       preflight; run before ANY further submission.
#
# It deliberately stops there.  Steps 3-8 (the GENESIS scan onward) are hours
# of compute and should be launched after reading this job's log.
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p slurm
set -u

export PYTHONPATH=src:${PYTHONPATH:-}
PY="${PY:-legacy/hpc/work/envs/abcd/bin/python}"
source genetic_analysis/config.sh          # PARC -> ABCD_CONFIG_PARC, PHENO_DIR
export ABCD_CONFIG="${ABCD_CONFIG_PARC:?}"
RUN="$("$PY" -m abcd.run_dir)"
PHENO_OUT="$PHENO_DIR"
echo "=== PARC=$PARC  run=$RUN  pheno_out=$PHENO_OUT"

echo "=== 1. verify the regenerated export  $(date -u +%FT%TZ)"
"$PY" genetic_analysis/step1_verify_export.py --parc "$PARC"

echo; echo "=== 2. align the export  $(date -u +%FT%TZ)"
"$PY" genetic_analysis/setup/align_export.py "$RUN/gcta_inputs" "$PHENO_OUT"

echo; echo "=== 3. family-id view of the dense GRM  $(date -u +%FT%TZ)"
bash genetic_analysis/setup/make_grm_famid.sh

echo; echo "=== 4. preflight  $(date -u +%FT%TZ)"
bash genetic_analysis/00_check_inputs.sh

echo; echo "=== done  $(date -u +%FT%TZ)"
echo "Steps 3-8 are NOT submitted.  Read this log, then run_all.sh 03 04."
