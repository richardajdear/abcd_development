#!/bin/bash
#SBATCH --job-name=sctwas_prep
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=slurm/%x_%j.log
# SCTWAS stage 1a for one cell type: sbatch --export=ALL,CT=Ext1 01_prep_eqtl.sbatch  (run from $W)
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development; W=$REPO/genetic_analysis/work/sctwas
[[ -e $W/raw/${CT}_eqtl_full_assoc.tsv.gz.done && -e $W/raw/${CT}_eqtl_top_assoc.tsv.gz.done ]] || { echo "FATAL: $CT download incomplete" >&2; exit 2; }
export TMPDIR=$W/tmp
$REPO/legacy/hpc/work/envs/abcd/bin/python $REPO/genetic_analysis/sctwas/01_prep_eqtl.py $CT
