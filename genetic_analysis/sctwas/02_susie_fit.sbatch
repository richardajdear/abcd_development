#!/bin/bash
#SBATCH --job-name=sctwas_fit
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=1-22
#SBATCH --output=slurm/%x_%A_%a.log
# SCTWAS stage 1b, array over chromosomes: sbatch --export=ALL,CT=Ext1,NSAMP=757[,TAG=Ext1_n983] 02_susie_fit.sbatch
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development; W=$REPO/genetic_analysis/work/sctwas
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK TMPDIR=$W/tmp
$REPO/legacy/hpc/work/envs/abcdR/bin/Rscript $REPO/genetic_analysis/sctwas/02_susie_fit.R $CT $SLURM_ARRAY_TASK_ID ${NSAMP:-757} ${TAG:-$CT}
