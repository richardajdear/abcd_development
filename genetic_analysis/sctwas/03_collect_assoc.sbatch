#!/bin/bash
#SBATCH --job-name=sctwas_assoc
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=slurm/%x_%j.log
# SCTWAS stage 1c + 2 for one tag: gather weights, then associations for all traits x {susie, top1}.
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development; W=$REPO/genetic_analysis/work/sctwas
RS=$REPO/legacy/hpc/work/envs/abcdR/bin/Rscript; export TMPDIR=$W/tmp
$RS $REPO/genetic_analysis/sctwas/02b_collect.R $TAG
$RS $REPO/genetic_analysis/sctwas/03_assoc.R $TAG
