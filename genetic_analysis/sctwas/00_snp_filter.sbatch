#!/bin/bash
#SBATCH --job-name=sctwas_snpfilter
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --output=slurm/%x_%j.log
# SCTWAS stage 0: SNP filter + GWAS alignment (run from $W).
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
export TMPDIR=$REPO/genetic_analysis/work/sctwas/tmp
$REPO/legacy/hpc/work/envs/abcd/bin/python $REPO/genetic_analysis/sctwas/00_snp_filter.py
