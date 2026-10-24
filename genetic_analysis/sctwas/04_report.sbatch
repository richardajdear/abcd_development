#!/bin/bash
#SBATCH --job-name=sctwas_report
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=slurm/%x_%j.log
# SCTWAS checks: report/*.tsv + README results section (run from $W after all tags finish).
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
$REPO/legacy/hpc/work/envs/abcd/bin/python $REPO/genetic_analysis/sctwas/04_report.py
