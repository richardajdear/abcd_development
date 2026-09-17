#!/bin/bash
#SBATCH --job-name=scz25_collect
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=slurm/%x_%A.log
set -eo pipefail
cd /home/rajd2/rds/hpc-work/abcd_development
legacy/hpc/work/envs/abcd/bin/python genetic_analysis/step9_scz2025_collect.py
