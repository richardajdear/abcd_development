#!/bin/bash
#SBATCH --job-name=geom_py
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
# Run a geom_gwas python helper off the login node:  sbatch py.sbatch <script.py> [args...]
set -euo pipefail
exec /home/rajd2/rds/hpc-work/abcd_development/legacy/hpc/work/envs/abcd/bin/python "$@"
