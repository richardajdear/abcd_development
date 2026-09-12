#!/bin/bash
#SBATCH --job-name=abcd_chr22
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=/home/rajd2/rds/hpc-work/abcd_development/hpc/work/logs/%x_%j.log
#
# Smoke-test genotypes: chromosome 22 only, all 5,678 subjects.
#
# Why a chromosome rather than a subject subsample: the GWAS step's cost is
# linear in SNPs and its correctness does not depend on which ones, whereas
# dropping subjects would change the relatedness structure and so make the
# sparse-GRM path untested.  Keeping every subject means the chr22 fastGWA run
# exercises exactly the code path the genome-wide run will.
#
# No --mem: cclake caps memory at 3410 MB/CPU, so `--cpus-per-task=8` already
# implies 26.6 GB.  Asking for more than cpus x 3410 makes SLURM silently widen
# the CPU request (or reject the job) -- see README_HPC.md.
set -euo pipefail

W=/home/rajd2/rds/hpc-work/abcd_development/hpc/work
SRC=/rds/user/rajd2/hpc-work/ABCD/genetics_copy/abcd_eur

"$W/bin/plink2" \
  --bfile "$SRC" \
  --chr 22 \
  --make-bed \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --out "$W/genotype/abcd_eur_chr22"

wc -l "$W/genotype/abcd_eur_chr22.bim" "$W/genotype/abcd_eur_chr22.fam"
