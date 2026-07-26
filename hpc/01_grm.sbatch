#!/bin/bash
#SBATCH -A VERTES-SL2-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --job-name=abcd-grm
#SBATCH -o slurm/abcd-grm-%j.out
#SBATCH -e slurm/abcd-grm-%j.err
source "$(dirname "$0")/config.sh"

# Full GRM. MAF 0.01 and genotyping rate 0.95 are the conventional floors;
# below them GRM entries are dominated by imputation error.
"$GCTA" --bfile "$GENO" --autosome --maf 0.01 --geno 0.05 \
        --make-grm --thread-num "$SLURM_CPUS_PER_TASK" --out "$GRM"

# Unrelated subset at 0.05. ABCD is family-based by design (~2000 sibling
# pairs, ~800 MZ), so the related and unrelated analyses answer different
# questions: the unrelated GRM gives SNP-h2 free of shared-environment
# confounding, while the full GRM is what fastGWA uses as its random effect.
"$GCTA" --grm "$GRM" --grm-cutoff 0.05 --make-grm \
        --thread-num "$SLURM_CPUS_PER_TASK" --out "${GRM}_unrel"

# Sparse GRM for fastGWA (thresholded at 0.05, stored sparse).
"$GCTA" --grm "$GRM" --make-bK-sparse 0.05 \
        --thread-num "$SLURM_CPUS_PER_TASK" --out "${GRM}_sparse"

wc -l "$GRM".grm.id "${GRM}_unrel".grm.id
