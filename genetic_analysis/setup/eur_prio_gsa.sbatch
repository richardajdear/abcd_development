#!/bin/bash
#SBATCH --job-name=v2_eurprio
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# EUR arm, MAGMA prioritised/high-confidence gene-set tests, reusing v1's own
# prioritised_gsa.sbatch against the v2 .genes.raw files.  Seconds of compute:
# the gene analysis is already done by the MAGMA step.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/genetic_analysis/config.sh"
W="$REPO/legacy/hpc/work"

MAGMA_DIR="$OUT_V2/magma_eur" \
OUT="$OUT_V2/magma_prio_eur" \
SETS="$W/genesets/scz_prioritised.txt" TAG=sczprio COND_SET=SCZ_locus_pool \
  bash "$W/prioritised_gsa.sbatch"

MAGMA_DIR="$OUT_V2/magma_eur" \
OUT="$OUT_V2/magma_prio_eur" \
SETS="$W/genesets/mdd_highconf.txt" TAG=mddhc COND_SET=MDD_pool \
  bash "$W/prioritised_gsa.sbatch"

echo "=== done  $(date) ==="
