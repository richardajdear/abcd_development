#!/bin/bash
#SBATCH --job-name=v3_prio
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#
# MAGMA prioritised / high-confidence gene-set tests for the four Experiment A
# phenotypes -- the SCZ_locus_pool and MDD_pool panel the slide draws.
#
# 02_ldsc_magma.sbatch produced the .genes.raw files in magma_eur_v3/ but not
# these set tests, which live in a separate harness
# (hpc/work/prioritised_gsa.sbatch, invoked for the settled five by
# hpc_v2/work/setup/eur_prio_gsa.sbatch).  Seconds of compute: the gene
# analysis is already done, only the set step runs.
#
# ONE PIPELINE MODIFICATION, documented per hpc_v3/README.md section 2:
# prioritised_gsa.sbatch had its phenotype list hard-coded.  It is now
# NAMES=${NAMES:-<the original list>}, so this script can pass the four new
# phenotypes and every existing caller keeps its previous behaviour verbatim.
#
# Fresh OUT directory for the same reason as 02_ldsc_magma.sbatch: the
# collector globs $OUT/*_$TAG.gsa.out and rewrites <TAG>_gsa_summary.tsv from
# whatever it finds, so pointing it at magma_prio_eur/ would rewrite a
# published v2 table.  make_results_v3.py reads both summaries.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"

V3_NAMES="slope_topDelta slope_topC3 slope_projDelta slope_projC3"
for n in $V3_NAMES; do
  [[ -f "$OUT_V2/magma_eur_v3/$n.genes.raw" ]] \
    || { echo "FATAL: $OUT_V2/magma_eur_v3/$n.genes.raw missing" >&2; exit 3; }
done

NAMES="$V3_NAMES" MAGMA_DIR="$OUT_V2/magma_eur_v3" \
OUT="$OUT_V2/magma_prio_eur_v3" \
SETS="$W/genesets/scz_prioritised.txt" TAG=sczprio COND_SET=SCZ_locus_pool \
  bash "$W/prioritised_gsa.sbatch"

NAMES="$V3_NAMES" MAGMA_DIR="$OUT_V2/magma_eur_v3" \
OUT="$OUT_V2/magma_prio_eur_v3" \
SETS="$W/genesets/mdd_highconf.txt" TAG=mddhc COND_SET=MDD_pool \
  bash "$W/prioritised_gsa.sbatch"

for f in sczprio mddhc; do
  s="$OUT_V2/magma_prio_eur_v3/${f}_gsa_summary.tsv"
  [[ -s "$s" ]] || { echo "FATAL: $s missing/empty" >&2; exit 4; }
  echo "  $(basename "$s"): $(( $(wc -l < "$s") - 1 )) rows"
done
echo "=== done  $(date) ==="
