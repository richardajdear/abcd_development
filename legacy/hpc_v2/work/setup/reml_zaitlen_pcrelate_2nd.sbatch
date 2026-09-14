#!/bin/bash
#SBATCH --job-name=v2_remlpc2
#SBATCH --partition=icelake
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-5
#
# Step 05, re-run with the PC-RELATE-derived second GRM.
#
# The shipped step 05 defines its second component with GCTA --make-bK on the
# POOLED multi-ancestry GRM.  On this data that kept 3,837,371 pairs (5.6% of
# all pairs) where PC-Relate finds 29,371 relatives -- so G2 was an ancestry
# component, and the fit degenerated: global_slope aborted with two components
# constrained, slope_PC1/PC2/PC3 drove V(G1) to zero, and baseline_thickness
# ran away to V(G1)=0.74 against v1's h2 of 0.246.  Same model, same GRM
# values, same threshold semantics -- only the SET of pairs called "close"
# changes, and it changes to the estimator hpc_v2 exists to use.
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/config.sh
OUTDIR="$OUT_V2/reml_zaitlen_pcrel2nd"
mgrm="$OUTDIR/mgrm.txt"
[[ -s "$mgrm" ]] || { echo "FATAL: $mgrm missing -- run build_bk_pcrelate.sbatch" >&2; exit 3; }
require_paths PHENO COVAR_QUANT COVAR_CAT MANIFEST GCTA

names=()
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)
idx=$((SLURM_ARRAY_TASK_ID - 1))
pheno="${names[$idx]}"
mpheno=$(mpheno_for "$pheno")
echo "=== reml_zaitlen_pcrel  $pheno  (--mpheno $mpheno)  $(date) ==="
cat "$mgrm"

"$GCTA" --reml --mgrm "$mgrm" \
  --pheno "$PHENO" --mpheno "$mpheno" \
  --qcovar "$COVAR_QUANT" --covar "$COVAR_CAT" \
  --reml-est-fix --thread-num "$THREADS" \
  --out "$OUTDIR/$pheno"
echo "=== done  $(date) ==="
