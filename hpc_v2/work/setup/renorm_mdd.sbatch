#!/bin/bash
#SBATCH --job-name=v2_mdd_neff
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=1:00:00
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/work/setup/paths.sh
"$PYTHON" "$V2ROOT/work/setup/normalise_gwas.py" "$GWAS_DIR" "$FINAL/gwas" \
  "$IN/ref/maf_union.tsv" MDD_pooled,MDD_eur
for t in MDD_pooled MDD_eur; do
  awk -v t="$t" 'NR>1 && $8+0>m{m=$8+0} END{printf "%s Neff max: %d\n", t, m}' "$FINAL/gwas/$t.ma"
done
