#!/bin/bash
#SBATCH --job-name=collect_1lmm
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=slurm/%x_%A.log
# Collect the single-LMM EUR-arm scan (rule 10: TMPDIR on rds -- run_all.sh's
# inline v2_collect uses the node's /tmp, which was full on 2026-09-18).
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
S=""; [[ "${PARC:-dsk}" == hcp ]] && S="_hcp"
ARMSUF="_eur"; [[ "${ARM:-eur}" == pooled ]] && ARMSUF=""   # ARM=pooled for the pooled arm (an EMPTY value is not distinguishable from unset via :-)
D="$REPO/genetic_analysis/work/results_70tab$S/scan_1lmm/assoc$ARMSUF"
export TMPDIR="$REPO/genetic_analysis/work/tmp/collect_$SLURM_JOB_ID"; mkdir -p "$TMPDIR"
R=/home/rajd2/rds/hpc-work/envs/genesis/bin/Rscript
for p in global_slope_1lmm baseline_thickness_1lmm; do
  n=$(ls "$D/${p}_chr"*.tsv.gz | wc -l); (( n == 22 )) || { echo "FATAL: $p has $n of 22" >&2; exit 3; }
  "$R" "$REPO/genetic_analysis/R/05_collect_assoc.R" --assoc-dir "$D" --pheno "$p" --out-dir "$D" --maf 0.01 --mac 20
done
rm -rf "$TMPDIR"; column -t "$D/gwas_summary.tsv"
