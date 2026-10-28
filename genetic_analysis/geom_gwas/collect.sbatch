#!/bin/bash
#SBATCH --job-name=geom_collect
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
# Final collection of a scan_geom arm (as collect_1lmm.sbatch): refuses a
# phenotype without all 22 chromosomes; TMPDIR on rds (rule 10).
#   sbatch collect.sbatch <assoc_dir> <pheno> [<pheno> ...]
set -euo pipefail
D="$1"; shift
REPO=/home/rajd2/rds/hpc-work/abcd_development
export TMPDIR="$REPO/genetic_analysis/work/tmp/geom_collect_${SLURM_JOB_ID:-$$}"; mkdir -p "$TMPDIR"
R=/home/rajd2/rds/hpc-work/envs/genesis/bin/Rscript
for p in "$@"; do
  n=$(ls "$D/${p}_chr"*.tsv.gz | wc -l); (( n == 22 )) || { echo "FATAL: $p has $n of 22" >&2; exit 3; }
  "$R" "$REPO/genetic_analysis/R/05_collect_assoc.R" --assoc-dir "$D" --pheno "$p" --out-dir "$D" --maf 0.01 --mac 20
done
rm -rf "$TMPDIR"; column -t "$D/gwas_summary.tsv"
