#!/bin/bash
#SBATCH --job-name=v3_paired
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#
# Experiment A -- the paired Delta-beta, the min-p permutation and the
# random-region null.  See hpc_v3/prs_paired_delta.py for the statistics.
#
# Primary arm reproduces prs_assoc.R's model exactly (which omits age -- see
# the module docstring); the second arm adds baseline_age and n_visits as a
# sensitivity check.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
PY="$REPO/hpc/work/envs/abcd/bin/python"
OUTD="$REPO/hpc_v3/work/prs"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "=== primary (prs_assoc.R's covariate set)  $(date) ==="
"$PY" "$REPO/hpc_v3/prs_paired_delta.py" --out-dir "$OUTD"

echo "=== sensitivity (+ baseline_age, n_visits)  $(date) ==="
"$PY" "$REPO/hpc_v3/prs_paired_delta.py" --out-dir "$OUTD" --with-age

for f in prs_paired_delta prs_minp_permutation prs_random_region_null; do
  for t in "" _withage; do
    [[ -s "$OUTD/$f$t.tsv" ]] || { echo "FATAL: $OUTD/$f$t.tsv missing" >&2; exit 4; }
    echo "  $f$t.tsv: $(( $(wc -l < "$OUTD/$f$t.tsv") - 1 )) rows"
  done
done
echo "=== done  $(date) ==="
