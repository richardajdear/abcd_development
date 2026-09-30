#!/bin/bash
# D7 on CSD3, from the repo root:  bash directions/d7_partitioned_prs/run_d7.sh
# Submits step 1-2 (scores), then step 3 (15-task array) afterok, then step 5.
# afterok on an array dies if any task fails (rule 10): check the logs and
# resubmit single tasks with  sbatch --array=<i> directions/d7_partitioned_prs/d7_partition.sbatch
set -euo pipefail
D=directions/d7_partitioned_prs
source "$D/config.sh"
"$D7PY" -m pytest -q "$D/code/test_partition_core.py"
[[ -s "$SYMPTOMS" ]] || echo "WARN: $SYMPTOMS missing: symptom tasks will fail until laptop step L1 is copied"
N=$(( ${#ARMS[@]} * 3 ))
j1=$(sbatch --parsable "$D/d7_scores.sbatch")
j2=$(sbatch --parsable --dependency=afterok:$j1 --array=1-$N "$D/d7_partition.sbatch")
j3=$(sbatch --parsable --dependency=afterok:$j2 --job-name=d7_collect --account=vertes-sl3-cpu \
     --partition=icelake --time=00:20:00 --mem=8G --output="$WORK/logs/%x_%A.log" \
     --wrap "cd $REPO && $D7PY $D/code/05_collect.py")
echo "scores $j1 -> partition $j2 (1-$N) -> collect $j3"
