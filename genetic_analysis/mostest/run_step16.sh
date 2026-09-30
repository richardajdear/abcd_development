#!/bin/bash
# Step 16 end to end: submits every stage with afterok dependencies.
#   bash genetic_analysis/mostest/run_step16.sh            (from the repo root, on a login node)
#   FROM=3 bash genetic_analysis/mostest/run_step16.sh     resume at stage 3 (earlier outputs must exist)
# Stages: 1 build | 2 REGENIE step 1 (x7) | 3 REGENIE step 2 (x154) |
#         4 MOSTest (x3) + engine check | 5 clump + MAGMA (x6) | 6 collect
# afterok on an array dies if ANY task fails (rule 10): that is intended --
# a failed chromosome must not flow into MOSTest silently.  Fix and resume with FROM.
set -euo pipefail
[[ -f genetic_analysis/mostest/paths.sh ]] || { echo "run from the repo root" >&2; exit 2; }
mkdir -p slurm
D=genetic_analysis/mostest; FROM="${FROM:-1}"; dep=""
sub() {  # sub <stage> <script>
  local st="$1" sc="$2" j
  if (( st < FROM )); then return 0; fi
  j=$(sbatch --parsable ${dep:+--dependency=afterok:$dep} "$D/$sc")
  echo "stage $st  $sc  job $j${dep:+  (after $dep)}"; dep="$j"
}
sub 1 01_build.sbatch
sub 2 02_regenie_step1.sbatch
sub 3 03_regenie_step2.sbatch
sub 4 04_mostest.sbatch
sub 5 05_clump_magma.sbatch
sub 6 06_collect.sbatch
