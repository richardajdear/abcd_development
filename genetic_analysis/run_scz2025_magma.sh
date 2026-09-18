#!/bin/bash
# Submit the step-10 MAGMA chain: prep -> genes (array 1-5) -> finish.
#   bash genetic_analysis/run_scz2025_magma.sh
#   START=finish bash genetic_analysis/run_scz2025_magma.sh   # resume
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development; cd "$REPO"; mkdir -p slurm
G=genetic_analysis; START="${START:-prep}"
sub() { sbatch --parsable "$@" | cut -d';' -f1; }
J0=""; J1=""
case "$START" in
  prep)   J0=$(sub $G/step10_scz2025_magma_prep.sbatch); echo "prep    $J0" ;;&
  prep|genes) J1=$(sub ${J0:+--dependency=afterok:$J0} $G/step10_scz2025_magma_genes.sbatch); echo "genes   $J1 (array 1-5)" ;;&
  prep|genes|finish) J2=$(sub ${J1:+--dependency=afterok:$J1} $G/step10_scz2025_magma_finish.sbatch); echo "finish  $J2" ;;
  *) echo "START must be prep|genes|finish" >&2; exit 2 ;;
esac
squeue -u "$USER" -h -o "%.14i %.20j %.8T %R" | grep magma || true
