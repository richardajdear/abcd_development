#!/bin/bash
# Submit the whole step-9 chain (Nature 2025 SCZ GWAS -> PRS -> both atlases)
# with afterok dependencies.  Run from anywhere; logs land in <repo>/slurm/.
#   bash genetic_analysis/run_scz2025.sh            # everything
#   START=gather bash genetic_analysis/run_scz2025.sh   # resume from a stage
# Stages: normalise -> score(heavy, chr) -> gather -> assoc(dsk, hcp) -> minp -> collect
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
cd "$REPO"; mkdir -p slurm
G=genetic_analysis
START="${START:-normalise}"
dep() { [[ -n "${1:-}" ]] && printf -- '--dependency=afterok:%s' "$1" || true; }
sub() { sbatch --parsable "$@" | cut -d';' -f1; }
J0=""; J1=""; J2=""; J3=""; J4=""; J5=""; J6=""
stage_ge() { local order=(normalise score gather assoc minp collect) i s=0 t=0
  for i in "${!order[@]}"; do [[ "${order[$i]}" == "$START" ]] && s=$i; [[ "${order[$i]}" == "$1" ]] && t=$i; done; (( t >= s )); }
if stage_ge normalise; then J0=$(sub $G/step9_scz2025_normalise.sbatch); echo "normalise      $J0"; fi
if stage_ge score; then
  J1=$(sub $(dep "$J0") --array=1-8 $G/step9_scz2025_score.sbatch heavy);  echo "score heavy    $J1 (array 1-8)"
  J2=$(sub $(dep "$J0") --array=1-66 --cpus-per-task=4 --mem=24G --time=08:00:00 $G/step9_scz2025_score.sbatch chr)
  echo "score chr      $J2 (array 1-66: PRS-CS EUR/META x22, PRS-CSx x22)"
fi
if stage_ge gather; then
  D=""; [[ -n "$J1" ]] && D="$J1"; [[ -n "$J2" ]] && D="${D:+$D:}$J2"
  J3=$(sub $(dep "$D") $G/step9_scz2025_gather.sbatch); echo "gather         $J3"
fi
if stage_ge assoc; then
  J4=$(sub $(dep "$J3") --export=ALL,PARC=dsk --job-name=scz25_assoc_dsk $G/step9_scz2025_assoc.sbatch); echo "assoc dsk      $J4"
  J5=$(sub $(dep "$J3") --export=ALL,PARC=hcp --job-name=scz25_assoc_hcp $G/step9_scz2025_assoc.sbatch); echo "assoc hcp      $J5"
fi
if stage_ge minp; then
  D=""; [[ -n "$J4" ]] && D="$J4"; [[ -n "$J5" ]] && D="${D:+$D:}$J5"
  J6=$(sub $(dep "$D") $G/step9_scz2025_minp.sbatch); echo "minp           $J6 (array 1-2: dsk, hcp)"
fi
if stage_ge collect; then
  # afterany: the tables are wanted even if a permutation task times out
  D=""; [[ -n "$J6" ]] && D="$J6"
  J7=$(sub ${D:+--dependency=afterany:$D} $G/step9_scz2025_collect.sbatch); echo "collect        $J7"
fi
echo; squeue -u "$USER" -o "%.10i %.22j %.8T %.10M %.6D %R" | grep -i "scz25\|JOBID" || true
