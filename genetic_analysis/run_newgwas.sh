#!/bin/bash
# Step 15 chain: new discovery GWAS (bipolar 2024 EUR + multi-ancestry, ADHD 2023,
# intelligence 2018) through the PRS layer, both atlases, both constructions.
#   bash genetic_analysis/run_newgwas.sh
# normalise -> score (heavy 1-12, chr 1-88) -> gather -> assoc (dsk, hcp; per-region)
#   + assoc_1lmm (dsk, hcp) -> minp -> collect
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development; cd "$REPO"; mkdir -p slurm
G=genetic_analysis; ROOT="$REPO/$G/work/scores_newgwas"
sub() { sbatch --parsable "$@" | cut -d';' -f1; }
PY=legacy/hpc/work/envs/abcd/bin/python
# NORM_ARMS: arms to (re)normalise; default all.  NORM_ARMS=ADHD re-runs one.
J0=$(sub --job-name=newgwas_norm --account=vertes-sl3-cpu --partition=icelake --time=02:00:00 --mem=48G --cpus-per-task=2 \
       --output=slurm/%x_%A.log --wrap "cd $REPO && $PY $G/step15_newgwas_normalise.py ${NORM_ARMS:-}")
J1=$(sub --dependency=afterok:$J0 --array=1-12 $G/step15_newgwas_score.sbatch heavy)
J2=$(sub --dependency=afterok:$J0 --array=1-88 --cpus-per-task=4 --mem=24G --time=08:00:00 $G/step15_newgwas_score.sbatch chr)
J3=$(sub --dependency=afterok:$J1:$J2 $G/step15_newgwas_gather.sbatch)
EXP="ALL,SCZ25_ROOT=$ROOT,PRS_TAG=prs_newgwas,ARM_GLOB=*"
J4=$(sub --dependency=afterok:$J3 --export=$EXP,PARC=dsk --job-name=newgwas_assoc_dsk $G/step9_scz2025_assoc.sbatch)
J5=$(sub --dependency=afterok:$J3 --export=$EXP,PARC=hcp --job-name=newgwas_assoc_hcp $G/step9_scz2025_assoc.sbatch)
J6=$(sub --dependency=afterok:$J3 --export=$EXP --job-name=newgwas_assoc_1lmm $G/step9_scz2025_assoc_1lmm.sbatch)
J7=$(sub --dependency=afterok:$J4:$J5 --job-name=newgwas_minp --account=vertes-sl3-cpu --partition=icelake --time=06:00:00 --mem=32G \
       --cpus-per-task=8 --output=slurm/%x_%A.log --wrap "cd $REPO && OMP_NUM_THREADS=8 $PY $G/step15_newgwas_minp.py")
J8=$(sub --dependency=afterany:$J6:$J7 --job-name=newgwas_collect --account=vertes-sl3-cpu --partition=icelake --time=00:30:00 --mem=8G \
       --output=slurm/%x_%A.log --wrap "cd $REPO && $PY $G/step15_newgwas_collect.py && $PY $G/fig1_prep_prs.py")
printf 'normalise %s\nscore heavy %s (1-12)\nscore chr %s (1-88)\ngather %s\nassoc dsk %s  hcp %s  1lmm %s (1-2)\nminp %s\ncollect %s\n' \
  "$J0" "$J1" "$J2" "$J3" "$J4" "$J5" "$J6" "$J7" "$J8"
