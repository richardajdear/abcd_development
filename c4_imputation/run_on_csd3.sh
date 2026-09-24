#!/bin/bash
# The cluster leg, driven from the laptop.  Needs a live CSD3 login first
# (ControlMaster is on in ~/.ssh/config, so one interactive ssh carries the rest):
#
#   ssh -o ControlPersist=4h login-q-1.hpc.cam.ac.uk     # authenticate, leave open
#   bash c4_imputation/run_on_csd3.sh setup    # git pull, fetch panel + jars, build the c4 env (login node, once)
#   bash c4_imputation/run_on_csd3.sh submit   # c4_impute -> (afterok) c4_assoc
#   bash c4_imputation/run_on_csd3.sh status
#   bash c4_imputation/run_on_csd3.sh pull     # summary tables back into c4_imputation/results/
#
# Nothing per-subject comes back: pull copies results/*.tsv|*.txt only.
set -euo pipefail
HOST="${CSD3_HOST:-login-q-1.hpc.cam.ac.uk}"
RREPO=rds/hpc-work/abcd_development
cd "$(dirname "$0")/.."
case "${1:-}" in
  setup)
    ssh "$HOST" "cd $RREPO && git pull --ff-only && bash c4_imputation/00_fetch_resources.sh && bash c4_imputation/env/setup_csd3_env.sh" ;;
  submit)
    ssh "$HOST" "cd $RREPO && mkdir -p slurm && j=\$(sbatch --parsable c4_imputation/c4_impute.sbatch) && \
      k=\$(sbatch --parsable --dependency=afterok:\$j c4_imputation/c4_assoc.sbatch) && \
      echo \"c4_impute \$j  c4_assoc \$k\" | tee -a c4_imputation/work/submitted_jobs.txt" ;;
  status)
    ssh "$HOST" "squeue -u \$USER -n c4_impute,c4_assoc; cd $RREPO && tail -n 5 slurm/c4_*.log 2>/dev/null" ;;
  pull)
    rsync -avh "$HOST:$RREPO/c4_imputation/results/" c4_imputation/results/ --include='*.tsv' --include='*.txt' --exclude='*' ;;
  *) echo "usage: $0 setup|submit|status|pull" >&2; exit 2 ;;
esac
