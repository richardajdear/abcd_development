#!/bin/bash
# H4 cluster leg, from the laptop.  Needs a live CSD3 login (ControlMaster):
#   ssh -o ControlPersist=4h login-q-1.hpc.cam.ac.uk
#   bash ahba_pls/h4_projection/run_on_csd3.sh push|status|pull
set -euo pipefail
HOST=login-q-1.hpc.cam.ac.uk; RREPO=rds/hpc-work/abcd_development; HERE=ahba_pls/h4_projection
cd "$(dirname "$0")/../.."
case "${1:-}" in
  push)
    # rsync the folder itself rather than relying on git pull: the cluster
    # checkout carries its own unpushed commits and may not fast-forward
    ssh "$HOST" "cd $RREPO && mkdir -p $HERE/work $HERE/results slurm"
    rsync -avh --exclude results/ --exclude __pycache__ "$HERE/" "$HOST:$RREPO/$HERE/"
    ssh "$HOST" "cd $RREPO && PARC=hcp sbatch --export=ALL,PARC=hcp $HERE/02_h4_reml.sbatch" ;;
  status) ssh "$HOST" "squeue -u rajd2 -n h4_reml; cd $RREPO && tail -n 4 slurm/h4_reml_*.log 2>/dev/null" ;;
  pull)   rsync -avh "$HOST:$RREPO/$HERE/results/reml/" "$HERE/results/reml/" --include '*.hsq' --include '*.tsv' --include '*.log' --exclude '*' ;;
  *) echo "usage: $0 push|status|pull" >&2; exit 2 ;;
esac
