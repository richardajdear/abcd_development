#!/bin/bash
# The cluster leg, from the laptop.  Needs a live CSD3 login: open one first
# (ControlMaster is on in ~/.ssh/config, so one interactive `ssh` lets the
# non-interactive calls below ride the same connection):
#
#   ssh -o ControlPersist=4h login-q-1.hpc.cam.ac.uk      # authenticate, then leave it open
#   bash prs_beta_map/prs_scz_beta_map_hcp/run_on_csd3.sh push      # slopes + scripts up, submit
#   bash prs_beta_map/prs_scz_beta_map_hcp/run_on_csd3.sh pull      # beta maps back
#   python prs_beta_map/prs_scz_beta_map_hcp/03_spin_test.py && python .../04_figure.py
set -euo pipefail
HOST=login-q-1.hpc.cam.ac.uk
RREPO=rds/hpc-work/abcd_development
HERE=prs_beta_map/prs_scz_beta_map_hcp
cd "$(dirname "$0")/../.."
case "${1:-}" in
  push)
    rsync -avh --progress "$HERE/" "$HOST:$RREPO/$HERE/" --exclude results/ --exclude '__pycache__'
    ssh "$HOST" "cd $RREPO && git pull --ff-only || true; PARC=hcp sbatch --export=ALL,PARC=hcp $HERE/02_parcel_prs_assoc.sbatch" ;;
  status)
    ssh "$HOST" "squeue -u rajd2 -n prs_beta_map; cd $RREPO && tail -n 3 slurm/prs_beta_map_*.log 2>/dev/null" ;;
  pull)
    rsync -avh "$HOST:$RREPO/$HERE/results/beta_map_*.tsv" "$HERE/results/" ;;
  *) echo "usage: $0 push|status|pull" >&2; exit 2 ;;
esac
