#!/usr/bin/env bash
# Re-run every locally fitted specification end to end (assemble -> fit ->
# phenotype, plus the GCTA export for the genetic configs).
#
# Written 2026-09-14 for the move from the 6.0-vintage tables to the true 7.0
# tabulation; it is the list of runs the report regenerators
# (tools/regen_*.py) and ahba_pls/ read.  Re-run it whenever the release
# tables change.  HCP-MMP (ct_70_hcp_noglobal_mv2) is deliberately not in the
# default list: run DK first, compare, then add it (see README "Status").
#
#   tools/rerun_local.sh                 # all configs below
#   tools/rerun_local.sh ct_70_baseline  # a subset
#
# PY / RSCRIPT / CORES override the interpreters and the lme4 worker count.
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=src

export PY="${PY:-$HOME/mambaforge/envs/abcd/bin/python}"
export RSCRIPT="${RSCRIPT:-$HOME/mambaforge/envs/abcd/bin/Rscript}"
CORES="${CORES:-8}"

DEFAULT=(
  ct_70_noglobal_mv2_genetic      # the settled specification
  ct_70_noglobal_mv2              # + family random effect (h2 contrast)
  ct_70_noglobal_mv3 ct_70_noglobal_mv3_genetic
  ct_70_noglobal_mv4 ct_70_noglobal_mv4_genetic
  ct_70_global_mv3_genetic        # global-covariate contrast
  ct_70_baseline ct_70_genetic    # global covariate, mv2 (fit_summary rows)
  t1t2_70_noglobal_mv2_genetic    # T1w/T2w maps for ahba_pls
  ct_60_noglobal_mv2_genetic      # the settled spec on the 6.0 tables (comparison only)
)
CONFIGS=("${@:-${DEFAULT[@]}}")

for cfg in "${CONFIGS[@]}"; do
  export ABCD_CONFIG="$cfg"
  echo "=================== $cfg  $(date -u +%FT%TZ)"
  "$PY" -m abcd.assemble
  "$RSCRIPT" R/fit_lmm.R --cores "$CORES"
  "$PY" -m abcd.phenotype
  if [[ "$cfg" == *genetic* ]]; then
    "$PY" -m abcd.gcta_export || echo "gcta_export skipped for $cfg (see message above)"
  fi
done
echo "=================== done $(date -u +%FT%TZ)"
