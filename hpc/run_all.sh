#!/bin/bash
# Single entry point for the genetics pipeline.
#
#   bash hpc/run_all.sh              # preflight, then run/submit 01-06
#   bash hpc/run_all.sh 01 02        # only those steps
#   DRY_RUN=1 bash hpc/run_all.sh    # print every command, execute nothing
#
# Two modes, chosen automatically:
#
#   sbatch present  -> submits each step as a SLURM job with afterok
#                      dependencies, so the chain runs unattended.  Prints job
#                      IDs and returns immediately.
#   sbatch absent   -> runs each step sequentially in this shell.  This is how
#                      the local test against the synthetic fixture exercises
#                      exactly the same script bodies the cluster will run.
#
# Configure paths with environment variables or hpc/config.local.sh; nothing in
# this file is machine-specific.  See hpc/README.md.
source "$(dirname "$0")/config.sh"

HPC_DIR="$REPO_ROOT/hpc"
# config.sh defines paths but deliberately creates nothing (so that sourcing it
# on a machine without ~/rds does not mkdir a CSD3 tree).  The runner is the
# thing that actually needs $LOG_DIR: without this, `tee` fails on the first
# step and the whole chain stops with a confusing "No such file or directory".
ensure_dirs "$LOG_DIR"
steps=("${@:-01 02 03 04 05 06}")
# Allow either `run_all.sh 01 02` or the default single string.
read -r -a steps <<< "${steps[*]}"

# A case statement rather than `declare -A`: associative arrays are bash 4+ and
# macOS ships bash 3.2, so the array form aborts on a developer laptop with
# "declare: -A: invalid option" while working on the cluster.  Same reason
# 02/03 avoid mapfile.
script_for() {
  case "$1" in
    01) echo "01_grm.sbatch" ;;
    02) echo "02_reml.sbatch" ;;
    03) echo "03_gwas.sbatch" ;;
    04) echo "04_magma.sbatch" ;;
    05) echo "05_ldsc_rg.sbatch" ;;
    06) echo "06_prs.sbatch" ;;
    *)  return 1 ;;
  esac
}

# ---------------------------------------------------------------------------
# Preflight always runs first: it is seconds, and it catches the ID-mismatch and
# stale-manifest failures that otherwise surface hours into a GRM build.
# ---------------------------------------------------------------------------
echo "### preflight"
if ! bash "$HPC_DIR/00_check_inputs.sh"; then
  echo "run_all: preflight failed, nothing submitted." >&2
  exit 1
fi
echo

n_pheno=$(phenotype_names | wc -l | tr -d ' ')

if command -v sbatch >/dev/null 2>&1; then
  echo "### SLURM mode: submitting with afterok dependencies"
  dep=""
  for s in "${steps[@]}"; do
    f=$(script_for "$s") || { echo "unknown step '$s'" >&2; exit 2; }

    # Array steps are sized from the manifest, so adding a phenotype to the
    # export does not require editing an --array range in the sbatch header.
    arr=""
    [[ "$s" == "02" || "$s" == "03" ]] && arr="--array=1-$n_pheno"

    # shellcheck disable=SC2086
    jid=$(sbatch --parsable \
            --account "$SLURM_ACCOUNT" --partition "$SLURM_PARTITION" \
            --chdir "$HPC_DIR" --output "$LOG_DIR/%x_%A_%a.log" \
            $arr $dep "$HPC_DIR/$f")
    # afterok on an array job waits for every task, which is what we want:
    # 03 must not start until the GRM exists, and 04-06 need all sumstats.
    #
    # 04, 05 and 06 are siblings, not a chain: each consumes 03's sumstats and
    # none reads another's output, so they all depend on the LAST GWAS job
    # rather than on each other.  Chaining them would serialise ~6 h of
    # independent work and make one failure block the other two.
    case "$s" in
      04|05|06) : ;;                       # keep $dep pointing at 03
      *) dep="--dependency=afterok:${jid%%;*}" ;;
    esac
    echo "  $s  $f  -> job ${jid%%;*} ${arr:+($arr)}"
  done
  echo
  echo "Monitor with: squeue -u \$USER    Logs in: $LOG_DIR"
else
  echo "### local mode: no sbatch, running sequentially"
  echo "    ($n_pheno phenotypes per step, looped in-process)"
  for s in "${steps[@]}"; do
    f=$(script_for "$s") || { echo "unknown step '$s'" >&2; exit 2; }
    echo
    echo "### $s  $f"
    # SLURM_ARRAY_TASK_ID deliberately unset: the steps loop over all
    # phenotypes when it is absent.
    bash "$HPC_DIR/$f" 2>&1 | tee "$LOG_DIR/${f%.sbatch}.log"
    status=${PIPESTATUS[0]}
    [[ "$status" == "0" ]] || { echo "run_all: $f failed (exit $status)" >&2; exit "$status"; }
  done
  echo
  echo "### all steps completed.  Logs in: $LOG_DIR"
fi
