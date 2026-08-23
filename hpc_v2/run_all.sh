#!/bin/bash
# Single entry point for the hpc_v2 pipeline.
#
#   bash hpc_v2/run_all.sh              # preflight, then run/submit 01-06
#   bash hpc_v2/run_all.sh 01 02        # only those steps
#   DRY_RUN=1 bash hpc_v2/run_all.sh    # print every command, execute nothing
#
# Two modes, chosen automatically (v1 convention):
#   sbatch present -> submits with afterok dependencies, returns immediately.
#   sbatch absent  -> runs each step sequentially in this shell; identical
#                     script bodies, which is what the local test exercises.
#
# DEPENDENCY SHAPE
#   01 (GDS) -> 02 (kinship) -> 03 (null models) -> 04 (assoc + collect)
#   05 (Zaitlen REML) depends only on the GRM + export: submitted immediately.
#   06 (within-family PRS) depends only on v1 score profiles + export:
#      submitted immediately.
#   So 05/06 run in parallel with the GENESIS chain, not after it.
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc_v2/config.sh" >&2; exit 2; }
source "$_cfg"

HPC_DIR="$REPO_ROOT/hpc_v2"
ensure_dirs "$LOG_DIR"
steps=("${@:-01 02 03 04 05 06}")
read -r -a steps <<< "${steps[*]}"

script_for() {   # case, not declare -A: bash 3.2 (macOS) compatibility
  case "$1" in
    01) echo "01_gds.sbatch" ;;
    02) echo "02_kinship.sbatch" ;;
    03) echo "03_null_model.sbatch" ;;
    04) echo "04_assoc.sbatch" ;;
    05) echo "05_reml_zaitlen.sbatch" ;;
    06) echo "06_prs_family.sbatch" ;;
    *)  return 1 ;;
  esac
}

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

    arr=""
    case "$s" in
      01) arr="--array=0-22" ;;
      03) arr="--array=1-$n_pheno" ;;
      04) arr="--array=1-$((n_pheno * 22))" ;;
      05) arr="--array=1-$n_pheno" ;;
    esac

    # 05 and 06 join the GENESIS chain nowhere: no --dependency.
    this_dep="$dep"
    case "$s" in 05|06) this_dep="" ;; esac

    # shellcheck disable=SC2086
    jid=$(sbatch --parsable \
            --account "$SLURM_ACCOUNT" --partition "$SLURM_PARTITION" \
            --chdir "$HPC_DIR" --output "$LOG_DIR/%x_%A_%a.log" \
            $arr $this_dep "$HPC_DIR/$f")
    jid="${jid%%;*}"
    case "$s" in
      01|02|03) dep="--dependency=afterok:$jid" ;;
      04)
        # Final collection pass AFTER the whole assoc array: the inline
        # chr==22 collect can see a partial set; this one cannot.
        col=$(sbatch --parsable \
                --account "$SLURM_ACCOUNT" --partition "$SLURM_PARTITION" \
                --chdir "$HPC_DIR" --output "$LOG_DIR/%x_%j.log" \
                --job-name=v2_collect --time=00:30:00 --mem=16G \
                --dependency="afterok:$jid" \
                --wrap "bash -c 'source $HPC_DIR/config.sh; for p in \$(phenotype_names); do \"\$RSCRIPT\" \"\$RSCRIPT_DIR/05_collect_assoc.R\" --assoc-dir \"\$ASSOC_DIR\" --pheno \"\$p\" --out-dir \"\$ASSOC_DIR\" --maf \"\$ASSOC_MAF\" --mac \"\$ASSOC_MAC\"; done'")
        echo "  04c v2_collect -> job ${col%%;*} (afterok:$jid)"
        ;;
    esac
    echo "  $s  $f  -> job $jid ${arr:+($arr)} ${this_dep:+[$this_dep]}"
  done
  echo
  echo "Monitor with: squeue -u \$USER    Logs in: $LOG_DIR"
  echo "REMEMBER (v1 lesson): verify *_summary.tsv outputs, not job state."
else
  echo "### local mode: no sbatch, running sequentially"
  for s in "${steps[@]}"; do
    f=$(script_for "$s") || { echo "unknown step '$s'" >&2; exit 2; }
    echo
    echo "### $s  $f"
    bash "$HPC_DIR/$f" 2>&1 | tee "$LOG_DIR/${f%.sbatch}.log"
    status=${PIPESTATUS[0]}
    [[ "$status" == "0" ]] || { echo "run_all: $f failed (exit $status)" >&2; exit "$status"; }
  done
  echo
  echo "### all steps completed.  Logs in: $LOG_DIR"
fi
