#!/bin/bash
#SBATCH --job-name=d7_partition
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=directions/d7_partitioned_prs/work/logs/%x_%A_%a.log
#
# D7 step 3: one task per (arm x outcome set).  Tasks = ${#ARMS[@]} x 3 (15 with
# the default 5 arms); set D7_C3AXIS=1 to add the optional c3axis outcome set.
#   sbatch --array=1-15 directions/d7_partitioned_prs/d7_partition.sbatch
# Extra arms (e.g. SCZ25_META_annot|...|pooled from step 4): D7_EXTRA_ARMS, space-separated.
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/directions/d7_partitioned_prs/config.sh
cd "$REPO"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
ALL=("${ARMS[@]}"); [[ -n "${D7_EXTRA_ARMS:-}" ]] && ALL+=(${D7_EXTRA_ARMS})
SETS=(imaging_hcp imaging_dk symptoms); [[ "${D7_C3AXIS:-0}" == 1 ]] && SETS+=(c3axis)
T=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
(( T < ${#ALL[@]} * ${#SETS[@]} )) || { echo "FATAL: task $((T+1)) > ${#ALL[@]}x${#SETS[@]}" >&2; exit 2; }
spec="${ALL[$(( T / ${#SETS[@]} ))]}"; OS="${SETS[$(( T % ${#SETS[@]} ))]}"
name=${spec%%|*}; cell=${spec##*|}
[[ -s "$WORK/gene_scores/$name.npy" ]] || { echo "FATAL: no gene scores for $name (step 2 gate?)" >&2; exit 3; }

base="${name%_annot}"
case "$OS" in
  imaging_hcp) PD="$PHENO_HCP"; R70="$REPO/genetic_analysis/work/results_70tab_hcp" ;;
  imaging_dk)  PD="$PHENO_DK";  R70="$REPO/genetic_analysis/work/results_70tab" ;;
  symptoms)    PD="$PHENO_HCP"; R70="" ;;
  c3axis)      PD="$C3PHENO";   R70="" ;;
esac
EXTRA=()
[[ "$OS" == symptoms ]] && { [[ -s "$SYMPTOMS" ]] || { echo "FATAL: $SYMPTOMS missing (laptop step L1)" >&2; exit 4; }; EXTRA+=(--symptoms "$SYMPTOMS"); }
if [[ -n "$R70" ]]; then
  z=""; [[ "$cell" == pooled ]] && z="_zanc"
  case "$base" in
    SCZ25_*) REF="$R70/prs_scz2025_1lmm/assoc_SBayesRC_${base}${z}.tsv" ;;
    *)       REF="$R70/prs_final_1lmm/assoc_SBayesRC_${base}${z}.tsv" ;;
  esac
  [[ -s "$REF" ]] && EXTRA+=(--assoc-ref "$REF") || echo "WARN: no gate reference $REF"
fi
echo "=== $name ($cell) x $OS  $(date +%T)"
"$D7PY" "$D7/code/03_partition_test.py" --arm "$name" --cell "$cell" --outcomes "$OS" \
  --work "$WORK" --pheno-dir "$PD" --strata "$STRATA" --eur-keep "$EUR_KEEP" ${EXTRA[@]+"${EXTRA[@]}"} \
  --out "$WORK/partition/${name}__${OS}.tsv"
echo "=== done $(date +%T)"
