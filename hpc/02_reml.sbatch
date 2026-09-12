#!/bin/bash
#SBATCH --job-name=abcd_reml
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-5
#
# Step 2: SNP heritability by GREML, one phenotype per array task.
#
# The array index selects a row of $MANIFEST, so the set of phenotypes is
# defined by the export -- not by an --array range that has to be kept in sync
# by hand.  The range above must cover the manifest; run_all.sh sets it from the
# file.  Outside SLURM, loops over all phenotypes instead.
#
# Runs on the UNRELATED subset (see 01_grm): with relatives retained, shared
# environment loads onto the additive term.
# Locate config.sh.  Under sbatch, $0 is a COPY of this script in
# /var/spool/slurm/slurmd/jobNNN/, so `dirname "$0"` does not contain
# config.sh -- run_all.sh's SLURM mode therefore sourced nothing, and because
# `set -euo pipefail` lives INSIDE config.sh it was never enabled either, so
# every step ran to completion with undefined functions and exited 0.  Search
# the submit directory too ($PWD, which run_all.sh sets via --chdir).
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc/config.sh (looked in $(dirname "$0"), $PWD, ${SLURM_SUBMIT_DIR:-unset})" >&2; exit 2; }
source "$_cfg"
require_paths GRM_UNREL PHENO COVAR_QUANT COVAR_CAT MANIFEST GCTA
ensure_dirs "$REML_DIR"

# Which phenotypes this invocation handles: one row under SLURM, all rows not.
# `while read` rather than `mapfile`: mapfile is bash 4+, and macOS ships bash
# 3.2, so mapfile fails with "command not found" on a developer laptop while
# working fine on the cluster.  Portability is free here.
names=()
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx=$((SLURM_ARRAY_TASK_ID - 1))
  if (( idx >= ${#names[@]} )); then
    echo "array index $SLURM_ARRAY_TASK_ID exceeds ${#names[@]} phenotypes; nothing to do"
    exit 0
  fi
  todo=("${names[$idx]}")
else
  todo=("${names[@]}")
fi

for pheno in "${todo[@]}"; do
  mpheno=$(mpheno_for "$pheno") || { echo "ERROR: '$pheno' not in $MANIFEST" >&2; exit 3; }
  echo "=== 02_reml  $pheno  (--mpheno $mpheno)  $(date) ==="

  # $GRM_UNREL is the pruned GRM written by 01_grm; its .grm.id is the
  # unrelated subject list, so no separate --keep is needed.
  run "$GCTA" \
    --reml \
    --grm "$GRM_UNREL" \
    --pheno "$PHENO" \
    --mpheno "$mpheno" \
    --qcovar "$COVAR_QUANT" \
    --covar "$COVAR_CAT" \
    --thread-num "$THREADS" \
    --out "$REML_DIR/$pheno"
done

# ---------------------------------------------------------------------------
# Collect.  GCTA writes one .hsq per phenotype; the numbers that matter are
# three lines inside it.  Without this, comparing five phenotypes means opening
# five files, and the array's logs are separate too.
#
# Collection is idempotent and rebuilt from whatever .hsq files exist, so it is
# safe to run after a partial array and re-run after the rest completes.
# ---------------------------------------------------------------------------
if [[ "$DRY_RUN" != "1" ]]; then
  summary="$REML_DIR/reml_summary.tsv"
  {
    echo -e "phenotype\th2\tse\tpval\tn"
    for pheno in "${names[@]}"; do
      f="$REML_DIR/$pheno.hsq"
      [[ -f "$f" ]] || continue
      # .hsq layout: "V(G)/Vp<TAB>estimate<TAB>SE"; Pval and n are their own rows.
      awk -v p="$pheno" '
        $1=="V(G)/Vp" {h2=$2; se=$3}
        $1=="Pval"    {pv=$2}
        $1=="n"       {n=$2}
        END {if (h2!="") printf "%s\t%.4f\t%.4f\t%s\t%s\n", p, h2, se, pv, n}
      ' "$f"
    done
  } > "$summary"
  echo
  column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
  echo
  echo "NOTE: n above is the unrelated subset, far below the phenotyped sample."
  echo "      A wide SE here is the design, not a failure -- see hpc/README_HPC.md."
fi

echo "=== 02_reml done  $(date) ==="
