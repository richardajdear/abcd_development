#!/bin/bash
#SBATCH --job-name=reml70
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-5
#
# README_HPC.md 3 step 5: SNP heritability by GCTA GREML, one phenotype per
# array task.  This is the CANONICAL h2 for the project (4 rule 11): the dense
# imputed GRM restricted to the PC-AiR unrelated set.  Do NOT report
# reml_zaitlen*.
#
# WHY THIS IS A NEW SCRIPT rather than `legacy/hpc/02_reml.sbatch` with
# overrides, as 6 Step 5 suggests.  That script sources v1's config, whose
# config.local.sh contains the PLAIN assignment
#     ABCD_HPC_ROOT="$REPO_ROOT/hpc/work"
# -- a path the 2.2 migration removed, and a plain assignment that would
# silently beat an exported override (the exact trap config.local.sh.example
# warns about).  Every v1-derived path would then point into a directory that
# no longer exists.  So the design is reproduced here against the 7.0 config
# instead; the GCTA call below is v1's, argument for argument, verified against
# legacy/hpc/work/results/reml_imp_pooled/*.log.
#
# The design, from that log:
#   --grm   results/grm_imp_pooled/abcd_imp.pcair   (8,178 ids: the dense
#           imputed GRM already restricted to the PC-AiR unrelated set, so no
#           separate --keep is needed and none was used)
#   --pheno the FID = IID export, because .grm.id spells FID = IID and GCTA
#           matches the FID+IID PAIR literally (rule 1).  Using the family-id
#           export here would silently return a zero-subject REML.
#   covariates: sex + site (discrete), baseline_age + n_visits + 10 PCs
#           (quantitative) -- carried by the export itself.
#
# Expect n around 6,010 (the overlap of the 8,177 unrelated ids with the 8,596
# analysis set), up from v1's 5,649.
#
#   sbatch genetic_analysis/step5_reml.sbatch
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate genetic_analysis/config.sh" >&2; exit 2; }
source "$_cfg"

# The FID = IID sibling of $PHENO_DIR (rule 2).
PHENO_FID="${PHENO_FID:-${PHENO_DIR}_fidiid}"
GRM_UNREL="${GRM_UNREL:-$V1_LEGACY/results/grm_imp_pooled/abcd_imp.pcair}"
REML_DIR="${REML_DIR:-$OUT_V2/reml_imp_pooled}"

PHENO="$PHENO_FID/phenotypes_gcta.txt"
COVAR_QUANT="$PHENO_FID/covar_quant.txt"
COVAR_CAT="$PHENO_FID/covar_categorical.txt"
MANIFEST="$PHENO_FID/phenotype_manifest.tsv"

require_paths PHENO COVAR_QUANT COVAR_CAT MANIFEST GCTA
[[ -f "$GRM_UNREL.grm.bin" && -f "$GRM_UNREL.grm.id" ]] \
  || { echo "FATAL: $GRM_UNREL incomplete" >&2; exit 2; }
ensure_dirs "$REML_DIR"

names=()
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx=$((SLURM_ARRAY_TASK_ID - 1))
  if (( idx >= ${#names[@]} )); then
    echo "array index $SLURM_ARRAY_TASK_ID exceeds ${#names[@]} phenotypes; nothing to do"; exit 0
  fi
  todo=("${names[$idx]}")
else
  todo=("${names[@]}")
fi

for pheno in "${todo[@]}"; do
  mpheno=$(mpheno_for "$pheno") || { echo "ERROR: '$pheno' not in $MANIFEST" >&2; exit 3; }
  echo "=== step5_reml  $pheno  (--mpheno $mpheno)  $(date) ==="
  run "$GCTA" \
    --reml \
    --grm "$GRM_UNREL" \
    --pheno "$PHENO" \
    --mpheno "$mpheno" \
    --qcovar "$COVAR_QUANT" \
    --covar "$COVAR_CAT" \
    --thread-num "${THREADS:-8}" \
    --out "$REML_DIR/$pheno"
done

# Collect: idempotent, rebuilt from whatever .hsq files exist (v1's collector).
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  summary="$REML_DIR/reml_summary.tsv"
  {
    printf 'phenotype\th2\tse\tpval\tn\n'
    for pheno in "${names[@]}"; do
      f="$REML_DIR/$pheno.hsq"
      [[ -f "$f" ]] || continue
      awk -v p="$pheno" '
        $1=="V(G)/Vp" {h2=$2; se=$3}
        $1=="Pval"    {pv=$2}
        $1=="n"       {n=$2}
        END {if (h2!="") printf "%s\t%.4f\t%.4f\t%s\t%s\n", p, h2, se, pv, n}
      ' "$f"
    done
  } > "$summary"
  echo; column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
fi
echo "=== step5_reml done  $(date) ==="
