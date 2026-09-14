#!/bin/bash
#SBATCH --job-name=v2_prsfam
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x_%j.log
#
# Step 6: within-family PRS association (Fulker between/within decomposition).
#
# WHY THIS IS THE STEP THAT MAKES THE PRS RESULT DEFENSIBLE
#
# The population PRS association (v1 06_prs: SCZ -> global_slope, full-sample
# p_adj = 0.034) is confounded by anything family-constant that correlates
# with both score and phenotype: residual stratification, assortative mating,
# genetic nurture.  Siblings differ only by Mendelian segregation, so the
# within-family coefficient b_W is immune to all three.  ABCD's sibling
# enrichment is exactly what powers this test -- the design reason hpc_v2
# exists.
#
# READ THE OUTPUT CAREFULLY (the script prints this too): at ABCD's pair
# count, b_W is underpowered for confirmation; the interpretable readout is
# beta_diff = b_B - b_W, which BOUNDS the confounding in the population
# estimate.  A null b_W alone is NOT evidence against the effect.
#
# INPUTS
#   - score_*.profile from v1's imputed 06_prs run ($PRS_PROFILE_DIR)
#   - the phenotype export, whose FID column carries REAL family IDs
#     (the genotype .fam has FID=IID; the script refuses to run on that)
#   - $EUR_IDS for the EUR-primary arm
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc_v2/config.sh (looked in $(dirname "$0"), $PWD, ${SLURM_SUBMIT_DIR:-unset})" >&2; exit 2; }
source "$_cfg"
require_paths PHENO COVAR_QUANT COVAR_CAT MANIFEST RSCRIPT PRS_PROFILE_DIR
ensure_dirs "$PRSFAM_DIR"

eur_args=()
if [[ -f "$EUR_IDS" ]]; then
  eur_args=(--eur-ids "$EUR_IDS")
else
  echo "WARN: EUR_IDS=$EUR_IDS not found; the script falls back to a PC-based"
  echo "      stand-in for the EUR stratum.  Fine on the fixture, NOT on CSD3."
fi

run "$RSCRIPT" "$RSCRIPT_DIR/06_prs_family.R" \
  --prs-dir "$PRS_PROFILE_DIR" \
  --pheno "$PHENO" \
  --covar-quant "$COVAR_QUANT" \
  --covar-cat "$COVAR_CAT" \
  --manifest "$MANIFEST" \
  ${eur_args[@]+"${eur_args[@]}"} \
  --out "$PRSFAM_DIR/prs_withinfamily.tsv"

if [[ "$DRY_RUN" != "1" ]]; then
  [[ -s "$PRSFAM_DIR/prs_withinfamily.tsv" ]] \
    || { echo "FATAL: $PRSFAM_DIR/prs_withinfamily.tsv missing/empty" >&2; exit 4; }
  n_rows=$(awk 'NR>1' "$PRSFAM_DIR/prs_withinfamily.tsv" | wc -l | tr -d ' ')
  echo "prs_withinfamily.tsv: $n_rows rows"
fi

echo "=== 06_prs_family done  $(date) ==="
