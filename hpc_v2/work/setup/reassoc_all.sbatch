#!/bin/bash
#SBATCH --job-name=v2_reassoc
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#
# Re-run ONLY the association step over the scores already computed, because the
# models were missing their age covariate.
#
# prs_assoc.R and 06_prs_family.R both asked for `age_c`; the covariate export
# calls that column `baseline_age`; and both scripts filter requested covariates
# down to those actually present -- silently.  So every PRS association in v1 and
# v2 was fitted as
#     phenotype ~ scale(PRS) + sex + PC1..PC10 + (1 | family_id)
# with NO age term, while the documented model claimed age was included.
#
# Scoring is unaffected, so nothing needs re-scoring: this walks every existing
# score directory and redoes just the regressions.
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/work/setup/paths.sh
shopt -s nullglob
for mdir in "$FINAL"/{CT,PRSCS,SBayesR,SBayesRC}; do
  [[ -d "$mdir" ]] || continue
  METH=$(basename "$mdir")
  for tdir in "$mdir"/*/; do
    TRAIT=$(basename "$tdir")
    ls "$tdir"score_*.profile >/dev/null 2>&1 || continue
    echo "--- $METH / $TRAIT"
    "$RSCRIPT" "$ASSOC_R" --prs-dir "$tdir" --pheno "$PHENO" \
      --covar-quant "$COVQ" --covar-cat "$COVC" --manifest "$MANIF" \
      --eur-ids "$EURKEEP" --out "$FINAL/assoc_${METH}_${TRAIT}.tsv" 2>&1 | grep -vE "^$" | tail -3
    "$RSCRIPT" "$FAM_R" --prs-dir "$tdir" --pheno "$PHENO" \
      --covar-quant "$COVQ" --covar-cat "$COVC" --manifest "$MANIF" \
      --eur-ids "$EURKEEP" --out "$FINAL/fam_${METH}_${TRAIT}.tsv" 2>&1 | tail -2
    # tdir ends in "/", so ${tdir}_zanc is <cell>/_zanc
    if [[ -d "${tdir}_zanc" ]]; then
      ZD="${tdir}_zanc"
      ls "$ZD"/score_*.profile >/dev/null 2>&1 || continue
      "$RSCRIPT" "$ASSOC_R" --prs-dir "$ZD" --pheno "$PHENO" \
        --covar-quant "$COVQ" --covar-cat "$COVC" --manifest "$MANIF" \
        --eur-ids "$EURKEEP" --out "$FINAL/assoc_${METH}_${TRAIT}_zanc.tsv" 2>&1 | tail -2
      "$RSCRIPT" "$FAM_R" --prs-dir "$ZD" --pheno "$PHENO" \
        --covar-quant "$COVQ" --covar-cat "$COVC" --manifest "$MANIF" \
        --eur-ids "$EURKEEP" --out "$FINAL/fam_${METH}_${TRAIT}_zanc.tsv" 2>&1 | tail -2
    fi
  done
done
"$PYTHON" "$V2ROOT/work/setup/collect_final.py" "$FINAL"
