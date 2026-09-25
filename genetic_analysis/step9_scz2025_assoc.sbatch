#!/bin/bash
#SBATCH --job-name=scz25_assoc
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm/%x_%A.log
#
# Step 9, association layer, ONE parcellation per job (PARC=dsk|hcp in the
# environment, exactly as every other step: config.local.sh turns it into
# $PHENO_DIR and $OUT_V2).  For every (method, trait-arm) score directory
# under $SCZ25_ROOT:
#   tools/prs_assoc.R        population association, EUR and full strata
#   R/06_prs_family.R        Fulker between/within-family decomposition
#   (both on the raw score and on the _zanc within-ancestry-standardised one)
#   R/09_prs_ancestry_strata.R   per-cluster fits + PRS x stratum heterogeneity
# Outputs: $OUT_V2/prs_scz2025/{assoc,fam,strata}_<METH>_<TRAIT>[_zanc].tsv
#
#   PARC=dsk sbatch --export=ALL,PARC=dsk genetic_analysis/step9_scz2025_assoc.sbatch
#   PARC=hcp sbatch --export=ALL,PARC=hcp genetic_analysis/step9_scz2025_assoc.sbatch
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/genetic_analysis/config.sh"
echo "PARC=$PARC  PHENO_DIR=$PHENO_DIR  OUT_V2=$OUT_V2"
ROOT="${SCZ25_ROOT:-$ABCD_HPC_ROOT/scores_scz2025}"
# Reused by step 15 (new GWAS traits): PRS_TAG names the output dir, ARM_GLOB the
# score-directory pattern; the defaults are step 9 exactly.
FINAL="$OUT_V2/${PRS_TAG:-prs_scz2025}"; ARM_GLOB="${ARM_GLOB:-SCZ25_*}"
ASSOC_R="$REPO/tools/prs_assoc.R"
FAM_R="$REPO/genetic_analysis/R/06_prs_family.R"
STRATA_R="$REPO/genetic_analysis/R/09_prs_ancestry_strata.R"
STRATA_TSV="${STRATA_TSV:-$V1_LEGACY/results/ancestry/strata_k4.tsv}"
PHENO="$PHENO_DIR/phenotypes_gcta.txt"; COVQ="$PHENO_DIR/covar_quant.txt"
COVC="$PHENO_DIR/covar_categorical.txt"; MANIF="$PHENO_DIR/phenotype_manifest.tsv"
require_paths PHENO COVQ COVC MANIF EUR_IDS RSCRIPT ASSOC_R FAM_R STRATA_R STRATA_TSV
ensure_dirs "$FINAL"

n_dirs=0
for d in "$ROOT"/{CT,PRSCS,SBayesR,SBayesRC,PRSCSX}/$ARM_GLOB; do
  [[ -d "$d" ]] || continue
  METH="$(basename "$(dirname "$d")")"; TRAIT="$(basename "$d")"
  n_prof=$(ls "$d"/score_*.profile 2>/dev/null | wc -l); (( n_prof > 0 )) || continue
  n_dirs=$((n_dirs+1))
  echo "=== $METH / $TRAIT : $n_prof profiles  $(date)"
  run "$RSCRIPT" "$ASSOC_R" --prs-dir "$d" --pheno "$PHENO" --covar-quant "$COVQ" \
    --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EUR_IDS" \
    --out "$FINAL/assoc_${METH}_${TRAIT}.tsv" || echo "WARN: assoc failed for $METH/$TRAIT"
  run "$RSCRIPT" "$FAM_R" --prs-dir "$d" --pheno "$PHENO" --covar-quant "$COVQ" \
    --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EUR_IDS" \
    --out "$FINAL/fam_${METH}_${TRAIT}.tsv" || echo "WARN: family failed for $METH/$TRAIT"
  if [[ -d "$d/_zanc" ]]; then
    run "$RSCRIPT" "$ASSOC_R" --prs-dir "$d/_zanc" --pheno "$PHENO" --covar-quant "$COVQ" \
      --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EUR_IDS" \
      --out "$FINAL/assoc_${METH}_${TRAIT}_zanc.tsv" || echo "WARN: assoc zanc failed for $METH/$TRAIT"
    run "$RSCRIPT" "$FAM_R" --prs-dir "$d/_zanc" --pheno "$PHENO" --covar-quant "$COVQ" \
      --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EUR_IDS" \
      --out "$FINAL/fam_${METH}_${TRAIT}_zanc.tsv" || echo "WARN: family zanc failed for $METH/$TRAIT"
  fi
  run "$RSCRIPT" "$STRATA_R" --prs-dir "$d" --pheno "$PHENO" --covar-quant "$COVQ" \
    --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EUR_IDS" --strata "$STRATA_TSV" \
    --out "$FINAL/strata_${METH}_${TRAIT}.tsv" || echo "WARN: strata failed for $METH/$TRAIT"
done
(( n_dirs > 0 )) || { echo "FATAL: no score directories under $ROOT" >&2; exit 3; }
echo "=== $n_dirs score directories done; $(ls "$FINAL"/*.tsv | wc -l) tables in $FINAL  $(date)"
