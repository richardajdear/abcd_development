#!/bin/bash
# D1 step 9: GENESIS single-LMM scan of the geometry slopes/intercepts, both arms.
# Same scripts and settings as scan_1lmm (03_null_model / 04_assoc / 05_collect,
# HCP sample set, sex + site + baseline_age + n_visits + 10 PC-AiR PCs); only the
# phenotype file and the output root differ.  Outputs under
# work/results_70tab_hcp/scan_geom/{nullmodel,assoc}{,_eur}.
#
#   bash genetic_analysis/geom_gwas/run_geom.sh eur    slopes   # 5 slopes, EUR arm (kinship_eur sample set)
#   bash genetic_analysis/geom_gwas/run_geom.sh pooled slopes   # 5 slopes, pooled arm (8,596)
#   bash genetic_analysis/geom_gwas/run_geom.sh eur    ints     # 5 intercepts (second priority)
#   bash genetic_analysis/geom_gwas/run_geom.sh reml            # GCTA GREML, 5 slopes (step5_reml.sbatch)
#
# Manifest rows 1-5 are the slopes, 6-10 the intercepts (00_build_pheno.py).
set -euo pipefail
ARM="${1:?eur|pooled|reml}"; SET="${2:-slopes}"
export PARC=hcp
REPO=/home/rajd2/rds/hpc-work/abcd_development; GA="$REPO/genetic_analysis"
source "$GA/config.sh"
W="$OUT_V2/scan_geom"; mkdir -p "$W/logs"
export PHENO_DIR="$W/pheno" LOG_DIR="$W/logs"
export PHENO="$PHENO_DIR/phenotypes_gcta.txt" COVAR_QUANT="$PHENO_DIR/covar_quant.txt" \
       COVAR_CAT="$PHENO_DIR/covar_categorical.txt" MANIFEST="$PHENO_DIR/phenotype_manifest.tsv"
SB=(sbatch --parsable --export=ALL,DRY_RUN=0 --account "$SLURM_ACCOUNT" --partition "$SLURM_PARTITION" --chdir "$GA")

if [[ "$ARM" == reml ]]; then
  export PHENO_FID="$W/pheno_fidiid" REML_DIR="$W/reml_imp_pooled"
  j=$("${SB[@]}" --array=1-5 --output "$LOG_DIR/%x_%A_%a.log" "$GA/step5_reml.sbatch")
  echo "reml -> $j"; exit 0
fi

case "$ARM" in
  eur)    export KIN_DIR="$KIN_DIR_EUR" NULL_DIR="$W/nullmodel_eur" ASSOC_DIR="$W/assoc_eur" ;;
  pooled) export KIN_DIR="$V2_LEGACY/kinship" NULL_DIR="$W/nullmodel" ASSOC_DIR="$W/assoc" ;;
  *) echo "ARM must be eur|pooled|reml" >&2; exit 2 ;;
esac
case "$SET" in
  slopes) lo=1; hi=5 ;;
  ints)   lo=6; hi=10 ;;
  *) echo "SET must be slopes|ints" >&2; exit 2 ;;
esac
names=$(awk -v lo=$lo -v hi=$hi 'NR>1 && $2>=lo && $2<=hi {print $1}' "$MANIFEST" | tr '\n' ' ')
echo "arm=$ARM set=$SET KIN_DIR=$KIN_DIR -> $NULL_DIR, $ASSOC_DIR"; echo "phenotypes: $names"

jn=$("${SB[@]}" --array=$lo-$hi --output "$LOG_DIR/%x_%A_%a.log" "$GA/03_null_model.sbatch")
ja=$("${SB[@]}" --array=$(( (lo-1)*22+1 ))-$(( hi*22 )) --dependency=afterok:$jn \
       --output "$LOG_DIR/%x_%A_%a.log" "$GA/04_assoc.sbatch")
# phenotype names go positionally (rule 10: --export splits on commas)
jc=$("${SB[@]}" --dependency=afterok:$ja --output "$LOG_DIR/%x_%j.log" \
       "$GA/geom_gwas/collect.sbatch" "$ASSOC_DIR" $names)
echo "null $jn -> assoc $ja -> collect $jc"
