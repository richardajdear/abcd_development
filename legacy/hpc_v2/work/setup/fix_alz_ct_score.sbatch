#!/bin/bash
#SBATCH --job-name=v2_alzct
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#
# Repair the APOE-INCLUSIVE Alzheimer's C+T score, and re-run both the
# population and the within-family association for the control disorders.
#
# WHY IT WAS BROKEN: PGC-ALZ2 carries rs283815 with z = inf.  PLINK loads it
# without complaint ("50,440 valid predictors loaded") and then returns
# SCORESUM = 0 for all 11,670 subjects, so the with-APOE arm silently produced
# no score at all and prs_assoc.R emitted rows for ALZnoAPOE only.  The
# APOE-EXCLUDED arm was unaffected because dropping the APOE region also drops
# the offending variant -- which is why the defect stayed invisible.
#
# The clumping is untouched (it reads SNP/P, not the effect), so only the
# weight files and the scoring need redoing.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
source "$REPO/hpc/config.sh"
W="$REPO/hpc/work"
OUT_V2="$REPO/hpc_v2/work/results_v2"
OUTD="$OUT_V2/control_ALZ"; PRSD="$OUTD/prs"
PRS_GENO="$W/genotype_imputed/prs/abcd_imp_prs"
PHENO="$W/pheno_allanc_prs/phenotypes_gcta.txt"
COVAR_QUANT="$W/pheno_allanc_prs/covar_quant.txt"
COVAR_CAT="$W/pheno_allanc_prs/covar_categorical.txt"
MANIFEST="$W/pheno_allanc_prs/phenotype_manifest.tsv"
EUR_IDS="$W/results/ancestry/eur_anchor.keep"
RSCRIPT=/rds/user/rajd2/hpc-work/envs/genesis/bin/Rscript

echo "=== rebuild ALZ weights with non-finite effects dropped  $(date) ==="
for thr in $PRS_THRESHOLDS; do
  tag="${thr//./p}"
  awk -v t="$thr" '
    NR==FNR {keep[$1]=1; next}
    FNR==1 {next}
    ($1 in keep) && ($5+0) < (t+0) {
      b=$4
      if (b=="" || b=="inf" || b=="-inf" || b=="nan" || b=="NA") next
      if (b+0 != b+0) next
      if (b+0 > 1e6 || b+0 < -1e6) next
      print $1"\t"$2"\t"b }' \
    "$PRSD/ALZ.index_snps" "$OUTD/ALZ.tsv" > "$PRSD/ALZ_$tag.weights"
  n=$(wc -l < "$PRSD/ALZ_$tag.weights")
  "$W/bin/plink" --bfile "$PRS_GENO" --score "$PRSD/ALZ_$tag.weights" 1 2 3 sum \
    --out "$PRSD/score_ALZ_$tag" >/dev/null 2>&1 || true
  sd=$(awk 'NR>1{s+=$NF; ss+=$NF*$NF; n++} END{m=s/n; printf "%.4g", sqrt(ss/n-m*m)}' "$PRSD/score_ALZ_$tag.profile")
  echo "  p<$thr: $n SNPs, score sd = $sd"
done

echo "=== re-associate (population + within-family), both ALZ arms  $(date) ==="
"$RSCRIPT" "$REPO/tools/prs_assoc.R" --prs-dir "$PRSD" --pheno "$PHENO" \
  --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
  --eur-ids "$EUR_IDS" --out "$OUTD/prs_association_ALZ.tsv"
"$RSCRIPT" "$REPO/hpc_v2/R/06_prs_family.R" --prs-dir "$PRSD" --pheno "$PHENO" \
  --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
  --eur-ids "$EUR_IDS" --out "$OUTD/prs_withinfamily_ALZ.tsv"
echo "=== done  $(date) ==="
