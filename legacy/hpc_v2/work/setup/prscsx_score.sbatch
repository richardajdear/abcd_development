#!/bin/bash
#SBATCH --job-name=v2_csxsc
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#
# Score and associate the PRS-CSx posteriors.
#
# PRS-CSx returns ONE posterior per discovery population, deliberately: the
# intended use is to score each and learn the combination weights in a
# validation set.  We have no independent validation set, so inventing weights
# by fitting them on the same sample we then test would be circular.  Three
# scores are therefore reported, all pre-specified:
#   EUR   the European posterior alone -- the direct like-for-like against
#         PRS-CS, differing only in that PRS-CSx coupled it to the EAS data
#   EAS   the East Asian posterior alone, for completeness
#   SUM   the unweighted sum of the two, the standard no-validation default
# Any of these beating the others is not evidence of anything by itself; the
# comparison that matters is EUR-here against PRS-CS's EUR-only score.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
OUTD="$OUT_V2/prscsx"
PRS_GENO="$W/genotype_imputed/prs/abcd_imp_prs"

# PRS-CSx output columns: CHR SNP BP A1 A2 BETA
for pop in EUR EAS; do
  cat "$OUTD"/SCZ_${pop}_pst_eff_a1_b0.5_phiauto_chr*.txt \
    | awk 'BEGIN{OFS="\t"} {print $2,$4,$6}' > "$OUTD/w_${pop}.txt"
  echo "$pop: $(wc -l < "$OUTD/w_${pop}.txt") SNPs"
done
# Unweighted sum: add posteriors on the shared effect allele.  A SNP present in
# only one population keeps that population's effect.
awk 'BEGIN{OFS="\t"}
     NR==FNR {a1[$1]=$2; b[$1]=$3; next}
     { if ($1 in a1) { if ($2==a1[$1]) b[$1]+=$3; else b[$1]-=$3 }
       else { a1[$1]=$2; b[$1]=$3 } }
     END { for (s in b) print s, a1[s], b[s] }' \
     "$OUTD/w_EUR.txt" "$OUTD/w_EAS.txt" > "$OUTD/w_SUM.txt"
echo "SUM: $(wc -l < "$OUTD/w_SUM.txt") SNPs"

for v in EUR EAS SUM; do
  d="$OUTD/csx_$v"; mkdir -p "$d"
  "$W/bin/plink" --bfile "$PRS_GENO" --score "$OUTD/w_${v}.txt" 1 2 3 sum \
    --out "$d/score_SCZcsx${v}_auto" >/dev/null 2>&1 || true
  sd=$(awk 'NR>1{s+=$NF; ss+=$NF*$NF; n++} END{m=s/n; printf "%.4g", sqrt(ss/n-m*m)}' "$d/score_SCZcsx${v}_auto.profile" 2>/dev/null)
  echo "  score sd ($v) = $sd"
  "$RSCRIPT" "$REPO/tools/prs_assoc.R" --prs-dir "$d" --pheno "$PHENO" \
    --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
    --eur-ids "$EUR_IDS" --out "$OUTD/prs_association_SCZcsx${v}.tsv" || true
  "$RSCRIPT" "$RSCRIPT_DIR/06_prs_family.R" --prs-dir "$d" --pheno "$PHENO" \
    --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
    --eur-ids "$EUR_IDS" --out "$OUTD/prs_withinfamily_SCZcsx${v}.tsv" || true
done
echo "=== done  $(date) ==="
