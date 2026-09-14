#!/bin/bash
#SBATCH --job-name=v2_prscsprep
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#
# Build the four PRS-CS input files ONCE, before the 88-task array runs.
#
# WHY SEPARATE: all 22 chromosome tasks for a disorder need the same input file.
# When the array built it itself they raced on one path, and each task read a
# partially-written file -- the SCZ run reported 391 SNPs instead of 7.5 million
# and would have produced a meaningless score without ever failing.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
OUTD="$OUT_V2/prscs"
mkdir -p "$OUTD"
for D in SCZ MDD ALZ ASD; do
echo "=== prepare $D  $(date) ==="
SS="$OUTD/${D}_prscs_input.txt"
# Input hygiene, learned the hard way on the first attempt:
#  - one all-tab blank row in SCZ_core.tsv made PRS-CS's parser (which splits on
#    whitespace and indexes ll[1]) raise IndexError on every chromosome;
#  - PGC-ALZ2 contains rs283815 with z = inf, which silently zeroed an entire
#    PLINK score (50,440 predictors loaded, every SCORESUM exactly 0).
# So every row is now required to have 6 non-empty fields, ATGC alleles, a
# finite effect and a p-value in (0,1].
clean () {  # clean <SNP> <A1> <A2> <EFFECT> <P>
  awk 'BEGIN{FS=OFS="\t"; print "SNP","A1","A2","BETA","P"}
       {  if (NF < 5) next
          snp=$1; a1=toupper($2); a2=toupper($3); b=$4; p=$5
          if (snp=="" || b=="" || p=="") next
          if (a1!~/^[ACGT]$/ || a2!~/^[ACGT]$/) next
          if (b+0 != b+0 || b=="inf" || b=="-inf" || b=="nan" || b=="NA") next
          if (b+0 > 1e6 || b+0 < -1e6) next
          if (p+0 <= 0 || p+0 > 1) next
          print snp,a1,a2,b,p }'
}
# Always regenerate.  The n_gwas assignments live in this case block, so a
# "skip if the file exists" guard would leave n_gwas unset whenever a stale
# input survived -- which silently fell back to 100,000 for MDD and ASD on the
# first attempt, mis-scaling PRS-CS's shrinkage.
case $D in
    SCZ)
      src="$OUT_V2/ldsc_eur_matched/SCZ_core.tsv"
      awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
           {print $c["ID"], $c["A1"], $c["A2"], $c["BETA"], $c["PVAL"]}' "$src" | clean > "$SS"
      NGWAS=130000 ;;
    MDD)
      src="$OUT_V2/mdd_eur/MDD_EUR.tsv"
      awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
           {print $c["rsid"], $c["effect_allele"], $c["other_allele"], $c["beta"], $c["p_value"]}' \
        "$src" | clean > "$SS"
      NGWAS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} FNR==2{print int($c["N"]); exit}' "$src") ;;
    ASD)
      awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2,$3,$4,$5}' "$OUT_V2/control_ASD/ASD.tsv" | clean > "$SS"
      NGWAS=$(awk -F'\t' 'NR==2{print int($6); exit}' "$OUT_V2/control_ASD/ASD.tsv") ;;
    ALZ)
      # PGC-ALZ2 gives Z, not BETA, and no allele frequency.  Feeding Z to
      # PRS-CS is worse than it is for P+T: PRS-CS models beta_hat ~ N(beta,
      # sigma^2/n), and Z = beta/se rescales every SNP by 1/se, which varies
      # with MAF.  Reconstruct beta with the standard approximation
      #   beta ~= z / sqrt(2 p (1-p) (n + z^2))
      # using MAF from the PRS-CS 1000G EUR reference itself, so the control is
      # weighted on the same scale as the other three disorders.
      awk 'BEGIN{FS=OFS="\t"}
           NR==FNR { if(FNR>1) maf[$2]=$6; next }
           FNR==1 { next }
           ($1 in maf) {
             z=$4+0; n=$6+0; p=maf[$1]+0
             if (p<=0 || p>=1 || n<=0) next
             denom = sqrt(2*p*(1-p)*(n + z*z))
             if (denom<=0) next
             print $1,$2,$3,z/denom,$5 }' \
        /rds/user/rajd2/hpc-work/prscs/ldblk_1kg_eur/snpinfo_1kg_hm3 \
        "$OUT_V2/control_ALZ/ALZ.tsv" | clean > "$SS"
      NGWAS=$(awk -F'\t' 'NR==2{print int($6); exit}' "$OUT_V2/control_ALZ/ALZ.tsv") ;;
esac

[[ -n "${NGWAS:-}" && "$NGWAS" -gt 1000 ]] || NGWAS=100000
echo "$NGWAS" > "$OUTD/${D}_n_gwas.txt"
echo "$D: $(( $(wc -l < "$SS") - 1 )) SNPs, n_gwas=$NGWAS"
unset NGWAS
done
echo "=== all inputs prepared  $(date) ==="
