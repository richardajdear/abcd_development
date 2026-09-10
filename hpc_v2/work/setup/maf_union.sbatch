#!/bin/bash
#SBATCH --job-name=v2_maf
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#
# ASD and ALZ ship no allele frequency, so the normaliser borrows one.  Using
# the PRS-CS UKB HapMap3 snpinfo alone (1,117,425 SNPs) threw away most of the
# ASD file: 378,844 of 3,085,577 SNPs survived.  The frequency we actually want
# for scoring is the frequency in OUR target sample, which covers all 7,068,032
# imputed SNPs -- and a SNP absent from the target cannot be scored anyway.
# So: target MAF first, HapMap3 panel only to fill SNPs the target lacks (those
# still matter to SBayesR/SBayesRC, whose LD matrices are HapMap3).
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/work/setup/paths.sh
"$PLINK" --bfile "$GENO" --freq --out "$IN/ref/target" --memory 28000
HM3="$IN/ref/ldblk_ukbb_eur/snpinfo_ukbb_hm3"
{ printf 'SNP\tMAF\n'
  awk 'NR>1 && $5>0 && $5<1 {print $2"\t"$5}' "$IN/ref/target.frq"
  awk 'NR==FNR{if(FNR>1 && $5>0 && $5<1) k[$2]=1; next}
       FNR>1 && !($2 in k) {print $2"\t"$6}' "$IN/ref/target.frq" "$HM3"
} > "$IN/ref/maf_union.tsv"
echo "target      : $(( $(wc -l < "$IN/ref/target.frq") - 1 ))"
echo "union panel : $(( $(wc -l < "$IN/ref/maf_union.tsv") - 1 ))"
