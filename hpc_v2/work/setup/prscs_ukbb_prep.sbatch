#!/bin/bash
#SBATCH --job-name=v2_csukprep
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#
# Inputs for the PRS-CS re-run against the UK BIOBANK LD reference.
#
# WHY RE-RUN.  PRS-CS was the one method that did not reproduce the SCZ ->
# global_slope association (p=0.195 against C+T 4.3e-03 and SBayesR 7.3e-04).
# The most likely reason is its LD reference: 1000 Genomes EUR is 503
# individuals, while SBayesR's shrunk-sparse panel is 50,000 unrelated UKB
# Europeans.  At ~1M SNPs, 503 individuals estimate LD poorly and posterior
# effects attenuate.  PRS-CS ships a UKB panel (ldblk_ukbb_eur, ~375k
# individuals) for exactly this reason, so the fix is to use it.
#
# The 1000G run is KEPT, in prscs/.  Holding everything else fixed and changing
# only the LD panel isolates the LD-reference effect -- which is the actual
# question, and is more informative than silently replacing the old numbers.
#
# TWO OTHER THINGS ALIGNED, so the three methods differ only in their prior:
#  1. SCZ now uses PGC3_SCZ_wave3.EUROPEAN, not .core.  The 1000G PRS-CS run
#     used core (EUR+EAS) -- see 13.7.  SBayesR already used .european.
#  2. ALZ's beta/se reconstruction takes MAF from the UKB panel's own snpinfo
#     rather than 1000G's, so the frequencies match the LD reference in use.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
G=/rds/user/rajd2/hpc-work/magma/gwas
REF=/rds/user/rajd2/hpc-work/prscs/ldblk_ukbb_eur
OUTD="$OUT_V2/prscs_ukbb"
mkdir -p "$OUTD"
SNPINFO=$(ls "$REF"/snpinfo_* 2>/dev/null | head -1)
[[ -s "$SNPINFO" ]] || { echo "FATAL: no snpinfo in $REF" >&2; exit 3; }
echo "UKB snpinfo: $SNPINFO ($(wc -l < "$SNPINFO") rows)"
head -1 "$SNPINFO"

clean () {
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

for D in SCZ MDD ASD ALZ; do
  SS="$OUTD/${D}_prscs_input.txt"
  case $D in
    SCZ)
      src="$G/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
      zcat "$src" | grep -v '^##' \
        | awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
               {print $c["ID"], $c["A1"], $c["A2"], $c["BETA"], $c["PVAL"]}' | clean > "$SS"
      NGWAS=$(zcat "$src" | grep -v '^##' | awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} NR==2{v=2*$c["NEFF"]} END{printf "%d", v}') ;;
    MDD)
      src="$OUT_V2/mdd_eur/MDD_EUR.tsv"
      awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
           {print $c["rsid"], $c["effect_allele"], $c["other_allele"], $c["beta"], $c["p_value"]}' "$src" | clean > "$SS"
      NGWAS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} NR==2{printf "%d", $c["N"]; exit}' "$src") ;;
    ASD)
      awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2,$3,$4,$5}' "$OUT_V2/control_ASD/ASD.tsv" | clean > "$SS"
      NGWAS=$(awk -F'\t' 'NR==2{printf "%d", $6; exit}' "$OUT_V2/control_ASD/ASD.tsv") ;;
    ALZ)
      # Z -> beta, with MAF from the UKB panel so frequency and LD agree.
      awk 'BEGIN{FS=OFS="\t"}
           NR==FNR { if(FNR>1) maf[$2]=$6; next }
           FNR==1 { next }
           ($1 in maf) { z=$4+0; n=$6+0; p=maf[$1]+0
             if (p<=0 || p>=1 || n<=0) next
             d=sqrt(2*p*(1-p)*(n + z*z)); if (d<=0) next
             print $1,$2,$3,z/d,$5 }' \
        "$SNPINFO" "$OUT_V2/control_ALZ/ALZ.tsv" | clean > "$SS"
      NGWAS=$(awk -F'\t' 'NR==2{printf "%d", $6; exit}' "$OUT_V2/control_ALZ/ALZ.tsv") ;;
  esac
  [[ -n "${NGWAS:-}" && "$NGWAS" -gt 1000 ]] || NGWAS=100000
  echo "$NGWAS" > "$OUTD/${D}_n_gwas.txt"
  echo "$D: $(( $(wc -l < "$SS") - 1 )) SNPs, n_gwas=$NGWAS"
  unset NGWAS
done
echo "=== done  $(date) ==="
