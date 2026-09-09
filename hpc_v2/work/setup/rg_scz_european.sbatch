#!/bin/bash
#SBATCH --job-name=v2_rgeur
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#
# CORRECTION to section 12.5: redo the "ancestry-matched" SCZ rg with the file
# that is actually European-only.
#
# 12.5 used PGC3_SCZ_wave3.CORE, describing it as the EUR-only release.  The
# release README says otherwise: "The core meta-analysis consists of cohorts
# primarily of east asian and european ancestry."  core is EUR+EAS
# (55,085 cases), and the EUR-only file is PGC3_SCZ_wave3.EUROPEAN
# (53,386 cases / 77,258 controls).  So 12.5 improved on v1's `primary`
# (EUR+EAS+AFR+LAT) but did not reach the matched pair it claimed.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
G=/rds/user/rajd2/hpc-work/magma/gwas
OUTD="$OUT_V2/ldsc_eur_matched"
mkdir -p "$OUTD"

SRC="$OUTD/SCZ_european.tsv"
if [[ ! -s "$SRC" ]]; then
  # N = 2 x NEFF, matching prepare_sumstats.sh's SCZ convention.
  zcat "$G/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz" | grep -v '^##' \
    | awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++)c[$i]=i; print $0,"N"; next}
           {print $0, 2*$c["NEFF"]}' > "$SRC"
  echo "rows: $(( $(wc -l < "$SRC") - 1 ))"
fi
if [[ ! -f "$OUTD/SCZ_european.sumstats.gz" ]]; then
  "$W/bin/munge_sumstats.py" --sumstats "$SRC" --snp ID --p PVAL --a1 A1 --a2 A2 \
    --signed-sumstats BETA,0 --N-col N \
    --merge-alleles "$W/ldsc_ref/eur_w_ld_chr/w_hm3.snplist" \
    --chunksize 500000 --out "$OUTD/SCZ_european"
fi
LD="$W/ldsc_ref/eur_w_ld_chr/"
for arm in v2 v1; do
  case $arm in
    v2) SS="$OUT_V2/ldsc_eur" ;;
    v1) SS="$W/results/ldsc_imp_eur" ;;
  esac
  list=""
  for p in $(phenotype_names); do
    [[ -f "$SS/$p.sumstats.gz" ]] && list="$list,$SS/$p.sumstats.gz"
  done
  "$W/bin/ldsc.py" --rg "$OUTD/SCZ_european.sumstats.gz${list}" \
    --ref-ld-chr "$LD" --w-ld-chr "$LD" --out "$OUTD/rg_SCZeuropean_${arm}" || true
done
echo "=== done  $(date) ==="
