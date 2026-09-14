#!/bin/bash
#SBATCH --job-name=v2_rgmatch
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#
# Genetic correlation with an ANCESTRY-MATCHED SCZ GWAS on both sides.
#
# WHY: v1 8.16.3 labels its EUR arm "properly specified", but that was true only
# of the ABCD side.  The disorder sumstats it used are
# PGC3_SCZ_wave3.PRIMARY -- the multi-ancestry release (74,776 cases, European
# plus East Asian).  The European-only release, PGC3_SCZ_wave3.CORE
# (55,085 cases / 78,957 controls), was already on disk and unused.  So no arm
# of v1 was ancestry-matched on both sides, and the same was true of this run
# until now.
#
# N convention: the core file carries NEFF where primary carries NEFFDIV2, and
# its magnitude (~60k against a theoretical Neff of ~130k for these case/control
# counts) says it is the same half-effective-N quantity under a different name.
# It is therefore doubled, matching prepare_sumstats.sh's SCZ handling.  rg is
# insensitive to a constant N scaling; the disorder h2 is not, so read that with
# the caveat.
#
# Runs the same rg for the v1 EUR sumstats too, so the v1/v2 comparison is made
# against an identical, matched disorder side.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
SRC=/rds/user/rajd2/hpc-work/magma/gwas
OUTD="$OUT_V2/ldsc_eur_matched"
mkdir -p "$OUTD"

SCZ_CORE="$OUTD/SCZ_core.tsv"
if [[ ! -s "$SCZ_CORE" ]]; then
  echo "=== prepare EUR-only SCZ (core)  $(date) ==="
  grep -v '^##' "$SRC/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv" \
    | awk 'BEGIN{FS=OFS="\t"}
           NR==1{for(i=1;i<=NF;i++) c[$i]=i; if(!c["NEFF"]){print "no NEFF" >"/dev/stderr"; exit 1}
                 print $0, "N"; next}
           {print $0, 2*$c["NEFF"]}' > "$SCZ_CORE"
  echo "rows: $(( $(wc -l < "$SCZ_CORE") - 1 ))"
fi

if [[ ! -f "$OUTD/SCZ_core.sumstats.gz" ]]; then
  "$W/bin/munge_sumstats.py" --sumstats "$SCZ_CORE" --snp ID --p PVAL \
    --a1 A1 --a2 A2 --signed-sumstats BETA,0 --N-col N \
    --merge-alleles "$W/ldsc_ref/eur_w_ld_chr/w_hm3.snplist" \
    --chunksize 500000 --out "$OUTD/SCZ_core"
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
  [[ -n "$list" ]] || { echo "SKIP $arm: no munged sumstats in $SS"; continue; }
  echo "=== rg SCZ-core x $arm EUR  $(date) ==="
  "$W/bin/ldsc.py" --rg "$OUTD/SCZ_core.sumstats.gz${list}" \
    --ref-ld-chr "$LD" --w-ld-chr "$LD" --out "$OUTD/rg_SCZcore_${arm}" || true
done
echo "=== done  $(date) ==="
