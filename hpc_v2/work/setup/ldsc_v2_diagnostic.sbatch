#!/bin/bash
#SBATCH --job-name=v2_ldsc
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#
# LDSC feasibility diagnostic on the v2 (GENESIS) sumstats.
#
# This does NOT produce a reportable rg.  The v2 scan is pooled multi-ancestry
# and both disorder GWAS are European, so no single set of LD scores is correct
# for both sides of a bivariate regression (hpc/README_HPC.md 8.16.6).  What it
# answers is two things that ARE well posed:
#
#   1. the univariate h2 z-score, which is the gate on whether an rg is
#      estimable at all (LDSC needs h2 z ~ 4; v1's EUR arm gave 3.89 for
#      baseline_thickness and 0.57 for global_slope);
#   2. whether the LDSC intercept still sits at ~1.10 with EUR LD scores now
#      that PC-AiR/PC-Relate demonstrably control structure better
#      (v2 lambda_GC 1.036 vs v1 1.107).  v1 8.16.2 attributed that intercept
#      to LD-score misspecification rather than confounding; if it is
#      unchanged under a better-controlled GWAS, that attribution is confirmed.
#
# Both LD-score panels are run: EUR (eur_w_ld_chr) and the v1 in-sample panel.
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/config.sh
V1="$V1_ROOT"
source "$REPO_ROOT/hpc/config.local.sh" 2>/dev/null || true
LDSC_MUNGE="$V1/bin/munge_sumstats.py"
LDSC_BIN="$V1/bin/ldsc.py"
HM3="$V1/ldsc_ref/eur_w_ld_chr/w_hm3.snplist"
LD_EUR="$V1/ldsc_ref/eur_w_ld_chr/"
LD_INS="$V1/results/ld_insample/"
OUTD="$OUT_V2/ldsc"
mkdir -p "$OUTD"

for pheno in $(phenotype_names); do
  f="$ASSOC_DIR/$pheno.sumstats.tsv.gz"
  [[ -f "$f" ]] || { echo "SKIP $pheno: $f missing"; continue; }
  if [[ ! -f "$OUTD/$pheno.sumstats.gz" ]]; then
    echo "=== munge $pheno  $(date) ==="
    "$LDSC_MUNGE" --sumstats "$f" --snp SNP --p P --signed-sumstats BETA,0 \
      --N-col N --merge-alleles "$HM3" --chunksize 500000 --out "$OUTD/$pheno"
  fi
done
cp -n "$V1/results/ldsc_imp/SCZ.sumstats.gz" "$V1/results/ldsc_imp/MDD.sumstats.gz" "$OUTD/" 2>/dev/null || true

for panel in eur insample; do
  case $panel in eur) LD="$LD_EUR";; insample) LD="$LD_INS";; esac
  for pheno in $(phenotype_names); do
    [[ -f "$OUTD/$pheno.sumstats.gz" ]] || continue
    echo "=== h2 $pheno ($panel LD)  $(date) ==="
    "$LDSC_BIN" --h2 "$OUTD/$pheno.sumstats.gz" \
      --ref-ld-chr "$LD" --w-ld-chr "$LD" --out "$OUTD/h2_${panel}_${pheno}" || true
  done
  for D in SCZ MDD; do
    [[ -f "$OUTD/$D.sumstats.gz" ]] || continue
    list=$(for p in $(phenotype_names); do [[ -f "$OUTD/$p.sumstats.gz" ]] && printf '%s,' "$OUTD/$p.sumstats.gz"; done)
    echo "=== rg $D ($panel LD)  $(date) ==="
    "$LDSC_BIN" --rg "$OUTD/$D.sumstats.gz,${list%,}" \
      --ref-ld-chr "$LD" --w-ld-chr "$LD" --out "$OUTD/rg_${panel}_${D}" || true
  done
done
echo "=== done  $(date) ==="
