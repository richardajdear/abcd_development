#!/bin/bash
#SBATCH --job-name=geom_ldsc
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#
# D1 step 9 (optional): LDSC h2 of the EUR-arm geometry slopes, munge and LD
# scores exactly as step12_ldsc_panel.sbatch.  Report only; writes scan_geom/ldsc_eur/.
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
L="$REPO/legacy/hpc/work"; LDSC="$L/bin/ldsc.py"; MUNGE="$L/bin/munge_sumstats.py"
LD="$L/ldsc_ref/eur_w_ld_chr/"; HM3="$L/ldsc_ref/eur_w_ld_chr/w_hm3.snplist"
W="$REPO/genetic_analysis/work/results_70tab_hcp/scan_geom"; OUT="$W/ldsc_eur"; mkdir -p "$OUT/munged"
export TMPDIR="$W/tmp"; mkdir -p "$TMPDIR"
printf 'phenotype\th2\tse\tz\tintercept\n' > "$OUT/ldsc_h2_summary.tsv"
for p in slope_thk slope_area slope_vol slope_thk_area slope_thk_vol_area; do
  out="$OUT/munged/$p"
  if [[ ! -s "$out.sumstats.gz" ]]; then
    tmp="$OUT/munged/$p.in.tsv"; zcat "$W/assoc_eur/$p.sumstats.tsv.gz" > "$tmp"
    "$MUNGE" --sumstats "$tmp" --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 --N-col N \
             --merge-alleles "$HM3" --chunksize 500000 --out "$out" > /dev/null 2>&1 || { echo "FATAL: munge $p" >&2; exit 3; }
    rm -f "$tmp"
  fi
  "$LDSC" --h2 "$out.sumstats.gz" --ref-ld-chr "$LD" --w-ld-chr "$LD" --out "$OUT/h2_$p" > /dev/null 2>&1 || { echo "FATAL: h2 $p" >&2; exit 4; }
  awk -v p="$p" '/Total Observed scale h2:/ {h=$5; se=$6; gsub(/[()]/,"",se)} /^Intercept:/ {i=$2} END {printf "%s\t%s\t%s\t%.2f\t%s\n", p, h, se, h/se, i}' \
      "$OUT/h2_$p.log" >> "$OUT/ldsc_h2_summary.tsv"
done
column -t "$OUT/ldsc_h2_summary.tsv"
