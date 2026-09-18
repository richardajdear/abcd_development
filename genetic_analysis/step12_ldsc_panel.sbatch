#!/bin/bash
#SBATCH --job-name=ldsc_panel
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm/%x_%A.log
#
# Step 12: LDSC SNP-heritability and genetic correlation, full disorder panel,
# both phenotype constructions, both atlases.  EUR-arm phenotype GWAS only
# (LDSC's EUR LD scores must match the GWAS sample).
#
# METHOD CHOICES, each a change from legacy 05_ldsc_rg.sbatch (README §8 2026-09-18):
#  * Disorder inputs are the normalised .ma files (SNP A1 A2 freq b se p N) --
#    one format, explicit alleles for --merge-alleles, per-SNP effective N.
#  * EUROPEAN discovery files only: SCZ25_EUR (2025), SCZ_eur (PGC3 european),
#    MDD_eur.  Step 7 correlated the ABCD EUR GWAS with the PGC3 *primary* and
#    MDD *div* (multi-ancestry) files against EUR LD scores; that is the same
#    LD mismatch rule 4 forbids for PRS, and it is dropped here.
#  * ALZ enters WITHOUT the APOE region (ALZ_noAPOE, ALZ_IGAP_noAPOE): one
#    locus of that size breaks LDSC's polygenic assumption and is excluded by
#    convention in AD LDSC work.
#  * h2 is read from the rg log's per-phenotype block (identical to --h2).
#  * rg is only read where the phenotype's own h2 z >= 4 (rule 13); the table
#    carries the flag.
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
W="$REPO/legacy/hpc/work"; LDSC="$W/bin/ldsc.py"; MUNGE="$W/bin/munge_sumstats.py"
LD="$W/ldsc_ref/eur_w_ld_chr/"; HM3="$W/ldsc_ref/eur_w_ld_chr/w_hm3.snplist"
GW25="$REPO/genetic_analysis/work/scores_scz2025/gwas"; GW3="$REPO/legacy/hpc_v2/work/results_v2/prs_final/gwas"
OUT="${LDSC_PANEL_DIR:-$REPO/genetic_analysis/work/ldsc_panel}"; mkdir -p "$OUT/munged"
export TMPDIR="$OUT/tmp"; mkdir -p "$TMPDIR"
[[ -s "$HM3" ]] || { echo "FATAL: $HM3 missing" >&2; exit 2; }

munge_ma() {  # <name> <.ma>
  local out="$OUT/munged/$1"; [[ -s "$out.sumstats.gz" ]] && { echo "  reuse $1"; return; }
  "$MUNGE" --sumstats "$2" --snp SNP --a1 A1 --a2 A2 --p p --signed-sumstats b,0 --N-col N \
           --merge-alleles "$HM3" --chunksize 500000 --out "$out" > /dev/null 2>&1 || { echo "FATAL: munge $1" >&2; exit 3; }
  echo "  munged $1: $(zcat "$out.sumstats.gz" | wc -l) rows"
}
munge_abcd() {  # <name> <sumstats.tsv.gz>: SNP CHR POS A1 A2 AF1 N BETA SE P
  local out="$OUT/munged/$1"; [[ -s "$out.sumstats.gz" ]] && { echo "  reuse $1"; return; }
  local tmp="$OUT/munged/$1.in.tsv"; zcat "$2" > "$tmp"
  "$MUNGE" --sumstats "$tmp" --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 --N-col N \
           --merge-alleles "$HM3" --chunksize 500000 --out "$out" > /dev/null 2>&1 || { echo "FATAL: munge $1" >&2; exit 3; }
  rm -f "$tmp"; echo "  munged $1: $(zcat "$out.sumstats.gz" | wc -l) rows"
}

echo "=== munge disorders  $(date)"
munge_ma SCZ25_EUR "$GW25/SCZ25_EUR.ma";  munge_ma PGC3_EUR "$GW3/SCZ_eur.ma"
munge_ma MDD_EUR "$GW3/MDD_eur.ma";       munge_ma ASD "$GW3/ASD.ma"
munge_ma ALZ_noAPOE "$GW3/ALZ_noAPOE.ma"; munge_ma ALZ_IGAP_noAPOE "$GW3/ALZ_IGAP_noAPOE.ma"
munge_ma EA "$GW3/EA.ma"
DIS=(SCZ25_EUR PGC3_EUR MDD_EUR ASD ALZ_noAPOE ALZ_IGAP_noAPOE EA)

echo "=== munge ABCD phenotypes (EUR arm)  $(date)"
PH=()
for parc in dk hcp; do
  s=""; [[ $parc == hcp ]] && s="_hcp"; R="$REPO/genetic_analysis/work/results_70tab$s"
  for p in global_slope baseline_thickness; do
    munge_abcd "${parc}__${p}__perregion" "$R/assoc_eur/$p.sumstats.tsv.gz";            PH+=("${parc}__${p}__perregion")
    munge_abcd "${parc}__${p}__1lmm"      "$R/scan_1lmm/assoc_eur/${p}_1lmm.sumstats.tsv.gz"; PH+=("${parc}__${p}__1lmm")
  done
done

echo "=== rg: each disorder against the 8 ABCD phenotype files  $(date)"
list=""; for p in "${PH[@]}"; do list="$list,$OUT/munged/$p.sumstats.gz"; done
for D in "${DIS[@]}"; do
  [[ -s "$OUT/rg_$D.log" ]] && grep -q "Summary of Genetic Correlation" "$OUT/rg_$D.log" && { echo "  reuse rg_$D"; continue; }
  "$LDSC" --rg "$OUT/munged/$D.sumstats.gz$list" --ref-ld-chr "$LD" --w-ld-chr "$LD" --out "$OUT/rg_$D" > /dev/null 2>&1 || echo "  WARN: rg $D failed"
  echo "  rg_$D done"
done
echo "=== collect  $(date)"
"$REPO/legacy/hpc/work/envs/abcd/bin/python" "$REPO/genetic_analysis/step12_ldsc_panel_collect.py"
echo "=== step 12 done  $(date)"
