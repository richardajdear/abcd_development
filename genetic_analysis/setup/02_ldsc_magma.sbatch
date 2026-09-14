#!/bin/bash
#SBATCH --job-name=v3_down
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#
# Experiment A, steps 5 and 7 -- LDSC and MAGMA on the EUR-arm sumstats for
# the four new phenotypes.
#
# Same construction as genetic_analysis/setup/eur_magma_ldsc.sbatch: it runs v1's
# OWN 04_magma.sbatch and 05_ldsc_rg.sbatch with the directory variables
# redirected, so any difference from the published v2 numbers is the GWAS and
# not the downstream code.
#
# EUR arm, not pooled, because that is where the published h2/rg panel lives
# (legacy/hpc_v2/README_HPC.md 12.4) and legacy/hpc_v3/h2_vs_scz_relevance.csv --- which
# Experiment A adds four rows to --- is built from ldsc_eur.
#
# TWO DELIBERATE ASYMMETRIES:
#
# * LDSC runs on all NINE phenotypes.  Munging is ~1 min each, and putting the
#   settled five through the same call reproduces the published rg table in
#   the same file as the new rows -- a reproduction check that costs minutes.
# * MAGMA runs on the FOUR new ones only.  Its gene analysis is the expensive
#   step (~1.5 h for five phenotypes) and the settled five are already done in
#   magma_eur/; nothing about them would change.
#
# Fresh output directories.  Writing into ldsc_eur/ and magma_eur/ would have
# overwritten ldsc_rg_summary.tsv and magma_gsa_summary.tsv with tables
# covering only the phenotypes in THIS manifest, silently deleting published
# v2 results.  LDSC's munged disorder sumstats ARE symlinked in, because
# munge_one explicitly reuses an existing .sumstats.gz; MAGMA's .genes.raw is
# deliberately not -- see the note beside the symlink loop below.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/genetic_analysis/config.sh"

STAGE="$OUT_V2/gwas_eur_staged"
LDSC_V3="$OUT_V2/ldsc_eur_v3"
MAGMA_V3="$OUT_V2/magma_eur_v3"
PDIR="$REPO/legacy/hpc_v2/work/pheno_v3"
ensure_dirs "$STAGE" "$LDSC_V3" "$MAGMA_V3"

# --- stage the four new EUR sumstats in the .fastGWA layout both tools read --
for pheno in slope_topDelta slope_topC3 slope_projDelta slope_projC3; do
  src="$OUT_V2/assoc_eur/$pheno.sumstats.tsv.gz"
  dst="$STAGE/$pheno.fastGWA"
  [[ -f "$src" ]] || { echo "FATAL: $src missing -- step 3e incomplete" >&2; exit 3; }
  if [[ ! -s "$dst" ]]; then
    echo "=== stage $pheno  $(date) ==="
    # v2 sumstats columns: SNP CHR POS A1 A2 AF1 N BETA SE P
    zcat "$src" | awk 'BEGIN{OFS="\t"} NR==1{print "SNP","A1","A2","AF1","N","BETA","SE","P"; next}
                       {print $1,$4,$5,$6,$7,$8,$9,$10}' > "$dst"
  fi
  echo "  $pheno: $(( $(wc -l < "$dst") - 1 )) SNPs"
done

# --- reuse the disorder-side work --------------------------------------------
for f in "$OUT_V2"/ldsc_eur/SCZ.sumstats.gz "$OUT_V2"/ldsc_eur/MDD.sumstats.gz; do
  [[ -f "$f" ]] && ln -sf "$f" "$LDSC_V3/$(basename "$f")"
done
# NOT symlinked for MAGMA.  hpc/04_magma.sbatch step 3 re-runs the disorder
# gene analysis unconditionally with --out "$MAGMA_DIR/$D", and MAGMA writing
# to a path that is a symlink would follow it and overwrite magma_eur/'s
# published SCZ.genes.raw.  It regenerates them inside magma_eur_v3/ instead:
# a few minutes of compute to keep v2's outputs untouchable.

export GWAS_DIR="$STAGE"
export DRY_RUN=0

echo "=== LDSC, all nine phenotypes  $(date) ==="
export LDSC_DIR="$LDSC_V3"
export MANIFEST="$PDIR/phenotype_manifest.tsv"
bash "$REPO/legacy/hpc/05_ldsc_rg.sbatch"

echo "=== MAGMA, the four new phenotypes  $(date) ==="
export MAGMA_DIR="$MAGMA_V3"
export MANIFEST="$PDIR/phenotype_manifest_new_only.tsv"
bash "$REPO/legacy/hpc/04_magma.sbatch"

for f in "$LDSC_V3/ldsc_rg_summary.tsv"; do
  [[ -s "$f" ]] || { echo "FATAL: $f missing/empty" >&2; exit 4; }
done
echo "=== done  $(date) ==="
