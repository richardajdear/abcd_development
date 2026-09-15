#!/bin/bash
#SBATCH --job-name=ldscmagma70
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=genetic_analysis/work/results_70tab/logs/%x_%j.log
#
# README_HPC.md 3 step 7: LDSC h2 + rg (EUR arm), MAGMA gene + gene-property,
# and the prioritised gene-set panel -- on the 7.0 EUR association sumstats.
#
# Same construction as setup/eur_magma_ldsc.sbatch + setup/eur_prio_gsa.sbatch:
# v1's OWN 04_magma.sbatch / 05_ldsc_rg.sbatch / work/prioritised_gsa.sbatch
# run verbatim with only their directory variables redirected, so any change
# from the published v2 numbers is the GWAS, not the downstream code.  Fresh
# output directories under results_70tab/ (rule 9): the collectors rewrite
# *_summary.tsv to cover exactly the manifest they are handed.
#
# Two migration facts this script depends on (verified 2026-09-15):
#  * legacy/hpc/config.sh derives REPO_ROOT from its own location, which after
#    the 2.2 move is legacy/, and its config.local.sh then sets
#    ABCD_HPC_ROOT="$REPO_ROOT/hpc/work" = legacy/hpc/work -- every v1 asset
#    resolves -- PROVIDED REPO_ROOT is not already set in the environment.
#    genetic_analysis/config.sh does not export it, but we unset it anyway.
#  * legacy/hpc/work/bin/{ldsc.py,munge_sumstats.py} were bash wrappers with
#    the pre-move path compiled in; fixed in place (originals *.pre70bak).
#  * prioritised_gsa.sbatch defaults MAGMA to $REPO/hpc/work/bin/magma, which
#    no longer exists, so MAGMA is passed explicitly.
#
# Rule 13: read h2_z FIRST.  rg for global_slope is uninformative unless its
# h2_z rose from ~1.2 to >= 4, which LDSC on ~4.2k EUR is unlikely to give.
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
source genetic_analysis/config.sh
STAGE="$OUT_V2/gwas_eur_staged"; MAG="$OUT_V2/magma_eur"; LD="$OUT_V2/ldsc_eur"; PRIO="$OUT_V2/magma_prio_eur"
MANI="$PHENO_DIR/phenotype_manifest.tsv"; MAGMA_BIN="$MAGMA"; W="$REPO_ROOT/legacy/hpc/work"
ensure_dirs "$STAGE" "$MAG" "$LD" "$PRIO"

echo "=== stage EUR sumstats -> .fastGWA layout  $(date) ==="
for pheno in $(phenotype_names); do
  src="$OUT_V2/assoc_eur/$pheno.sumstats.tsv.gz"; dst="$STAGE/$pheno.fastGWA"
  [[ -f "$src" ]] || { echo "FATAL: $src missing" >&2; exit 3; }
  if [[ ! -s "$dst" ]]; then
    # v2 sumstats columns: SNP CHR POS A1 A2 AF1 N BETA SE P
    zcat "$src" | awk 'BEGIN{OFS="\t"} NR==1{print "SNP","A1","A2","AF1","N","BETA","SE","P"; next}
                       {print $1,$4,$5,$6,$7,$8,$9,$10}' > "$dst"
  fi
  echo "  $pheno: $(( $(wc -l < "$dst") - 1 )) SNPs"
done

export GWAS_DIR="$STAGE" MAGMA_DIR="$MAG" LDSC_DIR="$LD" MANIFEST="$MANI" DRY_RUN=0
unset REPO_ROOT ABCD_HPC_ROOT     # let v1's config resolve to legacy/hpc/work

echo "=== MAGMA (v1 04_magma.sbatch)  $(date) ==="
bash legacy/hpc/04_magma.sbatch
echo "=== LDSC (v1 05_ldsc_rg.sbatch)  $(date) ==="
bash legacy/hpc/05_ldsc_rg.sbatch

echo "=== prioritised gene-set panel  $(date) ==="
MAGMA="$MAGMA_BIN" MAGMA_DIR="$MAG" OUT="$PRIO" SETS="$W/genesets/scz_prioritised.txt" TAG=sczprio COND_SET=SCZ_locus_pool \
  bash "$W/prioritised_gsa.sbatch"
MAGMA="$MAGMA_BIN" MAGMA_DIR="$MAG" OUT="$PRIO" SETS="$W/genesets/mdd_highconf.txt" TAG=mddhc COND_SET=MDD_pool \
  bash "$W/prioritised_gsa.sbatch"

for f in "$LD/ldsc_rg_summary.tsv" "$MAG/magma_summary.tsv" "$PRIO/sczprio_gsa_summary.tsv" "$PRIO/mddhc_gsa_summary.tsv"; do
  [[ -s "$f" ]] && echo "ok  $f ($(( $(wc -l < "$f") - 1 )) rows)" || { echo "FATAL: $f missing/empty" >&2; exit 4; }
done
echo "=== step 7 done  $(date) ==="
