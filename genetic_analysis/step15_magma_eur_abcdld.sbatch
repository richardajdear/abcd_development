#!/bin/bash
#SBATCH --job-name=magma_eur_abcdld
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-4
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 15a: EUR-arm phenotype MAGMA with the ABCD EUR analysis sample as its
# own LD reference -- the same design step 14 used for the pooled arm -- so the
# two arms differ only in sample, not in reference panel (step 13 used 1000G
# EUR, n = 503).  Reference = step-14 abcd_analysis (8,596) restricted to
# eur_anchor.keep; annotation = step-14 ABCD-SNP annotation (reused).
# Cells: {dk,hcp} x {global_slope_1lmm, baseline_thickness_1lmm}
#   -> results_70tab*/magma_panel_abcdld/genes/<pheno>.genes.{raw,out}
# Prerequisite: step14_magma_pooled_prep.sbatch has run.
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
PLINK="$REPO/legacy/hpc/work/bin/plink"; MAGMA="$REPO/legacy/hpc/work/bin/magma"
R14="$REPO/genetic_analysis/work/magma_scz2025/ref_abcd"
REF="$REPO/genetic_analysis/work/magma_scz2025/ref_abcd_eur"; mkdir -p "$REF"
EURKEEP="$REPO/legacy/hpc/work/results/ancestry/eur_anchor.keep"
export TMPDIR="$REPO/genetic_analysis/work/magma_scz2025/tmp"; mkdir -p "$TMPDIR"
[[ -s "$R14/abcd_analysis.bed" ]] || { echo "FATAL: run step14_magma_pooled_prep first" >&2; exit 2; }

# build the EUR reference once; the other array tasks wait on the lock
(
  flock 9
  if [[ ! -s "$REF/abcd_eur.bed" ]]; then
    "$PLINK" --bfile "$R14/abcd_analysis" --keep "$EURKEEP" --make-bed \
      --threads 8 --memory 40000 --out "$REF/abcd_eur" | grep -E "people|variants" | head -3
  fi
) 9>"$REF/.lock"
echo "reference: $(wc -l < "$REF/abcd_eur.fam") individuals"

CELLS=()
for parc in dk hcp; do s=""; [[ $parc == hcp ]] && s="_hcp"
  for p in global_slope_1lmm baseline_thickness_1lmm; do CELLS+=("$s:$p"); done; done
SPEC="${CELLS[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}"; S="${SPEC%%:*}"; P="${SPEC#*:}"
R="$REPO/genetic_analysis/work/results_70tab$S"; OUT="$R/magma_panel_abcdld/genes/$P"
SS="$R/scan_1lmm/assoc_eur/$P.sumstats.tsv.gz"
[[ -s "$SS" ]] || { echo "FATAL: $SS missing" >&2; exit 3; }
mkdir -p "$(dirname "$OUT")"
[[ -s "$OUT.genes.raw" ]] && { echo "exists $OUT"; exit 0; }
zcat "$SS" | awk 'BEGIN{OFS="\t"} NR==1{print "SNP","P","N"; next}{print $1,$10,$7}' > "$OUT.snp_p.tsv"
"$MAGMA" --bfile "$REF/abcd_eur" --gene-annot "$R14/ncbi37.window35-10.abcd.genes.annot" \
  --pval "$OUT.snp_p.tsv" use=SNP,P ncol=N --out "$OUT" | tail -3
rm -f "$OUT.snp_p.tsv"
echo "genes: $(( $(wc -l < "$OUT.genes.out") - 1 ))  $(date)"
