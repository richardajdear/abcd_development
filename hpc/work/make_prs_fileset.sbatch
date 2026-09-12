#!/bin/bash
#SBATCH --job-name=abcd_prsfs
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --output=%x_%A_%a.log
#
# Build the single merged fileset 06_prs needs, WITHOUT merging all 26M variants.
#
# 06_prs scores with PLINK 1.9 --score, which takes one --bfile and has no
# --mbfile, so unlike fastGWA it cannot read the 22 per-chromosome filesets as a
# list.  Merging them all would be 72 GB of duplicate genotypes.
#
# But a polygenic score can only ever use variants that appear in the discovery
# sumstats, and the clumped index SNPs are a subset of those.  The union of the
# SCZ and MDD rsIDs present in our data is 7,068,032 variants -- 20.6 GB rather
# than 72 GB, and it provably contains every SNP the score can reach.
#
# Array task N extracts chromosome N; task 0 (submitted separately, after) merges.
set -euo pipefail

if [[ -z "${REPO_ROOT:-}" || ! -f "${REPO_ROOT:-}/hpc/config.sh" ]]; then
  REPO_ROOT=""
  for _d in "$(dirname "${BASH_SOURCE[0]}")/../.." "${SLURM_SUBMIT_DIR:-}" "$PWD"; do
    [[ -n "$_d" ]] || continue
    _p="$(cd "$_d" 2>/dev/null && pwd)" || continue
    while [[ -n "$_p" && "$_p" != "/" ]]; do
      if [[ -f "$_p/hpc/config.sh" ]]; then REPO_ROOT="$_p"; break 2; fi
      _p="$(dirname "$_p")"
    done
  done
fi
[[ -n "$REPO_ROOT" ]] || { echo "FATAL: cannot find hpc/config.sh" >&2; exit 2; }
# shellcheck source=/dev/null
source "$REPO_ROOT/hpc/config.sh"

IMP_OUT="${IMP_OUT:-$ABCD_HPC_ROOT/genotype_imputed}"
SNPS="${PRS_SNPS:-$ABCD_HPC_ROOT/results/ancestry/prs_snps.txt}"
D="${PRS_FILESET_DIR:-$IMP_OUT/prs}"
ensure_dirs "$D"
CHR=${SLURM_ARRAY_TASK_ID:?run as an array job}

echo "=== chr$CHR  $(date) ==="
"$PLINK2" --bfile "$IMP_OUT/imp_union_chr${CHR}" \
          --extract "$SNPS" \
          --make-bed --threads "${SLURM_CPUS_PER_TASK:-4}" \
          --out "$D/prs_chr${CHR}"
bash "$REPO_ROOT/hpc/work/check_bfile_integrity.sh" "$D/prs_chr${CHR}" || exit 3
echo "  chr$CHR: $(wc -l < "$D/prs_chr${CHR}.bim") variants"
