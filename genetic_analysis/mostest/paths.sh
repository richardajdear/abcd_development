# Step 16 (MOSTest discovery arm) shared paths.  Sourced by every
# genetic_analysis/mostest/*.sbatch; never run directly.  Submit every script
# from the REPO ROOT (sbatch copies the script to a spool dir, so paths are
# resolved from SLURM_SUBMIT_DIR, not from $0).
#
# Every value is ${VAR:-default}: export a variable (or set it in
# genetic_analysis/config.local.sh) to relocate it.  The defaults are the
# CSD3 locations used by steps 14-15 (verified 2026-09).

REPO="${REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "$REPO/genetic_analysis/config.sh" ]] \
  || { echo "FATAL: submit from the abcd_development repo root (REPO=$REPO)" >&2; exit 2; }
# config.sh gives GENO_ARRAY (array genotypes, 11,670 x 515k, hg19) and the
# config.local.sh overrides; it sets -euo pipefail.
source "$REPO/genetic_analysis/config.sh"

# Output root (DK atlas).  table_*.tsv under it are tracked; everything else
# (pheno/, geno/, regenie/, zmat/, sumstats, genes.raw) is gitignored.
M16="${M16:-$REPO/genetic_analysis/work/results_70tab/mostest}"
PH="$M16/pheno"; GENO_DIR="$M16/geno"; RG="$M16/regenie"; ZM="$M16/zmat"
SS="$M16/sumstats"; MG="$M16/magma"; LOGS="$M16/logs"

# Inputs
DK_RUN="${DK_RUN:-$REPO/out/thickness_dsk_70_139406217085}"               # per-region LMM fits (step 1)
REF_1LMM="${REF_1LMM:-$REPO/genetic_analysis/work/results_70tab/prs_final_1lmm/pheno}"  # 8,596-child export
GENO_IMP_PRS="${GENO_IMP_PRS:-$REPO/genetic_analysis/work/inputs/geno/abcd_imp_prs}"    # 7.07 M SNPs, all chr
REF_ABCD="${REF_ABCD:-$REPO/genetic_analysis/work/magma_scz2025/ref_abcd}"             # step-14 MAGMA LD ref + annot
GENESIS_1LMM="${GENESIS_1LMM:-$REPO/genetic_analysis/work/results_70tab/scan_1lmm/assoc}"  # pooled GENESIS scan (engine check)
SETS_SCZ="${SETS_SCZ:-$REPO/genetic_analysis/work/magma_scz2025/genesets_scz2025.txt}"
SETS_MDD="${SETS_MDD:-$REPO/legacy/hpc/work/genesets/mdd_highconf.txt}"
SETS_TOP="${SETS_TOP:-$REPO/genetic_analysis/magma_gene_sets/genesets_topgenes.txt}"
SETS_WES="${SETS_WES:-$REPO/genetic_analysis/magma_gene_sets/genesets_wes.txt}"

# Binaries / environments
PLINK1="${PLINK1:-$REPO/legacy/hpc/work/bin/plink}"          # PLINK 1.9 (--clump)
MAGMA="${MAGMA:-$REPO/legacy/hpc/work/bin/magma}"
PY="${PY:-$REPO/legacy/hpc/work/envs/abcd/bin/python}"      # numpy/scipy/pandas
REGENIE="${REGENIE:-$HOME/rds/hpc-work/envs/regenie/bin/regenie}"   # built from envs/regenie_env.yml

# The REGENIE runs.  <family>[_perm]; the order fixes the array indices of
# 02_regenie_step1 (1-7) and 03_regenie_step2 ((run-1)*22 + chr, 1-154).
RUNS=(slope slope_perm ct ct_perm slopeols slopeols_perm global)
# MOSTest families (03/04): the 68-measure families.
FAMILIES=(slope ct slopeols)

# Rule 10: TMPDIR on rds, never node-local /tmp.
export TMPDIR="$M16/tmp"
mkdir -p "$PH" "$GENO_DIR" "$RG" "$ZM" "$SS" "$MG" "$LOGS" "$TMPDIR"

run_family() { local r="$1"; echo "${r%_perm}"; }
run_is_perm() { [[ "$1" == *_perm ]]; }
