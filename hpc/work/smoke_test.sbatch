#!/bin/bash
#SBATCH --job-name=abcd_smoke
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --output=/home/rajd2/rds/hpc-work/abcd_development/hpc/work/logs/%x_%j.log
#
# End-to-end smoke test of hpc/ steps 01-06 against the real ABCD genotypes,
# with the association scan restricted to chromosome 22.
#
# Scope of the restriction, stated precisely:
#
#   REML (02)   uses the real genome-wide GRM, so its h2 estimates are REAL
#               numbers on the real unrelated subset -- not fixture numbers.
#   GWAS (03)   scans chr22 only.  Correctness of the fastGWA path does not
#               depend on which SNPs it reads, and every subject is retained so
#               the sparse-GRM relatedness path is genuinely exercised.
#   MAGMA (04)  gene analysis of the ABCD side therefore covers chr22 genes
#               only; the SCZ/MDD side is genome-wide and real.
#   LDSC (05)   an ABCD sumstats file with only chr22 SNPs has far too few HM3
#               SNPs for a defined rg.  NA is the expected result; what is
#               under test is munging, the LD-score read, and the log parser.
#   PRS (06)    clumping is genome-wide in the reference panel and real; the
#               scores are then computed over chr22 SNPs only, so the effect
#               sizes mean nothing.  What is under test is the file plumbing
#               and that prs_assoc.R fits its models.
#
# Every step runs even if an earlier one fails: this is a defect hunt, and
# stopping at the first error would find one defect per submission.
set -uo pipefail

REPO=/home/rajd2/rds/hpc-work/abcd_development
W=$REPO/hpc/work
export GENO=$W/genotype/abcd_eur_chr22
export THREADS=${SLURM_CPUS_PER_TASK:-16}
# Unset so 02/03 loop over every phenotype rather than picking one array index.
unset SLURM_ARRAY_TASK_ID

cd "$REPO"
declare -a STEP_NAMES STEP_CODES

run_step() {                # run_step <label> <command...>
  local label="$1"; shift
  echo
  echo "################################################################"
  echo "### $label   $(date)"
  echo "################################################################"
  "$@"
  local code=$?
  STEP_NAMES+=("$label"); STEP_CODES+=("$code")
  echo "### $label exit=$code"
}

# --- 00 preflight ----------------------------------------------------------
run_step "00_check_inputs" bash hpc/00_check_inputs.sh

# --- 01 GRM ----------------------------------------------------------------
# $GRM and $GRM_SPARSE are the prebuilt read-only genome-wide GRMs, so running
# 01_grm as written would try to overwrite them.  Validate its command
# construction under DRY_RUN, then execute only the piece that is missing:
# the relatedness-pruned GRM that 02_reml needs.
run_step "01_grm (DRY_RUN)" env DRY_RUN=1 bash hpc/01_grm.sbatch

run_step "01_grm --grm-cutoff (real)" bash -c '
  source hpc/config.sh
  ensure_dirs "$(dirname "$GRM_UNREL")"
  if [[ -f "$GRM_UNREL.grm.bin" ]]; then
    echo "reuse $GRM_UNREL"; exit 0
  fi
  "$GCTA" --grm "$GRM" --grm-cutoff "$GRM_CUTOFF" --make-grm \
          --thread-num "$THREADS" --out "$GRM_UNREL"
  echo "unrelated at $GRM_CUTOFF: $(wc -l < "$GRM_UNREL.grm.id") of $(wc -l < "$GRM.grm.id")"
'

# --- 02-06 -----------------------------------------------------------------
run_step "02_reml"  bash hpc/02_reml.sbatch
run_step "03_gwas"  bash hpc/03_gwas.sbatch
run_step "04_magma" bash hpc/04_magma.sbatch
run_step "05_ldsc"  bash hpc/05_ldsc_rg.sbatch
run_step "06_prs"   bash hpc/06_prs.sbatch

# --- summary ---------------------------------------------------------------
echo
echo "################################################################"
echo "### SMOKE TEST SUMMARY   $(date)"
echo "################################################################"
fail=0
for i in "${!STEP_NAMES[@]}"; do
  status=ok; [[ "${STEP_CODES[$i]}" != "0" ]] && { status=FAIL; fail=1; }
  printf '%-28s %-5s exit=%s\n' "${STEP_NAMES[$i]}" "$status" "${STEP_CODES[$i]}"
done
exit $fail
