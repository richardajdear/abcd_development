#!/bin/bash
#SBATCH --job-name=c4_impute
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=slurm/%x_%j.log
#
# Steps 1-3 on CSD3: MHC extract of the 7.0 array fileset (all 11,670
# genotyped children; ancestry restriction happens in step 4) -> Beagle C4
# imputation against the Sekar 2016 panel -> per-child C4 copy numbers + GREx.
#
#   cd ~/rds/hpc-work/abcd_development && sbatch c4_imputation/c4_impute.sbatch
#
# Every path below is an overridable environment variable, so the same script
# runs on the laptop as a pipeline test (see README.md, "Smoke tests").
# Outputs: work/mhc.vcf.gz, work/c4_imputed.vcf.gz, work/c4_calls.tsv
# (per-subject: never committed) and results/c4_imputation_summary.tsv (commit).
set -euo pipefail
REPO="${REPO:-/home/rajd2/rds/hpc-work/abcd_development}"; cd "$REPO"
C4="$REPO/c4_imputation"; WORK="${C4_WORK:-$C4/work}"; RES="${C4_RESULTS:-$C4/results}"; mkdir -p "$WORK" "$RES" slurm
export TMPDIR="${TMPDIR_OVERRIDE:-$WORK/tmp}"; mkdir -p "$TMPDIR"     # rule 10: never node-local /tmp
# GENO_ARRAY comes from genetic_analysis/config.sh + config.local.sh (the CSD3 path authority)
if [[ -z "${GENO:-}" && -f genetic_analysis/config.sh ]]; then
  GENO="$(bash -c 'source genetic_analysis/config.sh >/dev/null 2>&1; echo "$GENO_ARRAY"')"; fi
: "${GENO:?set GENO to the 7.0 array PLINK prefix}"
C4_ENV="${C4_ENV:-$HOME/rds/hpc-work/envs/c4}"
export BCFTOOLS="${BCFTOOLS:-$C4_ENV/bin/bcftools}"
JAVA_DEF="$C4_ENV/bin/java"; [[ -x "$JAVA_DEF" ]] || JAVA_DEF="$C4_ENV/lib/jvm/bin/java"
export JAVA="${JAVA:-$JAVA_DEF}"
export PLINK="${PLINK:-$REPO/legacy/hpc/work/bin/plink}"
PY="${PY:-$REPO/legacy/hpc/work/envs/abcd/bin/python}"
EUR_KEEP="${EUR_KEEP:-$REPO/legacy/hpc/work/results/ancestry/eur_anchor.keep}"
export BEAGLE_JAR="$C4/resources/beagle.27Feb25.75f.jar" CONFORM_JAR="$C4/resources/conform-gt.24May16.cee.jar"
export PANEL="$C4/resources/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37.vcf.gz"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
for f in "$BCFTOOLS" "$JAVA" "$PLINK" "$PY" "$EUR_KEEP" "$BEAGLE_JAR" "$CONFORM_JAR" "$PANEL" "$GENO.bed"; do
  [[ -e "$f" ]] || command -v "$f" >/dev/null 2>&1 || { echo "FATAL: missing $f (env: bash c4_imputation/env/setup_csd3_env.sh; resources: bash c4_imputation/00_fetch_resources.sh)" >&2; exit 2; }
done
echo "c4_impute: GENO=$GENO  $(date)"
bash "$C4/01_extract_mhc.sh" "$GENO" "$WORK/mhc"
bash "$C4/02_impute_c4.sh" "$WORK/mhc.vcf.gz" "$WORK/c4_imputed" "$THREADS" 16
"$PY" "$C4/03_c4_grex.py" --vcf "$WORK/c4_imputed.vcf.gz" --out-calls "$WORK/c4_calls.tsv" \
  --out-summary "$RES/c4_imputation_summary.tsv" --keep-ids "$EUR_KEEP" --panel "$PANEL" --bcftools "$BCFTOOLS"
cp "$WORK/c4_imputed.conform.log" "$WORK/conform_marker_log.txt"
# marker fate tally (SUMMARY / INFO columns of the conform-gt log): record it in the README
awk -F'\t' 'NR>1{n[$9" "$10]++} END{for(k in n) print "conform-gt", k, n[k]}' "$WORK/conform_marker_log.txt" | sort | tee "$RES/conform_marker_tally.txt"
echo "done $(date)"
