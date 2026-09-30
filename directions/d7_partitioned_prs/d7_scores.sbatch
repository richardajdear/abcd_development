#!/bin/bash
#SBATCH --job-name=d7_scores
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=10:00:00
#SBATCH --output=directions/d7_partitioned_prs/work/logs/%x_%A.log
#
# D7 steps 1-2: assign SNPs to genes, then one pass over the genotypes that
# writes per-gene partial SBayesRC scores for every arm in config.sh (plus any
# extra arm strings given as arguments, e.g. the step-4 annotated refits).
# Run from the repo root:  sbatch directions/d7_partitioned_prs/d7_scores.sbatch
# Memory: 2 imputation modes x arms x 11.7k x ~19.5k float64 (~18 GB for 5 arms)
# plus the genotype chunk; 192G leaves headroom for 8 arms.
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/directions/d7_partitioned_prs/config.sh
cd "$REPO"; mkdir -p "$WORK/logs"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
[[ -s "$GENELOC" ]] || { echo "FATAL: $GENELOC missing" >&2; exit 2; }
[[ -s "$GENO.bed" ]] || { echo "FATAL: $GENO.bed missing" >&2; exit 2; }

ARGS=()
for spec in "${ARMS[@]}" "$@"; do
  name=${spec%%|*}; r=${spec#*|}; wf=${r%%|*}; r=${r#*|}; prof=${r%%|*}; cell=${r#*|}   # no here-strings (rule 10)
  if [[ -s "$wf" && -s "$prof" ]]; then ARGS+=(--arm "$spec"); echo "arm $name ($cell): $wf"
  else echo "WARN: arm $name skipped, missing $wf or $prof" >&2; fi
done
(( ${#ARGS[@]} )) || { echo "FATAL: no arm has both weights and a profile" >&2; exit 3; }

echo "=== step 1: SNP -> gene  $(date +%T)"
[[ -s "$WORK/snp2gene.tsv.gz" ]] || "$D7PY" "$D7/code/01_assign_snps.py" --bim "$GENO.bim" --gene-loc "$GENELOC" --out "$WORK"
echo "=== step 2: per-gene scores  $(date +%T)"
"$D7PY" "$D7/code/02_gene_scores.py" --bfile "$GENO" --work "$WORK" --threads "${SLURM_CPUS_PER_TASK:-8}" "${ARGS[@]}"
cp "$WORK/gene_scores/gate.tsv" "$RESULTS/table_d7_score_gate.tsv"
echo "=== done $(date +%T)"
