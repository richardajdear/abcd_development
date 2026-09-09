#!/bin/bash
#SBATCH --job-name=v3_prs
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#
# Experiment A, step 6 -- PRS association for the four new phenotypes.
#
# C+T ONLY at this stage (instruction of 2026-09-08).  PRS-CS and SBayesR
# exist for the settled five (hpc_v2/README_HPC.md 13.5) and are deliberately
# not run here.
#
# NOTHING IS RE-SCORED.  A polygenic score applies consortium weights to ABCD
# genotypes; it does not involve our phenotypes, so the score_*.profile files
# already on disk are exactly the ones a new phenotype needs.  $PRSD is a
# directory of symlinks gathering the three places C+T scores were written:
#   hpc/work/results/prs_imp                 SCZ, MDD   (v1 step 06)
#   hpc_v2/work/results_v2/control_ASD/prs   ASD        (13.1)
#   hpc_v2/work/results_v2/control_ALZ/prs   ALZ, ALZnoAPOE (13.6, after the
#                                            rs283815 z=inf repair)
#
# The manifest is the FULL nine, not phenotype_manifest_new_only.tsv: the
# primary readout is the new phenotype's beta COMPARED WITH global_slope's in
# the same subjects and the same model, so global_slope must be refitted here
# rather than quoted from a previous run.  Its reproduction of v2's published
# number is the check that this arm is wired correctly.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"

PRSD="$OUT_V2/prs_ct_v3"
OUTD="$REPO/hpc_v3/work/prs"
PHENO="$REPO/hpc_v2/work/pheno_v3/phenotypes_gcta.txt"
QUANT="$REPO/hpc_v2/work/pheno_v3/covar_quant.txt"
CAT="$REPO/hpc_v2/work/pheno_v3/covar_categorical.txt"
MAN="$REPO/hpc_v2/work/pheno_v3/phenotype_manifest.tsv"
# EUR = the 5,656-member anchor set, i.e. the same definition as the EUR GWAS
# arm (12.1) and the control-disorder PRS runs (13.1, 13.6); it yields 4,116
# phenotyped subjects.  prs_assoc.R's own PC-distance fallback gives a
# different 5,361 and would make these betas incomparable with the panel.
EURIDS="$REPO/hpc/work/results/ancestry/eur_anchor.keep"

mkdir -p "$OUTD"
for f in "$PHENO" "$QUANT" "$CAT" "$MAN" "$EURIDS"; do
  [[ -s "$f" ]] || { echo "FATAL: $f missing/empty" >&2; exit 3; }
done
n_prof=$(ls "$PRSD"/score_*.profile 2>/dev/null | wc -l)
(( n_prof == 40 )) || { echo "FATAL: expected 40 C+T profiles in $PRSD, found $n_prof" >&2; exit 3; }

echo "=== population PRS association (prs_assoc.R)  $(date) ==="
"$RSCRIPT" "$REPO/tools/prs_assoc.R" \
  --prs-dir "$PRSD" --pheno "$PHENO" \
  --covar-quant "$QUANT" --covar-cat "$CAT" --manifest "$MAN" \
  --eur-ids "$EURIDS" --out "$OUTD/prs_association_v3.tsv"

echo "=== within-family PRS association (06_prs_family.R)  $(date) ==="
"$RSCRIPT" "$REPO/hpc_v2/R/06_prs_family.R" \
  --prs-dir "$PRSD" --pheno "$PHENO" \
  --covar-quant "$QUANT" --covar-cat "$CAT" --manifest "$MAN" \
  --eur-ids "$EURIDS" --out "$OUTD/prs_withinfamily_v3.tsv"

for f in prs_association_v3.tsv prs_withinfamily_v3.tsv; do
  [[ -s "$OUTD/$f" ]] || { echo "FATAL: $OUTD/$f missing/empty" >&2; exit 4; }
  echo "  $f: $(( $(wc -l < "$OUTD/$f") - 1 )) rows"
done
echo "=== done  $(date) ==="
