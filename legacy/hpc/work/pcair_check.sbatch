#!/bin/bash
#SBATCH --job-name=abcd_pcair
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=/home/rajd2/rds/hpc-work/abcd_development/hpc/work/logs/%x_%j.log
#
# Does the h2 = 0.575 for baseline_thickness depend on the ancestry PCs being
# the wrong ones?
#
# The covariates in use are the release's `ab_g_stc__gen_pc__01..10`, computed
# on the full multi-ancestry ABCD cohort.  Restricted to the 4,126 EUR analysis
# subjects, PC1 has mean 0.00631 and SD 0.000523 -- the mean is 12x the spread.
# A PC that separates ancestry groups looks exactly like that inside one group,
# and it cannot absorb within-EUR structure (e.g. north-west vs south European
# gradients), which is what actually needs controlling here.  Uncorrected
# structure inflates a GRM-based h2.
#
# So: recompute PCs *within* the EUR sample from the existing GRM, and re-run
# REML on baseline_thickness with them.  The GRM already exists, so the PCA is
# minutes and the REML seconds.  If h2 barely moves, the release PCs were
# harmless here and 0.575 stands on its own; if it drops toward the published
# 0.20-0.40, the covariates were the problem.
#
# Either way the comparison is the answer, not an argument about it.
set -euo pipefail

REPO=/home/rajd2/rds/hpc-work/abcd_development
W=$REPO/hpc/work
GCTA="$W/bin/gcta64"
T="${SLURM_CPUS_PER_TASK:-16}"
OUT="$W/results/pcair_check"
mkdir -p "$OUT"

GRM=/rds/user/rajd2/hpc-work/ABCD/genetics/abcd_eur
GRM_UNREL="$W/results/grm/abcd_full.unrel"
PHENO="$W/pheno/phenotypes_gcta.txt"
QCOV="$W/pheno/covar_quant.txt"
CCOV="$W/pheno/covar_categorical.txt"

for f in "$GRM.grm.bin" "$GRM_UNREL.grm.bin" "$PHENO" "$QCOV" "$CCOV"; do
  [[ -f "$f" ]] || { echo "FATAL: $f not found" >&2; exit 2; }
done



# Build a covariate file with the in-sample PCs in place of the release ones,
# keeping baseline_age and n_visits exactly as they are so the only thing that
# differs between the two REML runs is the PC block.
echo
echo "=== rebuilding covariates with in-sample PCs  $(date) ==="
"$W/envs/abcd/bin/python" - "$QCOV" /home/rajd2/rds/rds-genetics_hpc-Nl99R8pHODQ/ABCD/abcd_eur.eigenvec "$OUT/covar_quant_inpc.txt" <<'PY'
import sys
import pandas as pd

qcov, eig, out = sys.argv[1:4]
q = pd.read_csv(qcov, sep=" ")
# GCTA writes eigenvec headerless: FID IID PC1..PC20
e = pd.read_csv(eig, sep=r"\s+", header=None)
e.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, e.shape[1] - 1)]

keep = [c for c in q.columns if not c.startswith("PC")]
pcs = [f"PC{i}" for i in range(1, 11)]
m = q[keep].merge(e[["FID", "IID"] + pcs], on=["FID", "IID"], how="inner")
m.to_csv(out, sep=" ", index=False, na_rep="NA")

print(f"  {len(q)} subjects in covar_quant, {len(e)} in eigenvec, {len(m)} merged")
print(f"  in-sample PC SDs: " + ", ".join(f"{c}={m[c].std():.4f}" for c in pcs[:4]))
print(f"  release   PC SDs: " + ", ".join(
    f"{c}={q[c].std():.6f}" for c in pcs[:4] if c in q.columns))
PY

run_reml () {   # run_reml <label> <qcovar-file>
  echo
  echo "=== REML baseline_thickness with $1 PCs  $(date) ==="
  "$GCTA" --reml --grm "$GRM_UNREL" \
          --pheno "$PHENO" --mpheno 1 \
          --qcovar "$2" --covar "$CCOV" \
          --thread-num "$T" --out "$OUT/reml_$1" || {
    echo "  REML failed for $1" >&2; return 0; }
  awk '/^V\(G\)\/Vp/{printf "  h2 = %s  SE %s\n", $2, $3}
       /^Pval/       {printf "  p  = %s\n", $2}
       /^n /         {printf "  n  = %s\n", $2}' "$OUT/reml_$1.hsq"
}

run_reml release  "$QCOV"
run_reml insample "$OUT/covar_quant_inpc.txt"

echo
echo "================= COMPARISON ================="
for lab in release insample; do
  f="$OUT/reml_$lab.hsq"
  [[ -f "$f" ]] && awk -v L="$lab" '/^V\(G\)\/Vp/{printf "  %-9s h2 = %-8s SE %-8s\n", L, $2, $3}' "$f"
done
echo "=============================================="
echo "=== done  $(date) ==="
