#!/bin/bash
#SBATCH --job-name=magma_set_tests
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=slurm/%x_%A.log
#
# Step 15b: gene-set tests on the single-LMM phenotypes for THREE versions of
# the phenotype gene results, so the LD-reference effect and the sample effect
# can be separated:
#   EUR_1000G   magma_panel/genes          (step 13; 1000G EUR LD)
#   EUR_ABCD    magma_panel_abcdld/genes   (step 15a; ABCD EUR LD)
#   pooled_ABCD magma_pooled/genes         (step 14; ABCD analysis-sample LD)
# Sets: the step-13 SCZ and MDD sets plus magma_gene_sets/genesets_topgenes.txt
# (MDD / SCZ25 Bonferroni and top-100/250/500 genes from the EUR discovery
# GWAS's own gene analysis, MHC excluded; built by build_topgene_sets.py).
# -> results_70tab*/magma_set_tests/table_magma_set_tests.tsv
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"; PY="$REPO/legacy/hpc/work/envs/abcd/bin/python"
ROOT="$REPO/genetic_analysis/work/magma_scz2025"
SETS=( "scz:$ROOT/genesets_scz2025.txt"
       "mdd:$REPO/legacy/hpc/work/genesets/mdd_highconf.txt"
       "top:$REPO/genetic_analysis/magma_gene_sets/genesets_topgenes.txt" )
for S in "" "_hcp"; do
  R="$REPO/genetic_analysis/work/results_70tab$S"; OUT="$R/magma_set_tests"; mkdir -p "$OUT"
  for V in "EUR_1000G:magma_panel" "EUR_ABCD:magma_panel_abcdld" "pooled_ABCD:magma_pooled"; do
    TAG="${V%%:*}"; DIR="$R/${V#*:}/genes"
    for p in global_slope_1lmm baseline_thickness_1lmm; do
      raw="$DIR/$p.genes.raw"; [[ -s "$raw" ]] || { echo "MISSING $raw (skipped)"; continue; }
      for ss in "${SETS[@]}"; do
        "$MAGMA" --gene-results "$raw" --set-annot "${ss#*:}" col=1,2 \
          --model direction=greater --out "$OUT/${TAG}__${p}__${ss%%:*}" > /dev/null
      done
    done
  done
  "$PY" - "$OUT" <<'PYEOF'
import sys, glob, os, pandas as pd
out = sys.argv[1]; rows = []
for f in sorted(glob.glob(f"{out}/*.gsa.out")):
    version, pheno, _ = os.path.basename(f)[:-8].split("__")
    d = pd.read_csv(f, sep=r"\s+", comment="#")
    d = d.rename(columns=str.lower).assign(version=version, phenotype=pheno.replace("_1lmm", ""))
    rows.append(d[["version", "phenotype", "variable", "ngenes", "beta", "se", "p"]])
t = pd.concat(rows).drop_duplicates(["version", "phenotype", "variable"])
t.to_csv(f"{out}/table_magma_set_tests.tsv", sep="\t", index=False)
print(t.pivot_table(index="variable", columns=["version", "phenotype"], values="p").round(4).to_string())
PYEOF
done
