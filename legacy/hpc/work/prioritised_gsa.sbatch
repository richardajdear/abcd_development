#!/bin/bash
#SBATCH --job-name=abcd_sczprio
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/home/rajd2/rds/hpc-work/abcd_development/hpc/work/logs/%x_%j.log
#
# Competitive MAGMA gene-set tests of every ABCD phenotype against the SCZ
# prioritised genes of Trubetskoy et al. 2022 (see make_prioritised_genesets.py
# for the source and the ID mapping).
#
# This reuses the .genes.raw files 04_magma already produced, so it re-runs no
# gene analysis -- only the set step, which takes seconds.
#
# The two disorder GWAS are included as controls, and they are controls of
# different kinds:
#
#   SCZ  is a POSITIVE CONTROL ONLY, and a circular one: these 120 genes were
#        prioritised *from this GWAS*, so a huge enrichment is guaranteed by
#        construction and means nothing biologically.  What it tests is the
#        plumbing -- that the Entrez mapping, the set file and --set-annot are
#        all being read as intended.  If this comes back null, everything below
#        it is uninterpretable.
#
#   MDD  is a genuine cross-disorder comparison: an independent GWAS, with the
#        same set.  It calibrates how much enrichment a well-powered psychiatric
#        GWAS shows in SCZ-prioritised genes.
#
# --model direction=greater: enrichment is directional.  The hypothesis is that
# ABCD association concentrates IN these genes, and a two-sided test would spend
# half its alpha on the uninteresting alternative that SCZ genes are unusually
# devoid of cortical signal.
set -euo pipefail

# REPO and MAGMA_DIR were hard-coded, which pinned this script to one account
# and to the 4.0/EUR MAGMA outputs.  Both are now overridable so the same harness
# runs against the imputed multi-ancestry gene results; the defaults reproduce
# the original behaviour exactly.
REPO=${REPO_ROOT:-/home/rajd2/rds/hpc-work/abcd_development}
W=$REPO/hpc/work
MAGMA=${MAGMA:-$W/bin/magma}
MAGMA_DIR=${MAGMA_DIR:-$W/results/magma}
# Parameterised so the same harness runs the SCZ and the MDD set files.
# TAG names the outputs and the conditioning set differs per disorder, so
# both are passed in rather than hard-coded.
SETS=${SETS:-$W/genesets/scz_prioritised.txt}
TAG=${TAG:-sczprio}
COND_SET=${COND_SET:-SCZ_locus_pool}
OUT=${OUT:-$W/results/magma_prio}
mkdir -p "$OUT"

[[ -f "$SETS" ]] || { echo "FATAL: $SETS missing; run the make_*_genesets.py step" >&2; exit 2; }

# NAMES was a hard-coded list, which pinned this harness to the settled five
# phenotypes.  Parameterised for hpc_v3 (Experiment A) in the same spirit as
# REPO/MAGMA_DIR/SETS/TAG/COND_SET above; the default is the original list
# verbatim, so every existing caller is unaffected.
NAMES=${NAMES:-"baseline_thickness global_slope slope_PC1 slope_PC2 slope_PC3 SCZ MDD"}

for name in $NAMES; do
  raw="$MAGMA_DIR/$name.genes.raw"
  if [[ ! -f "$raw" ]]; then
    echo "SKIP $name ($raw absent)" >&2
    continue
  fi
  echo "=== $name  $(date) ==="
  "$MAGMA" --gene-results "$raw" \
           --set-annot "$SETS" col=1,2 \
           --model direction=greater \
           --out "$OUT/${name}_$TAG"

  # Second pass, conditional on the full pool.  The marginal test above cannot
  # separate "enriched because these genes are causal for the disorder" from
  # "enriched because these genes are in the disorder's candidate pool at all":
  # the pool contains the prioritised set, so the two marginal tests overlap.
  # Conditioning asks whether the prioritised set explains anything the pool
  # does not.
  echo "--- $name conditional on $COND_SET ---"
  "$MAGMA" --gene-results "$raw" \
           --set-annot "$SETS" col=1,2 \
           --model direction=greater condition=$COND_SET \
           --out "$OUT/${name}_${TAG}_cond" || \
    echo "  (conditional pass failed for $name; marginal result above stands)" >&2
done

# ---------------------------------------------------------------------------
# Collect.  Column positions differ between the marginal and conditional
# outputs, so this parses by header name -- see the comment in the loop.
# ---------------------------------------------------------------------------
SUM="$OUT/${TAG}_gsa_summary.tsv"
printf 'trait\tmodel\tset\tsubmodel\tn_genes\tbeta\tbeta_std\tse\tp\n' > "$SUM"
for f in "$OUT"/*_$TAG.gsa.out "$OUT"/*_${TAG}_cond.gsa.out; do
  [[ -e "$f" ]] || continue
  b=$(basename "$f" .gsa.out)
  case "$b" in
    *_${TAG}_cond) trait=${b%_${TAG}_cond}; model=cond_on_pool ;;
    *)             trait=${b%_$TAG};          model=marginal ;;
  esac
  # Find columns by NAME, not position.  The conditional run inserts a MODEL
  # column that the marginal run does not have, so a fixed `$3 $4 $5...` layout
  # silently shifts every conditional row one column left -- which drops P
  # entirely and relabels SE as P.  Skip all `#` comment lines too: the
  # conditional output has five of them, one more than the marginal.
  awk -v t="$trait" -v m="$model" '
    /^#/ { next }
    $1 == "VARIABLE" { for (i = 1; i <= NF; i++) c[$i] = i; got = 1; next }
    got && NF {
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", t, m, $1,
             (c["MODEL"] ? $c["MODEL"] : "-"),
             $c["NGENES"], $c["BETA"], $c["BETA_STD"], $c["SE"], $c["P"]
    }' "$f" >> "$SUM"
done

echo
echo "=== summary ==="
column -t -s $'\t' "$SUM" 2>/dev/null || cat "$SUM"
echo
echo "wrote $SUM"
echo "=== done $(date) ==="
