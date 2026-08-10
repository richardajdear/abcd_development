#!/bin/bash
#SBATCH --job-name=abcd_grm
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=%x_%j.log
#
# Step 1: genetic relationship matrices.
#
# Produces three things, because the two downstream steps need different ones:
#
#   $GRM            full GRM over all genotyped subjects      -> input to both
#   $GRM_UNREL      relatedness-pruned GRM (its .grm.id is the subject list) -> 02_reml
#   $GRM_SPARSE     sparse GRM (off-diagonals below cutoff zeroed) -> 03_gwas
#
# Runs under sbatch on CSD3 or directly with bash (see hpc/run_all.sh).
source "$(dirname "$0")/config.sh"
require_paths GENO GCTA
ensure_dirs "$GRM_DIR"

echo "=== 01_grm  $(date) ==="
echo "genotypes: $GENO"
echo "output:    $GRM_DIR"

# ---------------------------------------------------------------------------
# Full GRM.
#
# Only --maf is applied here, because GCTA implements no other SNP filter:
# --geno and --hwe are PLINK options and GCTA exits with "invalid option" on
# both.  Per-SNP missingness and HWE filtering therefore have to happen upstream.
# The ABCD release used here is distributed already QC'd (hence the default
# GENO name), so this is normally nothing to do.  If you point GENO at an
# unfiltered fileset, filter it first:
#
#   plink --bfile RAW --autosome --geno 0.02 --hwe 1e-6 --maf 0.01 \
#         --make-bed --out RAW_qc
#
# and set GENO to the result.  00_check_inputs does not verify QC status --
# there is no way to tell from a .bed whether it has been filtered.
#
# For ~8k subjects this is a single-part GRM.  If it runs out of time or memory,
# split it with --part i n and concatenate; at this N that is not needed.
# ---------------------------------------------------------------------------
run "$GCTA" \
  --bfile "$GENO" \
  --autosome \
  --maf "$MAF" \
  --make-grm \
  --thread-num "$THREADS" \
  --out "$GRM"

# ---------------------------------------------------------------------------
# Unrelated subset for REML.
#
# This is the step that decides what REML can say.  ABCD is a family study
# built on twin and sibling pairs; removing everyone related above 0.05 discards
# a large share of the sample.  That is not optional -- with relatives left in,
# shared environment loads onto the additive genetic term and SNP-h2 is not
# SNP-h2 -- but it does mean the REML N is far below the phenotyped 8,192 and
# the SE will be wide.  Report the retained N with every h2.
#
# --grm-cutoff writes the pruned GRM; the .grm.id of that output IS the
# unrelated subject list, which is what 02_reml uses via --keep.
# ---------------------------------------------------------------------------
run "$GCTA" \
  --grm "$GRM" \
  --grm-cutoff "$GRM_CUTOFF" \
  --make-grm \
  --thread-num "$THREADS" \
  --out "$GRM_UNREL"

# ---------------------------------------------------------------------------
# Sparse GRM for fastGWA.
#
# fastGWA's mixed model needs a sparse GRM, not the dense one: it keeps
# close relatives' covariance and zeroes the rest, which is what makes it
# tractable.  Unlike REML, fastGWA *keeps* the relatives -- the sparse GRM is
# how it accounts for them.  So 03_gwas runs on the full sample while 02_reml
# runs on the unrelated subset, and their N differ by design, not by error.
# ---------------------------------------------------------------------------
run "$GCTA" \
  --grm "$GRM" \
  --make-bK-sparse "$SPARSE_CUTOFF" \
  --thread-num "$THREADS" \
  --out "$GRM_SPARSE"

# ---------------------------------------------------------------------------
# Report the three N values, since every downstream interpretation depends on
# them and they are otherwise buried in three separate GCTA logs.
# ---------------------------------------------------------------------------
if [[ "$DRY_RUN" != "1" ]]; then
  n_full=$(wc -l < "$GRM.grm.id" | tr -d ' ')
  n_unrel=$(wc -l < "$GRM_UNREL.grm.id" | tr -d ' ')
  {
    echo -e "quantity\tn"
    echo -e "genotyped_with_grm\t$n_full"
    echo -e "unrelated_at_${GRM_CUTOFF}\t$n_unrel"
    echo -e "removed_as_related\t$((n_full - n_unrel))"
  } | tee "$GRM_DIR/grm_sample_sizes.tsv"
fi

echo "=== 01_grm done  $(date) ==="
