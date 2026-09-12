#!/bin/bash
#SBATCH --job-name=v2_lava
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#
# LOCAL GENETIC CORRELATION (LAVA; Werme et al. 2022, Nat Genet).
#
# WHY THIS, AND WHY NOW.  Genome-wide rg averages the genetic covariance over
# the whole genome and came back null (-0.12 to -0.20, never significant),
# which at our h2 is what an underpowered test looks like regardless of the
# truth.  But MAGMA says something more specific: the SCZ signal in our data
# concentrates in the SCZ LOCUS POOL (462 genes; p = 1.1e-04 for
# baseline_thickness, 7.4e-03 for global_slope).  If the shared architecture is
# LOCALISED rather than spread evenly, a genome-wide average is exactly the
# wrong statistic -- it dilutes a real local signal with ~2,400 null blocks.
#
# LAVA estimates rg WITHIN each of 2,495 approximately-independent LD blocks.
# That converts the MAGMA gene-set enrichment from "these genes are enriched"
# into a testable statistical statement: "the genetic correlation between SCZ
# and cortical thinning is non-zero in these specific regions".  Two tests:
#
#   (a) TARGETED, the pre-specified one: only blocks containing SCZ locus-pool
#       genes.  ~100-200 tests instead of 2,495, so the multiple-testing burden
#       falls by an order of magnitude, and the hypothesis was set by MAGMA
#       BEFORE seeing any local rg -- it is not selected on the outcome.
#   (b) Genome-wide, for context and to see whether targeted blocks stand out.
#
# LAVA requires a univariate local h2 to be significant in BOTH traits before
# it will estimate a bivariate local rg; that filter is built into the tool and
# is why most blocks will return nothing at our sample size.  A null result
# here is therefore about power, not about absence of overlap -- the number to
# report is how many blocks passed the univariate filter at all.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
OUTD="$OUT_V2/lava"
mkdir -p "$OUTD"
RS="/rds/user/rajd2/hpc-work/envs/genesis/bin/Rscript"

# LAVA reads PLINK genotypes as its LD reference.  g1000_eur is GRCh37, matching
# the block definitions, and is the same panel MAGMA used -- so the two analyses
# describe LD identically.
REF="${MAGMA_REF:-/rds/user/rajd2/hpc-work/magma/reference_data/g1000_eur}"
LOCI="$REPO/hpc_v2/work/refs/lava_blocks_GRCh37.locfile"

"$RS" "$REPO/hpc_v2/work/setup/lava_local_rg.R" \
  --ref "$REF" --loci "$LOCI" \
  --gene-loc "/rds/user/rajd2/hpc-work/magma/gene_locations/NCBI37.3.gene.loc" \
  --geneset "$W/genesets/scz_prioritised.txt" \
  --set-name SCZ_locus_pool \
  --abcd-dir "$OUT_V2/ldsc_eur" \
  --scz "$OUT_V2/ldsc_eur_matched/SCZ_core.sumstats.gz" \
  --mdd "$OUT_V2/mdd_eur/MDD_EUR.sumstats.gz" \
  --out-dir "$OUTD"
echo "=== done  $(date) ==="
