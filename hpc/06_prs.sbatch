#!/bin/bash
#SBATCH --job-name=abcd_prs
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=%x_%j.log
#
# Step 6: polygenic score association -- the PRIMARY test of whether the ABCD
# developmental phenotypes are associated with SCZ / MDD genetic risk.
#
# WHY THIS RATHER THAN rg
#
# LDSC (05) needs adequate h2 power on BOTH traits, and the slope phenotype does
# not have it (reliability ~0.21 at N~8,000).  A polygenic score moves the
# statistical burden onto the DISCOVERY GWAS, which is large (SCZ N~160k, MDD
# N~674k).  The score is then a single well-estimated regressor, and only the
# target-side regression depends on our N.  With ~4-5k EUR subjects, a score
# explaining even 0.5% of phenotype variance is detectable -- so unlike rg, a
# null here is informative.
#
# Direction of inference: score (exposure) -> developmental phenotype (outcome).
# This is not a causal claim.  It is the same association rg estimates, measured
# with the variance concentrated where the data can support it.
#
# ANCESTRY (the decision that most affects interpretation)
#
# The discovery GWAS are predominantly European and polygenic scores lose
# accuracy across ancestry, so the PRIMARY analysis is restricted to the
# genetically-EUR subset and the full multi-ancestry sample is reported as a
# sensitivity analysis.  PC adjustment controls stratification but does NOT
# repair transferability: in the full sample a null is ambiguous between "no
# effect" and "the score does not transfer", which is exactly the ambiguity the
# EUR-primary design avoids.  Both are computed; the summary labels which is
# which, and tools/prs_assoc.R fits them separately rather than pooling.
source "$(dirname "$0")/config.sh"
require_paths GENO PHENO COVAR_QUANT COVAR_CAT MANIFEST PLINK
ensure_dirs "$PRS_DIR"

names=()   # not mapfile: bash 3.2 compatibility, see 02_reml.sbatch
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)

# ---------------------------------------------------------------------------
# Step 6a: clump each disorder GWAS to approximately independent index SNPs.
#
# Clumping (rather than using all SNPs) because a raw sum double-counts LD
# blocks: without it, one strongly-associated region with dense genotyping
# contributes as much as many independent loci.  Clumping is the conventional
# P+T approach and needs only the reference panel we already have for MAGMA.
#
# Clumping happens in the REFERENCE panel, not in ABCD: LD should come from a
# large well-powered panel, and using the target sample's own LD would let the
# target's noise choose the index SNPs.
# ---------------------------------------------------------------------------
for D in SCZ MDD; do
  case $D in
    SCZ) f="$SCZ_SUMSTATS"; snp="$SCZ_SNP_COL"; p="$SCZ_P_COL"; eff="$SCZ_EFFECT_COL"; a1="$SCZ_A1_COL" ;;
    MDD) f="$MDD_SUMSTATS"; snp="$MDD_SNP_COL"; p="$MDD_P_COL"; eff="$MDD_EFFECT_COL"; a1="$MDD_A1_COL" ;;
  esac
  if [[ ! -f "$f" && "$DRY_RUN" != "1" ]]; then
    echo "SKIP $D: $f not found (set ${D}_SUMSTATS)" >&2
    continue
  fi

  # PLINK --clump wants columns literally named SNP and P.  Write a small
  # normalised copy carrying the effect allele and effect size too, resolving
  # every column BY NAME so a consortium's layout change cannot silently
  # misalign the score weights.
  norm="$PRS_DIR/$D.assoc"
  if [[ -f "$f" ]]; then
    awk -v s="$snp" -v p="$p" -v e="$eff" -v a="$a1" '
      NR==1 {
        for (i=1; i<=NF; i++) c[$i]=i
        for (k in c) ;
        if (!c[s] || !c[p] || !c[e] || !c[a]) {
          printf("ERROR: %s missing one of %s/%s/%s/%s\n", FILENAME, s,p,e,a) > "/dev/stderr"
          exit 1
        }
        print "SNP\tP\tA1\tEFFECT"; next
      }
      # Drop rows with a non-numeric or zero p-value; PLINK refuses the file
      # otherwise and the error it gives does not name the offending row.
      ($c[p]+0) > 0 && ($c[p]+0) <= 1 {
        print $c[s] "\t" $c[p] "\t" $c[a] "\t" $c[e]
      }' "$f" > "$norm" \
      || { echo "ERROR: could not normalise $f" >&2; exit 3; }
    echo "  $D: $(( $(wc -l < "$norm") - 1 )) SNPs with usable p-values"
  fi

  echo "=== 06_prs  clump $D  $(date) ==="
  run "$PLINK" \
    --bfile "$MAGMA_REF" \
    --clump "$norm" \
    --clump-p1 "$CLUMP_P1" \
    --clump-p2 "$CLUMP_P2" \
    --clump-r2 "$CLUMP_R2" \
    --clump-kb "$CLUMP_KB" \
    --out "$PRS_DIR/$D"

  # ---------------------------------------------------------------------------
  # Step 6b: score ABCD subjects at several p-value thresholds.
  #
  # Multiple thresholds because the optimal one is unknown a priori and depends
  # on the disorder's polygenicity.  This is a multiple-testing liability, so
  # the summary reports every threshold and tools/prs_assoc.R applies an
  # explicit correction rather than letting a best-threshold be quoted as if it
  # were prespecified.
  # ---------------------------------------------------------------------------
  if [[ -f "$PRS_DIR/$D.clumped" ]]; then
    awk 'NR>1 && $3!="" {print $3}' "$PRS_DIR/$D.clumped" > "$PRS_DIR/$D.index_snps"
    echo "  $D: $(wc -l < "$PRS_DIR/$D.index_snps" | tr -d ' ') index SNPs after clumping"
  fi
  for thr in $PRS_THRESHOLDS; do
    tag="${thr//./p}"
    if [[ -f "$PRS_DIR/$D.index_snps" ]]; then
      # Restrict to clumped index SNPs at this threshold; --score reads
      # SNP/A1/EFFECT from the normalised file.
      awk -v t="$thr" 'NR==FNR {keep[$1]=1; next}
                       FNR==1 {next}
                       ($1 in keep) && ($2+0) < (t+0) {print $1"\t"$3"\t"$4}' \
        "$PRS_DIR/$D.index_snps" "$norm" > "$PRS_DIR/${D}_$tag.weights"
      echo "  $D p<$thr: $(wc -l < "$PRS_DIR/${D}_$tag.weights" | tr -d ' ') SNPs in score"
    fi
    run "$PLINK" \
      --bfile "$GENO" \
      --score "$PRS_DIR/${D}_$tag.weights" 1 2 3 sum \
      --out "$PRS_DIR/score_${D}_$tag"
  done
done

# ---------------------------------------------------------------------------
# Step 6c: associate each score with each phenotype, in R.
#
# In R rather than awk because the model needs a family random effect: ABCD
# contains twins and siblings, and treating related subjects as independent
# would understate the standard errors.  See tools/prs_assoc.R.
# ---------------------------------------------------------------------------
if [[ "$DRY_RUN" != "1" ]]; then
  if command -v "$RSCRIPT" >/dev/null 2>&1; then
    echo "=== 06_prs  association tests  $(date) ==="
    "$RSCRIPT" "$REPO_ROOT/tools/prs_assoc.R" \
      --prs-dir "$PRS_DIR" \
      --pheno "$PHENO" \
      --covar-quant "$COVAR_QUANT" \
      --covar-cat "$COVAR_CAT" \
      --manifest "$MANIFEST" \
      --out "$PRS_DIR/prs_association.tsv"
    echo
    column -t -s $'\t' "$PRS_DIR/prs_association.tsv" 2>/dev/null \
      || cat "$PRS_DIR/prs_association.tsv"
    echo
    echo "Read the EUR rows as primary and the full-sample rows as sensitivity."
    echo "p_adj is Bonferroni across thresholds within each phenotype-disorder"
    echo "pair; a result significant at one threshold only, and not after"
    echo "correction, is a threshold-selection artefact rather than a finding."
  else
    echo "Rscript not found; scores are written but not tested." >&2
    echo "Run: Rscript tools/prs_assoc.R --prs-dir $PRS_DIR ..." >&2
  fi
fi

echo "=== 06_prs done  $(date) ==="
