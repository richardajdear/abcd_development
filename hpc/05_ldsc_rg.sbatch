#!/bin/bash
#SBATCH --job-name=abcd_ldsc
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x_%j.log
#
# Step 5: genetic correlation (LDSC) between each ABCD developmental phenotype
# and the SCZ / MDD case-control GWAS.
#
# WHAT THIS TESTS, AND WHY IT IS THE SECONDARY ANALYSIS
#
# rg asks whether the two traits' genome-wide additive signals covary, using
# summary statistics only.  It is the conventional answer to "are these two GWAS
# associated", and it is symmetric and needs no individual genotypes.
#
# But LDSC is not well powered for OUR side of the comparison.  The published
# guidance is that rg is unreliable when the trait's own h2 z-score is below
# ~4.  With N~8,000 and slope reliability ~0.21, the slope h2 z-score is
# expected to be ~2 (see 02_reml output), so rg standard errors here will
# plausibly be +/-0.3-0.5.  That is wide enough that a null result would be
# uninformative -- it would not distinguish "no genetic overlap" from "not
# enough signal to tell".
#
# Run it anyway, for two reasons: the h2 and intercept it reports are a useful
# check on the GWAS itself, and the baseline_thickness phenotype (reliability
# ~0.9) IS adequately powered, so its rg is interpretable and acts as the
# positive control.  But read 06_prs as the primary test of this hypothesis.
#
# INTERPRETING THE OUTPUT
#   rg          genetic correlation, -1..1 (can exceed slightly from noise)
#   se, p       the quantity that decides whether anything was learned
#   h2_obs/h2   per-trait h2 on the observed scale; z = h2/se < 4 means the
#               rg for that trait should be treated as underpowered
#   gcov_int    genetic covariance intercept; should be ~0 for independent
#               samples.  Non-zero indicates sample overlap between ABCD and
#               the disorder GWAS, which would bias rg.
source "$(dirname "$0")/config.sh"
require_paths GWAS_DIR MANIFEST LDSC LDSC_MUNGE LD_REF LD_WEIGHTS HM3_SNPLIST
ensure_dirs "$LDSC_DIR"

names=()   # not mapfile: bash 3.2 compatibility, see 02_reml.sbatch
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)

# ---------------------------------------------------------------------------
# Munge: harmonise each GWAS to the HapMap3 SNP set with consistent alleles.
#
# --merge-alleles is not optional.  It restricts to HM3 and, critically, fixes
# the allele orientation against a reference; without it a strand flip in one
# input silently reverses the sign of rg.
# ---------------------------------------------------------------------------
munge_one() {   # munge_one <in.tsv> <snp-col> <p-col> <signed-spec> <n-spec> <out-prefix>
  local f="$1" snp="$2" p="$3" signed="$4" nspec="$5" out="$6"
  [[ -f "$out.sumstats.gz" ]] && { echo "  reuse $out.sumstats.gz"; return 0; }
  # shellcheck disable=SC2086
  run "$LDSC_MUNGE" \
    --sumstats "$f" \
    --snp "$snp" --p "$p" \
    --signed-sumstats "$signed" \
    $nspec \
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000 \
    --out "$out"
}

# ABCD phenotypes: fastGWA gives BETA/SE, so the signed statistic is BETA with
# a null of 0, and N is a per-SNP column.
for pheno in "${names[@]}"; do
  f="$GWAS_DIR/$pheno.fastGWA"
  if [[ ! -f "$f" && "$DRY_RUN" != "1" ]]; then
    echo "SKIP $pheno: $f not found (run 03_gwas first)" >&2
    continue
  fi
  echo "=== 05_ldsc  munge $pheno  $(date) ==="
  munge_one "$f" SNP P BETA,0 "--N-col N" "$LDSC_DIR/$pheno"
done

# Disorder GWAS: column names and the signed statistic differ per consortium,
# so they come from config rather than being assumed here.
for D in SCZ MDD; do
  case $D in
    SCZ) f="$SCZ_SUMSTATS"; snp="$SCZ_SNP_COL"; p="$SCZ_P_COL"
         signed="$SCZ_SIGNED"; nspec="$SCZ_N_SPEC" ;;
    MDD) f="$MDD_SUMSTATS"; snp="$MDD_SNP_COL"; p="$MDD_P_COL"
         signed="$MDD_SIGNED"; nspec="$MDD_N_SPEC" ;;
  esac
  if [[ ! -f "$f" && "$DRY_RUN" != "1" ]]; then
    echo "SKIP $D: $f not found (set ${D}_SUMSTATS)" >&2
    continue
  fi
  echo "=== 05_ldsc  munge $D  $(date) ==="
  munge_one "$f" "$snp" "$p" "$signed" "$nspec" "$LDSC_DIR/$D"
done

# ---------------------------------------------------------------------------
# rg: one LDSC call per disorder, with every ABCD phenotype in the same call.
#
# Batching this way is deliberate -- LDSC's --rg takes the first file as the
# reference and computes rg against each subsequent one, so one call per
# disorder gives all phenotypes against it while reading the LD scores once.
# ---------------------------------------------------------------------------
for D in SCZ MDD; do
  [[ -f "$LDSC_DIR/$D.sumstats.gz" || "$DRY_RUN" == "1" ]] || continue
  list="$LDSC_DIR/$D.sumstats.gz"
  for pheno in "${names[@]}"; do
    [[ -f "$LDSC_DIR/$pheno.sumstats.gz" || "$DRY_RUN" == "1" ]] || continue
    list="$list,$LDSC_DIR/$pheno.sumstats.gz"
  done
  case "$list" in
    *,*) : ;;
    *) echo "SKIP $D rg: no ABCD phenotypes munged" >&2; continue ;;
  esac
  echo "=== 05_ldsc  rg $D vs ABCD phenotypes  $(date) ==="
  run "$LDSC" \
    --rg "$list" \
    --ref-ld-chr "$LD_REF" \
    --w-ld-chr "$LD_WEIGHTS" \
    --out "$LDSC_DIR/rg_$D"
done

# ---------------------------------------------------------------------------
# Collect: LDSC writes the rg table as a whitespace block in the .log, after a
# "Summary of Genetic Correlation Results" line.  Parse from that marker rather
# than by line number, and resolve columns by name.
# ---------------------------------------------------------------------------
if [[ "$DRY_RUN" != "1" ]]; then
  summary="$LDSC_DIR/ldsc_rg_summary.tsv"
  {
    echo -e "disorder\tphenotype\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_z\tgcov_int\tunderpowered"
    for D in SCZ MDD; do
      log="$LDSC_DIR/rg_$D.log"
      [[ -f "$log" ]] || continue
      awk -v d="$D" '
        /Summary of Genetic Correlation Results/ {inblk=1; next}
        inblk && /^[[:space:]]*$/ {if (seen) inblk=0; next}
        inblk && /p2/ {for (i=1; i<=NF; i++) c[$i]=i; seen=1; next}
        inblk && seen && NF>3 {
          # p2 is the ABCD phenotype path; reduce to the phenotype name.
          ph=$c["p2"]; sub(/.*\//,"",ph); sub(/\.sumstats\.gz$/,"",ph)
          h2=$c["h2_obs"]; h2se=$c["h2_obs_se"]
          z=(h2se+0!=0 ? h2/h2se : "NA")
          # LDSC guidance: rg is unreliable below h2 z-score ~4.
          under=((z!="NA" && z<4) ? "yes" : "no")
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", d, ph,
            $c["rg"], $c["se"], $c["z"], $c["p"], h2, h2se,
            (z=="NA"?"NA":sprintf("%.2f",z)), $c["gcov_int"], under
        }' "$log"
    done
  } > "$summary"
  echo
  column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
  echo
  echo "underpowered=yes means that phenotype's own h2 z-score is below 4, so"
  echo "its rg is not interpretable regardless of the p-value.  Expect this for"
  echo "the slope phenotypes and NOT for baseline_thickness."
  echo "gcov_int far from 0 indicates sample overlap and would bias rg."
fi

echo "=== 05_ldsc done  $(date) ==="
