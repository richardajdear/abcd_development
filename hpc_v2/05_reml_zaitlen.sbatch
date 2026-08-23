#!/bin/bash
#SBATCH --job-name=v2_reml
#SBATCH --partition=icelake
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-5
#
# Step 5: Zaitlen two-GRM REML (Zaitlen et al. 2013, PLoS Genet), one
# phenotype per array task.
#
# WHY: v1's REML threw away every relative (--grm-cutoff 0.05 kept 5,649 of
# 8,082 phenotyped-genotyped subjects) because with relatives in a single-GRM
# fit, shared environment loads onto the additive term.  The Zaitlen method
# keeps ALL subjects and fits TWO variance components:
#
#   GRM1 = the full dense GRM (all pairs)
#   GRM2 = the same GRM with off-diagonals < BK_THRESH set to 0
#          (only close-relative pairs carry signal)
#
#   V(y) = s1*GRM1 + s2*GRM2 + se*I
#
#   h2_SNP  = s1 / Vp            (LD-tagged, comparable to v1's REML h2)
#   h2_ped  = (s1 + s2) / Vp     (pedigree h2, akin to twin h2)
#
# The GRM2 component absorbs the close-pair excess covariance (pedigree
# effects + shared environment), which is what protects h2_SNP from
# inflation while the relatives contribute to its precision.
#
# HONEST EXPECTATION: h2_SNP's SE improves less than raw N suggests, because
# relatives carry partly redundant genotype information.  The new quantity is
# h2_ped -- a family-design h2 from the SAME sample and model, bounding how
# much of the pedigree signal the SNPs tag.  Report both.
#
# INPUT: GRM_FULL must be the DENSE GRM over the FULL sample (relatives IN),
# built from the imputed genotypes -- v1's §8.15 GRM before any .unrel
# pruning.  The bK GRM is derived here.
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc_v2/config.sh (looked in $(dirname "$0"), $PWD, ${SLURM_SUBMIT_DIR:-unset})" >&2; exit 2; }
source "$_cfg"
require_paths GRM_FULL PHENO COVAR_QUANT COVAR_CAT MANIFEST GCTA
ensure_dirs "$REML2_DIR"

# --- 5a. Derive the thresholded (bK) GRM once; every array task checks. --------
# --make-bK zeroes off-diagonals below the threshold, keeping the diagonal --
# exactly Zaitlen's second component.  Cheap (one pass over the GRM), so task 1
# builds it and the rest wait via the file's existence; run_all.sh serialises
# by submitting task 1 alone first when the file is absent.
bk="$REML2_DIR/$(basename "$GRM_FULL").bK${BK_THRESH}"
mgrm="$REML2_DIR/mgrm.txt"
if [[ ! -f "$bk.grm.bin" ]]; then
  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" || "${SLURM_ARRAY_TASK_ID:-1}" == "1" ]]; then
    run "$GCTA" \
      --grm "$GRM_FULL" \
      --make-bK "$BK_THRESH" \
      --thread-num "$THREADS" \
      --out "$bk"
  else
    # Poll briefly rather than failing: task 1 is building it.
    for _ in $(seq 1 60); do [[ -f "$bk.grm.bin" ]] && break; sleep 10; done
    [[ -f "$bk.grm.bin" ]] || { echo "FATAL: $bk.grm.bin not built by task 1" >&2; exit 3; }
  fi
fi
if [[ "$DRY_RUN" != "1" ]]; then
  printf '%s\n%s\n' "$GRM_FULL" "$bk" > "$mgrm"
else
  echo "DRY_RUN: would write $mgrm with $GRM_FULL + $bk"
fi

# --- 5b. Two-component REML per phenotype ---------------------------------------
names=()   # not mapfile: bash 3.2 compatibility
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx=$((SLURM_ARRAY_TASK_ID - 1))
  if (( idx >= ${#names[@]} )); then
    echo "array index $SLURM_ARRAY_TASK_ID exceeds ${#names[@]} phenotypes; nothing to do"
    exit 0
  fi
  todo=("${names[$idx]}")
else
  todo=("${names[@]}")
fi

for pheno in "${todo[@]}"; do
  mpheno=$(mpheno_for "$pheno") || { echo "ERROR: '$pheno' not in $MANIFEST" >&2; exit 3; }
  echo "=== 05_reml_zaitlen  $pheno  (--mpheno $mpheno)  $(date) ==="
  # NO --grm-cutoff anywhere: the full sample, relatives included, is the point.
  # --reml-no-lrt off: we want the LRT-ready logL, default output includes it.
  run "$GCTA" \
    --reml \
    --mgrm "$mgrm" \
    --pheno "$PHENO" \
    --mpheno "$mpheno" \
    --qcovar "$COVAR_QUANT" \
    --covar "$COVAR_CAT" \
    --reml-est-fix \
    --thread-num "$THREADS" \
    --out "$REML2_DIR/$pheno"
done

# --- 5c. Collect ------------------------------------------------------------------
# Two-GRM .hsq layout differs from single: V(G1)/Vp, V(G2)/Vp and Sum of
# V(G)/Vp rows.  h2_SNP = V(G1)/Vp; h2_ped = Sum row.  Parse BY ROW NAME.
if [[ "$DRY_RUN" != "1" ]]; then
  summary="$REML2_DIR/reml_zaitlen_summary.tsv"
  {
    echo -e "phenotype\th2_snp\tse_snp\th2_ped\tse_ped\tpval\tn"
    for pheno in "${names[@]}"; do
      f="$REML2_DIR/$pheno.hsq"
      [[ -f "$f" ]] || continue
      awk -v p="$pheno" -F'\t' '
        $1=="V(G1)/Vp"        {h1=$2; s1=$3}
        $1=="Sum of V(G)/Vp"  {hs=$2; ss=$3}
        $1=="Pval"            {pv=$2}
        $1=="n"               {n=$2}
        END {if (h1!="") printf "%s\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\n", p, h1, s1, hs, ss, pv, n}
      ' "$f"
    done
  } > "$summary"
  echo
  column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
  echo
  echo "h2_snp = V(G1)/Vp (LD-tagged, compare with v1 reml_summary.tsv)."
  echo "h2_ped = V(G1)+V(G2) share (pedigree h2 from the same fit)."
  echo "n here should be the FULL phenotyped-genotyped sample -- if it matches"
  echo "v1's unrelated N instead, the wrong GRM (an .unrel one) was supplied."
fi

echo "=== 05_reml_zaitlen done  $(date) ==="
