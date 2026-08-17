#!/bin/bash
#SBATCH --job-name=abcd_strat_reml
#SBATCH --account=VERTES-SL3-CPU   # override at submit: sbatch -A <acct>
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --output=%x_%j.log
#
# Ancestry-stratified REML on the all-ancestry GRM.  README_HPC.md section
# 5.3.2 makes this non-optional: a single pooled GRM assumes a common allele
# frequency and LD structure across ancestry groups, which is exactly what does
# not hold, and in-sample PCs mitigate that without repairing it.  If the pooled
# h2 differs materially from the stratified ones, the pooled number is not the
# headline.
#
# THREE THINGS THIS DOES THAT A NAIVE --keep WOULD NOT
# ----------------------------------------------------
# 1. **PCs are recomputed within each stratum.**  Global PCs inside one stratum
#    are near-constant -- that is the defect insample_pcs.py was written for,
#    measured at PC1 mean 12x its own SD inside EUR.  A covariate that does not
#    vary absorbs nothing.  So each stratum gets `--pca 10` on the full GRM
#    restricted to its own members.
#
# 2. **The EUR stratum is defined by the abcd_eur fileset, not by our
#    clustering.**  That is the one label here we did not infer ourselves, and
#    it is what makes this h2 directly comparable to the published EUR h2 of
#    0.575 / 0.228: same ancestry definition, and -- because it runs on this
#    GRM -- the same 515k array variants.  A difference against the published
#    number is then attributable to the variant set alone, which is the
#    confound that would otherwise be entangled with ancestry.  See the header
#    of split_allanc_chroms.sbatch for why the variant set differs at all.
#
# 3. **Relatedness pruning is inherited, not redone.**  Strata are cut out of
#    $GRM_ALLANC_UNREL, which is already pruned at 0.05 genome-wide; a subset of
#    a mutually-unrelated set is mutually unrelated, so re-pruning would only
#    lose subjects.
#
# Phenotypes come from $MANIFEST, so this cannot drift out of sync with the
# pooled run.  h2 is a variance ratio, so the pooled standardisation of
# phenotypes_gcta.txt is harmless here -- no per-stratum rescaling is needed.
set -euo pipefail

# Locate the repository.  Under sbatch, $0 and ${BASH_SOURCE[0]} name a COPY of
# this script in /var/spool/slurm/slurmd/jobNNN/, so the obvious
# `dirname "${BASH_SOURCE[0]}"/../..` resolves to /var/spool/slurm and every
# array task dies in two seconds on
#   /var/spool/slurm/hpc/config.sh: No such file or directory
# 02_reml.sbatch and 06_prs.sbatch already carry this scar; these work/ scripts
# did not, only because they had always been run with REPO_ROOT already
# exported into the submitting environment, which hid it.  Observed here on
# 2026-08-17, job 33789921, all 22 tasks.  Walk up from the submit directory
# (and $PWD) looking for hpc/config.sh, which is what actually identifies the
# repo, rather than trusting any one path to be the right one.
if [[ -z "${REPO_ROOT:-}" || ! -f "${REPO_ROOT:-}/hpc/config.sh" ]]; then
  REPO_ROOT=""
  for _d in "$(dirname "${BASH_SOURCE[0]}")/../.." "${SLURM_SUBMIT_DIR:-}" "$PWD"; do
    [[ -n "$_d" ]] || continue
    _p="$(cd "$_d" 2>/dev/null && pwd)" || continue
    while [[ -n "$_p" && "$_p" != "/" ]]; do
      if [[ -f "$_p/hpc/config.sh" ]]; then REPO_ROOT="$_p"; break 2; fi
      _p="$(dirname "$_p")"
    done
  done
fi
[[ -n "$REPO_ROOT" ]] || {
  echo "FATAL: cannot find hpc/config.sh above $PWD or ${SLURM_SUBMIT_DIR:-unset}." >&2
  echo "       Submit from inside the repo, or export REPO_ROOT." >&2
  exit 2
}
# shellcheck source=/dev/null
source "$REPO_ROOT/hpc/config.sh"

require_paths GRM_ALLANC GRM_ALLANC_UNREL PHENO COVAR_QUANT COVAR_CAT MANIFEST GCTA

STRATA_TSV="${STRATA_TSV:-$ABCD_HPC_ROOT/results/ancestry/strata.tsv}"
OUTDIR="${STRAT_REML_DIR:-$OUT/reml_strat}"
PY="${PY:-$ABCD_HPC_ROOT/envs/abcd/bin/python}"
T="${SLURM_CPUS_PER_TASK:-16}"

[[ -f "$STRATA_TSV" ]] || { echo "FATAL: no strata table at $STRATA_TSV; run work/assign_ancestry.py first" >&2; exit 2; }
ensure_dirs "$OUTDIR" "$LOG_DIR"

# Minimum stratum size.  Below this a REML either fails to converge or returns
# an estimate whose SE spans the whole parameter space, and reporting such a
# number next to the others invites it to be read as evidence of a difference.
MIN_N="${STRAT_MIN_N:-300}"

names=()
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)
echo "phenotypes: ${names[*]}"

# The strata to fit: every k-means cluster, plus EUR_anchor (column anchor_eur).
# `while read` rather than `mapfile` for the same reason 02_reml gives: mapfile
# is bash 4+ and macOS ships 3.2.
clusters=()
while IFS= read -r line; do clusters+=("$line"); done \
  < <(awk -F'\t' 'NR>1 {print $3}' "$STRATA_TSV" | sort -u)
strata=("EUR_anchor" "${clusters[@]}")
echo "strata: ${strata[*]}"

# The GRM's own .grm.id is the universe: a keep-list naming subjects GCTA does
# not have is silently ignored, so intersecting here makes the printed N real.
declare -A INGRM=()
while read -r fid iid; do INGRM["$iid"]="$fid"; done < "$GRM_ALLANC_UNREL.grm.id"
echo "unrelated GRM: ${#INGRM[@]} subjects"

for s in "${strata[@]}"; do
  keep="$OUTDIR/keep_$s.txt"
  if [[ "$s" == "EUR_anchor" ]]; then
    awk -F'\t' -v OFS='\t' 'NR>1 && $4==1 {print $1}' "$STRATA_TSV" > "$OUTDIR/.ids_$s"
  else
    awk -F'\t' -v s="$s" 'NR>1 && $3==s {print $1}' "$STRATA_TSV" > "$OUTDIR/.ids_$s"
  fi
  : > "$keep"
  # `if`, not `[[ ... ]] && printf`: under `set -e` a failed [[ ]] at the head
  # of an && list is an unhandled failure, so the FIRST subject not in the GRM
  # would kill the whole job -- and would do it silently enough to look like a
  # short keep-list rather than an abort.
  while read -r iid; do
    if [[ -n "${INGRM[$iid]:-}" ]]; then
      printf '%s\t%s\n' "${INGRM[$iid]}" "$iid" >> "$keep"
    fi
  done < "$OUTDIR/.ids_$s"
  rm -f "$OUTDIR/.ids_$s"

  n=$(wc -l < "$keep")
  echo
  echo "=== stratum $s: $n unrelated subjects  $(date) ==="
  if (( n < MIN_N )); then
    echo "  SKIP: below STRAT_MIN_N=$MIN_N"
    continue
  fi

  # (1) within-stratum PCs, from the FULL GRM restricted to this stratum.
  "$GCTA" --grm "$GRM_ALLANC" --keep "$keep" --pca 10 --thread-num "$T" \
          --out "$OUTDIR/pca_$s"

  # (2) a covariate directory of its own, so the swap is the tested code path
  #     (insample_pcs.py) rather than an ad-hoc awk rewrite, and so each
  #     stratum keeps its own covar_quant.pcsource provenance file.
  pd="$OUTDIR/covar_$s"
  mkdir -p "$pd"
  cp "$COVAR_QUANT" "$pd/covar_quant.txt"
  rm -f "$pd/covar_quant.txt.releasepcs"
  "$PY" "$REPO_ROOT/hpc/work/insample_pcs.py" "$pd" "$OUTDIR/pca_$s.eigenvec" --allow-missing

  # (3) REML per phenotype.  --keep on the pruned GRM: a subset of a mutually
  #     unrelated set is still mutually unrelated.
  for pheno in "${names[@]}"; do
    mpheno=$(mpheno_for "$pheno") || { echo "ERROR: '$pheno' not in $MANIFEST" >&2; exit 3; }
    echo "--- $s / $pheno (--mpheno $mpheno) ---"
    # A stratum can fail to converge; that is information, not a reason to kill
    # the other 24 fits, so this one command is allowed to fail.
    "$GCTA" --reml \
      --grm "$GRM_ALLANC_UNREL" \
      --keep "$keep" \
      --pheno "$PHENO" \
      --mpheno "$mpheno" \
      --qcovar "$pd/covar_quant.txt" \
      --covar "$COVAR_CAT" \
      --thread-num "$T" \
      --out "$OUTDIR/${s}__${pheno}" || echo "  (REML did not complete for $s / $pheno)"
  done
done

# ---------------------------------------------------------------------------
# Collect.  Same .hsq parsing as 02_reml, one row per stratum x phenotype, and
# rebuilt from whatever exists so a partial run still yields a readable table.
# ---------------------------------------------------------------------------
summary="$OUTDIR/reml_strat_summary.tsv"
{
  echo -e "stratum\tphenotype\th2\tse\tpval\tn"
  for s in "${strata[@]}"; do
    for pheno in "${names[@]}"; do
      f="$OUTDIR/${s}__${pheno}.hsq"
      [[ -f "$f" ]] || continue
      awk -v s="$s" -v p="$pheno" '
        $1=="V(G)/Vp" {h2=$2; se=$3}
        $1=="Pval"    {pv=$2}
        $1=="n"       {n=$2}
        END {if (h2!="") printf "%s\t%s\t%.4f\t%.4f\t%s\t%s\n", s, p, h2, se, pv, n}
      ' "$f"
    done
  done
} > "$summary"

echo
column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
echo
echo "wrote $summary"
echo "=== strat_reml done  $(date) ==="
