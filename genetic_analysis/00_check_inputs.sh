#!/bin/bash
# Preflight for hpc_v2: seconds to run, catches the failures that otherwise
# surface hours into a queue.  Run before anything else; run_all.sh runs it
# automatically and refuses to submit if it fails.
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate genetic_analysis/config.sh" >&2; exit 2; }
source "$_cfg"

fail=0
note() { echo "  $*"; }
bad()  { echo "  FAIL  $*" >&2; fail=1; }

echo "=== hpc_v2 preflight  $(date) ==="

# --- required files -----------------------------------------------------------
require_paths GENO_ARRAY PHENO COVAR_QUANT COVAR_CAT MANIFEST RSCRIPT

# --- PLINK bed/bim/fam integrity (the 4.0 defect detector, inherited) ----------
# A PLINK1 .bed is exactly 3 + ceil(n/4)*m bytes.  A mismatch means the .bim or
# .fam was filtered without regenerating the .bed; DO NOT work around it.
check_bfile() {
  local prefix="$1"
  [[ -f "$prefix.bed" && -f "$prefix.bim" && -f "$prefix.fam" ]] || { bad "$prefix: missing member file"; return; }
  local n m size expect magic
  n=$(wc -l < "$prefix.fam" | tr -d ' ')
  m=$(wc -l < "$prefix.bim" | tr -d ' ')
  size=$(wc -c < "$prefix.bed" | tr -d ' ')
  expect=$((3 + ((n + 3) / 4) * m))
  magic=$(od -An -tx1 -N3 "$prefix.bed" | tr -d ' \n')
  if [[ "$magic" != "6c1b01" ]]; then
    bad "$prefix.bed magic=$magic (want 6c1b01; 6c1b00 = sample-major, fix with plink --make-bed)"
  elif [[ "$size" -ne "$expect" ]]; then
    bad "$prefix.bed is $size bytes, expected $expect (n=$n, m=$m) -- bim/fam/bed mismatch, see legacy/hpc/README_HPC.md §5.1"
  else
    note "ok  $prefix  (n=$n, m=$m)"
  fi
}
check_bfile "$GENO_ARRAY"

# Imputed per-chromosome filesets: check whichever exist, count them.
n_imp=0
for chr in $(seq 1 22); do
  p=$(geno_imp_prefix "$chr")
  [[ -f "$p.bed" ]] && { check_bfile "$p"; n_imp=$((n_imp + 1)); }
done
if [[ "$n_imp" -eq 0 ]]; then
  note "WARN: no imputed filesets under $GENO_IMP_DIR (template $GENO_IMP_TPL)."
  note "      Steps 04 (association) cannot run; 01-03 can."
elif [[ "$n_imp" -lt 22 ]]; then
  bad "only $n_imp/22 imputed filesets present -- an association scan on a partial genome is silently wrong"
fi

# --- ID sanity ------------------------------------------------------------------
# The genotype .fam on CSD3 has FID=IID; the export carries real family IDs.
# The 8-char NDAR token must intersect between the two, or every downstream
# join is empty (GCTA-style tools then report a zero-subject analysis with NO
# error -- v1's most expensive class of bug).
n_tok=$(awk 'NR>1{id=$2; sub(/^sub-/,"",id); gsub(/_/,"",id); print toupper(id)}' "$PHENO" \
        | sort -u \
        | join - <(awk '{id=$2; sub(/^sub-/,"",id); gsub(/_/,"",id); print toupper(id)}' "$GENO_ARRAY.fam" | sort -u) \
        | wc -l | tr -d ' ')
n_ph=$(awk 'NR>1' "$PHENO" | wc -l | tr -d ' ')
if [[ "$n_tok" -eq 0 ]]; then
  bad "phenotype IDs and $GENO_ARRAY.fam IDs do not intersect at all (n_pheno=$n_ph)"
else
  note "ok  ID token join: $n_tok of $n_ph phenotyped subjects present in the array .fam"
fi

# Family structure must exist in the EXPORT (step 06 needs it).
n_fam_multi=$(awk 'NR>1{print $1}' "$PHENO" | sort | uniq -c | awk '$1>1' | wc -l | tr -d ' ')
if [[ "$n_fam_multi" -eq 0 ]]; then
  bad "no multi-member families in $PHENO -- FID column has lost the family structure; re-export (src/abcd/gcta_export.py)"
else
  note "ok  $n_fam_multi multi-member families in the phenotype export"
fi

# --- manifest sanity ------------------------------------------------------------
n_pheno_cols=$(head -1 "$PHENO" | wc -w | tr -d ' ')
while read -r name mp; do
  col=$((mp + 2))
  if (( col > n_pheno_cols )); then
    bad "manifest says $name is --mpheno $mp but $PHENO has only $((n_pheno_cols - 2)) phenotype columns"
  fi
done < <(awk 'NR>1{print $1, $2}' "$MANIFEST")
note "ok  manifest covers $(phenotype_names | wc -l | tr -d ' ') phenotypes within $((n_pheno_cols - 2)) exported columns"

# --- R environment ----------------------------------------------------------------
if "$RSCRIPT" -e 'suppressMessages({library(GENESIS); library(SNPRelate); library(GWASTools)})' >/dev/null 2>&1; then
  note "ok  RSCRIPT has GENESIS/SNPRelate/GWASTools ($("$RSCRIPT" --version 2>&1 | head -1))"
else
  bad "RSCRIPT=$RSCRIPT cannot load GENESIS/SNPRelate/GWASTools -- build the env from genetic_analysis/envs/genesis_env.yml"
fi
if "$RSCRIPT" -e 'suppressMessages({library(lme4); library(lmerTest); library(data.table); library(optparse)})' >/dev/null 2>&1; then
  note "ok  RSCRIPT has lme4/lmerTest/data.table/optparse (step 06)"
else
  bad "RSCRIPT=$RSCRIPT lacks lme4/lmerTest/data.table/optparse (step 06 needs them)"
fi

# --- steps 05/06 inputs (warn, not fail: they can lag the GENESIS track) ---------
if [[ ! -f "$GRM_FULL.grm.bin" ]]; then
  note "WARN: GRM_FULL=$GRM_FULL not found -- step 05 (Zaitlen REML) cannot run yet."
fi
if ! ls "$PRS_PROFILE_DIR"/score_*.profile >/dev/null 2>&1; then
  note "WARN: no score_*.profile under $PRS_PROFILE_DIR -- step 06 cannot run yet (v1 06_prs writes them)."
fi

echo
if [[ "$fail" -eq 0 ]]; then
  echo "preflight PASSED"
else
  echo "preflight FAILED -- fix the FAIL lines above before submitting" >&2
  exit 1
fi
