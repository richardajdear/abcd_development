#!/bin/bash
# Preflight.  Verifies every input exists, that the binaries run, and that
# phenotype IDs actually intersect the .fam file -- before any expensive job is
# submitted.  Safe and quick; run it interactively on a login node.
#
#   bash hpc/00_check_inputs.sh
#
# Exit 0 means 01-03 have what they need.  MAGMA inputs are reported but do not
# fail the check, since 01-03 do not depend on them.
# Locate config.sh.  Under sbatch, $0 is a COPY of this script in
# /var/spool/slurm/slurmd/jobNNN/, so `dirname "$0"` does not contain
# config.sh -- run_all.sh's SLURM mode therefore sourced nothing, and because
# `set -euo pipefail` lives INSIDE config.sh it was never enabled either, so
# every step ran to completion with undefined functions and exited 0.  Search
# the submit directory too ($PWD, which run_all.sh sets via --chdir).
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc/config.sh (looked in $(dirname "$0"), $PWD, ${SLURM_SUBMIT_DIR:-unset})" >&2; exit 2; }
source "$_cfg"

fail=0
ok()   { echo "  ok       $*"; }
bad()  { echo "  MISSING  $*"; fail=1; }
warn() { echo "  absent   $*  (needed only by 04_magma)"; }

check_file() { [[ -f "$1" ]] && ok "$1" || bad "$1"; }
check_soft() { [[ -f "$1" ]] && ok "$1" || warn "$1"; }
check_exec() {
  if [[ ! -x "$1" ]]; then
    bad "$1 (not executable)"
    return
  fi
  # Executable is not sufficient: a Linux binary on a Mac, or a wrong-arch
  # build, is executable and fails only when run.  Actually starting it is the
  # only check that catches that.  GCTA exits non-zero on --help, so the test is
  # whether it produced output, not its exit status.
  if [[ -n "$("$1" --help 2>&1 | head -c 40)" ]]; then
    ok "$1"
  else
    bad "$1 (present but produced no output when run -- wrong architecture?)"
  fi
}

echo "genotypes:"
for e in bed bim fam; do check_file "$GENO.$e"; done

echo "binaries:"
check_exec "$GCTA"

echo "phenotypes:"
check_file "$PHENO"
check_file "$COVAR_QUANT"
check_file "$COVAR_CAT"
check_file "$MANIFEST"

echo "magma inputs:"
check_soft "$MAGMA_REF.bed"
check_soft "$MAGMA_GENE_LOC"

# ---------------------------------------------------------------------------
# The manifest must describe the file it ships with.  A stale manifest sends
# --mpheno at a column that has moved, and GCTA analyses the wrong phenotype
# without complaint.
# ---------------------------------------------------------------------------
if [[ -f "$MANIFEST" && -f "$PHENO" ]]; then
  echo "manifest consistency:"
  while read -r name mpheno _; do
    col=$(awk -v n="$name" 'NR==1{for(i=1;i<=NF;i++) if($i==n){print i; exit}}' "$PHENO")
    if [[ -z "$col" ]]; then
      bad "manifest lists '$name' but phenotypes_gcta.txt has no such column"
    elif [[ "$col" -ne $((mpheno + 2)) ]]; then
      bad "'$name': manifest says mpheno=$mpheno (column $((mpheno+2))) but it is column $col"
    else
      ok "$name -> --mpheno $mpheno"
    fi
  done < <(awk 'NR>1 {print $1, $2}' "$MANIFEST")
fi

# ---------------------------------------------------------------------------
# The failure mode that costs the most time: phenotype IDs that do not match
# the .fam file.  ABCD writes NDAR_INV... in genetics tables and
# sub-NDARINV... in imaging tables; abcd.gcta_export normalises to the genetics
# form, but verify rather than trust -- a mismatch yields an empty analysis, not
# an error.
# ---------------------------------------------------------------------------
if [[ -f "$PHENO" && -f "$GENO.fam" ]]; then
  echo "ID intersection:"
  n_pheno=$(awk 'NR>1{print $2}' "$PHENO" | sort -u | wc -l | tr -d ' ')
  n_fam=$(awk '{print $2}' "$GENO.fam" | sort -u | wc -l | tr -d ' ')
  n=$(comm -12 \
        <(awk 'NR>1{print $2}' "$PHENO" | sort -u) \
        <(awk '{print $2}' "$GENO.fam" | sort -u) | wc -l | tr -d ' ')
  echo "  phenotyped: $n_pheno   genotyped: $n_fam   intersection: $n"
  # Threshold as a fraction, not a fixed count: the same check then works on a
  # small test fixture, where a hard floor of 1000 would fail spuriously.
  min=$(( n_pheno / 2 ))
  if [[ "$n" -lt "$min" ]]; then
    bad "intersection $n is under half of $n_pheno -- check ID format (NDAR_INV vs sub-NDARINV)"
  else
    ok "intersection is $n of $n_pheno phenotyped subjects"
  fi
fi

echo
if [[ "$fail" == "0" ]]; then
  echo "PREFLIGHT PASSED -- 01_grm can be submitted."
else
  echo "PREFLIGHT FAILED -- fix the items marked MISSING above." >&2
fi
exit $fail
