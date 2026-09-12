#!/bin/bash
# Preflight: does each PLINK fileset's .bed actually match its .bim and .fam?
#
#   bash hpc/work/check_bfile_integrity.sh PREFIX [PREFIX ...]
#   bash hpc/work/check_bfile_integrity.sh "$GENO_ALLANC_DIR"/ABCD_chr{1..22}_hg19
#
# WHY THIS EXISTS
# ---------------
# The 4.0 all-ancestry cross-ancestry GRM build (job array, 22 tasks, 16 CPUs
# each) failed in seconds with GCTA's "Unexpected PLINK 1 .bed file size" -- on
# all 22 chromosomes.  The cause was that every .bim had been filtered to ~78 %
# of its variants without regenerating the .bed, so the two files disagreed
# about how many variants the fileset holds.
#
# That is worth catching BEFORE submitting an array job, and it is catchable
# from file sizes alone in about a second per chromosome.  A PLINK 1 .bed is
#
#     3 + ceil(n/4) * m  bytes        n = samples (.fam rows)
#                                     m = variants (.bim rows)
#
# so given n from the .fam, (size-3) must divide exactly by ceil(n/4), and the
# quotient must equal the .bim line count.  If the division is exact but the
# quotient disagrees with the .bim, the .bim was filtered without rebuilding the
# .bed -- exactly the 4.0 failure.  If the division is NOT exact, the .fam is
# the file that does not belong.
#
# The distinction matters because the two have different fixes: a stale .bim
# needs the original .bim recovered (the .bed column order is unrecoverable
# without it), whereas a mismatched .fam usually means the wrong .fam was
# copied alongside and the right one still exists.
#
# Exits non-zero if ANY fileset fails, so it can gate a submission:
#   bash hpc/work/check_bfile_integrity.sh ... && sbatch hpc/work/grm_allanc.sbatch
set -uo pipefail

[[ $# -ge 1 ]] || { echo "usage: $0 PREFIX [PREFIX ...]" >&2; exit 64; }

fail=0
pass=0

# stat(1) is not portable: GNU takes -c%s, BSD/macOS takes -f%z.
filesize() { stat -c%s "$1" 2>/dev/null || stat -f%z "$1"; }

for pfx in "$@"; do
  printf '%s\n' "$pfx"
  missing=""
  for ext in bed bim fam; do
    [[ -f "$pfx.$ext" ]] || missing="$missing .$ext"
  done
  if [[ -n "$missing" ]]; then
    echo "  FAIL  missing:$missing"
    fail=$((fail+1)); continue
  fi

  n=$(wc -l < "$pfx.fam" | tr -d ' ')
  m_bim=$(wc -l < "$pfx.bim" | tr -d ' ')
  bed=$(filesize "$pfx.bed")

  # Magic bytes: 0x6c 0x1b 0x01 is variant-major PLINK 1.  0x00 in byte 3 is
  # sample-major, which GCTA and PLINK 2 both refuse -- worth naming explicitly
  # rather than letting it surface as a size error downstream.
  magic=$(od -An -tx1 -N3 "$pfx.bed" | tr -d ' \n')
  if [[ "$magic" != "6c1b01" ]]; then
    if [[ "$magic" == "6c1b00" ]]; then
      echo "  FAIL  sample-major .bed (magic 6c1b00); convert with plink --make-bed"
    else
      echo "  FAIL  not a PLINK 1 .bed (magic $magic, expected 6c1b01)"
    fi
    fail=$((fail+1)); continue
  fi

  bytes_per_variant=$(( (n + 3) / 4 ))          # ceil(n/4)
  payload=$(( bed - 3 ))

  if (( payload % bytes_per_variant != 0 )); then
    implied=$(awk -v p="$payload" -v b="$bytes_per_variant" 'BEGIN{printf "%.4f", p/b}')
    echo "  FAIL  .bed size does not divide by ceil(n/4)"
    echo "        .fam n=$n -> $bytes_per_variant bytes/variant; (size-3)/that = $implied (not an integer)"
    echo "        => the .fam does not match this .bed (wrong .fam copied?)"
    fail=$((fail+1)); continue
  fi

  m_bed=$(( payload / bytes_per_variant ))
  if (( m_bed != m_bim )); then
    pct=$(awk -v a="$m_bim" -v b="$m_bed" 'BEGIN{printf "%.1f", 100*a/b}')
    echo "  FAIL  .bed holds $m_bed variants, .bim lists $m_bim (${pct}% of the .bed)"
    echo "        => the .bim was filtered without regenerating the .bed."
    echo "           The .bed column order cannot be recovered from the filtered"
    echo "           .bim alone; obtain the ORIGINAL .bim or re-derive the fileset."
    fail=$((fail+1)); continue
  fi

  echo "  ok    n=$n  m=$m_bim  bed=$bed bytes ($bytes_per_variant/variant)"
  pass=$((pass+1))
done

echo
echo "$pass passed, $fail failed"
(( fail == 0 )) || exit 1
