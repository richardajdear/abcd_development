#!/bin/bash
# Write a FID-matched id view of v1's dense pooled imputed GRM.
#
# WHY: hpc_v2 uses ONE phenotype file, and step 06 requires its FID column to
# be the real family id.  GCTA (step 05) matches an individual on the FID+IID
# PAIR, so the GRM it reads must spell FID the same way.  v1's .grm.id has
# FID = IID.  Rather than rewrite v1's GRM (read-only by policy: v1's results
# stay reproducible) this writes a new prefix whose .grm.bin/.grm.N.bin are
# SYMLINKS to v1's and whose .grm.id is the same file with column 1 replaced
# by family_id, in the SAME ROW ORDER -- the order is the matrix order and
# reordering it would silently permute the GRM.
#
# Subjects with no family_id (genotyped but not phenotyped) keep FID = IID;
# they are absent from $PHENO and GCTA drops them either way.
set -euo pipefail
source "$(cd "$(dirname "$0")/../.." && pwd)/config.sh"

SRC="${SRC:-$V1_ROOT/results/grm_imp_pooled/abcd_imp}"
FAMMAP="${FAMMAP:-$V1_ROOT/pheno_allanc/family_map.tsv}"
DST="$GRM_FULL"

[[ -f "$SRC.grm.bin" && -f "$SRC.grm.id" ]] || { echo "FATAL: $SRC incomplete" >&2; exit 2; }
[[ -f "$FAMMAP" ]] || { echo "FATAL: $FAMMAP not found" >&2; exit 2; }
ensure_dirs "$(dirname "$DST")"

ln -sf "$SRC.grm.bin"   "$DST.grm.bin"
ln -sf "$SRC.grm.N.bin" "$DST.grm.N.bin"

awk 'NR==FNR {if (FNR>1) fam[$1]=$2; next}
     {print (($2 in fam) ? fam[$2] : $1) "\t" $2}' "$FAMMAP" "$SRC.grm.id" > "$DST.grm.id"

n_src=$(wc -l < "$SRC.grm.id"); n_dst=$(wc -l < "$DST.grm.id")
[[ "$n_src" == "$n_dst" ]] || { echo "FATAL: row count changed $n_src -> $n_dst" >&2; exit 3; }
# IID column must be byte-identical to the source, in order.
diff <(awk '{print $2}' "$SRC.grm.id") <(awk '{print $2}' "$DST.grm.id") >/dev/null \
  || { echo "FATAL: IID column or its order changed" >&2; exit 3; }
n_re=$(awk '$1!=$2' "$DST.grm.id" | wc -l)
echo "ok  $DST.grm.id: $n_dst rows, $n_re FIDs replaced with family ids,"
echo "    $(awk '{print $1}' "$DST.grm.id" | sort -u | wc -l) distinct FIDs, IID order preserved."
