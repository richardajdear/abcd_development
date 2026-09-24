#!/bin/bash
# Step 2: impute C4 structural alleles into the MHC extract with Beagle.
#
#   JAVA=... BCFTOOLS=... BEAGLE_JAR=... CONFORM_JAR=... PANEL=... \
#     bash 02_impute_c4.sh <mhc.vcf.gz from step 1> <out prefix> [threads] [mem GB]
#
# 1. conform-gt aligns the target to the panel by POSITION: flips strand and
#    swaps REF/ALT where genotype correlation with the panel requires it, and
#    drops markers it cannot place (incl. ambiguous A/T, C/G it cannot resolve).
# 2. Beagle 5 (reference-based phasing + imputation) with the Sekar 2016 panel.
#    The genetic map is linear and compressed: cM = bp * MAP_CM_PER_BP, default
#    1e-6 (0.1 cM/Mb, ~10x below the genome average: long shared MHC
#    haplotypes).  The imputec4 protocol uses 1e-7; on the 5-fold leave-out
#    test (smoketest/, results/smoketest_map_scale_grid.tsv) 1e-6 imputed C4A
#    GREx best (r 0.89 vs 0.81 at 1e-7 with array-like density), so it is the
#    default.  Override with MAP_CM_PER_BP=... .
#    ap=true writes per-haplotype allele posteriors, which 03_c4_grex.py needs.
# Output: <out>.vcf.gz (+ .tbi), <out>.log, <out>.conform.log; marker overlap
# counts are echoed so the run log records how much of the panel the array hit.
set -euo pipefail
IN="$1"; OUT="$2"; THREADS="${3:-4}"; MEM="${4:-8}"
JAVA="${JAVA:-java}"; BCFTOOLS="${BCFTOOLS:-bcftools}"
: "${BEAGLE_JAR:?}" "${CONFORM_JAR:?}" "${PANEL:?}"
REGION="6:24894177-33890574"; MAP_CM_PER_BP="${MAP_CM_PER_BP:-1e-6}"
for f in "$IN" "$PANEL" "$BEAGLE_JAR" "$CONFORM_JAR"; do [[ -s "$f" ]] || { echo "missing $f" >&2; exit 2; }; done
mkdir -p "$(dirname "$OUT")"
# fresh scratch dir per run: conform-gt refuses to overwrite its output
W=$(mktemp -d "$OUT.tmp.XXXXXX")

"$JAVA" -jar "$CONFORM_JAR" ref="$PANEL" gt="$IN" chrom="$REGION" match=POS out="$W/conform" > "$W/conform.stdout"
cp "$W/conform.log" "$OUT.conform.log"
"$BCFTOOLS" index -f -t "$W/conform.vcf.gz"
n_in=$("$BCFTOOLS" index -n "$IN.tbi" 2>/dev/null || "$BCFTOOLS" view -H "$IN" | wc -l)
n_conf=$("$BCFTOOLS" index -n "$W/conform.vcf.gz"); n_panel=$("$BCFTOOLS" view -H "$PANEL" | wc -l | tr -d ' ')
echo "conform-gt: $n_in target markers -> $n_conf aligned to the panel ($n_panel panel markers)"
[[ "$n_conf" -gt 300 ]] || { echo "only $n_conf markers overlap the panel -- build/position mismatch?" >&2; exit 3; }

"$BCFTOOLS" query -f '%CHROM\t%POS\n' "$PANEL" | awk -v s="$MAP_CM_PER_BP" '{print $1"\t.\t"$2*s"\t"$2}' > "$W/flat.map"
"$JAVA" -Xmx"${MEM}"g -jar "$BEAGLE_JAR" gt="$W/conform.vcf.gz" ref="$PANEL" map="$W/flat.map" \
  chrom="$REGION" out="$W/beagle" nthreads="$THREADS" ap=true gp=true > "$W/beagle.stdout"
mv -f "$W/beagle.vcf.gz" "$OUT.vcf.gz"; mv -f "$W/beagle.log" "$OUT.log"
"$BCFTOOLS" index -f -t "$OUT.vcf.gz"
"$BCFTOOLS" view -H -i 'ID=="C4"' "$OUT.vcf.gz" | cut -f1-8 | cut -c1-300
echo "impute_c4: $(${BCFTOOLS} query -l "$OUT.vcf.gz" | wc -l | tr -d ' ') samples -> $OUT.vcf.gz"
