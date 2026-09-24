#!/bin/bash
# Reproduces results/smoketest_map_scale_grid.tsv: leave-fold-out C4 accuracy
# for five genetic-map scales (cM per bp) at two SNP densities.  Run after
# run_smoketest.sh (it reuses the conformed fold targets in smoketest/out/).
#   JAVA=... BCFTOOLS=... PY=python bash c4_imputation/smoketest/map_scale_grid.sh
# Chose the default MAP_CM_PER_BP=1e-6 in 02_impute_c4.sh (best C4A GREx r).
# Caveat: the scale was chosen on the same 111 people it is scored on, so the
# accuracy at the chosen value is mildly optimistic (5 candidate values).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; C4="$(dirname "$HERE")"; O="$HERE/out"; G="$O/mapgrid"; mkdir -p "$G"
JAVA="${JAVA:-java}"; BCFTOOLS="${BCFTOOLS:-bcftools}"; PY="${PY:-python}"
for THIN in 1 3; do for SC in 1e-7 3e-7 1e-6 3e-6 1e-5; do
  pre=()
  for F in 1 2 3 4 5; do
    P="$O/fold${F}_thin${THIN}"; C=$(ls -d "$P".imputed.tmp.*/ | tail -1); Q="$G/f${F}_t${THIN}_$SC"
    "$BCFTOOLS" query -f '%CHROM\t%POS\n' "$P.ref.vcf.gz" | awk -v s="$SC" '{print $1"\t.\t"$2*s"\t"$2}' > "$G/m.map"
    "$JAVA" -jar "$C4/resources/beagle.27Feb25.75f.jar" gt="$C/conform.vcf.gz" ref="$P.ref.vcf.gz" map="$G/m.map" \
      chrom=6:24894177-33890574 out="$Q" ap=true nthreads=4 > /dev/null
    "$PY" "$C4/03_c4_grex.py" --vcf "$Q.vcf.gz" --out-calls "$Q.calls.tsv" --out-summary "$G/s.tsv" --bcftools "$BCFTOOLS" > /dev/null
    cp "$P.truth.tsv" "$Q.truth.tsv"; pre+=("$Q")
  done
  "$PY" "$HERE/evaluate.py" --pairs "${pre[@]}" --label "t${THIN}_map$SC" --out "$G/acc.tsv" > /dev/null
done; done
"$PY" - "$G/acc.tsv" "$C4/results/smoketest_map_scale_grid.tsv" <<'EOF'
import sys, pandas as pd
t = pd.read_csv(sys.argv[1], sep="\t")
t.insert(1, "thin", t.label.str.extract(r"t(\d)_")[0].astype(int)); t.insert(2, "map_cM_per_bp", t.label.str.extract(r"map(.*)$")[0])
t["_s"] = t.map_cM_per_bp.astype(float); t = t.sort_values(["thin", "_s"]).drop(columns=["label", "_s"])
num = t.select_dtypes("float").columns; t[num] = t[num].round(4); t.to_csv(sys.argv[2], sep="\t", index=False)
EOF
echo "wrote $C4/results/smoketest_map_scale_grid.tsv"
