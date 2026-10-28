"""D1 step 9 report tables, from whatever scan_geom outputs exist.

Writes under scan_geom/:
  table_gwas.tsv   arm x phenotype: n, lambda_GC, min P, #P<5e-8, #P<1e-6
  table_magma.tsv  arm x phenotype: genes, lambda of gene P, gene-Z Spearman with
                   slope_thk, and (slope_thk only) gene-Z Pearson with the
                   global_slope_1lmm reference; top-10 genes
  plus the GREML table if present.

  python summarise.py <results_70tab_hcp dir> <NCBI37.3.gene.loc>
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats

R, LOC = sys.argv[1:3]; W = f"{R}/scan_geom"
sym = pd.read_csv(LOC, sep="\t", header=None, usecols=[0, 5], names=["GENE", "SYMBOL"]).set_index("GENE").SYMBOL
TAGS = ["thk", "area", "vol", "thk_area", "thk_vol_area"]
PHENOS = [f"slope_{t}" for t in TAGS] + [f"int_{t}" for t in TAGS]
ARMS = {"eur": ("assoc_eur", "magma_eur/genes", f"{R}/magma_panel/genes"),
        "pooled": ("assoc", "magma_pooled/genes", f"{R}/magma_pooled/genes")}
CHI = stats.chi2.ppf(0.5, 1)

rows = []
for arm, (ad, _, _) in ARMS.items():
    for p in PHENOS:
        f = f"{W}/{ad}/{p}.sumstats.tsv.gz"
        if not os.path.exists(f): continue
        d = pd.read_csv(f, sep="\t", usecols=["N", "BETA", "SE", "P"])
        rows.append(dict(arm=arm, phenotype=p, n_snps=len(d), n_mean=round(d.N.mean()),
                         lambda_gc=round(np.median((d.BETA / d.SE) ** 2) / CHI, 4), min_p=d.P.min(),
                         n_p5e8=int((d.P < 5e-8).sum()), n_p1e6=int((d.P < 1e-6).sum())))
if rows:
    t = pd.DataFrame(rows); t.to_csv(f"{W}/table_gwas.tsv", sep="\t", index=False); print(t.to_string(index=False))

rows = []
for arm, (_, gd, refd) in ARMS.items():
    rd = lambda p: pd.read_csv(f"{W}/{gd}/{p}.genes.out", sep=r"\s+").set_index("GENE")
    if not os.path.exists(f"{W}/{gd}/slope_thk.genes.out"): continue
    thk = rd("slope_thk")
    for p in PHENOS:
        if not os.path.exists(f"{W}/{gd}/{p}.genes.out"): continue
        g = rd(p); j = g.join(thk[["ZSTAT"]], rsuffix="_thk", how="inner")
        row = dict(arm=arm, phenotype=p, genes=len(g),
                   lambda_geneP=round(np.median(stats.chi2.isf(g.P.clip(lower=1e-300), 1)) / CHI, 4),
                   spearman_Z_vs_slope_thk=round(stats.spearmanr(j.ZSTAT, j.ZSTAT_thk)[0], 4),
                   top10=",".join(f"{sym.get(i, i)}({pv:.1e})" for i, pv in g.P.nsmallest(10).items()))
        if p in ("slope_thk", "int_thk"):
            ref = pd.read_csv(f"{refd}/{'global_slope_1lmm' if p == 'slope_thk' else 'baseline_thickness_1lmm'}.genes.out",
                              sep=r"\s+").set_index("GENE")
            jr = g.join(ref[["ZSTAT"]], rsuffix="_ref", how="inner")
            row["pearson_Z_vs_1lmm_ref"] = round(np.corrcoef(jr.ZSTAT, jr.ZSTAT_ref)[0, 1], 4)
        rows.append(row)
if rows:
    t = pd.DataFrame(rows); t.to_csv(f"{W}/table_magma.tsv", sep="\t", index=False)
    print(t.drop(columns="top10").to_string(index=False))
    for r in rows: print(f"{r['arm']:6s} {r['phenotype']:20s} {r['top10']}")
f = f"{W}/reml_imp_pooled/reml_summary.tsv"
if os.path.exists(f): print(pd.read_csv(f, sep="\t").to_string(index=False))
