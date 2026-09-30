#!/usr/bin/env python
"""Step 16.7: tracked summary tables and the gate verdicts.

Writes under <root> (only table_*.tsv are tracked by .gitignore):
  table_mostest_build.tsv        sample and phenotype-construction checks (01)
  table_mostest_engine_check.tsv REGENIE vs GENESIS (gate G1)
  table_mostest_summary.tsv      one row per family x scan: R spectrum, null fits,
                                 lambda, tail counts, loci, gene-level lambda
  table_mostest_univariate.tsv   per-measure lambda and hits
  table_mostest_qq.tsv           binned QQ points (figure input)
  table_mostest_loci.tsv         clumped index SNPs (p < 1e-6), family, locus genes
  table_mostest_genes.tsv        genes with MAGMA p < 1e-4 per scan (+ Bonferroni flag)
  table_mostest_set_tests.tsv    step-15b gene sets on each scan's genes.raw
  table_mostest_gates.tsv        G0-G4 with value, threshold, pass
All SNP/gene-level; nothing per-subject (rule 16).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

FAMS = ["slope", "ct", "slopeols"]


def read_clumped(f: Path) -> pd.DataFrame:
    if not f.exists() or f.stat().st_size == 0:
        return pd.DataFrame(columns=["CHR", "BP", "SNP", "P", "TOTAL", "SP2"])
    return pd.read_csv(f, sep=r"\s+", usecols=["CHR", "BP", "SNP", "P", "TOTAL", "SP2"])


def read_genes(f: Path) -> pd.DataFrame:
    return pd.read_csv(f, sep=r"\s+")


def gene_lambda(g: pd.DataFrame) -> float:
    return float(np.median(stats.chi2.isf(g.P.clip(1e-300), 1)) / stats.chi2.ppf(0.5, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = Path(a.root)
    gates = []

    # ---- build (G0)
    b = json.loads((root / "pheno" / "build_summary.json").read_text())
    counts = dict(pd.read_csv(root / "geno" / "step1_counts.tsv", sep="\t", header=None).values)
    bt = pd.DataFrame([dict(item="n_children", value=b["n"]),
                       dict(item="n_analysis_set", value=b["n_analysis_set"]),
                       dict(item="n_step1_array_snps", value=counts["n_step1_snps"]),
                       dict(item="r_mean_slope_blup_vs_global_slope_1lmm",
                            value=round(b["r_mean_slope_blup_vs_global_slope_1lmm"], 4)),
                       dict(item="r_mean_ct_blup_vs_baseline_thickness_1lmm",
                            value=round(b["r_mean_ct_blup_vs_baseline_thickness_1lmm"], 4))]
                      + [dict(item=f"median_pheno_r_{k}", value=round(v["median_pheno_r"], 4))
                         for k, v in b["families"].items() if v["median_pheno_r"] is not None])
    bt.to_csv(root / "table_mostest_build.tsv", sep="\t", index=False)
    gates.append(dict(gate="G0_sample", value=b["n"], threshold=">= 8500", passed=b["n"] >= 8500))

    # ---- engine (G1)
    ec = pd.read_csv(root / "engine_check.tsv", sep="\t")
    ec.to_csv(root / "table_mostest_engine_check.tsv", sep="\t", index=False)
    gates.append(dict(gate="G1_engine_|r_z|_min", value=round(ec.r_z.abs().min(), 4),
                      threshold=">= 0.90 and |dlambda| <= 0.05", passed=bool(ec["pass"].all())))

    # ---- per family
    rows, uni, qq, loci, genes = [], [], [], [], []
    for fam in FAMS:
        s = json.loads((root / f"{fam}_mostest.json").read_text())
        uni.append(pd.read_csv(root / f"univariate_{fam}.tsv", sep="\t"))
        qq.append(pd.read_csv(root / f"qq_{fam}.tsv", sep="\t"))
        for scan in ("real", "perm"):
            sc = fam if scan == "real" else f"{fam}_perm"
            r = dict(family=fam, scan=scan, k=s["k"], n_snps=s["n_snps"], n_median=s["n_median"],
                     R_median_offdiag=round(s["R"]["median_offdiag"], 4),
                     R_eff_dim_pr=round(s["R"]["eff_dim_pr"], 2), R_cond=round(s["R"]["cond"], 1),
                     R_n_floored=s["R"]["n_floored"],
                     R_maxdiff_real_vs_perm=round(s["R"]["max_abs_diff_real_vs_perm"], 4),
                     gamma_shape=round(s["gamma_fit"]["shape"], 4),
                     gamma_scale=round(s["gamma_fit"]["scale"], 4),
                     minp_m_eff=round(s["minp_fit"]["m_eff"], 2))
            r.update({kk: (round(v, 4) if isinstance(v, float) else v) for kk, v in s[scan].items()})
            for thr in ("5e-8", "1e-6"):
                c = read_clumped(root / "clump" / f"{sc}_p{thr}.clumped")
                r[f"n_loci_p{thr}"] = len(c) if (root / "clump" / f"{sc}_p{thr}.clumped").exists() else np.nan
            gf = root / "magma" / "genes" / f"mostest_{sc}.genes.out"
            if gf.exists():
                g = read_genes(gf)
                bonf = 0.05 / len(g)
                r.update(n_genes=len(g), gene_lambda=round(gene_lambda(g), 4),
                         n_genes_bonf=int((g.P < bonf).sum()), n_genes_p1e4=int((g.P < 1e-4).sum()))
                top = g[g.P < 1e-4].assign(family=fam, scan=scan, bonferroni=lambda d: d.P < bonf)
                genes.append(top[["family", "scan", "GENE", "CHR", "START", "STOP", "NSNPS",
                                  "ZSTAT", "P", "bonferroni"]])
            rows.append(r)
            if scan == "real":
                c = read_clumped(root / "clump" / f"{sc}_p1e-6.clumped")
                if len(c):
                    loci.append(c.assign(family=fam, genome_wide=c.P < 5e-8))
    summ = pd.DataFrame(rows)
    summ.to_csv(root / "table_mostest_summary.tsv", sep="\t", index=False)
    pd.concat(uni).to_csv(root / "table_mostest_univariate.tsv", sep="\t", index=False)
    pd.concat(qq).to_csv(root / "table_mostest_qq.tsv", sep="\t", index=False)
    (pd.concat(loci) if loci else pd.DataFrame(columns=["family", "CHR", "BP", "SNP", "P", "TOTAL", "SP2",
                                                        "genome_wide"])) \
        .drop(columns=["SP2"], errors="ignore") \
        .to_csv(root / "table_mostest_loci.tsv", sep="\t", index=False)
    (pd.concat(genes) if genes else pd.DataFrame()).to_csv(root / "table_mostest_genes.tsv",
                                                           sep="\t", index=False)

    # ---- gene sets
    sets = []
    for f in sorted((root / "magma" / "sets").glob("*.gsa.out")):
        scan, setfile = f.name[:-8].split("__")
        d = pd.read_csv(f, sep=r"\s+", comment="#").rename(columns=str.lower)
        sets.append(d.assign(scan=scan, setfile=setfile)[["scan", "setfile", "variable", "ngenes",
                                                          "beta", "se", "p"]])
    if sets:
        pd.concat(sets).to_csv(root / "table_mostest_set_tests.tsv", sep="\t", index=False)

    # ---- gene-level concordance across scans: is the slope signal its own, or
    # the CT signal again?  (Spearman of MAGMA gene Z; perm rows are the floor.)
    gz = {}
    for sc in ["slope", "ct", "slopeols", "slope_perm", "ct_perm"]:
        f = root / "magma" / "genes" / f"mostest_{sc}.genes.out"
        if f.exists():
            gz[sc] = read_genes(f).set_index("GENE").ZSTAT
    pairs = [("slope", "slopeols"), ("slope", "ct"), ("slopeols", "ct"),
             ("slope", "slope_perm"), ("ct", "ct_perm")]
    zc = [dict(scan_a=x, scan_b=y, n_genes=int(pd.concat([gz[x], gz[y]], axis=1).dropna().shape[0]),
               spearman=round(float(pd.concat([gz[x], gz[y]], axis=1).corr("spearman").iloc[0, 1]), 4))
          for x, y in pairs if x in gz and y in gz]
    pd.DataFrame(zc).to_csv(root / "table_mostest_gene_z_corr.tsv", sep="\t", index=False)

    # ---- gates G2-G4
    perm = summ[summ.scan == "perm"].set_index("family"); real = summ[summ.scan == "real"].set_index("family")
    for fam in FAMS:
        lam = perm.loc[fam, "lambda_mostest"]
        gates.append(dict(gate=f"G2a_perm_lambda_{fam}", value=lam, threshold="0.95-1.05",
                          passed=bool(0.95 <= lam <= 1.05)))
        nl = perm.loc[fam, "n_loci_p5e-8"]
        gates.append(dict(gate=f"G2b_perm_loci_{fam}", value=nl, threshold="<= 1",
                          passed=bool(nl <= 1) if pd.notna(nl) else False))
        d = real.loc[fam, "R_maxdiff_real_vs_perm"]
        gates.append(dict(gate=f"G3_R_real_vs_perm_{fam}", value=d, threshold="<= 0.05",
                          passed=bool(d <= 0.05)))
        gl = perm.loc[fam, "gene_lambda"] if "gene_lambda" in perm.columns else np.nan
        gates.append(dict(gate=f"G2c_perm_gene_lambda_{fam}", value=gl, threshold="0.90-1.10",
                          passed=bool(0.90 <= gl <= 1.10) if pd.notna(gl) else False))
    nct = real.loc["ct", "n_loci_p5e-8"]
    gates.append(dict(gate="G4_positive_control_ct_loci", value=nct, threshold=">= 1",
                      passed=bool(nct >= 1) if pd.notna(nct) else False))
    gt = pd.DataFrame(gates)
    gt.to_csv(root / "table_mostest_gates.tsv", sep="\t", index=False)
    print(gt.to_string(index=False))
    cols = ["family", "scan", "lambda_mostest", "lambda_minp", "n_loci_p5e-8", "n_loci_p1e-6",
            "gene_lambda", "n_genes_bonf"]
    print(summ[[c for c in cols if c in summ.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
