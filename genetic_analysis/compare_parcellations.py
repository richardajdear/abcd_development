#!/usr/bin/env python
"""Side-by-side DK vs HCP-MMP readouts from the two results roots.

Reads whatever summary tables exist under work/results_70tab (DK) and
work/results_70tab_hcp (HCP) and writes one long table, one row per
(readout, phenotype, arm, method) with a DK and an HCP column, to
work/results_70tab_hcp/compare/table_dk_vs_hcp.tsv (tracked).  Missing tables
are reported as blank, so this can be run at any stage.

Only `global_slope` and `baseline_thickness` are comparable across atlases:
the slope PCs are atlas-specific decompositions (DK<->HCP |r| 0.88 / 0.79 /
0.51 for PC1-3, PC1 with a sign flip), so their rows are kept but flagged
`comparable=no`.  See README_HPC.md 8, 2026-09-15.

    python genetic_analysis/compare_parcellations.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ROOTS = {"dk": REPO / "genetic_analysis/work/results_70tab",
         "hcp": REPO / "genetic_analysis/work/results_70tab_hcp"}
PHENOS = ["global_slope", "baseline_thickness", "slope_PC1", "slope_PC2", "slope_PC3"]
COMPARABLE = {"global_slope", "baseline_thickness"}
MATCHED = {"SCZ_pooled": "full", "MDD_pooled": "full", "SCZ_eur": "EUR", "MDD_eur": "EUR",
           "ASD": "EUR", "ALZ": "EUR", "ALZ_noAPOE": "EUR", "ALZ_IGAP": "EUR",
           "ALZ_IGAP_noAPOE": "EUR", "EA": "EUR"}


def rd(p: Path, **kw) -> pd.DataFrame | None:
    return pd.read_csv(p, sep="\t", **kw) if p.exists() else None


def collect(root: Path) -> list[dict]:
    rows: list[dict] = []
    def add(readout, pheno, arm, method, value, se=None, p=None, n=None):
        rows.append(dict(readout=readout, phenotype=pheno, arm=arm, method=method,
                         value=value, se=se, p=p, n=n))
    # null models
    for arm, d in (("pooled", "nullmodel"), ("EUR", "nullmodel_eur")):
        for ph in PHENOS:
            t = rd(root / d / f"{ph}_null_summary.tsv")
            if t is not None:
                r = t.iloc[0]; add("null prop_kin", ph, arm, "GENESIS", r.prop_kin, n=r.n)
    # GWAS
    for arm, d in (("pooled", "assoc"), ("EUR", "assoc_eur")):
        t = rd(root / d / "gwas_summary.tsv")
        if t is not None:
            for _, r in t.iterrows():
                if r.n_chr == 22:
                    add("GWAS lambda_GC", r.phenotype, arm, "GENESIS", r.lambda_gc, n=r.n_mean)
                    add("GWAS hits p<5e-8", r.phenotype, arm, "GENESIS", r.n_p5e8, n=r.n_mean)
    # GREML
    t = rd(root / "reml_imp_pooled" / "reml_summary.tsv")
    if t is not None:
        for _, r in t.iterrows():
            add("GREML h2", r.phenotype, "unrelated", "GCTA", r.h2, se=r.se, p=r.pval, n=r.n)
    # PRS, matched cells; C+T at its best threshold; pooled arm uses zanc
    t = rd(root / "prs_final" / "table_main.tsv")
    if t is not None:
        t = t[t.matched == "yes"]
        for (arm, meth, ph), g in t.groupby(["trait_arm", "method", "phenotype"]):
            g = g[g.score == ("zanc" if MATCHED[arm] == "full" else "raw")]
            if g.empty: continue
            r = g.sort_values("p").iloc[0]
            add(f"PRS beta {arm}", ph, MATCHED[arm], meth, r.beta, se=r.se, p=r.p_adj, n=r.n)
    # min-p permutation (C+T), matched design
    t = rd(root / "prs_final" / "table_minp_permutation.tsv")
    if t is not None:
        for _, r in t[t.design == "matched"].iterrows():
            add(f"PRS min-p perm {r.trait_arm}", r.phenotype, r.stratum, "CT", r.p_perm, n=r.n)
    # LDSC
    t = rd(root / "ldsc_eur" / "ldsc_rg_summary.tsv")
    if t is not None:
        for _, r in t.iterrows():
            add("LDSC h2_z", r.phenotype, "EUR", "LDSC", r.h2_z)
            add(f"LDSC rg {r.disorder}", r.phenotype, "EUR", "LDSC", r.rg, se=r.se, p=r.p)
    # prioritised gene sets (marginal)
    for tag, setname in (("sczprio", "SCZ_locus_pool"), ("mddhc", "MDD_pool")):
        t = rd(root / "magma_prio_eur" / f"{tag}_gsa_summary.tsv")
        if t is not None:
            g = t[(t.model == "marginal") & (t["set"] == setname)]
            for _, r in g.iterrows():
                add(f"MAGMA {setname}", r.trait, "EUR", "MAGMA", r.beta, se=r.se, p=r.p, n=r.n_genes)
    # ahba_pls H3 (step 8): signature weights conditioned on AHBA C3
    t = rd(root / "magma_ahba_pls_h3" / "table_h3.tsv")
    if t is not None:
        for _, r in t[(t.model == "condC3") & (t.variable != "AHBA_C3")].iterrows():
            add(f"MAGMA H3 {r.variable} | AHBA_C3 (beta_std)", r.phenotype, "EUR", "MAGMA", r.beta_std, se=r.se, p=r.p, n=r.n_genes)
        for _, r in t[(t.model == "marginal")].iterrows():
            add(f"MAGMA H3 {r.variable} marginal (beta_std)", r.phenotype, "EUR", "MAGMA", r.beta_std, se=r.se, p=r.p, n=r.n_genes)
    return rows


def main() -> int:
    frames = []
    for tag, root in ROOTS.items():
        f = pd.DataFrame(collect(root))
        if f.empty:
            print(f"{tag}: no summary tables under {root}"); continue
        f["parc"] = tag; frames.append(f)
    if not frames:
        return 1
    long = pd.concat(frames, ignore_index=True)
    key = ["readout", "phenotype", "arm", "method"]
    wide = long.pivot_table(index=key, columns="parc", values=["value", "se", "p", "n"], aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    wide["comparable"] = wide.phenotype.map(lambda x: "yes" if x in COMPARABLE else "no (atlas-specific PC)")
    order = [c for c in ["readout", "phenotype", "arm", "method", "comparable",
                         "value_dk", "value_hcp", "se_dk", "se_hcp", "p_dk", "p_hcp", "n_dk", "n_hcp"] if c in wide.columns]
    wide = wide[order].sort_values(["readout", "phenotype", "arm", "method"])
    out = ROOTS["hcp"] / "compare" / "table_dk_vs_hcp.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(out, sep="\t", index=False, float_format="%.5g")
    print(f"wrote {out} ({len(wide)} rows; parcellations present: {sorted(long.parc.unique())})")
    show = wide[wide.phenotype.isin(COMPARABLE)]
    pd.set_option("display.width", 220); pd.set_option("display.max_rows", 200)
    print(show.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
