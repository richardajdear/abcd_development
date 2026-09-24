"""One table of the current results, read from the committed cluster tables.

Usage (repo root):  python genetic_analysis/build_current_results.py
Output:             genetic_analysis/current_results.tsv   (tracked; summary-level only)

Every number quoted in genetic_analysis/README_HPC.md comes from this table.
Rows cover the primary specification (HCP-MMP, single-LMM slope, 2025 SCZ GWAS)
and the two sensitivity axes that are kept live: the DK parcellation and the
per-region-BLUP-mean slope construction.  `source` names the table each row
was read from, so any row can be traced back to the cluster job that wrote it.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
ROOTS = {"hcp": WORK / "results_70tab_hcp", "dk": WORK / "results_70tab"}
OUT = HERE / "current_results.tsv"
rows: list[dict] = []


def add(section, readout, atlas, construction, trait, method, est, se, p, n,
        source, note=""):
    rows.append(dict(section=section, readout=readout, atlas=atlas,
                     construction=construction, trait=trait, method=method,
                     estimate=est, se=se, p=p, n=n, note=note,
                     source=str(Path(source).relative_to(HERE))))


for atlas, root in ROOTS.items():
    # --- heritability (GCTA GREML, dense imputed GRM, 6,011 unrelated) ----------
    for construction, f in (("1lmm", root / "scan_1lmm/reml_imp_pooled/reml_summary.tsv"),
                            ("perregion", root / "reml_imp_pooled/reml_summary.tsv")):
        t = pd.read_csv(f, sep="\t")
        for r in t.itertuples():
            ph = r.phenotype.replace("_1lmm", "")
            if ph not in ("global_slope", "baseline_thickness"):
                continue
            add("heritability", "GREML h2", atlas, construction, ph, "GCTA",
                r.h2, r.se, r.pval, r.n, f)

    # --- GWAS scan summaries ---------------------------------------------------
    for construction, arm, f in (("1lmm", "EUR", root / "scan_1lmm/assoc_eur/gwas_summary.tsv"),
                                 ("perregion", "EUR", root / "assoc_eur/gwas_summary.tsv"),
                                 ("perregion", "pooled", root / "assoc/gwas_summary.tsv")):
        t = pd.read_csv(f, sep="\t")
        for r in t.itertuples():
            ph = r.phenotype.replace("_1lmm", "")
            add("gwas", f"GWAS {arm}: lambda_GC (estimate), n hits p<5e-8 (note)",
                atlas, construction, ph, "GENESIS", r.lambda_gc, None, r.min_p,
                r.n_mean, f, note=f"hits={r.n_p5e8}; p<1e-5={r.n_p1e5}")

    # --- LDSC h2 z and rg (EUR) --------------------------------------------------
    f = root / "ldsc_1lmm/table_ldsc_panel.tsv"
    t = pd.read_csv(f, sep="\t")
    t = t[t.atlas == atlas]
    for (ph, con), g in t.groupby(["phenotype", "construction"]):
        r0 = g.iloc[0]
        add("ldsc", "LDSC h2 (EUR); note = h2 z", atlas, con, ph, "LDSC",
            r0.h2_obs, r0.h2_obs_se, None, None, f, note=f"h2_z={r0.h2_z}")
        for r in g.itertuples():
            add("ldsc", f"LDSC rg with {r.disorder}", atlas, con, ph, "LDSC",
                r.rg, r.se, r.p, None, f,
                note="uninformative (phenotype h2 z<4)" if r.underpowered == "yes" else "")

    # --- SCZ 2025 PRS, every cell, both constructions -------------------------------
    f = root / "prs_scz2025/table_order_of_operations_scz2025.tsv"
    t = pd.read_csv(f, sep="\t")
    for r in t.itertuples():
        for con in ("1lmm", "perregion"):
            add("prs_scz2025", f"PRS beta (SD/SD): {r.cell}", atlas, con,
                r.phenotype, r.method, getattr(r, f"beta_{con}"),
                getattr(r, f"se_{con}"), getattr(r, f"p_adj_{con}"), r.n, f,
                note=f"score={r.score}; p = threshold-adjusted for C+T")

    # --- min-p permutation (C+T, per-region construction) ------------------------
    f = root / "prs_scz2025/table_minp_permutation.tsv"
    t = pd.read_csv(f, sep="\t")
    for r in t[t.phenotype.isin(["global_slope", "baseline_thickness"])].itertuples():
        add("prs_scz2025", f"C+T min-p permutation ({r.design}, {r.stratum})", atlas,
            "perregion", r.phenotype, "CT", None, None, r.p_perm, r.n, f,
            note=f"{r.n_perm} family-block permutations; trait={r.trait_arm}")

    # --- ancestry strata and PRS x stratum heterogeneity (SBayesRC, per-region) --
    f = root / "prs_scz2025/table_scz2025_strata.tsv"
    t = pd.read_csv(f, sep="\t")
    t = t[(t.phenotype == "global_slope") & (t.method == "SBayesRC")
          & t.trait_arm.isin(["SCZ25_META", "SCZ25_EUR"])]
    for r in t.itertuples():
        add("prs_scz2025_strata", f"{r.analysis} ({r.stratum})", atlas, "perregion",
            "global_slope", f"SBayesRC {r.trait_arm}", r.beta, r.se, r.p, r.n, f)

    # --- within-family (Fulker), SCZ 2025, per-region construction ---------------
    f = root / "prs_scz2025/table_scz2025_family.tsv"
    t = pd.read_csv(f, sep="\t")
    t = t[(t.phenotype == "global_slope") & t.trait_arm.isin(["SCZ25_EUR", "SCZ25_META"])
          & ((t.method != "CT") | (t.threshold == "0p5"))]
    t = t.drop_duplicates(["method", "trait_arm", "stratum", "n_pairs"])
    for r in t.itertuples():
        add("prs_scz2025_family", f"within-family beta_W ({r.stratum}, {r.n_pairs} pairs)",
            atlas, "perregion", "global_slope", f"{r.method} {r.trait_arm}",
            r.beta_within, r.se_within, r.p_within, r.n, f,
            note=f"beta_B={r.beta_between:.3f} (p {r.p_between:.3g}); p_diff={r.p_diff:.2f}")

    # --- control panel (PGC3-era SCZ, MDD, ASD, ALZ +/- APOE, EA), both constructions
    f = root / "prs_final_1lmm/table_order_of_operations_all.tsv"
    t = pd.read_csv(f, sep="\t")
    for r in t.itertuples():
        for con in ("1lmm", "perregion"):
            add("prs_panel", f"PRS beta (SD/SD): {r.trait_arm} -> {r.stratum}", atlas,
                con, r.phenotype, r.method, getattr(r, f"beta_{con}"),
                getattr(r, f"se_{con}"), getattr(r, f"p_adj_{con}"), r.n, f,
                note=f"score={r.score}")

    # --- MAGMA panel: reverse / forward gene-property and gene sets (EUR) ---------
    f = root / "magma_panel/table_magma_panel.tsv"
    t = pd.read_csv(f, sep="\t")
    for r in t.itertuples():
        add("magma", f"MAGMA {r.kind}: {r.disorder_result} / {r.variable}", atlas,
            r.construction, r.phenotype, "MAGMA", r.beta, r.se, r.p, r.n_genes, f,
            note="n = genes")

    # --- AHBA / ahba_pls H3 gene-property (per-region construction) --------------
    f = root / "magma_ahba_pls_h3/table_h3.tsv"
    t = pd.read_csv(f, sep="\t")
    for r in t.itertuples():
        add("ahba_h3", f"MAGMA gene-property ({r.model}): {r.variable}", atlas,
            "perregion", r.phenotype, "MAGMA", r.beta, r.se, r.p, r.n_genes, f,
            note="n = genes")

out = pd.DataFrame(rows)
out.to_csv(OUT, sep="\t", index=False, float_format="%.5g")
print(f"{OUT.relative_to(HERE.parent)}: {len(out)} rows")
