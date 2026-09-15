#!/usr/bin/env python
"""README_HPC.md 7 item 1: one row per (readout, phenotype, arm, method) with
estimate, SE, p, n -- and the 6.0-vintage legacy value alongside, so every
row of 5 has a partner.  Writes work/results_70tab<suffix>/summary_70tab.tsv
(tracked).  --parc hcp writes the HCP root's table; the legacy column is DK
6.0 in both cases (there is no HCP legacy), flagged in `legacy_note`.

    python genetic_analysis/build_summary_70tab.py [--parc dsk|hcp]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
LEG = REPO / "legacy"
V2 = LEG / "hpc_v2/work/results_v2"
V1 = LEG / "hpc/work/results"
V3 = LEG / "hpc_v3/prs_tables"
MATCHED = {"SCZ_pooled": "full", "MDD_pooled": "full", "SCZ_eur": "EUR", "MDD_eur": "EUR",
           "ASD": "EUR", "ALZ": "EUR", "ALZ_noAPOE": "EUR", "ALZ_IGAP": "EUR",
           "ALZ_IGAP_noAPOE": "EUR", "EA": "EUR"}
PHENOS = ["global_slope", "baseline_thickness", "slope_PC1", "slope_PC2", "slope_PC3"]


def rd(p: Path) -> pd.DataFrame | None:
    return pd.read_csv(p, sep="\t") if p.exists() else None


def rows_from(root: Path, legacy: bool) -> list[dict]:
    """Same extraction for the new root and for the legacy roots."""
    out: list[dict] = []
    def add(readout, pheno, arm, method, est, se=None, p=None, n=None):
        out.append(dict(readout=readout, phenotype=pheno, arm=arm, method=method,
                        estimate=est, se=se, p=p, n=n))
    # null models
    for arm, d in (("pooled", "nullmodel"), ("EUR", "nullmodel_eur")):
        for ph in PHENOS:
            t = rd(root / d / f"{ph}_null_summary.tsv")
            if t is not None:
                r = t.iloc[0]; add("null model prop_kin", ph, arm, "GENESIS", r.prop_kin, n=r.n)
    # GWAS
    for arm, d in (("pooled", "assoc"), ("EUR", "assoc_eur")):
        t = rd(root / d / "gwas_summary.tsv")
        if t is not None:
            for _, r in t[t.n_chr == 22].iterrows():
                add("GWAS lambda_GC", r.phenotype, arm, "GENESIS", r.lambda_gc, n=r.n_mean)
                add("GWAS hits p<5e-8", r.phenotype, arm, "GENESIS", r.n_p5e8, n=r.n_mean)
    # GREML: legacy lives in v1's results/reml_imp_pooled
    t = rd((V1 / "reml_imp_pooled" if legacy else root / "reml_imp_pooled") / "reml_summary.tsv")
    if t is not None:
        for _, r in t.iterrows():
            add("GREML h2 (dense imputed GRM, PC-AiR unrelated)", r.phenotype, "unrelated", "GCTA",
                r.h2, se=r.se, p=r.pval, n=r.n)
    # PRS matched cells
    t = rd(root / "prs_final" / "table_main.tsv")
    if t is not None:
        t = t[t.matched == "yes"]
        for (arm, meth, ph), g in t.groupby(["trait_arm", "method", "phenotype"]):
            g = g[g.score == ("zanc" if MATCHED[arm] == "full" else "raw")]
            if g.empty: continue
            r = g.sort_values("p").iloc[0]
            add(f"PRS {arm} -> {MATCHED[arm]} (beta SD/SD, p_adj)", ph, MATCHED[arm], meth,
                r.beta, se=r.se, p=r.p_adj, n=r.n)
    # within-family (Fulker), C+T best + single-score methods
    t = rd(root / "prs_final" / "table_family.tsv")
    if t is not None:
        t = t[(t.score == "raw") & t.trait_arm.isin(MATCHED)]
        t = t[t.apply(lambda r: MATCHED[r.trait_arm] == r.stratum, axis=1)]
        for (arm, meth, ph), g in t.groupby(["trait_arm", "method", "phenotype"]):
            r = g.sort_values("p_between").iloc[0]
            add(f"PRS within-family beta_W {arm}", ph, r.stratum, meth, r.beta_within, se=r.se_within, p=r.p_within, n=r.n_pairs)
            add(f"PRS between-family beta_B {arm}", ph, r.stratum, meth, r.beta_between, se=r.se_between, p=r.p_between, n=r.n_pairs)
    # min-p permutation: new = matched design; legacy = hpc_v3 table (legacy-like design)
    if legacy:
        t = rd(V3 / "prs_minp_permutation.tsv")
        if t is not None:
            for _, r in t.iterrows():
                add(f"PRS C+T min-p permutation {r.disorder} (legacy-like design)", r.phenotype, r.stratum, "CT", r.p_perm, n=r.n)
    else:
        t = rd(root / "prs_final" / "table_minp_permutation.tsv")
        if t is not None:
            for _, r in t.iterrows():
                add(f"PRS C+T min-p permutation {r.trait_arm} ({r.design} design)", r.phenotype, r.stratum, "CT", r.p_perm, n=r.n)
    # LDSC
    t = rd(root / "ldsc_eur" / "ldsc_rg_summary.tsv")
    if t is not None:
        # h2 is estimated on the SNPs shared with each disorder file, so it
        # differs slightly by disorder row; keep the disorder in the readout
        # name or the legacy merge cross-multiplies the rows.
        for _, r in t.iterrows():
            add(f"LDSC h2_obs (SNPs shared with {r.disorder})", r.phenotype, "EUR", "LDSC", r.h2_obs, se=r.h2_obs_se)
            add(f"LDSC h2_z (SNPs shared with {r.disorder})", r.phenotype, "EUR", "LDSC", r.h2_z)
            add(f"LDSC rg with {r.disorder}", r.phenotype, "EUR", "LDSC", r.rg, se=r.se, p=r.p)
    # prioritised gene sets, marginal
    for tag, setname in (("sczprio", "SCZ_locus_pool"), ("mddhc", "MDD_pool")):
        t = rd(root / "magma_prio_eur" / f"{tag}_gsa_summary.tsv")
        if t is not None:
            g = t[(t.model == "marginal") & (t["set"] == setname)]
            for _, r in g.iterrows():
                add(f"MAGMA gene-set {setname} (marginal beta)", r.trait, "EUR", "MAGMA", r.beta, se=r.se, p=r.p, n=r.n_genes)
    # ahba_pls H3 (step 8): signature weights conditioned on AHBA C3
    t = rd(root / "magma_ahba_pls_h3" / "table_h3.tsv")
    if t is not None:
        for _, r in t[(t.model == "condC3") & (t.variable != "AHBA_C3")].iterrows():
            add(f"MAGMA H3 {r.variable} | AHBA_C3 (beta_std)", r.phenotype, "EUR", "MAGMA", r.beta_std, se=r.se, p=r.p, n=r.n_genes)
        for _, r in t[(t.model == "marginal")].iterrows():
            add(f"MAGMA H3 {r.variable} marginal (beta_std)", r.phenotype, "EUR", "MAGMA", r.beta_std, se=r.se, p=r.p, n=r.n_genes)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parc", default="dsk", choices=["dsk", "hcp"])
    a = ap.parse_args()
    root = REPO / "genetic_analysis/work" / ("results_70tab" if a.parc == "dsk" else "results_70tab_hcp")
    new = pd.DataFrame(rows_from(root, legacy=False))
    leg = pd.DataFrame(rows_from(V2, legacy=True))
    if new.empty:
        print(f"no summaries under {root}"); return 1
    key = ["readout", "phenotype", "arm", "method"]
    leg = leg.rename(columns={"estimate": "legacy_estimate", "se": "legacy_se", "p": "legacy_p", "n": "legacy_n"})
    # min-p readout names differ by design on purpose; pair the matched design with the legacy-like benchmark
    leg["readout"] = leg.readout.str.replace(" (legacy-like design)", " (matched design)", regex=False)
    tab = new.merge(leg, on=key, how="left")
    tab["legacy_note"] = "6.0-vintage DK (n 8,082 / EUR 4,116)" + ("; no HCP legacy exists" if a.parc == "hcp" else "")
    tab.loc[tab.readout.str.contains("min-p"), "legacy_note"] += "; legacy min-p was the ancestry-MISMATCHED design"
    tab = tab.sort_values(key)
    out = root / "summary_70tab.tsv"
    tab.to_csv(out, sep="\t", index=False, float_format="%.5g")
    unpaired = tab.legacy_estimate.isna().sum()
    print(f"wrote {out}: {len(tab)} rows, {len(tab) - unpaired} with a legacy partner, {unpaired} without")
    print(tab[tab.phenotype == "global_slope"][key + ["estimate", "se", "p", "n", "legacy_estimate", "legacy_p"]]
          .to_string(index=False, max_rows=60))
    return 0


if __name__ == "__main__":
    sys.exit(main())
