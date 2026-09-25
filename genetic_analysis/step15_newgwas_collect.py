#!/usr/bin/env python
"""Step 15 collector: new discovery GWAS (BIP_eur, BIP_pooled, ADHD, INT) x four
PRS methods, both atlases, both phenotype constructions.

Per results root (results_70tab = DK, results_70tab_hcp = HCP-MMP), from
prs_newgwas/ (per-region-BLUP-mean phenotypes) and prs_newgwas_1lmm/ (single LMM):
  prs_newgwas/table_newgwas_main.tsv      every method x arm x stratum x score,
                                          best C+T threshold with its Bonferroni
                                          p_adj (as collect_final.py), matched flag
  prs_newgwas/table_newgwas_family.tsv    Fulker between/within rows
  prs_newgwas/table_newgwas_strata.tsv    per-cluster fits + PRS x stratum LRT
  prs_newgwas/table_order_of_operations_newgwas.tsv
        matched cells, per-region vs single LMM, global_slope and
        baseline_thickness, the COLUMNS OF prs_final_1lmm/
        table_order_of_operations_all.tsv -- which is what
        fig1_prep_prs.py globs for (trait_arm BIP_eur / BIP_pooled / ADHD).

Matching (README_HPC.md 4 rule 4): BIP_pooled is a multi-ancestry meta-analysis
-> pooled target with the within-cluster-standardised score; every other arm is
European -> the 4,308-child EUR anchor, raw score.
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ROOTS = {"dk": REPO / "genetic_analysis/work/results_70tab", "hcp": REPO / "genetic_analysis/work/results_70tab_hcp"}
METHODS = ["CT", "PRSCS", "SBayesR", "SBayesRC"]
ARMS = {  # arm -> (display trait, discovery ancestry, description)
    "BIP_eur":    ("BIP",  "EUR",   "PGC bip2024 EUR (O'Connell 2025)"),
    "BIP_pooled": ("BIP",  "multi", "PGC bip2024 multi-ancestry (EUR+AFR+EAS+LAT)"),
    "ADHD":       ("ADHD", "EUR",   "Demontis 2023 iPSYCH+deCODE+PGC"),
    "INT":        ("INT",  "EUR",   "Savage & Jansen 2018 intelligence"),
}
MATCHED = {"multi": "full", "EUR": "EUR"}
HEADLINE = ["global_slope", "baseline_thickness"]


def rd(p: Path) -> pd.DataFrame | None:
    return pd.read_csv(p, sep="\t") if p.exists() else None


def best(t: pd.DataFrame, stratum: str) -> pd.DataFrame:
    s = t[t.stratum == stratum].dropna(subset=["p"])
    return s.sort_values("p").groupby("phenotype", as_index=False).head(1)


def main_tables(root: Path) -> None:
    d = root / "prs_newgwas"
    main, fam, strata = [], [], []
    for meth in METHODS:
        for arm, (trait, anc, desc) in ARMS.items():
            for ver in ("raw", "zanc"):
                t = rd(d / f"assoc_{meth}_{arm}{'' if ver == 'raw' else '_zanc'}.tsv")
                if t is not None:
                    for st in ("full", "EUR"):
                        if ver == "zanc" and st == "EUR":
                            continue
                        for _, r in best(t, st).iterrows():
                            main.append(dict(method=meth, trait=trait, trait_arm=arm, discovery_ancestry=anc,
                                             target_stratum=st, matched="yes" if MATCHED[anc] == st else "no",
                                             score=ver, discovery_gwas=desc, phenotype=r.phenotype, threshold=r.threshold,
                                             n=r.n, n_families=r.n_families, n_snps=r.get("n_snps", ""),
                                             beta=r.beta, se=r.se, p=r.p, p_adj=r.p_adj))
                f = rd(d / f"fam_{meth}_{arm}{'' if ver == 'raw' else '_zanc'}.tsv")
                if f is not None:
                    f.insert(0, "method", meth); f.insert(1, "trait_arm", arm); f.insert(2, "score", ver); fam.append(f)
            s = rd(d / f"strata_{meth}_{arm}.tsv")
            if s is not None:
                s.insert(0, "method", meth); s.insert(1, "trait_arm", arm); strata.append(s)
    m = pd.DataFrame(main)
    m.to_csv(d / "table_newgwas_main.tsv", sep="\t", index=False)
    if fam: pd.concat(fam).to_csv(d / "table_newgwas_family.tsv", sep="\t", index=False)
    if strata: pd.concat(strata).to_csv(d / "table_newgwas_strata.tsv", sep="\t", index=False)
    return m


def orderops(atlas: str, root: Path, m: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for meth in METHODS:
        for arm, (_, anc, _) in ARMS.items():
            st = MATCHED[anc]; sv = "zanc" if st == "full" else "raw"
            t = rd(root / "prs_newgwas_1lmm" / f"assoc_{meth}_{arm}{'' if sv == 'raw' else '_zanc'}.tsv")
            if t is None:
                continue
            for ph in HEADLINE:
                one = t[(t.phenotype == f"{ph}_1lmm") & (t.stratum == st)].sort_values("p")
                b = m[(m.trait_arm == arm) & (m.method == meth) & (m.score == sv) & (m.target_stratum == st) & (m.phenotype == ph)]
                if one.empty or b.empty:
                    continue
                one, b = one.iloc[0], b.iloc[0]
                rows.append(dict(atlas=atlas, phenotype=ph, trait_arm=arm, stratum=st, method=meth, score=sv,
                                 beta_perregion=b.beta, se_perregion=b.se, p_adj_perregion=b.p_adj,
                                 beta_1lmm=one.beta, se_1lmm=one.se, p_adj_1lmm=one.p_adj, n=int(one.n)))
    o = pd.DataFrame(rows)
    o.to_csv(root / "prs_newgwas" / "table_order_of_operations_newgwas.tsv", sep="\t", index=False)
    return o


def main() -> int:
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
    for atlas, root in ROOTS.items():
        if not (root / "prs_newgwas").exists():
            print(f"[{atlas}] no prs_newgwas/"); continue
        m = main_tables(root); o = orderops(atlas, root, m)
        print(f"\n[{atlas}] main {len(m)} rows; order-of-operations {len(o)} rows")
        if len(o):
            f = o.copy()
            for c in ["beta_perregion", "se_perregion", "beta_1lmm", "se_1lmm"]: f[c] = f[c].round(4)
            for c in ["p_adj_perregion", "p_adj_1lmm"]: f[c] = f[c].map(lambda v: f"{v:.3g}")
            print(f.drop(columns=["atlas", "score", "n"]).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
