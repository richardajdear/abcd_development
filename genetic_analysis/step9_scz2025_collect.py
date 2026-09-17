#!/usr/bin/env python
"""Assemble the step-9 tables: the Nature 2025 SCZ GWAS, four methods + PRS-CSx,
both ancestry arms, both parcellations, read against the PGC3 result.

Per results root (work/results_70tab = DK, work/results_70tab_hcp = HCP-MMP),
from prs_scz2025/{assoc,fam,strata}_*.tsv:
  table_scz2025_main.tsv     headline phenotypes; one row per method x trait-arm
                             x target stratum x score version (raw / zanc), best
                             C+T threshold by p with p_adj already Bonferroni'd
                             over the 8 thresholds (as collect_final.py does)
  table_scz2025_all.tsv      every phenotype
  table_scz2025_family.tsv   Fulker between/within rows
  table_scz2025_strata.tsv   per-cluster fits, pooled-with-stratum-FE, and the
                             PRS x stratum heterogeneity LRT
Then, in results_70tab_hcp/compare/:
  table_scz2025_vs_pgc3_dk_vs_hcp.tsv   global_slope, the DESIGN cells below,
                             PGC3 and 2025 side by side, DK and HCP side by side.

DESIGN LABELS (README_HPC.md 4 rule 4 extended to multiple discovery
ancestries).  A cell is one (discovery arm -> target stratum):
  primary          rule-4 matched: EUR GWAS -> EUR target; multi-ancestry GWAS
                   -> pooled target (score standardised within cluster, zanc);
                   ancestry-matched composite -> pooled target
  secondary        multi-ancestry GWAS -> EUR target: less power, not confounded
  confounded       EUR (or any single-ancestry) GWAS -> pooled target: the score
                   tracks ancestry and so does the phenotype; emitted, flagged,
                   never read as a result (the false ASD "hit" of rule 4)
  strata_only      AFR / EAS single-ancestry weights: meaningful only inside
                   the matching cluster (table_scz2025_strata.tsv), reported
                   here for completeness
"""
from __future__ import annotations

import sys
from collections import OrderedDict
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ROOTS = {"dk": REPO / "genetic_analysis/work/results_70tab",
         "hcp": REPO / "genetic_analysis/work/results_70tab_hcp"}
METHODS = ["CT", "PRSCS", "SBayesR", "SBayesRC", "PRSCSX"]
HEADLINE = ["global_slope", "baseline_thickness"]
# trait-arm -> (discovery ancestry, description)
ARMS = OrderedDict([
    ("SCZ25_EUR",     ("EUR",     "2025 SCZ EUR (MVP+PGC3+AoU+FinnGen)")),
    ("SCZ25_AFR",     ("AFR",     "2025 SCZ AFR (MVP+AoU+GPC+MGS)")),
    ("SCZ25_EAS",     ("EAS",     "2025 SCZ EAS (PGC3 asian)")),
    ("SCZ25_META",    ("multi",   "2025 SCZ AFR+EUR+EAS meta")),
    ("SCZ25_MATCHED", ("matched", "per-child ancestry-matched (EURlike<-EUR, cluster1<-AFR, cluster3<-EAS, cluster2<-META), z within cluster")),
])


def design(anc: str, stratum: str, meth: str, arm: str, score: str = "raw") -> str:
    if anc == "EUR":
        if stratum == "EUR":
            return "primary"
        # EUR weights on the pooled target: the RAW score tracks ancestry and is
        # confounded (rule 4); z-scored within cluster the between-cluster
        # component is gone and it is a legitimate secondary pooled design
        return "secondary (EUR weights, z within cluster)" if score == "zanc" else "confounded"
    if anc == "multi":
        return "primary" if stratum == "full" else "secondary"
    if anc == "matched":
        # PRS-CSx per-population posteriors live under PRSCSX/SCZ25_<POP>; the
        # composite is z within cluster, so the pooled stratum is its home
        return "primary" if stratum == "full" else "n/a (composite read in EUR only = EUR weights)"
    return "strata_only"


def rd(p: Path) -> pd.DataFrame | None:
    return pd.read_csv(p, sep="\t") if p.exists() else None


def best_per_pheno(t: pd.DataFrame, stratum: str) -> pd.DataFrame:
    s = t[t.stratum == stratum].dropna(subset=["p"])
    if s.empty:
        return s
    return s.sort_values("p").groupby("phenotype", as_index=False).head(1)


def collect(root: Path) -> dict[str, pd.DataFrame]:
    fin = root / "prs_scz2025"
    main, allrows, fam, strata = [], [], [], []
    for meth in METHODS:
        for arm, (anc, desc) in ARMS.items():
            for ver in ("raw", "zanc"):
                suf = "" if ver == "raw" else "_zanc"
                t = rd(fin / f"assoc_{meth}_{arm}{suf}.tsv")
                if t is None:
                    continue
                for stratum in ("full", "EUR"):
                    if ver == "zanc" and stratum == "EUR":
                        continue          # one cluster: z-score == raw
                    for _, r in best_per_pheno(t, stratum).iterrows():
                        rec = OrderedDict([
                            ("method", meth), ("trait_arm", arm), ("discovery_ancestry", anc),
                            ("target_stratum", stratum), ("design", design(anc, stratum, meth, arm, ver)),
                            ("score", ver), ("discovery_gwas", desc), ("phenotype", r.phenotype),
                            ("threshold", r.threshold), ("n", r.n), ("n_families", r.n_families),
                            ("n_snps", r.get("n_snps", "")), ("beta", r.beta), ("se", r.se),
                            ("p", r.p), ("p_adj", r.get("p_adj", r.p)),
                        ])
                        allrows.append(rec)
                        if r.phenotype in HEADLINE:
                            main.append(rec)
                f = rd(fin / f"fam_{meth}_{arm}{suf}.tsv")
                if f is not None:
                    f.insert(0, "method", meth); f.insert(1, "trait_arm", arm); f.insert(2, "score", ver)
                    fam.append(f)
            s = rd(fin / f"strata_{meth}_{arm}.tsv")
            if s is not None:
                s.insert(0, "method", meth); s.insert(1, "trait_arm", arm)
                s.insert(2, "discovery_ancestry", anc)
                strata.append(s)
    out = {"main": pd.DataFrame(main), "all": pd.DataFrame(allrows),
           "family": pd.concat(fam, ignore_index=True) if fam else pd.DataFrame(),
           "strata": pd.concat(strata, ignore_index=True) if strata else pd.DataFrame()}
    for k, t in out.items():
        if not t.empty:
            t.to_csv(fin / f"table_scz2025_{k}.tsv", sep="\t", index=False)
    return out


def pgc3_rows(root: Path) -> pd.DataFrame:
    """The PGC3 counterparts from step 6: SCZ_eur -> EUR raw; SCZ_pooled -> full zanc."""
    t = rd(root / "prs_final" / "table_main.tsv")
    if t is None:
        return pd.DataFrame()
    t = t[(t.phenotype == "global_slope") & (t.trait.eq("SCZ"))]
    keep = t[((t.trait_arm == "SCZ_eur") & (t.target_stratum == "EUR") & (t.score == "raw")) |
             ((t.trait_arm == "SCZ_pooled") & (t.target_stratum == "full") & (t.score == "zanc")) |
             ((t.trait_arm == "SCZ_pooled") & (t.target_stratum == "EUR") & (t.score == "raw"))].copy()
    def cell(r):
        if r.trait_arm == "SCZ_eur": return "EUR GWAS -> EUR"
        return "multi GWAS -> pooled (zanc)" if r.target_stratum == "full" else "multi GWAS -> EUR (secondary)"
    keep["cell"] = keep.apply(cell, axis=1)
    return keep[["cell", "method", "beta", "se", "p_adj", "n"]]


def new_rows(main: pd.DataFrame) -> pd.DataFrame:
    if main.empty:
        return main
    t = main[main.phenotype == "global_slope"].copy()
    sel = t[((t.trait_arm == "SCZ25_EUR") & (t.target_stratum == "EUR") & (t.score == "raw")) |
            ((t.trait_arm == "SCZ25_META") & (t.target_stratum == "full") & (t.score == "zanc")) |
            ((t.trait_arm == "SCZ25_META") & (t.target_stratum == "EUR") & (t.score == "raw")) |
            ((t.trait_arm == "SCZ25_MATCHED") & (t.target_stratum == "full") & (t.score == "raw"))].copy()
    def cell(r):
        if r.trait_arm == "SCZ25_EUR": return "EUR GWAS -> EUR"
        if r.trait_arm == "SCZ25_META": return "multi GWAS -> pooled (zanc)" if r.target_stratum == "full" else "multi GWAS -> EUR (secondary)"
        return "ancestry-matched composite -> pooled"
    sel["cell"] = sel.apply(cell, axis=1)
    return sel[["cell", "method", "beta", "se", "p_adj", "n"]]


def compare(mains: dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = []
    for parc, root in ROOTS.items():
        a = pgc3_rows(root); a["gwas"] = "PGC3"; a["parc"] = parc
        b = new_rows(mains[parc]); b["gwas"] = "2025"; b["parc"] = parc
        frames += [a, b]
    long = pd.concat(frames, ignore_index=True)
    if long.empty:
        return long
    wide = long.pivot_table(index=["cell", "method"], columns=["gwas", "parc"],
                            values=["beta", "se", "p_adj", "n"], aggfunc="first")
    wide.columns = [f"{v}_{g}_{p}" for v, g, p in wide.columns]
    order = [c for g in ("PGC3", "2025") for p in ("dk", "hcp") for v in ("beta", "se", "p_adj", "n")
             if (c := f"{v}_{g}_{p}") in wide.columns]
    wide = wide[order].reset_index()
    cell_order = ["EUR GWAS -> EUR", "multi GWAS -> pooled (zanc)", "ancestry-matched composite -> pooled",
                  "multi GWAS -> EUR (secondary)"]
    wide["_c"] = wide.cell.map({c: i for i, c in enumerate(cell_order)})
    wide["_m"] = wide.method.map({m: i for i, m in enumerate(METHODS)})
    return wide.sort_values(["_c", "_m"]).drop(columns=["_c", "_m"])


def main() -> int:
    mains = {}
    for parc, root in ROOTS.items():
        out = collect(root)
        mains[parc] = out["main"]
        print(f"[{parc}] main {len(out['main'])} rows, all {len(out['all'])}, family {len(out['family'])}, strata {len(out['strata'])}")
    cmp_ = compare(mains)
    if not cmp_.empty:
        p = ROOTS["hcp"] / "compare" / "table_scz2025_vs_pgc3_dk_vs_hcp.tsv"
        p.parent.mkdir(parents=True, exist_ok=True)
        cmp_.to_csv(p, sep="\t", index=False, float_format="%.4g")
        print(f"\nwrote {p}\n")
        pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
        show = cmp_[["cell", "method"] + [c for c in cmp_.columns if c.startswith(("beta_", "p_adj_"))]]
        print(show.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    for parc in ROOTS:
        s = mains[parc]
        if s.empty:
            continue
        print(f"\n[{parc}] strata (global_slope, within-stratum, single-score methods + C+T best):")
        st = rd(ROOTS[parc] / "prs_scz2025" / "table_scz2025_strata.tsv")
        if st is not None:
            g = st[(st.phenotype == "global_slope") & (st.analysis.isin(["within_stratum", "heterogeneity_PRSxStratum"]))]
            g = g[(g.method != "CT") | (g.threshold == "0p05")]
            print(g[["method", "trait_arm", "threshold", "analysis", "stratum", "n", "beta", "se", "p"]]
                  .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
