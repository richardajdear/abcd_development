"""Assemble the Experiment A results table (hpc_v3/results_v3_summary.csv).

Run as ``python hpc_v3/make_results_v3.py`` from the repo root.

One long table, one row per (panel, phenotype, disorder, stratum), drawn from
the committed result files listed in SOURCES below.  For each of the six
phenotypes the script looks the value up in the source table; a phenotype not
yet present (the four hpc_v3 phenotypes, until the cluster run lands) gets a
clearly-marked ``status=pending`` row with dummy values so the slide layout
can be designed now.  **Re-running this script after the new result files are
pulled from CSD3 flips those rows to observed automatically** — no edits here,
no hand-typed numbers (the project rule: a figure reads statistics from a
committed table, never from prose).

Choices fixed here, so the slide cannot drift from them:

* h²: Zaitlen two-GRM ``h2_snp`` on the PC-Relate GRM, relatives kept
  (n = 8,082) — one estimate per phenotype as agreed. The kinship-threshold
  sensitivity (0.183 → 0.115 for global_slope) is a caveat for text, not a
  second error bar.
* PRS: population association at C+T p < 0.5 (the project's headline
  threshold), pooled ("full") and EUR strata.  Threshold multiplicity is
  corrected with the effective-test count, p_corr = min(1, p × M_EFF) with
  M_EFF = 4.6 — the most conservative of the three estimators in
  hpc/README_HPC.md (range 2.7–4.6) — not Bonferroni × 8.
* MAGMA: the full locus pools only (SCZ_locus_pool, MDD_pool), marginal
  model, EUR.  Unsigned test — say so wherever these are shown.
* rg: LDSC vs SCZ/MDD, EUR, with the source table's ``underpowered`` flag
  (either trait h² z < 4) carried through.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

PHENOS = ["baseline_thickness", "global_slope", "slope_topDelta",
          "slope_topC3", "slope_projDelta", "slope_projC3"]

M_EFF = 4.6          # conservative end of hpc/README_HPC.md's 2.7-4.6
PRS_THRESHOLD = "0p5"

SOURCES = {
    "h2": "hpc_v2/work/results_v2/reml_zaitlen_pcrel/reml_zaitlen_summary.tsv",
    "prs": "hpc/work/results/prs/prs_association.tsv",
    "magma_scz": "hpc_v2/work/results_v2/magma_prio_eur/sczprio_gsa_summary.tsv",
    "magma_mdd": "hpc_v2/work/results_v2/magma_prio_eur/mddhc_gsa_summary.tsv",
    "rg": "hpc_v2/work/results_v2/ldsc_eur/ldsc_rg_summary.tsv",
}

#: dummy values for pending rows -- chosen to be visually plausible so the
#: layout is realistic, and flagged status=pending so the slide greys them out.
DUMMY = {
    "h2": dict(estimate=0.15, se=0.10),
    "prs": dict(estimate=0.0, se=0.016),
    "magma": dict(estimate=0.0, se=0.05),
    "rg": dict(estimate=0.0, se=0.11),
}


def rows() -> list[dict]:
    out: list[dict] = []

    def add(panel, pheno, disorder, stratum, hit, source, **kw):
        base = dict(panel=panel, phenotype=pheno, disorder=disorder,
                    stratum=stratum, source=source)
        if hit is None:
            out.append(base | DUMMY[panel if panel in DUMMY else "prs"]
                       | dict(p=np.nan, n=np.nan, status="pending") | kw)
        else:
            out.append(base | hit | dict(status="observed") | kw)

    # ---- h2 -----------------------------------------------------------------
    h2 = pd.read_csv(REPO / SOURCES["h2"], sep="\t")
    for ph in PHENOS:
        m = h2[h2.phenotype == ph]
        hit = (dict(estimate=float(m.h2_snp.iloc[0]),
                    se=float(m.se_snp.iloc[0]), p=float(m.pval.iloc[0]),
                    n=int(m.n.iloc[0])) if len(m) else None)
        add("h2", ph, "", "", hit, SOURCES["h2"])

    # ---- PRS ----------------------------------------------------------------
    prs = pd.read_csv(REPO / SOURCES["prs"], sep="\t")
    prs = prs[prs.threshold == PRS_THRESHOLD]
    for ph in PHENOS:
        for dis in ("SCZ", "MDD"):
            for stratum in ("EUR", "full"):
                m = prs[(prs.phenotype == ph) & (prs.disorder == dis)
                        & (prs.stratum == stratum)]
                hit = (dict(estimate=float(m.beta.iloc[0]),
                            se=float(m.se.iloc[0]), p=float(m.p.iloc[0]),
                            n=int(m.n.iloc[0])) if len(m) else None)
                add("prs", ph, dis, stratum, hit, SOURCES["prs"],
                    m_eff=M_EFF, threshold=PRS_THRESHOLD)

    # ---- MAGMA locus pools --------------------------------------------------
    for dis, key, gset in (("SCZ", "magma_scz", "SCZ_locus_pool"),
                           ("MDD", "magma_mdd", "MDD_pool")):
        g = pd.read_csv(REPO / SOURCES[key], sep="\t")
        g = g[(g.model == "marginal") & (g["set"] == gset)].drop_duplicates(
            subset=["trait", "set"])
        for ph in PHENOS:
            m = g[g.trait == ph]
            hit = (dict(estimate=float(m.beta.iloc[0]),
                        se=float(m.se.iloc[0]), p=float(m.p.iloc[0]),
                        n=int(m.n_genes.iloc[0])) if len(m) else None)
            add("magma", ph, dis, "EUR", hit, SOURCES[key], gene_set=gset)

    # ---- LDSC rg ------------------------------------------------------------
    rg = pd.read_csv(REPO / SOURCES["rg"], sep="\t")
    for ph in PHENOS:
        for dis in ("SCZ", "MDD"):
            m = rg[(rg.phenotype == ph) & (rg.disorder == dis)]
            if len(m) and pd.notna(m.rg.iloc[0]):
                hit = dict(estimate=float(m.rg.iloc[0]),
                           se=float(m.se.iloc[0]), p=float(m.p.iloc[0]),
                           n=np.nan)
                add("rg", ph, dis, "EUR", hit, SOURCES["rg"],
                    h2_z=float(m.h2_z.iloc[0]),
                    underpowered=str(m.underpowered.iloc[0]))
            else:
                add("rg", ph, dis, "EUR", None, SOURCES["rg"],
                    h2_z=np.nan, underpowered="pending")
    return out


def main() -> int:
    df = pd.DataFrame(rows())
    n_obs = int((df.status == "observed").sum())
    out = HERE / "results_v3_summary.csv"
    df.to_csv(out, index=False)
    print(f"{out}: {len(df)} rows, {n_obs} observed / "
          f"{len(df) - n_obs} pending")
    print(df[df.status == "observed"]
          .groupby("panel").size().rename("observed_rows").to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
