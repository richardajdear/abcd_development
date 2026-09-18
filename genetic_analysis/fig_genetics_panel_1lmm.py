#!/usr/bin/env python
"""Figure: the genetics of the single-LMM phenotypes (7.0 tabulation) in one view.

Rows: global slope (thinning rate) and baseline thickness, both the single-LMM
construction (one mixed model on the per-scan whole-cortex mean).  Columns:
  (a) SNP heritability -- GCTA GREML on the dense imputed GRM, PC-AiR unrelated
      set (n 6,011), and LDSC on the EUR-arm GWAS (n ~4,300);
  (b) polygenic score association, beta per SD of phenotype per SD of score with
      95 % CI, SBayesRC, EUR-arm matched cells (EUR discovery GWAS -> 4,308 EUR
      children) so the three methods share one ancestry frame;
  (c) LDSC genetic correlation with 95 % CI, EUR discovery GWAS vs the EUR-arm
      phenotype GWAS; drawn hollow-and-faint where the phenotype's own h2 z < 4
      (rule 13: not interpretable);
  (d) MAGMA gene-property, disorder gene-level Z ~ phenotype gene-level Z
      (beta with 95 % CI), the LD-matched disorder gene result for each trait.
DK = circle, HCP-MMP = square; filled = p < 0.05 (PRS: threshold-adjusted).
One hue throughout (emphasis form: significance is the story, not identity).

Reads committed tables only; run from anywhere:
    python genetic_analysis/fig_genetics_panel_1lmm.py
Output: docs/figures/fig_genetics_panel_1lmm.png (+ .tsv twin of every value)
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[1]
ROOTS = {"DK": REPO / "genetic_analysis/work/results_70tab", "HCP": REPO / "genetic_analysis/work/results_70tab_hcp"}
OUT = REPO / "docs/figures/fig_genetics_panel_1lmm"
PHENOS = [("global_slope", "Global slope (thinning rate), single LMM"),
          ("baseline_thickness", "Baseline thickness, single LMM")]
# trait rows: label, PRS arm, LDSC disorder file, MAGMA disorder result
TRAITS = [("Schizophrenia (2025 GWAS)", "SCZ25_EUR", "SCZ25_EUR", "SCZ25_META"),
          ("Depression (MDD 2025)",     "MDD_eur",   "MDD_EUR",   "MDD_EUR"),
          ("Autism",                    "ASD",       "ASD",       "ASD"),
          ("Alzheimer's (Wightman)",    "ALZ",       None,        "ALZ"),
          ("Alzheimer's, APOE excluded","ALZ_noAPOE","ALZ_noAPOE","ALZ_noAPOE"),
          ("Educational attainment",    "EA",        "EA",        "EA")]
PRS_METHOD = "SBayesRC"
MAGMA_KIND = "reverse"          # disorder genes ~ phenotype gene Z
BLUE, INK, MUTED, GRID = "#2a78d6", "#0b0b0b", "#52514e", "#e6e5e1"
MARK = {"DK": "o", "HCP": "s"}

def rd(p):
    return pd.read_csv(p, sep="\t") if Path(p).exists() else None

def h2_rows():
    rows = []
    for atlas, root in ROOTS.items():
        g = rd(root / "scan_1lmm/reml_imp_pooled/reml_summary.tsv")
        if g is not None:
            for _, r in g.iterrows():
                rows.append(dict(phenotype=r.phenotype.replace("_1lmm", ""), atlas=atlas, estimator="GREML", h2=r.h2, se=r.se, n=r.n))
        l = rd(root / "ldsc_1lmm/table_ldsc_panel.tsv")
        if l is not None:
            l = l[(l.atlas == atlas.lower()) & (l.construction == "1lmm")].drop_duplicates(["phenotype"])
            for _, r in l.iterrows():
                rows.append(dict(phenotype=r.phenotype, atlas=atlas, estimator="LDSC", h2=r.h2_obs, se=r.h2_obs_se, n=np.nan))
    return pd.DataFrame(rows)

def prs_rows():
    rows = []
    for atlas, root in ROOTS.items():
        a = rd(root / "prs_final_1lmm/table_order_of_operations_all.tsv")
        b = rd(root / "prs_scz2025/table_order_of_operations_scz2025.tsv")
        for label, arm, _, _ in TRAITS:
            t = b if arm.startswith("SCZ25") else a
            if t is None: continue
            s = t[(t.trait_arm == arm) & (t.method == PRS_METHOD) & (t.stratum == "EUR")]
            for ph, _ in PHENOS:
                r = s[s.phenotype == ph]
                if r.empty: continue
                r = r.iloc[0]
                rows.append(dict(trait=label, phenotype=ph, atlas=atlas, est=r.beta_1lmm, se=r.se_1lmm, p=r.p_adj_1lmm, n=r.n))
    return pd.DataFrame(rows)

def rg_rows():
    rows = []
    t = rd(ROOTS["DK"] / "ldsc_1lmm/table_ldsc_panel.tsv")
    if t is None: return pd.DataFrame(rows)
    t = t[t.construction == "1lmm"]
    for label, _, dis, _ in TRAITS:
        if dis is None: continue
        for ph, _ in PHENOS:
            for atlas in ROOTS:
                r = t[(t.disorder == dis) & (t.phenotype == ph) & (t.atlas == atlas.lower())]
                if r.empty: continue
                r = r.iloc[0]
                rows.append(dict(trait=label, phenotype=ph, atlas=atlas, est=r.rg, se=r.se, p=r.p, h2_z=r.h2_z, underpowered=r.underpowered == "yes"))
    return pd.DataFrame(rows)

def magma_rows():
    rows = []
    for atlas, root in ROOTS.items():
        t = rd(root / "magma_panel/table_magma_panel.tsv")
        if t is None: continue
        t = t[(t.kind == MAGMA_KIND) & (t.construction == "1lmm")]
        for label, _, _, dis in TRAITS:
            for ph, _ in PHENOS:
                r = t[(t.disorder_result == dis) & (t.phenotype == ph)]
                if r.empty: continue
                r = r.iloc[0]
                rows.append(dict(trait=label, phenotype=ph, atlas=atlas, est=r.beta, se=r.se, p=r.p, n_genes=r.n_genes))
    return pd.DataFrame(rows)

def forest(ax, t, ph, xlabel, title, zero=True, faint_col=None):
    labels = [x[0] for x in TRAITS]; ny = len(labels)
    ax.set_yticks(range(ny)); ax.set_yticklabels(labels[::-1] if False else labels, fontsize=8.5)
    ax.invert_yaxis(); ax.set_ylim(ny - 0.5, -0.5)
    for i in range(ny):
        if i % 2: ax.axhspan(i - 0.5, i + 0.5, color="#f6f5f2", zorder=0, lw=0)
    if zero: ax.axvline(0, color=MUTED, lw=1, zorder=1)
    for i, lab in enumerate(labels):
        for atlas, dy in (("DK", -0.18), ("HCP", 0.18)):
            r = t[(t.trait == lab) & (t.phenotype == ph) & (t.atlas == atlas)] if not t.empty else t
            if r.empty:
                ax.text(0, i + dy, "n/a", fontsize=6.5, color=MUTED, ha="center", va="center"); continue
            r = r.iloc[0]
            faint = bool(r[faint_col]) if faint_col and faint_col in r else False
            sig = (r.p < 0.05) and not faint
            col = BLUE if not faint else "#9ec5f4"
            ax.errorbar(r.est, i + dy, xerr=1.96 * r.se, fmt="none", ecolor=col, elinewidth=1.4, capsize=0, zorder=2)
            ax.plot(r.est, i + dy, marker=MARK[atlas], ms=7, mfc=col if sig else "white", mec=col, mew=1.4, zorder=3)
    ax.set_xlabel(xlabel, fontsize=8); ax.set_title(title, fontsize=9.5, loc="left", pad=6, color=INK)
    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID); ax.tick_params(axis="x", labelsize=8, colors=MUTED); ax.tick_params(axis="y", length=0)

def h2_panel(ax, h2, ph, title):
    ests = ["GREML", "LDSC"]; ax.set_yticks(range(2)); ax.set_yticklabels(["GREML (n 6,011)", "LDSC (EUR arm)"], fontsize=8.5)
    ax.set_ylim(1.5, -0.5); ax.axvline(0, color=MUTED, lw=1)
    for i, e in enumerate(ests):
        for atlas, dy in (("DK", -0.18), ("HCP", 0.18)):
            r = h2[(h2.phenotype == ph) & (h2.estimator == e) & (h2.atlas == atlas)] if not h2.empty else h2
            if r.empty: ax.text(0.02, i + dy, "n/a", fontsize=6.5, color=MUTED, va="center"); continue
            r = r.iloc[0]; z = r.h2 / r.se if r.se else 0
            ax.errorbar(r.h2, i + dy, xerr=1.96 * r.se, fmt="none", ecolor=BLUE, elinewidth=1.4, capsize=0)
            ax.plot(r.h2, i + dy, marker=MARK[atlas], ms=7, mfc=BLUE if z >= 1.96 else "white", mec=BLUE, mew=1.4)
            ax.annotate(f"{r.h2:.2f}", (r.h2, i + dy), xytext=(0, 7 if dy < 0 else -9), textcoords="offset points", fontsize=6.5, ha="center", color=MUTED)
    ax.set_xlim(-0.15, 0.75); ax.set_xlabel("SNP heritability h² (95% CI)", fontsize=8); ax.set_title(title, fontsize=9.5, loc="left", pad=6)
    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID); ax.tick_params(axis="x", labelsize=8, colors=MUTED); ax.tick_params(axis="y", length=0)

def main():
    h2, prs, rg, mg = h2_rows(), prs_rows(), rg_rows(), magma_rows()
    for name, t in (("h2", h2), ("prs", prs), ("rg", rg), ("magma", mg)):
        t.to_csv(f"{OUT}_{name}.tsv", sep="\t", index=False)
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "figure.facecolor": "white", "savefig.facecolor": "white"})
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.4), gridspec_kw=dict(width_ratios=[1.0, 1.25, 1.25, 1.25], hspace=0.55, wspace=0.55))
    for row, (ph, phlab) in enumerate(PHENOS):
        h2_panel(axes[row, 0], h2, ph, f"{phlab}\n(a) heritability")
        forest(axes[row, 1], prs, ph, "β (SD phenotype per SD score), 95% CI", f"(b) PRS β, {PRS_METHOD}, EUR arm")
        forest(axes[row, 2], rg, ph, "genetic correlation rg, 95% CI", "(c) LDSC rg (faint: h² z < 4)", faint_col="underpowered")
        forest(axes[row, 3], mg, ph, "β (disorder gene Z ~ phenotype gene Z), 95% CI", "(d) MAGMA gene-property")
        for c in (2, 3): axes[row, c].set_yticklabels([])
        axes[row, 1].set_yticklabels([x[0] for x in TRAITS], fontsize=8.5)
    handles = [Line2D([], [], marker="o", color=BLUE, ls="", ms=7, label="DK (68 regions)"),
               Line2D([], [], marker="s", color=BLUE, ls="", ms=7, label="HCP-MMP (358 regions)"),
               Line2D([], [], marker="o", color=BLUE, mfc=BLUE, ls="", ms=7, label="p < 0.05"),
               Line2D([], [], marker="o", color=BLUE, mfc="white", ls="", ms=7, label="p ≥ 0.05")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("ABCD 7.0 cortical development × psychiatric polygenic risk — single-LMM phenotypes, 8,596 children", fontsize=11.5, x=0.02, ha="left", y=0.995)
    fig.text(0.02, 0.955, "Rows: the two phenotypes. Columns: heritability; polygenic-score association; genetic correlation; MAGMA gene-level association. "
             "European-ancestry discovery GWAS throughout; PRS shown for the European arm so the three methods share one ancestry frame "
             "(pooled-arm PRS results are in the run log).", fontsize=8, color=MUTED, ha="left")
    fig.savefig(f"{OUT}.png", dpi=200, bbox_inches="tight")
    print(f"wrote {OUT}.png; rows h2 {len(h2)}, prs {len(prs)}, rg {len(rg)}, magma {len(mg)}")

if __name__ == "__main__":
    sys.exit(main())
