"""Figure 1 DRAFT -- the paper's pitch as data, HCP-MMP throughout.

Usage (repo root):  python genetic_analysis/fig1_draft.py

Panels (all HCP-MMP, release 7.0, thickness_hcp_70_aa6e91efba82 phenotypes):
  a  design: age at scan by visit, scans per child
  b  the phenotype: per-scan cortical mean thickness vs age; one line per
     child (thin), the population trend; the child-level slope is the trait
  c  how well the slope is measured: model-based reliability of a parcel
     slope by scans per child (358 parcels), and the whole-cortex mean
  d  polygenic scores -> thinning rate (single-LMM slope): SCZ 2025, MDD,
     ALZ without APOE, education, autism; four scoring methods; EUR and
     pooled arms where the discovery GWAS allows
  e  GWAS of the thinning rate -- **MOCK** until the EUR-arm summary
     statistics are pulled (see NEEDED below); the real gwas_summary numbers
     are annotated
  f  MAGMA gene-set enrichment of the thinning-rate and baseline GWAS in
     SCZ gene sets

Reads ONLY committed tables in genetic_analysis/fig1_inputs/ (built from the
run directory and the results trees; individual-level rows are gitignored
and the script falls back to binned aggregates when they are absent).

NEEDED from CSD3 for the real panel e (any one of):
  results_70tab_hcp/assoc_eur/global_slope.fastGWA[.gz]      full sumstats
  or a thinned copy (all p < 1e-3 + 2 % of the rest, as docs/results/
  gwas_thinned/ did for the legacy run): columns CHR POS SNP P at minimum.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
IN = HERE / "fig1_inputs"
OUT = REPO / "docs/figures/fig1_draft.png"

# -- sizes: Nature Neuroscience double column 180 mm; 3-step ladder 7/6/5 pt
W_IN, H_IN = 180 / 25.4, 125 / 25.4
BASE, ANN, TICK = 7, 6, 5.5

METHOD_ORDER = ["CT", "PRSCS", "SBayesR", "SBayesRC"]
METHOD_LABEL = {"CT": "C+T", "PRSCS": "PRS-CS", "SBayesR": "SBayesR",
                "SBayesRC": "SBayesRC"}
METHOD_COLOR = {"CT": "0.30", "PRSCS": "#56B4E9", "SBayesR": "#D55E00",
                "SBayesRC": "#CC79A7"}
VISIT_COLOR = {"v0": "#cfe0f2", "v2": "#8fb8de", "v4": "#4f86c6",
               "v6": "#1f4e8c"}
VISIT_LABEL = {"v0": "baseline", "v2": "2-year", "v4": "4-year", "v6": "6-year"}
SLOPE_C, BASE_C = "#B2182B", "0.55"      # thinning rate vs baseline thickness


def _style():
    mpl.rcParams.update({
        "font.size": BASE, "axes.titlesize": BASE, "axes.labelsize": BASE,
        "xtick.labelsize": TICK, "ytick.labelsize": TICK,
        "legend.fontsize": ANN, "axes.linewidth": 0.6,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "axes.spines.top": False, "axes.spines.right": False,
    })


def letter(ax, s, dx=-0.18, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=BASE + 2,
            fontweight="bold", va="bottom", ha="left")


# ---------------------------------------------------------------- panel a --
def panel_a(ax):
    age = pd.read_csv(IN / "hcp70_age_by_visit.csv", index_col=0)
    per = pd.read_csv(IN / "hcp70_scans_per_child.csv", index_col=0)
    scans = IN / "hcp70_scans.csv"
    if scans.exists():
        sc = pd.read_csv(scans)
        bins = np.arange(8, 18.6, 0.25)
        for v in ["v0", "v2", "v4", "v6"]:
            ax.hist(sc.age[sc.visit == v], bins=bins, color=VISIT_COLOR[v],
                    lw=0)
    else:  # fallback: normal approximation from the committed summary
        x = np.linspace(8, 18.5, 300)
        for v, r in age.iterrows():
            y = r["count"] * 0.25 * np.exp(-0.5 * ((x - r["mean"]) / r["std"]) ** 2) \
                / (r["std"] * np.sqrt(2 * np.pi))
            ax.fill_between(x, y, color=VISIT_COLOR[v], lw=0)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.55)
    for v, r in age.iterrows():
        ax.text(r["mean"], ax.get_ylim()[1] * (0.62 if v in ("v0", "v4")
                                                else 0.54), VISIT_LABEL[v],
                ha="center", va="bottom", fontsize=ANN,
                color=VISIT_COLOR[v] if v != "v0" else "0.45")
    ax.set_xlabel("age at scan (years)")
    ax.set_ylabel("scans")
    ax.set_xlim(8, 18.6)
    tot = int(per.n_children.sum())
    txt = "scans per child\n" + "\n".join(
        f"{k}: {int(v):,} ({v / tot:.0%})" for k, v in per.n_children.items())
    ax.text(0.98, 0.98, txt, transform=ax.transAxes, ha="right", va="top",
            fontsize=ANN, color="0.25")
    ax.set_title(f"{tot:,} children, {int(age['count'].sum()):,} scans",
                 loc="left")


# ---------------------------------------------------------------- panel b --
def panel_b(ax, rng):
    summ = pd.read_csv(IN / "hcp70_scan_summary.csv", index_col=0).value
    scans = IN / "hcp70_scans.csv"
    if scans.exists():
        sc = pd.read_csv(scans)
        kids = rng.choice(sc.sid.unique(), 250, replace=False)
        sub = sc[sc.sid.isin(kids)].sort_values(["sid", "age"])
        for _, g in sub.groupby("sid"):
            ax.plot(g.age, g.mean_ct, color="0.6", lw=0.35, alpha=0.5,
                    zorder=1)
        # three children, highlighted, with their own OLS lines
        hi = rng.choice(sc[sc.n_visits == 4].sid.unique(), 3, replace=False)
        for k, sid in enumerate(hi):
            g = sc[sc.sid == sid].sort_values("age")
            c = ["#1f4e8c", "#D55E00", "#009E73"][k]
            ax.plot(g.age, g.mean_ct, "o-", color=c, lw=0.9, ms=2.2,
                    zorder=3)
    else:
        h = pd.read_csv(IN / "hcp70_scans_hist2d.csv")
        ax.scatter(h.age_bin + 0.125, h.ct_bin + 0.01, s=h.n / 6, color="0.6",
                   alpha=0.5, lw=0, zorder=1)
    x = np.array([8.3, 18.1])
    ax.plot(x, summ.ols_intercept_mm + summ.ols_slope_mm_per_yr * x,
            color="black", lw=1.4, zorder=4)
    ax.text(18.2, summ.ols_intercept_mm + summ.ols_slope_mm_per_yr * 18.1,
            f"population\n{summ.ols_slope_mm_per_yr * 1000:.0f} µm / year",
            fontsize=ANN, va="center", ha="left")
    ax.set_xlim(8, 20.6)
    ax.set_ylim(2.25, 3.0)
    ax.set_xticks([8, 10, 12, 14, 16, 18])
    ax.set_xlabel("age (years)")
    ax.set_ylabel("mean cortical thickness (mm)")
    ax.set_title("the trait: each child's thinning slope", loc="left")


# ---------------------------------------------------------------- panel c --
def panel_c(ax):
    rel = pd.read_csv(IN / "hcp70_regional_slope_reliability.csv")
    sh = pd.read_csv(IN / "hcp70_global_slope_splithalf.csv")
    for k in (2, 3, 4):
        v = rel.reliability[rel.n_visits == k].to_numpy()
        ax.boxplot(v, positions=[k], widths=0.5, showfliers=False,
                   medianprops=dict(color=SLOPE_C, lw=1.1),
                   boxprops=dict(color="0.35", lw=0.6),
                   whiskerprops=dict(color="0.35", lw=0.6),
                   capprops=dict(color="0.35", lw=0.6))
    ax.plot([2, 3, 4], [sh.spearman_brown[sh.n_visits == k].iloc[0]
                        for k in (2, 3, 4)], "s-", color="black", ms=3,
            lw=0.9)
    ax.text(4.15, sh.spearman_brown[sh.n_visits == 4].iloc[0],
            "cortex-wide\nmean (split-half)", fontsize=ANN, va="center")
    ax.text(4.15, rel.reliability[rel.n_visits == 4].median(),
            "one parcel\n(358 parcels)", fontsize=ANN, va="center",
            color=SLOPE_C)
    ax.set_xticks([2, 3, 4])
    ax.set_xlim(1.5, 6.4)
    ax.set_ylim(0, 1)
    ax.set_xlabel("scans per child")
    ax.set_ylabel("slope reliability")
    ax.set_title("slope reliability by scans per child", loc="left")


# ---------------------------------------------------------------- panel d --
ROWS_D = [
    ("SCZ25_EUR", "schizophrenia (2025) · EUR", "EUR"),
    ("SCZ25_META", "schizophrenia (2025) · pooled", "pooled"),
    ("MDD_eur", "depression · EUR", "EUR"),
    ("MDD_pooled", "depression · pooled", "pooled"),
    ("ALZ_noAPOE", "Alzheimer's, APOE excluded", "EUR"),
    ("EA", "education (years)", "EUR"),
    ("ASD", "autism", "EUR"),
]
DODGE = {"CT": 0.27, "PRSCS": 0.09, "SBayesR": -0.09, "SBayesRC": -0.27}


def panel_d(ax):
    t = pd.read_csv(IN / "hcp70_prs_key_arms.tsv", sep="\t")
    t = t[t.phenotype == "global_slope"]
    n = len(ROWS_D)
    for i, (arm, lab, stratum) in enumerate(ROWS_D):
        y0 = n - 1 - i
        if i % 2 == 1:
            ax.axhspan(y0 - 0.5, y0 + 0.5, color="0.95", zorder=0, lw=0)
        for m in METHOD_ORDER:
            r = t[(t.trait_arm == arm) & (t.method == m)]
            if r.empty:
                continue
            r = r.iloc[0]
            y = y0 + DODGE[m]
            c = METHOD_COLOR[m]
            sig = r.p_adj_1lmm < 0.05
            ax.errorbar(r.beta_1lmm, y, xerr=1.96 * r.se_1lmm, fmt="none",
                        ecolor=c, elinewidth=0.7, capsize=0, zorder=2)
            ax.plot(r.beta_1lmm, y, "o" if stratum == "EUR" else "D",
                    ms=2.6 if stratum == "EUR" else 2.3,
                    mfc=c if sig else "white", mec=c, mew=0.7, zorder=3)
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    ax.set_yticks(range(n))
    ax.set_yticklabels([lab for _, lab, _ in ROWS_D][::-1], fontsize=ANN)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.085, 0.145)
    ax.set_xticks([-0.05, 0, 0.05])
    ax.set_xlabel("β on thinning rate, per SD of score (95% CI)")
    ax.set_ylim(-0.55, n - 0.45)
    ax.set_title("polygenic scores → thinning rate", loc="left")
    # method key, direct
    for k, m in enumerate(METHOD_ORDER):
        ax.text(0.99, 0.97 - 0.065 * k, METHOD_LABEL[m],
                transform=ax.transAxes, ha="right", va="top", fontsize=ANN,
                color=METHOD_COLOR[m])
    ax.text(0.99, 0.30, "○ EUR arm\n   n = 4,308\n◇ pooled arm\n   n = 8,596\n"
            "   (within-ancestry z)\nfilled = p̃ < .05", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=ANN, color="0.3")
    ax.text(0.35, -0.2, "← faster thinning", transform=ax.transAxes,
            ha="center", va="top", fontsize=ANN, color="0.35")


# ---------------------------------------------------------------- panel e --
def panel_e(ax, rng):
    gw = pd.read_csv(IN / "hcp70_gwas_summary_eur.tsv", sep="\t")
    g = gw[gw.phenotype == "global_slope"].iloc[0]
    # ---- MOCK: uniform null p-values, thinned like a real Manhattan -------
    chr_len = np.array([248, 242, 198, 190, 181, 171, 159, 145, 138, 134, 135,
                        133, 114, 107, 102, 90, 83, 80, 59, 64, 47, 51])
    offs = np.concatenate([[0], np.cumsum(chr_len)[:-1]])
    xs, ys, cs = [], [], []
    for c in range(22):
        m = int(chr_len[c] * 25)
        p = rng.uniform(size=m)
        keep = (p < 1e-2) | (rng.uniform(size=m) < 0.02)
        xs.append(offs[c] + rng.uniform(0, chr_len[c], size=keep.sum()))
        ys.append(-np.log10(p[keep]))
        cs.append(np.full(keep.sum(), "0.55" if c % 2 else "0.30"))
    ax.scatter(np.concatenate(xs), np.concatenate(ys), s=1.2,
               c=np.concatenate(cs), lw=0, rasterized=True)
    ax.axhline(-np.log10(5e-8), color=SLOPE_C, lw=0.7, ls="--")
    ax.axhline(5, color="0.6", lw=0.5, ls=":")
    ax.set_xlim(0, offs[-1] + chr_len[-1])
    ax.set_xticks(offs + chr_len / 2)
    ax.set_xticklabels([str(c + 1) if (c < 12 or c % 2) else ""
                        for c in range(22)], fontsize=TICK - 1)
    ax.tick_params(axis="x", length=0)
    ax.set_xlabel("chromosome")
    ax.set_ylabel("−log10 p")
    ax.set_title("GWAS of the thinning rate", loc="left")
    ax.text(0.5, 0.42, "MOCK POINTS — awaiting sumstats",
            transform=ax.transAxes, ha="center", va="center", fontsize=ANN,
            color=SLOPE_C, alpha=0.8, fontweight="bold")
    ax.text(0.02, 0.98, f"n = {int(g.n_mean):,} · {g.n_snps / 1e6:.1f} M SNPs "
            f"· λGC {g.lambda_gc:.2f}\n{int(g.n_p5e8)} loci at p < 5×10⁻⁸ · "
            f"min p {g.min_p:.1e}", transform=ax.transAxes, ha="left",
            va="top", fontsize=ANN, color="0.3")
    ax.set_ylim(0, 11)
    ax.set_yticks([0, 2, 4, 6, 8])


# ---------------------------------------------------------------- panel f --
SETS_F = [
    ("SCZ_locus_pool", "SCZ loci, curated"),
    ("SCZ_pool_not_prio", "  ↳ non-prioritised genes"),
    ("SCZ_prioritised", "  ↳ prioritised genes"),
    ("PGC3_locus_pool", "under PGC3 peaks"),
    ("SCZ25_locus_pool", "under 2025 peaks"),
    ("SCZ25_genesig", "2025 gene-level sig."),
]


def panel_f(ax):
    m = pd.read_csv(IN / "hcp70_magma_genesets.tsv", sep="\t")
    n = len(SETS_F)
    for i, (key, lab) in enumerate(SETS_F):
        y0 = n - 1 - i
        for ph, c, dy, mk in (("global_slope", SLOPE_C, 0.16, "o"),
                              ("baseline_thickness", BASE_C, -0.16, "s")):
            r = m[(m.variable == key) & (m.phenotype == ph)].iloc[0]
            ax.errorbar(r.beta, y0 + dy, xerr=1.96 * r.se, fmt="none",
                        ecolor=c, elinewidth=0.7, capsize=0, zorder=2)
            ax.plot(r.beta, y0 + dy, mk, ms=2.6, mfc=c if r.p < 0.05
                    else "white", mec=c, mew=0.7, zorder=3)
    labs = [f"{lab} ({int(m[m.variable == k].n_genes.iloc[0])})"
            for k, lab in SETS_F]
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labs[::-1], fontsize=ANN)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.12, 0.62)
    ax.set_xticks([0, 0.2, 0.4])
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("enrichment β (95% CI)")
    ax.set_title("SCZ gene sets", loc="left")
    ax.text(0.99, 0.5, "● thinning\n   rate\n■ baseline\n   thickness\n"
            "filled =\np < .05", transform=ax.transAxes, ha="right",
            va="center", fontsize=ANN, color="0.3")


def draw() -> Path:
    _style()
    rng = np.random.default_rng(7)
    fig = plt.figure(figsize=(W_IN, H_IN))
    # [x0, y0, w, h] in figure fraction; row 1 y 0.60-0.92, row 2 y 0.10-0.45
    axa = fig.add_axes([0.060, 0.600, 0.200, 0.320])
    axb = fig.add_axes([0.345, 0.600, 0.260, 0.320])
    axc = fig.add_axes([0.700, 0.600, 0.190, 0.320])
    axd = fig.add_axes([0.190, 0.105, 0.245, 0.360])
    axe = fig.add_axes([0.495, 0.105, 0.185, 0.360])
    axf = fig.add_axes([0.858, 0.105, 0.125, 0.360])
    panel_a(axa); letter(axa, "a", dx=-0.25)
    panel_b(axb, rng); letter(axb, "b", dx=-0.18)
    panel_c(axc); letter(axc, "c", dx=-0.22)
    panel_d(axd); letter(axd, "d", dx=-0.70)
    panel_e(axe, rng); letter(axe, "e", dx=-0.20)
    panel_f(axf); letter(axf, "f", dx=-1.25)
    fig.savefig(OUT, dpi=300)
    plt.close(fig)
    print(OUT)
    return OUT


if __name__ == "__main__":
    draw()
    sys.exit(0)
