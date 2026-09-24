"""Figure 1 DRAFT -- the paper's pitch as data, HCP-MMP throughout.

Usage (repo root):  python genetic_analysis/fig1_draft.py

Panels (HCP-MMP1.0, release 7.0, run thickness_hcp_70_aa6e91efba82; the
trait is the single-LMM slope of the per-scan cortical mean):
  a  group maps: baseline thickness (visit 0) and the mean thinning rate per
     parcel (bilateral average, drawn on the left hemisphere)
  b  design: age at scan by visit, scans per child
  c  the trait: per-scan cortical mean vs age, one line per child
  d  reliability of a parcel slope vs the cortex-wide slope, by scans
  e  polygenic scores -> thinning rate      (C+T, PRS-CS, SBayesRC)
  f  polygenic scores -> baseline thickness (same scores, same children)
  g  MAGMA: SCZ and MDD locus gene sets, thinning rate and baseline thickness,
     EUR arm (1000G EUR LD reference); pooled arm (ABCD-sample LD reference,
     step 14) drawn as n.d. until its table is committed

Reads ONLY committed tables in genetic_analysis/fig1_inputs/; individual-level
rows are gitignored and the script falls back to binned aggregates.
Writes the figure and a caption whose numbers are read from the same tables.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
IN = HERE / "fig1_inputs"
OUT = REPO / "docs/figures/fig1_draft.png"
CAPTION = REPO / "docs/figures/fig1_draft_caption.md"
POLY = REPO / "ahba_pls/data/hcp_polygons.csv"
POOLED_MAGMA = HERE / "work/results_70tab_hcp/magma_pooled/table_magma_pooled.tsv"

# -- Nature Neuroscience double column 180 mm; 3-step ladder 7/6/5.5 pt
W_IN, H_IN = 180 / 25.4, 118 / 25.4
BASE, ANN, TICK = 7, 6, 5.5

#: SBayesR dropped: redundant with SBayesRC (same prior family, r > .9)
METHOD_ORDER = ["CT", "PRSCS", "SBayesRC"]
METHOD_LABEL = {"CT": "C+T", "PRSCS": "PRS-CS", "SBayesRC": "SBayesRC"}
METHOD_COLOR = {"CT": "0.30", "PRSCS": "#56B4E9", "SBayesRC": "#CC79A7"}
VISIT_COLOR = {"v0": "#cfe0f2", "v2": "#8fb8de", "v4": "#4f86c6",
               "v6": "#1f4e8c"}
VISIT_LABEL = {"v0": "baseline", "v2": "2-year", "v4": "4-year", "v6": "6-year"}
#: colour threading: thinning = red family, baseline thickness = grey family
SLOPE_C, BASE_C = "#B2182B", "0.45"


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


# -------------------------------------------------------------- panel a --
def _polys() -> pd.DataFrame:
    p = pd.read_csv(POLY)
    out = []
    for view, off in (("lateral", 0.0), ("medial", 1.0)):
        v = p[p.view == view].copy()
        v["x"] = (v.x - v.x.min()) / (v.x.max() - v.x.min())
        v["y"] = (v.y - v.y.min()) / (v.y.max() - v.y.min()) * 0.62 - off * 0.72
        out.append(v)
    return pd.concat(out)


def panel_map(ax, fig, values: pd.Series, cmap, vmin, vmax, title, cbar_label,
              ticks, ticklabels):
    polys = _polys()
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    verts, cols = [], []
    for (label, view, g, sg), q in polys.groupby(["label", "view", "group",
                                                  "subgroup"], sort=False):
        verts.append(q[["x", "y"]].to_numpy())
        v = values.get(label, np.nan)
        cols.append(cmap(norm(v)) if np.isfinite(v) else (0.93, 0.93, 0.92, 1))
    ax.add_collection(PolyCollection(verts, facecolors=cols,
                                     edgecolors="white", linewidths=0.15))
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-1.0, 0.66)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=BASE, loc="left", pad=1, x=0.02)
    ax.text(1.0, 0.63, "lateral", fontsize=TICK, color="0.5", ha="right",
            va="top")
    ax.text(1.0, -0.09, "medial", fontsize=TICK, color="0.5", ha="right",
            va="top")
    cax = ax.inset_axes([0.10, 0.085, 0.80, 0.030])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                      orientation="horizontal")
    cb.outline.set_visible(False)
    cb.set_ticks(ticks); cb.set_ticklabels(ticklabels)
    cb.ax.tick_params(labelsize=TICK, length=2, width=0.5, pad=1)
    cb.set_label(cbar_label, fontsize=ANN, labelpad=1)


def panel_maps(ax_top, ax_bot, fig):
    m = pd.read_csv(IN / "hcp70_group_maps.csv").set_index("label")
    lo, hi = np.percentile(m.baseline_ct, [2, 98])
    lo, hi = np.floor(lo * 10) / 10, np.ceil(hi * 10) / 10
    panel_map(ax_top, fig, m.baseline_ct, mpl.colormaps["Greys"], lo, hi,
              "baseline thickness", "mm (age ~10)", [lo, (lo + hi) / 2, hi],
              [f"{lo:.1f}", f"{(lo + hi) / 2:.2f}", f"{hi:.1f}"])
    um = m.slope_mm_per_yr * 1000
    vmin = np.floor(np.percentile(um, 2) / 5) * 5
    panel_map(ax_bot, fig, um, mpl.colormaps["Reds_r"], vmin, 0,
              "change in thickness", "µm / year", [vmin, vmin / 2, 0],
              [f"{vmin:.0f}", f"{vmin / 2:.0f}", "0"])


def panel_design(ax):
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
    ax.set_xticks([8, 10, 12, 14, 16, 18])
    tot = int(per.n_children.sum())
    txt = f"{int(age['count'].sum()):,} scans\nscans per child\n" + "\n".join(
        f"{k}: {int(v):,} ({v / tot:.0%})" for k, v in per.n_children.items())
    ax.text(0.98, 0.98, txt, transform=ax.transAxes, ha="right", va="top",
            fontsize=ANN, color="0.25")
    ax.set_title(f"{tot:,} children", loc="left")


def panel_traj(ax, rng):
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
            f"{summ.ols_slope_mm_per_yr * 1000:.0f} µm\n/ year",
            fontsize=ANN, va="center", ha="left")
    handles = [
        mpl.lines.Line2D([], [], color="0.6", lw=0.6, alpha=0.8,
                         label="one child (250 of 8,716 shown)"),
        mpl.lines.Line2D([], [], color="#1f4e8c", lw=0.9, marker="o", ms=2.2,
                         label="three children with 4 scans"),
        mpl.lines.Line2D([], [], color="black", lw=1.4,
                         label="population trend (OLS)"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False,
              fontsize=ANN, handlelength=1.6, handletextpad=0.5,
              borderaxespad=0.2, labelspacing=0.3)
    ax.set_xlim(8, 19.9)
    ax.set_ylim(2.25, 3.0)
    ax.set_xticks([8, 10, 12, 14, 16, 18])
    ax.set_xlabel("age (years)")
    ax.set_ylabel("mean cortical thickness (mm)")
    ax.set_title("each child's slope is the trait", loc="left")


def panel_rel(ax):
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
    ax.text(4.45, sh.spearman_brown[sh.n_visits == 4].iloc[0],
            "cortex-wide\nmean (split-half)", fontsize=ANN, va="center")
    ax.text(4.45, rel.reliability[rel.n_visits == 4].median(),
            "one parcel\n(358 parcels)", fontsize=ANN, va="center",
            color=SLOPE_C)
    ax.set_xticks([2, 3, 4])
    ax.set_xlim(1.5, 6.4)
    ax.set_ylim(0, 1)
    ax.set_xlabel("scans per child")
    ax.set_ylabel("slope reliability")
    ax.set_title("slope reliability", loc="left")


# ---------------------------------------------------------- panels e, f --
ROWS_PRS = [
    ("SCZ25_EUR", "schizophrenia · EUR", "EUR"),
    ("SCZ25_META", "schizophrenia · pooled", "pooled"),
    ("MDD_eur", "depression · EUR", "EUR"),
    ("MDD_pooled", "depression · pooled", "pooled"),
    ("ALZ", "Alzheimer's", "EUR"),
    ("ALZ_noAPOE", "  ↳ APOE excluded", "EUR"),
    ("EA", "education", "EUR"),
    ("ASD", "autism", "EUR"),
]
DODGE = {"CT": 0.24, "PRSCS": 0.0, "SBayesRC": -0.24}


def panel_prs(ax, phenotype, title, show_ylab, colour):
    t = pd.read_csv(IN / "hcp70_prs_key_arms.tsv", sep="\t")
    t = t[(t.phenotype == phenotype) & t.method.isin(METHOD_ORDER)]
    n = len(ROWS_PRS)
    for i, (arm, lab, stratum) in enumerate(ROWS_PRS):
        y0 = n - 1 - i
        if i % 2 == 1:
            ax.axhspan(y0 - 0.5, y0 + 0.5, color="0.95", zorder=0, lw=0)
        for m in METHOD_ORDER:
            r = t[(t.trait_arm == arm) & (t.method == m)]
            if r.empty:
                continue
            r = r.iloc[0]
            y, c = y0 + DODGE[m], METHOD_COLOR[m]
            sig = r.p_adj_1lmm < 0.05
            ax.errorbar(r.beta_1lmm, y, xerr=1.96 * r.se_1lmm, fmt="none",
                        ecolor=c, elinewidth=0.7, capsize=0, zorder=2)
            ax.plot(r.beta_1lmm, y, "o" if stratum == "EUR" else "D",
                    ms=2.7 if stratum == "EUR" else 2.4,
                    mfc=c if sig else "white", mec=c, mew=0.7, zorder=3)
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    ax.set_yticks(range(n))
    ax.set_yticklabels([lab for _, lab, _ in ROWS_PRS][::-1] if show_ylab
                       else [], fontsize=ANN)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.09, 0.09)
    ax.set_xticks([-0.05, 0, 0.05])
    ax.set_ylim(-0.55, n - 0.45)
    ax.set_title(title, loc="left", color=colour)


# -------------------------------------------------------------- panel g --
SETS_G = [("SCZ_locus_pool", "SCZ loci"), ("MDD_pool", "MDD loci")]
SLOTS_G = [("global_slope", "EUR", SLOPE_C, 0.30),
           ("global_slope", "pooled", SLOPE_C, 0.10),
           ("baseline_thickness", "EUR", BASE_C, -0.10),
           ("baseline_thickness", "pooled", BASE_C, -0.30)]


def _magma():
    eur = pd.read_csv(IN / "hcp70_magma_locus_sets.tsv", sep="\t")
    if POOLED_MAGMA.exists():                 # step 14, once committed
        pl = pd.read_csv(POOLED_MAGMA, sep="\t")
        pl = pl[(pl.kind == "gene-set") & pl.variable.isin(dict(SETS_G))]
        if "construction" in pl:
            pl = pl[pl.construction == "1lmm"]
        pl = pl.assign(arm="pooled")
        eur = pd.concat([eur, pl[eur.columns.intersection(pl.columns)]])
    return eur


def panel_magma(ax):
    m = _magma()
    n = len(SETS_G)
    for i, (key, lab) in enumerate(SETS_G):
        y0 = n - 1 - i
        if i % 2 == 1:
            ax.axhspan(y0 - 0.5, y0 + 0.5, color="0.95", zorder=0, lw=0)
        for ph, arm, c, dy in SLOTS_G:
            r = m[(m.variable == key) & (m.phenotype == ph) & (m.arm == arm)]
            mk = "o" if arm == "EUR" else "D"
            if r.empty:
                ax.text(0.0, y0 + dy, "n.d.", fontsize=TICK - 0.5,
                        color="0.6", ha="center", va="center")
                continue
            r = r.iloc[0]
            ax.errorbar(r.beta, y0 + dy, xerr=1.96 * r.se, fmt="none",
                        ecolor=c, elinewidth=0.7, capsize=0, zorder=2)
            ax.plot(r.beta, y0 + dy, mk, ms=2.7, mfc=c if r.p < 0.05
                    else "white", mec=c, mew=0.7, zorder=3)
    ng = {k: int(m[m.variable == k].n_genes.iloc[0]) for k, _ in SETS_G}
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{lab}\n({ng[k]:,} genes)" for k, lab in SETS_G][::-1],
                       fontsize=ANN)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_xlim(-0.12, 0.36)
    ax.set_xticks([0, 0.15, 0.3])
    ax.set_ylim(-0.55, n - 0.45)
    ax.set_xlabel("enrichment β (95% CI)")
    ax.set_title("MAGMA gene sets", loc="left")
    ax.text(0.98, 0.10, "thinning rate", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=ANN, color=SLOPE_C)
    ax.text(0.98, 0.03, "baseline", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=ANN, color=BASE_C)


# ----------------------------------------------------------------- draw --
def draw() -> Path:
    _style()
    rng = np.random.default_rng(7)
    fig = plt.figure(figsize=(W_IN, H_IN))
    # left column: the two maps span both rows
    ama = fig.add_axes([0.006, 0.505, 0.160, 0.405])
    amb = fig.add_axes([0.006, 0.040, 0.160, 0.405])
    # top row
    axb = fig.add_axes([0.250, 0.615, 0.170, 0.300])
    axc = fig.add_axes([0.495, 0.615, 0.205, 0.300])
    axd = fig.add_axes([0.805, 0.615, 0.105, 0.300])
    # bottom row
    axe = fig.add_axes([0.330, 0.125, 0.165, 0.345])
    axf = fig.add_axes([0.525, 0.125, 0.165, 0.345])
    axg = fig.add_axes([0.830, 0.125, 0.150, 0.345])

    panel_maps(ama, amb, fig)
    fig.text(0.006, 0.955, "a", fontsize=BASE + 2, fontweight="bold", va="bottom")
    panel_design(axb); letter(axb, "b", dx=-0.30)
    panel_traj(axc, rng); letter(axc, "c", dx=-0.22)
    panel_rel(axd); letter(axd, "d", dx=-0.45)
    panel_prs(axe, "global_slope", "PRS → thinning rate", True, SLOPE_C)
    letter(axe, "e", dx=-0.95)
    panel_prs(axf, "baseline_thickness", "PRS → baseline thickness", False,
              BASE_C)
    letter(axf, "f", dx=-0.10)
    fig.text((0.330 + 0.690) / 2, 0.062, "β per SD of polygenic score (95% CI)",
             ha="center", va="bottom", fontsize=BASE)
    panel_magma(axg); letter(axg, "g", dx=-0.62)

    # shared key for e-g, under the bottom row
    keys = [mpl.lines.Line2D([], [], color=METHOD_COLOR[m], lw=0.8,
                             marker="o", ms=2.7, label=METHOD_LABEL[m])
            for m in METHOD_ORDER]
    keys += [mpl.lines.Line2D([], [], color="0.3", ls="", marker="o", ms=2.7,
                              mfc="white", label="EUR arm (n = 4,308)"),
             mpl.lines.Line2D([], [], color="0.3", ls="", marker="D", ms=2.4,
                              mfc="white",
                              label="pooled, within-ancestry z (n = 8,596)"),
             mpl.lines.Line2D([], [], color="0.3", ls="", marker="o", ms=2.7,
                              label="filled = p < .05")]
    fig.legend(handles=keys, loc="lower left", bbox_to_anchor=(0.215, -0.008),
               ncol=6, frameon=False, fontsize=ANN, handlelength=1.3,
               handletextpad=0.4, columnspacing=1.1)
    fig.savefig(OUT, dpi=300)
    plt.close(fig)
    print(OUT)
    return OUT


def write_caption() -> Path:
    """Caption with every number read from fig1_inputs/ at write time."""
    per = pd.read_csv(IN / "hcp70_scans_per_child.csv", index_col=0)
    age = pd.read_csv(IN / "hcp70_age_by_visit.csv", index_col=0)
    summ = pd.read_csv(IN / "hcp70_scan_summary.csv", index_col=0).value
    maps = pd.read_csv(IN / "hcp70_group_maps.csv")
    rel = pd.read_csv(IN / "hcp70_regional_slope_reliability.csv")
    sh = pd.read_csv(IN / "hcp70_global_slope_splithalf.csv")
    prs = pd.read_csv(IN / "hcp70_prs_key_arms.tsv", sep="\t")
    prs = prs[prs.method.isin(METHOD_ORDER)]
    mg = _magma()

    n_kids, n_scans = int(per.n_children.sum()), int(age["count"].sum())
    med = {k: rel.reliability[rel.n_visits == k].median() for k in (2, 3, 4)}
    sb = {int(r.n_visits): r.spearman_brown for _, r in sh.iterrows()}
    um = maps.slope_mm_per_yr * 1000

    def k(arm, ph):
        d = prs[(prs.trait_arm == arm) & (prs.phenotype == ph)]
        return f"{int((d.p_adj_1lmm < 0.05).sum())}/{len(d)}"

    def br(arm, ph="global_slope"):
        d = prs[(prs.trait_arm == arm) & (prs.phenotype == ph)]
        return f"{d.beta_1lmm.min():.3f} to {d.beta_1lmm.max():.3f}"

    def mp(key, ph):
        r = mg[(mg.variable == key) & (mg.phenotype == ph) & (mg.arm == "EUR")]
        return f"{r.p.iloc[0]:.1g}"

    nsz = {kk: int(mg[mg.variable == kk].n_genes.iloc[0]) for kk, _ in SETS_G}
    pooled_state = ("pooled-arm results use the ABCD analysis sample itself as "
                    "the LD reference" if (mg.arm == "pooled").any() else
                    "pooled-arm results (ABCD analysis sample as LD reference) "
                    "are pending and marked n.d.")

    txt = f"""**Figure 1 | Polygenic risk for schizophrenia predicts the rate, not the
baseline level, of adolescent cortical thinning.** ABCD release 7.0; all panels
use the HCP-MMP1.0 parcellation (358 parcels).

**a**, Group maps: mean thickness at the baseline visit (top) and mean
per-child thinning rate (bottom) per parcel, averaged over hemispheres and
drawn on the left hemisphere; every parcel thins (range {um.min():.0f} to
{um.max():.1f} µm per year). **b**, Age at scan by visit for the {n_kids:,}
children with at least two usable scans ({n_scans:,} scans;
{int(per.n_children[2]):,} / {int(per.n_children[3]):,} / {int(per.n_children[4]):,} children with 2 / 3 / 4 scans). **c**, Mean cortical
thickness per scan against age: 250 randomly drawn children (grey), three
children with four scans (colour), and the population trend (black; {summ.ols_slope_mm_per_yr * 1000:.0f} µm per
year). The trait analysed in **e–g** is each child's age slope from one linear
mixed model on the per-scan cortical mean. **d**, Reliability of the slope:
model-based reliability of a single parcel's slope across 358 parcels (boxes;
median {med[2]:.2f} / {med[3]:.2f} / {med[4]:.2f} for 2 / 3 / 4 scans) and split-half consistency of the
cortex-wide slope (squares; left vs right hemisphere, Spearman–Brown corrected,
{sb[2]:.2f} / {sb[3]:.2f} / {sb[4]:.2f}). **e, f**, Association of polygenic scores with the thinning rate
(**e**) and with baseline thickness (**f**) in the same children, β per SD of
score with 95 % CI, three scoring methods. Circles, European-ancestry arm
scored with a European discovery GWAS; diamonds, all ancestries scored with a
multi-ancestry discovery GWAS and standardised within genetic-ancestry cluster.
Filled, p < 0.05 (C+T corrected over its thresholds). Schizophrenia (2025
GWAS) is significant for thinning under {k("SCZ25_EUR", "global_slope")} methods in the European arm (β
{br("SCZ25_EUR")}) and {k("SCZ25_META", "global_slope")} pooled, and for baseline thickness under
{k("SCZ25_EUR", "baseline_thickness")} and {k("SCZ25_META", "baseline_thickness")}. Alzheimer's disease scores are associated with thinning
({k("ALZ", "global_slope")}) only when the APOE region is retained ({k("ALZ_noAPOE", "global_slope")} without it); depression
{k("MDD_eur", "global_slope")} / {k("MDD_pooled", "global_slope")}; educational attainment {k("EA", "global_slope")}, opposite in sign; autism
{k("ASD", "global_slope")}. **g**, MAGMA competitive gene-set enrichment of the thinning rate
(red) and baseline thickness (grey) in genes of schizophrenia loci (curated,
Trubetskoy et al. 2022 Supplementary Table 12; {nsz["SCZ_locus_pool"]} genes) and depression loci
({nsz["MDD_pool"]:,} genes). European arm, 1000 Genomes European LD reference; {pooled_state}.
The schizophrenia set is enriched for both traits (p = {mp("SCZ_locus_pool", "global_slope")} thinning,
{mp("SCZ_locus_pool", "baseline_thickness")} baseline); the depression set for neither (p = {mp("MDD_pool", "global_slope")}, {mp("MDD_pool", "baseline_thickness")}).

Methods notes
- Thickness parsed from release FreeSurfer surfaces into HCP-MMP parcels; scans passing the release QC code; children with ≥ 2 scans.
- Slope model: thickness ~ age + sex + (1 + age | child) + (1 | site), fitted per parcel (a, d) and to the per-scan cortical mean (c, e–g).
- Polygenic scores: C+T, PRS-CS, SBayesRC (SBayesR in Supplementary Information); association model score + age + sex + 10 PCs + (1 | family).
- Sensitivity analyses (DK parcellation, mean of per-parcel slopes as the trait, PGC3 2022 discovery GWAS, per-ancestry strata, SBayesR) in Supplementary Information.
"""
    CAPTION.write_text(txt)
    print(CAPTION)
    return CAPTION


if __name__ == "__main__":
    draw()
    write_caption()
    sys.exit(0)
