"""Step 4 (laptop): the integrating figure -- SCZ polygenic score -> parcel
thinning slope, against the transcriptomic maps, on HCP-MMP.

  a  the beta map: per-parcel beta of the thinning slope on the primary SCZ
     score (SD per SD; negative = higher score, faster thinning), left
     hemisphere, lateral and medial
  b  AHBA C3 (Dear et al. 2024) on the same parcels (137 AHBA-covered)
  c  the ahba_pls lead component (thinning-oriented) on the same parcels
  d  beta vs C3, Spearman rho and spin p (5,000 rotations of the complete
     beta map)
  e  beta vs the lead component, likewise
  f  left vs right hemisphere beta per Glasser area: how reproducible the
     beta map itself is

Reads results/beta_map_<cell>.tsv, results/spin_tests.tsv, work/spin_null_<cell>.npz,
ahba_pls/data/hcp_polygons.csv (ggsegGlasser flat polygons, lh).  Writes
results/fig_prs_beta_map_<cell>.png and results/fig_prs_beta_map_<cell>.md (caption).

Usage (repo root):
  python prs_beta_map/prs_scz_beta_map_hcp/04_figure.py [--cell SCZ25_META_SBayesRC_zanc] [--variant lh]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy import stats

mpl.use("Agg")
REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from importlib import import_module  # noqa: E402
spin = import_module("03_spin_test")

POLY = REPO / "ahba_pls" / "data" / "hcp_polygons.csv"

# -- dataviz palette (references/palette.md): diverging blue <-> red, gray midpoint
BLUE, RED, MID = "#2a78d6", "#e34948", "#f0efec"
BLUE_DEEP, RED_DEEP = "#184f95", "#a82f2e"
INK, INK2, INK3, HAIR = "#0b0b0b", "#52514e", "#8a8983", "#dcdbd6"
DIV = LinearSegmentedColormap.from_list("div", [BLUE_DEEP, BLUE, MID, RED, RED_DEEP])

W_IN, H_IN = 180 / 25.4, 118 / 25.4
BASE, ANN, TICK = 7, 6, 5.5


def _style():
    mpl.rcParams.update({
        "font.size": BASE, "axes.titlesize": BASE, "axes.labelsize": BASE,
        "xtick.labelsize": TICK, "ytick.labelsize": TICK, "legend.fontsize": ANN,
        "axes.linewidth": 0.5, "axes.edgecolor": INK3, "xtick.color": INK3, "ytick.color": INK3,
        "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "axes.labelcolor": INK,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.5,
        "ytick.major.size": 2.5, "figure.facecolor": "white", "savefig.facecolor": "white",
        "axes.spines.top": False, "axes.spines.right": False, "font.family": "sans-serif",
    })


def letter(ax, s, dx=-0.08, dy=1.02):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=BASE + 2, fontweight="bold",
            va="bottom", ha="left", color=INK)


# ---------------------------------------------------------------- flat map --
def load_polys() -> pd.DataFrame:
    p = pd.read_csv(POLY)
    # lay the medial view under the lateral one, both centred
    out = []
    for view, off in (("lateral", 0.0), ("medial", 1.0)):
        v = p[p.view == view].copy()
        v["x"] = (v.x - v.x.min()) / (v.x.max() - v.x.min())
        v["y"] = (v.y - v.y.min()) / (v.y.max() - v.y.min()) * 0.62 - off * 0.72
        out.append(v)
    return pd.concat(out)


def draw_map(ax, values: pd.Series, polys: pd.DataFrame, vlim: float, title: str,
             cbar_label: str, fig, missing_note: str | None = None):
    norm = TwoSlopeNorm(vmin=-vlim, vcenter=0, vmax=vlim)
    verts, cols = [], []
    for (label, view, g, sg), q in polys.groupby(["label", "view", "group", "subgroup"], sort=False):
        verts.append(q[["x", "y"]].to_numpy())
        v = values.get(label, np.nan)
        cols.append(DIV(norm(v)) if np.isfinite(v) else (0.90, 0.90, 0.89, 1.0))
    pc = PolyCollection(verts, facecolors=cols, edgecolors="white", linewidths=0.25)
    ax.add_collection(pc)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-1.0, 0.66); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=BASE, color=INK, pad=3, loc="left", x=0.04)
    ax.text(0.98, 0.64, "lateral", fontsize=TICK, color=INK3, ha="right", va="top")
    ax.text(0.98, -0.10, "medial", fontsize=TICK, color=INK3, ha="right", va="top")
    cax = ax.inset_axes([0.06, 0.0, 0.40, 0.03])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=DIV), cax=cax, orientation="horizontal")
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=TICK, length=2, width=0.5, color=INK3, pad=1)
    cb.set_ticks([-vlim, 0, vlim]); cb.set_ticklabels([f"−{vlim:.2g}", "0", f"+{vlim:.2g}"])
    cb.set_label(cbar_label, fontsize=TICK, color=INK2, labelpad=1)
    if missing_note:
        ax.text(0.52, -0.975, missing_note, fontsize=TICK, color=INK3, va="center", ha="left")


# ----------------------------------------------------------------- scatter --
def scatter(ax, x: pd.Series, y: pd.Series, row: pd.Series, xlabel: str, ylabel: str,
            null: np.ndarray | None, cond: pd.Series | None = None):
    common = x.index.intersection(y.index)
    xv, yv = x.loc[common].to_numpy(float), y.loc[common].to_numpy(float)
    ax.axhline(0, color=HAIR, lw=0.5, zorder=0); ax.axvline(0, color=HAIR, lw=0.5, zorder=0)
    ax.scatter(xv, yv, s=9, color=INK2, edgecolor="white", linewidth=0.4, zorder=2)
    b, a0 = np.polyfit(xv, yv, 1)
    xx = np.linspace(xv.min(), xv.max(), 20)
    ax.plot(xx, a0 + b * xx, color=BLUE, lw=1.2, zorder=3)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.30 * (hi - lo))
    p_txt = f"p_spin = {row.p_spin:.3f}" if row.p_spin >= 0.001 else "p_spin < 0.001"
    ax.text(0.03, 0.97, f"Spearman ρ = {row.rho:+.2f}\n{p_txt}\n{int(row.n_parcels)} parcels",
            transform=ax.transAxes, va="top", ha="left", fontsize=ANN, color=INK)
    if cond is not None:
        pc = f"p_spin {cond.p_spin:.3f}" if cond.p_spin >= 0.001 else "p_spin < 0.001"
        ax.text(0.97, 0.03, f"whole-cortex slope held fixed:\nρ = {cond.rho:+.2f}, {pc}",
                transform=ax.transAxes, va="bottom", ha="right", fontsize=TICK, color=INK2)
    if null is not None:   # the spin null of rho, as a strip along the top edge
        ins = ax.inset_axes([0.62, 0.84, 0.36, 0.09])
        ins.hist(null, bins=40, color=HAIR, lw=0)
        ins.axvline(row.rho, color=BLUE, lw=1.0)
        ins.set_xlim(-1, 1); ins.set_yticks([]); ins.set_xticks([-1, 0, 1])
        ins.tick_params(labelsize=TICK - 1, length=1.5, pad=1, colors=INK3)
        for s in ins.spines.values(): s.set_visible(False)
        ins.text(0.5, 1.15, "spin null of ρ", transform=ins.transAxes, fontsize=TICK - 0.5,
                 color=INK3, ha="center")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default="SCZ25_META_SBayesRC_zanc")
    ap.add_argument("--variant", default="lh", choices=["lh", "bilateral_mean"])
    ap.add_argument("--results", type=Path, default=HERE / "results")
    ap.add_argument("--work", type=Path, default=HERE / "work")
    a = ap.parse_args(argv)
    _style()

    bm = pd.read_csv(a.results / f"beta_map_{a.cell}.tsv", sep="\t")
    st = pd.read_csv(a.results / "spin_tests.tsv", sep="\t")
    st_all = st[(st.cell == a.cell) & (st.beta_variant == a.variant)]
    st = st_all[st_all["map"] == "beta"].set_index("reference")
    st_cond = st_all[st_all["map"] == "beta_cond"].set_index("reference")
    st_load = st_all[st_all["map"] == "r_global"].set_index("reference")
    nulls = np.load(a.work / f"spin_null_{a.cell}.npz")
    refs = spin.references()
    parc = bm[bm.label != "cortex_mean"]
    maps = spin.split_map(parc)
    beta = maps[a.variant]
    cortex = bm[bm.label == "cortex_mean"].iloc[0]
    polys = load_polys()
    covered = refs["C3"].dropna().index

    fig = plt.figure(figsize=(W_IN, H_IN))
    gs = fig.add_gridspec(2, 3, left=0.075, right=0.985, top=0.92, bottom=0.09,
                          hspace=0.40, wspace=0.40, height_ratios=[1.0, 0.95])
    axm = [fig.add_subplot(gs[0, i]) for i in range(3)]
    axs = [fig.add_subplot(gs[1, i]) for i in range(3)]

    vb = float(np.nanpercentile(np.abs(beta), 98))
    arm = "pooled arm" if cortex.stratum == "full" else "European arm"
    draw_map(axm[0], beta, polys, vb,
             f"SCZ score → thinning slope, per parcel\n({arm}, n = {int(cortex.n):,})",
             "β, SD slope per SD score", fig)
    axm[0].text(0.56, -0.90, f"whole cortex β = {cortex.beta:+.3f}\n(SE {cortex.se:.3f}, p = {cortex.p:.2g})",
                fontsize=TICK, color=INK2, va="center")
    c3 = refs["C3"].reindex(covered); lead = refs["PLS_lead"].reindex(covered)
    draw_map(axm[1], c3, polys, float(np.nanpercentile(np.abs(c3), 98)),
             "AHBA C3\n(Dear et al. 2024)", "component score", fig,
             missing_note="grey: no AHBA coverage")
    draw_map(axm[2], lead, polys, float(np.nanpercentile(np.abs(lead), 98)),
             "ABCD lead component\n(ahba_pls, thinning-oriented)", "regional score", fig,
             missing_note="grey: no AHBA coverage")
    for ax, L in zip(axm, "abc"):
        letter(ax, L, dx=-0.04, dy=1.06)

    ylab = "β (SD slope / SD score)" if a.variant == "lh" else "β, lh–rh mean"
    scatter(axs[0], c3, beta, st.loc["C3"], "AHBA C3 score", ylab,
            nulls[f"beta__{a.variant}__C3"], st_cond.loc["C3"])
    scatter(axs[1], lead, beta, st.loc["PLS_lead"], "lead component score", "",
            nulls[f"beta__{a.variant}__PLS_lead"], st_cond.loc["PLS_lead"])
    lh, rh = maps["lh"], maps["rh_as_lh"].reindex(maps["lh"].index)
    rlr = stats.spearmanr(lh, rh).statistic
    ax = axs[2]
    ax.axhline(0, color=HAIR, lw=0.5, zorder=0); ax.axvline(0, color=HAIR, lw=0.5, zorder=0)
    lim = float(np.nanmax(np.abs(np.r_[lh, rh]))) * 1.08
    ax.plot([-lim, lim], [-lim, lim], color=HAIR, lw=0.6, zorder=1)
    ax.scatter(lh, rh, s=9, color=INK2, edgecolor="white", linewidth=0.4, zorder=2)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.set_xlabel("β, left-hemisphere parcel"); ax.set_ylabel("β, homologous right parcel")
    ax.text(0.03, 0.97, f"Spearman ρ = {rlr:+.2f}\n(179 areas)", transform=ax.transAxes,
            va="top", ha="left", fontsize=ANN, color=INK)
    for ax, L in zip(axs, "def"):
        letter(ax, L, dx=-0.28, dy=1.0)

    if a.cell.upper().startswith("SYNTHETIC"):
        fig.text(0.5, 0.5, "SYNTHETIC SCORE — PIPELINE TEST, NOT A RESULT", fontsize=22,
                 color=RED, alpha=0.35, ha="center", va="center", rotation=20, zorder=10)
    out = a.results / f"fig_prs_beta_map_{a.cell}.png"
    fig.savefig(out, dpi=300)
    # -- caption from the numbers
    r = st
    cap = (f"**SCZ polygenic score → parcel thinning slope, against the transcriptomic maps "
           f"(HCP-MMP, {arm}).** (a) β of each parcel's thinning slope (SD) on the score (SD), from "
           f"slope ~ score + sex + age + 10 PCs + (1 | family); n = {int(cortex.n):,} children, whole-cortex "
           f"β {cortex.beta:+.3f} (SE {cortex.se:.3f}, p = {cortex.p:.2g}); left hemisphere shown"
           f"{'' if a.variant == 'lh' else ' (lh–rh mean)'}; {int(r.n_parcels_p05.iloc[0])} of 358 parcels p < 0.05, "
           f"{int(r.n_parcels_beta_neg.iloc[0])} negative. (b) AHBA C3, (c) the ahba_pls lead component, on the "
           f"{int(r.loc['C3', 'n_parcels'])} AHBA-covered parcels. (d) β vs C3: ρ = {r.loc['C3', 'rho']:+.2f}, "
           f"p_spin = {r.loc['C3', 'p_spin']:.3f}; (e) β vs the lead component: ρ = {r.loc['PLS_lead', 'rho']:+.2f}, "
           f"p_spin = {r.loc['PLS_lead', 'p_spin']:.3f} ({int(r.n_perm.iloc[0]):,} rotations of the complete "
           f"179-parcel β map, bijective assignment). Context: β vs the group thinning map ρ = "
           f"{r.loc['dCT', 'rho']:+.2f} (p_spin {r.loc['dCT', 'p_spin']:.3f}). With the whole-cortex mean slope as a "
           f"covariate (parcel-specific β): vs C3 ρ = {st_cond.loc['C3', 'rho']:+.2f} (p_spin "
           f"{st_cond.loc['C3', 'p_spin']:.3f}), vs the lead component ρ = {st_cond.loc['PLS_lead', 'rho']:+.2f} "
           f"(p_spin {st_cond.loc['PLS_lead', 'p_spin']:.3f}); each parcel's loading on the cortex mean alone vs C3: "
           f"ρ = {st_load.loc['C3', 'rho']:+.2f} (p_spin {st_load.loc['C3', 'p_spin']:.3f}). (f) left vs right β per area, "
           f"ρ = {rlr:+.2f}. Score: {cortex.score}.")
    (a.results / f"fig_prs_beta_map_{a.cell}.md").write_text(cap + "\n")
    print(f"wrote {out}\n{cap}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
