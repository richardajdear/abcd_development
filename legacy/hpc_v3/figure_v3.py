"""One figure for the hpc_v3 phenotypes: brain maps + comparison matrices.

Usage (repo root):
    PYTHONPATH=src python hpc_v3/figure_v3.py --recompute   # rebuild CSVs from run dir
    PYTHONPATH=src python hpc_v3/figure_v3.py               # draw from committed CSVs

Follows the repo's provenance rule: the figure is drawn only from committed
summary CSVs in hpc_v3/, never from a run directory.  ``--recompute``
regenerates those CSVs (needs the settled run under out/):

  phenotypes_v3_maps.csv          68 labels x: loading map per phenotype
                                  (Pearson r between the subject-level score
                                  and each region's slope BLUP, n=8192),
                                  plus group_delta, C3 (mirrored), membership
  phenotypes_v3_spatial_corr.csv  pairwise spatial Spearman rho + spin p
                                  (5000 rotations) among the 9 loading maps

The phenotypic (subject-space) matrix comes from
phenotypes_v3_characterization.csv, written by make_phenotypes_v3.py.

Output: hpc_v3/phenotypes_v3_figure.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

PHENOS = ["baseline_thickness", "global_slope", "slope_PC1", "slope_PC2",
          "slope_PC3", "slope_topDelta", "slope_topC3", "slope_projDelta",
          "slope_projC3"]
SHORT = ["baseline", "global", "PC1", "PC2", "PC3",
         "topΔCT", "topC3", "projΔCT", "projC3"]
DISPLAY = dict(zip(PHENOS, [
    "baseline thickness", "global slope", "slope PC1", "slope PC2",
    "slope PC3", "top-ΔCT mean", "top-C3 mean", "ΔCT projection",
    "C3 projection"]))

# membership encoding in phenotypes_v3_maps.csv
MEMBER_COLORS = {1: "#0072B2", 2: "#CC79A7", 3: "#E69F00"}   # ΔCT / both / C3
MEMBER_NAMES = {1: "top-ΔCT only", 2: "both sets", 3: "top-C3 only"}


def recompute() -> None:
    """Rebuild the committed CSVs from the settled run directory."""
    import make_phenotypes_v3 as mp
    from abcd import spatial
    from abcd.config import active_run_dir
    from abcd.gcta_export import subject_phenotypes

    run_dir = active_run_dir()
    sl = mp.slope_matrix(run_dir)                       # 8192 x 68
    delta, c3 = mp.region_rankings()
    new, _regions = mp.build_new_columns(run_dir)
    scores = subject_phenotypes(run_dir).join(new)[PHENOS]
    assert scores.index.equals(sl.index)

    # loading maps: corr(score, regional slope) per region, all 9 phenotypes
    slz = (sl - sl.mean()) / sl.std()
    scz = (scores - scores.mean()) / scores.std()
    load = pd.DataFrame(slz.to_numpy().T @ scz.to_numpy() / (len(sl) - 1),
                        index=sl.columns, columns=PHENOS)

    maps = load.copy()
    maps["group_delta"] = pd.Series(
        {f"{h}_{b}": v for b, v in delta.items() for h in ("lh", "rh")})
    maps["C3"] = pd.Series(
        {f"{h}_{b}": v for b, v in c3.items() for h in ("lh", "rh")})
    set_d = set(delta.index[:mp.K_BILATERAL])
    set_c = set(c3.sort_values(ascending=False).index[:mp.K_BILATERAL])
    member = {}
    for lab in sl.columns:
        base = lab.split("_", 1)[1]
        code = (2 if base in set_d & set_c else
                1 if base in set_d else 3 if base in set_c else np.nan)
        member[lab] = code
    maps["membership"] = pd.Series(member)
    maps.index.name = "label"
    maps.round(6).to_csv(HERE / "phenotypes_v3_maps.csv")

    # spatial correlations among the loading maps, spin null
    geom = spatial.load_dk_geometry()
    rows = []
    for i, a in enumerate(PHENOS):
        for b in PHENOS[i + 1:]:
            t = spatial.spatial_corr(load[a], load[b], geom=geom,
                                     n_perm=5000, seed=42)
            rows.append(dict(map_x=a, map_y=b, rho=t.r, p_spin=t.p))
    pd.DataFrame(rows).round(6).to_csv(
        HERE / "phenotypes_v3_spatial_corr.csv", index=False)
    print("recomputed phenotypes_v3_maps.csv and phenotypes_v3_spatial_corr.csv")


# ---------------------------------------------------------------- plotting
# Slide-deck dimensions (project convention: 13.333 x 7.5 in, dpi 200).
# Layout: panel a = single column of definition maps (left); panel b = the
# nine loading maps in two columns; right block = k-sensitivity (c) over the
# two comparison matrices (d, e).

BASE, MID, SMALL = 9.5, 8.5, 7.0


def _style() -> None:
    mpl.rcParams.update({
        "font.size": BASE, "axes.titlesize": BASE, "axes.labelsize": BASE,
        "legend.fontsize": MID, "xtick.labelsize": SMALL,
        "ytick.labelsize": SMALL, "axes.titlelocation": "left",
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })


def _panel_letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=BASE + 3, fontweight="bold", va="top")


def _matrix(ax, M, title, star=None, annot_size=SMALL - 1.5):
    """Annotated symmetric correlation matrix on a +/-1 diverging scale."""
    n = len(M)
    im = ax.imshow(M.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(n), SHORT, rotation=45, ha="right", fontsize=SMALL - 1)
    ax.set_yticks(range(n), SHORT, fontsize=SMALL - 1)
    ax.set_title(title, pad=4, fontsize=MID)
    for i in range(n):
        for j in range(n):
            v = M.iloc[i, j]
            if i == j:
                continue
            txt = f"{v:.2f}".replace("0.", ".").replace("-.",  "\u2212.")
            if star is not None and star.iloc[i, j]:
                txt += "*"
            ax.text(j, i, txt, ha="center", va="center", fontsize=annot_size,
                    color="white" if abs(v) > 0.6 else "black")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def draw() -> Path:
    from abcd import brainplot

    maps = pd.read_csv(HERE / "phenotypes_v3_maps.csv").set_index("label")
    char = pd.read_csv(HERE / "phenotypes_v3_characterization.csv",
                       index_col=0).loc[PHENOS, PHENOS]
    sp = pd.read_csv(HERE / "phenotypes_v3_spatial_corr.csv")
    S = pd.DataFrame(np.eye(len(PHENOS)), index=PHENOS, columns=PHENOS)
    P = S.copy()
    for _, r in sp.iterrows():
        S.loc[r.map_x, r.map_y] = S.loc[r.map_y, r.map_x] = r.rho
        P.loc[r.map_x, r.map_y] = P.loc[r.map_y, r.map_x] = r.p_spin
    star = P < 0.05

    cons = pd.read_csv(HERE / "phenotypes_v3_consistency.csv")

    _style()
    fig = plt.figure(figsize=(13.333, 7.5))

    def brain(rect, series, title, **kw):
        ax = fig.add_axes(rect)
        _, _, m = brainplot.plot_dk(series, ax=ax, colorbar=False,
                                    fontsize=SMALL, **kw)
        ax.set_title(title, fontsize=MID, pad=1.5)
        return ax, m

    def hcbar(rect, mappable, label):
        cax = fig.add_axes(rect)
        cb = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        cb.set_label(label, fontsize=SMALL, labelpad=1.5)
        cb.ax.tick_params(labelsize=SMALL - 1, length=2)
        cb.outline.set_visible(False)

    # --- panel a: definitions, single column on the left ---------------------
    AX_W, AX_H = 0.155, 0.215                     # one brain-map slot
    _panel_letter(fig, 0.006, 0.985, "a")
    fig.text(0.022, 0.972, "definitions (top-8 bilateral sets)",
             fontsize=BASE)
    _, m1 = brain([0.012, 0.715, AX_W, AX_H], maps["group_delta"],
                  "group-mean thinning rate", diverging=True, center=0.0)
    hcbar([0.037, 0.700, 0.105, 0.010], m1, "ΔCT (mm/year)")
    _, m2 = brain([0.012, 0.415, AX_W, AX_H], maps["C3"],
                  "AHBA C3 (lh mirrored)", diverging=True, center=0.0)
    hcbar([0.037, 0.400, 0.105, 0.010], m2, "C3 score (z)")
    cmap = mpl.colors.ListedColormap([MEMBER_COLORS[k] for k in (1, 2, 3)])
    ax3, _ = brain([0.012, 0.125, AX_W, AX_H], maps["membership"],
                   "selection sets (3 of 8 shared)", cmap=cmap,
                   vminmax=(0.5, 3.5))
    handles = [mpl.patches.Patch(facecolor=MEMBER_COLORS[k],
                                 label=MEMBER_NAMES[k]) for k in (1, 2, 3)]
    ax3.legend(handles=handles, loc="upper center", ncol=1, frameon=False,
               bbox_to_anchor=(0.5, -0.02), handlelength=1.1,
               labelspacing=0.25, fontsize=SMALL)

    # --- panel b: loading maps, two columns -----------------------------------
    load = maps[PHENOS]
    lim = float(np.abs(load.to_numpy()).max())
    _panel_letter(fig, 0.185, 0.985, "b")
    fig.text(0.202, 0.972, "loading maps: r(score, regional thinning rate), "
             "n = 8,192", fontsize=BASE)
    X0, PITCH_X, PITCH_Y = 0.195, 0.170, 0.185
    BW, BH = 0.160, 0.158                       # map slot; title uses the rest
    # column-wise: settled five down the left column, the four new on the right
    mappable = None
    for k, ph in enumerate(PHENOS):
        col, row = k // 5, k % 5
        rect = [X0 + col * PITCH_X, 0.790 - row * PITCH_Y, BW, BH]
        _, mappable = brain(rect, load[ph], DISPLAY[ph], diverging=True,
                            center=0.0, vminmax=(-lim, lim))
    # shared colorbar in the empty 10th slot
    hcbar([X0 + PITCH_X + 0.025, 0.095, 0.110, 0.010],
          mappable, "r (score vs regional thinning rate)")

    # --- panel c: k-sensitivity ------------------------------------------------
    axk = fig.add_axes([0.635, 0.615, 0.348, 0.30])
    _panel_letter(fig, 0.578, 0.985, "c")
    axk.set_title("choosing k: reliability falls before\ncollinearity does",
                  fontsize=MID, pad=3)
    kcol = {"topDelta": MEMBER_COLORS[1], "topC3": MEMBER_COLORS[3]}
    for ph in ("topDelta", "topC3"):
        d = cons[cons.phenotype == ph].sort_values("k_bilateral")
        axk.plot(d.k_bilateral, d.r_with_global, "o-", color=kcol[ph],
                 lw=1.6, ms=3.5)
        axk.plot(d.k_bilateral, d.spearman_brown, "o--", color=kcol[ph],
                 lw=1.3, ms=3.5, alpha=0.85)
    g_sb = float(cons.loc[cons.phenotype == "global_slope",
                          "spearman_brown"].iloc[0])
    axk.axhline(g_sb, color="0.55", lw=0.8, ls=":")
    axk.axvspan(7.2, 8.8, color="0.90", zorder=0)
    axk.set_ylim(0.695, 1.0)
    axk.set_xlim(3.2, 17.8)
    axk.text(8, 0.703, "k = 8", ha="center", fontsize=SMALL, color="0.35")
    axk.text(16.8, g_sb + 0.006, f"global-slope consistency ({g_sb:.2f})",
             ha="right", va="bottom", fontsize=SMALL, color="0.35")
    axk.text(0.03, 0.97, "— r with global slope (collinearity; 1 at k = 34)",
             transform=axk.transAxes, ha="left", va="top", fontsize=SMALL)
    axk.text(0.03, 0.90, "- - lh–rh consistency (reliability)",
             transform=axk.transAxes, ha="left", va="top", fontsize=SMALL)
    axk.set_xlabel("regions per set (k, bilateral)", fontsize=SMALL,
                   labelpad=1.5)
    axk.set_ylabel("correlation", fontsize=SMALL, labelpad=1.5)
    axk.set_xticks([4, 6, 8, 10, 12, 17])
    axk.tick_params(labelsize=SMALL - 1)
    axk.spines[["top", "right"]].set_visible(False)
    leg = [mpl.lines.Line2D([], [], color=kcol["topDelta"], lw=1.6,
                            label="top-ΔCT"),
           mpl.lines.Line2D([], [], color=kcol["topC3"], lw=1.6,
                            label="top-C3")]
    axk.legend(handles=leg, frameon=False, fontsize=SMALL,
               loc="lower right", handlelength=1.3, borderaxespad=0.2)

    # --- panels d, e: matrices --------------------------------------------------
    axc = fig.add_axes([0.635, 0.115, 0.155, 0.36])
    axd = fig.add_axes([0.828, 0.115, 0.155, 0.36])
    _matrix(axc, char, "subject space:\ncorrelation of scores")
    _matrix(axd, S, "map space: spatial corr.\nof loading maps*", star=star)
    _panel_letter(fig, 0.578, 0.53, "d")
    _panel_letter(fig, 0.795, 0.53, "e")
    fig.text(0.828, 0.022, "* spin-test p < 0.05 (5,000 rotations)",
             fontsize=SMALL - 1)

    out = HERE / "phenotypes_v3_figure.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(out)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--recompute", action="store_true",
                    help="rebuild the committed CSVs from the run directory")
    a = ap.parse_args(argv)
    if a.recompute:
        recompute()
    draw()
    return 0


if __name__ == "__main__":
    sys.exit(main())
