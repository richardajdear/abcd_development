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


def _matrix(ax, M, title, star=None):
    """Annotated symmetric correlation matrix on a +/-1 diverging scale."""
    n = len(M)
    im = ax.imshow(M.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(n), SHORT, rotation=45, ha="right")
    ax.set_yticks(range(n), SHORT)
    ax.set_title(title, pad=8, fontsize=BASE)
    for i in range(n):
        for j in range(n):
            v = M.iloc[i, j]
            if i == j:
                continue
            txt = f"{v:.2f}".replace("0.", ".").replace("-.",  "\u2212.")
            if star is not None and star.iloc[i, j]:
                txt += "*"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=SMALL - 0.5,
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

    _style()
    fig = plt.figure(figsize=(13.2, 15.5))
    # maps block: row a (definitions) + 3 rows of loading maps
    gs = fig.add_gridspec(4, 3, hspace=0.42, wspace=0.05,
                          left=0.03, right=0.985, top=0.925, bottom=0.36)

    # --- row a: definitions -------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    brainplot.plot_dk(maps["group_delta"], ax=ax, diverging=True, center=0.0,
                      label="mm/year", fontsize=SMALL)
    ax.set_title("group-mean thinning rate (ΔCT)", fontsize=MID)

    ax = fig.add_subplot(gs[0, 1])
    brainplot.plot_dk(maps["C3"], ax=ax, diverging=True, center=0.0,
                      label="AHBA C3 (z)", fontsize=SMALL)
    ax.set_title("AHBA C3 score (lh mirrored)", fontsize=MID)

    ax = fig.add_subplot(gs[0, 2])
    memb = maps["membership"]
    cmap = mpl.colors.ListedColormap([MEMBER_COLORS[k] for k in (1, 2, 3)])
    brainplot.plot_dk(memb, ax=ax, cmap=cmap, vminmax=(0.5, 3.5),
                      colorbar=False, fontsize=SMALL)
    ax.set_title("selected regions (top-8 bilateral each)", fontsize=MID)
    handles = [mpl.patches.Patch(facecolor=MEMBER_COLORS[k],
                                 label=MEMBER_NAMES[k]) for k in (1, 2, 3)]
    ax.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
              bbox_to_anchor=(0.5, -0.18), handlelength=1.2,
              columnspacing=0.9, fontsize=SMALL)

    fig.text(0.03, 0.952,
             "Fast-thinning and high-C3 cortex overlap in 3 of 8 regions",
             fontsize=BASE)
    _panel_letter(fig, 0.008, 0.962, "a")

    # --- rows b: loading maps, shared scale ---------------------------------
    load = maps[PHENOS]
    lim = float(np.abs(load.to_numpy()).max())
    mappable, axes_b = None, []
    for k, ph in enumerate(PHENOS):
        ax = fig.add_subplot(gs[1 + k // 3, k % 3])
        _, _, mappable = brainplot.plot_dk(
            load[ph], ax=ax, diverging=True, center=0.0,
            vminmax=(-lim, lim), colorbar=False, fontsize=SMALL)
        ax.set_title(DISPLAY[ph], fontsize=MID)
        axes_b.append(ax)

    # header sits in the gap above row 1, derived from real axes geometry
    b_top = max(ax.get_position().y1 for ax in axes_b[:3])
    fig.text(0.03, b_top + 0.030, "What each phenotype measures: correlation "
             "of the score with each region's thinning rate (n = 8,192)",
             fontsize=BASE)
    _panel_letter(fig, 0.008, b_top + 0.040, "b")

    # one shared colorbar just below the block, clear of the matrices
    b_bot = min(ax.get_position().y0 for ax in axes_b[-3:])
    cax = fig.add_axes([0.40, b_bot - 0.035, 0.22, 0.008])
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal")
    cb.set_label("r (phenotype vs regional thinning rate)", fontsize=SMALL)
    cb.ax.tick_params(labelsize=SMALL)

    # --- row c/d: matrices ---------------------------------------------------
    axc = fig.add_axes([0.10, 0.035, 0.33, 0.20])
    axd = fig.add_axes([0.585, 0.035, 0.33, 0.20])
    _matrix(axc, char,
            "subject space: subset means are near-collinear with global slope")
    _matrix(axd, S, "map space: spatial correlation of loading maps",
            star=star)
    top_c = axc.get_position().y1
    _panel_letter(fig, 0.008, top_c + 0.035, "c")
    _panel_letter(fig, 0.50, top_c + 0.035, "d")
    fig.text(0.585, 0.006, "* spin-test p < 0.05 (5,000 rotations)",
             fontsize=SMALL)

    out = HERE / "phenotypes_v3_figure.png"
    fig.savefig(out, dpi=300)
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
