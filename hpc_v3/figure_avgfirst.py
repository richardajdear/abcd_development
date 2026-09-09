"""Figure: avg-first vs mean-of-BLUPs construction for the slope phenotypes.

Usage (repo root):  python hpc_v3/figure_avgfirst.py

Draws ONLY from the committed summary tables written by make_avgfirst.py
(avgfirst_density.csv, avgfirst_comparison.csv) — no per-subject data.
Panels a-c: binned density of the two constructions per phenotype, with r
annotated from the table.  Panel d: lh/rh Spearman-Brown consistency of the
two constructions side by side — the direct answer to "would averaging first
raise reliability (and so heritability)?": no, it is marginally lower.

Output: hpc_v3/figure_avgfirst.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

BASE, MID, SMALL = 10.0, 9.0, 7.5
DISPLAY = {"global": "global slope", "topDelta": "top-ΔCT mean",
           "topC3": "top-C3 mean"}
PHENO_COLOR = {"global": "0.25", "topDelta": "#0072B2", "topC3": "#E69F00"}


def draw() -> Path:
    dens = pd.read_csv(HERE / "avgfirst_density.csv")
    cmp_ = pd.read_csv(HERE / "avgfirst_comparison.csv").set_index("phenotype")

    mpl.rcParams.update({
        "font.size": BASE, "axes.titlesize": MID, "axes.labelsize": SMALL,
        "xtick.labelsize": SMALL - 1, "ytick.labelsize": SMALL - 1,
        "axes.titlelocation": "left", "figure.facecolor": "white",
        "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(12.5, 3.7))
    fig.text(0.005, 0.955, "Averaging regions before the LMM does not beat "
             "averaging the per-region BLUPs — and is slightly less reliable",
             fontsize=BASE + 1)

    for k, nm in enumerate(("global", "topDelta", "topC3")):
        ax = fig.add_axes([0.068 + k * 0.213, 0.16, 0.178, 0.62])
        d = dens[dens.phenotype == nm]
        # draw binned counts as a quadmesh on the bin rectangles
        x = np.sort(d.x_lo.unique())
        for _, r in d.iterrows():
            ax.add_patch(mpl.patches.Rectangle(
                (r.x_lo, r.y_lo), r.x_hi - r.x_lo, r.y_hi - r.y_lo,
                facecolor=PHENO_COLOR[nm],
                alpha=min(1.0, 0.12 + 0.88 * (r["count"] / d["count"].max())**0.5),
                edgecolor="none"))
        lo = min(d.x_lo.min(), d.y_lo.min())
        hi = max(d.x_hi.max(), d.y_hi.max())
        ax.plot([lo, hi], [lo, hi], color="0.55", lw=0.7, ls=":")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        row = cmp_.loc[nm]
        ax.set_title(f"{DISPLAY[nm]}   r = {row.r_methods:.3f}", fontsize=MID)
        ax.set_xlabel("mean of per-region BLUPs (mm/yr)")
        if k == 0:
            ax.set_ylabel("avg-first: one LMM on\nregion-mean CT (mm/yr)",
                          fontsize=SMALL)
        ax.tick_params(length=2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        fig.text(0.068 + k * 0.213 - 0.030, 0.90, "abc"[k],
                 fontsize=BASE + 1, fontweight="bold")

    # panel d: consistency dumbbells
    ax = fig.add_axes([0.76, 0.16, 0.225, 0.62])
    names = ["global", "topDelta", "topC3"]
    for i, nm in enumerate(names):
        row = cmp_.loc[nm]
        y = len(names) - 1 - i
        ax.plot([row.sb_blupmean, row.sb_avgfirst], [y, y], color="0.7",
                lw=1.2, zorder=1)
        ax.plot(row.sb_blupmean, y, "o", color=PHENO_COLOR[nm], ms=6,
                zorder=3, label="mean of BLUPs" if i == 0 else None)
        ax.plot(row.sb_avgfirst, y, "o", mfc="white", mec=PHENO_COLOR[nm],
                mew=1.4, ms=6, zorder=3,
                label="avg-first" if i == 0 else None)
        for v, va in ((row.sb_blupmean, "bottom"), (row.sb_avgfirst, "top")):
            ax.text(v, y + (0.13 if va == "bottom" else -0.13), f"{v:.3f}",
                    fontsize=SMALL - 1, ha="center", va=va, color="0.25")
    ax.set_yticks(range(len(names)), [DISPLAY[n] for n in names[::-1]],
                  fontsize=SMALL)
    ax.set_xlabel("lh–rh consistency (Spearman-Brown)")
    ax.set_title("reliability: filled = settled construction", fontsize=MID)
    ax.set_ylim(-0.5, len(names) - 0.2)
    ax.tick_params(length=2)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    fig.text(0.76 - 0.028, 0.90, "d", fontsize=BASE + 1, fontweight="bold")

    out = HERE / "figure_avgfirst.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(out)
    return out


if __name__ == "__main__":
    draw()
    sys.exit(0)
