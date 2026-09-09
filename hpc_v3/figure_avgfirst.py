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


def _density_panel(fig, rect, d, color, title, xlabel, ylabel):
    ax = fig.add_axes(rect)
    for _, r in d.iterrows():
        ax.add_patch(mpl.patches.Rectangle(
            (r.x_lo, r.y_lo), r.x_hi - r.x_lo, r.y_hi - r.y_lo,
            facecolor=color, edgecolor="none",
            alpha=min(1.0, 0.12 + 0.88 * (r["count"] / d["count"].max())**0.5)))
    lo = min(d.x_lo.min(), d.y_lo.min())
    hi = max(d.x_hi.max(), d.y_hi.max())
    ax.plot([lo, hi], [lo, hi], color="0.55", lw=0.7, ls=":")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=MID)
    ax.set_xlabel(xlabel, fontsize=SMALL, labelpad=1.5)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=SMALL)
    ax.tick_params(length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return ax


def _dumbbell_panel(fig, rect, cmp_, col_filled, col_open, title,
                    lab_filled, lab_open):
    ax = fig.add_axes(rect)
    names = ["global", "topDelta", "topC3"]
    for i, nm in enumerate(names):
        row = cmp_.loc[nm]
        y = len(names) - 1 - i
        v1, v2 = row[col_filled], row[col_open]
        ax.plot([v1, v2], [y, y], color="0.7", lw=1.2, zorder=1)
        ax.plot(v1, y, "o", color=PHENO_COLOR[nm], ms=6, zorder=3)
        ax.plot(v2, y, "o", mfc="white", mec=PHENO_COLOR[nm], mew=1.4, ms=6,
                zorder=3)
        for v, va in ((v1, "bottom"), (v2, "top")):
            ax.text(v, y + (0.14 if va == "bottom" else -0.14), f"{v:.3f}",
                    fontsize=SMALL - 1, ha="center", va=va, color="0.25")
    ax.set_yticks(range(len(names)), [DISPLAY[n] for n in names[::-1]],
                  fontsize=SMALL)
    ax.set_xlabel("lh–rh consistency (Spearman-Brown)", fontsize=SMALL,
                  labelpad=1.5)
    ax.set_title(f"{title}\nfilled = {lab_filled} · open = {lab_open}",
                 fontsize=MID)
    ax.set_ylim(-0.5, len(names) - 0.2)
    ax.tick_params(length=2)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    return ax


def draw() -> Path:
    dens = pd.read_csv(HERE / "avgfirst_density.csv")
    cmp_ = pd.read_csv(HERE / "avgfirst_comparison.csv").set_index("phenotype")
    dens_pc = pd.read_csv(HERE / "avgfirst_pc_density.csv")
    cmp_pc = (pd.read_csv(HERE / "avgfirst_pc_comparison.csv")
              .set_index("phenotype"))

    mpl.rcParams.update({
        "font.size": BASE, "axes.titlesize": MID, "axes.labelsize": SMALL,
        "xtick.labelsize": SMALL - 1, "ytick.labelsize": SMALL - 1,
        "axes.titlelocation": "left", "figure.facecolor": "white",
        "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(12.5, 7.6))

    # ---- row 1: avg-first vs mean-of-BLUPs ---------------------------------
    fig.text(0.005, 0.975, "Construction: averaging regions before the LMM "
             "does not beat averaging the per-region BLUPs — and is slightly "
             "less reliable", fontsize=BASE + 1)
    Y1, H = 0.585, 0.30
    for k, nm in enumerate(("global", "topDelta", "topC3")):
        row = cmp_.loc[nm]
        _density_panel(
            fig, [0.068 + k * 0.213, Y1, 0.178, H],
            dens[dens.phenotype == nm], PHENO_COLOR[nm],
            f"{DISPLAY[nm]}   r = {row.r_methods:.3f}",
            "mean of per-region BLUPs (mm/yr)",
            "avg-first: one LMM on\nregion-mean CT (mm/yr)" if k == 0 else "")
        fig.text(0.068 + k * 0.213 - 0.030, Y1 + H + 0.035, "abc"[k],
                 fontsize=BASE + 1, fontweight="bold")
    _dumbbell_panel(fig, [0.76, Y1, 0.225, H], cmp_,
                    "sb_blupmean", "sb_avgfirst", "reliability",
                    "mean of BLUPs", "avg-first")
    fig.text(0.76 - 0.028, Y1 + H + 0.035, "d", fontsize=BASE + 1,
             fontweight="bold")

    # ---- row 2: ancestry PCs in the LMM vs not ------------------------------
    n_pc = int(cmp_pc.n.iloc[0])
    fig.text(0.005, 0.470, "Covariates: adding ancestry PC1–10 (+ PC × age) "
             f"to the LMM leaves the slope phenotype unchanged (n = {n_pc:,})",
             fontsize=BASE + 1)
    Y2 = 0.075
    for k, nm in enumerate(("global", "topDelta", "topC3")):
        row = cmp_pc.loc[nm]
        _density_panel(
            fig, [0.068 + k * 0.213, Y2, 0.178, H],
            dens_pc[dens_pc.phenotype == nm], PHENO_COLOR[nm],
            f"{DISPLAY[nm]}   r = {row.r_pc_vs_nopc:.3f}",
            "avg-first slope, no PCs (mm/yr)",
            "avg-first slope,\nPC1–10 × age in LMM (mm/yr)" if k == 0 else "")
        fig.text(0.068 + k * 0.213 - 0.030, Y2 + H + 0.035, "efg"[k],
                 fontsize=BASE + 1, fontweight="bold")
    _dumbbell_panel(fig, [0.76, Y2, 0.225, H], cmp_pc,
                    "sb_nopc", "sb_pc", "reliability",
                    "no PCs (settled)", "with PCs")
    fig.text(0.76 - 0.028, Y2 + H + 0.035, "h", fontsize=BASE + 1,
             fontweight="bold")

    out = HERE / "figure_avgfirst.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(out)
    return out


if __name__ == "__main__":
    draw()
    sys.exit(0)
