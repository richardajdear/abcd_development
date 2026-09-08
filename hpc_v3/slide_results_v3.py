"""Experiment A collection slide: h2, PRS, MAGMA, rg for the six phenotypes.

Usage (repo root):  PYTHONPATH=src python hpc_v3/slide_results_v3.py

Reads ONLY hpc_v3/results_v3_summary.csv (built by make_results_v3.py from
the committed result tables).  Rows with status=pending — the four new
phenotypes until the CSD3 run lands — are drawn as grey open squares with
dashed intervals; re-run make_results_v3.py after pulling the new result
files and this slide updates itself.  No statistic is hardcoded here: p
annotations, the m_eff correction, and the underpowered daggers are all
computed from the table at plot time.

Output: hpc_v3/slide_results_v3.png  (13.333 x 7.5 in, dpi 200)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

BASE, MID, SMALL = 10.5, 9.0, 7.5

PHENOS = ["baseline_thickness", "global_slope", "slope_topDelta",
          "slope_topC3", "slope_projDelta", "slope_projC3"]
YLAB = ["baseline thickness\n(control)", "global slope", "top-ΔCT mean",
        "top-C3 mean", "ΔCT projection", "C3 projection"]

DIS_COLOR = {"SCZ": "#7570b3", "MDD": "#1b9e77", "": "0.25"}
PENDING_GREY = "0.62"
STRATUM_MARKER = {"EUR": "o", "full": "D", "": "o"}   # EUR filled, pooled open


def _style():
    mpl.rcParams.update({
        "font.size": BASE, "axes.titlesize": MID, "axes.labelsize": SMALL,
        "legend.fontsize": SMALL, "xtick.labelsize": SMALL,
        "ytick.labelsize": MID, "axes.titlelocation": "left",
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })


def _fmt_p(p: float) -> str:
    return f"p={p:.1e}" if p < 1e-3 else f"p={p:.3f}".rstrip("0")


def forest(ax, sub: pd.DataFrame, dodges: dict, annotate: str | None,
           xlabel: str, zero_line: bool = True):
    """One forest panel. `dodges` maps (disorder, stratum) -> y offset."""
    for _, r in sub.iterrows():
        y = PHENOS.index(r.phenotype)
        y = (len(PHENOS) - 1 - y) + dodges[(r.disorder, r.stratum)]
        pending = r.status == "pending"
        color = PENDING_GREY if pending else DIS_COLOR[r.disorder]
        marker = "s" if pending else STRATUM_MARKER[r.stratum]
        filled = (not pending) and r.stratum != "full"
        ax.errorbar(r.estimate, y, xerr=r.se, fmt="none", ecolor=color,
                    elinewidth=1.2, capsize=2, capthick=1.0,
                    ls="--" if pending else "-", alpha=0.55 if pending else 1)
        ax.plot(r.estimate, y, marker=marker, ms=4.5, mew=1.1,
                mfc=color if filled else "white", mec=color)
        # dagger for rg rows LDSC itself flags as underpowered
        if (not pending and "underpowered" in sub.columns
                and str(r.get("underpowered")) == "yes"):
            ax.text(r.estimate + r.se + 0.02, y, "†", fontsize=SMALL,
                    va="center", color="0.35")
        if annotate and not pending:
            p_show = (min(1.0, r.p * r.m_eff) if annotate == "meff"
                      else r.p)
            if p_show < 0.05 or (annotate == "meff" and r.p < 0.05):
                lbl = _fmt_p(p_show)
                if annotate == "meff":
                    lbl = lbl.replace("p=", "p̃=")
                sgn = 1 if r.estimate >= 0 else -1   # away from zero
                ax.text(r.estimate + sgn * r.se * 1.45, y, lbl,
                        fontsize=SMALL - 1, ha="left" if sgn > 0 else "right",
                        va="center", color="0.15")
    n = len(PHENOS)
    ax.set_ylim(-0.55, n - 0.45)
    for i in range(n):                       # alternating phenotype bands
        if i % 2 == 1:
            ax.axhspan(i - 0.5, i + 0.5, color="0.955", zorder=0)
    if zero_line:
        ax.axvline(0, color="0.4", lw=0.8, zorder=1)
    ax.set_xlabel(xlabel, fontsize=SMALL, labelpad=2)
    ax.set_yticks(range(n))
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.margins(x=0.12)


def draw() -> Path:
    df = pd.read_csv(HERE / "results_v3_summary.csv").fillna({"disorder": "",
                                                              "stratum": ""})
    _style()
    fig = plt.figure(figsize=(13.333, 7.5))
    n_pend = df[df.status == "pending"].phenotype.nunique()
    fig.text(0.008, 0.955, "Experiment A — genetics of the regional-subset "
             "phenotypes vs the two anchors", fontsize=BASE + 1.5)
    fig.text(0.008, 0.915, f"grey squares = dummy placeholders for the "
             f"{n_pend} phenotypes awaiting the CSD3 run "
             "(re-run make_results_v3.py to fill)", fontsize=SMALL,
             color="0.35")

    rects = {"h2":   [0.115, 0.115, 0.14, 0.71],
             "prs":  [0.315, 0.115, 0.27, 0.71],
             "magma": [0.645, 0.115, 0.155, 0.71],
             "rg":   [0.855, 0.115, 0.135, 0.71]}
    axes = {}

    # --- h2 ---
    ax = axes["h2"] = fig.add_axes(rects["h2"])
    forest(ax, df[df.panel == "h2"], {("", ""): 0.0}, annotate=None,
           xlabel="SNP h²", zero_line=False)
    ax.set_xlim(0, 0.72)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticklabels(YLAB[::-1], fontsize=MID)
    ax.set_title("heritability\n(Zaitlen two-GRM, n = 8,082)", pad=6)

    dodge4 = {("SCZ", "EUR"): 0.28, ("SCZ", "full"): 0.10,
              ("MDD", "EUR"): -0.10, ("MDD", "full"): -0.28}
    dodge2 = {("SCZ", "EUR"): 0.17, ("MDD", "EUR"): -0.17}

    # --- PRS ---
    ax = axes["prs"] = fig.add_axes(rects["prs"])
    forest(ax, df[df.panel == "prs"], dodge4, annotate="meff",
           xlabel="β per SD of score (± SE)")
    ax.set_yticklabels([])
    ax.set_title("PRS association, population, C+T p < 0.5\n"
                 "(p̃ = m_eff-corrected across 8 thresholds)", pad=6)

    # --- MAGMA ---
    ax = axes["magma"] = fig.add_axes(rects["magma"])
    forest(ax, df[df.panel == "magma"], dodge2, annotate="raw",
           xlabel="gene-set β (± SE)")
    ax.set_yticklabels([])
    ax.set_title("MAGMA locus-pool enrichment\n(unsigned test)", pad=6)

    # --- rg ---
    ax = axes["rg"] = fig.add_axes(rects["rg"])
    forest(ax, df[df.panel == "rg"], dodge2, annotate="raw",
           xlabel="LDSC rg (± SE)")
    ax.set_yticklabels([])
    ax.set_title("genetic correlation\n(† h² z < 4: uninformative)", pad=6)

    # legend (panel b's palette applies to b-d)
    handles = [
        mpl.lines.Line2D([], [], color=DIS_COLOR["SCZ"], marker="o", ls="",
                         mfc=DIS_COLOR["SCZ"], label="SCZ"),
        mpl.lines.Line2D([], [], color=DIS_COLOR["MDD"], marker="o", ls="",
                         mfc=DIS_COLOR["MDD"], label="MDD"),
        mpl.lines.Line2D([], [], color="0.25", marker="o", ls="",
                         mfc="0.25", label="EUR stratum"),
        mpl.lines.Line2D([], [], color="0.25", marker="D", ls="",
                         mfc="white", label="pooled"),
        mpl.lines.Line2D([], [], color=PENDING_GREY, marker="s", ls="--",
                         mfc="white", label="pending (dummy)"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.99),
               ncol=5, frameon=False, fontsize=SMALL, handletextpad=0.4,
               columnspacing=1.0)

    srcs = sorted({"/".join(Path(s).parts[-2:]) for s in df.source.unique()})
    fig.text(0.008, 0.012, "sources: " + "  ·  ".join(srcs),
             fontsize=SMALL - 1.5, color="0.45")

    out = HERE / "slide_results_v3.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(out)
    return out


if __name__ == "__main__":
    draw()
    sys.exit(0)
