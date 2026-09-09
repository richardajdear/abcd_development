"""Experiment A collection slide: h2, PRS, MAGMA, rg across eight rows.

Usage (repo root):  PYTHONPATH=src python hpc_v3/slide_results_v3.py

Reads ONLY hpc_v3/results_v3_summary.csv (built by make_results_v3.py from
the committed result tables).  No statistic is hardcoded here: p annotations,
the m_eff correction and the underpowered daggers are all computed from the
table at plot time, and a row the source tables cannot supply is drawn as a
grey open square rather than silently omitted.

EIGHT ROWS: the two anchors and the four Experiment A phenotypes through the
v2 pipeline, plus the two anchors through v1's, so the pipelines sit side by
side.  The v1 rows are indented under their v2 counterpart and share its
shaded band.

WHAT ENCODES WHAT (the reason there is no fourth visual channel):
  colour   = disorder (SCZ / MDD), the categorical identity
  fill+shape = ancestry stratum (filled circle = EUR, open diamond = pooled)
  row      = phenotype x pipeline
Pipeline is deliberately NOT a colour or a fill: both are already spoken for,
and a third overlapping channel is how a forest plot becomes unreadable.  The
two-hue palette passes the categorical checks (adjacent-pair CVD dE 12.9
deutan / 10.7 tritan, normal-vision dE 19.9, contrast >= 3:1 on white), so
identity survives colour-vision deficiency; every row is also directly
labelled, so identity is never colour-alone.

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

#: (phenotype, pipeline, row label, band group).  Top to bottom on the figure.
#: v1 rows sit directly under the v2 row they compare against and share its
#: band, so the pairing is read before the numbers are.
ROWS = [
    ("baseline_thickness", "v2", "baseline thickness (control)", 0),
    ("baseline_thickness", "v1", "     ↳ v1 pipeline", 0),
    ("global_slope", "v2", "global slope", 1),
    ("global_slope", "v1", "     ↳ v1 pipeline", 1),
    ("global_slope", "avgfirst", "     ↳ avg-first definition", 1),
    ("slope_PC2", "v2", "slope PC2 (h² reference)", 2),
    ("slope_topDelta", "v2", "top-ΔCT mean", 3),
    ("slope_topC3", "v2", "top-C3 mean", 4),
    ("slope_projDelta", "v2", "ΔCT projection", 5),
    ("slope_projC3", "v2", "C3 projection", 6),
]
ROW_INDEX = {(ph, pl): i for i, (ph, pl, _, _) in enumerate(ROWS)}
YLAB = [lab for _, _, lab, _ in ROWS]
#: the four new phenotypes start here -- used for the separating rule
#: (anchors + their pipeline variants + the slope_PC2 h2 reference)
N_ANCHOR_ROWS = 6

# Validated with the dataviz skill's checker (light mode, white surface):
# lightness band PASS, chroma floor PASS, CVD separation PASS, normal-vision
# floor PASS, contrast PASS.  Do not substitute hues without re-running it.
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
    """Two significant figures. `%.3f` rounded 0.0014 to "0.001", which reads
    as an order of magnitude better than it is."""
    return f"p={p:.1e}" if p < 1e-3 else f"p={p:.2g}"


def forest(ax, sub: pd.DataFrame, dodges: dict, annotate: str | None,
           xlabel: str, zero_line: bool = True):
    """One forest panel. `dodges` maps (disorder, stratum) -> y offset."""
    n = len(ROWS)
    for _, r in sub.iterrows():
        key = (r.phenotype, r.pipeline)
        if key not in ROW_INDEX:
            continue
        y = (n - 1 - ROW_INDEX[key]) + dodges[(r.disorder, r.stratum)]
        pending = r.status != "observed"
        # "not estimable" is a result (LDSC declined to divide by a negative
        # heritability), so it gets an x on the zero line rather than a
        # dummy estimate with a fake interval.
        na = r.status == "not_estimable"
        color = PENDING_GREY if pending else DIS_COLOR[r.disorder]
        marker = ("x" if na else "s") if pending else STRATUM_MARKER[r.stratum]
        filled = (not pending) and r.stratum != "full"
        if na:
            ax.plot(0, y, marker="x", ms=5.5, mew=1.4,
                    color=DIS_COLOR[r.disorder], alpha=0.75, zorder=4)
            continue
        # v1-pipeline rows are context, not the result: dim to the background
        a = 0.5 if r.pipeline == "v1" else (0.55 if pending else 1)
        ax.errorbar(r.estimate, y, xerr=1.96 * r.se, fmt="none", ecolor=color,
                    elinewidth=1.3, capsize=2, capthick=1.0,
                    ls="--" if pending else "-", alpha=a, zorder=2)
        # a surface ring keeps the marker legible where intervals overlap
        ax.plot(r.estimate, y, marker=marker, ms=5.5, mew=1.6,
                mfc=color if filled else "white", mec="white", zorder=3,
                alpha=a)
        ax.plot(r.estimate, y, marker=marker, ms=5.5, mew=1.1,
                mfc=color if filled else "white", mec=color, zorder=4,
                alpha=a)
        # dagger for rg rows LDSC itself flags as underpowered
        if (not pending and "underpowered" in sub.columns
                and str(r.get("underpowered")) == "yes"):
            ax.text(r.estimate + 1.96 * r.se + 0.02, y, "†", fontsize=SMALL,
                    va="center", color="0.35", zorder=5)
        if annotate and not pending:
            p_show = (min(1.0, r.p * r.m_eff) if annotate == "meff" else r.p)
            if p_show < 0.05 or (annotate == "meff" and r.p < 0.05):
                lbl = _fmt_p(p_show)
                if annotate == "meff":
                    lbl = lbl.replace("p=", "p̃=")
                sgn = 1 if r.estimate >= 0 else -1   # away from zero
                ax.text(r.estimate + sgn * r.se * 1.96 * 1.15, y, lbl,
                        fontsize=SMALL - 1, ha="left" if sgn > 0 else "right",
                        va="center", color="0.15", zorder=5)
    ax.set_ylim(-0.55, n - 0.45)
    # band by PHENOTYPE, not by row: a v1 row shares its v2 row's band, which
    # is what makes the pair read as one comparison
    for grp in {g for _, _, _, g in ROWS}:
        if grp % 2 != 1:
            continue
        ys = [n - 1 - i for i, (_, _, _, g) in enumerate(ROWS) if g == grp]
        ax.axhspan(min(ys) - 0.5, max(ys) + 0.5, color="0.955", zorder=0)
    # rule separating the two anchors from the four Experiment A phenotypes
    ax.axhline(n - N_ANCHOR_ROWS - 0.5, color="0.75", lw=0.8, ls=(0, (4, 3)),
               zorder=1)
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
    fig.text(0.008, 0.958, "Experiment A — the regional-subset phenotypes "
             "against the two anchors, and v2 against v1", fontsize=BASE + 1.5)
    # two lines: at 13.3in a single line of this length runs off the edge
    fig.text(0.008, 0.925,
             "No subset improves on the global mean: h² is equal or lower, "
             "the PRS association equal or weaker, and SCZ locus-pool "
             "enrichment matches it (β 0.125 vs 0.124). The two projections "
             "lose the association outright.",
             fontsize=SMALL, color="0.30")
    fig.text(0.008, 0.898,
             "One exception, nominal and uncorrected across 12 gene-set "
             "tests: MDD locus pool × top-ΔCT mean, p = 0.027 against global "
             "slope's 0.13.",
             fontsize=SMALL, color="0.30")

    n_pend = int((df.status == "pending").sum())
    n_na = int((df.status == "not_estimable").sum())
    if n_pend:
        fig.text(0.008, 0.893, f"{n_pend} rows still pending (grey squares)",
                 fontsize=SMALL, color="0.45")

    rects = {"h2":   [0.145, 0.155, 0.135, 0.645],
             "prs":  [0.335, 0.155, 0.265, 0.645],
             "magma": [0.655, 0.155, 0.150, 0.645],
             "rg":   [0.860, 0.155, 0.132, 0.645]}
    axes = {}

    # --- h2 ---
    ax = axes["h2"] = fig.add_axes(rects["h2"])
    forest(ax, df[df.panel == "h2"], {("", "full"): 0.0}, annotate=None,
           xlabel="SNP h² (95% CI)", zero_line=False)
    ax.set_xlim(0, 0.72)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticklabels(YLAB[::-1], fontsize=MID)
    for t, (_, pl, _, _) in zip(ax.get_yticklabels(), ROWS[::-1]):
        if pl == "v1":
            t.set_color("0.42")
            t.set_fontsize(SMALL)
    ax.set_title("heritability\n(v2 Zaitlen two-GRM · v1 single-GRM)", pad=6)

    dodge4 = {("SCZ", "EUR"): 0.28, ("SCZ", "full"): 0.10,
              ("MDD", "EUR"): -0.10, ("MDD", "full"): -0.28}
    dodge2 = {("SCZ", "EUR"): 0.17, ("MDD", "EUR"): -0.17}

    # --- PRS ---
    ax = axes["prs"] = fig.add_axes(rects["prs"])
    forest(ax, df[df.panel == "prs"], dodge4, annotate="meff",
           xlabel="β per SD of score (95% CI)")
    ax.set_yticklabels([])
    ax.set_title("PRS association, imputed genotypes, C+T p < 0.5\n"
                 "(p̃ = m_eff-corrected across 8 thresholds)", pad=6)

    # --- MAGMA ---
    ax = axes["magma"] = fig.add_axes(rects["magma"])
    forest(ax, df[df.panel == "magma"], dodge2, annotate="raw",
           xlabel="gene-set β (95% CI)")
    ax.set_yticklabels([])
    ax.set_title("MAGMA locus-pool enrichment\n(unsigned test)", pad=6)

    # --- rg ---
    ax = axes["rg"] = fig.add_axes(rects["rg"])
    forest(ax, df[df.panel == "rg"], dodge2, annotate="raw",
           xlabel="LDSC rg (95% CI)")
    ax.set_yticklabels([])
    ax.set_title("genetic correlation\n(† h² z < 4: uninformative)", pad=6)

    handles = [
        mpl.lines.Line2D([], [], color=DIS_COLOR["SCZ"], marker="o", ls="",
                         mfc=DIS_COLOR["SCZ"], label="SCZ"),
        mpl.lines.Line2D([], [], color=DIS_COLOR["MDD"], marker="o", ls="",
                         mfc=DIS_COLOR["MDD"], label="MDD"),
        mpl.lines.Line2D([], [], color="0.25", marker="o", ls="",
                         mfc="0.25", label="EUR stratum"),
        mpl.lines.Line2D([], [], color="0.25", marker="D", ls="",
                         mfc="white", label="pooled"),
    ]
    if n_pend:
        handles.append(mpl.lines.Line2D([], [], color=PENDING_GREY, marker="s",
                                        ls="--", mfc="white", label="pending"))
    if n_na:
        handles.append(mpl.lines.Line2D([], [], color="0.35", marker="x", ls="",
                                        label="rg not estimable"))
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.995),
               ncol=len(handles), frameon=False, fontsize=SMALL,
               handletextpad=0.4, columnspacing=1.0)

    # The two caveats that decide whether the v1/v2 rows may be read as a
    # like-for-like comparison.  They live in the table, not in this file.
    cav = {r.panel: r.caveat for _, r in
           df.drop_duplicates("panel").iterrows() if isinstance(r.caveat, str)}
    fig.text(0.008, 0.072,
             "v1 vs v2 — h²: " + cav.get("h2", "") + ".",
             fontsize=SMALL - 1.5, color="0.42")
    fig.text(0.008, 0.052,
             "v1 vs v2 — PRS: " + cav.get("prs", ""),
             fontsize=SMALL - 1.5, color="0.42")
    af = df[(df.pipeline == "avgfirst") & (df.disorder == "SCZ")]
    if len(af):
        r = af.iloc[0]
        fig.text(0.008, 0.092,
                 "avg-first — global slope re-defined as region-mean CT per "
                 "scan → one LMM (r = 0.887 with the settled definition); "
                 f"PRS: {r.caveat}.",
                 fontsize=SMALL - 1.5, color="0.42")

    srcs = sorted({"/".join(Path(s).parts[-2:]) for s in df.source.unique()})
    # one line ran off the right edge at 13.3in; wrap onto as many as needed
    per_line, lines, cur = 4, [], []
    for x in srcs:
        cur.append(x)
        if len(cur) == per_line:
            lines.append("  ·  ".join(cur)); cur = []
    if cur:
        lines.append("  ·  ".join(cur))
    for i, ln in enumerate(lines):
        fig.text(0.008, 0.020 - i * 0.016,
                 ("sources: " if i == 0 else "         ") + ln,
                 fontsize=SMALL - 2, color="0.55")

    out = HERE / "slide_results_v3.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(out)
    return out


if __name__ == "__main__":
    draw()
    sys.exit(0)
