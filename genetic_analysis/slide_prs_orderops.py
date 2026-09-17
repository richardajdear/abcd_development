"""Slide: does the order of operations change the PRS result?

Usage (repo root):  python genetic_analysis/slide_prs_orderops.py

Two constructions of the same whole-cortex thinning rate, on the same 8,596
children:
  * mean of per-region slope BLUPs  -- the pipeline's registered primary
    `global_slope` (68 DK regions / 358 HCP-MMP parcels, each BLUP shrunk by
    its OWN variance ratio, so the average is over-shrunk);
  * single LMM on the per-scan whole-cortex mean -- mean first, fit once.

Reads ONLY committed tables:
  work/results_70tab{,_hcp}/prs_final_1lmm/table_order_of_operations.tsv
  orderops/table_construction_stats.tsv   (transcribed from the run log)

SCOPE -- this run covered the 16 matched SCZ/MDD cells per atlas ONLY.  There
are no single-LMM results for the controls (ASD, ALZ x4, EA) or for
baseline_thickness, so the specificity verdict still rests entirely on the
mean-of-BLUPs construction; see docs/figures/slide_prs_methods.png.

Output: docs/figures/slide_prs_orderops.png  (13.333 x 7.5 in, dpi 200)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT = REPO / "docs/figures/slide_prs_orderops.png"

ATLASES = [("dk", "DK (68 regions)",
            HERE / "work/results_70tab/prs_final_1lmm"),
           ("hcp", "HCP-MMP (358 parcels)",
            HERE / "work/results_70tab_hcp/prs_final_1lmm")]

#: rows, top to bottom -- strongest cell first
ROWS = [("SCZ_pooled", "SCZ · pooled arm\nn = 8,596"),
        ("SCZ_eur", "SCZ · EUR arm\nn = 4,308"),
        ("MDD_pooled", "MDD · pooled arm\nn = 8,596"),
        ("MDD_eur", "MDD · EUR arm\nn = 4,308")]

METHOD_ORDER = ["CT", "PRSCS", "SBayesR", "SBayesRC"]
METHOD_LABEL = {"CT": "C+T", "PRSCS": "PRS-CS", "SBayesR": "SBayesR",
                "SBayesRC": "SBayesRC"}
METHOD_COLOR = {"CT": "0.30", "PRSCS": "#56B4E9", "SBayesR": "#D55E00",
                "SBayesRC": "#CC79A7"}
#: offsets within a row band (band height 1.0)
MDODGE = {"CT": 0.33, "PRSCS": 0.11, "SBayesR": -0.11, "SBayesRC": -0.33}
CDODGE = {"perregion": 0.055, "1lmm": -0.055}
#: marker encodes the ARM, matching slide_prs_methods.py -- the construction
#: is carried by alpha and the arrow, not by the shape.
ARM_MARKER = {"EUR": "o", "pooled": "D"}
ARM_OF = {"SCZ_eur": "EUR", "SCZ_pooled": "pooled",
          "MDD_eur": "EUR", "MDD_pooled": "pooled"}

BASE, MID, SMALL = 10.5, 9.0, 7.5
XLIM = (-0.072, 0.012)


def _fmt_p(p: float) -> str:
    return f"{p:.4f}"[1:] if p < 0.001 else f"{p:.3f}"[1:]


def load(d: Path) -> pd.DataFrame:
    t = pd.read_csv(d / "table_order_of_operations.tsv", sep="\t")
    long = []
    for _, r in t.iterrows():
        for c in ("perregion", "1lmm"):
            long.append(dict(trait_arm=r.trait_arm, method=r.method,
                             construction=c, beta=r[f"beta_{c}"],
                             se=r[f"se_{c}"], p_adj=r[f"p_adj_{c}"], n=r.n))
    out = pd.DataFrame(long)
    assert len(out) == 32, len(out)
    return out


def panel(fig, rect, d, title, show_ylab):
    ax = fig.add_axes(rect)
    nrow = len(ROWS)
    for _, r in d.iterrows():
        ri = [k for k, _ in ROWS].index(r.trait_arm)
        y = (nrow - 1 - ri) + MDODGE[r.method] + CDODGE[r.construction]
        c = METHOD_COLOR[r.method]
        sig = r.p_adj < 0.05
        one = r.construction == "1lmm"
        ax.errorbar(r.beta, y, xerr=1.96 * r.se, fmt="none", ecolor=c,
                    elinewidth=1.0 if one else 0.8, capsize=1.2,
                    alpha=1.0 if one else 0.40, zorder=2)
        ax.plot(r.beta, y, ARM_MARKER[ARM_OF[r.trait_arm]],
                ms=4.4 if one else 3.4,
                mfc=c if sig else "white", mec=c, mew=1.0,
                alpha=1.0 if one else 0.50, zorder=3)
        if sig:
            ax.text(r.beta - 1.96 * r.se - 0.0025, y, _fmt_p(r.p_adj),
                    fontsize=SMALL - 2.5, ha="right", va="center",
                    color="0.15", zorder=4)
    # the shift, drawn as a connector per method
    piv = d.pivot_table(index=["trait_arm", "method"], columns="construction",
                        values="beta")
    for (ta, m), row in piv.iterrows():
        ri = [k for k, _ in ROWS].index(ta)
        y0 = (nrow - 1 - ri) + MDODGE[m] + CDODGE["perregion"]
        y1 = (nrow - 1 - ri) + MDODGE[m] + CDODGE["1lmm"]
        ax.annotate("", xy=(row["1lmm"], y1), xytext=(row["perregion"], y0),
                    arrowprops=dict(arrowstyle="-|>,head_width=0.09,"
                                    "head_length=0.22",
                                    color=METHOD_COLOR[m], lw=0.7,
                                    alpha=0.55, shrinkA=2.5, shrinkB=2.5),
                    zorder=1.5)
    for i in range(nrow):
        if i % 2 == 1:
            ax.axhspan(i - 0.5, i + 0.5, color="0.955", zorder=0)
    ax.axvline(0, color="0.4", lw=0.8, zorder=1)
    ax.set_xlim(*XLIM)
    ax.set_ylim(-0.5, nrow - 0.5)
    ax.set_xticks([-0.06, -0.04, -0.02, 0.0])
    ax.set_yticks(range(nrow))
    ax.set_yticklabels([lab for _, lab in ROWS][::-1] if show_ylab else [],
                       fontsize=SMALL)
    ax.set_xlabel("β per SD of score (95% CI)", fontsize=SMALL - 1,
                  labelpad=2)
    ax.set_title(title, fontsize=MID, pad=4)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=SMALL - 2)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return ax


def counts(d: pd.DataFrame) -> dict:
    g = d.groupby(["trait_arm", "construction"]).p_adj.apply(
        lambda s: int((s < 0.05).sum()))
    return g.to_dict()


def draw() -> Path:
    mpl.rcParams.update({"font.size": BASE, "figure.facecolor": "white",
                         "savefig.facecolor": "white"})
    st = pd.read_csv(HERE / "orderops/table_construction_stats.tsv", sep="\t")
    S = {(r.metric, r.atlas): r.value for _, r in st.iterrows()}
    data = {k: load(p) for k, _, p in ATLASES}
    C = {k: counts(v) for k, v in data.items()}

    fig = plt.figure(figsize=(13.333, 7.5))
    fig.text(0.006, 0.955, "Order of operations: mean of per-region slope "
             "BLUPs vs one LMM on the per-scan mean — SCZ pooled invariant, "
             "MDD strengthens throughout", fontsize=BASE + 0.5)
    fig.text(0.006, 0.921, "same 8,596 children, scores and model · ○ EUR "
             "· ◇ pooled · faded = mean of per-region BLUPs (primary), solid = "
             "single LMM, arrow = shift · filled = p̃ < .05",
             fontsize=SMALL - 0.5, color="0.35")

    x0, w, gap, y0, h = 0.098, 0.245, 0.035, 0.115, 0.735
    for k, (key, label, path) in enumerate(ATLASES):
        panel(fig, [x0 + k * (w + gap), y0, w, h], data[key], label,
              show_ylab=(k == 0))

    handles = [mpl.lines.Line2D([], [], color=METHOD_COLOR[m], marker="s",
                                ls="", mfc=METHOD_COLOR[m],
                                label=METHOD_LABEL[m]) for m in METHOD_ORDER]
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(0.995, 0.940), ncol=4, frameon=False,
               fontsize=SMALL, handletextpad=0.3, columnspacing=0.9)

    def cnt(atlas, ta):
        return f"{C[atlas][(ta, 'perregion')]}/4 → {C[atlas][(ta, '1lmm')]}/4"

    notes = [
        ("Two constructions, one measure", [
            f"r(mean of BLUPs, single-LMM slope) = "
            f"{S[('r_slope_constructions','dk')]:.3f} DK / "
            f"{S[('r_slope_constructions','hcp')]:.3f} HCP-MMP; intercepts",
            f"agree at {S[('r_intercept_constructions','dk')]:.4f}.  Each "
            "regional BLUP is shrunk by its OWN",
            "variance ratio, so averaging 68 (358) of them over-shrinks child",
            f"by child: the mean-of-BLUPs slope keeps "
            f"{S[('sd_slope_perregion_mm_per_yr','dk')] / S[('sd_slope_1lmm_mm_per_yr','dk')]:.0%}"
            f" / "
            f"{S[('sd_slope_perregion_mm_per_yr','hcp')] / S[('sd_slope_1lmm_mm_per_yr','hcp')]:.0%}"
            " of the",
            f"single-LMM SD, agreement rises with scans "
            f"({S[('r_slope_visits2','dk')]:.2f}/{S[('r_slope_visits3','dk')]:.2f}/"
            f"{S[('r_slope_visits4','dk')]:.2f} for 2/3/4),",
            f"and mean-first is more atlas-invariant "
            f"({S[('cross_atlas_r_slope','perregion')]:.3f} → "
            f"{S[('cross_atlas_r_slope','1lmm')]:.3f}).",
        ]),
        ("β moves, precision does not", [
            "SEs identical to the third decimal (ratio 0.996) — set by n and",
            "the score–phenotype correlation, not the construction.  What",
            "moves is β (mean −0.004, mean-first more negative, max |Δ| .020):",
            "recovered signal, not a precision gain.",
        ]),
        ("Cells significant (p̃ < .05), mean-of-BLUPs → single LMM", [
            f"SCZ pooled   {cnt('dk','SCZ_pooled')} (DK)   "
            f"{cnt('hcp','SCZ_pooled')} (HCP-MMP)",
            f"SCZ EUR      {cnt('dk','SCZ_eur')} (DK)   "
            f"{cnt('hcp','SCZ_eur')} (HCP-MMP)",
            f"MDD pooled   {cnt('dk','MDD_pooled')} (DK)   "
            f"{cnt('hcp','MDD_pooled')} (HCP-MMP)",
            f"MDD EUR      {cnt('dk','MDD_eur')} (DK)   "
            f"{cnt('hcp','MDD_eur')} (HCP-MMP)",
        ]),
        ("Reading", [
            "SCZ pooled is invariant (16/16 cells across method × atlas ×",
            "construction, βs within .003).  MDD strengthens in all 16 cells",
            "and MDD-EUR turns nominally significant under 2–3 methods — a",
            "one-directional shift ≈ 1.5–2 SE of the difference, consistent",
            "with a weak, spatially spread signal being less attenuated in the",
            "less-shrunk phenotype (MDD loads on slope PC3 as much as on the",
            "global mean; SCZ is carried almost entirely by PC1).  The SCZ-EUR",
            "'weakening' is ≈ 0.5 SE and absent from SCZ pooled: noise.",
        ]),
        ("Scope and the SCZ-EUR question since resolved", [
            "This run covered the 16 matched SCZ/MDD cells per atlas with the",
            "PGC3 SCZ GWAS only — no controls, no baseline thickness, and not",
            "yet the 2025 SCZ GWAS, under which SCZ-EUR is 4/4 on both atlases",
            "(slide_prs_methods).  Mean-of-BLUPs remains the registered",
            "primary; the single LMM is a documented sensitivity phenotype",
            "until the full grid is re-run on it.",
        ]),
    ]

    xn, y = 0.635, 0.905
    for head, lines in notes:
        fig.text(xn, y, head, fontsize=SMALL - 1, fontweight="bold",
                 color="0.20", va="top")
        y -= 0.0205
        for ln in lines:
            fig.text(xn, y, ln, fontsize=SMALL - 2.2, color="0.30", va="top")
            y -= 0.0158
        y -= 0.007

    fig.savefig(OUT, dpi=200)
    plt.close(fig)
    print(OUT)
    return OUT


if __name__ == "__main__":
    draw()
    sys.exit(0)
