"""Slide: the four-method PRS grid (C+T, PRS-CS, SBayesR, SBayesRC).

Usage (repo root):  python hpc_v3/slide_prs_methods.py

Reads ONLY hpc_v2/work/results_v2/prs_final/table_main.tsv -- the canonical
grid built by collect_final.py (commit 9c64dc6: 40 cells complete, age
covariate restored, discovery GWAS matched to target arm).  Layout: GWAS
traits as facet rows, phenotypes as columns, a reading-guide column on the
right carrying the caveats a reader seeing only this slide needs.

Cell selection (the honest subset, per README 14.8):
  * matched cells only -- EUR-only discovery GWAS are never read against the
    pooled target (that mismatch produced the ASD false positive);
  * pooled cells are shown WITHIN-ANCESTRY STANDARDISED (score='zanc') --
    raw pooled Bayesian betas inflate 2-3x with their SEs because EUR-panel
    re-weighting makes score variance ancestry-dependent.

Superseded versions of this slide (pre-ff63dd8) compared methods whose
discovery GWAS differed; nothing from them should be quoted.

Output: hpc_v3/slide_prs_methods.png  (13.333 x 7.5 in, dpi 200)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
TABLE = REPO / "hpc_v2/work/results_v2/prs_final/table_main.tsv"

BASE, MID, SMALL = 10.5, 9.0, 7.5

PHENOS = ["global_slope", "baseline_thickness"]
PHENO_TITLE = {"global_slope": "global slope (thinning rate)",
               "baseline_thickness": "baseline thickness (control)"}

#: facet rows, top to bottom: hypothesis traits, then controls.  The height
#: field gives two-arm bands proportionally more vertical space (8 markers vs
#: 4) -- the ggplot space='free_y' behaviour.
BANDS = [
    ("SCZ", "schizophrenia", ("SCZ_eur", "SCZ_pooled"), 1.8),
    ("MDD", "depression", ("MDD_eur", "MDD_pooled"), 1.8),
    ("ASD", "autism (control)", ("ASD",), 1.0),
    ("ALZ", "Alzheimer's — Wightman", ("ALZ",), 1.0),
    ("ALZ_noAPOE", "  ↳ APOE excluded", ("ALZ_noAPOE",), 1.0),
    ("ALZ_IGAP", "Alzheimer's — Kunkle, no proxy", ("ALZ_IGAP",), 1.0),
    ("ALZ_IGAP_noAPOE", "  ↳ APOE excluded", ("ALZ_IGAP_noAPOE",), 1.0),
    ("EA", "education (Okbay, + control)", ("EA",), 1.0),
]

#: band geometry, top to bottom (y decreasing)
_TOTAL_H = sum(b[3] for b in BANDS)
BAND_GEOM = {}
_y = _TOTAL_H
for _key, _lab, _tas, _h in BANDS:
    BAND_GEOM[_key] = dict(top=_y, bottom=_y - _h, center=_y - _h / 2, h=_h)
    _y -= _h

# Okabe-Ito, one hue per method (grey = the non-Bayesian baseline)
METHOD_ORDER = ["CT", "PRSCS", "SBayesR", "SBayesRC"]
METHOD_LABEL = {"CT": "C+T", "PRSCS": "PRS-CS", "SBayesR": "SBayesR",
                "SBayesRC": "SBayesRC"}
METHOD_COLOR = {"CT": "0.30", "PRSCS": "#56B4E9", "SBayesR": "#D55E00",
                "SBayesRC": "#CC79A7"}
DODGE = {"CT": 0.30, "PRSCS": 0.10, "SBayesR": -0.10, "SBayesRC": -0.30}
ARM_MARKER = {"EUR": "o", "pooled": "D"}
ARM_DODGE = {"EUR": 0.045, "pooled": -0.045}


def _fmt_p(p: float) -> str:
    return f"{p:.0e}".replace("e-0", "e-") if p < 1e-3 else f"{p:.3f}"[1:]


def load() -> pd.DataFrame:
    t = pd.read_csv(TABLE, sep="\t")
    t = t[t.phenotype.isin(PHENOS) & (t.matched == "yes")]
    # pooled arm: standardized score only; EUR arm: raw
    keep = ((t.target_stratum == "EUR") & (t.score == "raw")) | \
           ((t.target_stratum == "full") & (t.score == "zanc"))
    t = t[keep].copy()
    t["arm"] = t.target_stratum.map({"EUR": "EUR", "full": "pooled"})
    band_of = {ta: name for name, _, tas, _h in BANDS for ta in tas}
    t["band"] = t.trait_arm.map(band_of)
    assert t.band.notna().all(), sorted(t[t.band.isna()].trait_arm.unique())
    # one row per cell (C+T rows are already the selected best threshold)
    dup = t.duplicated(["band", "arm", "method", "phenotype"], keep=False)
    assert not dup.any(), t[dup][["trait_arm", "method", "threshold"]]
    return t


def panel(fig, rect, d, title, show_ylab):
    ax = fig.add_axes(rect)
    for _, r in d.iterrows():
        g = BAND_GEOM[r.band]
        # dodge in units of THIS band's height, so markers fill it evenly
        y = g["center"] + (DODGE[r.method] + ARM_DODGE[r.arm]) * g["h"]
        c = METHOD_COLOR[r.method]
        sig = r.p_adj < 0.05
        ax.errorbar(r.beta, y, xerr=1.96 * r.se, fmt="none", ecolor=c,
                    elinewidth=1.1, capsize=1.5, zorder=2)
        ax.plot(r.beta, y, ARM_MARKER[r.arm], ms=4.0,
                mfc=c if sig else "white", mec=c, mew=1.0, zorder=3)
        if sig:
            sgn = 1 if r.beta >= 0 else -1
            ax.text(r.beta + sgn * (1.96 * r.se + 0.004), y, _fmt_p(r.p_adj),
                    fontsize=SMALL - 2, ha="left" if sgn > 0 else "right",
                    va="center", color="0.15", zorder=4)
    ax.set_ylim(0, _TOTAL_H)
    for i, (key, *_rest) in enumerate(BANDS):
        if i % 2 == 1:
            g = BAND_GEOM[key]
            ax.axhspan(g["bottom"], g["top"], color="0.955", zorder=0)
    ax.axvline(0, color="0.4", lw=0.8, zorder=1)
    ax.set_xlim(-0.105, 0.105)
    ax.set_yticks([BAND_GEOM[b[0]]["center"] for b in BANDS])
    if show_ylab:
        ax.set_yticklabels([b[1] for b in BANDS], fontsize=SMALL + 0.5)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("β per SD of score (95% CI)", fontsize=SMALL, labelpad=2)
    ax.set_title(title, fontsize=MID, loc="left")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=SMALL - 1)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return ax


NOTES = [
    ("Model & samples", [
        "phenotype ~ score + age + sex + PC1–10 + (1 | family), score",
        "standardised.  EUR arm n = 4,116; pooled arm n = 8,082 (all",
        "ancestries).  β per SD of score; 95% CI.",
    ]),
    ("Discovery GWAS matched to arm", [
        "Pooled arm ← multi-ancestry releases (SCZ: PGC3 primary",
        "EUR+EAS+AFR+LAT; MDD: MDD2025 trans-ancestry).  EUR arm ←",
        "European-only releases.  ASD, both ALZ releases and EA exist",
        "only as European GWAS, so they appear in the EUR arm ONLY:",
        "read against the pooled target they are confounded (score and",
        "phenotype both track ancestry) — doing so had made ASD appear",
        "significant (+0.036, p̃=.015); it is null in all 4 methods here.",
    ]),
    ("Pooled cells are standardised", [
        "Pooled points use the within-ancestry-standardised score (z-",
        "scored inside each genetic-ancestry stratum).  Raw pooled",
        "Bayesian βs run 2–3× larger with matching SEs — an artefact of",
        "EUR-panel re-weighting making score variance ancestry-",
        "dependent — so raw pooled magnitudes are never quoted.",
    ]),
    ("LD references", [
        "PRS-CS / SBayesR / SBayesRC all use UK Biobank EUROPEAN LD",
        "(no multi-ancestry panel is distributed), so pooled Bayesian",
        "cells stay LD-mismatched on the discovery side.  C+T clumps on",
        "the target genotypes themselves and is the exception.",
    ]),
    ("Multiplicity & controls", [
        "C+T: best of 8 p-thresholds, p̃ threshold-adjusted.  Bayesian",
        "methods: one score, raw p.  APOE exclusion: chr19 windows",
        "removed on BOTH builds (GRCh37 44.4–46.5 ∪ GRCh38 43.5–46.5",
        "Mb), verified 0 residual genome-wide-significant SNPs.  Kunkle",
        "= clinically diagnosed cases only (no UKB by-proxy), but Neff",
        "57.7k vs Wightman's 763k — its nulls are suggestive, not",
        "decisive.  EA is a positive control for SES confounding: it",
        "runs OPPOSITE in sign (+) to the disorders (−).",
    ]),
    ("Verdict (cross-method agreement)", [
        "SCZ: 3/4 methods, BOTH arms (PRS-CS same sign, n.s.).",
        "MDD: pooled only — 2/4 as shown (standardised); 3/4 on the",
        "raw pooled score.  EUR arm n.s. at n = 4,116.",
        "ALZ without APOE: C+T alone; null under joint modelling",
        "(SBayesR/RC) and null in Kunkle under all 4 — diffuse LD",
        "accumulation, not polygenic AD signal.  ASD: 0/4.",
    ]),
]


def draw() -> Path:
    d = load()
    mpl.rcParams.update({
        "font.size": BASE, "axes.titlesize": MID,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(13.333, 7.5))
    fig.text(0.008, 0.955, "PRS → cortical development, four methods on one "
             "grid — SCZ is the only association robust across methods and "
             "arms", fontsize=BASE + 1.5)
    fig.text(0.008, 0.917, "ancestry-matched cells from "
             "results_v2/prs_final/table_main.tsv (commit 9c64dc6, "
             "age-adjusted) · filled = p̃ < .05 · ○ EUR arm · ◇ pooled, "
             "within-ancestry standardised", fontsize=SMALL, color="0.35")

    panel(fig, [0.145, 0.085, 0.26, 0.77], d[d.phenotype == "global_slope"],
          PHENO_TITLE["global_slope"], show_ylab=True)
    panel(fig, [0.435, 0.085, 0.26, 0.77],
          d[d.phenotype == "baseline_thickness"],
          PHENO_TITLE["baseline_thickness"], show_ylab=False)

    handles = [mpl.lines.Line2D([], [], color=METHOD_COLOR[m], marker="o",
                                ls="", mfc=METHOD_COLOR[m],
                                label=METHOD_LABEL[m]) for m in METHOD_ORDER]
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(0.995, 0.995), ncol=4, frameon=False,
               fontsize=SMALL, handletextpad=0.3, columnspacing=0.9)

    # reading-guide column
    x0, y = 0.725, 0.875
    for head, lines in NOTES:
        fig.text(x0, y, head, fontsize=SMALL - 0.5, fontweight="bold",
                 color="0.20", va="top")
        y -= 0.023
        for ln in lines:
            fig.text(x0, y, ln, fontsize=SMALL - 1.5, color="0.30", va="top")
            y -= 0.0178
        y -= 0.009

    out = HERE / "slide_prs_methods.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(out)
    return out


if __name__ == "__main__":
    draw()
    sys.exit(0)
