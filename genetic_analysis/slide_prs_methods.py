"""Slide: the four-method PRS grid in BOTH parcellations (DK and HCP-MMP).

Usage (repo root):  python genetic_analysis/slide_prs_methods.py

Reads ONLY committed grid tables, per parcellation root
(work/results_70tab = DK 68 regions; work/results_70tab_hcp = HCP-MMP 358):
    prs_final/table_main.tsv            PGC3 SCZ, MDD2025, controls (step 6)
    prs_scz2025/table_scz2025_main.tsv  the 2025 multi-ancestry SCZ GWAS
                                        re-scored (step 9), primary cells
Traits are facet ROWS; the four COLUMNS are phenotype x parcellation, so the
atlas comparison sits side by side within each phenotype.

Cell selection (the honest subset, README_HPC.md sec 14.8):
  * matched cells only -- a European-only discovery GWAS is never read
    against the pooled target (that mismatch produced the ASD false positive);
  * pooled cells are shown WITHIN-ANCESTRY STANDARDISED (score='zanc'),
    because raw pooled Bayesian betas inflate 2-3x with their SEs.

Predecessor: legacy/hpc_v3/slide_prs_methods.py, DK-only and on the
6.0-vintage tables.  Nothing from it should be quoted.

Output: docs/figures/slide_prs_methods.png  (13.333 x 7.5 in, dpi 200)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT = REPO / "docs/figures/slide_prs_methods.png"

#: (phenotype, parcellation label, results root).  Each root holds the PGC3
#: grid (prs_final/table_main.tsv) and the step-9 2025 SCZ GWAS re-scoring
#: (prs_scz2025/table_scz2025_main.tsv); load() reads both.
COLUMNS = [
    ("global_slope", "DK (68)", HERE / "work/results_70tab"),
    ("global_slope", "HCP-MMP (358)", HERE / "work/results_70tab_hcp"),
    ("baseline_thickness", "DK (68)", HERE / "work/results_70tab"),
    ("baseline_thickness", "HCP-MMP (358)", HERE / "work/results_70tab_hcp"),
]
GROUPS = [("global slope (thinning rate)", 0, 2),
          ("baseline thickness (control phenotype)", 2, 4)]

BASE, MID, SMALL = 10.5, 9.0, 7.5

#: facet rows, top to bottom; height gives the two-arm disorder bands more
#: room (8 markers vs 4) -- the ggplot space='free_y' behaviour.
BANDS = [
    ("SCZ25", "schizophrenia (2025 GWAS)", ("SCZ25_EUR", "SCZ25_META"), 1.8),
    ("SCZ", "  ↳ PGC3 2022 GWAS", ("SCZ_eur", "SCZ_pooled"), 1.8),
    ("MDD", "depression", ("MDD_eur", "MDD_pooled"), 1.8),
    ("ASD", "autism (control)", ("ASD",), 1.0),
    ("ALZ", "Alzheimer's — Wightman", ("ALZ",), 1.0),
    ("ALZ_noAPOE", "  ↳ APOE excluded", ("ALZ_noAPOE",), 1.0),
    ("ALZ_IGAP", "Alzheimer's — Kunkle", ("ALZ_IGAP",), 1.0),
    ("ALZ_IGAP_noAPOE", "  ↳ APOE excluded", ("ALZ_IGAP_noAPOE",), 1.0),
    ("EA", "education (Okbay)", ("EA",), 1.0),
]
_TOTAL_H = sum(b[3] for b in BANDS)
BAND_GEOM, _y = {}, _TOTAL_H
for _key, _lab, _tas, _h in BANDS:
    BAND_GEOM[_key] = dict(top=_y, bottom=_y - _h, center=_y - _h / 2, h=_h)
    _y -= _h

METHOD_ORDER = ["CT", "PRSCS", "SBayesR", "SBayesRC"]
METHOD_LABEL = {"CT": "C+T", "PRSCS": "PRS-CS", "SBayesR": "SBayesR",
                "SBayesRC": "SBayesRC"}
METHOD_COLOR = {"CT": "0.30", "PRSCS": "#56B4E9", "SBayesR": "#D55E00",
                "SBayesRC": "#CC79A7"}
DODGE = {"CT": 0.30, "PRSCS": 0.10, "SBayesR": -0.10, "SBayesRC": -0.30}
ARM_MARKER = {"EUR": "o", "pooled": "D"}
ARM_DODGE = {"EUR": 0.045, "pooled": -0.045}
XLIM = (-0.088, 0.088)


def _fmt_p(p: float) -> str:
    return f"{p:.0e}".replace("e-0", "e-") if p < 1e-3 else f"{p:.3f}"[1:]


def load(root: Path, phenotype: str) -> pd.DataFrame:
    g = pd.read_csv(root / "prs_final/table_main.tsv", sep="\t")
    g = g[g.matched == "yes"]
    # step 9: the primary cells only (EUR GWAS -> EUR anchor; AFR+EUR+EAS
    # meta -> pooled).  The four-method set is kept so the rows align with
    # PGC3; PRS-CSx and the ancestry-matched composite live in the notes.
    s = pd.read_csv(root / "prs_scz2025/table_scz2025_main.tsv", sep="\t")
    s = s[(s.design == "primary") & s.method.isin(METHOD_ORDER)
          & s.trait_arm.isin(["SCZ25_EUR", "SCZ25_META"])]
    t = pd.concat([g, s], ignore_index=True)
    t = t[t.phenotype == phenotype]
    keep = ((t.target_stratum == "EUR") & (t.score == "raw")) | \
           ((t.target_stratum == "full") & (t.score == "zanc"))
    t = t[keep].copy()
    t["arm"] = t.target_stratum.map({"EUR": "EUR", "full": "pooled"})
    band_of = {ta: name for name, _, tas, _h in BANDS for ta in tas}
    t["band"] = t.trait_arm.map(band_of)
    assert t.band.notna().all(), sorted(t[t.band.isna()].trait_arm.unique())
    dup = t.duplicated(["band", "arm", "method"], keep=False)
    assert not dup.any(), t[dup][["trait_arm", "method", "threshold"]]
    return t


def panel(fig, rect, d, title, show_ylab):
    ax = fig.add_axes(rect)
    for _, r in d.iterrows():
        g = BAND_GEOM[r.band]
        y = g["center"] + (DODGE[r.method] + ARM_DODGE[r.arm]) * g["h"]
        c = METHOD_COLOR[r.method]
        sig = r.p_adj < 0.05
        ax.errorbar(r.beta, y, xerr=1.96 * r.se, fmt="none", ecolor=c,
                    elinewidth=1.0, capsize=1.2, zorder=2)
        ax.plot(r.beta, y, ARM_MARKER[r.arm], ms=3.6,
                mfc=c if sig else "white", mec=c, mew=0.9, zorder=3)
        if sig:
            sgn = 1 if r.beta >= 0 else -1
            ax.text(r.beta + sgn * (1.96 * r.se + 0.003), y, _fmt_p(r.p_adj),
                    fontsize=SMALL - 2.5, ha="left" if sgn > 0 else "right",
                    va="center", color="0.15", zorder=4)
    ax.set_ylim(0, _TOTAL_H)
    for i, (key, *_rest) in enumerate(BANDS):
        if i % 2 == 1:
            g = BAND_GEOM[key]
            ax.axhspan(g["bottom"], g["top"], color="0.955", zorder=0)
    ax.axvline(0, color="0.4", lw=0.8, zorder=1)
    ax.set_xlim(*XLIM)
    ax.set_xticks([-0.05, 0.0, 0.05])
    ax.set_yticks([BAND_GEOM[b[0]]["center"] for b in BANDS])
    ax.set_yticklabels([b[1] for b in BANDS] if show_ylab else [],
                       fontsize=SMALL)
    ax.set_xlabel("β per SD (95% CI)", fontsize=SMALL - 1, labelpad=1.5)
    ax.set_title(title, fontsize=SMALL + 0.5, loc="center", pad=3,
                 color="0.25")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=SMALL - 2)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return ax


NOTES = [
    ("Model & samples", [
        "phenotype ~ score + age + sex + PC1–10 + (1 | family); score",
        "standardised.  EUR arm n = 4,308; pooled n = 8,596 (all ancestries),",
        "score z-scored within ancestry cluster.  Release 7.0 tables.",
    ]),
    ("Discovery GWAS matched to arm", [
        "EUR arm ← European-only GWAS; pooled ← multi-ancestry GWAS.  ASD,",
        "ALZ and EA exist only as European GWAS so appear EUR-only: read",
        "against the pooled target they are confounded (score and phenotype",
        "both track ancestry; that made ASD look significant, +0.036).",
        "SCZ 2025 = Nature multi-ancestry release (EUR Neff 117k → 174k,",
        "no Latino cohort); PGC3 kept as the row below for comparison.",
    ]),
    ("What the second row (2025 GWAS) settles", [
        "The EUR arm moves from borderline (PGC3: 0/4 DK, 3/4 HCP-MMP) to",
        "4/4 on both atlases, β −0.035 to −0.051, and the C+T permutation",
        "p 0.053 → 0.0085.  Nothing was selected — the discovery GWAS grew",
        "1.5×.  Pooled stays 4/4, slightly weaker (the meta dropped the",
        "Latino cohort).  PRS × ancestry-cluster LRT p ≥ .10 in every cell,",
        "so one pooled β is licensed; within clusters the meta weights",
        "carry same-sign effects into the Hispanic- and mixed/Asian-like",
        "clusters (β −.03 to −.09, n 928 / 380) but ~0 into the African-",
        "American-like cluster under the Bayesian methods (C+T −.03); AFR-",
        "and EAS-derived weights predict nothing (Neff 35k / 29k).",
    ]),
    ("Atlases and scores", [
        "DK and HCP-MMP are one whole-cortex measure formed two ways on the",
        "same children (GREML h² 0.166 vs 0.156) — a robustness check, not",
        "replication; where they differ (MDD pooled 3/4 vs 1/4) they differ",
        "at p ≈ .05 with near-identical βs.  Raw pooled Bayesian βs inflate",
        "2–3× with their SEs, so pooled cells are within-cluster z-scored.",
        "All three Bayesian methods use UKB EUROPEAN LD; C+T clumps on the",
        "target.  C+T: best of 8 thresholds, p̃ adjusted; others: one score,",
        "raw p; no correction across the grid (FDR reported separately).",
    ]),
    ("Controls", [
        "APOE removed on both builds, 0 residual GWS SNPs.  Kunkle: no UKB",
        "by-proxy but Neff 57.7k vs Wightman 763k — nulls suggestive.  EA is",
        "POSITIVE (opposite sign; cannot manufacture the disorders' βs) yet",
        "p̃ < .05 in 3/4 methods on DK, 0/4 HCP-MMP — not a clean null; the",
        "conditional-on-EA test matters.  Baseline thickness is null in every",
        "SCZ cell, both GWAS: the association is slope-specific.",
    ]),
    ("Verdict", [
        "SCZ: significant under every method, both arms, both atlases with",
        "the 2025 GWAS; the EUR estimate is now the stronger one.",
        "MDD: pooled only — 3/4 (DK) / 1/4 (HCP-MMP); EUR arm null.",
        "ALZ without APOE ≤ 1/4 either atlas, 0/4 Kunkle.  ASD 0/4.",
    ]),
]


def draw() -> Path:
    mpl.rcParams.update({
        "font.size": BASE, "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(13.333, 7.5))
    fig.text(0.006, 0.958, "PRS → cortical development: four methods × two "
             "parcellations × two SCZ GWAS — SCZ is the only association "
             "robust across all of them", fontsize=BASE + 0.5)
    fig.text(0.006, 0.925, "ancestry-matched cells from "
             "results_70tab{,_hcp}/prs_final/table_main.tsv and "
             "prs_scz2025/table_scz2025_main.tsv · filled = p̃ < .05 · "
             "○ EUR arm · ◇ pooled, within-ancestry standardised",
             fontsize=SMALL - 0.5, color="0.35")

    x0, w, gap, y0, h = 0.128, 0.127, 0.0140, 0.075, 0.755
    for k, (pheno, parc, table) in enumerate(COLUMNS):
        rect = [x0 + k * (w + gap), y0, w, h]
        panel(fig, rect, load(table, pheno), parc, show_ylab=(k == 0))
    # phenotype group headers spanning their two columns
    for label, a, b in GROUPS:
        xc = x0 + (a * (w + gap)) + ((b - a) * (w + gap) - gap) / 2
        fig.text(xc, y0 + h + 0.045, label, fontsize=MID, ha="center")
        xa = x0 + a * (w + gap)
        xb = xa + (b - a) * (w + gap) - gap
        fig.add_artist(mpl.lines.Line2D([xa, xb], [y0 + h + 0.038] * 2,
                                        color="0.75", lw=0.8,
                                        transform=fig.transFigure))

    handles = [mpl.lines.Line2D([], [], color=METHOD_COLOR[m], marker="o",
                                ls="", mfc=METHOD_COLOR[m],
                                label=METHOD_LABEL[m]) for m in METHOD_ORDER]
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(0.995, 0.955), ncol=4, frameon=False,
               fontsize=SMALL, handletextpad=0.3, columnspacing=0.9)

    xn, y = 0.695, 0.905
    for head, lines in NOTES:
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
