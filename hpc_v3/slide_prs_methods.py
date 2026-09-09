"""Slide: the three PRS methods (C+T, PRS-CS, SBayesR) side by side.

Usage (repo root):  python hpc_v3/slide_prs_methods.py

Context (hpc_v2/README_HPC.md §13.5 → §13.8): §13.5 concluded the SCZ →
global_slope association "does not survive a method without threshold
selection" from PRS-CS alone; SBayesR — a second untuned method, with a
50,000-individual LD reference against PRS-CS's 503 — reverses that reading.
This slide puts the three methods on one axis for every phenotype each has
been scored on, POOLED stratum (the only stratum the SBayesR record covers).

Sources (all committed):
  C+T      hpc_v3/prs_tables/prs_association_v3.tsv  (p<0.5 score; imputed;
           p̃ = p x m_eff, m_eff = 4.6, as on slide_results_v3)
  PRS-CS   hpc_v2/work/results_v2/prscs/prs_association_<DIS>_prscs.tsv
           (threshold "auto"; single score, no correction needed)
  SBayesR  hpc_v3/prs_tables/prs_sbayesr_readme.tsv — transcribed from
           hpc_v2/README_HPC.md §13.8 (the association tables live on CSD3
           and are not yet committed).  SEs are RECONSTRUCTED from beta and p
           via the normal quantile, se = |beta| / Phi^-1(1 - p/2) — exact for
           a Wald z test, so intervals are faithful to the reported p.

The four Experiment A phenotypes exist under C+T only — PRS-CS and SBayesR
have not been scored against them (single-score files exist on CSD3; a
one-line assoc job would fill them).

Output: hpc_v3/slide_prs_methods.png  (13.333 x 7.5 in, dpi 200)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

BASE, MID, SMALL = 10.5, 9.0, 7.5
M_EFF = 4.6

PHENOS = ["baseline_thickness", "global_slope", "slope_PC2",
          "slope_topDelta", "slope_topC3", "slope_projDelta", "slope_projC3"]
YLAB = ["baseline thickness\n(control)", "global slope",
        "slope PC2 (h² reference)", "top-ΔCT mean", "top-C3 mean",
        "ΔCT projection", "C3 projection"]
N_SETTLED_ROWS = 3          # rule between settled phenotypes and Experiment A

#: the 1000G PRS-CS run stays in results_v2/prscs/ for the record but is not
#: drawn -- the UKB-LD re-run supersedes it (user decision 2026-09-09; the
#: 1000G-vs-UKB contrast lives in the footnote and README §13.8 discussion)
METHODS = ["ct", "prscs_ukbb", "sbayesr"]
METHOD_LABEL = {"ct": "C+T", "prscs_ukbb": "PRS-CS (UKB LD)",
                "sbayesr": "SBayesR"}
METHOD_COLOR = {"ct": "0.30", "prscs_ukbb": "#0072B2", "sbayesr": "#D55E00"}
DODGE = {"ct": 0.24, "prscs_ukbb": 0.0, "sbayesr": -0.24}
#: the UKB-LD re-run (hpc_v2 commit 8fe2312's prscs_ukbb.sbatch) -- drawn
#: automatically once its association tables are pulled from CSD3
UKBB_DIR = "hpc_v2/work/results_v2/prscs_ukbb"


def load() -> pd.DataFrame:
    rows = []
    ct = pd.read_csv(REPO / "hpc_v3/prs_tables/prs_association_v3.tsv",
                     sep="\t")
    ct = ct[ct.threshold == "0p5"]
    for _, r in ct.iterrows():
        rows.append(dict(method="ct", disorder=r.disorder, stratum=r.stratum,
                         phenotype=r.phenotype, beta=r.beta, se=r.se,
                         p=r.p, p_corr=min(1.0, r.p * M_EFF), n=r.n))
    # UKB-LD PRS-CS; note the cluster's filename inconsistency for ALZnoAPOE,
    # whose file also tags the disorder column 'ALZnoAPOEukbb'
    ukbb_files = {"SCZ": "prs_association_SCZ_prscs_ukbb.tsv",
                  "MDD": "prs_association_MDD_prscs_ukbb.tsv",
                  "ASD": "prs_association_ASD_prscs_ukbb.tsv",
                  "ALZ": "prs_association_ALZ_prscs_ukbb.tsv",
                  "ALZnoAPOE": "prs_association_ALZnoAPOE_ukbb.tsv"}
    for dis, fname in ukbb_files.items():
        f = REPO / UKBB_DIR / fname
        if not f.exists():
            continue
        cs = pd.read_csv(f, sep="\t")
        for _, r in cs.iterrows():
            rows.append(dict(method="prscs_ukbb", disorder=dis,
                             stratum=r.stratum, phenotype=r.phenotype,
                             beta=r.beta, se=r.se, p=r.p, p_corr=r.p, n=r.n))
    # SBayesR: transcription of README 13.8 (pooled) + the cluster agent's
    # 2026-09-09 report (EUR arm, ALZnoAPOE) -- see the tsv's source column;
    # replace with results_v2/sbayesr/ tables when committed
    sb = pd.read_csv(REPO / "hpc_v3/prs_tables/prs_sbayesr_readme.tsv",
                     sep="\t")
    for _, r in sb.iterrows():
        se = abs(r.beta) / norm.ppf(1 - r.p / 2)
        rows.append(dict(method="sbayesr", disorder=r.disorder,
                         stratum=r.stratum, phenotype=r.phenotype,
                         beta=r.beta, se=se, p=r.p, p_corr=r.p, n=np.nan))
    return pd.DataFrame(rows)


def _fmt_p(p: float) -> str:
    return f"{p:.1e}" if p < 1e-3 else f"{p:.2g}"


STRATUM_MARKER = {"EUR": "o", "full": "D"}    # matches slide_results_v3
STRATUM_DODGE = {"EUR": 0.055, "full": -0.055}


def panel(fig, rect, df, phenos, ylabels, title, show_ylab, xlim=None):
    ax = fig.add_axes(rect)
    n = len(phenos)
    for _, r in df.iterrows():
        if r.phenotype not in phenos:
            continue
        y = (n - 1 - phenos.index(r.phenotype)) + DODGE[r.method] \
            + STRATUM_DODGE[r.stratum]
        c = METHOD_COLOR[r.method]
        ax.errorbar(r.beta, y, xerr=1.96 * r.se, fmt="none", ecolor=c,
                    elinewidth=1.2, capsize=1.5, zorder=2)
        ax.plot(r.beta, y, STRATUM_MARKER[r.stratum], ms=4.2,
                mfc=c if r.p_corr < 0.05 else "white",
                mec=c, mew=1.1, zorder=3)
        if r.p_corr < 0.05:
            sgn = 1 if r.beta >= 0 else -1
            ax.text(r.beta + sgn * r.se * 1.96 * 1.12, y, _fmt_p(r.p_corr),
                    fontsize=SMALL - 1.5, ha="left" if sgn > 0 else "right",
                    va="center", color="0.15", zorder=4)
    for i in range(n):
        if i % 2 == 1:
            ax.axhspan(i - 0.5, i + 0.5, color="0.955", zorder=0)
    if n > N_SETTLED_ROWS and phenos is PHENOS:
        ax.axhline(n - N_SETTLED_ROWS - 0.5, color="0.75", lw=0.8,
                   ls=(0, (4, 3)), zorder=1)
    ax.axvline(0, color="0.4", lw=0.8, zorder=1)
    ax.set_ylim(-0.55, n - 0.45)
    ax.set_yticks(range(n))
    ax.set_yticklabels(ylabels[::-1] if show_ylab else [],
                       fontsize=MID if show_ylab else 0)
    ax.set_xlabel("β per SD of score (95% CI)", fontsize=SMALL, labelpad=2)
    ax.set_title(title, pad=6)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    if xlim:
        ax.set_xlim(*xlim)
    else:
        ax.margins(x=0.22)
    return ax


def draw() -> Path:
    df = load()
    mpl.rcParams.update({
        "font.size": BASE, "axes.titlesize": MID, "axes.labelsize": SMALL,
        "legend.fontsize": SMALL, "xtick.labelsize": SMALL,
        "ytick.labelsize": MID, "axes.titlelocation": "left",
        "figure.facecolor": "white", "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(13.333, 7.5))
    fig.text(0.008, 0.955, "PRS methods compared — SCZ → global slope stands "
             "under C+T and SBayesR; PRS-CS attenuates under BOTH LD panels",
             fontsize=BASE + 1.5)
    fig.text(0.008, 0.915, "filled marker = p < 0.05 after its method's "
             "threshold correction (C+T: m_eff = 4.6; single-score methods: "
             "raw p) · ○ EUR stratum · ◇ pooled", fontsize=SMALL,
             color="0.35")

    panel(fig, [0.135, 0.115, 0.225, 0.71], df[df.disorder == "SCZ"],
          PHENOS, YLAB, "SCZ", True)
    panel(fig, [0.395, 0.115, 0.225, 0.71], df[df.disorder == "MDD"],
          PHENOS, YLAB, "MDD", False)
    # controls share the same 7-row grid so rows align across all panels;
    # points exist only where the disorder was scored on that phenotype.
    # ALZ appears twice: with APOE (nominally significant, ~one locus) and
    # with APOE excluded (null) -- the specificity question is the contrast.
    panel(fig, [0.655, 0.115, 0.095, 0.71], df[df.disorder == "ASD"],
          PHENOS, YLAB, "ASD (control)", False)
    panel(fig, [0.775, 0.115, 0.095, 0.71], df[df.disorder == "ALZ"],
          PHENOS, YLAB, "ALZ (control)", False, xlim=(-0.105, 0.15))
    panel(fig, [0.895, 0.115, 0.095, 0.71], df[df.disorder == "ALZnoAPOE"],
          PHENOS, YLAB, "ALZ, no APOE", False, xlim=(-0.105, 0.15))

    present = [m for m in METHODS if (df.method == m).any()]
    handles = [mpl.lines.Line2D([], [], color=METHOD_COLOR[m], marker="o",
                                ls="", mfc=METHOD_COLOR[m],
                                label=METHOD_LABEL[m]) for m in present]
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(0.995, 0.995), ncol=3, frameon=False,
               fontsize=SMALL, handletextpad=0.4, columnspacing=1.0)

    fig.text(0.008, 0.062,
             "PRS-CS on the UKB LD panel (~375k EUR — 7× LARGER than "
             "SBayesR's 50k) does NOT recover the SCZ association in the "
             "EUR arm (p = 0.12 vs 0.20 under 1000G): §13.8's "
             "LD-panel-size account fails a fortiori — the attenuation is "
             "PRS-CS's prior.",
             fontsize=SMALL - 1.5, color="0.42")
    fig.text(0.008, 0.042,
             "Pooled-arm SE inflation worsens under UKB LD (SCZ SE 0.0267 "
             "pooled vs 0.0157 EUR): expected — no EUR panel of any size "
             "fits the 32% non-EUR subjects; the EUR arm is the inference "
             "arm.",
             fontsize=SMALL - 1.5, color="0.42")
    fig.text(0.008, 0.022,
             "SBayesR: transcribed (README §13.8 + agent report; SEs from β "
             "and p, normal quantile; CSD3 tables pending).  ALZ contrast: "
             "with APOE, nominal under both Bayesian methods (92.7% of "
             "SBayesR's squared weight is APOE); APOE excluded, null — the "
             "polygenic late-onset control is clean.",
             fontsize=SMALL - 1.5, color="0.42")

    out = HERE / "slide_prs_methods.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(out)
    return out


if __name__ == "__main__":
    draw()
    sys.exit(0)
