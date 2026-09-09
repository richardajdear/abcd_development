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

METHODS = ["ct", "prscs", "prscs_ukbb", "sbayesr"]
METHOD_LABEL = {"ct": "C+T (p<0.5, p̃ = p×4.6)", "prscs": "PRS-CS (1000G LD)",
                "prscs_ukbb": "PRS-CS (UKB LD)", "sbayesr": "SBayesR"}
METHOD_COLOR = {"ct": "0.30", "prscs": "#56B4E9", "prscs_ukbb": "#0072B2",
                "sbayesr": "#D55E00"}
DODGE = {"ct": 0.27, "prscs": 0.09, "prscs_ukbb": -0.09, "sbayesr": -0.27}
#: the UKB-LD re-run (hpc_v2 commit 8fe2312's prscs_ukbb.sbatch) -- drawn
#: automatically once its association tables are pulled from CSD3
UKBB_DIR = "hpc_v2/work/results_v2/prscs_ukbb"


def load() -> pd.DataFrame:
    rows = []
    ct = pd.read_csv(REPO / "hpc_v3/prs_tables/prs_association_v3.tsv",
                     sep="\t")
    ct = ct[(ct.threshold == "0p5") & (ct.stratum == "full")]
    for _, r in ct.iterrows():
        rows.append(dict(method="ct", disorder=r.disorder,
                         phenotype=r.phenotype, beta=r.beta, se=r.se,
                         p=r.p, p_corr=min(1.0, r.p * M_EFF), n=r.n))
    for method, d in (("prscs", REPO / "hpc_v2/work/results_v2/prscs"),
                      ("prscs_ukbb", REPO / UKBB_DIR)):
        for dis in ("SCZ", "MDD", "ASD", "ALZ"):
            f = d / f"prs_association_{dis}_prscs.tsv"
            if not f.exists():      # UKB-LD tables not pulled from CSD3 yet
                continue
            cs = pd.read_csv(f, sep="\t")
            cs = cs[(cs.stratum == "full") & (cs.disorder == dis)]
            for _, r in cs.iterrows():
                rows.append(dict(method=method, disorder=dis,
                                 phenotype=r.phenotype, beta=r.beta, se=r.se,
                                 p=r.p, p_corr=r.p, n=r.n))
    sb = pd.read_csv(REPO / "hpc_v3/prs_tables/prs_sbayesr_readme.tsv",
                     sep="\t")
    for _, r in sb.iterrows():
        se = abs(r.beta) / norm.ppf(1 - r.p / 2)
        rows.append(dict(method="sbayesr", disorder=r.disorder,
                         phenotype=r.phenotype, beta=r.beta, se=se,
                         p=r.p, p_corr=r.p, n=np.nan))
    return pd.DataFrame(rows)


def _fmt_p(p: float) -> str:
    return f"{p:.1e}" if p < 1e-3 else f"{p:.2g}"


def panel(fig, rect, df, phenos, ylabels, title, show_ylab, xlim=None):
    ax = fig.add_axes(rect)
    n = len(phenos)
    for _, r in df.iterrows():
        if r.phenotype not in phenos:
            continue
        y = (n - 1 - phenos.index(r.phenotype)) + DODGE[r.method]
        c = METHOD_COLOR[r.method]
        ax.errorbar(r.beta, y, xerr=1.96 * r.se, fmt="none", ecolor=c,
                    elinewidth=1.4, capsize=2, zorder=2)
        ax.plot(r.beta, y, "o", ms=5, mfc=c if r.p_corr < 0.05 else "white",
                mec=c, mew=1.2, zorder=3)
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
        ax.margins(x=0.14)
    return ax


def draw() -> Path:
    df = load()
    mpl.rcParams.update({
        "font.size": BASE, "axes.titlesize": MID, "axes.labelsize": SMALL,
        "legend.fontsize": SMALL, "xtick.labelsize": SMALL,
        "ytick.labelsize": MID, "axes.titlelocation": "left",
        "figure.facecolor": "white", "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(13.333, 7.5))
    fig.text(0.008, 0.955, "PRS methods compared, pooled arm — SCZ → global "
             "slope stands under C+T and SBayesR; PRS-CS is the outlier",
             fontsize=BASE + 1.5)
    fig.text(0.008, 0.915, "filled marker = p < 0.05 after its method's "
             "threshold correction (C+T: m_eff = 4.6; single-score methods: "
             "raw p)", fontsize=SMALL, color="0.35")

    panel(fig, [0.135, 0.115, 0.255, 0.71], df[df.disorder == "SCZ"],
          PHENOS, YLAB, "SCZ", True)
    panel(fig, [0.425, 0.115, 0.255, 0.71], df[df.disorder == "MDD"],
          PHENOS, YLAB, "MDD", False)
    # controls share the same 7-row grid so rows align across all panels;
    # points exist only where the disorder was scored on that phenotype
    panel(fig, [0.72, 0.115, 0.125, 0.71], df[df.disorder == "ASD"],
          PHENOS, YLAB, "ASD (control)", False)
    panel(fig, [0.875, 0.115, 0.115, 0.71], df[df.disorder == "ALZ"],
          PHENOS, YLAB, "ALZ (control)", False)

    present = [m for m in METHODS if (df.method == m).any()]
    handles = [mpl.lines.Line2D([], [], color=METHOD_COLOR[m], marker="o",
                                ls="", mfc=METHOD_COLOR[m],
                                label=METHOD_LABEL[m]) for m in present]
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(0.995, 0.995), ncol=3, frameon=False,
               fontsize=SMALL, handletextpad=0.4, columnspacing=1.0)

    fig.text(0.008, 0.052,
             "SBayesR: transcribed from hpc_v2/README_HPC.md §13.8 (pooled; "
             "CSD3 tables not yet committed); SEs reconstructed from β and p "
             "(normal quantile), exact for a Wald test.  LD references: "
             "SBayesR 50k UKB EUR · PRS-CS 1000G EUR (503) — the attenuation "
             "ordering SBayesR ≳ C+T > PRS-CS tracks LD-panel size.",
             fontsize=SMALL - 1.5, color="0.42")
    fig.text(0.008, 0.032,
             "ALZ caution: 92.7% of the SBayesR ALZ score's squared weight "
             "is APOE-region variants — its p = 0.039 on global slope is one "
             "locus, not polygenic AD risk.  Experiment A phenotypes are "
             "C+T-only so far; PRS-CS/SBayesR scores exist on CSD3, one "
             "assoc job fills the gaps.",
             fontsize=SMALL - 1.5, color="0.42")

    out = HERE / "slide_prs_methods.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(out)
    return out


if __name__ == "__main__":
    draw()
    sys.exit(0)
