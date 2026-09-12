"""Figures for the multi-ancestry pre-GRM diagnostics report.

Palette: the validated categorical order (blue, orange, aqua, yellow, magenta,
violet) run through the data-viz validator -- all six checks pass on the
adjacent pairlist, which is the one stacked bars and grouped bars use.  Three
slots sit below 3:1 contrast on a light surface, so the relief rule applies and
every figure ships direct labels plus a data table in the report.

Categorical hues are assigned to a FIXED category order and never cycled, so a
category keeps its colour across every facet and figure.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, numpy as np
from pathlib import Path

FIG = Path("docs/multianc/fig"); FIG.mkdir(parents=True, exist_ok=True)
SURF = "#fcfcfb"; INK = "#0b0b0b"; INK2 = "#52514e"; GRID = "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": GRID, "grid.linewidth": 0.6,
})

def strat_order(strata):
    """EURlike first, then clusters by number -- stable across facets."""
    return ["EURlike"] + sorted([s for s in strata if s != "EURlike"])

# ---------------------------------------------------------------- Figure 1
# Self-reported race/ethnicity composition of each genetic cluster, per k.
# 100% stacked columns: the question is COMPOSITION, so the part-to-whole form
# is right and the shared 0-100 axis makes facets comparable.
sr = pd.read_csv("docs/multianc/selfreport_by_cluster_by_k.tsv", sep="\t")
races = ["White", "Black", "Asian", "AmInd", "Multiple", "Other/NA"]
ks = sorted(sr.k.unique())
fig, axes = plt.subplots(1, len(ks), figsize=(3.5 * len(ks), 4.3), sharey=True)
for ax, k in zip(np.atleast_1d(axes), ks):
    d = sr[(sr.k == k) & (sr.race.isin(races))]
    order = strat_order(d.stratum.unique())
    piv = d.pivot(index="stratum", columns="race", values="pct").reindex(order)[races].fillna(0)
    ns = {r.stratum: r.n for r in sr[sr.k == k].itertuples()}
    bottom = np.zeros(len(piv))
    x = np.arange(len(piv))
    for j, r in enumerate(races):
        vals = piv[r].values
        ax.bar(x, vals, bottom=bottom, color=SLOT[j], width=0.68,
               edgecolor=SURF, linewidth=1.6, label=r if k == ks[0] else None)
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 8:   # direct-label only substantial segments
                ax.text(xi, b + v / 2, f"{v:.0f}", ha="center", va="center",
                        fontsize=8, color="white", fontweight="bold")
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\nn={ns.get(s,0)}" for s in piv.index], fontsize=8)
    ax.set_title(f"k = {k}", fontsize=11, color=INK, pad=8)
    ax.set_ylim(0, 100); ax.yaxis.grid(True); ax.set_axisbelow(True)
np.atleast_1d(axes)[0].set_ylabel("% of cluster (self-reported)")
handles, labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False,
           bbox_to_anchor=(0.5, -0.02), fontsize=9)
fig.suptitle("Genetic clusters recover self-reported groups at every k",
             fontsize=12.5, y=0.99, color=INK)
fig.tight_layout(rect=[0, 0.06, 1, 0.96])
fig.savefig(FIG / "fig1_selfreport_by_k.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("wrote fig1")

# ---------------------------------------------------------------- Figure 2
# Retention through each relatedness-pruning rule, by stratum, per k.
# Grouped bars: two methods compared within each stratum; the comparison is
# magnitude against a common 0-100 scale, so one axis, no dual scale.
gd = pd.read_csv("docs/multianc/grm_diagnostics_by_k.tsv", sep="\t")
fig, axes = plt.subplots(1, len(ks), figsize=(3.5 * len(ks), 4.0), sharey=True)
for ax, k in zip(np.atleast_1d(axes), ks):
    d = gd[(gd.k == k) & (gd.stratum != "ALL")]
    piv = d.pivot(index="stratum", columns="metric", values="value")
    order = strat_order(piv.index)
    piv = piv.reindex(order)
    x = np.arange(len(piv)); w = 0.36
    for off, col, slot, lab in ((-w/2, "pct_kept_gctacutoff", SLOT[1], "--grm-cutoff 0.05"),
                                (+w/2, "pct_kept_pcair",      SLOT[0], "PC-AiR")):
        vals = piv[col].values
        ax.bar(x + off, vals, width=w, color=slot, edgecolor=SURF, linewidth=1.4,
               label=lab if k == ks[0] else None)
        for xi, v in zip(x + off, vals):
            ax.text(xi, v + 2, f"{v:.0f}", ha="center", va="bottom", fontsize=7.5, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\nn={int(piv.loc[s,'n_in_grm'])}" for s in piv.index], fontsize=8)
    ax.set_title(f"k = {k}", fontsize=11, color=INK, pad=8)
    ax.set_ylim(0, 100); ax.yaxis.grid(True); ax.set_axisbelow(True)
np.atleast_1d(axes)[0].set_ylabel("% of stratum retained as 'unrelated'")
handles, labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, -0.02), fontsize=9)
fig.suptitle("--grm-cutoff removes non-European strata, not relatives — at every k",
             fontsize=12.5, y=0.99, color=INK)
fig.tight_layout(rect=[0, 0.07, 1, 0.96])
fig.savefig(FIG / "fig2_retention_by_k.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("wrote fig2")

# ---------------------------------------------------------------- Figures 3&4
ks_df = pd.read_csv("hpc/work/results/impqual/ksweep_chr22.tsv", sep="\t")
KS = sorted(ks_df.k.unique())

# Fig 3: where the pooled-MAF filter's losses fall.
# Stacked columns of the LOST variants attributed to the stratum each is
# commonest in -- the question is composition of a loss, so part-to-whole.
fig, axes = plt.subplots(1, len(KS), figsize=(3.3 * len(KS), 4.0), sharey=True)
for ax, k in zip(np.atleast_1d(axes), KS):
    d = ks_df[(ks_df.k == k) & (ks_df.metric == "lost_attributed")]
    order = strat_order(d.stratum.unique())
    d = d.set_index("stratum").reindex(order)
    tot = int(ks_df[(ks_df.k == k) & (ks_df.metric == "n_lost_by_pooled")].value.iloc[0])
    com = int(ks_df[(ks_df.k == k) & (ks_df.metric == "n_stratum_common")].value.iloc[0])
    bottom = 0.0
    for j, s in enumerate(order):
        v = float(d.loc[s, "value"])
        ax.bar([0], [v], bottom=[bottom], color=SLOT[j], width=0.5,
               edgecolor=SURF, linewidth=1.6, label=s if k == KS[0] else None)
        # Direct-label only segments tall enough to hold two lines of text.
        # At 4% the label overflowed the segment and was clipped by the axis --
        # caught by looking at the rendered PNG, which is why that step exists.
        # Smaller segments are carried by the legend and the table instead.
        if v / max(tot, 1) > 0.10:
            ax.text(0, bottom + v / 2, f"{s}\n{v:.0f} ({100*v/tot:.0f}%)",
                    ha="center", va="center", fontsize=8.5, color="white", fontweight="bold")
        bottom += v
    ax.set_title(f"k = {k}\n{tot} of {com} lost ({100*tot/com:.0f}%)", fontsize=10, color=INK, pad=8)
    ax.set_xticks([]); ax.yaxis.grid(True); ax.set_axisbelow(True)
    ax.set_ylim(0, 500)
    if k == KS[0]:
        ax.legend(loc="upper left", frameon=False, fontsize=8, ncol=1)
np.atleast_1d(axes)[0].set_ylabel("variants common in a stratum but\nfailing the POOLED MAF filter")
fig.suptitle("A pooled MAF filter discards mostly non-European variation",
             fontsize=12.5, y=0.99, color=INK)
fig.tight_layout(rect=[0, 0.02, 1, 0.93])
fig.savefig(FIG / "fig3_maf_loss_by_k.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("wrote fig3")

# Fig 4: imputation uncertainty among variants common in each stratum.
fig, axes = plt.subplots(1, len(KS), figsize=(3.3 * len(KS), 4.0), sharey=True)
for ax, k in zip(np.atleast_1d(axes), KS):
    d = ks_df[(ks_df.k == k) & (ks_df.metric == "pct_uncertain_ge05")]
    order = [s for s in strat_order(d.stratum.unique())]
    d = d.set_index("stratum").reindex(order)
    x = np.arange(len(order))
    ax.bar(x, d.value.values, color=[SLOT[j] for j in range(len(order))],
           width=0.62, edgecolor=SURF, linewidth=1.4)
    for xi, v in zip(x, d.value.values):
        ax.text(xi, v + 0.4, f"{v:.1f}", ha="center", va="bottom", fontsize=8, color=INK2)
    ax.axhline(float(d.loc["EURlike", "value"]), color=INK2, lw=1, ls=(0, (4, 3)))
    ax.set_xticks(x); ax.set_xticklabels(order, fontsize=8, rotation=20, ha="right")
    ax.set_title(f"k = {k}", fontsize=11, color=INK, pad=8)
    ax.yaxis.grid(True); ax.set_axisbelow(True); ax.set_ylim(0, 20)
np.atleast_1d(axes)[0].set_ylabel("% of stratum-common variants\nsubstantially uncertain (ambiguity ≥ 0.05)")
fig.suptitle("Imputation quality by stratum — dashed line = European-like reference",
             fontsize=12.5, y=0.99, color=INK)
fig.tight_layout(rect=[0, 0.02, 1, 0.94])
fig.savefig(FIG / "fig4_impqual_by_k.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("wrote fig4")
