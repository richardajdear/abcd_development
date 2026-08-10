#!/usr/bin/env python
"""Regenerate the figures in ``docs/figures/`` that are computed from tables.

Companion to ``regen_h2_tables.py`` and ``regen_report_tables.py``.  Every
figure here is drawn from a committed CSV in ``docs/``, never from a run
directory, so a figure can never disagree with the table the report cites
beside it -- which is exactly what happened when these were drawn in ad-hoc
cells: three figures survived a correction to the numbers underneath them, and
two cited figures did not exist on disk at all.

Brain-surface maps are NOT here.  They need the ggseg geometry and are written
by ``abcd.brainplot``; see ``tools/regen_brain_maps.py``.

Usage
-----
    PYTHONPATH=src python tools/regen_report_figures.py
    PYTHONPATH=src python tools/regen_report_figures.py --only reliability_design_grid

Every figure is saved at 300 dpi with a geometric overlap check (see
``_verify``); a figure with overlapping text raises rather than being written.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DPI = 300
# Colour is bound to release, and reused for every mark representing it.
C_70 = "#1f6f8b"
C_51 = "#b0b7bd"
C_FOCAL = "#1f6f8b"
C_MUTED = "#9aa5ad"
C_ALARM = "#c1443c"


def _verify(fig: plt.Figure, name: str) -> None:
    """Fail loudly on overlapping visible text rather than shipping it."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    bad = [(a.get_text(), b.get_text())
           for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:]
           if ba.overlaps(bb)]
    if bad:
        raise AssertionError(f"{name}: overlapping text {bad[:4]}")


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

def reliability_design_grid(docs: Path) -> plt.Figure:
    """Effective N and reliability across the design grid.

    The point of the panel is that reliability and effective N move in
    OPPOSITE directions as the visit filter tightens, so plotting them on one
    axis pair is the whole argument.
    """
    g = pd.read_csv(docs / "reliability_grid.csv")
    g = g[~g.family_effect].copy()          # family effect is a separate figure
    g["spec"] = np.where(g.global_cov == "none", "no global covariate",
                         "global-adjusted")
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.7))
    for ax, (col, lab) in zip(axes, [("mean_reliability", "mean slope reliability"),
                                     ("effective_N", "effective N")]):
        for spec, sub in g.groupby("spec"):
            sub = sub.sort_values("min_visits")
            focal = spec == "no global covariate"
            ax.plot(sub.min_visits, sub[col], "-o",
                    color=C_FOCAL if focal else C_MUTED,
                    lw=2.0 if focal else 1.2, ms=5 if focal else 4,
                    zorder=3 if focal else 2, label=spec)
        ax.set_xlabel("minimum visits required")
        ax.set_ylabel(lab)
        ax.set_xticks([2, 3, 4])
        ax.margins(0.12)
    axes[0].set_title("Longer follow-up measures each slope better", loc="left")
    axes[1].set_title("...but costs more subjects than it gains", loc="left")
    peak = (g[g.global_cov == "none"].set_index("min_visits").effective_N.loc[2])
    axes[1].annotate("power peaks at ≥2 visits", xy=(2, peak),
                     xytext=(18, -6), textcoords="offset points", color=C_FOCAL,
                     fontsize=6, va="center",
                     arrowprops=dict(arrowstyle="-", color=C_FOCAL, lw=0.8))
    axes[0].legend(frameon=False, loc="lower right")
    fig.tight_layout()
    return fig


def release_power_tradeoff(docs: Path) -> plt.Figure:
    """Cross-release effective N, decomposed into subjects x reliability."""
    r = pd.read_csv(docs / "handoff_release_comparison.csv")
    # Three panels at this width cannot carry the full run names; the release
    # is the categorical variable and the filter the ordinal one, so label with
    # the filter and key release to colour (threaded through all three panels).
    lab = [("baseline" if v == 2 and rel == 5.1 else f"≥{v}" if v < 4 else "4")
           for rel, v in zip(r.release.astype(float), r.min_visits)]
    lab = [f"5.1\n≥2" if rel == 5.1 else f"7.0\n{l}"
           for rel, l in zip(r.release.astype(float), ["≥2", "≥2", "≥3", "4"])]
    col = [C_51 if x == 5.1 else C_70 for x in r.release.astype(float)]
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.6))
    for ax, (c, t) in zip(axes, [("subjects", "subjects"),
                                 ("mean_reliability", "mean slope reliability"),
                                 ("effective_N", "effective N")]):
        ax.bar(range(len(r)), r[c], color=col, width=0.66)
        ax.set_xticks(range(len(r)))
        ax.set_xticklabels(lab)
        ax.set_ylabel(t)
        ax.margins(y=0.16)
        ax.tick_params(axis="x", length=0)
    for i, v in enumerate(r.effective_N):
        axes[2].annotate(f"{r['vs_5.1'].iloc[i]:.2f}x", (i, v),
                         ha="center", va="bottom", xytext=(0, 2),
                         textcoords="offset points")
    handles = [mpl.patches.Patch(color=C_51, label="release 5.1"),
               mpl.patches.Patch(color=C_70, label="release 7.0")]
    # Upper right is the only empty corner of the counts panel (bars descend
    # left to right); lower left sits on top of the tallest bar.
    axes[0].legend(handles=handles, frameon=False, loc="upper right", fontsize=6)
    axes[0].set_title("Subject counts are nearly equal", loc="left")
    axes[1].set_title("7.0 measures slopes better", loc="left")
    axes[2].set_title("so 7.0 wins on power at every filter", loc="left")
    fig.tight_layout()
    return fig


def family_effect_heritability(docs: Path) -> plt.Figure:
    """The family random effect makes the positive control impossible.

    h2 > 1 is not a large estimate, it is out of the parameter space, which is
    why this figure draws the h2 = 1 bound as an alarm line.
    """
    fe = pd.read_csv(docs / "h2_family_effect_contrast.csv")
    order = ["global mean slope", "baseline thickness (control)"]
    fig, ax = plt.subplots(figsize=(5.0, 2.9))
    for j, (on, colr, lab) in enumerate([(False, C_FOCAL, "family effect off"),
                                         (True, C_ALARM, "family effect on")]):
        sub = fe[fe.family_effect == on].set_index("phenotype").loc[order]
        x = np.arange(len(order)) + (j - 0.5) * 0.3
        ax.errorbar(x, sub.h2,
                    yerr=[sub.h2 - sub.ci_lo, sub.ci_hi - sub.h2],
                    fmt="o", color=colr, ms=6, lw=1.4, capsize=3, label=lab)
    ax.axhline(1.0, color=C_ALARM, ls=":", lw=1.2, zorder=1)
    ax.annotate("h² = 1: upper bound of the parameter space", (len(order) - 0.5, 1.0),
                ha="right", va="bottom", color=C_ALARM, xytext=(0, 3),
                textcoords="offset points")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["global mean\nslope", "baseline thickness\n(control)"])
    ax.set_ylabel("Falconer h²")
    ax.set_title("The family effect returns an impossible h² on the control",
                 loc="left")
    ax.legend(frameon=False, loc="upper left")
    ax.margins(x=0.22, y=0.14)
    fig.tight_layout()
    return fig


def heritability_noglobal(docs: Path) -> plt.Figure:
    """Candidate phenotype heritabilities, with the DZ-class correction shown."""
    h = pd.read_csv(docs / "h2_candidate_phenotypes.csv")
    d = pd.read_csv(docs / "h2_dz_class_comparison.csv").set_index("phenotype")
    h = h.sort_values("h2", ascending=True).reset_index(drop=True)
    y = np.arange(len(h))
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.9),
                             gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    ctrl = h.phenotype.str.contains("control")
    ax.errorbar(h.h2[~ctrl], y[~ctrl],
                xerr=[(h.h2 - h.ci_lo)[~ctrl], (h.ci_hi - h.h2)[~ctrl]],
                fmt="o", color=C_FOCAL, ms=6, lw=1.4, capsize=3)
    ax.errorbar(h.h2[ctrl], y[ctrl],
                xerr=[(h.h2 - h.ci_lo)[ctrl], (h.ci_hi - h.h2)[ctrl]],
                fmt="s", mfc="white", mec=C_MUTED, ecolor=C_MUTED,
                ms=6, lw=1.4, capsize=3)
    ax.axvline(0, color="0.7", lw=0.8, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([p.replace(" (control)", "\n(control)") for p in h.phenotype])
    ax.set_xlabel("Falconer h²  (95% CI)")
    ax.set_title("Every slope phenotype is heritable", loc="left")
    ax.margins(y=0.16)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4, prune="upper"))

    ax = axes[1]
    # Same phenotypes, same order as the left panel -- the two panels share a
    # y axis by construction, so the order is taken from ``h`` rather than
    # rebuilt (an earlier version dropped the control here and silently
    # mis-aligned every row against its label).
    ph = list(h.phenotype)
    w = 0.34
    for j, (col, colr, lab) in enumerate([("r_MZ", C_FOCAL, "MZ twins"),
                                          ("r_DZ_twin", C_MUTED, "DZ twins"),
                                          ("r_full_sib", C_ALARM, "non-twin sibs")]):
        ax.barh(np.arange(len(ph)) + (j - 1) * w / 1.2, d.loc[ph, col],
                height=w / 1.2, color=colr, label=lab)
    ax.set_yticks(range(len(ph)))
    ax.set_yticklabels([])
    ax.set_xlabel("within-pair correlation")
    ax.set_title("DZ twins ≠ non-twin siblings", loc="left")
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4, prune="lower"))
    ax.legend(frameon=False, loc="lower right")
    ax.margins(y=0.16)
    fig.tight_layout()
    return fig


def gwas_phenotype_priority(docs: Path) -> plt.Figure:
    """Held-out h2 against transcriptional coupling, the actual tradeoff."""
    p = pd.read_csv(docs / "gwas_phenotype_priority.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.9))
    p = p.sort_values("h2_heldout", ascending=True).reset_index(drop=True)
    y = np.arange(len(p))
    # ``definition_uses_data`` marks a phenotype whose DEFINITION is fitted on
    # the sample; only that flag distinguishes the rejected candidate, so it is
    # read from the table rather than pattern-matched on the name.
    unstable = p.definition_uses_data.astype(bool).values
    ax = axes[0]
    ax.barh(y[~unstable], p.h2_heldout[~unstable], color=C_FOCAL, height=0.62)
    ax.barh(y[unstable], p.h2_heldout[unstable], color="white", height=0.62,
            edgecolor=C_ALARM, hatch="///", lw=1.1)
    ax.errorbar(p.h2_heldout, y, xerr=p.h2_heldout_sd, fmt="none",
                ecolor="0.35", lw=1.0, capsize=2.5)
    ax.set_yticks(y)
    ax.set_yticklabels(p.phenotype)
    ax.set_xlabel("held-out h²  (SD over 40 splits)")
    ax.set_title("Heritability barely separates the candidates", loc="left")
    ax.margins(y=0.14)
    # Adjacent panels' end ticks collide at this width; cap the count.
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4, prune="upper"))

    ax = axes[1]
    ahba = p.ahba_rho.abs()
    absent = ahba.isna().values          # no spatial map -> not measurable
    ax.barh(y[~unstable & ~absent], ahba[~unstable & ~absent],
            color=C_FOCAL, height=0.62)
    ax.barh(y[unstable & ~absent], ahba[unstable & ~absent], color="white",
            height=0.62, edgecolor=C_ALARM, hatch="///", lw=1.1)
    for yi in y[absent]:
        ax.annotate("no regional map", (0.02, yi), va="center", fontsize=6,
                    color="0.45")
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel("|ρ| with best AHBA component")
    ax.set_title("Transcriptional coupling does", loc="left")
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4, prune="lower"))
    ax.legend(handles=[mpl.patches.Patch(facecolor="white", edgecolor=C_ALARM,
                                         hatch="///",
                                         label="definition fitted in-sample")],
              frameon=False, loc="lower right", fontsize=6)
    ax.margins(y=0.14)
    fig.tight_layout()
    return fig


def ahba_vs_maps_noglobal(docs: Path) -> plt.Figure:
    """Spin-tested correlation of every developmental map with AHBA C1-C3.

    Grouped bars, one group per map, with spin-test significance marked.  The
    spin null SD is drawn as a shaded band rather than an error bar: it is the
    width of the *null*, not the uncertainty of rho, and drawing it as an error
    bar would invite reading it as a confidence interval.
    """
    a = pd.read_csv(docs / "ahba_vs_maps_noglobal.csv")
    order = ["slope_total", "tau_slope", "slopePC1", "slopePC2", "slopePC3", "h2"]
    a = a[a.map_col.isin(order)].copy()
    a["map_col"] = pd.Categorical(a.map_col, order, ordered=True)
    a = a.sort_values(["map_col", "component"])
    names = a.drop_duplicates("map_col").set_index("map_col")["map"]

    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    comps = ["C1", "C2", "C3"]
    colors = {"C1": C_MUTED, "C2": "#7d8f99", "C3": C_FOCAL}
    w = 0.26
    x = np.arange(len(order))
    for j, c in enumerate(comps):
        sub = a[a.component == c].set_index("map_col").reindex(order)
        xs = x + (j - 1) * w
        ax.bar(xs, sub.rho, width=w, color=colors[c], label=f"AHBA {c}")
        # Null width as a band centred on zero, drawn once per bar position.
        ax.bar(xs, 2 * sub.null_sd_spin, bottom=-sub.null_sd_spin, width=w,
               color="0.5", alpha=0.18, zorder=0)
        for xi, rho, p in zip(xs, sub.rho, sub.p_spin):
            if p < 0.05:
                ax.annotate("*", (xi, rho + (0.04 if rho >= 0 else -0.11)),
                            ha="center", fontsize=9, color="0.2")
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([names[o].replace(" of ", "\nof ").replace("Between-subject", "Between-subj.")
                        for o in order], fontsize=6)
    ax.set_ylabel("Spearman ρ")
    # Title states each map's STRONGEST partner, not its only one: the
    # thinning rate is also significantly coupled to C2 (rho = -0.33,
    # p_spin = 0.003), so "only C3" would misread the table.
    ax.set_title("C3 is strongest for the thinning rate, C2 for slope PC3, "
                 "C1 for slope PC2", loc="left")
    ax.legend(frameon=False, ncols=3, fontsize=6, loc="lower left")
    ax.annotate("* p_spin < 0.05;  grey band = spin-null SD", (0.99, 0.02),
                xycoords="axes fraction", ha="right", fontsize=6, color="0.45")
    fig.tight_layout()
    return fig


def structural_covariance(docs: Path) -> plt.Figure:
    """Slope structural-covariance matrix, lobe-ordered, with its PC spectrum.

    The matrix is plotted on a diverging scale centred on zero so the 13
    negative pairs are visible as such rather than being absorbed into the low
    end of a sequential ramp.
    """
    m = pd.read_csv(docs / "sc_matrix_lobe_ordered.csv").set_index("label")
    m = m.loc[:, m.index]                      # square, same order as rows
    s = pd.read_csv(docs / "sc_matrix_summary.csv").set_index("quantity")["value"]

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    v = float(np.nanmax(np.abs(m.values)))
    im = axes[0].imshow(m.values, cmap="RdBu_r", vmin=-v, vmax=v)
    axes[0].set_xticks([]); axes[0].set_yticks([])
    axes[0].set_title("Slope structural covariance (lobe-ordered)", loc="left")
    cb = fig.colorbar(im, ax=axes[0], fraction=0.046)
    cb.set_label("Pearson r between regional slopes", fontsize=7)

    ax = axes[1]
    lobes = [q.replace("within-lobe r ", "") for q in s.index
             if q.startswith("within-lobe r ")]
    vals = [s[f"within-lobe r {l}"] for l in lobes]
    y = np.arange(len(lobes))
    ax.barh(y, vals, color=C_FOCAL, height=0.6)
    ax.axvline(s["between-lobe mean r"], color=C_ALARM, ls=":", lw=1.2)
    ax.annotate(f"between-lobe mean\n{s['between-lobe mean r']:.2f}",
                (s["between-lobe mean r"], len(lobes) - 0.4), fontsize=6,
                color=C_ALARM, ha="center", va="bottom")
    ax.set_yticks(y); ax.set_yticklabels(lobes, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("mean within-lobe r")
    ax.set_title(f"PC1 explains {s['var_explained PC1']:.0%} of variance "
                 f"({s['var_explained_sq PC1']:.0%} of r²)", loc="left")
    fig.tight_layout()
    return fig


def site_scanner_supplement(docs: Path) -> plt.Figure:
    """Site and scanner variance in the slope phenotype, as a 2x2 panel.

    Top-left is the per-region site ICC on the cortical surface; the remaining
    three are the ICC of each grouping against the between-subject variance it
    would have to rival to matter, the per-region distribution of both ICCs, and
    the variance inflation in subjects who switched scanner manufacturer.

    The map was previously a separate figure.  Combining them puts the spatial
    and summary views of the same quantity side by side, which is how a reader
    checks that the small overall ICC is not hiding a few high regions.
    """
    from abcd import brainplot

    icc = pd.read_csv(docs / "site_scanner_icc.csv")
    byreg = pd.read_csv(docs / "site_scanner_icc_by_region.csv")
    sw = pd.read_csv(docs / "scanner_switching_summary.csv").set_index("quantity")["value"]

    # 2x2 with the per-region map in the top-left: the map is the spatial view
    # of the same quantity the boxplot summarises, so keeping them in one figure
    # lets a reader check that no region is an outlier without changing figures.
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.0))
    axm = axes[0, 0]
    brainplot.plot_dk(byreg.set_index("label")["site_icc"], ax=axm,
                      diverging=False, colorbar=True, label="site ICC",
                      fontsize=6, vminmax=(0, float(byreg.site_icc.max())))
    axm.set_title("Per-region site ICC of regional slope", loc="left")

    axes = np.array([axes[0, 1], axes[1, 0], axes[1, 1]])
    ax = axes[0]
    ax.bar(np.arange(len(icc)), icc.icc, color=C_FOCAL, width=0.55)
    ax.set_xticks(np.arange(len(icc)))
    ax.set_xticklabels([f"{g}\n(k={k})" for g, k in zip(icc.grouping, icc.k)],
                       fontsize=6)
    ax.set_ylabel("ICC of global slope")
    for xi, (v, p) in enumerate(zip(icc.icc, icc.p)):
        ax.annotate(f"{v:.4f}\np = {p:.3g}", (xi, v), ha="center", va="bottom",
                    fontsize=6, xytext=(0, 2), textcoords="offset points")
    ax.set_ylim(0, max(icc.icc) * 1.6)
    # Bars of 2 categories on a half-width axis render as slabs; pad the x range
    # so the bar width carries no accidental visual weight.
    ax.set_xlim(-0.75, len(icc) - 0.25)
    ax.set_title(f"Both groupings explain < {max(icc.icc)*100:.0f}%", loc="left")

    ax = axes[1]
    ax.boxplot([byreg.site_icc, byreg.manufacturer_icc], widths=0.55,
               tick_labels=["site", "scanner\nmanufacturer"],
               medianprops=dict(color=C_FOCAL), showfliers=True,
               flierprops=dict(marker=".", markersize=3, mfc="0.5", mec="0.5"))
    ax.set_ylabel("per-region ICC")
    ax.set_title(f"No region exceeds "
                 f"{max(byreg.site_icc.max(), byreg.manufacturer_icc.max()):.2f}",
                 loc="left")

    ax = axes[2]
    ax.bar([0], [sw["var ratio switch/same (global slope)"]], color=C_ALARM,
           width=0.5)
    ax.axhline(1.0, color="0.3", lw=0.8, ls=":")
    ax.set_xticks([0]); ax.set_xticklabels(["switched\nmanufacturer"], fontsize=7)
    ax.set_ylabel("slope variance ratio vs non-switchers")
    ax.set_ylim(0, 1.7)
    ax.set_xlim(-1.0, 1.0)   # single bar: keep it from spanning the panel
    ax.annotate(f"{sw['var ratio switch/same (global slope)']:.2f}×\n"
                f"Cohen d = {sw['Cohen d switch vs same']:.2f}\n"
                f"{sw['pct switched manufacturer']:.1%} of subjects switched",
                (0, sw["var ratio switch/same (global slope)"]), ha="center",
                va="bottom", fontsize=6, xytext=(0, 3),
                textcoords="offset points")
    ax.set_title("Switchers are noisier, but few", loc="left")
    fig.tight_layout()
    return fig


FIGURES = {
    "reliability_design_grid": reliability_design_grid,
    "ahba_vs_maps_noglobal": ahba_vs_maps_noglobal,
    "structural_covariance": structural_covariance,
    "site_scanner_supplement": site_scanner_supplement,
    "release_power_tradeoff": release_power_tradeoff,
    "family_effect_heritability": family_effect_heritability,
    "heritability_noglobal": heritability_noglobal,
    "gwas_phenotype_priority": gwas_phenotype_priority,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", nargs="*", choices=sorted(FIGURES), default=None)
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    docs, figs = root / "docs", root / "docs" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    try:
        from abcd.figstyle import apply_figure_style
        apply_figure_style(sizes=(8, 7, 6))
    except Exception:
        plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                             "axes.spines.right": False, "figure.dpi": DPI,
                             "axes.titlesize": 8, "legend.fontsize": 7,
                             "xtick.labelsize": 6, "ytick.labelsize": 6})
    for name in (a.only or sorted(FIGURES)):
        fig = FIGURES[name](docs)
        _verify(fig, name)
        out = figs / f"{name}.png"
        fig.savefig(out, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out.relative_to(root)}")


if __name__ == "__main__":
    main()
