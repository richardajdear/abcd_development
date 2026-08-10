#!/usr/bin/env python
"""Regenerate the DK brain-surface figures in ``docs/figures/``.

Companion to ``regen_report_figures.py``, split off from it because these need
the ggseg geometry via :mod:`abcd.brainplot` while that script needs only
matplotlib.  The rule is the same and it is the point of both scripts: every
map is drawn from a committed CSV in ``docs/``, never from a run directory, so
a figure cannot disagree with the table the report cites beside it.  Before
these scripts existed the maps were drawn in ad-hoc cells, and several
survived corrections to the numbers underneath them.

Usage
-----
    PYTHONPATH=src python tools/regen_brain_maps.py
    PYTHONPATH=src python tools/regen_brain_maps.py --only map_reliability

Colour scales
-------------
Signed quantities (thinning rate, PC loadings) use a diverging scale centred on
the semantic zero with symmetric limits, so hue sign carries meaning.
Non-negative quantities (reliability, h2, ICC) use a sequential scale anchored
at zero, so lightness is comparable across panels rather than being rescaled
per figure.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from abcd import brainplot

DPI = 300

#: Cited in the report as the primary regional maps.  Each entry is
#: (column in developmental_maps_noglobal.csv, colourbar label, diverging?).
_MAPS = {
    "map_slope_total": ("slope_total", "mm/year", True),
    "map_tau_slope": ("tau_slope", "between-subject SD (mm/year)", False),
    "map_slopePC1": ("slopePC1", "slope PC1 loading", True),
    "map_reliability": ("reliability", "slope reliability", False),
    "map_h2": ("h2", "Falconer h²", False),
}


def _single(docs: Path, name: str) -> plt.Figure:
    col, label, diverging = _MAPS[name]
    m = pd.read_csv(docs / "developmental_maps_noglobal.csv").set_index("label")
    vals = m[col]
    # Reliability is undefined (not zero) in boundary-fit regions; plotting a
    # dropped NaN as "no data" grey is the honest rendering.
    fig, ax, _ = brainplot.plot_dk(vals, diverging=diverging, label=label,
                                   center=0.0 if diverging else None,
                                   vminmax=None if diverging else (0, float(vals.max())))
    n_missing = int(vals.isna().sum())
    if n_missing:
        ax.set_title(f"{n_missing} region(s) undefined, shown in grey",
                     loc="right", fontsize=6, color="0.45")
    return fig


def site_icc_map(docs: Path) -> plt.Figure:
    """Per-region site ICC as a standalone map.

    No longer embedded in the report: §9 now shows this as the top-left panel of
    ``site_scanner_supplement``, beside the boxplot of the same quantity.  Kept
    because it is cheap to render and useful on its own (slides, or checking a
    single region), but if you change it, change the panel too.
    """
    t = pd.read_csv(docs / "site_scanner_icc_by_region.csv").set_index("label")
    fig, ax, _ = brainplot.plot_dk(t["site_icc"], diverging=False,
                                   label="site ICC of regional slope",
                                   vminmax=(0, float(t["site_icc"].max())))
    return fig


def _grid(values: dict, labels: dict, diverging: bool, ncols: int = 2) -> plt.Figure:
    n = len(values)
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.8 * ncols, 3.1 * nrows))
    axes = axes.ravel()
    for ax, (k, v) in zip(axes, values.items()):
        brainplot.plot_dk(v, ax=ax, diverging=diverging, colorbar=True,
                          label=labels[k], fontsize=6)
        ax.set_title(k, loc="left", fontsize=8)
    for ax in axes[n:]:
        ax.set_visible(False)
    fig.tight_layout()
    return fig


def developmental_maps_noglobal(docs: Path) -> plt.Figure:
    """The six primary developmental maps, 2 rows x 3 columns.

    Row 1 is what the phenotype *is* -- where cortex starts and how fast it
    thins, in mm and mm/yr.  Row 2 is what can be done with it -- whether the
    per-subject slope is measurable, heritable, and structured.  That grouping
    is the reason for the layout: the two rows answer different questions and
    a reader comparing within a row is comparing like with like.
    """
    m = pd.read_csv(docs / "developmental_maps_noglobal.csv").set_index("label")
    panels = [
        ("Baseline thickness (mm)",        "baseline_thickness", False),
        ("Absolute thinning rate (mm/yr)", "slope_total",        True),
        ("SD of thinning rate (mm/yr)",    "tau_slope",          False),
        ("Slope reliability",              "reliability",        False),
        ("Slope h² (Falconer)",            "h2",                 False),
        ("Slope PC1 loading",              "slopePC1",           True),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 6.0))
    for ax, (title, col, div) in zip(axes.ravel(), panels):
        # Anchoring a sequential scale at 0 is right when 0 is the meaningful
        # floor (reliability, SD, h2: "none of this quantity"), and wrong for
        # baseline thickness, where every value sits between 1.8 and 3.8 mm and
        # a 0-anchored scale spends its whole range on empty space -- the panel
        # renders as one flat colour and the regional pattern disappears.
        # Thinning rate is negative everywhere and h2 dips slightly below 0 from
        # estimation noise, so those must autoscale too or real values clip.
        vminmax = None
        if not div and float(m[col].min()) >= 0 and col != "baseline_thickness":
            vminmax = (0, float(m[col].max()))
        brainplot.plot_dk(m[col], ax=ax, diverging=div, colorbar=True,
                          label="", fontsize=6, vminmax=vminmax)
        ax.set_title(title, loc="left", fontsize=9)
    fig.tight_layout()
    return fig


def ahba_maps_grid(docs: Path) -> plt.Figure:
    """Each developmental map beside the AHBA component it best matches.

    The pairing is read from ``ahba_vs_maps_noglobal.csv`` rather than being
    hard-coded, so the panel follows the table if a correlation changes sign
    or a different component becomes the best match.
    """
    from abcd import genemaps
    m = pd.read_csv(docs / "developmental_maps_noglobal.csv").set_index("label")
    a = pd.read_csv(docs / "ahba_vs_maps_noglobal.csv")
    best = (a.assign(absrho=a.rho.abs())
             .sort_values("absrho", ascending=False)
             .groupby("map_col", as_index=False).first())
    comps = genemaps.ahba_components("dsk")
    order = [c for c in ["slope_total", "tau_slope", "slopePC1", "slopePC2",
                         "slopePC3", "h2"] if c in set(best.map_col)]
    best = best.set_index("map_col").loc[order]

    # Three columns per row: ABCD map, AHBA component, and the scatter behind
    # the rho.  The scatter is what makes the rho auditable -- a coefficient
    # printed in a title cannot show whether it rests on a couple of extreme
    # regions or on a consistent gradient, and the maps alone cannot either.
    # width_ratios keeps the scatter square-ish rather than letting it inherit
    # the wide aspect the surface panels need.
    fig, axes = plt.subplots(len(order), 3, figsize=(16.5, 3.1 * len(order)),
                             gridspec_kw=dict(width_ratios=[1.0, 1.0, 0.52]))
    for i, (col, r) in enumerate(best.iterrows()):
        brainplot.plot_dk(m[col], ax=axes[i, 0], diverging=True, label="",
                          fontsize=6)
        axes[i, 0].set_title(r["map"], loc="left", fontsize=8)
        # AHBA components are bilateral (34 regions); mirror onto both
        # hemispheres so the two columns are visually comparable.  The rho
        # beside it is computed on bilateral averages of BOTH maps -- hence
        # n_regions = 34 in the table, not 68.
        #
        # Mirroring shows the same 34 values twice and so implies no information
        # that is not there, but a reader cannot tell that by looking, which is
        # why the figure footnote says so explicitly.  The alternative -- passing
        # lh_ labels only -- renders the right hemisphere in the no-data grey,
        # which reads as missing rather than as deliberate.
        c = comps[r["component"]]
        mirrored = pd.concat([c.add_prefix("lh_"), c.add_prefix("rh_")])
        brainplot.plot_dk(mirrored, ax=axes[i, 1], diverging=True, label="",
                          fontsize=6)
        star = "*" if r["p_spin"] < 0.05 else ""
        axes[i, 1].set_title(
            f"AHBA {r['component']}   ρ = {r['rho']:+.2f}{star} "
            f"(p_spin = {r['p_spin']:.3f})", loc="left", fontsize=8)

        # --- scatter -------------------------------------------------------
        # Bilateral averages of BOTH maps, matching how rho in the table was
        # computed (n_regions = 34).  Plotting the 68 unilateral regions here
        # would show a different, more optimistic-looking cloud than the
        # coefficient it sits beside.
        ax = axes[i, 2]
        bil = (m[[col]].assign(region=m.index.str.replace(r"^[lr]h_", "", regex=True))
                       .groupby("region")[col].mean())
        j = pd.concat([bil.rename("abcd"), c.rename("ahba")], axis=1).dropna()
        # z within each map: the two axes are in incommensurable units (mm/yr or
        # a PC loading vs an expression component), and z-scoring is monotone so
        # Spearman rho is unchanged.
        z = (j - j.mean()) / j.std()
        ax.axhline(0, color="0.85", lw=0.6, zorder=0)
        ax.axvline(0, color="0.85", lw=0.6, zorder=0)
        ax.scatter(z.abcd, z.ahba, s=14, color="#3b6e8f", alpha=0.75,
                   linewidths=0)
        # Monotone fit on ranks, to match Spearman rather than implying a
        # least-squares relationship the coefficient does not describe.
        rk = z.rank()
        b, a0 = np.polyfit(rk.abcd, rk.ahba, 1)
        xs = np.array([rk.abcd.min(), rk.abcd.max()])
        # Map the rank-space line back onto the z axes for display.
        ax.plot(np.interp(xs, sorted(rk.abcd), sorted(z.abcd)),
                np.interp(a0 + b * xs, sorted(rk.ahba), sorted(z.ahba)),
                color="0.35", lw=1.0, ls="--", zorder=1)
        ax.set_xlabel(f"{r['map']} (z)", fontsize=6)
        ax.set_ylabel(f"AHBA {r['component']} (z)", fontsize=6)
        ax.tick_params(labelsize=6)
        ax.set_title(f"n = {len(z)} bilateral regions", loc="left", fontsize=7)
    fig.tight_layout()
    # Footnote rather than a per-panel note: it applies to every row, and the
    # two facts a reader needs are that the hemispheres are not independent
    # evidence and that a few regions are unrenderable in this atlas.
    fig.text(0.005, 0.002,
             "Both hemispheres show the same 34 bilateral values (mirrored), so "
             "left and right are not independent. ρ and the scatters use the 34 "
             "bilateral regions. Frontal/temporal pole are unrenderable and "
             "entorhinal, fusiform and parahippocampal draw as slivers in this "
             "atlas; see docs/developmental_maps_noglobal.csv for values.",
             fontsize=6, color="0.45", ha="left", va="bottom")
    return fig


FIGURES = {
    **{k: (lambda docs, _k=k: _single(docs, _k)) for k in _MAPS},
    "site_icc_map": site_icc_map,
    "developmental_maps_noglobal": developmental_maps_noglobal,
    "ahba_maps_grid": ahba_maps_grid,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", nargs="*", choices=sorted(FIGURES))
    ap.add_argument("--docs", type=Path,
                    default=Path(__file__).resolve().parents[1] / "docs")
    args = ap.parse_args()

    outdir = args.docs / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    for name in (args.only or sorted(FIGURES)):
        fig = FIGURES[name](args.docs)
        fig.savefig(outdir / f"{name}.png", dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote docs/figures/{name}.png")


if __name__ == "__main__":
    main()
