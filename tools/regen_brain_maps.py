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
    """Per-region site ICC -- the scanner-confound supplement's headline map."""
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
    """The five primary developmental maps as one grid."""
    m = pd.read_csv(docs / "developmental_maps_noglobal.csv").set_index("label")
    cols = {"Absolute thinning rate": ("slope_total", True),
            "Between-subject SD of rate": ("tau_slope", False),
            "Slope PC1": ("slopePC1", True),
            "Slope reliability": ("reliability", False),
            "Regional h²": ("h2", False)}
    fig, axes = plt.subplots(len(cols), 1, figsize=(6.8, 3.1 * len(cols)))
    for ax, (title, (col, div)) in zip(axes, cols.items()):
        brainplot.plot_dk(m[col], ax=ax, diverging=div, colorbar=True,
                          label="", fontsize=6,
                          vminmax=None if div else (0, float(m[col].max())))
        ax.set_title(title, loc="left", fontsize=8)
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

    fig, axes = plt.subplots(len(order), 2, figsize=(13.6, 3.1 * len(order)))
    for i, (col, r) in enumerate(best.iterrows()):
        brainplot.plot_dk(m[col], ax=axes[i, 0], diverging=True, label="",
                          fontsize=6)
        axes[i, 0].set_title(r["map"], loc="left", fontsize=8)
        # AHBA components are bilateral (34 regions); mirror onto both
        # hemispheres so the two columns are visually comparable.  The rho
        # beside it is computed on bilateral averages of BOTH maps -- hence
        # n_regions = 34 in the table, not 68.
        c = comps[r["component"]]
        mirrored = pd.concat([c.add_prefix("lh_"), c.add_prefix("rh_")])
        brainplot.plot_dk(mirrored, ax=axes[i, 1], diverging=True, label="",
                          fontsize=6)
        star = "*" if r["p_spin"] < 0.05 else ""
        axes[i, 1].set_title(
            f"AHBA {r['component']}   ρ = {r['rho']:+.2f}{star} "
            f"(p_spin = {r['p_spin']:.3f})", loc="left", fontsize=8)
    fig.tight_layout()
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
