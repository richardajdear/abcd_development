"""Cortical surface maps for DK-parcellated region values.

Rendering uses the polygon outlines bundled with the ``ggseg`` Python package
(a port of the R package of the same name). Those outlines are 2-D schematic
projections -- lateral and medial views of each hemisphere -- not true surface
meshes, so this module needs no FreeSurfer ``fsaverage`` files and no network
access. The trade-off is that it is a schematic, not an inflated surface.

Known limitations of the bundled schematic geometry, all verified against the
polygon files rather than assumed:

* Two DK regions have no polygon at all -- ``frontalpole`` and
  ``temporalpole``. They are reported by :func:`missing_regions` and drawn in
  the "no data" colour rather than silently dropped.
* Ventral-surface regions are visible only as slivers, because neither a
  lateral nor a medial view faces them. ``entorhinal`` (~170 units^2) and
  ``fusiform`` (~320 units^2) are the extreme cases -- 20-40x smaller than a
  typical parcel -- so their colour is essentially unreadable. Never make a
  claim about these regions from this figure alone; read the value off the
  table. :data:`SLIVER_REGIONS` lists them.
* Three regions (``superiorfrontal``, ``superiorparietal``,
  ``inferiortemporal``) legitimately appear in both the lateral and the medial
  view as separate sub-paths of one polygon; both are filled with the same
  value, which is correct, not double-counting.

Typical use::

    from abcd import brainplot
    fig = brainplot.plot_dk(values, diverging=True, label="Age effect")

where ``values`` is a Series indexed by ``lh_``/``rh_``-prefixed DK labels
(the ``label`` column of ``data/region_labels.csv``).
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

#: Regions in the DK atlas with no polygon in the bundled ggseg geometry.
UNSUPPORTED_REGIONS = ("frontalpole", "temporalpole")

#: Ventral regions whose rendered polygon is too small to read a colour from.
#: Measured drawn areas: entorhinal ~1.7e2, fusiform ~3.2e2 units^2, against
#: ~1e4 for a typical parcel.
SLIVER_REGIONS = ("entorhinal", "fusiform", "parahippocampal")

#: Outline files that bound the canvas rather than naming a cortical region.
_OUTLINE_FILES = (
    "lateral_left",
    "medial_left",
    "lateral_right",
    "medial_right",
)

#: Non-cortical / filler polygons that should never be filled with data.
_NON_REGION_FILES = ("NA_left", "NA_right", "corpuscallosum_left",
                     "corpuscallosum_right")


class GeometryMissing(RuntimeError):
    """Raised when the ggseg polygon geometry cannot be located."""


@dataclass(frozen=True)
class DKGeometry:
    """Parsed DK polygon geometry, keyed by ``<region>_<left|right>``."""

    paths: dict          # name -> matplotlib.path.Path
    outlines: dict       # name -> matplotlib.path.Path
    extent: tuple        # (xmin, xmax, ymin, ymax)

    def region_names(self) -> set:
        """DK region stems (no hemisphere suffix) that have geometry."""
        return {n.rsplit("_", 1)[0] for n in self.paths}


def _geometry_dir() -> str:
    try:
        import ggseg
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise GeometryMissing(
            "The 'ggseg' Python package is required for cortical surface maps. "
            "Install it with: pip install ggseg"
        ) from exc
    d = os.path.join(os.path.dirname(ggseg.__file__), "data", "dk")
    if not os.path.isdir(d):
        raise GeometryMissing(f"ggseg is installed but {d} is missing")
    return d


def load_dk_geometry() -> DKGeometry:
    """Parse the bundled DK polygons into matplotlib Paths.

    Returns a :class:`DKGeometry` whose ``paths`` excludes the canvas outlines
    and the non-cortical filler polygons, so every entry is fillable with data.
    """
    from matplotlib.path import Path as MplPath
    from ggseg import _svg_parse_

    d = _geometry_dir()
    paths, outlines = {}, {}
    for name in sorted(os.listdir(d)):
        codes, verts = _svg_parse_(open(os.path.join(d, name)).read())
        p = MplPath(verts, codes)
        if name in _OUTLINE_FILES:
            outlines[name] = p
        elif name in _NON_REGION_FILES:
            continue
        else:
            paths[name] = p

    allv = np.vstack([p.vertices for p in outlines.values()])
    extent = (allv[:, 0].min(), allv[:, 0].max(),
              allv[:, 1].min(), allv[:, 1].max())
    return DKGeometry(paths=paths, outlines=outlines, extent=extent)


def _ggseg_key(label: str) -> str:
    """Map a repo region label (``lh_bankssts``) to a ggseg key."""
    hemi, region = label.split("_", 1)
    side = {"lh": "left", "rh": "right"}[hemi]
    return f"{region}_{side}"


def missing_regions(values: pd.Series) -> list:
    """Labels in ``values`` that have no polygon in the bundled geometry."""
    geom = load_dk_geometry()
    return [lab for lab in values.index
            if _ggseg_key(lab) not in geom.paths]


def plot_dk(values, *, ax=None, cmap=None, diverging=False, center=0.0,
            vminmax=None, label="", title="", nodata_color="0.82",
            edgecolor="white", linewidth=0.4, colorbar=True,
            fontsize=8, figsize=(6.8, 3.4), symmetric=True):
    """Draw a DK-parcellated map of ``values``.

    Parameters
    ----------
    values : pandas.Series or dict
        Indexed by ``lh_``/``rh_``-prefixed DK region labels.
    diverging : bool
        Use a diverging colormap centred on ``center``. When ``symmetric``,
        the limits are set to +/- max|value - center| so the colour scale's
        neutral point is the semantic zero, not the data midpoint.
    center : float
        The semantically meaningful zero for a diverging scale.
    vminmax : (float, float), optional
        Explicit colour limits, overriding the automatic choice.
    label : str
        Colourbar label (include units).
    nodata_color : matplotlib color
        Fill for regions with geometry but no value, and for the two DK
        regions absent from the bundled atlas.

    Returns
    -------
    (fig, ax, mappable)
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    values = pd.Series(values).dropna().astype(float)
    geom = load_dk_geometry()

    if cmap is None:
        cmap = "RdBu_r" if diverging else "viridis"
    cmap = mpl.colormaps[cmap] if isinstance(cmap, str) else cmap

    if vminmax is not None:
        vmin, vmax = vminmax
    elif diverging and symmetric:
        r = float(np.abs(values - center).max())
        vmin, vmax = center - r, center + r
    else:
        vmin, vmax = float(values.min()), float(values.max())
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    xmin, xmax, ymin, ymax = geom.extent
    ax.set_xlim(xmin - 1, xmax + 1)
    ax.set_ylim(ymax + 1, ymin - 1)   # SVG y grows downward
    ax.set_aspect(1)
    ax.set_axis_off()

    wanted = {}
    for lab, v in values.items():
        k = _ggseg_key(lab)
        if k in geom.paths:
            wanted[k] = v

    # every region first in the no-data colour, then overpaint those with values
    for k, p in geom.paths.items():
        ax.add_patch(patches.PathPatch(
            p, facecolor=nodata_color, edgecolor=edgecolor, lw=linewidth))
    for k, v in wanted.items():
        ax.add_patch(patches.PathPatch(
            geom.paths[k], facecolor=cmap(norm(v)),
            edgecolor=edgecolor, lw=linewidth))

    if title:
        ax.set_title(title, fontsize=fontsize, loc="left")

    mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    if colorbar:
        cb = fig.colorbar(mappable, ax=ax, fraction=0.022, pad=0.01,
                          shrink=0.62)
        cb.set_label(label, fontsize=fontsize - 1)
        cb.ax.tick_params(labelsize=fontsize - 2)
        cb.outline.set_visible(False)
    return fig, ax, mappable


def view_labels(ax, geom=None, fontsize=6, color="0.35"):
    """Annotate the four anatomical views. Call after :func:`plot_dk`."""
    geom = geom or load_dk_geometry()
    xmin, xmax, ymin, ymax = geom.extent
    xmid, ymid = (xmin + xmax) / 2, (ymin + ymax) / 2
    for txt, x, y in [("left lateral", xmin + 0.02 * (xmax - xmin), ymid - 0.02 * (ymax - ymin)),
                      ("left medial", xmid + 0.02 * (xmax - xmin), ymid - 0.02 * (ymax - ymin)),
                      ("right lateral", xmin + 0.02 * (xmax - xmin), ymax),
                      ("right medial", xmid + 0.02 * (xmax - xmin), ymax)]:
        ax.text(x, y, txt, fontsize=fontsize, color=color, va="bottom")
