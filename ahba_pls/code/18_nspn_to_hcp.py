"""
18_nspn_to_hcp.py -- put the Whitaker & Vertes 2016 (NSPN) maps into HCP-MMP
space, so the HCP-MMP arm can be compared with NSPN PLS2 directly instead of
leaving that cell of the figure empty.

Method (as in ~/Git/AHBA/notebooks/MT_whitakervertes.ipynb, cell 24): the NSPN
maps are defined on the 308-region subdivision of Desikan-Killiany
(`500.aparc`), and both that parcellation and HCP-MMP1 have fsaverage
(164k-vertex) annot files, so the transfer is parcels -> vertices -> parcels:
each 308-parcel's value is painted onto its vertices, then averaged within every
HCP-MMP parcel.  neuromaps' parcels_to_vertices / vertices_to_parcels do exactly
this; they are not installed here and they match the 308 values to annot parcels
BY POSITION, so this script does the same two steps with nibabel and matches
them BY NAME instead -- `<region>_part<N>` in the annot against
(hemi, region, n_sub_regions) in the published table.  The name mapping is
asserted to be a bijection (308 <-> 308; `corpuscallosum_part1` carries no value
and is treated as background, like `unknown` and the medial wall).

Two things the vertex route makes visible and that the figure must respect:
  * coverage -- an HCP parcel can straddle the medial wall or the corpus
    callosum, so the fraction of its vertices carrying a 308-value is recorded
    per parcel and parcels below COVER_MIN are set to NaN;
  * it is a resampling, not a measurement: HCP parcels are smaller than
    308-parcels in places, so neighbouring HCP parcels can inherit the same
    308-value and the map is smoother than a native HCP fit would be.

Validation: the same route is run 308 -> Desikan-Killiany (lh/rh.aparc.annot)
and the result is correlated with the DK map this project already uses
(data/reference/nspn_dk_maps_bilateral_34.csv, built from the published
DK-level table).  A high correlation there is the evidence that the vertex
transfer is doing what it should.

Outputs
  data/reference/nspn_hcp_maps.csv        bilateral HCP-MMP maps (label = lh_*)
  data/reference/nspn_hcp_coverage.csv    per-parcel vertex coverage, both hemis
  data/reference/NSPN_HCP.md             provenance note
"""
from __future__ import annotations
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
import nibabel as nib
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "data" / "reference"
AHBA = Path.home() / "Git" / "AHBA"
PARC = AHBA / "data" / "parcellations"
MEASURES = ["PLS2", "CT", "CT_delta", "MT", "MT_delta"]
COVER_MIN = 0.50          # a parcel needs half its vertices covered to be kept
BACKGROUND = ("unknown", "???", "corpuscallosum", "Background")


def annot(path: Path):
    lab, _, names = nib.freesurfer.read_annot(str(path))
    return lab, [n.decode() if isinstance(n, bytes) else n for n in names]


def is_bg(name: str) -> bool:
    return name.startswith(BACKGROUND)


# ---------------------------------------------------------------- 308 values --
tab = pd.read_csv(AHBA / "data" / "whitakervertes2016_complete_308.csv", index_col=0)
# -99 is a missing-value sentinel in the published 308-region table, not a score
# (PLS2 there otherwise spans -0.16..0.11). Left unhandled it survives the
# vertex averaging and wrecks any Pearson correlation while leaving Spearman
# almost intact -- which is exactly what the DK validation showed before this.
SENTINEL = -99
n_sentinel = int((tab[MEASURES] == SENTINEL).sum().sum())
sentinel_rows = tab.loc[(tab[MEASURES] == SENTINEL).any(axis=1), ["hemi", "region", "n_sub_regions"]]
tab[MEASURES] = tab[MEASURES].replace(SENTINEL, np.nan)
val = {(h, r, int(pp)): row for (h, r, pp), row in
       tab.set_index(["hemi", "region", "n_sub_regions"])[MEASURES].iterrows()}
assert len(val) == len(tab) == 308, (len(val), len(tab))


def paint_308(hemi: str) -> tuple[int, dict[str, np.ndarray]]:
    """Per-vertex value of each measure on the 308 parcellation (NaN outside).

    `in_308` marks vertices inside a cortical 308-parcel regardless of whether
    the published table has a value there, so parcel COVERAGE (a geometric
    property) stays separate from a measure being missing."""
    lab, names = annot(PARC / f"{hemi}h.500.aparc.annot")
    out = {m: np.full(lab.shape, np.nan) for m in MEASURES}
    out["in_308"] = np.zeros(lab.shape)
    matched = 0
    for i, nm in enumerate(names):
        if is_bg(nm):
            continue
        m = re.match(r"^(.*)_part(\d+)$", nm)
        assert m, f"unparsed 308 label: {nm}"
        key = (hemi[0], m.group(1), int(m.group(2)))
        assert key in val, f"no published value for {hemi} {nm}"
        matched += 1
        sel = lab == i
        out["in_308"][sel] = 1.0
        for meas in MEASURES:
            out[meas][sel] = val[key][meas]
    return matched, out


def to_target(vert: dict[str, np.ndarray], hemi: str, target: str):
    """Average each measure within every parcel of the target annot."""
    lab, names = annot(PARC / f"{hemi}h.{target}.annot")
    rows = []
    for i, nm in enumerate(names):
        if is_bg(nm):
            continue
        sel = lab == i
        n_v = int(sel.sum())
        if n_v == 0:
            continue
        r = dict(hemi=hemi, parcel=nm, n_vertices=n_v,
                 coverage=float(vert["in_308"][sel].mean()))
        for meas in MEASURES:
            v = vert[meas][sel]
            r[meas] = float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
            r[f"valid_{meas}"] = float(np.isfinite(v).mean())
        rows.append(r)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------- transfer --
parts = {}
frames = {"HCPMMP1": [], "aparc": []}
for hemi in ("l", "r"):
    n_matched, vert = paint_308(hemi)
    parts[hemi] = n_matched
    for target in frames:
        frames[target].append(to_target(vert, hemi, target))
assert sum(parts.values()) == 308, parts
HCP = pd.concat(frames["HCPMMP1"], ignore_index=True)
DKv = pd.concat(frames["aparc"], ignore_index=True)

HCP["region"] = HCP.parcel.str.replace(r"^[LR]_", "", regex=True).str.replace("_ROI$", "", regex=True)
HCP.to_csv(REF / "nspn_hcp_coverage.csv", index=False, float_format="%.5g")

kept = HCP[HCP.coverage >= COVER_MIN]
bilat = (kept.groupby("region")[MEASURES].mean()
         .join(kept.groupby("region").coverage.min().rename("min_coverage"))
         .join(kept.groupby("region").size().rename("n_hemispheres")))
bilat.index = "lh_" + bilat.index
bilat.index.name = "label"
bilat.to_csv(REF / "nspn_hcp_maps.csv", float_format="%.6g")

# ----------------------------------------------------------------- validation --
DKv["region"] = DKv.parcel
dk_bilat = (DKv[DKv.coverage >= COVER_MIN].groupby("region")[MEASURES].mean())
dk_ref = pd.read_csv(REF / "nspn_dk_maps_bilateral_34.csv", index_col=0)
dk_ref.index = dk_ref.index.str.replace("^lh_", "", regex=True)
shared = dk_bilat.index.intersection(dk_ref.index)
vrows = []
for meas, refcol in (("PLS2", "PLS2"), ("CT", "CT"), ("CT_delta", "CT_delta"),
                     ("MT", "MT"), ("MT_delta", "MT_delta")):
    if refcol not in dk_ref.columns:
        continue
    a, b = dk_bilat.loc[shared, meas], dk_ref.loc[shared, refcol]
    ok = a.notna() & b.notna()
    vrows.append(dict(measure=meas, n=int(ok.sum()),
                      pearson_r=stats.pearsonr(a[ok], b[ok]).statistic,
                      spearman_rho=stats.spearmanr(a[ok], b[ok]).statistic))
V = pd.DataFrame(vrows)

note = f"""# NSPN maps in HCP-MMP space

`code/18_nspn_to_hcp.py` -> `nspn_hcp_maps.csv`, `nspn_hcp_coverage.csv`.

Route: 308-region (`500.aparc`) parcel values -> fsaverage 164k vertices ->
HCP-MMP1 parcels, both parcellations read from their fsaverage annot files in
`~/Git/AHBA/data/parcellations/`. This is the method used in
`~/Git/AHBA/notebooks/MT_whitakervertes.ipynb` (neuromaps
`parcels_to_vertices` + `vertices_to_parcels`), reimplemented with nibabel so
the 308 values are matched to annot parcels **by name**
(`<region>_part<N>` against the published table's hemi / region /
n_sub_regions) rather than by row position.

- 308 of 308 published parcels matched ({parts['l']} left, {parts['r']} right);
  {n_sentinel} cell(s) held the -99 missing-value sentinel and were set to NaN
  ({', '.join(f"{r.hemi}h {r.region} part {r.n_sub_regions}" for r in sentinel_rows.itertuples())});
  `corpuscallosum_part1` carries no published value and is treated as
  background, as are `unknown` / `???` / the medial wall.
- Parcels keep a value only if at least {COVER_MIN:.0%} of their vertices carry
  one: {len(HCP[HCP.coverage >= COVER_MIN])} of {len(HCP)} HCP-MMP parcels across both hemispheres
  ({len(bilat)} regions after bilateral averaging, of 180).
- Bilateral average of the two hemispheres, labelled `lh_<region>`, to match
  this project's other HCP maps.

**Validation** — the same vertex route run 308 -> Desikan-Killiany, against the
published DK-level table this project already uses
(`nspn_dk_maps_bilateral_34.csv`), over {int(V.n.iloc[0])} regions:

{V.round(3).to_string(index=False)}

**Caveat.** This is a resampling, not a measurement in HCP-MMP space: where an
HCP parcel is smaller than the 308-parcel containing it, neighbouring HCP
parcels inherit the same value, so the map is smoother than a native HCP fit.
Spatial statistics against it (spin tests) are therefore conservative in the
sense that the effective spatial resolution is the 308 parcellation, not 180.
"""
(REF / "NSPN_HCP.md").write_text(note)

print(V.round(4).to_string(index=False))
print(f"\nHCP parcels kept: {len(HCP[HCP.coverage >= COVER_MIN])}/{len(HCP)} "
      f"(bilateral regions: {len(bilat)})")
print("lowest coverage kept:", HCP[HCP.coverage >= COVER_MIN].coverage.min().round(3))
print("dropped parcels:", HCP.loc[HCP.coverage < COVER_MIN, "parcel"].tolist())
print(bilat[MEASURES].describe().round(3).to_string())
