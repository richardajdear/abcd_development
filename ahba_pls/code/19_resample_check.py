"""
19_resample_check.py -- does the 308 -> HCP-MMP resampling in 18_nspn_to_hcp.py
lose spatial resolution, and is the NSPN map's weak agreement with the HCP fits
explained by that?

Written to check a claim I had made loosely: that the resampled NSPN map's
"effective resolution is the 308 parcellation, not 180, because neighbouring HCP
parcels inherit one 308-value".  That reasoning is wrong twice over -- 308 is a
BILATERAL count (152 left + 156 right) against HCP-MMP's 180 per hemisphere, so
the target is the FINER grid, and each HCP parcel is a vertex-weighted average
of however many 308-parcels it overlaps, not a copy of one.  This script
measures what actually happens.

Four measurements:

1. MIXING -- for every HCP parcel, how many 308-parcels contribute vertices and
   what fraction comes from the largest contributor.  If the dominant fraction
   were ~1.0 the "inherits one value" story would be right.

2. ROUND TRIP -- take a map that is native to HCP-MMP (the ABCD dCT+CT scores,
   AHBA C3), push it through the same vertex route to 308 and straight back to
   HCP-MMP, and correlate with the original.  This is the decisive number: it is
   how much of an HCP-resolution map survives passing through the 308 grid, and
   it bounds how much of the NSPN attenuation resampling can explain.

3. BACK-PROJECTION -- 308 -> HCP -> DK against 308 -> DK directly, i.e. whether
   detouring through HCP-MMP damages the map at DK resolution.

4. LABEL ALIGNMENT -- that the labels of nspn_hcp_maps.csv are the same objects
   as the HCP arm's parcels, so the correlations are not computed across a
   mislabelled index.

Outputs
  results/resample_check.tsv   the numbers behind each measurement
  results/resample_mixing.csv  per-HCP-parcel mixing detail
"""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np, pandas as pd
import nibabel as nib
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RES, REF, DATA = ROOT / "results", ROOT / "data" / "reference", ROOT / "data"
PARC = Path.home() / "Git" / "AHBA" / "data" / "parcellations"
BACKGROUND = ("unknown", "???", "corpuscallosum", "Background")


def annot(name: str):
    lab, _, names = nib.freesurfer.read_annot(str(PARC / name))
    return lab, [n.decode() if isinstance(n, bytes) else n for n in names]


def strip_hemi(n: str) -> str:
    # lh.aparc.annot labels its parcels `lh_bankssts`, rh.aparc.annot the same
    # parcel `bankssts`, HCPMMP1 `L_V1_ROI` -- normalise before any matching.
    return re.sub(r"_ROI$", "", re.sub(r"^(lh_|rh_|L_|R_)", "", n))


def is_bg(n: str) -> bool:
    return strip_hemi(n).startswith(BACKGROUND)


lab308, nm308 = annot("lh.500.aparc.annot")
labhcp, nmhcp = annot("lh.HCPMMP1.annot")
labdk, nmdk = annot("lh.aparc.annot")
assert lab308.shape == labhcp.shape == labdk.shape, "annots must share the fsaverage mesh"

hcp_region = {i: strip_hemi(n) for i, n in enumerate(nmhcp) if not is_bg(n)}
n_308_lh = sum(not is_bg(n) for n in nm308)
n_hcp_lh = len(hcp_region)
n_dk_lh = sum(not is_bg(n) for n in nmdk)

rows = []

# ---------------------------------------------------------------- 1. mixing ---
mix = []
for i, reg in hcp_region.items():
    sel = labhcp == i
    src = lab308[sel]
    src = src[np.isin(src, [j for j, n in enumerate(nm308) if not is_bg(n)])]
    if src.size == 0:
        mix.append(dict(parcel=reg, n_vertices=int(sel.sum()), n_sources=0,
                        dominant_frac=np.nan, n_sources_5pct=0))
        continue
    _, cnt = np.unique(src, return_counts=True)
    frac = cnt / cnt.sum()
    mix.append(dict(parcel=reg, n_vertices=int(sel.sum()), n_sources=len(cnt),
                    dominant_frac=float(frac.max()),
                    n_sources_5pct=int((frac >= 0.05).sum())))
MIX = pd.DataFrame(mix)
MIX.to_csv(RES / "resample_mixing.csv", index=False, float_format="%.4g")
ok = MIX[MIX.n_sources > 0]
rows += [
    dict(measurement="parcel counts (left hemisphere)", value=f"{n_308_lh} in the 308 scheme, "
         f"{n_hcp_lh} in HCP-MMP, {n_dk_lh} in DK", note="the target grid is FINER than the source"),
    dict(measurement="308-parcels contributing per HCP parcel (median)",
         value=f"{ok.n_sources.median():.0f} (>=5% of vertices: {ok.n_sources_5pct.median():.0f})",
         note="an HCP parcel is a vertex-weighted average of several 308-parcels"),
    dict(measurement="dominant 308-parcel share per HCP parcel (median, IQR)",
         value=f"{ok.dominant_frac.median():.2f} ({ok.dominant_frac.quantile(.25):.2f}-{ok.dominant_frac.quantile(.75):.2f})",
         note="1.00 would mean 'inherits one value'"),
    dict(measurement="HCP parcels whose dominant share > 0.90",
         value=f"{int((ok.dominant_frac > 0.90).sum())} of {len(ok)}", note=""),
]


# ------------------------------------------------------------ 2. round trip ---
def paint(values: pd.Series, lab: np.ndarray, names: list[str], key) -> np.ndarray:
    """Parcel values -> per-vertex values (NaN where no value)."""
    out = np.full(lab.shape, np.nan)
    for i, n in enumerate(names):
        if is_bg(n):
            continue
        k = key(n)
        if k in values.index and np.isfinite(values[k]):
            out[lab == i] = values[k]
    return out


def collect(vert: np.ndarray, lab: np.ndarray, names: list[str], key) -> pd.Series:
    """Per-vertex values -> parcel means."""
    out = {}
    for i, n in enumerate(names):
        if is_bg(n):
            continue
        v = vert[lab == i]
        if np.isfinite(v).any():
            out[key(n)] = float(np.nanmean(v))
    return pd.Series(out)


hcp_key = strip_hemi
p308_key = lambda n: n          # `<region>_part<N>`, no hemisphere prefix
dk_key = strip_hemi

SPTS = pd.read_csv(RES / "design_grid_points_scores.csv")
native = {}
for v in ("dCT + CT", "AHBA C3", "dCT rate"):
    x = SPTS[(SPTS.parcellation == "HCP") & (SPTS.variable == v)]
    native[v] = pd.Series(x.value.to_numpy(), index=x.label.str.replace("^lh_", "", regex=True))

for v, s in native.items():
    back = collect(paint(collect(paint(s, labhcp, nmhcp, hcp_key), lab308, nm308, p308_key),
                         lab308, nm308, p308_key), labhcp, nmhcp, hcp_key)
    sh = s.index.intersection(back.dropna().index)
    rows.append(dict(measurement=f"round trip HCP -> 308 -> HCP: {v}",
                     value=f"r = {stats.pearsonr(s[sh], back[sh]).statistic:.3f}, "
                           f"rho = {stats.spearmanr(s[sh], back[sh]).statistic:.3f} (n = {len(sh)})",
                     note="how much of an HCP-resolution map survives the 308 grid"))

# ------------------------------------------------------- 3. back-projection ---
tab = pd.read_csv(Path.home() / "Git" / "AHBA" / "data" / "whitakervertes2016_complete_308.csv",
                  index_col=0)
tab["PLS2"] = tab["PLS2"].replace(-99, np.nan)
lh = tab[tab.hemi == "l"].copy()
lh["key"] = lh.region + "_part" + lh.n_sub_regions.astype(int).astype(str)
nspn308 = lh.set_index("key")["PLS2"]

v308 = paint(nspn308, lab308, nm308, p308_key)
dk_direct = collect(v308, labdk, nmdk, dk_key)
hcp_mid = collect(v308, labhcp, nmhcp, hcp_key)
dk_via_hcp = collect(paint(hcp_mid, labhcp, nmhcp, hcp_key), labdk, nmdk, dk_key)
sh = dk_direct.dropna().index.intersection(dk_via_hcp.dropna().index)
rows.append(dict(measurement="NSPN PLS2: 308 -> DK direct vs 308 -> HCP -> DK",
                 value=f"r = {stats.pearsonr(dk_direct[sh], dk_via_hcp[sh]).statistic:.3f} (n = {len(sh)})",
                 note="detouring through HCP-MMP costs nothing at DK resolution"))

# ------------------------------------------------------ 4. label alignment ----
nspnH = pd.read_csv(REF / "nspn_hcp_maps.csv", index_col=0)
hcp_arm = SPTS[(SPTS.parcellation == "HCP") & (SPTS.variable == "dCT + CT")].label
shared = set(nspnH.index) & set(hcp_arm)
rows.append(dict(measurement="label alignment, nspn_hcp_maps vs the HCP arm",
                 value=f"{len(shared)} shared of {len(nspnH)} and {hcp_arm.nunique()}",
                 note="labels are Glasser region names in both, no positional matching"))

# ------------------------------------------- 5. attenuation is not specific ---
PS = pd.read_csv(RES / "design_grid_pairs_scores.tsv", sep="\t")
piv = PS.pivot_table(index=["var_x", "var_y"], columns="parcellation", values="rho")
piv["ratio"] = piv.HCP.abs() / piv.DK.abs()
with_nspn = piv[[("NSPN PLS2" in a) or ("NSPN PLS2" in b) for a, b in piv.index]]
without = piv[[("NSPN PLS2" not in a) and ("NSPN PLS2" not in b) for a, b in piv.index]]
rows.append(dict(measurement="|rho| HCP / |rho| DK, pairs WITH NSPN PLS2",
                 value=f"median {with_nspn.ratio.median():.2f} "
                       f"({with_nspn.ratio.min():.2f}-{with_nspn.ratio.max():.2f}, n = {len(with_nspn)})",
                 note=""))
rows.append(dict(measurement="|rho| HCP / |rho| DK, pairs WITHOUT NSPN PLS2",
                 value=f"median {without.ratio.median():.2f} "
                       f"({without.ratio.min():.2f}-{without.ratio.max():.2f}, n = {len(without)})",
                 note="the NSPN-specific drop is not the general parcellation effect"))

# ------------------------------------------- 6. at which spatial scale do the
#                                               two maps actually agree? --------
# Aggregate BOTH maps to DK regions through the same vertex route and correlate
# there, then compare with the correlation at HCP resolution. If the coarse-scale
# correlation is high while the fine-scale one is not, the two maps agree on the
# lobar pattern and disagree on sub-DK detail -- which is a statement about the
# maps, not about the resampling.
nspn_dk_v = collect(v308, labdk, nmdk, dk_key)
for v, s_ in native.items():
    hcp_v = paint(s_, labhcp, nmhcp, hcp_key)
    dk_v = collect(hcp_v, labdk, nmdk, dk_key)
    sh_dk = dk_v.dropna().index.intersection(nspn_dk_v.dropna().index)
    fine = nspnH["PLS2"].rename(lambda x: re.sub("^lh_", "", x))
    sh_hcp = s_.index.intersection(fine.dropna().index)
    rows.append(dict(
        measurement=f"{v} vs NSPN PLS2, both via vertices",
        value=f"DK regions: rho = {stats.spearmanr(dk_v[sh_dk], nspn_dk_v[sh_dk]).statistic:+.2f} "
              f"(n = {len(sh_dk)})   |   HCP parcels: rho = "
              f"{stats.spearmanr(s_[sh_hcp], fine[sh_hcp]).statistic:+.2f} (n = {len(sh_hcp)})",
        note="same two maps, two spatial scales"))

# ------------------------------------ 7. is it the NSPN map or the HCP fit? ----
# Push each HCP-space map down to DK regions and compare with the DK-space
# version of the SAME quantity. This separates three candidate sources of the
# weak HCP-NSPN agreement: the imaging (dCT rate), the component (dCT + CT
# scores) and the reference (C3).
dk_native = {}
for v in ("dCT + CT", "AHBA C3", "dCT rate"):
    x = SPTS[(SPTS.parcellation == "DK") & (SPTS.variable == v)]
    dk_native[v] = pd.Series(x.value.to_numpy(), index=x.label.str.replace("^lh_", "", regex=True))
for v, s_hcp in native.items():
    down = collect(paint(s_hcp, labhcp, nmhcp, hcp_key), labdk, nmdk, dk_key)
    sh = down.dropna().index.intersection(dk_native[v].dropna().index)
    rows.append(dict(measurement=f"HCP {v} pushed to DK vs the native DK {v}",
                     value=f"rho = {stats.spearmanr(down[sh], dk_native[v][sh]).statistic:+.2f} (n = {len(sh)})",
                     note="same quantity, two parcellations"))

T = pd.DataFrame(rows)
T.to_csv(RES / "resample_check.tsv", sep="\t", index=False)
pd.set_option("display.max_colwidth", 90, "display.width", 220)
print(T.to_string(index=False))
