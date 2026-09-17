# NSPN maps in HCP-MMP space

`code/18_nspn_to_hcp.py` -> `nspn_hcp_maps.csv`, `nspn_hcp_coverage.csv`.

Route: 308-region (`500.aparc`) parcel values -> fsaverage 164k vertices ->
HCP-MMP1 parcels, both parcellations read from their fsaverage annot files in
`~/Git/AHBA/data/parcellations/`. This is the method used in
`~/Git/AHBA/notebooks/MT_whitakervertes.ipynb` (neuromaps
`parcels_to_vertices` + `vertices_to_parcels`), reimplemented with nibabel so
the 308 values are matched to annot parcels **by name**
(`<region>_part<N>` against the published table's hemi / region /
n_sub_regions) rather than by row position.

- 308 of 308 published parcels matched (152 left, 156 right);
  2 cell(s) held the -99 missing-value sentinel and were set to NaN
  (lh lateraloccipital part 8, rh parahippocampal part 1);
  `corpuscallosum_part1` carries no published value and is treated as
  background, as are `unknown` / `???` / the medial wall.
- Parcels keep a value only if at least 50% of their vertices carry
  one: 357 of 360 HCP-MMP parcels across both hemispheres
  (179 regions after bilateral averaging, of 180).
- Bilateral average of the two hemispheres, labelled `lh_<region>`, to match
  this project's other HCP maps.

**Validation** — the same vertex route run 308 -> Desikan-Killiany, against the
published DK-level table this project already uses
(`nspn_dk_maps_bilateral_34.csv`), over 34 regions:

 measure  n  pearson_r  spearman_rho
    PLS2 34      0.970         0.981
      CT 34      0.996         0.981
CT_delta 34      0.981         0.974
      MT 34      0.993         0.992
MT_delta 34      0.976         0.981

**Caveat.** This is a resampling, not a measurement in HCP-MMP space: where an
HCP parcel is smaller than the 308-parcel containing it, neighbouring HCP
parcels inherit the same value, so the map is smoother than a native HCP fit.
Spatial statistics against it (spin tests) are therefore conservative in the
sense that the effective spatial resolution is the 308 parcellation, not 180.
