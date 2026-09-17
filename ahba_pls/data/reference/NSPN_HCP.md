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

**Validation** — the same vertex route run 308 -> Desikan-Killiany and averaged
over both hemispheres, against the published DK-level table this project already
uses (`nspn_dk_maps_bilateral_34.csv`), over 34 regions:

 measure  n  pearson_r  spearman_rho
    PLS2 34      0.985         0.971
      CT 34      1.000         0.997
CT_delta 34      0.999         0.998
      MT 34      0.998         0.995
MT_delta 34      1.000         0.999

**How lossy is it?** Measured, not assumed, in `code/19_resample_check.py`:

- The target grid is **finer** than the source, not coarser: 152 parcels per
  hemisphere in the 308 scheme (308 is the bilateral count) against 180 in
  HCP-MMP.
- Each HCP parcel is a vertex-weighted average of a median of **4** source
  parcels (3 contributing >=5% of its vertices), and the largest contributor
  supplies a median of 52% of the vertices (IQR 40-67%). Only 12 of 180 parcels
  take >90% from one source parcel, so "neighbouring parcels inherit one value"
  is not what happens.
- An HCP-native map pushed through the 308 grid and straight back returns at
  r = 0.90-0.91 (ABCD dCT+CT scores, AHBA C3, the thinning rate), which bounds
  what the resampling alone can attenuate.
- Routing 308 -> HCP -> DK instead of 308 -> DK directly costs r = 0.96.

So the weak agreement between this map and the HCP-MMP fits is **not** a
resampling artefact. It is the ABCD component that differs between
parcellations: its DK and HCP score maps correlate rho = 0.74 (against 0.91 for
AHBA C3 and 0.85 for the thinning rate across the same two parcellations), and
NSPN PLS2 tracks the DK version.

**Earlier bug, fixed.** `lh.aparc.annot` names its parcels `lh_bankssts` while
`rh.aparc.annot` names the same parcel `bankssts`. The first version of the DK
validation grouped on the raw annot name, so the two hemispheres became separate
keys and the reported correlation compared RIGHT-hemisphere values against the
published bilateral map (and `lh_unknown` / `lh_corpuscallosum` escaped the
background filter). With the hemisphere prefix normalised the validation is
bilateral-vs-bilateral, and PLS2 agreement rises from r = 0.97 to r = 0.985.
