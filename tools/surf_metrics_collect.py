"""
Assemble the abcd.surf_metrics part files into release-convention wide TSVs (D1).

Reads ``<release>/processed/hcp/surf_metrics_parts/part_*.parquet`` (+ ``*_qc.tsv``) and
writes, for every measure ``m`` sampled:

- ``processed/hcp/mr_y_smri__<m>__hcp.tsv``: ``participant_id, session_id``,
  ``mr_y_smri__<m>__hcp__<region>__<hemi>_mean`` (180 x 2, region tokens as in
  ``mr_y_smri__thk__hcp.tsv``), ``..__lh_mean``, ``..__rh_mean``, ``.._mean``;
- ``processed/dsk_local/mr_y_smri__<m>__dsk.tsv``: the same with the release DK codes;
- ``processed/hcp/surf_metrics_session_qc.tsv``: one row per session x measure;
- ``processed/hcp/PROVENANCE_surf_metrics.txt``.

Parcel values are vertex-area-weighted means; hemisphere / whole-cortex values are
area-weighted over all non-medial-wall vertices (see abcd.surf_metrics).

    PYTHONPATH=src python tools/surf_metrics_collect.py
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import pyarrow.dataset as ds

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abcd import paths  # noqa: E402
from abcd.surf_metrics import MEASURES  # noqa: E402

KEY = ["participant_id", "session_id"]


def wide(long: pd.DataFrame, stem: str) -> pd.DataFrame:
    parc = long[long.region != ""]
    piv = parc.pivot_table(index=KEY, columns=["region", "hemi"], values="mean_w",
                           aggfunc="first", observed=True)
    piv.columns = [f"{stem}__{r}__{h}_mean" for r, h in piv.columns]
    summ = long[long.region == ""].pivot_table(index=KEY, columns="hemi", values="mean_w",
                                               aggfunc="first", observed=True)
    for h in ("lh", "rh"):
        piv[f"{stem}__{h}_mean"] = summ.get(h)
    piv[f"{stem}_mean"] = summ.get("all")
    return piv.reset_index()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--release-dir", type=Path, default=None)
    ap.add_argument("--pattern", default="part_*", help="part-file glob (without suffix)")
    a = ap.parse_args(argv)
    rel = a.release_dir or paths.release_dir("7.0")
    hcp_dir, dsk_dir = rel / "processed" / "hcp", rel / "processed" / "dsk_local"
    parts_dir = hcp_dir / "surf_metrics_parts"
    files = sorted(parts_dir.glob(f"{a.pattern}.parquet"))
    if not files:
        print(f"no part files in {parts_dir}")
        return 1
    print(f"{len(files)} part files", flush=True)

    qc = pd.concat([pd.read_csv(f.with_name(f.stem + "_qc.tsv"), sep="\t") for f in files],
                   ignore_index=True)
    qc.to_csv(hcp_dir / "surf_metrics_session_qc.tsv", sep="\t", index=False)
    print("status by measure:\n" + qc.assign(status=qc.status.str.split(":").str[0])
          .groupby(["measure", "status"]).size().to_string(), flush=True)

    dset = ds.dataset([str(f) for f in files], format="parquet")
    written = []
    for atlas, out in (("hcp", hcp_dir), ("dsk", dsk_dir)):
        for m in sorted(qc.measure.unique()):
            t = dset.to_table(filter=(ds.field("atlas") == atlas) & (ds.field("measure") == m),
                              columns=[*KEY, "hemi", "region", "mean_w"]).to_pandas()
            for c in (*KEY, "hemi", "region"):
                t[c] = t[c].astype(str)
            stem = f"mr_y_smri__{m}__{atlas}"
            w = wide(t, stem).sort_values(KEY)
            path = out / f"{stem}.tsv"
            w.to_csv(path, sep="\t", index=False, float_format="%.4f")
            written.append((path, w.shape))
            print(f"wrote {path}: {w.shape[0]} rows x {w.shape[1]} cols", flush=True)

    import nibabel, numpy, scipy
    ok = qc[qc.status == "ok"].groupby("measure").size().to_dict()
    lines = [
        "Per-parcel surface-sampled cortical measures (D1) -- abcd.surf_metrics",
        f"generated {time.strftime('%Y-%m-%d %H:%M')} by tools/surf_metrics_collect.py",
        "freesurfer tree : /rds/project/rds-CeXlNYOYMxw/derivatives/freesurfer (DAIRC FreeSurfer 7.1.1)",
        "HCP annotation  : <fs>/<sub>/<ses>/label/?h.HCP.fsaverage.aparc.annot (rr480, parcellate.sh), "
        "or legacy/hpc/work/parcellations_backfill/T1/<sub>/<ses>/HCP.fsaverage.aparc/ for backfilled "
        "sessions -- same resolution as abcd.hcp_stats / mr_y_smri__thk__hcp.tsv",
        "DK annotation   : <fs>/<sub>/<ses>/label/?h.aparc.annot",
        "measures        :",
        *[f"  {k}: {v.kind}, {v.source}" + (f", depth {v.depth}" if v.depth else "")
          for k, v in MEASURES.items() if k in ok],
        "depths          : gm50 = white + 0.5*(pial - white); wmD = white - D mm * outward unit "
        "normal of the white surface (area-weighted face normals, oriented white->pial)",
        "intensity       : T1 = mri/nu.mgz (N4/nu-corrected, FreeSurfer-rescaled; conformed 1 mm), "
        "trilinear interpolation (scipy.ndimage.map_coordinates order 1), tkr-RAS -> voxel via inv(vox2ras-tkr)",
        "T2w             : ABCC fMRIPrep v25.1.4 desc-preproc_T2w (in fMRIPrep's T1w space; raw scanner "
        "units, NOT bias-corrected by fMRIPrep). Mapped to DAIRC orig.mgz by mri_coreg (--sep 4 --ftol 1e-4 "
        "--linmintol 1e-2) + bbregister --t1 of desc-preproc_T1w (tools/t2reg.sbatch; "
        "legacy/hpc/work/d1_t2reg/<sub>/<ses>/bbr.lta); divided by an ANTs N4 field (negatives clipped "
        "to 1, fMRIPrep brain mask, -s 4 -b [200]; kept at 8 mm, normalised to median 1 in the mask). "
        "T2 level differs by scanner -- not comparable across manufacturers without harmonisation",
        "T1w/T2w         : fMRIPrep desc-preproc_T1w (N4) / N4-corrected desc-preproc_T2w, both "
        "sampled at gm50 through the same bbr.lta; NaN where T2w <= 0",
        "GWC             : surf/?h.w-g.pct.mgh as written by recon-all pctsurfcon (FreeSurfer 7.1.1)",
        "summary         : parcel = vertex-area-weighted mean (surf/?h.area) over finite samples; "
        "hemisphere/whole cortex = area-weighted over all non-medial-wall vertices "
        "(HCP '???', DK unknown/corpuscallosum, label -1 excluded)",
        f"sessions ok     : {ok}",
        f"software        : python {sys.version.split()[0]}, nibabel {nibabel.__version__}, "
        f"numpy {numpy.__version__}, scipy {scipy.__version__}, pandas {pd.__version__}",
        "files           :",
        *[f"  {p}  {s[0]} x {s[1]}" for p, s in written],
    ]
    (hcp_dir / "PROVENANCE_surf_metrics.txt").write_text("\n".join(lines) + "\n")
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
