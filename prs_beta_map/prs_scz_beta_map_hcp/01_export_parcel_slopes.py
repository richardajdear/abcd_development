"""Step 1 (laptop): the per-child, per-parcel thinning slopes, wide, for CSD3.

Reads the settled HCP-MMP run (thickness_hcp_70_aa6e91efba82, 8,716 children,
358 parcels; the medial-wall parcel H is excluded by the config) and writes
one row per child with the 358 parcel slopes in native units (mm/yr):

    work/parcel_slopes_hcp70.tsv.gz      IID <lh_1> <lh_10d> ... <rh_v23ab>

The slope is the same quantity the pipeline's `global_slope` averages
("per-region" construction: fixed age effect + the child's random slope for
that parcel), so the per-parcel betas below are the decomposition of the
tabled whole-cortex beta.  IID is the export spelling (NDAR_INV...), and the
CSD3 side re-joins on the 8-char NDAR token exactly as setup/align_export.py
does, so the spelling does not matter.

Usage (repo root):  python prs_beta_map/prs_scz_beta_map_hcp/01_export_parcel_slopes.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RUN = REPO / "out" / "thickness_hcp_70_aa6e91efba82"
OUT = HERE / "work" / "parcel_slopes_hcp70.tsv.gz"


def main() -> int:
    p = pd.read_parquet(RUN / "phenotypes" / "phenotypes.parquet",
                        columns=["subject", "label", "phenotype", "value"])
    s = p[p.phenotype == "slope"]
    wide = s.pivot(index="subject", columns="label", values="value")
    assert wide.shape[1] == 358, wide.shape
    # sub-NDARINV003RTV85 -> NDAR_INV003RTV85 (the gcta export spelling)
    wide.index = "NDAR_" + wide.index.str.replace(r"^sub-NDAR", "", regex=True)
    wide.index.name = "IID"
    OUT.parent.mkdir(exist_ok=True)
    wide.to_csv(OUT, sep="\t", float_format="%.7g")
    # the whole-cortex mean must be the pipeline's global_slope; check against
    # the run's own phenotype table for the same children
    g = p[(p.phenotype == "slope")].groupby("subject").value.mean()
    print(f"wrote {OUT.relative_to(REPO)}: {wide.shape[0]} children x {wide.shape[1]} parcels; "
          f"whole-cortex mean of the rows = {g.mean():.6f} mm/yr (tabled HCP global_slope "
          f"mean -0.019930)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
