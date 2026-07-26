"""Shared fixtures.

Fixtures that need a fitted run skip rather than fail when the run directory
is absent, so the suite is still meaningful on a fresh clone that has not run
the pipeline yet.
"""
from __future__ import annotations

import pandas as pd
import pytest

from abcd import paths, spatial

HCP_RUN = "thickness_hcp_51_f408a620fde4"


@pytest.fixture(scope="session")
def hcp_geom():
    p = paths.DATA_DIR / "hcp_centroids.csv"
    if not p.exists():
        pytest.skip("hcp_centroids.csv not present")
    return spatial.load_dk_geometry(p)


@pytest.fixture(scope="session")
def hcp_geom_lh(hcp_geom):
    """Left-hemisphere-only geometry with hemi prefixes stripped.

    AHBA expression is left-hemisphere only (donor sampling), so gene-level
    tests spin within one hemisphere and match on bare region names.
    """
    lh = [lab for lab in hcp_geom.labels if lab.startswith("lh_")]
    idx = [hcp_geom.labels.index(lab) for lab in lh]
    return spatial.ParcelGeometry(
        labels=tuple(lab[3:] for lab in lh),
        coords=hcp_geom.coords[idx],
        space="sphere",
        hemi=tuple(["lh"] * len(lh)),
    )


@pytest.fixture(scope="session")
def hcp_fits():
    p = paths.OUT_DIR / HCP_RUN / "fits" / "fixed.parquet"
    if not p.exists():
        pytest.skip(f"fitted run {HCP_RUN} not present; run the pipeline first")
    return pd.read_parquet(p)


@pytest.fixture(scope="session")
def dev_map_hcp(hcp_fits):
    """Group-mean absolute thickness change per region (no global covariate)."""
    return hcp_fits[hcp_fits.term == "age_c"].set_index("label").estimate
