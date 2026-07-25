"""
Release-agnostic access to ABCD tabulated data.

The pipeline never opens a release CSV directly.  It asks a
:class:`ReleaseAdapter` for a *concept* -- imaging, demographics, longitudinal
tracking -- and the adapter knows which file and which column names that
release uses.  Swapping 5.1 for 7.0 is then a one-line config change rather
than a search-and-replace through the analysis.

Two release-specific irregularities the 5.1 adapter absorbs:

1. **Region abbreviations differ between metric tables.**  ABCD names the same
   DK region differently in the T2 tables than in the thickness/area/volume
   tables -- ``smri_thick_cdk_cdacatelh`` but ``smri_t2wg02_cdk_cdatcgatelh``
   for left caudal anterior cingulate; 19 of 34 regions per hemisphere differ.
   ``data/region_labels.csv`` carries both spellings (``stem_a``/``stem_b``)
   keyed to one canonical label, verified to resolve for all 71 columns of all
   10 Desikan tables.

2. **Global columns are named by aggregation type.**  Area and volume use
   ``total*``; every other metric uses ``mean*``.

All returned frames use the canonical schema: ``subject`` (BIDS-style
``sub-NDARINV...``), ``visit`` (abstract ``v0``/``v2``/``v4``/``v6``, not the
release's event strings), and long-format ``region``/``hemi``/``value``.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from . import paths

# --------------------------------------------------------------------------
# Canonical region labels (shared by Python and R; R reads the same CSV)
# --------------------------------------------------------------------------

REGION_LABEL_FILE = paths.DATA_DIR / "region_labels.csv"


@lru_cache(maxsize=None)
def region_labels(parcellation: str = "dsk") -> pd.DataFrame:
    """Canonical region table: one row per (hemi, region) plus global rows."""
    df = pd.read_csv(REGION_LABEL_FILE)
    out = df[df.parcellation == parcellation].reset_index(drop=True)
    if out.empty:
        raise ValueError(
            f"no region labels for parcellation={parcellation!r} in "
            f"{REGION_LABEL_FILE}; available: {sorted(df.parcellation.unique())}"
        )
    return out


class ReleaseAdapter(ABC):
    """Interface every ABCD release must satisfy."""

    release: str
    #: Abstract visit code -> release event string.
    VISITS: dict[str, str]

    @property
    def root(self) -> Path:
        return paths.release_dir(self.release)

    @abstractmethod
    def imaging(self, metric: str, parcellation: str = "dsk") -> pd.DataFrame:
        """Long format: subject, visit, hemi, region, value, is_global."""

    @abstractmethod
    def demographics(self) -> pd.DataFrame:
        """One row per subject: subject, sex."""

    @abstractmethod
    def longitudinal(self) -> pd.DataFrame:
        """One row per (subject, visit): age, site, family_id."""

    @abstractmethod
    def genetics(self) -> dict[str, Path]:
        """Named paths to genetic resources (relatedness, genotypes, ...)."""

    # -- shared helpers -------------------------------------------------
    @staticmethod
    def _to_bids(subject: pd.Series) -> pd.Series:
        """``NDAR_INV00CY2MDM`` -> ``sub-NDARINV00CY2MDM`` (idempotent)."""
        return subject.str.replace("NDAR_", "sub-NDAR", regex=False)

    def _map_visits(self, eventname: pd.Series) -> pd.Series:
        inverse = {v: k for k, v in self.VISITS.items()}
        return eventname.map(inverse)


class Release51Adapter(ReleaseAdapter):
    """ABCD data release 5.1 (structural imaging only in the local copy)."""

    release = "5.1"

    VISITS = {
        "v0": "baseline_year_1_arm_1",
        "v2": "2_year_follow_up_y_arm_1",
        "v4": "4_year_follow_up_y_arm_1",
    }

    #: metric -> (table stem, column prefix, stem family, global aggregation)
    METRIC_TABLES = {
        "thickness":   ("thk",      "smri_thick_cdk_",      "a", "mean"),
        "area":        ("area",     "smri_area_cdk_",       "a", "total"),
        "volume":      ("vol",      "smri_vol_cdk_",        "a", "total"),
        "sulc":        ("sulc",     "smri_sulc_cdk_",       "a", "mean"),
        "t1_gray":     ("t1_gray",  "smri_t1wgray02_cdk_",  "a", "mean"),
        "t1_white":    ("t1_white", "smri_t1ww02_cdk_",     "a", "mean"),
        "t1_contrast": ("t1_contr", "smri_t1wcnt_cdk_",     "a", "mean"),
        "t2_gray":     ("t2_gray",  "smri_t2wg02_cdk_",     "b", "mean"),
        "t2_white":    ("t2_white", "smri_t2ww02_cdk_",     "b", "mean"),
        "t2_contrast": ("t2_contr", "smri_t2wcnt_cdk_",     "b", "mean"),
    }

    @property
    def imaging_dir(self) -> Path:
        return self.root / "core" / "imaging"

    @property
    def general_dir(self) -> Path:
        return self.root / "core" / "abcd-general"

    # ------------------------------------------------------------------
    def imaging(self, metric: str, parcellation: str = "dsk") -> pd.DataFrame:
        if metric not in self.METRIC_TABLES:
            raise KeyError(
                f"metric {metric!r} is not a raw 5.1 table. Available: "
                f"{sorted(self.METRIC_TABLES)}. Derived metrics such as "
                "'t1t2_ratio' are built in assemble.py, not here."
            )
        if parcellation != "dsk":
            raise NotImplementedError(
                f"parcellation {parcellation!r} not yet wired; only 'dsk' has "
                "a verified label mapping in data/region_labels.csv"
            )
        table, prefix, family, agg = self.METRIC_TABLES[metric]
        labels = region_labels(parcellation).copy()
        stem_col = f"stem_{family}"

        # globals are spelled total*/mean* depending on the metric
        labels["stem"] = labels[stem_col]
        if agg == "total":
            labels.loc[labels.is_global, "stem"] = labels.loc[
                labels.is_global, "stem"
            ].str.replace("mean", "total", regex=False)
        labels["column"] = prefix + labels["stem"]

        path = self.imaging_dir / f"mri_y_smr_{table}_dsk.csv"
        raw = pd.read_csv(path, low_memory=False)
        missing = sorted(set(labels.column) - set(raw.columns))
        if missing:
            raise KeyError(f"{path.name} is missing expected columns: {missing}")

        keep = ["src_subject_id", "eventname"] + labels.column.tolist()
        raw = raw[keep]
        long = raw.melt(
            id_vars=["src_subject_id", "eventname"],
            var_name="column",
            value_name="value",
        )
        long = long.merge(
            labels[["column", "hemi", "region", "label", "is_global"]],
            on="column",
            how="left",
        )
        long["subject"] = self._to_bids(long.src_subject_id)
        long["visit"] = self._map_visits(long.eventname)
        long["metric"] = metric
        long = long.dropna(subset=["visit"])
        return long[
            ["subject", "visit", "metric", "hemi", "region", "label",
             "is_global", "value"]
        ].reset_index(drop=True)

    # ------------------------------------------------------------------
    def demographics(self) -> pd.DataFrame:
        """Sex from the baseline parent demographics form.

        ``demo_sex_v2``: 1=male, 2=female, 3=intersex-male.  The 5.1 code
        dropped 3 outright; we keep the raw code alongside a two-level factor
        so the exclusion is visible downstream rather than silent.
        """
        path = self.general_dir / "abcd_p_demo.csv"
        demo = pd.read_csv(
            path,
            usecols=["src_subject_id", "eventname", "demo_sex_v2"],
            low_memory=False,
        )
        demo = demo[demo.eventname == self.VISITS["v0"]]
        demo["subject"] = self._to_bids(demo.src_subject_id)
        demo["sex_code"] = demo.demo_sex_v2
        demo["sex"] = demo.demo_sex_v2.map({1: "M", 2: "F"})
        return (
            demo[["subject", "sex", "sex_code"]]
            .dropna(subset=["subject"])
            .drop_duplicates("subject")
            .reset_index(drop=True)
        )

    # ------------------------------------------------------------------
    def longitudinal(self) -> pd.DataFrame:
        """Age, site and family from the longitudinal tracking table."""
        path = self.general_dir / "abcd_y_lt.csv"
        lt = pd.read_csv(
            path,
            usecols=[
                "src_subject_id", "eventname", "interview_age",
                "site_id_l", "rel_family_id",
            ],
            low_memory=False,
        )
        lt["subject"] = self._to_bids(lt.src_subject_id)
        lt["visit"] = self._map_visits(lt.eventname)
        lt = lt.dropna(subset=["visit", "interview_age"])
        lt["age"] = lt.interview_age / 12.0
        lt = lt.rename(columns={"site_id_l": "site", "rel_family_id": "family_id"})

        # family_id is recorded per visit but is a subject-level attribute;
        # carry the baseline value forward so a subject cannot change family
        # mid-study (which would silently split a family random effect).
        base = (
            lt[lt.visit == "v0"][["subject", "family_id"]]
            .dropna()
            .drop_duplicates("subject")
            .rename(columns={"family_id": "family_id_baseline"})
        )
        lt = lt.merge(base, on="subject", how="left")
        lt["family_id"] = lt.family_id_baseline.fillna(lt.family_id)
        return lt[
            ["subject", "visit", "age", "site", "family_id"]
        ].reset_index(drop=True)

    # ------------------------------------------------------------------
    def genetics(self) -> dict[str, Path]:
        d = self.root / "core" / "genetics"
        return {
            "pihat": d / "gen_y_pihat.csv",
            "zygosity": d / "gen_y_zygrat.csv",
            # NB: the local 5.1 copy contains no genotype/WGS data. GCTA and
            # GWAS steps therefore cannot run against 5.1 here; see hpc/README.
        }

    def relatedness(self) -> pd.DataFrame:
        """One row per subject: family, sibling group, and genetic relatedness.

        Subject ids are normalised to the same BIDS-style form the imaging
        tables use, so this joins directly to phenotype tables.  The release
        ships raw ids as ``NDAR_INV...`` here but ``sub-NDARINV...`` elsewhere;
        joining without normalising silently yields zero matched rows, which
        is exactly the kind of failure that looks like a null result.

        ``pair_type`` classifies subjects by genetic relatedness where it is
        available: MZ twins (pi-hat > 0.9), DZ twins or full siblings
        (0.35-0.65), and everything else.  ~3,670 of 11,868 subjects have a
        pi-hat value in 5.1; the rest are singletons or ungenotyped.
        """
        p = pd.read_csv(self.genetics()["pihat"], low_memory=False)
        p = p[p.eventname == self.VISITS["v0"]].copy()
        out = pd.DataFrame({
            "subject": self._to_bids(p.src_subject_id),
            "family_id": p.rel_family_id,
            "group_id": p.rel_group_id,
            "relationship": p.rel_relationship,
            "pi_hat": p.genetic_pi_hat_1,
        })
        out["pair_type"] = np.select(
            [out.pi_hat > 0.9, out.pi_hat.between(0.35, 0.65)],
            ["MZ", "DZ_or_sib"], default="other",
        )
        out.loc[out.pi_hat.isna(), "pair_type"] = "unknown"
        return out.reset_index(drop=True)


class Release70Adapter(Release51Adapter):
    """ABCD data release 7.0 -- **stub, not yet verified against real data**.

    Inherits 5.1 behaviour so the pipeline is runnable the moment the release
    is linked, but every assumption below must be checked before any 7.0
    result is trusted.  Run ``python -m abcd.io --verify 7.0`` (see
    :func:`verify_adapter`) once ``abcd-data-release-7.0/`` exists.

    TO VERIFY ON LINKING
    --------------------
    1. **Event names.** ``VISITS`` below assumes the 6-year visit is
       ``6_year_follow_up_y_arm_1``.  Confirm against the release notes; also
       confirm the baseline/2y/4y strings are unchanged.
    2. **Column naming.** Confirm the ``smri_*_cdk_*`` prefixes and the
       stem_a/stem_b split still hold. If ABCD harmonised the T2 spellings,
       ``METRIC_TABLES`` families need updating (the label CSV already carries
       both, so this is a one-character change per metric).
    3. **New modalities.** The local 5.1 copy is structural-only. 7.0 is
       expected to carry resting-state and diffusion tables; add them to
       ``METRIC_TABLES`` with their own prefixes and, for connectivity
       matrices, a separate reader -- they are edge-level, not region-level,
       so they need a different long schema (``region_i``/``region_j``).
    4. **Genotypes.** Populate :meth:`genetics` with the genotype/WGS paths.
       This is the blocker for every GCTA/GWAS step.
    5. **QC variables.** Confirm the FreeSurfer QC and Euler columns used by
       ``qc.py`` still exist under the same names.
    """

    release = "7.0"

    VISITS = {
        "v0": "baseline_year_1_arm_1",
        "v2": "2_year_follow_up_y_arm_1",
        "v4": "4_year_follow_up_y_arm_1",
        "v6": "6_year_follow_up_y_arm_1",  # ASSUMED - verify
    }

    VERIFIED = False


ADAPTERS: dict[str, type[ReleaseAdapter]] = {
    "5.1": Release51Adapter,
    "7.0": Release70Adapter,
}


def get_adapter(release: str) -> ReleaseAdapter:
    """Instantiate the adapter for a release, warning if it is unverified."""
    if release not in ADAPTERS:
        raise KeyError(f"no adapter for release {release!r}; have {sorted(ADAPTERS)}")
    adapter = ADAPTERS[release]()
    if not getattr(adapter, "VERIFIED", True):
        import warnings

        warnings.warn(
            f"Release {release} adapter is an unverified stub. Its visit names, "
            "column prefixes and genetics paths are assumptions inherited from "
            "5.1. Run verify_adapter() before trusting any result.",
            stacklevel=2,
        )
    return adapter


def verify_adapter(release: str, parcellation: str = "dsk") -> pd.DataFrame:
    """Smoke-test an adapter against real files on disk.

    Returns a check table (one row per assertion) rather than raising, so a
    partial pass is visible.  This is the first thing to run when 7.0 lands.
    """
    adapter = ADAPTERS[release]()
    checks: list[dict] = []

    def record(name, ok, detail=""):
        checks.append({"check": name, "ok": bool(ok), "detail": str(detail)})

    try:
        root = adapter.root
        record("release_dir_exists", root.exists(), root)
    except Exception as exc:  # noqa: BLE001
        record("release_dir_exists", False, exc)
        return pd.DataFrame(checks)

    for metric in adapter.METRIC_TABLES:
        table = adapter.METRIC_TABLES[metric][0]
        p = adapter.imaging_dir / f"mri_y_smr_{table}_{parcellation}.csv"
        record(f"table_exists:{metric}", p.exists(), p.name)

    for metric in adapter.METRIC_TABLES:
        try:
            df = adapter.imaging(metric, parcellation)
            n_reg = df[~df.is_global].region.nunique()
            record(f"imaging_readable:{metric}", n_reg == 34,
                   f"{n_reg} regions, {df.visit.nunique()} visits")
        except Exception as exc:  # noqa: BLE001
            record(f"imaging_readable:{metric}", False, exc)

    for name, fn in [("demographics", adapter.demographics),
                     ("longitudinal", adapter.longitudinal)]:
        try:
            df = fn()
            record(f"{name}_readable", len(df) > 0, f"{len(df)} rows")
        except Exception as exc:  # noqa: BLE001
            record(f"{name}_readable", False, exc)

    # visit coverage: does every declared visit actually appear?
    try:
        lt = adapter.longitudinal()
        seen = set(lt.visit.unique())
        for v in adapter.VISITS:
            record(f"visit_present:{v}", v in seen,
                   f"{(lt.visit == v).sum()} rows")
    except Exception as exc:  # noqa: BLE001
        record("visit_coverage", False, exc)

    for name, p in adapter.genetics().items():
        record(f"genetics:{name}", p.exists(), p)

    return pd.DataFrame(checks)


if __name__ == "__main__":  # pragma: no cover
    import argparse

    ap = argparse.ArgumentParser(description="Verify an ABCD release adapter.")
    ap.add_argument("--verify", default="5.1", help="release, e.g. 5.1 or 7.0")
    args = ap.parse_args()
    result = verify_adapter(args.verify)
    print(result.to_string(index=False))
    print(f"\n{result.ok.sum()}/{len(result)} checks passed")
