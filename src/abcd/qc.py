"""
Quality control as composable, logged predicates.

The 5.1 pipeline applied QC by inner-joining a static file
(``SUBJECTS_MRI_ONLY_post_EULER_no_phillips.csv``) whose construction is not
recorded anywhere in the repository.  Every exclusion was therefore invisible:
you could see that 20,197 scans survived, but not how many were dropped by
which criterion, and the criteria could not be varied without regenerating a
file by unknown means.

Here each criterion is a :class:`Predicate` that returns a boolean mask over
scan rows and records how many it removed.  Applying a stack yields a
:class:`QCResult` carrying the surviving rows *and* the attrition ledger, which
is what the sample-flow figure and the run manifest are built from.

Unavailable criteria are explicit
---------------------------------
A predicate whose source data is missing reports ``status='unavailable'`` and
passes every row through.  It does **not** silently succeed: the ledger shows
it was skipped, and :meth:`QCResult.unavailable` lists them so a caller can
refuse to proceed.  This matters for release portability -- the local 5.1 copy
has no ``mri_y_qc_incl`` table and no Euler numbers, so the release-QC and
Euler-threshold predicates are inert here but should activate on 7.0 without
a code change.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Literal

import numpy as np
import pandas as pd

from . import paths

Status = Literal["applied", "unavailable"]


@dataclass
class Predicate:
    """One QC criterion.

    ``fn`` maps the scan table to a boolean Series (True = keep).  It may raise
    :class:`SourceUnavailable`, in which case the predicate is recorded as
    unavailable and no rows are dropped.
    """

    name: str
    description: str
    fn: Callable[[pd.DataFrame], pd.Series]

    def __call__(self, df: pd.DataFrame) -> tuple[pd.Series, Status, str]:
        try:
            mask = self.fn(df)
        except SourceUnavailable as exc:
            return pd.Series(True, index=df.index), "unavailable", str(exc)
        mask = mask.reindex(df.index).fillna(False).astype(bool)
        return mask, "applied", ""


class SourceUnavailable(RuntimeError):
    """Raised by a predicate whose input data is not present in this release."""


@dataclass
class QCResult:
    """Surviving scans plus the attrition ledger that produced them."""

    scans: pd.DataFrame
    ledger: pd.DataFrame
    policy: str
    _initial: int = 0

    @property
    def unavailable(self) -> list[str]:
        return self.ledger.loc[
            self.ledger.status == "unavailable", "criterion"
        ].tolist()

    @property
    def n_scans(self) -> int:
        return len(self.scans)

    @property
    def n_subjects(self) -> int:
        return self.scans.subject.nunique()

    def require_available(self, *names: str) -> None:
        """Raise if any named criterion could not actually be applied."""
        missing = [n for n in names if n in self.unavailable]
        if missing:
            raise SourceUnavailable(
                f"QC criteria {missing} were requested but their source data is "
                f"absent for this release; refusing to proceed silently."
            )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"QCResult(policy={self.policy!r}, {self._initial} -> {self.n_scans} scans, "
            f"{self.n_subjects} subjects, unavailable={self.unavailable})"
        )


# --------------------------------------------------------------------------
# Individual criteria
# --------------------------------------------------------------------------

LEGACY_QC_FILE = "SUBJECTS_MRI_ONLY_post_EULER_no_phillips.csv"

#: The legacy list keys sessions by BIDS-style names.
LEGACY_SESSION_TO_VISIT = {
    "ses-baseline-year1": "v0",
    "ses-2YearFollowUpYArm1": "v2",
    "ses-4YearFollowUpYArm1": "v4",
    "ses-6YearFollowUpYArm1": "v6",
}


def load_legacy_list() -> pd.DataFrame:
    """Read the 5.1 static QC list (subject x session that passed)."""
    p = paths.abcd_root() / LEGACY_QC_FILE
    if not p.exists():
        raise SourceUnavailable(f"legacy QC list not found at {p}")
    q = pd.read_csv(p)
    q = q.rename(columns={"Subject": "subject"})
    q["visit"] = q.session.map(LEGACY_SESSION_TO_VISIT)
    unmapped = q[q.visit.isna()].session.unique()
    if len(unmapped):
        raise ValueError(f"unmapped sessions in legacy QC list: {unmapped}")
    return q[["subject", "visit"]].drop_duplicates()


def legacy_euler_predicate() -> Predicate:
    """Scan is in the pre-computed Euler + non-Philips list.

    Opaque but reproducible: it is the exact criterion behind the published
    5.1 results, retained so the new pipeline can be compared against them.
    """

    def fn(df: pd.DataFrame) -> pd.Series:
        keep = load_legacy_list()
        idx = pd.MultiIndex.from_frame(df[["subject", "visit"]])
        allowed = pd.MultiIndex.from_frame(keep)
        return pd.Series(idx.isin(allowed), index=df.index)

    return Predicate(
        "legacy_euler_list",
        f"scan present in {LEGACY_QC_FILE} (Euler QC + Philips scanners removed)",
        fn,
    )


def complete_regions_predicate(n_expected: int | None = None) -> Predicate:
    """Scan has a non-missing value for every region of the parcellation.

    Applied to the *wide* scan table produced by :func:`scan_table`, where
    ``n_regions_present`` has already been counted.
    """

    def fn(df: pd.DataFrame) -> pd.Series:
        if "n_regions_present" not in df.columns:
            raise SourceUnavailable(
                "scan table lacks n_regions_present; build it with scan_table()"
            )
        target = n_expected if n_expected is not None else df.n_regions_present.max()
        return df.n_regions_present >= target

    return Predicate(
        "complete_regions",
        "no missing regional values in the scan",
        fn,
    )


def release_qc_include_predicate(adapter) -> Predicate:
    """ABCD's own recommended structural-imaging inclusion flag.

    Reads ``mri_y_qc_incl.csv`` (``imgincl_t1w_include``).  Absent from the
    local 5.1 copy; expected to be present in 7.0.
    """

    def fn(df: pd.DataFrame) -> pd.Series:
        p = adapter.imaging_dir / "mri_y_qc_incl.csv"
        if not p.exists():
            raise SourceUnavailable(
                f"{p.name} not present in release {adapter.release} at {p.parent}"
            )
        inc = pd.read_csv(p, low_memory=False)
        col = "imgincl_t1w_include"
        if col not in inc.columns:
            raise SourceUnavailable(f"{p.name} has no column {col}")
        inc["subject"] = adapter._to_bids(inc.src_subject_id)
        inc["visit"] = adapter._map_visits(inc.eventname)
        inc = inc[inc[col] == 1][["subject", "visit"]].drop_duplicates()
        idx = pd.MultiIndex.from_frame(df[["subject", "visit"]])
        return pd.Series(idx.isin(pd.MultiIndex.from_frame(inc)), index=df.index)

    return Predicate(
        "release_qc_include",
        "ABCD imgincl_t1w_include == 1",
        fn,
    )


def euler_threshold_predicate(adapter, min_euler: float = -200.0) -> Predicate:
    """FreeSurfer Euler number above a threshold (surface-defect count).

    Absent from the local 5.1 copy. When 7.0 exposes an Euler column, point
    ``EULER_SOURCES`` at it; the threshold then becomes a config knob rather
    than something baked into a pre-made subject list.
    """

    EULER_SOURCES = [
        ("mri_y_qc_motion.csv", "iqc_t1_euler_total"),
    ]

    def fn(df: pd.DataFrame) -> pd.Series:
        for fname, col in EULER_SOURCES:
            p = adapter.imaging_dir / fname
            if p.exists():
                tab = pd.read_csv(p, low_memory=False)
                if col in tab.columns:
                    tab["subject"] = adapter._to_bids(tab.src_subject_id)
                    tab["visit"] = adapter._map_visits(tab.eventname)
                    m = tab[["subject", "visit", col]].dropna()
                    merged = df[["subject", "visit"]].merge(
                        m, on=["subject", "visit"], how="left"
                    )
                    return pd.Series(
                        (merged[col] >= min_euler).values, index=df.index
                    )
        raise SourceUnavailable(
            f"no Euler-number source found for release {adapter.release} "
            f"(looked for {[f for f, _ in EULER_SOURCES]})"
        )

    return Predicate(
        "euler_threshold",
        f"FreeSurfer Euler number >= {min_euler}",
        fn,
    )


def scanner_predicate(adapter, exclude: tuple[str, ...] = ("Philips",)) -> Predicate:
    """Exclude named scanner manufacturers.

    Philips sites were excluded in the 5.1 work because of a known
    distortion-correction issue affecting early ABCD Philips data.
    """

    def fn(df: pd.DataFrame) -> pd.Series:
        p = adapter.imaging_dir / "mri_y_adm_info.csv"
        if not p.exists():
            raise SourceUnavailable(f"{p.name} not present in release {adapter.release}")
        info = pd.read_csv(p, low_memory=False)
        col = next(
            (c for c in ("mri_info_manufacturer", "mri_info_deviceserialnumber")
             if c in info.columns),
            None,
        )
        if col is None:
            raise SourceUnavailable(f"{p.name} has no manufacturer column")
        info["subject"] = adapter._to_bids(info.src_subject_id)
        info["visit"] = adapter._map_visits(info.eventname)
        merged = df[["subject", "visit"]].merge(
            info[["subject", "visit", col]], on=["subject", "visit"], how="left"
        )
        bad = merged[col].astype(str).str.contains("|".join(exclude), case=False, na=False)
        return pd.Series(~bad.values, index=df.index)

    return Predicate(
        "scanner_manufacturer",
        f"scanner manufacturer not in {exclude}",
        fn,
    )


# --------------------------------------------------------------------------
# Policies
# --------------------------------------------------------------------------

def build_policy(policy: str, adapter, n_expected_regions: int | None = None
                 ) -> list[Predicate]:
    """Return the predicate stack for a named QC policy (see config.QC_POLICIES)."""
    if policy == "none":
        return []
    if policy == "legacy_euler":
        return [complete_regions_predicate(n_expected_regions),
                legacy_euler_predicate()]
    if policy == "coded":
        # Ordered cheapest-first so the ledger reads as a funnel.
        return [
            complete_regions_predicate(n_expected_regions),
            release_qc_include_predicate(adapter),
            euler_threshold_predicate(adapter),
            scanner_predicate(adapter),
        ]
    raise ValueError(f"unknown QC policy {policy!r}")


def scan_table(imaging_long: pd.DataFrame) -> pd.DataFrame:
    """Collapse long imaging data to one row per scan, with completeness counts.

    QC operates on scans, not on region-level rows, so this is the unit the
    predicates see.
    """
    regional = imaging_long[~imaging_long.is_global]
    grouped = regional.groupby(["subject", "visit"], observed=True)
    return (
        grouped.agg(
            n_regions_present=("value", "count"),
            n_regions_total=("value", "size"),
        )
        .reset_index()
    )


def apply_qc(scans: pd.DataFrame, predicates: list[Predicate],
             policy: str = "custom") -> QCResult:
    """Apply a predicate stack sequentially, recording attrition at each step."""
    initial = len(scans)
    rows = [{
        "step": 0, "criterion": "initial", "description": "all scans with imaging data",
        "status": "applied", "n_scans_in": initial, "n_scans_out": initial,
        "n_dropped": 0, "n_subjects_out": scans.subject.nunique(), "detail": "",
    }]
    current = scans
    for i, pred in enumerate(predicates, start=1):
        mask, status, detail = pred(current)
        kept = current[mask]
        rows.append({
            "step": i,
            "criterion": pred.name,
            "description": pred.description,
            "status": status,
            "n_scans_in": len(current),
            "n_scans_out": len(kept),
            "n_dropped": len(current) - len(kept),
            "n_subjects_out": kept.subject.nunique(),
            "detail": detail,
        })
        current = kept
    return QCResult(scans=current.reset_index(drop=True),
                    ledger=pd.DataFrame(rows), policy=policy, _initial=initial)


# --------------------------------------------------------------------------
# Sample definition (which *subjects*, given which scans passed)
# --------------------------------------------------------------------------

def visit_pattern(scans: pd.DataFrame, visit_order: tuple[str, ...]) -> pd.Series:
    """Per subject, a '+'-joined string of visits present, in temporal order."""
    rank = {v: i for i, v in enumerate(visit_order)}
    return (
        scans.groupby("subject").visit
        .apply(lambda s: "+".join(sorted(set(s), key=lambda v: rank.get(v, 99))))
    )


def select_subjects(scans: pd.DataFrame, rule: str, min_visits: int,
                    visit_order: tuple[str, ...]) -> tuple[pd.Index, dict]:
    """Apply a sample rule; returns (subjects, info dict for the ledger)."""
    counts = scans.groupby("subject").visit.nunique()
    if rule == "all":
        keep = counts.index
    elif rule == "min_visits":
        keep = counts[counts >= min_visits].index
    elif rule == "has_last":
        last = visit_order[-1]
        keep = scans.loc[scans.visit == last, "subject"].unique()
        keep = pd.Index(keep)
    else:
        raise ValueError(f"unknown sample rule {rule!r}")
    info = {
        "rule": rule,
        "min_visits": min_visits if rule == "min_visits" else None,
        "n_subjects_before": len(counts),
        "n_subjects_after": len(keep),
    }
    return pd.Index(keep), info
