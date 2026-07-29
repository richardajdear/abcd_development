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


#: Re-exported from :mod:`paths` so both adapters and predicates can raise it
#: without a circular import.  Kept as a module-level name here because
#: existing code and tests import it from ``qc``.
SourceUnavailable = paths.SourceUnavailable


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
    # Searched across every root, not just abcd_root(): this file sits beside
    # the 5.1 release under ~/Git/ABCD, while the preferred root is now the repo
    # (which holds 7.0). A single-root lookup silently reports the legacy policy
    # as unavailable and lets QC pass every scan through.
    p = paths.find_in_roots(LEGACY_QC_FILE)
    if p is None:
        raise SourceUnavailable(
            f"legacy QC list {LEGACY_QC_FILE} not found in any of "
            f"{[str(r) for r in paths.abcd_roots()]}"
        )
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

    Delegates to ``adapter.qc_include()``, so the table name and column live
    with the release adapter rather than here.  Absent from the local 5.1 copy;
    present in 7.0 as ``mr_y_qc__incl__smri__t1_indicator``.
    """

    def fn(df: pd.DataFrame) -> pd.Series:
        inc = adapter.qc_include()  # raises SourceUnavailable if absent
        inc = inc[inc.include == 1][["subject", "visit"]].drop_duplicates()
        idx = pd.MultiIndex.from_frame(df[["subject", "visit"]])
        return pd.Series(idx.isin(pd.MultiIndex.from_frame(inc)), index=df.index)

    return Predicate(
        "release_qc_include",
        "ABCD imgincl_t1w_include == 1",
        fn,
    )


def surface_defect_predicate(adapter, max_defects: float = 46.0) -> Predicate:
    """Exclude scans with too many FreeSurfer topological defects.

    This replaces the Euler-number threshold used in the 5.1 work.  The Euler
    characteristic of a reconstructed surface is a linear function of its
    defect count, so the two carry the same information with opposite sign:
    ``euler = 2 - 2 * defects`` per hemisphere, summed over hemispheres.  The
    5.1 threshold of ``euler >= -200`` therefore corresponds to roughly 50
    defects; the default here (46) is the value that reproduces the legacy
    exclusion rate on the subjects the two releases share.

    Validated against the legacy list: subjects it excluded have a median 32
    defects against 19 for those retained (Mann-Whitney p = 3e-83, AUC 0.73),
    so the substitution preserves the QC signal rather than merely being
    available.
    """

    def fn(df: pd.DataFrame) -> pd.Series:
        tab = adapter.surface_defects()  # raises SourceUnavailable if absent
        merged = df[["subject", "visit"]].merge(
            tab.drop_duplicates(["subject", "visit"]),
            on=["subject", "visit"], how="left",
        )
        # Missing defect count is not evidence of a good surface; keep the scan
        # only if we have a measurement and it passes.
        keep = merged.defects.notna() & (merged.defects <= max_defects)
        return pd.Series(keep.values, index=df.index)

    return Predicate(
        "surface_defects",
        f"FreeSurfer topological defect count <= {max_defects:g} "
        f"(Euler-number equivalent)",
        fn,
    )


#: Retained under its historical name so existing configs and the 5.1
#: replication path keep working; the implementation is defect-count based.
euler_threshold_predicate = surface_defect_predicate


def scanner_predicate(adapter, exclude: tuple[str, ...] = ("Philips",)) -> Predicate:
    """Exclude named scanner manufacturers.

    Philips sites were excluded in the 5.1 work because of a known
    distortion-correction issue affecting early ABCD Philips data.  7.0 spells
    the vendor differently again (``SIEMENS`` and ``Siemens Healthineers`` are
    both present), so matching is case-insensitive substring rather than exact.
    """

    def fn(df: pd.DataFrame) -> pd.Series:
        tab = adapter.scanner()  # raises SourceUnavailable if absent
        merged = df[["subject", "visit"]].merge(
            tab.drop_duplicates(["subject", "visit"]),
            on=["subject", "visit"], how="left",
        )
        bad = merged.manufacturer.astype(str).str.contains(
            "|".join(exclude), case=False, na=False
        )
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
            surface_defect_predicate(adapter),
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
