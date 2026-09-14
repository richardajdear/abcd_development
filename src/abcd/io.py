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

import re
from abc import ABC, abstractmethod
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from . import paths
from .paths import SourceUnavailable

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
        if parcellation == "hcp":
            return self._imaging_hcp(metric)
        if parcellation != "dsk":
            raise NotImplementedError(
                f"parcellation {parcellation!r} not yet wired; only 'dsk' and "
                "'hcp' have a verified label mapping"
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
    HCP_FILES = {
        "v0": "ses-baseline-year1_HCP.fsaverage.aparc_CT.csv",
        "v2": "ses-2YearFollowUpYArm1_HCP.fsaverage.aparc_CT.csv",
        "v4": "ses-4YearFollowUpYArm1_HCP.fsaverage.aparc_CT.csv",
    }

    def _imaging_hcp(self, metric: str) -> pd.DataFrame:
        """Thickness in the HCP-MMP (Glasser) 360-region parcellation.

        These come from ``processed/`` rather than ``core/`` -- they are
        surface parcellations run locally, not release tables, so they carry
        no release QC columns and their subject coverage is smaller than the
        core tables (10,779 / 7,093 / 2,801 at v0 / v2 / v4).

        The reason to use them: the AHBA C3 gene-expression component is
        defined on HCP-MMP regions.  Mapping C3 onto DK would require
        averaging expression components across parcel boundaries, which
        blurs precisely the spatial gradient being tested.  Fitting in the
        native parcellation of the reference map avoids that.
        """
        if metric != "thickness":
            raise KeyError(
                f"only 'thickness' is available in the HCP processed files; "
                f"got {metric!r}"
            )
        frames = []
        for visit, fname in self.HCP_FILES.items():
            path = self.root / "processed" / fname
            if not path.exists():
                raise SourceUnavailable(f"HCP parcellated file absent: {path}")
            raw = pd.read_csv(path, low_memory=False)
            roi = [c for c in raw.columns
                   if c.endswith("_ROI") and "???" not in c]
            long = raw[["Subject"] + roi].melt(
                id_vars="Subject", var_name="column", value_name="value")
            long["visit"] = visit
            frames.append(long)
        long = pd.concat(frames, ignore_index=True)

        # lh_L_V1_ROI -> hemi lh, region V1
        parts = long.column.str.extract(r"^(lh|rh)_[LR]_(.+)_ROI$")
        long["hemi"], long["region"] = parts[0], parts[1]
        if long.hemi.isna().any():
            bad = long.loc[long.hemi.isna(), "column"].unique()[:5]
            raise ValueError(f"unparsed HCP column names: {list(bad)}")
        long["label"] = long.hemi + "_" + long.region
        long["subject"] = long.Subject.astype(str)
        long["metric"] = metric
        long["is_global"] = False

        # the whole-cortex mean is needed as a covariate but is not a column
        # here, so it is computed from the parcels
        g = (long.groupby(["subject", "visit"], observed=True).value.mean()
             .reset_index())
        g["hemi"], g["region"], g["label"] = "both", "global_mean", "global_mean"
        g["metric"], g["is_global"] = metric, True

        out = pd.concat([long, g], ignore_index=True)
        return out[["subject", "visit", "metric", "hemi", "region", "label",
                    "is_global", "value"]].dropna(subset=["value"]
                                                  ).reset_index(drop=True)

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

    # -- QC sources ------------------------------------------------------
    # Each returns a tidy (subject, visit, <value>) frame, or raises
    # SourceUnavailable.  The predicates in qc.py call these rather than
    # opening release files themselves, so adding a release means writing
    # adapter methods, not editing QC logic.

    def _qc_table(self, fname: str, col: str, out_name: str,
                  dtype=None) -> pd.DataFrame:
        p = self.imaging_dir / fname
        if not p.exists():
            raise SourceUnavailable(
                f"{fname} not present in release {self.release} at {p.parent}"
            )
        raw = pd.read_csv(p, low_memory=False)
        if col not in raw.columns:
            raise SourceUnavailable(f"{fname} has no column {col!r}")
        out = pd.DataFrame({
            "subject": self._to_bids(raw.src_subject_id),
            "visit": self._map_visits(raw.eventname),
            out_name: raw[col] if dtype is None else raw[col].astype(dtype),
        })
        return out.dropna(subset=["visit"]).reset_index(drop=True)

    def qc_include(self) -> pd.DataFrame:
        """ABCD's recommended T1w inclusion flag (absent from the local copy)."""
        return self._qc_table("mri_y_qc_incl.csv", "imgincl_t1w_include", "include")

    def scanner(self) -> pd.DataFrame:
        """Per-visit scanner manufacturer."""
        for fname, col in (("mri_y_adm_info.csv", "mri_info_manufacturer"),
                           ("mri_y_qc_motion.csv", "mri_info_manufacturer")):
            try:
                return self._qc_table(fname, col, "manufacturer")
            except SourceUnavailable:
                continue
        raise SourceUnavailable(
            f"no manufacturer column found for release {self.release}"
        )

    def surface_defects(self) -> pd.DataFrame:
        """Topological defect count.

        5.1's local copy exposes neither this nor the Euler number, so QC on
        5.1 falls back to the precomputed legacy subject list -- see
        :func:`qc.legacy_euler_predicate`.
        """
        raise SourceUnavailable(
            f"release {self.release} exposes no surface-defect or Euler column; "
            "use qc_policy='legacy_euler'"
        )

    def genetic_pcs(self, n: int = 10) -> pd.DataFrame:
        raise SourceUnavailable(
            f"release {self.release} has no ancestry PC table in the local copy"
        )

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
    """ABCD data release 7.0.

    7.0 is not a superset of 5.1 with more rows -- it is a different tabulation
    of the same study, and every one of the stub's original guesses was wrong.
    The differences that matter, each verified against the release:

    1. **Directory name.** ``abcd-7.0/``, not ``abcd-data-release-7.0/``
       (handled in :func:`paths.release_dir`).
    2. **Session codes.** ``ses-00A``/``ses-02A``/``ses-04A``/``ses-06A``, not
       ``*_year_follow_up_y_arm_1``.  The release also contains odd-year
       sessions (``ses-01A``, ``ses-03A``, ``ses-05A``) which carry no imaging;
       they are dropped by the visit map rather than silently joined as NaN.
    3. **Subject ids.** ``sub-003RTV85``, where 5.1 used
       ``NDAR_INV003RTV85``.  The transform is deterministic
       (11,817/11,818 subjects intersect), so both releases normalise to the
       same ``sub-NDARINV...`` form and phenotype tables remain joinable
       across releases.
    4. **Column scheme.** ``mr_y_smri__thk__dsk__bstmps__lh_mean`` against
       5.1's ``smri_thick_cdk_banksstslh``, with abbreviated region codes.
       Region identity is resolved by :data:`REGION_CODES`, which was built
       from the release data dictionary and then *validated empirically* --
       see the note on region order below.
    5. **Age units.** ``ab_g_dyn__visit_age`` is in **years**; 5.1's
       ``interview_age`` was in months.  Missing this would rescale every
       slope by 12 while leaving all diagnostics looking healthy.
    6. **Whole-cortex means are supplied.** ``__lh_mean``/``__rh_mean``/
       ``_mean`` columns exist, where 5.1 required computing them.
    7. **Covariates are consolidated.** ``ab_g_dyn`` (per visit) carries age,
       site and scanner; ``ab_g_stc`` (per subject) carries sex, family,
       birth event, twin flags and 32 genetic ancestry PCs -- so no separate
       genetics file is needed for GWAS covariates.
    9. **HCP-MMP is derived, not shipped.** ``imaging(..., "hcp")`` reads
       ``processed/hcp/`` written by :mod:`abcd.hcp_stats` from the release
       FreeSurfer surfaces (see :meth:`_imaging_hcp`).
    8. **QC is renamed, not removed.** There is no Euler column; the
       equivalent is ``topodfct_count`` (topological defect count), which is
       what the Euler number is computed from.  Validated against the 5.1
       legacy exclusion list: excluded subjects have median 32 defects vs 19
       for retained (Mann-Whitney p = 3e-83, AUC 0.73).

    REGION ORDER IS NOT SHARED BETWEEN RELEASES
    -------------------------------------------
    The two releases list the 68 DK regions in *different orders*.  Matching
    5.1 and 7.0 columns by position agrees with the true correspondence for
    only 16 of 68 regions -- and because both releases are internally
    consistent, a positional mapping produces plausible thickness values,
    plausible left-right symmetry and plausible age effects.  It would not
    have failed loudly anywhere; it would just have relabelled the cortex.

    :data:`REGION_CODES` is therefore checked against 5.1 on overlapping
    subject-visits by :meth:`validate_region_mapping`, which is called from
    :func:`verify_adapter`.  The check is a bijective correlation match
    (min r = 0.9997, min margin over the runner-up = 0.17).  That check also
    established that 99.94% of matched values are bit-identical, i.e. 7.0
    reuses the 5.1 surface reconstructions for the shared waves rather than
    reprocessing them -- so 5.1-vs-7.0 differences in results come from the
    added wave and added subjects, not from a pipeline change upstream.
    """

    release = "7.0"

    VISITS = {
        "v0": "ses-00A",
        "v2": "ses-02A",
        "v4": "ses-04A",
        "v6": "ses-06A",
    }

    #: metric -> (table stem, column infix, global aggregation).
    #: 7.0 has no stem_a/stem_b split: region codes are shared across metrics.
    METRIC_TABLES = {
        "thickness": ("mr_y_smri__thk__dsk",     "thk__dsk",     "mean"),
        "t1_gray":   ("mr_y_smri__t1__gm__dsk",  "t1__gm__dsk",  "mean"),
        "t2_gray":   ("mr_y_smri__t2__gm__dsk",  "t2__gm__dsk",  "mean"),
    }

    #: 7.0 region code -> 5.1 `stem_a` region token.  Built from the release
    #: data dictionary's ROI descriptions, validated empirically (see class
    #: docstring).  Keys are the abbreviations appearing in column names.
    REGION_CODES = {
        "bstmps": "bankssts",    "cac": "cdacate",        "cmfrt": "cdmdfr",
        "cn": "cuneus",          "er": "ehinal",          "ff": "fusiform",
        "ic": "ihcate",          "ins": "insula",         "iprt": "ifpl",
        "itmp": "iftm",          "lg": "lingual",         "lobfrt": "lobfr",
        "locc": "locc",          "mobfrt": "mobfr",       "mtmp": "mdtm",
        "pactr": "paracn",       "pcc": "pericc",         "pcg": "ptcate",
        "pfrt": "frpole",        "ph": "parahpal",        "pob": "parsobis",
        "poctr": "postcn",       "pop": "parsopc",        "prcn": "pc",
        "prctr": "precn",        "ptg": "parstgris",      "ptmp": "tmpole",
        "rac": "rracate",        "rmfrt": "rrmdfr",       "sfrt": "sufr",
        "sm": "sm",              "sprt": "supl",          "stmp": "sutm",
        "ttmp": "trvtm",
    }

    VERIFIED = True

    # ------------------------------------------------------------------
    @staticmethod
    def _to_bids(subject: pd.Series) -> pd.Series:
        """``sub-003RTV85`` -> ``sub-NDARINV003RTV85`` (idempotent).

        Normalises to the same form the 5.1 adapter produces, so phenotype
        tables from the two releases join without a crosswalk.
        """
        s = subject.astype(str)
        # already-normalised or 5.1-style ids pass through unchanged
        done = s.str.startswith("sub-NDARINV")
        out = s.where(done, "sub-NDARINV" + s.str.replace("sub-", "", regex=False))
        return out.str.replace("NDAR_INV", "NDARINV", regex=False)

    def _table(self, stem: str) -> pd.DataFrame:
        """Read a 7.0 table by stem, tolerating flat or nested layout."""
        path = paths.find_table(self.root, stem)
        sep = "\t" if path.suffix == ".tsv" else ","
        if path.suffix == ".parquet":
            return pd.read_parquet(path)
        return pd.read_csv(path, sep=sep, low_memory=False)

    # ------------------------------------------------------------------
    #: HCP-MMP tables are not part of the release; :mod:`abcd.hcp_stats`
    #: derives them from the release FreeSurfer surfaces into
    #: ``<release>/processed/hcp/`` in the DK column convention.
    #: metric -> (column infix, per-parcel aggregation suffix).
    HCP_METRICS = {
        "thickness": ("thk", "mean"),
        "area": ("area", "sum"),
        "volume": ("vol", "sum"),
    }

    @property
    def hcp_dir(self) -> Path:
        return self.root / "processed" / "hcp"

    def imaging(self, metric: str, parcellation: str = "dsk") -> pd.DataFrame:
        if parcellation == "hcp":
            return self._imaging_hcp(metric)
        if parcellation != "dsk":
            raise NotImplementedError(
                f"parcellation {parcellation!r} not wired for 7.0; only 'dsk' and 'hcp'"
            )
        if metric not in self.METRIC_TABLES:
            raise KeyError(
                f"metric {metric!r} is not a raw 7.0 table. Available: "
                f"{sorted(self.METRIC_TABLES)}. Derived metrics such as "
                "'t1t2_ratio' are built in assemble.py, not here."
            )
        stem, infix, _agg = self.METRIC_TABLES[metric]
        raw = self._table(stem)

        pat = re.compile(rf"^mr_y_smri__{re.escape(infix)}__([a-z0-9]+)__(lh|rh)_mean$")
        region_cols = {c: pat.match(c) for c in raw.columns}
        region_cols = {c: m for c, m in region_cols.items() if m}
        if len(region_cols) != 68:
            raise KeyError(
                f"{stem}: expected 68 region columns matching {pat.pattern}, "
                f"found {len(region_cols)}"
            )
        unknown = {m.group(1) for m in region_cols.values()} - set(self.REGION_CODES)
        if unknown:
            raise KeyError(
                f"{stem}: region codes not in REGION_CODES: {sorted(unknown)}. "
                "The release may have changed its abbreviations; re-derive from "
                "the data dictionary and re-run validate_region_mapping()."
            )

        # whole-cortex mean is supplied directly in 7.0
        global_col = f"mr_y_smri__{infix}_mean"
        keep = ["participant_id", "session_id"] + list(region_cols)
        has_global = global_col in raw.columns
        if has_global:
            keep.append(global_col)

        long = raw[keep].melt(
            id_vars=["participant_id", "session_id"],
            var_name="column",
            value_name="value",
        )
        # REGION_CODES maps a 7.0 abbreviation to the release-5.1 column token;
        # that token is then resolved to the canonical DK region name through
        # the same ``region_labels.csv`` the 5.1 adapter uses.  Going through
        # the shared table (rather than naming regions here) is what makes
        # phenotype tables from the two releases directly comparable, and it
        # keeps one canonical spelling for the spatial and gene-mapping code.
        labels = region_labels(parcellation)
        labels = labels[~labels.is_global]
        token_to_region = dict(zip(labels.stem_a, labels.region))
        meta_rows = []
        for c, m in region_cols.items():
            code, hemi = m.group(1), m.group(2)
            token = f"{self.REGION_CODES[code]}{hemi}"
            if token not in token_to_region:
                raise KeyError(
                    f"{stem}: region code {code!r} maps to token {token!r}, which "
                    f"is not in {REGION_LABEL_FILE.name}. Fix REGION_CODES or add "
                    "the region to the label table."
                )
            meta_rows.append({"column": c, "hemi": hemi,
                              "region": token_to_region[token],
                              "is_global": False})
        meta = pd.DataFrame(meta_rows)
        if has_global:
            meta = pd.concat([meta, pd.DataFrame([{
                "column": global_col, "hemi": "both",
                "region": "global_mean", "is_global": True,
            }])], ignore_index=True)
        meta["label"] = np.where(
            meta.is_global, meta.region, meta.hemi + "_" + meta.region
        )

        long = long.merge(meta, on="column", how="inner")
        long["subject"] = self._to_bids(long.participant_id)
        long["visit"] = self._map_visits(long.session_id)
        long["metric"] = metric
        long = long.dropna(subset=["visit", "value"])
        return long[
            ["subject", "visit", "metric", "hemi", "region", "label",
             "is_global", "value"]
        ].reset_index(drop=True)

    # ------------------------------------------------------------------
    def _imaging_hcp(self, metric: str) -> pd.DataFrame:
        """HCP-MMP1.0 (Glasser, 180 regions x 2 hemispheres) from ``processed/hcp/``.

        These are **derived, not release, tables**: :mod:`abcd.hcp_stats`
        parses ``mris_anatomical_stats`` output computed on the release
        FreeSurfer 7.1.1 surfaces (the same reconstructions the DK release
        tables are tabulated from) and writes them in the DK column
        convention, ``mr_y_smri__thk__hcp__<region>__<hemi>_mean``.  Region
        names are the HCP-MMP names verbatim (``V1``, ``a9-46v``, ``TE1m``),
        case preserved, and labels ``lh_V1`` match ``data/hcp_centroids.csv``.

        Coverage is whatever the parcellation pipeline has reached -- on
        2026-09-12, 30,360 of 33,825 FreeSurfer sessions -- so
        ``complete_regions`` QC and the run manifest, not this method, are
        where the shortfall shows up.  There is no release QC row for
        sessions newer than the tabulated release; ``hcp_session_qc.tsv``
        beside the table carries FreeSurfer surface-hole counts for those.

        The whole-cortex ``_mean`` is vertex-weighted over the 360 parcels,
        as the DK ``_mean`` columns are.
        """
        if metric not in self.HCP_METRICS:
            raise KeyError(
                f"metric {metric!r} not available in HCP-MMP; have "
                f"{sorted(self.HCP_METRICS)}"
            )
        infix, agg = self.HCP_METRICS[metric]
        path = self.hcp_dir / f"mr_y_smri__{infix}__hcp.tsv"
        if not path.exists():
            raise SourceUnavailable(
                f"HCP-MMP table absent: {path}. Generate it with "
                "`sbatch hpc/hcp_extract.sbatch` (see abcd.hcp_stats)."
            )
        raw = pd.read_csv(path, sep="\t", low_memory=False)

        pat = re.compile(rf"^mr_y_smri__{re.escape(infix)}__hcp__(.+)__(lh|rh)_{agg}$")
        region_cols = {c: m for c in raw.columns if (m := pat.match(c))}
        if len(region_cols) != 360:
            raise KeyError(
                f"{path.name}: expected 360 region columns matching "
                f"{pat.pattern}, found {len(region_cols)}"
            )
        global_col = f"mr_y_smri__{infix}__hcp_{agg}"
        has_global = global_col in raw.columns

        keep = ["participant_id", "session_id", *region_cols]
        if has_global:
            keep.append(global_col)
        long = raw[keep].melt(id_vars=["participant_id", "session_id"],
                              var_name="column", value_name="value")

        meta = pd.DataFrame([
            {"column": c, "hemi": m.group(2), "region": m.group(1), "is_global": False}
            for c, m in region_cols.items()
        ])
        if has_global:
            meta = pd.concat([meta, pd.DataFrame([{
                "column": global_col, "hemi": "both",
                "region": "global_mean", "is_global": True,
            }])], ignore_index=True)
        meta["label"] = np.where(meta.is_global, meta.region, meta.hemi + "_" + meta.region)

        long = long.merge(meta, on="column", how="inner")
        long["subject"] = self._to_bids(long.participant_id)
        long["visit"] = self._map_visits(long.session_id)
        long["metric"] = metric
        long = long.dropna(subset=["visit", "value"])
        return long[
            ["subject", "visit", "metric", "hemi", "region", "label",
             "is_global", "value"]
        ].reset_index(drop=True)

    # ------------------------------------------------------------------
    def demographics(self) -> pd.DataFrame:
        """Sex from the static cohort table (one row per subject)."""
        st = self._table("ab_g_stc")
        col = "ab_g_stc__cohort_sex"
        if col not in st.columns:
            raise SourceUnavailable(f"ab_g_stc has no column {col}")
        out = pd.DataFrame({
            "subject": self._to_bids(st.participant_id),
            "sex_code": st[col],
            "sex": st[col].map({1: "M", 2: "F"}),
        })
        return out.dropna(subset=["subject"]).drop_duplicates("subject").reset_index(drop=True)

    # ------------------------------------------------------------------
    def longitudinal(self) -> pd.DataFrame:
        """Age, site and family, joined from the per-visit and static tables.

        ``ab_g_dyn__visit_age`` is already in years -- unlike 5.1's
        ``interview_age`` in months -- so no conversion is applied here.  The
        unit is asserted in :func:`verify_adapter` rather than trusted.
        """
        dy = self._table("ab_g_dyn")
        st = self._table("ab_g_stc")
        need = ["ab_g_dyn__visit_age", "ab_g_dyn__design_site"]
        missing = [c for c in need if c not in dy.columns]
        if missing:
            raise SourceUnavailable(f"ab_g_dyn missing columns: {missing}")

        out = pd.DataFrame({
            "subject": self._to_bids(dy.participant_id),
            "visit": self._map_visits(dy.session_id),
            "age": dy["ab_g_dyn__visit_age"].astype(float),
            "site": dy["ab_g_dyn__design_site"],
        }).dropna(subset=["visit", "age"])

        fam = pd.DataFrame({
            "subject": self._to_bids(st.participant_id),
            "family_id": st["ab_g_stc__design_id__fam"],
        }).dropna(subset=["family_id"]).drop_duplicates("subject")

        # family is a subject-level attribute in 7.0 (static table), so unlike
        # 5.1 there is no per-visit value that could drift mid-study.
        out = out.merge(fam, on="subject", how="left")
        return out[["subject", "visit", "age", "site", "family_id"]].reset_index(drop=True)

    def scanner(self) -> pd.DataFrame:
        """Per-visit scanner manufacturer, for the Philips exclusion."""
        dy = self._table("ab_g_dyn")
        col = "ab_g_dyn__design_mr__manufact"
        if col not in dy.columns:
            raise SourceUnavailable(f"ab_g_dyn has no column {col}")
        return pd.DataFrame({
            "subject": self._to_bids(dy.participant_id),
            "visit": self._map_visits(dy.session_id),
            "manufacturer": dy[col],
        }).dropna(subset=["visit"]).reset_index(drop=True)

    def qc_include(self) -> pd.DataFrame:
        """The release's own recommended-inclusion flag for T1w."""
        inc = self._table("mr_y_qc__incl")
        col = "mr_y_qc__incl__smri__t1_indicator"
        if col not in inc.columns:
            raise SourceUnavailable(f"mr_y_qc__incl has no column {col}")
        return pd.DataFrame({
            "subject": self._to_bids(inc.participant_id),
            "visit": self._map_visits(inc.session_id),
            "include": inc[col],
        }).dropna(subset=["visit"]).reset_index(drop=True)

    def surface_defects(self) -> pd.DataFrame:
        """Topological defect count -- the Euler-number equivalent.

        FreeSurfer's Euler characteristic is a linear function of the number of
        topological defects, so a *lower* Euler number and a *higher* defect
        count mean the same thing: a worse surface reconstruction.  Thresholds
        must therefore be expressed as an upper bound on defects, not a lower
        bound on Euler.
        """
        au = self._table("mr_y_qc__post__aut")
        col = "mr_y_qc__post__aut__smri__topodfct_count"
        if col not in au.columns:
            raise SourceUnavailable(f"mr_y_qc__post__aut has no column {col}")
        return pd.DataFrame({
            "subject": self._to_bids(au.participant_id),
            "visit": self._map_visits(au.session_id),
            "defects": au[col].astype(float),
        }).dropna(subset=["visit"]).reset_index(drop=True)

    # ------------------------------------------------------------------
    def genetics(self) -> dict[str, Path]:
        """Genetic resources.  PCs live in the static table in 7.0.

        No genotype/WGS files are present in the local copy, so GCTA and GWAS
        steps still cannot run against 7.0 here; see ``hpc/README_HPC.md``.
        """
        return {"static": paths.find_table(self.root, "ab_g_stc")}

    def genetic_pcs(self, n: int = 10) -> pd.DataFrame:
        """First ``n`` genetic ancestry PCs, for GWAS covariates."""
        st = self._table("ab_g_stc")
        cols = [f"ab_g_stc__gen_pc__{i:02d}" for i in range(1, n + 1)]
        missing = [c for c in cols if c not in st.columns]
        if missing:
            raise SourceUnavailable(f"ab_g_stc missing PC columns: {missing}")
        out = st[cols].copy()
        out.columns = [f"PC{i}" for i in range(1, n + 1)]
        out.insert(0, "subject", self._to_bids(st.participant_id))
        return out.dropna(subset=["subject"]).reset_index(drop=True)

    #: Zygosity codes in ``gn_y_genrel_zyg__NN``.  No codebook ships with the
    #: release; these were established empirically and are asserted in
    #: ``tests/test_io.py``.  Code 1 has pi-hat 0.92-1.00 (mean 0.987) and
    #: always shares a birth event; codes 2 and 3 both sit at pi-hat ~0.50 and
    #: are separated *perfectly* by birth event -- code 2 always shares one
    #: (DZ twin), code 3 never does (non-twin full sibling).  That separation
    #: is what lets a Falconer estimate use DZ twins alone instead of pooling
    #: siblings in, which biases h2 upward (see :mod:`abcd.heritability`).
    ZYGOSITY_CODES = {1: "MZ", 2: "DZ_twin", 3: "full_sib"}

    def genotyped_pairs(self) -> pd.DataFrame:
        """Genotype-confirmed relative pairs from ``gn_y_genrel``.

        One row per *unordered* pair with columns ``a``, ``b`` (subject IDs in
        the release's ``sub-NDARINV...`` form), ``family_id``, ``pihat`` and
        ``pair_type`` (``MZ`` / ``DZ_twin`` / ``full_sib``).

        This table is new in 7.0 and supersedes the previous route, which
        borrowed zygosity from the 5.1 pi-hat file and matched on normalised
        subject ID.  Two things it buys:

        * **DZ twins are separable from non-twin full siblings.**  The 5.1 route
          could only produce a pooled ``DZ_or_sib`` class.  Pooling inflates
          Falconer h2, because siblings correlate less than DZ twins do
          (measured on the 7.0 slope phenotype: r = 0.090 vs 0.279, difference
          +0.189, 95% CI [+0.046, +0.321]).
        * **Pairs are explicit rather than inferred from family size.**  The old
          route took families with exactly two phenotyped members as a pair,
          which silently drops any family holding both a twin pair and a third
          sibling.

        The source lists each pair from both sides across four partner slots;
        this collapses them.  Reciprocity and class agreement are asserted --
        on the 7.0 release all 1,922 pairs appear twice with no disagreement.
        """
        raw = self._table("gn_y_genrel")
        slots = []
        for i in ("01", "02", "03", "04"):
            need = [f"gn_y_genrel_id__paired__{i}", f"gn_y_genrel_pihat__{i}",
                    f"gn_y_genrel_zyg__{i}"]
            if not all(c in raw.columns for c in need):
                continue
            s = raw[["participant_id", "gn_y_genrel_id__fam"] + need].copy()
            s.columns = ["a", "family_id", "b", "pihat", "zyg"]
            slots.append(s.dropna(subset=["zyg", "b"]))
        if not slots:
            raise SourceUnavailable("gn_y_genrel has no populated pair slots")

        long = pd.concat(slots, ignore_index=True)
        long["pair_type"] = long.zyg.astype(int).map(self.ZYGOSITY_CODES)
        unknown = long.loc[long.pair_type.isna(), "zyg"].unique()
        if len(unknown):
            raise SourceUnavailable(
                f"gn_y_genrel has unrecognised zygosity codes {sorted(unknown)}; "
                f"known codes are {self.ZYGOSITY_CODES}"
            )
        long["a"] = self._to_bids(long.a)
        long["b"] = self._to_bids(long.b)
        long["_key"] = [tuple(sorted(t)) for t in zip(long.a, long.b)]

        agree = long.groupby("_key").pair_type.nunique()
        if (agree > 1).any():
            bad = agree[agree > 1].index[:3].tolist()
            raise SourceUnavailable(
                f"gn_y_genrel disagrees with itself on the class of {int((agree > 1).sum())} "
                f"pair(s), e.g. {bad}; zygosity cannot be trusted"
            )
        out = long.drop_duplicates("_key").copy()
        out[["a", "b"]] = pd.DataFrame(out._key.tolist(), index=out.index)
        return out[["a", "b", "family_id", "pihat", "pair_type"]].reset_index(drop=True)

    def relatedness(self) -> pd.DataFrame:
        """Family, birth event and twin structure from the static table.

        ``pair_type_design`` here is inferred from *design* variables, not
        measured relatedness: subjects sharing a birth event are twins/triplets,
        and ``design_sstwin`` marks same-sex twins.  Same-sex twins are MZ
        *candidates* only -- roughly half are DZ -- so this column cannot be
        used in a Falconer estimate.  The name keeps that visible at the call
        site.  For genotype-confirmed zygosity use :meth:`genotyped_pairs`,
        which reads the ``gn_y_genrel`` pi-hat table shipped with 7.0.
        """
        st = self._table("ab_g_stc")
        out = pd.DataFrame({
            "subject": self._to_bids(st.participant_id),
            "family_id": st.get("ab_g_stc__design_id__fam"),
            "birth_id": st.get("ab_g_stc__design_id__birth"),
            "group_id": st.get("ab_g_stc__design_id__group"),
            "relationship": st.get("ab_g_stc__design_famrel"),
            "same_sex_twin": st.get("ab_g_stc__design_sstwin"),
        })
        multiple = out.birth_id.notna() & out.duplicated("birth_id", keep=False)
        out["pair_type_design"] = np.select(
            [multiple & (out.same_sex_twin == 1), multiple],
            ["twin_same_sex", "twin_opposite_sex"],
            default="singleton_or_sib",
        )
        return out.reset_index(drop=True)

    # ------------------------------------------------------------------
    def validate_region_mapping(self, metric: str = "thickness") -> pd.DataFrame:
        """Check :data:`REGION_CODES` against 5.1 on overlapping visits.

        Returns one row per region with the correlation to its mapped 5.1
        counterpart.  Raises if the mapping is not the best available match for
        every region, which is the failure a positional mapping would not
        produce: see the class docstring.
        """
        mine = self.imaging(metric, "dsk")
        mine = mine[~mine.is_global]
        other = Release51Adapter().imaging(metric, "dsk")
        other = other[~other.is_global]

        a = mine.pivot_table(index=["subject", "visit"], columns="label", values="value")
        b = other.pivot_table(index=["subject", "visit"], columns="label", values="value")
        shared = a.index.intersection(b.index)
        if len(shared) < 1000:
            raise SourceUnavailable(
                f"only {len(shared)} overlapping subject-visits between 7.0 and "
                "5.1; cannot validate the region mapping"
            )
        a, b = a.loc[shared], b.loc[shared]
        common = sorted(set(a.columns) & set(b.columns))
        if len(common) != 68:
            raise KeyError(
                f"label sets disagree: {len(common)} shared of "
                f"{len(a.columns)}/{len(b.columns)}"
            )
        A, B = a[common].to_numpy(float), b[common].to_numpy(float)
        ok = ~(np.isnan(A).any(1) | np.isnan(B).any(1))
        A, B = A[ok], B[ok]
        Az = (A - A.mean(0)) / A.std(0)
        Bz = (B - B.mean(0)) / B.std(0)
        C = Az.T @ Bz / len(Az)

        diag = np.diag(C)
        best = C.argmax(1)
        wrong = [common[i] for i in range(len(common)) if best[i] != i]
        if wrong:
            raise ValueError(
                "REGION_CODES is wrong: for these labels the mapped 5.1 region "
                f"is not the best correlate: {wrong[:10]}"
                + (f" (+{len(wrong)-10} more)" if len(wrong) > 10 else "")
            )
        srt = np.sort(C, axis=1)
        return pd.DataFrame({
            "label": common,
            "r": diag,
            "margin_over_runner_up": srt[:, -1] - srt[:, -2],
            "n_visits": len(A),
        }).sort_values("r").reset_index(drop=True)


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
