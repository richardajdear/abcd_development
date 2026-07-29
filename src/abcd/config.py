"""
Run configuration for the ABCD longitudinal pipeline.

Every analysis decision that used to be encoded in a directory name
(``CT_y4_published_QC_lhmean_nosingular`` and its 35 siblings) is a field on
:class:`RunConfig`.  A run writes ``out/<hash>/config.yaml`` next to its
outputs, so the mapping from result to decisions is recorded rather than
inferred from a filename.

The hash is content-derived and stable across interpreter sessions: it is a
SHA-256 over the canonical JSON of the field dict, truncated to 12 hex chars.
Field *order* does not affect it; field *values* do.

    >>> cfg = RunConfig(metric="thickness")
    >>> cfg.hash == RunConfig.from_dict(cfg.to_dict()).hash
    True
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Literal

import yaml

# --------------------------------------------------------------------------
# Controlled vocabularies.  Keeping these as module constants (rather than
# free strings) means a typo fails at construction instead of silently
# producing a new config hash and a new output directory.
# --------------------------------------------------------------------------

RELEASES = ("5.1", "7.0")

#: Imaging metrics.  ``t1t2_ratio`` is derived (see ``assemble.DERIVED_METRICS``);
#: the rest map onto release tables via the adapter's ``METRIC_TABLES``.
METRICS = (
    "thickness",
    "area",
    "volume",
    "sulc",
    "t1_gray",
    "t2_gray",
    "t1_white",
    "t2_white",
    "t1_contrast",
    "t2_contrast",
    "t1t2_ratio",
)

PARCELLATIONS = ("dsk", "dst", "fzy", "hcp")  # Desikan-Killiany, Destrieux, fuzzy-cluster, HCP-MMP
HEMISPHERES = ("lh", "rh", "both")

#: How to decide which scans enter the model.
#:   ``legacy_euler``  reproduce the static SUBJECTS_MRI_ONLY_post_EULER list
#:   ``coded``         apply the predicate stack in ``qc.py``
#:   ``none``          no QC (diagnostic use only)
QC_POLICIES = ("legacy_euler", "coded", "none")

#: How to decide which *subjects* enter the model.
#:   ``min_visits``    subjects with >= ``min_visits`` QC-passing scans
#:   ``has_last``      subjects present at the final visit (the 5.1 thesis rule)
#:   ``all``           every subject with >=1 scan
SAMPLE_RULES = ("min_visits", "has_last", "all")

#: Random-effects structure to attempt first.  ``fallback_ladder`` in the R
#: fitter degrades from here when a fit is singular, and records what it used.
RE_STRUCTURES = ("slope_correlated", "slope_uncorrelated", "intercept_only")

#: Fixed-effect age basis.
AGE_BASES = ("linear", "spline")

#: Global-signal covariate.  ``observed_mean`` uses the per-visit hemisphere
#: mean measured at that visit.  ``none`` omits it.  Note there is deliberately
#: no ``predicted_mean`` option: the 5.1 pipeline used a *fitted, extrapolated*
#: hemisphere mean as a per-region covariate, which propagates first-stage
#: model error into the second stage without accounting for it.
GLOBAL_COVARIATES = ("none", "observed_mean")


class ConfigError(ValueError):
    """Raised when a RunConfig field is outside its controlled vocabulary."""


@dataclass(frozen=True)
class RunConfig:
    """A fully-specified pipeline run.

    Frozen so a config cannot drift after its hash is taken.  Use
    :meth:`replace` to derive a variant.
    """

    # --- data selection -------------------------------------------------
    release: str = "5.1"
    metric: str = "thickness"
    parcellation: str = "dsk"
    hemisphere: str = "both"

    # --- sample construction --------------------------------------------
    qc_policy: str = "coded"
    sample_rule: str = "min_visits"
    min_visits: int = 2
    #: Visits to include, in temporal order.  Release-specific event names are
    #: resolved by the adapter, so these stay abstract ("v0", "v2", "v4", "v6").
    visits: tuple[str, ...] = ("v0", "v2", "v4")

    # --- model ------------------------------------------------------------
    age_basis: str = "linear"
    spline_df: int = 3
    #: Age is centred at this value.  ``None`` means the QC-passing sample mean,
    #: which is resolved at fit time and written back into the run manifest.
    age_centre: float | None = None
    re_structure: str = "slope_correlated"
    global_covariate: str = "observed_mean"
    covariates: tuple[str, ...] = ("sex",)
    family_effect: bool = True
    site_effect: bool = True

    # --- bookkeeping ------------------------------------------------------
    #: Free-text note; excluded from the hash so annotating a run does not
    #: invalidate its outputs.
    note: str = ""

    #: Fields that do not participate in the identity hash.
    HASH_EXCLUDE: tuple[str, ...] = field(
        default=("note",), repr=False, compare=False
    )

    # ------------------------------------------------------------------
    def __post_init__(self) -> None:
        checks = {
            "release": RELEASES,
            "metric": METRICS,
            "parcellation": PARCELLATIONS,
            "hemisphere": HEMISPHERES,
            "qc_policy": QC_POLICIES,
            "sample_rule": SAMPLE_RULES,
            "age_basis": AGE_BASES,
            "re_structure": RE_STRUCTURES,
            "global_covariate": GLOBAL_COVARIATES,
        }
        for name, allowed in checks.items():
            value = getattr(self, name)
            if value not in allowed:
                raise ConfigError(
                    f"{name}={value!r} is not one of {allowed}"
                )
        if self.min_visits < 1:
            raise ConfigError(f"min_visits must be >= 1, got {self.min_visits}")
        if self.age_basis == "spline" and self.spline_df < 2:
            raise ConfigError(
                f"spline_df must be >= 2 for a spline basis, got {self.spline_df}"
            )
        if not self.visits:
            raise ConfigError("visits must not be empty")
        # tuples, not lists, so the dataclass stays hashable after a YAML round-trip
        object.__setattr__(self, "visits", tuple(self.visits))
        object.__setattr__(self, "covariates", tuple(self.covariates))

    # ------------------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        """Plain-Python dict, tuples flattened to lists (YAML/JSON friendly)."""
        d = asdict(self)
        d.pop("HASH_EXCLUDE", None)
        for k, v in d.items():
            if isinstance(v, tuple):
                d[k] = list(v)
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "RunConfig":
        known = {f.name for f in fields(cls)} - {"HASH_EXCLUDE"}
        unknown = set(d) - known
        if unknown:
            raise ConfigError(f"unknown config keys: {sorted(unknown)}")
        return cls(**d)

    # ------------------------------------------------------------------
    @property
    def hash(self) -> str:
        """Stable 12-char content hash of the identity-bearing fields."""
        d = self.to_dict()
        for k in self.HASH_EXCLUDE:
            d.pop(k, None)
        payload = json.dumps(d, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()[:12]

    @property
    def slug(self) -> str:
        """Human-readable prefix; the hash carries the actual identity."""
        return f"{self.metric}_{self.parcellation}_{self.release.replace('.', '')}"

    @property
    def run_id(self) -> str:
        return f"{self.slug}_{self.hash}"

    # ------------------------------------------------------------------
    def replace(self, **kwargs: Any) -> "RunConfig":
        """Derive a variant config (returns a new frozen instance)."""
        d = self.to_dict()
        d.update(kwargs)
        return RunConfig.from_dict(d)

    # ------------------------------------------------------------------
    def to_yaml(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            yaml.safe_dump(
                {"run_id": self.run_id, "hash": self.hash, **self.to_dict()},
                fh,
                sort_keys=False,
            )
        return path

    @classmethod
    def from_yaml(cls, path: str | Path) -> "RunConfig":
        with open(path) as fh:
            d = yaml.safe_load(fh)
        # run_id/hash are written for human readability; they are derived, not inputs
        recorded_hash = d.pop("hash", None)
        d.pop("run_id", None)
        cfg = cls.from_dict(d)
        if recorded_hash is not None and recorded_hash != cfg.hash:
            raise ConfigError(
                f"{path}: recorded hash {recorded_hash} != recomputed {cfg.hash}. "
                "The file was hand-edited, or the config schema changed."
            )
        return cfg


def find_run(out_dir: str | Path, *, metric: str | None = None,
             parcellation: str | None = None, release: str | None = None,
             require_fits: bool = True, **config_filters: Any) -> Path:
    """Resolve a run directory by *what it is*, not by its hash.

    ``run_id`` embeds a hash of the config, so any hardcoded id is specific to
    the machine whose config produced it.  Notebooks that hardcode one break for
    every other user (and for the same user after a config change).  Prefer::

        RUN = find_run(ROOT / "out", release="7.0", family_effect=False)

    Any keyword not named above is matched against the run's ``config.yaml``,
    so ``family_effect=False`` or ``min_visits=2`` work as filters.

    Raises ConfigError if zero or several runs match -- ambiguity here would
    silently analyse the wrong phenotype, so it is never resolved by guessing.
    """
    out_dir = Path(out_dir)
    if not out_dir.is_dir():
        raise ConfigError(f"{out_dir} does not exist; run `python -m abcd.assemble <config>` first")

    matches = []
    for d in sorted(out_dir.iterdir()):
        cfg_path = d / "config.yaml"
        if not d.is_dir() or not cfg_path.exists():
            continue
        with open(cfg_path) as fh:
            cfg = yaml.safe_load(fh) or {}
        want = {"metric": metric, "parcellation": parcellation, "release": release,
                **config_filters}
        if any(v is not None and cfg.get(k) != v for k, v in want.items()):
            continue
        if require_fits and not (d / "fits").is_dir():
            continue
        matches.append(d)

    criteria = {k: v for k, v in
                {"metric": metric, "parcellation": parcellation, "release": release,
                 **config_filters}.items() if v is not None}
    if not matches:
        avail = [p.name for p in out_dir.iterdir() if (p / "config.yaml").exists()]
        raise ConfigError(
            f"no run in {out_dir} matches {criteria}"
            + (f" (with fits/)" if require_fits else "")
            + f". Available: {avail or 'none'}"
        )
    if len(matches) > 1:
        raise ConfigError(
            f"{len(matches)} runs match {criteria}: {[m.name for m in matches]}. "
            "Add filters (e.g. family_effect=, min_visits=) to disambiguate -- "
            "picking one arbitrarily would risk analysing the wrong phenotype."
        )
    return matches[0]
