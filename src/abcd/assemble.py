"""
Assemble a modelling table from raw release data.

This is the one place that turns "a config" into "a tidy table R can fit".
Everything downstream -- the R model scripts, the phenotype step, the spatial
analyses -- consumes the Parquet this module writes and never touches a
release CSV.

Metric-agnostic by construction
-------------------------------
A metric is either *raw* (a table the adapter can read directly) or *derived*
(a function of two or more raw metrics, registered in :data:`DERIVED_METRICS`).
Adding the T1w/T2w ratio, a myelin proxy, or a future resting-state measure
means registering a function, not editing the pipeline.  The long schema

    subject | visit | hemi | region | value

carries any region-level metric unchanged, which is what makes
multi-metric joint modelling a downstream concern rather than a rewrite.

The global covariate
--------------------
The 5.1 pipeline used a *predicted* hemisphere mean from a separate model as a
regressor, which propagates that model's uncertainty into every regional fit
while treating it as known.  Here the global covariate is the **observed**
hemisphere mean for the same scan (``global_covariate='observed_mean'``), or
omitted entirely (``'none'``).  It is centred within-subject and
between-subject separately so the regional models can distinguish a subject
whose whole cortex is thin from a subject whose cortex thinned at this visit.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from . import io, paths, qc
from .config import RunConfig

# --------------------------------------------------------------------------
# Derived metrics
# --------------------------------------------------------------------------

DerivedFn = Callable[[dict[str, pd.DataFrame]], pd.DataFrame]

#: derived metric -> (raw metrics it needs, combining function)
DERIVED_METRICS: dict[str, tuple[tuple[str, ...], DerivedFn]] = {}


def register_derived(name: str, requires: tuple[str, ...]):
    """Decorator registering a derived metric built from raw metric frames."""

    def wrap(fn: DerivedFn) -> DerivedFn:
        DERIVED_METRICS[name] = (requires, fn)
        return fn

    return wrap


def _align_on_key(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Inner-join raw long frames on the region key, one value column each."""
    key = ["subject", "visit", "hemi", "region", "label", "is_global"]
    out = None
    for name, df in frames.items():
        part = df[key + ["value"]].rename(columns={"value": name})
        out = part if out is None else out.merge(part, on=key, how="inner")
    return out


@register_derived("t1t2_ratio", ("t1_gray", "t2_gray"))
def _t1t2_ratio(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """T1w/T2w grey-matter intensity ratio, a common myelin proxy.

    Both terms come from the same scan, so scanner intensity scaling largely
    cancels -- which is the whole point of the ratio.
    """
    j = _align_on_key(frames)
    j["value"] = j["t1_gray"] / j["t2_gray"].replace(0, np.nan)
    j["metric"] = "t1t2_ratio"
    return j.drop(columns=["t1_gray", "t2_gray"])


@register_derived("t1t2_contrast_ratio", ("t1_contrast", "t2_contrast"))
def _t1t2_contrast_ratio(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Ratio of grey/white contrast in T1w vs T2w."""
    j = _align_on_key(frames)
    j["value"] = j["t1_contrast"] / j["t2_contrast"].replace(0, np.nan)
    j["metric"] = "t1t2_contrast_ratio"
    return j.drop(columns=["t1_contrast", "t2_contrast"])


def load_metric(adapter, metric: str, parcellation: str) -> pd.DataFrame:
    """Load a raw or derived metric in the canonical long schema."""
    if metric in DERIVED_METRICS:
        requires, fn = DERIVED_METRICS[metric]
        frames = {m: adapter.imaging(m, parcellation) for m in requires}
        return fn(frames)
    return adapter.imaging(metric, parcellation)


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------

def assemble(cfg: RunConfig, adapter=None, verbose: bool = True
             ) -> tuple[pd.DataFrame, dict]:
    """Build the modelling table for a config.

    Returns ``(table, manifest)``.  The table is long over regions with the
    covariates joined on; the manifest records every sample decision and the
    QC ledger, so a fitted result can always be traced back to the rows it saw.
    """
    adapter = adapter or io.get_adapter(cfg.release)
    manifest: dict = {
        "config": cfg.to_dict(),
        "config_hash": cfg.hash,
        "run_id": cfg.run_id,
        "assembled_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "release": cfg.release,
        "stages": [],
    }
    # Which tables this run actually saw.  The 6.0 and 7.0 tabulations are
    # column-identical; only the six-year row count distinguishes them, and it
    # was not recorded anywhere for two months.  ``assert_vintage`` refuses to
    # assemble a 7.0 run from 6.0-sized tables.
    if hasattr(adapter, "assert_vintage"):
        base_metric = cfg.metric if cfg.metric in getattr(adapter, "METRIC_TABLES", {}) else "thickness"
        manifest["data_vintage"] = adapter.assert_vintage(base_metric)
        if verbose:
            six = manifest["data_vintage"]["imaging_rows_by_session"]
            print(f"  data vintage: {manifest['data_vintage']['release_dir']}  rows/session={six}")

    def stage(name, df, **extra):
        rec = {"stage": name, "n_rows": len(df),
               "n_subjects": int(df.subject.nunique()) if "subject" in df else None,
               **extra}
        manifest["stages"].append(rec)
        if verbose:
            print(f"  {name:34s} rows={len(df):>9,}  subjects={rec['n_subjects']}")
        return df

    if verbose:
        print(f"assemble: {cfg.metric} / {cfg.release} / {cfg.run_id}")

    # 1. imaging -----------------------------------------------------------
    img = load_metric(adapter, cfg.metric, cfg.parcellation)
    img = img[img.visit.isin(cfg.visits)]
    if cfg.exclude_regions:
        drop = img.region.isin(cfg.exclude_regions) & ~img.is_global
        n_lab = img.loc[drop, "label"].nunique()
        img = img[~drop]
        if verbose:
            print(f"  excluded regions {list(cfg.exclude_regions)}: {n_lab} labels dropped, "
                  f"{img.label.nunique()} remain")
    stage("imaging loaded", img)

    n_regions = int(img.loc[~img.is_global, ["hemi", "region"]].drop_duplicates().shape[0])

    # 2. QC on scans -------------------------------------------------------
    scans = qc.scan_table(img)
    preds = qc.build_policy(cfg.qc_policy, adapter, n_expected_regions=n_regions)
    qc_res = qc.apply_qc(scans, preds, policy=cfg.qc_policy)
    manifest["qc_ledger"] = qc_res.ledger.to_dict("records")
    manifest["qc_unavailable"] = qc_res.unavailable
    if verbose and qc_res.unavailable:
        print(f"  ! QC criteria unavailable for this release: {qc_res.unavailable}")

    img = img.merge(qc_res.scans[["subject", "visit"]], on=["subject", "visit"])
    stage("after QC", img, policy=cfg.qc_policy)

    # 3. covariates --------------------------------------------------------
    lt = adapter.longitudinal()
    demo = adapter.demographics()
    img = img.merge(lt[["subject", "visit", "age", "site", "family_id"]],
                    on=["subject", "visit"], how="inner")
    img = img.merge(demo[["subject", "sex"]].dropna(subset=["sex"]),
                    on="subject", how="inner")
    stage("covariates joined", img)

    # 4. global covariate --------------------------------------------------
    glob = img[img.is_global].copy()
    reg = img[~img.is_global].copy()

    if cfg.global_covariate == "observed_mean":
        # per (subject, visit, hemi) global value; fall back to the regional
        # mean when the release's global column is missing for that scan
        g = (glob.groupby(["subject", "visit", "hemi"], observed=True)
             .value.mean().rename("global_value").reset_index())
        g_both = g[g.hemi == "both"][["subject", "visit", "global_value"]]
        if cfg.hemisphere == "both" and len(g_both):
            reg = reg.merge(g_both, on=["subject", "visit"], how="left")
        else:
            reg = reg.merge(g[g.hemi.isin(["lh", "rh"])],
                            on=["subject", "visit", "hemi"], how="left")
        fallback = reg.global_value.isna()
        if fallback.any():
            rm = (reg.groupby(["subject", "visit"], observed=True)
                  .value.mean().rename("_rm").reset_index())
            reg = reg.merge(rm, on=["subject", "visit"], how="left")
            reg.loc[fallback, "global_value"] = reg.loc[fallback, "_rm"]
            reg = reg.drop(columns="_rm")
        manifest["global_covariate_fallback_rows"] = int(fallback.sum())
    elif cfg.global_covariate == "none":
        reg["global_value"] = np.nan
    else:
        raise ValueError(f"unsupported global_covariate {cfg.global_covariate!r}")

    # split the global covariate into between- and within-subject parts, so a
    # regional model can separate "this cortex is thin" from "this cortex
    # thinned at this visit"
    if cfg.global_covariate != "none":
        subj_mean = (reg.drop_duplicates(["subject", "visit"])
                     .groupby("subject").global_value.mean().rename("global_between"))
        reg = reg.merge(subj_mean, on="subject", how="left")
        reg["global_within"] = reg.global_value - reg.global_between
        reg["global_between_c"] = reg.global_between - reg.global_between.mean()
    else:
        reg["global_between"] = np.nan
        reg["global_within"] = np.nan
        reg["global_between_c"] = np.nan

    # 5. hemisphere selection ---------------------------------------------
    if cfg.hemisphere in ("lh", "rh"):
        reg = reg[reg.hemi == cfg.hemisphere]
    stage("hemisphere selected", reg, hemisphere=cfg.hemisphere)

    # 6. sample rule -------------------------------------------------------
    scans_ok = reg[["subject", "visit"]].drop_duplicates()
    subs, sample_info = qc.select_subjects(
        scans_ok, cfg.sample_rule, cfg.min_visits, cfg.visits
    )
    manifest["sample"] = sample_info
    reg = reg[reg.subject.isin(subs)]
    stage("sample rule applied", reg, **sample_info)

    # 7. age coding --------------------------------------------------------
    centre = cfg.age_centre if cfg.age_centre is not None else float(reg.age.mean())
    manifest["age_centre"] = round(centre, 4)
    reg["age_c"] = reg.age - centre

    # per-subject age summary: needed to flag which slopes are interpolated
    span = (reg.drop_duplicates(["subject", "visit"])
            .groupby("subject").age.agg(["min", "max", "count"])
            .rename(columns={"min": "age_first", "max": "age_last",
                             "count": "n_visits"}))
    span["age_span"] = span.age_last - span.age_first
    reg = reg.merge(span, on="subject", how="left")

    # 8. tidy up -----------------------------------------------------------
    reg["sex"] = pd.Categorical(reg.sex, categories=["F", "M"])
    reg["site"] = reg.site.astype("category")
    reg["family_id"] = reg.family_id.astype("Int64").astype("string")
    reg["metric"] = cfg.metric
    reg["release"] = cfg.release

    cols = [
        "subject", "visit", "metric", "release", "hemi", "region", "label",
        "value", "age", "age_c", "sex", "site", "family_id",
        "global_value", "global_between", "global_between_c", "global_within",
        "age_first", "age_last", "age_span", "n_visits",
    ]
    out = reg[cols].dropna(subset=["value"]).reset_index(drop=True)

    manifest["final"] = {
        "n_rows": len(out),
        "n_subjects": int(out.subject.nunique()),
        "n_scans": int(out[["subject", "visit"]].drop_duplicates().shape[0]),
        "n_regions": int(out[["hemi", "region"]].drop_duplicates().shape[0]),
        "n_sites": int(out.site.nunique()),
        "n_families": int(out.family_id.nunique()),
        "age_range": [round(float(out.age.min()), 3), round(float(out.age.max()), 3)],
        "visits_present": sorted(out.visit.unique().tolist()),
        "visit_counts": out.drop_duplicates(["subject", "visit"])
                           .visit.value_counts().sort_index().to_dict(),
        "n_visits_distribution": out.drop_duplicates("subject")
                                    .n_visits.value_counts().sort_index().to_dict(),
    }
    if verbose:
        f = manifest["final"]
        print(f"  final: {f['n_rows']:,} rows, {f['n_subjects']:,} subjects, "
              f"{f['n_regions']} regions, {f['n_sites']} sites, {f['n_families']:,} families")
    return out, manifest


def write(table: pd.DataFrame, manifest: dict, run_dir: Path | None = None) -> Path:
    """Write the table and its manifest to the run directory."""
    run_dir = run_dir or paths.run_dir(manifest["run_id"])
    run_dir.mkdir(parents=True, exist_ok=True)
    table.to_parquet(run_dir / "model_table.parquet", index=False)
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    (run_dir / "config.yaml").write_text(
        __import__("yaml").safe_dump(manifest["config"], sort_keys=False)
    )
    pd.DataFrame(manifest["qc_ledger"]).to_csv(run_dir / "qc_ledger.csv", index=False)
    return run_dir


def build(cfg: RunConfig, verbose: bool = True) -> tuple[pd.DataFrame, dict, Path]:
    """Assemble and write in one call. Returns (table, manifest, run_dir)."""
    table, manifest = assemble(cfg, verbose=verbose)
    run_dir = write(table, manifest)
    if verbose:
        print(f"  written -> {run_dir}")
    return table, manifest, run_dir


if __name__ == "__main__":  # pragma: no cover
    import argparse

    from .config import resolve_config

    ap = argparse.ArgumentParser(
        description="Assemble an ABCD modelling table.",
        epilog="config defaults to $ABCD_CONFIG.",
    )
    ap.add_argument("config", nargs="?", default=None,
                    help="config name or path; omit to use $ABCD_CONFIG")
    args = ap.parse_args()
    build(RunConfig.from_yaml(resolve_config(args.config)))
