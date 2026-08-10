#!/usr/bin/env python
"""Regenerate the non-heritability tables in ``docs/`` from runs on disk.

Companion to ``tools/regen_h2_tables.py``, which owns every Falconer number.
Between the two, every ``docs/*.csv`` cited by ``docs/REPORT_7.0.md`` has a
generator, so a reviewer can re-run rather than trust.

Tables written here
-------------------
``fit_summary``            one row per fitted 7.0 run: design, n, convergence.
``reliability_grid``       BLUP slope reliability by min_visits x family_effect.
``site_scanner_icc``       ICC of the global slope by site and by scanner.
``developmental_maps_noglobal``  the 68-region map table used by sections 5-7.
``regional_slope_pc_loadings``   slope-PC loadings for the settled run.

A note on scanner
-----------------
Earlier versions of ``site_scanner_icc.csv`` carried a "baseline scanner serial"
row with k=29 levels.  Neither the local 7.0 tree nor the local 5.1 tree ships a
device-serial column (``mri_info_deviceserial`` is absent from both), so that row
cannot be reproduced and is not written.  7.0 exposes scanner *manufacturer*
only, via ``io.Release70Adapter.scanner``, which is 3 levels and removes
correspondingly less variance.  Report section 9 states this limitation.

Usage
-----
    PYTHONPATH=src python tools/regen_report_tables.py [--only TABLE ...]

Deterministic given ``SEED``.  Every table is rewritten in place.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import covariance as cov  # noqa: E402
from abcd import heritability as her  # noqa: E402
from abcd import io  # noqa: E402
from abcd import maps as M  # noqa: E402

#: The settled design (report section 2).  Must match regen_h2_tables.SETTLED.
SETTLED = "thickness_dsk_70_139406217085"
SEED = 0

#: Lobe assignment for the 34 DK bilateral regions, used by section 6.
LOBES = {
    "frontal": ("superiorfrontal rostralmiddlefrontal caudalmiddlefrontal parsopercularis "
                "parstriangularis parsorbitalis lateralorbitofrontal medialorbitofrontal "
                "precentral paracentral frontalpole"),
    "parietal": ("superiorparietal inferiorparietal supramarginal postcentral precuneus"),
    "temporal": ("superiortemporal middletemporal inferiortemporal bankssts fusiform "
                 "transversetemporal entorhinal temporalpole parahippocampal"),
    "occipital": ("lateraloccipital lingual cuneus pericalcarine"),
    "cingulate": ("rostralanteriorcingulate caudalanteriorcingulate posteriorcingulate "
                  "isthmuscingulate"),
    "insula": ("insula"),
}
_LOBE_OF = {r: lobe for lobe, rs in LOBES.items() for r in rs.split()}


def _runs(root: Path) -> pd.DataFrame:
    """Every 7.0 run with phenotypes, keyed by its design."""
    rows = []
    for d in sorted((root / "out").glob("thickness_dsk_70_*")):
        if not (d / "phenotypes" / "phenotypes.parquet").exists():
            continue
        c = yaml.safe_load((d / "config.yaml").read_text())
        rows.append({"run": d.name, "dir": str(d),
                     "global_cov": c.get("global_covariate"),
                     "min_visits": c.get("min_visits"),
                     "family_effect": bool(c.get("family_effect"))})
    if not rows:
        raise SystemExit("no fitted 7.0 runs found under out/ -- run the pipeline first")
    return pd.DataFrame(rows)


def _key(r) -> str:
    # global_covariate is the *string* "none" when absent, which is truthy --
    # compare explicitly or every run is labelled as global-adjusted.
    g = "ng" if r.global_cov in (None, "none", "") else "g"
    return f"{g}_mv{r.min_visits}{'_fam' if r.family_effect else ''}"


def _mean_visits(mt: pd.DataFrame) -> float:
    """Mean visits per subject.

    ``model_table`` is long over region x visit, so ``groupby(subject).size()``
    counts rows (68 x visits), not visits.  Count distinct visits.
    """
    return float(mt.groupby("subject")["visit"].nunique().mean())


def _global_slope(run_dir: str) -> pd.Series:
    ph = pd.read_parquet(Path(run_dir) / "phenotypes" / "phenotypes.parquet")
    return ph[ph.phenotype == "slope"].groupby("subject")["value"].mean()


def fit_summary(root: Path) -> pd.DataFrame:
    rows = []
    for r in _runs(root).itertuples():
        mt = pd.read_parquet(Path(r.dir) / "model_table.parquet")
        diag = Path(r.dir) / "fits" / "diagnostics.parquet"
        d = pd.read_parquet(diag) if diag.exists() else pd.DataFrame()
        rows.append({
            "key": _key(r), "run": r.run, "global_cov": r.global_cov,
            "min_visits": r.min_visits, "family_effect": r.family_effect,
            "subjects": mt.subject.nunique(),
            "mean_visits": round(_mean_visits(mt), 4),
            "singular": int(d.get("singular", pd.Series(dtype=bool)).sum()),
            "nonconverged": int(d.get("nonconverged", pd.Series(dtype=bool)).sum()),
        })
    return pd.DataFrame(rows).sort_values("key")


def reliability_grid(root: Path) -> pd.DataFrame:
    """BLUP slope reliability across the design grid.

    Reliability is 1 - (posterior SD / prior SD)^2 per subject, as written by
    the phenotype step.  ``effective_N`` is n x mean reliability -- the sample
    size a perfectly measured phenotype would need to carry the same
    information, which is what makes the min_visits tradeoff legible.
    """
    rows = []
    for r in _runs(root).itertuples():
        ph = pd.read_parquet(Path(r.dir) / "phenotypes" / "phenotypes.parquet")
        sl = ph[ph.phenotype == "slope"]
        if "reliability" not in sl.columns:
            continue
        rel = sl.groupby("subject")["reliability"].mean()
        mt = pd.read_parquet(Path(r.dir) / "model_table.parquet")
        rows.append({
            "key": _key(r), "global_cov": r.global_cov, "min_visits": r.min_visits,
            "family_effect": r.family_effect, "subjects": len(rel),
            "mean_visits": round(_mean_visits(mt), 4),
            "mean_reliability": round(rel.mean(), 4),
            "median_reliability": round(rel.median(), 4),
            "effective_N": int(round(len(rel) * rel.mean())),
        })
    return pd.DataFrame(rows).sort_values("key")


def _icc(y: pd.Series, g: pd.Series) -> dict:
    """One-way ICC of ``y`` by grouping ``g``, with its F test."""
    from scipy import stats
    df = pd.DataFrame({"y": y, "g": g}).dropna()
    grps = [v.values for _, v in df.groupby("g")["y"] if len(v) > 1]
    k, n = len(grps), len(df)
    F, p = stats.f_oneway(*grps)
    # ICC(1) from the one-way F: (F - 1) / (F + n0 - 1), with n0 the mean group
    # size.  Equivalent to (MSB - MSW) / (MSB + (n0 - 1) MSW) and less prone to
    # arithmetic slips than assembling the mean squares by hand.
    n0 = n / k
    icc = max(0.0, (F - 1) / (F + n0 - 1))
    return {"k": k, "n": n, "icc": round(icc, 5), "F": round(F, 5), "p": p}


def site_scanner_icc(root: Path) -> pd.DataFrame:
    run = str(root / "out" / SETTLED)
    y = _global_slope(run)
    mt = pd.read_parquet(Path(run) / "model_table.parquet")
    rows = []
    site = mt.groupby("subject")["site"].first()
    rows.append({"grouping": "site", **_icc(y, site.reindex(y.index))})
    try:
        sc = io.get_adapter("7.0").scanner()
        col = [c for c in sc.columns if c not in ("subject", "visit")][0]
        base = sc.sort_values("visit").groupby("subject")[col].first()
        rows.append({"grouping": f"scanner {col} (baseline)",
                     **_icc(y, base.reindex(y.index))})
    except Exception as exc:  # SourceUnavailable, or a schema change
        print(f"  (scanner ICC skipped: {type(exc).__name__}: {exc})")
    return pd.DataFrame(rows)


def _baseline_manufacturer(y_index) -> pd.Series:
    """Each subject's manufacturer at their first imaged visit.

    Manufacturer is the ONLY scanner grouping either release ships locally:
    ``mri_info_deviceserial`` is absent from both trees (verified against the
    release CSV headers, not a docstring), so serial-level scanner analyses are
    not reproducible here and this is the coarser substitute.  See the
    limitation stated in REPORT_7.0 section 9.
    """
    sc = io.get_adapter("7.0").scanner()
    col = [c for c in sc.columns if c not in ("subject", "visit")][0]
    return sc.sort_values("visit").groupby("subject")[col].first().reindex(y_index)


def site_scanner_icc_by_region(root: Path) -> pd.DataFrame:
    """Per-region ICC of the regional slope by site and by manufacturer.

    Previously this table's ``scanner_icc`` column was computed at serial
    level and its ``var_ratio_switch`` from serial switching; neither is
    reproducible (no serial column exists), so both are recomputed at
    manufacturer level here and the column is renamed to say so.
    """
    run = str(root / "out" / SETTLED)
    ph = pd.read_parquet(Path(run) / "phenotypes" / "phenotypes.parquet")
    sl = ph[ph.phenotype == "slope"]
    wide = sl.pivot_table(index="subject", columns="label", values="value")
    mt = pd.read_parquet(Path(run) / "model_table.parquet")
    site = mt.groupby("subject")["site"].first().reindex(wide.index)
    manu = _baseline_manufacturer(wide.index)

    rows = []
    for lab in wide.columns:
        y = wide[lab].dropna()
        rows.append({"label": lab,
                     "site_icc": _icc(y, site.reindex(y.index))["icc"],
                     "manufacturer_icc": _icc(y, manu.reindex(y.index))["icc"]})
    out = pd.DataFrame(rows)
    out["region"] = out.label.str.replace(r"^(lh|rh)_", "", regex=True)
    return out


def scanner_switching_summary(root: Path) -> pd.DataFrame:
    """How often manufacturer changes mid-study, and whether it adds variance.

    Serial-level switching (previously reported as 26.5% of subjects) is not
    reproducible; only manufacturer switching is.
    """
    run = str(root / "out" / SETTLED)
    y = _global_slope(run)
    mt = pd.read_parquet(Path(run) / "model_table.parquet")
    sc = io.get_adapter("7.0").scanner()
    col = [c for c in sc.columns if c not in ("subject", "visit")][0]
    sc = sc[sc.subject.isin(y.index)]
    nuniq = sc.groupby("subject")[col].nunique()
    switched = (nuniq > 1).reindex(y.index).fillna(False)

    v_sw = y[switched].var(ddof=1)
    v_no = y[~switched].var(ddof=1)
    d = ((y[switched].mean() - y[~switched].mean())
         / np.sqrt((v_sw + v_no) / 2))
    rows = [("subjects", float(len(y))),
            ("sites", float(mt.groupby("subject")["site"].first().nunique())),
            ("distinct manufacturers", float(sc[col].nunique())),
            ("manufacturer metadata coverage",
             float(nuniq.reindex(y.index).notna().mean())),
            ("pct switched manufacturer", float(switched.mean())),
            ("n switchers", float(switched.sum())),
            ("var ratio switch/same (global slope)", float(v_sw / v_no)),
            ("Cohen d switch vs same", float(d))]
    return pd.DataFrame(rows, columns=["quantity", "value"]).round(4)


def developmental_maps_noglobal(root: Path) -> pd.DataFrame:
    """The 68-region map table behind sections 5-7.

    ``hidden`` flags regions that the bundled ggseg schematic either cannot draw
    at all or draws as an unreadable sliver; it is sourced from
    ``brainplot.UNSUPPORTED_REGIONS`` / ``SLIVER_REGIONS`` rather than
    hardcoded, so the figure code and the table cannot drift apart.
    """
    from abcd import brainplot as bp

    run = str(root / "out" / SETTLED)
    mp = M.regional_maps(run).copy()

    # Baseline thickness: the model intercept, i.e. mm at the centring age (not
    # at the first visit).  Section 5 shows it beside the rate maps, and it is
    # the positive-control phenotype for the GWAS, so it belongs in the same
    # table rather than being recomputed per figure.
    fx = M.load_fits(run, "fits")["fixed"]
    icept = (fx[fx.term == "(Intercept)"].set_index("label")["estimate"]
             .rename("baseline_thickness"))
    mp = mp.join(icept, how="left")
    if mp.baseline_thickness.isna().any():
        raise SystemExit("missing intercept for some regions")
    W = cov.slope_matrix(run)
    # slope_pcs names its components PC1..PC3; sc_pcs names them SC1..SC3.
    # Published column names are slopePC1.. and scPC1.., so normalise both.
    for prefix, load in (("slope", cov.slope_pcs(W, 3).loadings),
                         ("sc", cov.sc_pcs(cov.sc_matrix(W), 3).loadings)):
        load = load.rename(columns=lambda c: f"{prefix}PC{c[-1]}")
        mp = mp.join(load, how="left")

    ph = pd.read_parquet(Path(run) / "phenotypes" / "phenotypes.parquet")
    rel = ph[ph.phenotype == "slope"].groupby("label")["reliability"].mean()
    mp["reliability"] = rel.reindex(mp.index)

    mp["region"] = mp.index.str.replace(r"^[lr]h_", "", regex=True)
    mp["lobe"] = mp.region.map(_LOBE_OF)
    if mp.lobe.isna().any():
        raise SystemExit(f"unmapped regions: {sorted(mp.loc[mp.lobe.isna(), 'region'].unique())}")
    mp["hidden"] = mp.region.isin(bp.UNSUPPORTED_REGIONS + bp.SLIVER_REGIONS)

    h2 = her.falconer_by_region(ph, her.pair_table("7.0"))
    mp = mp.join(h2[["h2", "r_MZ", "r_DZ"]].rename(
        columns={"r_MZ": "h2_rMZ", "r_DZ": "h2_rDZ"}), how="left")
    return mp.reset_index()


def ahba_vs_maps_noglobal(root: Path) -> pd.DataFrame:
    """Spin-tested correlation of every developmental map with AHBA C1-C3.

    This table previously had NO generator: it was written once and then read by
    the section 7 figure and by ``gwas_phenotype_priority``.  It went stale --
    its h2 row held rho = +0.342 against C1, while every h2 map in the repo
    gives +0.313, so the value came from a superseded h2 definition and nothing
    could detect the drift.  Regenerating it from the maps table closes that.

    The spin test is the expensive part (1000 rotations x 3 components x 6 maps),
    which is presumably why it was cached; it is a few minutes, not hours, and
    correctness is worth more than the cache.
    """
    from abcd import genemaps, spatial

    mp = pd.read_csv(root / "docs" / "developmental_maps_noglobal.csv").set_index("label")
    geom = spatial.load_dk_geometry()
    comps = genemaps.ahba_components("dsk")

    # Labels are the published column names; keep the section-5/7 ordering.
    wanted = [("Absolute thinning rate", "slope_total"),
              ("Between-subject SD of rate", "tau_slope"),
              ("Slope PC1", "slopePC1"),
              ("Slope PC2", "slopePC2"),
              ("Slope PC3", "slopePC3"),
              ("Regional h²", "h2")]
    rows = []
    for label, col in wanted:
        if col not in mp.columns:
            continue
        r = genemaps.map_vs_components(mp[col].dropna(), geom, components=comps)
        r.insert(0, "map_col", col)
        r.insert(0, "map", label)
        rows.append(r)
    return pd.concat(rows, ignore_index=True)


def handoff_release_comparison(root: Path) -> pd.DataFrame:
    """Cross-release power comparison, matched on design.

    Every row must come from a run with the SAME global-covariate, family-effect
    and site-effect settings, differing only in release and visit filter -- the
    previous hand-maintained version of this table mixed a global-adjusted 5.1
    row with no-global 7.0 rows, which made the >=3-visit filter appear to lose
    to 5.1 when matched runs show it winning.  Rows whose design does not match
    the settled 7.0 specification are refused rather than silently reported.
    """
    want = {"global_cov": "none", "family_effect": False}
    picks = [("5.1", 2), ("7.0", 2), ("7.0", 3), ("7.0", 4)]
    rows = []
    for rel, mv in picks:
        hit = None
        for d in sorted((root / "out").glob(f"thickness_dsk_{rel.replace('.', '')}_*")):
            if not (d / "phenotypes" / "phenotypes.parquet").exists():
                continue
            c = yaml.safe_load((d / "config.yaml").read_text())
            g = c.get("global_covariate")
            if (g in (None, "none", "") ) != (want["global_cov"] == "none"):
                continue
            if bool(c.get("family_effect")) != want["family_effect"]:
                continue
            if int(c.get("min_visits")) != mv:
                continue
            hit = d
            break
        if hit is None:
            raise SystemExit(
                f"no fitted {rel} run with min_visits={mv}, global_covariate=none, "
                f"family_effect=false -- this table must be matched on design; "
                f"fit it (configs/ct_51_noglobal_mv2_matched.yaml for 5.1) first"
            )
        ph = pd.read_parquet(hit / "phenotypes" / "phenotypes.parquet")
        sl = ph[ph.phenotype == "slope"]
        mt = pd.read_parquet(hit / "model_table.parquet")
        rel_mean = float(sl.groupby("subject")["reliability"].mean().mean())
        n = int(sl.subject.nunique())
        rows.append({"run": f"{rel} (>={mv} visits)" if mv < 4 else f"{rel} (4 visits)",
                     "release": rel, "min_visits": mv, "run_id": hit.name,
                     "subjects": n, "mean_visits": round(_mean_visits(mt), 4),
                     "mean_reliability": round(rel_mean, 5),
                     "effective_N": round(rel_mean * n, 1)})
    df = pd.DataFrame(rows)
    df["vs_5.1"] = (df.effective_N / df.effective_N.iloc[0]).round(4)
    return df


def regional_slope_pc_loadings(root: Path) -> pd.DataFrame:
    W = cov.slope_matrix(str(root / "out" / SETTLED))
    return cov.slope_pcs(W, 3).loadings.rename_axis("label").reset_index()


TABLES = {
    "fit_summary": fit_summary,
    "reliability_grid": reliability_grid,
    "site_scanner_icc": site_scanner_icc,
    "site_scanner_icc_by_region": site_scanner_icc_by_region,
    "scanner_switching_summary": scanner_switching_summary,
    "developmental_maps_noglobal": developmental_maps_noglobal,
    "ahba_vs_maps_noglobal": ahba_vs_maps_noglobal,
    "regional_slope_pc_loadings": regional_slope_pc_loadings,
    "handoff_release_comparison": handoff_release_comparison,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="docs")
    ap.add_argument("--only", nargs="*", choices=sorted(TABLES), default=None)
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = root / a.outdir
    for name in (a.only or sorted(TABLES)):
        df = TABLES[name](root)
        p = out / f"{name}.csv"
        df.to_csv(p, index=False, float_format="%.6g")
        print(f"wrote {p.relative_to(root)}  ({len(df)} rows)")


if __name__ == "__main__":
    main()
