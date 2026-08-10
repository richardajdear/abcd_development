#!/usr/bin/env python
"""Regenerate every Falconer heritability table in ``docs/`` from runs on disk.

Why this script exists
----------------------
The ``docs/h2_*.csv`` and ``docs/*heritability*.csv`` tables were originally
built ad hoc in a session, which meant (a) nobody could reproduce them and (b)
when the pairing was corrected, there was no single place to re-run.  Every
Falconer number quoted in ``docs/REPORT_7.0.md`` now comes from here.

What changed in the numbers
---------------------------
The tables written before 2026-07-30 used zygosity borrowed from the 5.1 pi-hat
file, pairs inferred from family size, and DZ twins pooled with non-twin full
siblings.  7.0 ships its own pi-hat table (``y/gn_y_genrel.tsv``), which is read
by ``io.Release70Adapter.genotyped_pairs``.  Three consequences:

* Pairs are explicit, so the usable counts rise (MZ 260 -> 266,
  DZ+sibling 871 -> 967 on the settled run).
* DZ twins separate from non-twin siblings, and they correlate substantially
  more (0.279 vs 0.090 on the global slope).  Pooling therefore inflated h2.
* Headline global-slope h2 falls from 0.58 to 0.44.  Rank order across
  candidate phenotypes is unchanged; see ``docs/REPORT_7.0.md`` section 9.

Usage
-----
    PYTHONPATH=src python tools/regen_h2_tables.py [--outdir docs] [--boot 2000]

Every table is rewritten in place.  The script is deterministic given a seed.
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

#: The settled design (report section 2): no global covariate, >= 2 visits, no
#: family random effect.  Everything else is a sensitivity analysis.
SETTLED = "thickness_dsk_70_139406217085"
SEED = 0


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


def _subject_maps(run_dir: str) -> dict[str, pd.Series]:
    """Candidate subject-level phenotypes from one run.

    Delegates to :func:`abcd.gcta_export.subject_phenotypes` so that the vector
    whose heritability this table reports is the identical vector exported to
    the cluster for GWAS.  The definitions used to be written out twice, here
    and in the export module, which is one refactor away from publishing an h2
    for a phenotype nobody ran a GWAS on.  Display names are this table's
    convention and are mapped from the export's column names.
    """
    from abcd.gcta_export import PHENOTYPES, subject_phenotypes

    frame = subject_phenotypes(run_dir)
    return {p["display"]: frame[p["name"]] for p in PHENOTYPES}


def _falc_row(y: pd.Series, pairs: pd.DataFrame, boot: int, **extra) -> dict:
    r = her.falconer(y, pairs, n_boot=boot, seed=SEED)
    return {**extra, "h2": r["h2"], "ci_lo": r.get("ci_lo"), "ci_hi": r.get("ci_hi"),
            "r_MZ": r["r_MZ"], "r_DZ": r["r_DZ"], "n_MZ": r["n_MZ"], "n_DZ": r["n_DZ"],
            "dz_class": r["dz_class"], "r_sib_excluded": r.get("r_sib_excluded")}


def candidate_phenotypes(root: Path, pairs, boot: int) -> pd.DataFrame:
    ys = _subject_maps(str(root / "out" / SETTLED))
    rows = [_falc_row(y, pairs, boot, phenotype=k) for k, y in ys.items()]
    return pd.DataFrame(rows)


def dz_class_comparison(root: Path, pairs, boot: int) -> pd.DataFrame:
    """The correction itself, as a table: pooled vs twins-only, both reported.

    This is the evidence for defaulting to DZ twins.  ``delta_ci`` is a paired
    bootstrap on the *difference* of the two DZ correlations, which is the
    quantity that decides whether pooling is defensible.
    """
    ys = _subject_maps(str(root / "out" / SETTLED))
    rng = np.random.default_rng(SEED)
    rows = []
    for name, y in ys.items():
        y = y[~y.index.duplicated()]
        mats = {c: her._pair_matrix(y, pairs, c)
                for c in (her.MZ, her.DZ_TWIN, her.SIB)}
        from scipy.stats import pearsonr
        r = {c: pearsonr(A[:, 0], A[:, 1]).statistic for c, A in mats.items()}
        A, B = mats[her.DZ_TWIN], mats[her.SIB]
        draws = np.empty(boot or 2000)
        for b in range(len(draws)):
            a = A[rng.integers(0, len(A), len(A))]
            bb = B[rng.integers(0, len(B), len(B))]
            draws[b] = (pearsonr(a[:, 0], a[:, 1]).statistic
                        - pearsonr(bb[:, 0], bb[:, 1]).statistic)
        rows.append({
            "phenotype": name, "r_MZ": r[her.MZ], "r_DZ_twin": r[her.DZ_TWIN],
            "r_full_sib": r[her.SIB], "n_MZ": len(mats[her.MZ]),
            "n_DZ_twin": len(A), "n_full_sib": len(B),
            "delta_DZtwin_minus_sib": r[her.DZ_TWIN] - r[her.SIB],
            "delta_ci_lo": np.percentile(draws, 2.5),
            "delta_ci_hi": np.percentile(draws, 97.5),
            "p_delta_le_0": float((draws <= 0).mean()),
            "h2_twins_only": her.falconer(y, pairs, dz_class="twins_only")["h2"],
            "h2_pooled": her.falconer(y, pairs, dz_class="pooled")["h2"],
        })
    return pd.DataFrame(rows)


def global_slope_by_design(root: Path, pairs, boot: int) -> pd.DataFrame:
    """Sensitivity of the headline estimate to the four design switches."""
    rows = []
    for _, r in _runs(root).iterrows():
        ph = pd.read_parquet(Path(r["dir"]) / "phenotypes" / "phenotypes.parquet")
        y = ph[ph.phenotype == "slope"].groupby("subject")["value"].mean()
        rows.append(_falc_row(y, pairs, boot, run=r["run"], global_cov=r["global_cov"],
                              min_visits=r["min_visits"], family_effect=r["family_effect"]))
    return pd.DataFrame(rows).sort_values(["global_cov", "min_visits", "family_effect"])


def family_effect_contrast(root: Path, pairs, boot: int) -> pd.DataFrame:
    """The positive control that exposed the family-random-effect artefact.

    Baseline thickness is strongly heritable, so a design that returns h2 > 1
    for it is broken.  ``family_effect: true`` does exactly that, because the
    random effect absorbs between-family variance and siblings' BLUPs then
    anti-correlate.
    """
    runs = _runs(root)
    sel = runs[(runs.global_cov == "none") & (runs.min_visits == 2)]
    rows = []
    for _, r in sel.iterrows():
        for name, y in _subject_maps(r["dir"]).items():
            if name not in ("global mean slope", "baseline thickness (control)"):
                continue
            rows.append(_falc_row(y, pairs, boot, family_effect=r["family_effect"],
                                  phenotype=name))
    return pd.DataFrame(rows).sort_values(["phenotype", "family_effect"])


def site_scanner_sensitivity(root: Path, pairs, boot: int) -> pd.DataFrame:
    """Does removing site/scanner structure change h2?

    If the heritable signal were really site or scanner artefact shared within
    families (twins are scanned at the same site, usually the same session), then
    residualising it out would collapse h2.  It does not.
    """
    run_dir = root / "out" / SETTLED
    ph = pd.read_parquet(run_dir / "phenotypes" / "phenotypes.parquet")
    mt = pd.read_parquet(run_dir / "model_table.parquet")
    y = ph[ph.phenotype == "slope"].groupby("subject")["value"].mean()

    keys = {}
    if "site" in mt.columns:
        keys["site-mean removed"] = mt.drop_duplicates("subject").set_index("subject")["site"]
    # Scanner: the local 7.0 tree ships manufacturer only, not a device serial
    # (5.1 had ``mri_info_deviceserial``).  Manufacturer is a 2-3 level grouping
    # rather than ~29 scanners, so it removes much less variance -- the row is
    # labelled to make that visible instead of implying a serial-level control.
    try:
        sc = io.get_adapter("7.0").scanner()
        col = next((c for c in ("scanner_serial", "device_serial", "scanner_id",
                                "manufacturer") if c in sc.columns), None)
        if col is None:
            print(f"  (scanner sensitivity skipped: no usable column in {list(sc.columns)})")
        else:
            first = (sc.sort_values("visit").drop_duplicates("subject")
                       .set_index("subject")[col])
            label = ("baseline-scanner-mean removed" if col != "manufacturer"
                     else f"baseline-manufacturer-mean removed (k={first.nunique()})")
            keys[label] = first
    except Exception as exc:  # SourceUnavailable, or a schema change
        print(f"  (scanner sensitivity skipped: {type(exc).__name__}: {exc})")

    rows = [_falc_row(y, pairs, boot, phenotype="global slope, raw", n=len(y))]
    for label, key in keys.items():
        g = y.to_frame("y").join(key.rename("k")).dropna()
        resid = g.y - g.groupby("k").y.transform("mean")
        rows.append(_falc_row(resid, pairs, boot, phenotype=label, n=len(resid)))
    return pd.DataFrame(rows)


def top_h2_selection_optimism(root: Path, pairs, boot: int,
                              n_splits: int = 40, k: int = 10) -> pd.DataFrame:
    """In-sample vs held-out h2 for a "top-k most heritable regions" composite.

    Selecting the k regions with highest h2 and then reporting the composite's h2
    in the same sample is circular.  This quantifies the gap, and how unstable
    the selected set is across splits.
    """
    run_dir = str(root / "out" / SETTLED)
    ph = pd.read_parquet(Path(run_dir) / "phenotypes" / "phenotypes.parquet")
    W = cov.slope_matrix(run_dir)
    full_top = (her.falconer_by_region(ph, pairs, phenotype="slope")
                .h2.nlargest(k).index)

    fams = pairs.family_id.dropna().unique()
    rng = np.random.default_rng(SEED)
    rows = []
    for s in range(n_splits):
        tr_fam = set(rng.choice(fams, size=len(fams) // 2, replace=False))
        tr, te = pairs[pairs.family_id.isin(tr_fam)], pairs[~pairs.family_id.isin(tr_fam)]
        tr_subj = set(tr.a) | set(tr.b)
        h2_tr = her.falconer_by_region(ph[ph.subject.isin(tr_subj)], tr,
                                       phenotype="slope", min_pairs=20)
        top = h2_tr.h2.nlargest(k).index
        y = W[list(top)].mean(axis=1)
        try:
            rows.append({"split": s,
                         "in_sample": her.falconer(y, tr)["h2"],
                         "held_out": her.falconer(y, te)["h2"],
                         f"n_overlap_with_full_top{k}": len(set(top) & set(full_top))})
        except ValueError:
            continue
    return pd.DataFrame(rows)


def regional(root: Path, pairs, boot: int) -> pd.DataFrame:
    ph = pd.read_parquet(Path(root / "out" / SETTLED) / "phenotypes" / "phenotypes.parquet")
    slope = her.falconer_by_region(ph, pairs, phenotype="slope")
    inter = her.falconer_by_region(ph, pairs, phenotype="intercept")
    out = slope.rename(columns={"h2": "h2_slope"}).join(
        inter[["h2"]].rename(columns={"h2": "h2_intercept"}))
    rel = ph[ph.phenotype == "slope"].groupby("label")["reliability"].mean()
    return out.join(rel.rename("reliability")).reset_index()


def _held_out_splits(root: Path, pairs: pd.DataFrame, n_splits: int = 40,
                     seed: int = SEED) -> pd.DataFrame:
    """Held-out h2 for phenotypes whose *definition* is data-dependent.

    Slope PCs and "top-k h2 regions" are chosen using the data, so their
    in-sample h2 is optimistic.  Each split defines the phenotype on a training
    half and estimates h2 on the held-out half.

    Splitting is by **pair**, not by subject: a pair straddling the two halves
    would leak the very covariance Falconer measures, and would also be dropped
    from both halves' pair matrices, silently shrinking n.  Families are kept
    whole for the same reason.
    """
    run_dir = str(root / "out" / SETTLED)
    ph = pd.read_parquet(Path(run_dir) / "phenotypes" / "phenotypes.parquet")
    W = cov.slope_matrix(run_dir)
    gs = ph[ph.phenotype == "slope"].groupby("subject")["value"].mean()

    fams = pairs.family_id.dropna().unique()
    rng = np.random.default_rng(seed)
    rows = []
    for s in range(n_splits):
        tr_fam = set(rng.choice(fams, size=len(fams) // 2, replace=False))
        tr = pairs[pairs.family_id.isin(tr_fam)]
        te = pairs[~pairs.family_id.isin(tr_fam)]
        tr_subj = set(tr.a) | set(tr.b)

        # phenotypes defined on train only
        Wtr = W.loc[W.index.intersection(sorted(tr_subj))]
        if len(Wtr) < 50:
            continue
        load = cov.slope_pcs(Wtr, 3).loadings
        scores = cov.subject_scores(W, load)          # applied to everyone
        h2_tr = her.falconer_by_region(ph[ph.subject.isin(tr_subj)], tr,
                                       phenotype="slope", min_pairs=20)
        top10 = h2_tr.h2_slope.nlargest(10).index if "h2_slope" in h2_tr else \
            h2_tr.h2.nlargest(10).index
        topmean = W[list(top10)].mean(axis=1)

        cand = {"global mean slope": ("none", gs),
                "slope PC1": ("train (PCA)", scores.PC1),
                "slope PC2": ("train (PCA)", scores.PC2),
                "slope PC3": ("train (PCA)", scores.PC3),
                "top-10 h2 regions (mean)": ("train (h2)", topmean)}
        for name, (defined_on, y) in cand.items():
            try:
                r = her.falconer(y, te)
            except ValueError:
                continue
            rows.append({"split": s, "phenotype": name, "defined_on": defined_on,
                         "h2": r["h2"], "r_MZ": r["r_MZ"], "r_DZ": r["r_DZ"],
                         "n_MZ": r["n_MZ"], "n_DZ": r["n_DZ"]})
    return pd.DataFrame(rows)


def held_out_h2(root: Path, pairs, boot: int) -> pd.DataFrame:
    return _held_out_splits(root, pairs)


def heldout_h2_paired_tests(root: Path, pairs, boot: int) -> pd.DataFrame:
    """Each candidate against the global mean slope, on identical splits.

    Reported under both DZ classes, because the choice reverses one conclusion.
    The pre-correction tables ranked a "top-10 most heritable regions" composite
    above the global mean slope.  It does not survive: with DZ twins alone the
    composite is *worse* than the global slope, and even with siblings pooled
    back in it only ties.  The composite also selects regions on h2 estimated in
    the training half, so its apparent advantage was partly selection optimism
    that pair-level splitting removes.
    """
    from scipy import stats
    ref = "global mean slope"
    out = []
    for label, kw in (("twins_only", {}), ("pooled", {"dz_class": "pooled"})):
        orig = her.falconer
        if kw:
            her.falconer = lambda y, p, **k: orig(y, p, **{**k, **kw})  # noqa: E731
        try:
            piv = (_held_out_splits(root, pairs)
                   .pivot_table(index="split", columns="phenotype", values="h2"))
        finally:
            her.falconer = orig
        for c in piv.columns:
            if c == ref:
                continue
            d = (piv[c] - piv[ref]).dropna()
            t = stats.ttest_rel(piv.loc[d.index, c], piv.loc[d.index, ref])
            out.append({"dz_class": label, "comparison": f"{c} - {ref}",
                        "delta": d.mean(), "sd": d.std(), "t": t.statistic,
                        "p": t.pvalue, "n_splits": len(d)})
    return pd.DataFrame(out).sort_values(["dz_class", "delta"], ascending=[True, False])


def gwas_phenotype_priority(root: Path, pairs, boot: int) -> pd.DataFrame:
    """Ranked shortlist of phenotypes to take to the cluster.

    Columns are the three things that decide cluster time: how heritable the
    phenotype is when its definition cannot cheat (``h2_heldout``), how well it
    matches the transcriptional axes we predicted (``ahba_*``), and whether the
    definition is reproducible outside this sample (``definition_stable``).
    """
    hs = _held_out_splits(root, pairs)
    agg = (hs.groupby(["phenotype", "defined_on"])
             .agg(h2_heldout=("h2", "mean"), h2_heldout_sd=("h2", "std"),
                  n_splits=("h2", "size")).reset_index())
    full = candidate_phenotypes(root, pairs, boot=0).set_index("phenotype")["h2"]
    agg["h2_full"] = agg.phenotype.map(full)

    ahba_p = root / "docs" / "ahba_vs_maps_noglobal.csv"
    if ahba_p.exists():
        A = pd.read_csv(ahba_p)
        key = {"global mean slope": "slope_total", "slope PC1": "slopePC1",
               "slope PC2": "slopePC2", "slope PC3": "slopePC3"}
        best = []
        for nm in agg.phenotype:
            sub = A[A.map_col == key.get(nm, "")]
            if len(sub):
                b = sub.loc[sub.rho.abs().idxmax()]
                best.append((b.component, b.rho, b.p_spin))
            else:
                best.append((None, np.nan, np.nan))
        agg[["ahba_best", "ahba_rho", "ahba_p_spin"]] = pd.DataFrame(best, index=agg.index)

    agg["definition_uses_data"] = agg.defined_on != "none"
    return agg.sort_values("h2_heldout", ascending=False)[
        ["phenotype", "defined_on", "h2_full", "h2_heldout", "h2_heldout_sd", "n_splits",
         "ahba_best", "ahba_rho", "ahba_p_spin", "definition_uses_data"]]


TABLES = {
    "h2_candidate_phenotypes": candidate_phenotypes,
    "h2_dz_class_comparison": dz_class_comparison,
    "h2_global_slope_by_design": global_slope_by_design,
    "h2_family_effect_contrast": family_effect_contrast,
    "h2_site_scanner_sensitivity": site_scanner_sensitivity,
    "topH2_selection_optimism": top_h2_selection_optimism,
    "regional_h2": regional,
    "heldout_h2_by_split": held_out_h2,
    "heldout_h2_paired_tests": heldout_h2_paired_tests,
    "gwas_phenotype_priority": gwas_phenotype_priority,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="repo root (default: cwd)")
    ap.add_argument("--outdir", default="docs")
    ap.add_argument("--boot", type=int, default=2000, help="bootstrap resamples")
    ap.add_argument("--only", nargs="*", choices=sorted(TABLES), default=None)
    a = ap.parse_args()

    root = Path(a.root).resolve()
    outdir = root / a.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    pairs = her.pair_table("7.0")
    counts = pairs.pair_type.value_counts().to_dict()
    print(f"pair table: {len(pairs)} pairs {counts}")

    for name in (a.only or TABLES):
        df = TABLES[name](root, pairs, a.boot)
        p = outdir / f"{name}.csv"
        df.to_csv(p, index=False, float_format="%.6g")
        print(f"wrote {p.relative_to(root)}  ({len(df)} rows)")


if __name__ == "__main__":
    main()
