#!/usr/bin/env python
"""Compare the settled specification fitted on two tabulations of ABCD.

Written for the 2026-09-14 move from the 6.0 tables (analysed under a 7.0
label) to the true 7.0 tabulation.  Reads two run directories with the same
config hash and writes one long table, ``docs/vintage_comparison.csv``, that
the README and docs/RERUN_7.0_TABULATED.md quote from.

    PYTHONPATH=src python tools/compare_vintage.py \
        --old out/legacy_6.0_tabulated/thickness_dsk_70_139406217085 \
        --new out/thickness_dsk_70_139406217085

Everything here is descriptive; no number is hand-entered anywhere else.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]


def _load(run: Path) -> dict:
    m = json.loads((run / "manifest.json").read_text())
    mt = pd.read_parquet(run / "model_table.parquet", columns=["subject", "visit", "family_id"])
    ph = pd.read_parquet(run / "phenotypes" / "phenotypes.parquet")
    fx = pd.read_parquet(run / "fits" / "fixed.parquet")
    rel = pd.read_csv(run / "phenotypes" / "phenotypes_reliability.csv")
    ledger = pd.read_csv(run / "qc_ledger.csv")
    return dict(m=m, mt=mt, ph=ph, fx=fx, rel=rel, ledger=ledger)


def _global_slope(ph: pd.DataFrame) -> pd.Series:
    s = ph[ph.phenotype == "slope"]
    return s.groupby("subject")["value"].mean()


def _rows(tag: str, d: dict) -> list[dict]:
    m, mt, ph, rel, led = d["m"], d["mt"], d["ph"], d["rel"], d["ledger"]
    fin = m["final"]
    sv = mt.drop_duplicates(["subject", "visit"])
    nvis = sv.groupby("subject").size()
    slope_rel = ph[ph.phenotype == "slope"]
    rows = [
        ("scans_with_imaging", int(led.loc[led.criterion == "initial", "n_scans_in"].iloc[0])),
        ("scans_after_qc", int(led.n_scans_out.iloc[-1])),
        ("scans_dropped_release_flag", int(led.loc[led.criterion == "release_qc_include", "n_dropped"].iloc[0])),
        ("scans_dropped_defects", int(led.loc[led.criterion == "surface_defects", "n_dropped"].iloc[0])),
        ("scans_dropped_philips", int(led.loc[led.criterion == "scanner_manufacturer", "n_dropped"].iloc[0])),
        ("subjects_ge2_visits", fin["n_subjects"]),
        ("scans_in_model", fin["n_scans"]),
        ("families", fin["n_families"]),
        ("mean_visits_per_subject", round(float(nvis.mean()), 4)),
        ("subjects_2_visits", int((nvis == 2).sum())),
        ("subjects_3_visits", int((nvis == 3).sum())),
        ("subjects_4_visits", int((nvis == 4).sum())),
        ("scans_v0", fin["visit_counts"].get("v0")),
        ("scans_v2", fin["visit_counts"].get("v2")),
        ("scans_v4", fin["visit_counts"].get("v4")),
        ("scans_v6", fin["visit_counts"].get("v6")),
        ("age_centre_years", m.get("age_centre")),
        ("age_max_years", fin["age_range"][1]),
        ("slope_reliability_median", round(float(slope_rel.reliability.median()), 4)),
        ("slope_reliability_mean", round(float(slope_rel.reliability.mean()), 4)),
        ("effective_n_slope", round(float(slope_rel.groupby("subject").reliability.mean().sum()), 1)),
        ("group_mean_thinning_mm_per_yr",
         round(float(d["fx"].loc[d["fx"].term == "age_c", "estimate"].mean()), 5)),
    ]
    return [{"quantity": k, "vintage": tag, "value": v} for k, v in rows]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--old", required=True, type=Path)
    ap.add_argument("--new", required=True, type=Path)
    ap.add_argument("--out", default=ROOT / "docs" / "vintage_comparison.csv", type=Path)
    a = ap.parse_args()
    old, new = _load(a.old), _load(a.new)
    assert old["m"]["config_hash"] == new["m"]["config_hash"], "different specifications"

    rows = _rows("6.0_tables", old) + _rows("7.0_tables", new)

    # group-level maps: are the 68 regional thinning rates the same map?
    fo = old["fx"].query("term == 'age_c'").set_index("label").estimate
    fn = new["fx"].query("term == 'age_c'").set_index("label").estimate
    common = fo.index.intersection(fn.index)
    rows += [
        {"quantity": "map_spearman_old_vs_new_slope", "vintage": "both", "value": round(spearmanr(fo[common], fn[common])[0], 4)},
        {"quantity": "map_pearson_old_vs_new_slope", "vintage": "both", "value": round(pearsonr(fo[common], fn[common])[0], 4)},
    ]
    io_ = old["fx"].query("term == '(Intercept)'").set_index("label").estimate
    in_ = new["fx"].query("term == '(Intercept)'").set_index("label").estimate
    rows.append({"quantity": "map_spearman_old_vs_new_intercept", "vintage": "both",
                 "value": round(spearmanr(io_[common], in_[common])[0], 4)})

    # subject-level global slope on the shared subjects
    go, gn = _global_slope(old["ph"]), _global_slope(new["ph"])
    shared = go.index.intersection(gn.index)
    rows += [
        {"quantity": "subjects_shared", "vintage": "both", "value": int(len(shared))},
        {"quantity": "subjects_new_only", "vintage": "both", "value": int(len(gn.index.difference(go.index)))},
        {"quantity": "subjects_old_only", "vintage": "both", "value": int(len(go.index.difference(gn.index)))},
        {"quantity": "global_slope_pearson_shared_subjects", "vintage": "both",
         "value": round(pearsonr(go[shared], gn[shared])[0], 4)},
    ]
    # how many shared subjects gained a visit
    vo = old["mt"].drop_duplicates(["subject", "visit"]).groupby("subject").size()
    vn = new["mt"].drop_duplicates(["subject", "visit"]).groupby("subject").size()
    gained = (vn.reindex(shared) - vo.reindex(shared))
    rows.append({"quantity": "shared_subjects_gaining_a_visit", "vintage": "both", "value": int((gained > 0).sum())})

    out = pd.DataFrame(rows)
    out["old_run"], out["new_run"] = str(a.old.resolve().relative_to(ROOT)), str(a.new.resolve().relative_to(ROOT))
    a.out.parent.mkdir(exist_ok=True)
    out.to_csv(a.out, index=False)
    wide = out[out.vintage != "both"].pivot(index="quantity", columns="vintage", values="value")
    print(wide.to_string())
    print(out[out.vintage == "both"][["quantity", "value"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
