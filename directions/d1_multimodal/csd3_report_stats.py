"""
Numbers for CSD3_REPORT.md, from the collected surf_metrics tables (run after tools/surf_metrics_collect.py).

Per measure: sessions processed (overall and per visit), failures by reason, whole-cortex
mean by visit (median, IQR), and the share of variance in the whole-cortex mean explained
by site and by scanner manufacturer (one-way R^2, ab_g_dyn), with thickness as reference.
Also lists every output file with its size.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from abcd import paths

REL = paths.release_dir("7.0")
HCP, DSK = REL / "processed" / "hcp", REL / "processed" / "dsk_local"
KEY = ["participant_id", "session_id"]


def md(df: pd.DataFrame, index: bool = True, floatfmt: str = ".3f") -> str:
    """Markdown table without the tabulate dependency."""
    d = df.reset_index() if index else df
    fmt = lambda x: format(x, floatfmt) if isinstance(x, float) else str(x)
    lines = ["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * d.shape[1]]
    lines += ["| " + " | ".join(fmt(x) for x in row) + " |" for row in d.itertuples(index=False)]
    return "\n".join(lines)


def r2(y: pd.Series, g: pd.Series) -> float:
    ok = y.notna() & g.notna()
    y, g = y[ok], g[ok]
    fit = y.groupby(g).transform("mean")
    return float(1 - ((y - fit) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def main() -> int:
    qc = pd.read_csv(HCP / "surf_metrics_session_qc.tsv", sep="\t")
    qc["reason"] = qc.status.str.split(":").str[0]
    print("## sessions by measure and status\n")
    print(md(qc.pivot_table(index="measure", columns="reason", values="session_id",
                            aggfunc="count", fill_value=0)))
    bad = qc[qc.status != "ok"]
    if len(bad):
        print("\nfailing sessions:\n")
        print(md(bad[[*KEY, "measure", "status", "missing_input"]], index=False))
    print("\n## ok sessions by visit\n")
    print(md(qc[qc.status == "ok"].pivot_table(index="measure", columns="session_id",
                                               values="participant_id", aggfunc="count")))
    print(f"\nannotation source: {qc.drop_duplicates(KEY).annot_source.value_counts().to_dict()}")
    agree = qc[["normal_agree_lh", "normal_agree_rh"]].min(axis=1)
    print(f"white-normal/white->pial agreement: median {agree.median():.4f}, min {agree.min():.4f}")

    dyn = pd.read_csv(REL / "g" / "dyn" / "ab_g_dyn.tsv", sep="\t",
                      usecols=[*KEY, "ab_g_dyn__design_site", "ab_g_dyn__design_mr__manufact"])
    stems = {f"{m}__{a}": (HCP if a == "hcp" else DSK) / f"mr_y_smri__{m}__{a}.tsv"
             for m in sorted(qc.measure.unique()) for a in ("hcp", "dsk")}
    stems["thk__hcp (reference)"] = HCP / "mr_y_smri__thk__hcp.tsv"
    rows, vis = [], []
    for name, path in stems.items():
        col = f"mr_y_smri__{name.split(' ')[0]}_mean"
        t = pd.read_csv(path, sep="\t", usecols=[*KEY, col]).merge(dyn, on=KEY, how="left")
        rows.append({"measure": name, "n": int(t[col].notna().sum()),
                     "R2 site": r2(t[col], t.ab_g_dyn__design_site),
                     "R2 manufacturer": r2(t[col], t.ab_g_dyn__design_mr__manufact),
                     "n missing site": int(t.ab_g_dyn__design_site.isna().sum())})
        g = t.groupby("session_id")[col]
        vis.append(pd.DataFrame({"measure": name, "median": g.median(), "q25": g.quantile(.25),
                                 "q75": g.quantile(.75), "n": g.count()}).reset_index())
    print("\n## whole-cortex mean by visit (median [IQR])\n")
    v = pd.concat(vis)
    v["value"] = v.apply(lambda r: f"{r['median']:.2f} [{r.q25:.2f}-{r.q75:.2f}]", axis=1)
    print(md(v.pivot(index="measure", columns="session_id", values="value")))
    print("\n## variance in whole-cortex mean explained by site / manufacturer\n")
    print(md(pd.DataFrame(rows), index=False))

    print("\n## output files\n")
    files = sorted([*HCP.glob("mr_y_smri__*__hcp.tsv"), *DSK.glob("mr_y_smri__*__dsk.tsv"),
                    HCP / "surf_metrics_session_qc.tsv", HCP / "PROVENANCE_surf_metrics.txt",
                    HCP / "surf_metrics_inventory.tsv"])
    files = [f for f in files if not f.name.startswith("mr_y")
             or f.name.split("__")[1] in set(qc.measure.str.split("__").str[0])]
    for f in files:
        print(f"- `{f}` ({f.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
