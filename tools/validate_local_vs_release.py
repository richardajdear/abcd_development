"""
Are the FreeSurfer derivatives on rds the surfaces the release tabulated?

Three comparisons, written as a Markdown report plus TSVs to
``docs/hcp_census/``:

1. **DK local vs release.**  ``abcd-data-release-7.0/processed/dsk_local/mr_y_smri__thk__dsk.tsv``
   (parsed from ``stats/?h.aparc.stats`` by :mod:`abcd.dk_stats`) against the
   release ``mr_y_smri__thk__dsk`` table, joined on participant x session.
   Reports the fraction of values equal at the release's printed precision,
   per-region max |delta|, the hemisphere means, and sessions present on one
   side only.
2. **HCP backfill vs July run** on the deliberately re-run overlap sample
   (``legacy/hpc/work/parcellations_overlap/T1``): same-session, same-surface,
   different execution -- must be identical.
3. **HCP vs DK on the same surfaces**: vertex-weighted whole-cortex means.

Usage (from the repo root, with the release dir resolvable)::

    python tools/validate_local_vs_release.py \
        --release-root /rds/project/rds-CeXlNYOYMxw/derivatives/tabulated \
        --release-label 6.0
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abcd import paths  # noqa: E402
from abcd import hcp_stats  # noqa: E402

KEY = ["participant_id", "session_id"]


def _md(df: pd.DataFrame, floatfmt: str = ".5f", index: bool = False) -> str:
    """Minimal Markdown table (avoids the optional ``tabulate`` dependency)."""
    if index:
        df = df.reset_index()
    def fmt(v):
        if isinstance(v, float):
            return format(v, floatfmt)
        return str(v)
    head = "| " + " | ".join(map(str, df.columns)) + " |"
    sep = "|" + "|".join("---" for _ in df.columns) + "|"
    rows = ["| " + " | ".join(fmt(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join([head, sep, *rows])

def _read(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t", low_memory=False)


def compare_dk(local: pd.DataFrame, release: pd.DataFrame, out: Path, label: str) -> str:
    cols = [c for c in local.columns if c.startswith("mr_y_smri__thk__dsk") and c in release.columns]
    region_cols = [c for c in cols if c.count("__") == 4]
    lk = set(map(tuple, local[KEY].values)); rk = set(map(tuple, release[KEY].values))
    only_local = sorted(lk - rk); only_release = sorted(rk - lk)
    m = local.merge(release, on=KEY, suffixes=("_loc", "_rel"))
    # release prints 3 decimals; compare after rounding local to the same
    rows = []
    for c in region_cols:
        a = m[f"{c}_loc"].round(3); b = m[f"{c}_rel"]
        d = (a - b).abs()
        rows.append({"column": c, "n": int(d.notna().sum()),
                     "frac_equal": float((d <= 1e-9).mean()),
                     "frac_within_0.001": float((d <= 0.0011).mean()),
                     "max_abs_diff": float(d.max()), "n_diff_gt_0.01": int((d > 0.01).sum())})
    per = pd.DataFrame(rows)
    per.to_csv(out / f"dk_local_vs_release_{label}_by_region.tsv", sep="\t", index=False)
    glob = []
    for c in [c for c in cols if c.count("__") < 4]:
        a = m[f"{c}_loc"]; b = m[f"{c}_rel"]
        glob.append({"column": c, "corr": float(a.corr(b)), "mean_abs_diff": float((a - b).abs().mean()),
                     "max_abs_diff": float((a - b).abs().max())})
    glob = pd.DataFrame(glob)
    pd.DataFrame(only_local, columns=KEY).to_csv(out / f"dk_sessions_only_local_vs_{label}.tsv", sep="\t", index=False)
    pd.DataFrame(only_release, columns=KEY).to_csv(out / f"dk_sessions_only_release_{label}.tsv", sep="\t", index=False)

    ol = pd.DataFrame(only_local, columns=KEY).session_id.value_counts().to_dict() if only_local else {}
    orl = pd.DataFrame(only_release, columns=KEY).session_id.value_counts().to_dict() if only_release else {}
    bad = per[per.frac_equal < 1]
    txt = [
        f"### 1. DK local (FreeSurfer `aparc.stats`) vs release {label} tables",
        "",
        f"- sessions local / release / shared: {len(lk):,} / {len(rk):,} / {len(m):,}",
        f"- local only (no release row): {len(only_local):,} by wave {ol}",
        f"- release only (no FreeSurfer directory on rds): {len(only_release):,} by wave {orl}",
        f"- region columns compared: {len(region_cols)}; values equal at 3 dp in "
        f"{per.frac_equal.mean():.5f} of cells overall; columns not 100 % equal: {len(bad)}",
        f"- largest absolute difference in any region cell: {per.max_abs_diff.max():.4f} mm",
        "",
        "Hemisphere / whole-cortex means (local FreeSurfer `Cortex MeanThickness` vs release):",
        "",
        _md(glob),
        "",
    ]
    if len(bad):
        txt += ["Columns with any unequal cell:", "", _md(bad, ".4f"), ""]
    return "\n".join(txt)


def compare_hcp_overlap(primary_root: Path, overlap_root: Path, fs_root: Path, out: Path) -> str:
    sessions = [s for s in hcp_stats.discover_sessions(overlap_root)]
    rows = []
    for s in sessions:
        a, qa = hcp_stats.extract_session(s, overlap_root, fs_root)
        b, qb = hcp_stats.extract_session(s, primary_root, fs_root)
        if a is None or b is None:
            rows.append({"session": s.key, "status": f"{qa['status']} / {qb['status']}"}); continue
        mm = a.merge(b, on=["label"], suffixes=("_new", "_old"))
        for c in ("nverts", "thickness_mm", "area_mm2", "gmv_mm3"):
            rows.append({"session": s.key, "measure": c,
                         "max_abs_diff": float((mm[f"{c}_new"] - mm[f"{c}_old"]).abs().max())})
    df = pd.DataFrame(rows)
    df.to_csv(out / "hcp_overlap_rerun_vs_july.tsv", sep="\t", index=False)
    ok = df.dropna(subset=["measure"]) if "measure" in df else df
    txt = ["### 2. HCP backfill re-run vs July output, same sessions", "",
           f"- sessions re-run: {len(sessions)}; comparable: {ok.session.nunique() if len(ok) else 0}"]
    if len(ok):
        summ = ok.groupby("measure").max_abs_diff.agg(["max", lambda x: (x == 0).mean()])
        summ.columns = ["max_abs_diff", "frac_sessions_identical"]
        txt += ["", _md(summ, ".4f", index=True), ""]
    if "status" in df and df.status.notna().any():
        txt += ["Sessions not comparable:", _md(df[df.status.notna()]), ""]
    return "\n".join(txt)


def compare_hcp_dk(hcp: pd.DataFrame, dk: pd.DataFrame) -> str:
    m = hcp[KEY + ["mr_y_smri__thk__hcp_mean"]].merge(dk[KEY + ["mr_y_smri__thk__dsk_mean"]], on=KEY)
    d = m["mr_y_smri__thk__hcp_mean"] - m["mr_y_smri__thk__dsk_mean"]
    return "\n".join([
        "### 3. HCP vs DK whole-cortex mean on the same surfaces", "",
        f"- shared sessions: {len(m):,}; r = {m.iloc[:, 2].corr(m.iloc[:, 3]):.5f}; "
        f"mean diff {d.mean():+.4f} mm; max |diff| {d.abs().max():.4f} mm",
        "- the two parcellations cover slightly different vertex sets (HCP's medial-wall "
        "`???` vs FreeSurfer's cortex label), so exact equality is not expected.", "",
    ])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release-root", type=Path, required=True,
                    help="directory holding the release tables (y/..., g/...)")
    ap.add_argument("--release-label", default="6.0")
    ap.add_argument("--processed", type=Path, default=None,
                    help="default: <release 7.0 dir>/processed")
    ap.add_argument("--overlap-root", type=Path, default=ROOT / "legacy/hpc/work/parcellations_overlap/T1")
    ap.add_argument("--out", type=Path, default=ROOT / "docs/hcp_census")
    a = ap.parse_args(argv)
    processed = a.processed or paths.release_dir("7.0") / "processed"
    a.out.mkdir(parents=True, exist_ok=True)

    parts = [f"# Local FreeSurfer-derived tables vs release {a.release_label}",
             "", f"Generated by `tools/validate_local_vs_release.py` on {pd.Timestamp.now():%Y-%m-%d}.", ""]
    dk_local = _read(processed / "dsk_local" / "mr_y_smri__thk__dsk.tsv")
    dk_rel = _read(paths.find_table(a.release_root, "mr_y_smri__thk__dsk"))
    parts.append(compare_dk(dk_local, dk_rel, a.out, a.release_label))
    if a.overlap_root.exists():
        parts.append(compare_hcp_overlap(hcp_stats.DEFAULT_PARC_ROOT, a.overlap_root,
                                         hcp_stats.DEFAULT_FS_ROOT, a.out))
    hcp = _read(processed / "hcp" / "mr_y_smri__thk__hcp.tsv")
    parts.append(compare_hcp_dk(hcp, dk_local))
    report = a.out / f"local_vs_release_{a.release_label}.md"
    report.write_text("\n".join(parts))
    print("\n".join(parts))
    print(f"\n-> {report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
