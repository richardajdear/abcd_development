"""
Desikan-Killiany cortical statistics parsed directly from FreeSurfer.

Why this exists when ABCD already tabulates DK: (1) to check that the release
tables really are these surfaces, value for value, and (2) to have DK
thickness for the sessions newer than the tabulated release on hand (the
FreeSurfer derivatives on rds run to July 2025; the 6.0 tables stop in early
2024).  The release table remains the canonical DK source for the pipeline;
this module writes a *parallel* table with the **same column names** under
``<release_dir>/processed/dsk_local/`` so the two can be joined directly.

No FreeSurfer command is run.  ``recon-all`` already wrote
``stats/{lh,rh}.aparc.stats`` for every session; this is a parse.  The
``aparc.stats`` layout differs from the ``mris_anatomical_stats -b`` table
that :mod:`abcd.hcp_stats` reads -- the structure name comes *first*, and the
header carries hemisphere totals (``Cortex, NumVert`` / ``MeanThickness`` /
``WhiteSurfArea``) which are what the release ``__lh_mean`` columns should
equal.

    sbatch tools/dk_extract.sbatch
    python -m abcd.dk_stats --limit 5 --out-dir /tmp/dk_test
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from .hcp_stats import DEFAULT_FS_ROOT, STATS_COLUMNS, parse_aseg_holes

#: aparc.stats ColHeaders -> our names (same measures as the HCP long table)
APARC_COLUMNS = {
    "StructName": "region", "NumVert": "nverts", "SurfArea": "area_mm2",
    "GrayVol": "gmv_mm3", "ThickAvg": "thickness_mm", "ThickStd": "thickness_sd",
    "MeanCurv": "mean_curv", "GausCurv": "gauss_curv", "FoldInd": "folding_index",
    "CurvInd": "intrinsic_curv",
}
N_REGIONS_PER_HEMI = 34
_MEASURE_RE = re.compile(r"^# Measure Cortex, (\w+), .*?, ([-\d.]+), (\S+)")

#: metric -> (7.0 column infix, per-region aggregation suffix), as in the release
WIDE_METRICS = {
    "thickness_mm": ("thk", "mean"),
    "area_mm2": ("area", "sum"),
    "gmv_mm3": ("vol", "sum"),
}


def _region_to_code() -> dict[str, str]:
    """Canonical DK region name -> 7.0 column abbreviation (``bankssts`` -> ``bstmps``).

    Built from the two mappings the adapter already validates against the
    release: ``region_labels.csv`` (name -> 5.1 token) and
    ``Release70Adapter.REGION_CODES`` (7.0 code -> 5.1 token).  Asserted
    complete for all 34 regions.
    """
    from .io import Release70Adapter, region_labels
    lab = region_labels("dsk")
    lab = lab[~lab.is_global]
    tok = lab.stem_a.str.replace(r"(lh|rh)$", "", regex=True)
    name_to_tok = dict(zip(lab.region, tok))
    tok_to_code = {v: k for k, v in Release70Adapter.REGION_CODES.items()}
    out = {}
    for name, t in name_to_tok.items():
        if t not in tok_to_code:
            raise KeyError(f"DK region {name!r} (token {t!r}) has no 7.0 column code")
        out[name] = tok_to_code[t]
    if len(out) != N_REGIONS_PER_HEMI:
        raise ValueError(f"expected 34 DK regions, mapped {len(out)}")
    return out


def parse_aparc_stats(path: Path) -> tuple[pd.DataFrame, dict[str, float]]:
    """One ``?h.aparc.stats`` -> (per-region table, hemisphere measures).

    The header ``# Measure Cortex, <name>, ..., <value>, <unit>`` lines give
    ``NumVert``, ``WhiteSurfArea`` and ``MeanThickness`` for the hemisphere.
    """
    measures: dict[str, float] = {}
    header: list[str] | None = None
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                m = _MEASURE_RE.match(line)
                if m:
                    measures[m.group(1)] = float(m.group(2))
                elif line.startswith("# ColHeaders"):
                    header = line.split()[2:]
                continue
            if not line.strip():
                continue
            if header is None:
                raise ValueError(f"{path}: data before ColHeaders")
            parts = line.split()
            if len(parts) != len(header):
                raise ValueError(f"{path}: bad row {line.rstrip()!r}")
            rows.append(dict(zip(header, parts)))
    if not rows:
        raise ValueError(f"{path}: no region rows")
    df = pd.DataFrame(rows).rename(columns=APARC_COLUMNS)
    for c in STATS_COLUMNS:
        df[c] = pd.to_numeric(df[c])
    df["nverts"] = df.nverts.astype(int)
    return df[["region", *STATS_COLUMNS]], measures


def extract_session(sub: str, ses: str, fs_root: Path
                    ) -> tuple[pd.DataFrame | None, dict]:
    qc: dict = {"participant_id": sub, "session_id": ses, "status": "ok"}
    sdir = fs_root / sub / ses
    qc.update(parse_aseg_holes(sdir / "stats" / "aseg.stats"))
    qc["recon_done"] = (sdir / "scripts" / "recon-all.done").exists()
    frames = []
    for hemi in ("lh", "rh"):
        p = sdir / "stats" / f"{hemi}.aparc.stats"
        if not p.exists() or p.stat().st_size == 0:
            qc["status"] = f"missing_{hemi}"
            return None, qc
        try:
            df, meas = parse_aparc_stats(p)
        except ValueError as exc:
            qc["status"] = f"parse_error_{hemi}: {exc}"
            return None, qc
        qc[f"n_regions_{hemi}"] = len(df)
        qc[f"cortex_nverts_{hemi}"] = meas.get("NumVert")
        qc[f"cortex_mean_thickness_{hemi}"] = meas.get("MeanThickness")
        qc[f"cortex_white_area_{hemi}"] = meas.get("WhiteSurfArea")
        df["hemi"] = hemi
        frames.append(df)
    if any(qc[f"n_regions_{h}"] != N_REGIONS_PER_HEMI for h in ("lh", "rh")):
        qc["status"] = f"region_count_{qc['n_regions_lh']}_{qc['n_regions_rh']}"
        return None, qc
    long = pd.concat(frames, ignore_index=True)
    long.insert(0, "session_id", ses)
    long.insert(0, "participant_id", sub)
    long["label"] = long.hemi + "_" + long.region
    return long[["participant_id", "session_id", "hemi", "region", "label",
                 *STATS_COLUMNS]], qc


def _worker(args):
    sub, ses, fs_root = args
    return extract_session(sub, ses, Path(fs_root))


def discover_sessions(fs_root: Path) -> list[tuple[str, str]]:
    out = []
    with os.scandir(fs_root) as subs:
        for sub in subs:
            if not (sub.is_dir(follow_symlinks=False) and sub.name.startswith("sub-")):
                continue
            with os.scandir(sub.path) as sess:
                for s in sess:
                    if s.is_dir(follow_symlinks=False) and s.name.startswith("ses-"):
                        out.append((sub.name, s.name))
    return sorted(out)


def wide_table(long: pd.DataFrame, measure: str, qc: pd.DataFrame,
               codes: dict[str, str]) -> pd.DataFrame:
    """Release-shaped table: ``mr_y_smri__thk__dsk__bstmps__lh_mean`` etc.

    Hemisphere and whole-cortex summaries follow the release convention,
    established empirically on ``sub-003RTV85/ses-00A``: the release
    ``__lh_mean`` (2.71228) is the **surface-area-weighted mean over the 34
    DK regions** -- not FreeSurfer's ``Cortex MeanThickness`` (2.71093,
    vertex-weighted over the cortex label) and not the unweighted region mean
    (2.75106).  Area and volume summaries are sums over the regions.
    FreeSurfer's cortex mean is kept in the QC table for reference.
    """
    infix, agg = WIDE_METRICS[measure]
    prefix = f"mr_y_smri__{infix}__dsk"
    piv = long.pivot_table(index=["participant_id", "session_id"],
                           columns=["region", "hemi"], values=measure, aggfunc="first")
    piv.columns = [f"{prefix}__{codes[r]}__{h}_{agg}" for r, h in piv.columns]

    key = ["participant_id", "session_id"]
    if measure == "thickness_mm":
        g = long.assign(w=long.area_mm2 * long.thickness_mm)
        per_hemi = g.groupby(key + ["hemi"], sort=False)
        hemi = (per_hemi.w.sum() / per_hemi.area_mm2.sum()).unstack("hemi")
        for h in ("lh", "rh"):
            piv[f"{prefix}__{h}_mean"] = hemi[h]
        tot = g.groupby(key, sort=False)
        piv[f"{prefix}_mean"] = tot.w.sum() / tot.area_mm2.sum()
    else:
        per_hemi = long.groupby(key + ["hemi"], sort=False)[measure].sum().unstack("hemi")
        for h in ("lh", "rh"):
            piv[f"{prefix}__{h}_sum"] = per_hemi[h]
        piv[f"{prefix}_sum"] = per_hemi.sum(axis=1)
    return piv.reset_index()


def run(fs_root: Path, out_dir: Path, workers: int = 8, limit: int | None = None,
        log=print) -> None:
    t0 = time.time()
    codes = _region_to_code()
    sessions = discover_sessions(fs_root)
    log(f"{len(sessions)} FreeSurfer sessions under {fs_root}")
    if limit:
        sessions = sessions[:limit]
    longs, qcs = [], []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, (long, qc) in enumerate(
                ex.map(_worker, [(s, e, str(fs_root)) for s, e in sessions], chunksize=64), 1):
            qcs.append(qc)
            if long is not None:
                longs.append(long)
            if i % 4000 == 0:
                log(f"  {i}/{len(sessions)}  {time.time() - t0:.0f}s")
    qc = pd.DataFrame(qcs)
    long = pd.concat(longs, ignore_index=True)
    ok = qc[qc.status == "ok"]
    log(f"parsed {len(ok)}/{len(qc)} sessions; status counts:\n"
        f"{qc.status.str.split(':').str[0].value_counts().to_string()}")

    out_dir.mkdir(parents=True, exist_ok=True)
    long.to_parquet(out_dir / "dsk_aparc_stats_long.parquet", index=False)
    qc.to_csv(out_dir / "dsk_session_qc.tsv", sep="\t", index=False)
    for measure, (infix, _agg) in WIDE_METRICS.items():
        w = wide_table(long, measure, ok, codes)
        name = f"mr_y_smri__{infix}__dsk.tsv"
        w.to_csv(out_dir / name, sep="\t", index=False, float_format="%.5f")
        log(f"wrote {name}: {w.shape[0]} rows x {w.shape[1]} cols")
    (out_dir / "PROVENANCE.txt").write_text(
        "Desikan-Killiany statistics parsed from FreeSurfer stats/?h.aparc.stats\n"
        f"generated {time.strftime('%Y-%m-%d %H:%M')} by abcd.dk_stats\n"
        f"freesurfer tree : {fs_root}\nsessions found  : {len(qc)}\nsessions parsed : {len(ok)}\n"
        "columns named exactly as the release mr_y_smri__{thk,area,vol}__dsk tables; "
        "hemisphere/whole-cortex thickness means are surface-area-weighted over regions (release convention).\n"
        "NOT the canonical DK source for the pipeline -- a comparison object; see README.\n")
    log(f"done in {time.time() - t0:.0f}s -> {out_dir}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--fs-root", type=Path, default=DEFAULT_FS_ROOT)
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="default: <release 7.0 dir>/processed/dsk_local")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args(argv)
    out_dir = a.out_dir
    if out_dir is None:
        from . import paths
        out_dir = paths.release_dir("7.0") / "processed" / "dsk_local"
    run(a.fs_root, out_dir, workers=a.workers, limit=a.limit,
        log=lambda m: print(m, flush=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
