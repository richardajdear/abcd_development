"""
HCP-MMP1.0 (Glasser) cortical statistics from FreeSurfer surfaces.

ABCD does not tabulate the HCP-MMP parcellation, so it has to be produced from
the release FreeSurfer reconstructions.  On CSD3 those live at::

    /rds/project/rds-CeXlNYOYMxw/derivatives/freesurfer/<sub>/<ses>/     (DAIRC FreeSurfer 7.1.1)
    /rds/project/rds-CeXlNYOYMxw/derivatives/parcellations/T1/<sub>/<ses>/HCP.fsaverage.aparc/
        {lh,rh}.HCP.fsaverage.aparc.log                                   (mris_anatomical_stats)

The second tree is the output of the group's parcellation pipeline
(``Code/parcellation_T1/parcellate.sh``, R. Romero-Garcia): the fsaverage
HCP-MMP1.0 annotation is carried to each subject's sphere with
``mri_surf2surf`` and summarised with ``mris_anatomical_stats -a <annot> -b``.
Each ``.log`` is a fixed-width table with one row per parcel (180 ``*_ROI``
rows plus a ``???`` medial-wall row).

This module parses those tables into

1. a **long** parquet with every measure `mris_anatomical_stats` reports
   (vertex count, area, GM volume, thickness mean/SD, curvatures, folding);
2. **wide TSVs in the release-7.0 column convention** --
   ``mr_y_smri__thk__hcp__<region>__<hemi>_mean`` and so on -- so the
   :class:`abcd.io.Release70Adapter` can read them exactly as it reads the
   Desikan release tables; and
3. a **per-session QC table** with parcel counts, vertex counts and the
   surface-hole counts from ``aseg.stats`` (the Euler-number equivalent,
   ``euler = 2 - 2*holes`` per hemisphere).

Outputs land in ``<release_dir>/processed/hcp/``, next to the release's
tabulated tables, mirroring the ``processed/`` convention the 5.1 HCP files
already use.  They are derived data, not release tables, and stay gitignored
with the rest of ``abcd-7.0/``.

Run it on the cluster (30k sessions is not a login-node job)::

    sbatch hpc/hcp_extract.sbatch

or, for a handful of sessions::

    python -m abcd.hcp_stats --limit 5 --out-dir /tmp/hcp_test
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------
# Locations (CSD3).  Override on the command line.
# --------------------------------------------------------------------------

RDS_ROOT = Path("/rds/project/rds-CeXlNYOYMxw")
DEFAULT_PARC_ROOT = RDS_ROOT / "derivatives" / "parcellations" / "T1"
DEFAULT_FS_ROOT = RDS_ROOT / "derivatives" / "freesurfer"
ATLAS = "HCP.fsaverage.aparc"

#: Columns of the ``mris_anatomical_stats`` table, in order.
STATS_COLUMNS = (
    "nverts", "area_mm2", "gmv_mm3", "thickness_mm", "thickness_sd",
    "mean_curv", "gauss_curv", "folding_index", "intrinsic_curv",
)
N_PARCELS_PER_HEMI = 180
MEDIAL_WALL = "???"

#: metric -> (7.0 column infix, per-parcel aggregation suffix)
#: ``mean`` measures are averaged over vertices, ``sum`` measures are totals,
#: matching the release DK tables (``thk__dsk__*__lh_mean``,
#: ``area__dsk__*__lh_sum``, ``vol__dsk__*__lh_sum``).
WIDE_METRICS = {
    "thickness_mm": ("thk", "mean"),
    "area_mm2": ("area", "sum"),
    "gmv_mm3": ("vol", "sum"),
}

_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+"                # nverts area gmv
    r"([-\d.]+)\s+([-\d.]+)\s+"                    # thickness mean, sd
    r"([-\d.]+)\s+([-\d.]+)\s+"                    # mean curv, gauss curv
    r"([-\d.]+)\s+([-\d.]+)\s+"                    # folding index, ici
    r"(\S+)\s*$"                                   # structure name
)
_ROI_RE = re.compile(r"^[LR]_(.+)_ROI$")
_HOLES_RE = re.compile(r"^# Measure (lh|rh)?SurfaceHoles, .*?, (\d+), unitless")


# --------------------------------------------------------------------------
# Parsers
# --------------------------------------------------------------------------

def parse_anatomical_stats(path: Path) -> pd.DataFrame:
    """Parse one ``mris_anatomical_stats`` table.

    Returns one row per structure with :data:`STATS_COLUMNS` plus
    ``structure`` (the raw name, e.g. ``L_V1_ROI`` or ``???``).  Header and
    diagnostic lines are skipped; a line that looks numeric but does not have
    exactly nine values is an error rather than silently dropped.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            m = _ROW_RE.match(line)
            if m:
                vals = m.groups()
                rows.append({
                    "nverts": int(vals[0]),
                    "area_mm2": float(vals[1]),
                    "gmv_mm3": float(vals[2]),
                    "thickness_mm": float(vals[3]),
                    "thickness_sd": float(vals[4]),
                    "mean_curv": float(vals[5]),
                    "gauss_curv": float(vals[6]),
                    "folding_index": float(vals[7]),
                    "intrinsic_curv": float(vals[8]),
                    "structure": vals[9],
                })
            elif line.strip() and line.strip()[0].isdigit():
                raise ValueError(f"{path}: unparsed numeric line: {line.rstrip()!r}")
    if not rows:
        raise ValueError(f"{path}: no parcel rows found")
    return pd.DataFrame(rows)


def parse_aseg_holes(path: Path) -> dict[str, float]:
    """Surface-hole counts from ``aseg.stats`` (NaN if the file is missing)."""
    out = {"lh_holes": float("nan"), "rh_holes": float("nan"),
           "total_holes": float("nan")}
    if not path.exists():
        return out
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = _HOLES_RE.match(line)
            if m:
                key = {"lh": "lh_holes", "rh": "rh_holes", None: "total_holes"}[m.group(1)]
                out[key] = float(m.group(2))
    return out


def region_from_structure(structure: str) -> str | None:
    """``L_V1_ROI`` -> ``V1``; the medial wall and anything unexpected -> None."""
    m = _ROI_RE.match(structure)
    return m.group(1) if m else None


# --------------------------------------------------------------------------
# Per-session extraction
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Session:
    subject: str   # sub-003RTV85
    session: str   # ses-00A

    @property
    def key(self) -> str:
        return f"{self.subject}/{self.session}"


def discover_sessions(parc_root: Path) -> list[Session]:
    """Every ``<sub>/<ses>`` under the parcellation tree, sorted.

    Uses ``os.scandir`` rather than ``glob`` so a 30k-session tree is two
    directory listings deep, not a stat per file.
    """
    out = []
    with os.scandir(parc_root) as subs:
        for sub in subs:
            if not (sub.is_dir() and sub.name.startswith("sub-")):
                continue
            with os.scandir(sub.path) as sess:
                for s in sess:
                    if s.is_dir() and s.name.startswith("ses-"):
                        out.append(Session(sub.name, s.name))
    return sorted(out, key=lambda s: s.key)


def extract_session(sess: Session, parc_root: Path, fs_root: Path
                    ) -> tuple[pd.DataFrame | None, dict]:
    """Parse both hemispheres of one session.

    Returns ``(long, qc)``.  ``long`` is None when either hemisphere table is
    missing or malformed; the reason is recorded in ``qc['status']`` so the
    session is accounted for rather than vanishing.
    """
    qc: dict = {
        "participant_id": sess.subject, "session_id": sess.session,
        "status": "ok", "n_parcels_lh": 0, "n_parcels_rh": 0,
        "nverts_lh": 0, "nverts_rh": 0, "medial_wall_verts": 0,
    }
    qc.update(parse_aseg_holes(fs_root / sess.subject / sess.session / "stats" / "aseg.stats"))
    qc["recon_done"] = (fs_root / sess.subject / sess.session / "scripts" / "recon-all.done").exists()

    frames = []
    for hemi in ("lh", "rh"):
        p = parc_root / sess.subject / sess.session / ATLAS / f"{hemi}.{ATLAS}.log"
        if not p.exists() or p.stat().st_size == 0:
            qc["status"] = f"missing_{hemi}"
            return None, qc
        try:
            df = parse_anatomical_stats(p)
        except ValueError as exc:
            qc["status"] = f"parse_error_{hemi}: {exc}"
            return None, qc
        df["hemi"] = hemi
        df["region"] = df.structure.map(region_from_structure)
        wall = df.structure == MEDIAL_WALL
        qc["medial_wall_verts"] += int(df.loc[wall, "nverts"].sum())
        df = df[~wall]
        if df.region.isna().any():
            bad = df.loc[df.region.isna(), "structure"].tolist()[:3]
            qc["status"] = f"unexpected_structure_{hemi}: {bad}"
            return None, qc
        qc[f"n_parcels_{hemi}"] = int(len(df))
        qc[f"nverts_{hemi}"] = int(df.nverts.sum())
        frames.append(df)

    if qc["n_parcels_lh"] != N_PARCELS_PER_HEMI or qc["n_parcels_rh"] != N_PARCELS_PER_HEMI:
        qc["status"] = f"parcel_count_{qc['n_parcels_lh']}_{qc['n_parcels_rh']}"
        return None, qc

    long = pd.concat(frames, ignore_index=True)
    long.insert(0, "session_id", sess.session)
    long.insert(0, "participant_id", sess.subject)
    long["label"] = long.hemi + "_" + long.region
    return long[["participant_id", "session_id", "hemi", "region", "label",
                 *STATS_COLUMNS]], qc


def _worker(args):
    sess, parc_root, fs_root = args
    return extract_session(sess, Path(parc_root), Path(fs_root))


# --------------------------------------------------------------------------
# Wide tables in the release-7.0 convention
# --------------------------------------------------------------------------

def wide_table(long: pd.DataFrame, measure: str) -> pd.DataFrame:
    """Pivot one measure to ``participant_id, session_id, mr_y_smri__<infix>__hcp__<region>__<hemi>_<agg>``.

    Whole-cortex and per-hemisphere summaries follow the DK tables:
    ``..._mean`` columns are **vertex-weighted** means over the 360 (or 180)
    parcels, ``..._sum`` columns are totals.  The DK release tables are
    computed the same way (a mean over all cortical vertices), which is why the
    weighted rather than the parcel-average form is used here.
    """
    infix, agg = WIDE_METRICS[measure]
    prefix = f"mr_y_smri__{infix}__hcp"

    piv = long.pivot_table(index=["participant_id", "session_id"],
                           columns=["region", "hemi"], values=measure,
                           aggfunc="first")
    piv.columns = [f"{prefix}__{r}__{h}_{agg}" for r, h in piv.columns]

    g = long.assign(w=long.nverts * long[measure] if agg == "mean" else long[measure])
    grp = g.groupby(["participant_id", "session_id"], sort=False)
    per_hemi = g.groupby(["participant_id", "session_id", "hemi"], sort=False)
    if agg == "mean":
        tot = grp.w.sum() / grp.nverts.sum()
        hemi = per_hemi.w.sum() / per_hemi.nverts.sum()
    else:
        tot = grp.w.sum()
        hemi = per_hemi.w.sum()
    hemi = hemi.unstack("hemi")
    for h in ("lh", "rh"):
        piv[f"{prefix}__{h}_{agg}"] = hemi[h]
    piv[f"{prefix}_{agg}"] = tot
    return piv.reset_index()


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def run(parc_root: Path, fs_root: Path, out_dir: Path, workers: int = 8,
        limit: int | None = None, log=print) -> dict:
    t0 = time.time()
    sessions = discover_sessions(parc_root)
    log(f"{len(sessions)} sessions under {parc_root}")
    if limit:
        sessions = sessions[:limit]
        log(f"limiting to first {limit}")

    longs, qcs = [], []
    tasks = [(s, str(parc_root), str(fs_root)) for s in sessions]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, (long, qc) in enumerate(ex.map(_worker, tasks, chunksize=64), 1):
            qcs.append(qc)
            if long is not None:
                longs.append(long)
            if i % 2000 == 0:
                log(f"  {i}/{len(sessions)}  {time.time() - t0:.0f}s")

    qc = pd.DataFrame(qcs)
    long = pd.concat(longs, ignore_index=True)
    ok = (qc.status == "ok").sum()
    log(f"parsed {ok}/{len(qc)} sessions; status counts:\n{qc.status.str.split(':').str[0].value_counts().to_string()}")

    out_dir.mkdir(parents=True, exist_ok=True)
    long.to_parquet(out_dir / "hcp_anatomical_stats_long.parquet", index=False)
    qc.to_csv(out_dir / "hcp_session_qc.tsv", sep="\t", index=False)
    written = {"long": len(long), "qc": len(qc)}
    for measure, (infix, _agg) in WIDE_METRICS.items():
        w = wide_table(long, measure)
        name = f"mr_y_smri__{infix}__hcp.tsv"
        w.to_csv(out_dir / name, sep="\t", index=False, float_format="%.4f")
        written[name] = w.shape
        log(f"wrote {name}: {w.shape[0]} rows x {w.shape[1]} cols")

    # provenance next to the data
    stamp = out_dir / "PROVENANCE.txt"
    stamp.write_text(
        "HCP-MMP1.0 cortical statistics parsed from mris_anatomical_stats tables\n"
        f"generated {time.strftime('%Y-%m-%d %H:%M')} by abcd.hcp_stats\n"
        f"parcellation tree : {parc_root}\n"
        f"freesurfer tree   : {fs_root}\n"
        f"sessions found    : {len(qc)}\n"
        f"sessions parsed   : {ok}\n"
        "column convention : mr_y_smri__{thk,area,vol}__hcp__<region>__<hemi>_{mean,sum};"
        " whole-cortex/hemisphere means are vertex-weighted\n"
        "medial wall (???) excluded; its vertex count is in hcp_session_qc.tsv\n"
    )
    log(f"done in {time.time() - t0:.0f}s -> {out_dir}")
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--parc-root", type=Path, default=DEFAULT_PARC_ROOT)
    ap.add_argument("--fs-root", type=Path, default=DEFAULT_FS_ROOT)
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="default: <release 7.0 dir>/processed/hcp")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None,
                    help="parse only the first N sessions (testing)")
    a = ap.parse_args(argv)

    out_dir = a.out_dir
    if out_dir is None:
        from . import paths
        out_dir = paths.release_dir("7.0") / "processed" / "hcp"

    def log(msg):
        print(msg, flush=True)

    run(a.parc_root, a.fs_root, out_dir, workers=a.workers, limit=a.limit, log=log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
