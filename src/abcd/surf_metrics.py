"""
Per-parcel cortical intensity and contrast measures sampled on the FreeSurfer surfaces (D1).

The thickness table ``processed/hcp/mr_y_smri__thk__hcp.tsv`` (:mod:`abcd.hcp_stats`)
comes from the DAIRC FreeSurfer 7.1.1 surfaces and rr480's subject-space HCP-MMP1.0
annotation.  This module samples further measures on **the same surfaces and the same
annotation**, so they are parcel-for-parcel comparable with thickness:

=============  =============================================  ==============================
measure        what                                           source
=============  =============================================  ==============================
``gwc``        grey/white contrast, FreeSurfer ``pctsurfcon``  ``surf/?h.w-g.pct.mgh``
``t1__gm50``   T1w intensity at GM mid-depth                   ``mri/nu.mgz``, white + 0.5 (pial - white)
``t1__wm1``    T1w intensity 1 mm below the white surface      ``mri/nu.mgz``, white - 1 mm x outward normal
=============  =============================================  ==============================

Volume measures are sampled with trilinear interpolation
(``scipy.ndimage.map_coordinates``, order 1) at tkr-RAS vertex positions mapped to
voxels with the inverse of the volume's ``vox2ras-tkr``.  The white-surface normal is the
normalised sum of the area-weighted normals of the faces around each vertex, oriented
towards the pial surface.  Samples outside the volume are NaN.

Per parcel the summary is the **vertex-area-weighted mean** (weights ``surf/?h.area``)
over the finite samples of the parcel's vertices; the unweighted vertex mean is kept too
(FreeSurfer's ``mri_segstats`` reports that one, so the GWC gate uses it).  Hemisphere and
whole-cortex means are area-weighted over all vertices outside the medial wall
(``???`` for HCP; ``unknown``/``corpuscallosum`` for DK; label -1 for both).

Annotation per session, resolved exactly as :func:`abcd.hcp_stats.discover_sessions`
resolves the ``.log`` tables: a session in the backfill tree reads
``<backfill>/<sub>/<ses>/HCP.fsaverage.aparc/?h.HCP.fsaverage.aparc.annot``; every other
session reads the annotation rr480's pipeline wrote next to the surfaces,
``<fs>/<sub>/<ses>/label/?h.HCP.fsaverage.aparc.annot``.  DK is ``label/?h.aparc.annot``.

Work is split in chunks (one array task each); each chunk writes a long parquet and a QC
table under ``<out>/surf_metrics_parts/``.  ``tools/surf_metrics_collect.py`` assembles the
release-convention wide TSVs.

    sbatch --array=0-169 tools/surf_metrics.sbatch           # all sessions, 200 per task
    python -m abcd.surf_metrics --sessions pilot.tsv --part-name pilot --out-dir /tmp/sm
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.freesurfer import read_annot, read_geometry, read_morph_data
from scipy.ndimage import map_coordinates

from .hcp_stats import (ATLAS, DEFAULT_FS_ROOT, DEFAULT_PARC_ROOT, Session,
                        discover_sessions, region_from_structure)

DEFAULT_BACKFILL_ROOT = Path(
    "/rds/user/rajd2/hpc-work/abcd_development/legacy/hpc/work/parcellations_backfill/T1")
HEMIS = ("lh", "rh")
DK_EXCLUDE = {"unknown", "corpuscallosum"}
HCP_MEDIAL_WALL = "???"


@dataclass(frozen=True)
class Measure:
    name: str               # column infix after mr_y_smri__, e.g. "t1__gm50"
    kind: str               # "surface" (per-vertex file) or "volume" (sampled)
    source: str             # path under the FreeSurfer session; {hemi} filled in
    depth: str | None = None  # gm50 / wm1 / wm2 for volume measures


MEASURES = {
    m.name: m for m in (
        Measure("gwc", "surface", "surf/{hemi}.w-g.pct.mgh"),
        Measure("t1__gm50", "volume", "mri/nu.mgz", "gm50"),
        Measure("t1__wm1", "volume", "mri/nu.mgz", "wm1"),
    )
}
DEFAULT_MEASURES = tuple(MEASURES)


# --------------------------------------------------------------------------
# Geometry and sampling
# --------------------------------------------------------------------------

def vertex_normals(coords: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Unit vertex normals: the sum of the adjacent faces' area-weighted normals.

    ``cross(v1 - v0, v2 - v0)`` has length twice the face area, so summing it is
    area-weighting.  The sign follows FreeSurfer's face winding; callers orient it with
    :func:`orient_outward`.
    """
    v0, v1, v2 = (coords[faces[:, k]] for k in range(3))
    fn = np.cross(v1 - v0, v2 - v0)
    n = len(coords)
    vn = np.empty_like(coords)
    for d in range(3):
        vn[:, d] = sum(np.bincount(faces[:, k], weights=fn[:, d], minlength=n) for k in range(3))
    norm = np.linalg.norm(vn, axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return vn / norm


def orient_outward(normals: np.ndarray, white: np.ndarray, pial: np.ndarray
                   ) -> tuple[np.ndarray, float]:
    """Flip ``normals`` so they point from white towards pial for most vertices.

    Returns the oriented normals and the fraction of vertices (with thickness > 0.1 mm)
    whose normal agrees with the white->pial direction, a sanity figure for the QC table.
    """
    d = pial - white
    ok = np.linalg.norm(d, axis=1) > 0.1
    dots = np.einsum("ij,ij->i", normals[ok], d[ok])
    if np.nanmedian(dots) < 0:
        normals = -normals
        dots = -dots
    return normals, float(np.mean(dots > 0))


def depth_coords(depth: str, white: np.ndarray, pial: np.ndarray, normals: np.ndarray
                 ) -> np.ndarray:
    if depth == "gm50":
        return white + 0.5 * (pial - white)
    if depth.startswith("wm"):
        return white - float(depth[2:]) * normals
    raise ValueError(f"unknown depth {depth!r}")


def sample_volume(data: np.ndarray, tkr2vox: np.ndarray, xyz: np.ndarray) -> np.ndarray:
    """Trilinear samples of ``data`` at tkr-RAS points ``xyz`` (n x 3); NaN outside."""
    ijk = xyz @ tkr2vox[:3, :3].T + tkr2vox[:3, 3]
    return map_coordinates(data, ijk.T, order=1, mode="constant", cval=np.nan)


def load_volume(path: Path) -> tuple[np.ndarray, np.ndarray]:
    img = nib.load(str(path))
    data = np.asarray(img.dataobj, dtype=np.float32)
    return data, np.linalg.inv(img.header.get_vox2ras_tkr())


def sample_hemi(sdir: Path, hemi: str, measures: tuple[str, ...],
                volumes: dict[str, tuple[np.ndarray, np.ndarray]]) -> tuple[dict, dict]:
    """Per-vertex values of every measure for one hemisphere.

    Returns ``(values, info)``: ``values[measure]`` is an n-vertex array (missing input ->
    all NaN), ``info`` holds the vertex area and diagnostics.
    """
    white, faces = read_geometry(str(sdir / "surf" / f"{hemi}.white"))
    pial, _ = read_geometry(str(sdir / "surf" / f"{hemi}.pial"))
    area = read_morph_data(str(sdir / "surf" / f"{hemi}.area"))
    if not (len(white) == len(pial) == len(area)):
        raise ValueError(f"{hemi}: white/pial/area vertex counts differ")
    normals, agree = orient_outward(vertex_normals(white, faces), white, pial)
    out, info = {}, {"area": area, "normal_agree": agree, "n": len(white)}
    for name in measures:
        m = MEASURES[name]
        if m.kind == "surface":
            p = sdir / m.source.format(hemi=hemi)
            out[name] = (np.asarray(nib.load(str(p)).dataobj, dtype=np.float64).reshape(-1)
                         if p.exists() else None)
            if out[name] is not None and len(out[name]) != len(white):
                raise ValueError(f"{p.name}: {len(out[name])} values for {len(white)} vertices")
        else:
            vol = volumes.get(m.source)
            out[name] = (None if vol is None
                         else sample_volume(vol[0], vol[1], depth_coords(m.depth, white, pial, normals)))
    return out, info


# --------------------------------------------------------------------------
# Parcel summaries
# --------------------------------------------------------------------------

def annot_regions(path: Path, atlas: str, dk_codes: dict[str, str]
                  ) -> tuple[np.ndarray, list[str | None]]:
    """Vertex label indices and the region token per index (None = medial wall/excluded).

    HCP tokens are as in :mod:`abcd.hcp_stats` (``L_V1_ROI`` -> ``V1``); DK tokens are the
    release 7.0 codes (``bankssts`` -> ``bstmps``).
    """
    labels, _ctab, names = read_annot(str(path))
    names = [n.decode() if isinstance(n, bytes) else n for n in names]
    if atlas == "hcp":
        tokens = [None if n == HCP_MEDIAL_WALL else region_from_structure(n) for n in names]
        bad = [n for n, t in zip(names, tokens) if t is None and n != HCP_MEDIAL_WALL]
    else:
        tokens = [None if n in DK_EXCLUDE else dk_codes.get(n) for n in names]
        bad = [n for n, t in zip(names, tokens) if t is None and n not in DK_EXCLUDE]
    if bad:
        raise ValueError(f"{path}: unexpected structures {bad[:3]}")
    return labels, tokens


def summarise(values: np.ndarray, area: np.ndarray, labels: np.ndarray,
              tokens: list[str | None]) -> tuple[list[dict], dict]:
    """Parcel rows plus the hemisphere summary over all non-medial vertices."""
    keep_tok = np.array([t is not None for t in tokens] + [False])   # index -1 -> last
    in_cortex = keep_tok[labels]
    finite = np.isfinite(values) & in_cortex
    lab = np.where(finite, labels, len(tokens))
    k = len(tokens) + 1
    n_all = np.bincount(np.where(in_cortex, labels, len(tokens)), minlength=k)
    n_fin = np.bincount(lab, minlength=k)
    v = np.where(finite, values, 0.0)
    a = np.where(finite, area, 0.0)
    sw = np.bincount(lab, weights=v * a, minlength=k)
    sa = np.bincount(lab, weights=a, minlength=k)
    su = np.bincount(lab, weights=v, minlength=k)
    rows = []
    with np.errstate(invalid="ignore", divide="ignore"):
        for i, tok in enumerate(tokens):
            if tok is None or n_all[i] == 0:
                continue
            rows.append({"region": tok, "nverts": int(n_all[i]), "n_finite": int(n_fin[i]),
                         "area": float(sa[i]), "mean_w": sw[i] / sa[i] if sa[i] else np.nan,
                         "mean_u": su[i] / n_fin[i] if n_fin[i] else np.nan})
    hemi = {"nverts": int(in_cortex.sum()), "n_finite": int(finite.sum()),
            "area": float(a.sum()), "wsum": float((v * a).sum()), "usum": float(v.sum())}
    return rows, hemi


# --------------------------------------------------------------------------
# Per session
# --------------------------------------------------------------------------

def annot_path(sess: Session, hemi: str, fs_root: Path, backfill_root: Path) -> Path:
    if sess.source == "backfill":
        return backfill_root / sess.subject / sess.session / ATLAS / f"{hemi}.{ATLAS}.annot"
    return fs_root / sess.subject / sess.session / "label" / f"{hemi}.{ATLAS}.annot"


def process_session(sess: Session, fs_root: Path, backfill_root: Path,
                    measures: tuple[str, ...], dk_codes: dict[str, str]
                    ) -> tuple[pd.DataFrame | None, list[dict]]:
    """All measures x {hcp, dsk} for one session: (long rows, QC rows per measure)."""
    t0 = time.time()
    sdir = fs_root / sess.subject / sess.session
    base = {"participant_id": sess.subject, "session_id": sess.session,
            "annot_source": "backfill" if sess.source == "backfill" else "fs_label",
            "annot_regenerated": False}
    qc = {m: {**base, "measure": m, "status": "ok", "missing_input": ""} for m in measures}

    def fail(status: str):
        for q in qc.values():
            q["status"] = status
        return None, list(qc.values())

    volumes = {}
    for src in sorted({MEASURES[m].source for m in measures if MEASURES[m].kind == "volume"}):
        p = sdir / src
        if p.exists():
            try:
                volumes[src] = load_volume(p)
            except Exception as exc:          # corrupt file: account for it, carry on
                for m in measures:
                    if MEASURES[m].source == src:
                        qc[m]["missing_input"] = f"{src} unreadable: {exc}"
        else:
            for m in measures:
                if MEASURES[m].source == src:
                    qc[m]["missing_input"] = src

    rows, summ = [], {}
    for hemi in HEMIS:
        try:
            vals, info = sample_hemi(sdir, hemi, measures, volumes)
        except Exception as exc:
            return fail(f"surface_error_{hemi}: {exc}")
        for m in measures:
            qc[m][f"normal_agree_{hemi}"] = round(info["normal_agree"], 4)
        for atlas in ("hcp", "dsk"):
            p = (annot_path(sess, hemi, fs_root, backfill_root) if atlas == "hcp"
                 else sdir / "label" / f"{hemi}.aparc.annot")
            if not p.exists():
                return fail(f"missing_annot_{atlas}_{hemi}")
            try:
                labels, tokens = annot_regions(p, atlas, dk_codes)
            except Exception as exc:
                return fail(f"annot_error_{atlas}_{hemi}: {exc}")
            if len(labels) != info["n"]:
                return fail(f"annot_size_{atlas}_{hemi}: {len(labels)} vs {info['n']}")
            for m in measures:
                if vals[m] is None:
                    if not qc[m]["missing_input"]:
                        qc[m]["missing_input"] = MEASURES[m].source.format(hemi=hemi)
                    continue
                prow, h = summarise(vals[m], info["area"], labels, tokens)
                rows += [{"atlas": atlas, "measure": m, "hemi": hemi, **r} for r in prow]
                summ[(atlas, m, hemi)] = h
                if atlas == "hcp":
                    qc[m][f"n_vertices_{hemi}"] = h["nverts"]
                    qc[m][f"n_nonfinite_{hemi}"] = h["nverts"] - h["n_finite"]

    for atlas, m in sorted({k[:2] for k in summ}):
        hs = {h: summ[(atlas, m, h)] for h in HEMIS if (atlas, m, h) in summ}
        if len(hs) == len(HEMIS):      # whole cortex only when both hemispheres sampled
            hs["all"] = {k: sum(x[k] for x in hs.values())
                         for k in ("nverts", "n_finite", "area", "wsum", "usum")}
        for hemi, h in hs.items():
            rows.append({"atlas": atlas, "measure": m, "hemi": hemi, "region": "",
                         "nverts": h["nverts"], "n_finite": h["n_finite"], "area": h["area"],
                         "mean_w": h["wsum"] / h["area"] if h["area"] else np.nan,
                         "mean_u": h["usum"] / h["n_finite"] if h["n_finite"] else np.nan})
    for m in measures:
        if qc[m]["missing_input"]:
            qc[m]["status"] = "missing_input"
        qc[m]["seconds"] = round(time.time() - t0, 2)
    if not rows:
        return None, list(qc.values())
    long = pd.DataFrame(rows)
    long.insert(0, "session_id", sess.session)
    long.insert(0, "participant_id", sess.subject)
    return long, list(qc.values())


def _worker(args):
    sess, fs_root, backfill_root, measures, dk_codes = args
    try:
        return process_session(sess, Path(fs_root), Path(backfill_root), measures, dk_codes)
    except Exception as exc:   # never lose a session silently
        return None, [{"participant_id": sess.subject, "session_id": sess.session,
                       "measure": m, "status": f"error: {exc!r}"} for m in measures]


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def all_sessions(parc_root: Path, backfill_root: Path) -> list[Session]:
    """The :func:`abcd.hcp_stats.discover_sessions` list the thickness table was built from."""
    return discover_sessions(parc_root, {"backfill": backfill_root})


def run_chunk(sessions: list[Session], part_name: str, out_dir: Path, fs_root: Path,
              backfill_root: Path, measures: tuple[str, ...], workers: int,
              log=print) -> Path:
    from .dk_stats import _region_to_code
    parts = out_dir / "surf_metrics_parts"
    parts.mkdir(parents=True, exist_ok=True)
    target = parts / f"{part_name}.parquet"
    if target.exists():
        log(f"{target} exists; skipping")
        return target
    t0 = time.time()
    codes = _region_to_code()
    tasks = [(s, str(fs_root), str(backfill_root), measures, codes) for s in sessions]
    longs, qcs = [], []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, (long, qc) in enumerate(ex.map(_worker, tasks, chunksize=1), 1):
            qcs += qc
            if long is not None:
                longs.append(long)
            if i % 50 == 0:
                log(f"  {i}/{len(sessions)}  {time.time() - t0:.0f}s")
    qc = pd.DataFrame(qcs)
    qc.to_csv(parts / f"{part_name}_qc.tsv", sep="\t", index=False)
    long = (pd.concat(longs, ignore_index=True) if longs
            else pd.DataFrame(columns=["participant_id", "session_id", "atlas", "measure",
                                       "hemi", "region", "nverts", "n_finite", "area",
                                       "mean_w", "mean_u"]))
    for c in ("participant_id", "session_id", "atlas", "measure", "hemi", "region"):
        long[c] = long[c].astype("category")
    tmp = target.with_suffix(".tmp")
    long.to_parquet(tmp, index=False)
    tmp.rename(target)        # the part file appears only once complete (restartable)
    log(f"{part_name}: {len(sessions)} sessions, status "
        f"{qc.status.str.split(':').str[0].value_counts().to_dict()}, "
        f"{time.time() - t0:.0f}s on {workers} workers")
    return target


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--parc-root", type=Path, default=DEFAULT_PARC_ROOT)
    ap.add_argument("--backfill-root", type=Path, default=DEFAULT_BACKFILL_ROOT)
    ap.add_argument("--fs-root", type=Path, default=DEFAULT_FS_ROOT)
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="default: <release 7.0 dir>/processed/hcp")
    ap.add_argument("--measures", default=",".join(DEFAULT_MEASURES))
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--chunk", type=int, default=None, help="chunk index (array task id)")
    ap.add_argument("--chunk-size", type=int, default=200)
    ap.add_argument("--sessions", type=Path, default=None,
                    help="TSV with participant_id, session_id: process only these")
    ap.add_argument("--part-name", default=None)
    ap.add_argument("--list-chunks", action="store_true",
                    help="print the number of chunks and exit")
    a = ap.parse_args(argv)

    out_dir = a.out_dir
    if out_dir is None:
        from . import paths
        out_dir = paths.release_dir("7.0") / "processed" / "hcp"
    measures = tuple(a.measures.split(","))
    unknown = set(measures) - set(MEASURES)
    if unknown:
        ap.error(f"unknown measures {sorted(unknown)}")

    sessions = all_sessions(a.parc_root, a.backfill_root)
    if a.list_chunks:
        print(-(-len(sessions) // a.chunk_size))
        return 0
    if a.sessions is not None:
        want = pd.read_csv(a.sessions, sep="\t")
        want = set(zip(want.participant_id, want.session_id))
        sessions = [s for s in sessions if (s.subject, s.session) in want]
        part = a.part_name or a.sessions.stem
    elif a.chunk is not None:
        sessions = sessions[a.chunk * a.chunk_size:(a.chunk + 1) * a.chunk_size]
        part = a.part_name or f"part_{a.chunk:04d}"
    else:
        ap.error("give --chunk or --sessions")
    if not sessions:
        print("no sessions in this chunk")
        return 0
    run_chunk(sessions, part, out_dir, a.fs_root, a.backfill_root, measures, a.workers,
              log=lambda m: print(m, flush=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
