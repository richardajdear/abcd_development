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
``t2__gm50``   T2w intensity at GM mid-depth                   fMRIPrep ``desc-preproc_T2w`` (registered)
``t2__wm1``    T2w intensity 1 mm below the white surface      fMRIPrep ``desc-preproc_T2w`` (registered)
``t1t2__gm50`` T1w/T2w ratio at GM mid-depth, per vertex       fMRIPrep ``desc-preproc_T1w`` / ``_T2w``
=============  =============================================  ==============================

Volume measures are sampled with trilinear interpolation
(``scipy.ndimage.map_coordinates``, order 1) at tkr-RAS vertex positions mapped to
voxels with the inverse of the volume's ``vox2ras-tkr``.  fMRIPrep volumes (``fmriprep:``
sources) live in fMRIPrep's own T1w space; they are reached through a per-session
``bbregister`` LTA from the fMRIPrep T1w to the DAIRC ``orig.mgz`` (``<reg_root>/<sub>/<ses>/
bbr.lta``, made by ``tools/t2reg.sbatch``).  The fMRIPrep T2w is already resampled into that
T1w space, so one transform serves both.  fMRIPrep N4-corrects its T1w but not its T2w, so
the T2w samples are divided by an N4 field computed the same way (``t2_biasfield_8mm.nii.gz``,
also from ``tools/t2reg.sbatch``).  The ratio is formed per vertex before averaging, from the
two (both N4-corrected) fMRIPrep images.  The white-surface normal is the
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

    sbatch --array=0-168 tools/surf_metrics.sbatch           # all sessions, 200 per task
    sbatch --array=0-168 tools/surf_metrics.sbatch --measures t2__gm50,t2__wm1,t1t2__gm50 \
        --part-prefix part_t2                                # after tools/t2reg.sbatch
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
DEFAULT_FMRIPREP_ROOT = Path("/rds/project/rds-CeXlNYOYMxw/derivatives/abcc-fmriprep_v25.1.4")
DEFAULT_REG_ROOT = Path("/rds/user/rajd2/hpc-work/abcd_development/legacy/hpc/work/d1_t2reg")
#: fMRIPrep images that are NOT bias-corrected by fMRIPrep -> N4 field written by
#: tools/t2reg.sbatch next to bbr.lta (fMRIPrep's T1w is already N4-corrected)
BIASFIELD = {"T2w": "t2_biasfield_8mm.nii.gz"}
#: bbregister --t1 costs are bimodal over the 28.8k sessions: 99 % < 0.065, then a gap from
#: 0.3 to 0.8 and ~75 sessions at 0.85-1.0 that fail from every initialisation tried (mri_coreg
#: default/fast, --init-header).  Above this the fMRIPrep volumes are treated as unregistered.
MAX_BBR_COST = 0.5
#: fMRIPrep resamples each image into its T1w grid with zero fill, so cortex outside an image's
#: field of view samples ~0 (seen for the T2w: up to 11 % of vertices in a session, ratios up to
#: 1e8).  Samples of an fMRIPrep source below this fraction of the hemisphere's median gm50
#: sample are set to NaN -- far below any tissue (CSF is brighter, WM ~0.8 x GM in T2w).
FOV_FLOOR = 0.2
HEMIS = ("lh", "rh")
DK_EXCLUDE = {"unknown", "corpuscallosum"}
HCP_MEDIAL_WALL = "???"


@dataclass(frozen=True)
class Measure:
    name: str               # column infix after mr_y_smri__, e.g. "t1__gm50"
    kind: str               # "surface" (per-vertex file), "volume" (sampled) or "ratio"
    source: str             # path under the FreeSurfer session ({hemi} filled in), or
                            # "fmriprep:<suffix>"; a ratio is "<numerator>/<denominator>"
    depth: str | None = None  # gm50 / wm1 / wm2 for volume measures

    @property
    def volume_sources(self) -> list[str]:
        return {"surface": [], "volume": [self.source], "ratio": self.source.split("/")}[self.kind]


MEASURES = {
    m.name: m for m in (
        Measure("gwc", "surface", "surf/{hemi}.w-g.pct.mgh"),
        Measure("t1__gm50", "volume", "mri/nu.mgz", "gm50"),
        Measure("t1__wm1", "volume", "mri/nu.mgz", "wm1"),
        Measure("t2__gm50", "volume", "fmriprep:T2w", "gm50"),
        Measure("t2__wm1", "volume", "fmriprep:T2w", "wm1"),
        Measure("t1t2__gm50", "ratio", "fmriprep:T1w/fmriprep:T2w", "gm50"),
    )
}
#: FreeSurfer-only measures (the first pass); the fMRIPrep ones need tools/t2reg.sbatch first
DEFAULT_MEASURES = ("gwc", "t1__gm50", "t1__wm1")
T2_MEASURES = ("t2__gm50", "t2__wm1", "t1t2__gm50")


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


def sample_volume(data: np.ndarray, tkr2vox: np.ndarray, xyz: np.ndarray,
                  mode: str = "constant") -> np.ndarray:
    """Trilinear samples of ``data`` at tkr-RAS points ``xyz`` (n x 3); NaN outside."""
    ijk = xyz @ tkr2vox[:3, :3].T + tkr2vox[:3, 3]
    return map_coordinates(data, ijk.T, order=1, mode=mode, cval=np.nan)


def sample_source(vol: tuple, xyz: np.ndarray) -> np.ndarray:
    """Sample a loaded source: ``(data, tkr2vox)`` or, bias-corrected,
    ``(data, tkr2vox, field, field_tkr2vox)`` -> data / field (N4 convention)."""
    v = sample_volume(vol[0], vol[1], xyz)
    if len(vol) == 4:          # smooth, coarse field: extend at the grid edge
        with np.errstate(invalid="ignore", divide="ignore"):
            v = v / sample_volume(vol[2], vol[3], xyz, mode="nearest")
    return v


def load_volume(path: Path) -> tuple[np.ndarray, np.ndarray]:
    img = nib.load(str(path))
    data = np.asarray(img.dataobj, dtype=np.float32)
    return data, np.linalg.inv(img.header.get_vox2ras_tkr())


def read_lta(path: Path) -> tuple[int, np.ndarray]:
    """(type, 4x4 matrix) of a single-transform FreeSurfer LTA (0 = vox2vox, 1 = ras2ras)."""
    lines = open(path).read().splitlines()
    ltype = next(int(l.split("=")[1].split()[0]) for l in lines if l.strip().startswith("type"))
    i = next(k for k, l in enumerate(lines) if l.split() == ["1", "4", "4"])
    return ltype, np.array([[float(x) for x in lines[i + 1 + r].split()] for r in range(4)])


def fmriprep_image(fmriprep_root: Path, sub: str, ses: str, suffix: str) -> Path | None:
    """``*_desc-preproc_<suffix>.nii.gz``; an empty suffix gives the ``desc-brain_mask``."""
    pat = f"*_desc-preproc_{suffix}.nii.gz" if suffix else "*_desc-brain_mask.nii.gz"
    hits = sorted((fmriprep_root / sub / ses / "anat").glob(pat))
    return hits[0] if hits else None


def load_biasfield(field: Path, mov: Path, mask: Path | None, tkr2vox: np.ndarray
                   ) -> tuple[np.ndarray, np.ndarray]:
    """N4 field (coarse grid) and its tkr-RAS(orig) -> field-voxel matrix.

    N4 does not normalise its field, so dividing by it would rescale the image by an
    arbitrary per-session factor.  The field is divided by its median inside the fMRIPrep
    brain mask, so the corrected image keeps the original intensity scale and only the
    field's spatial shape is removed.
    """
    f = nib.load(str(field))
    data = np.asarray(f.dataobj, dtype=np.float32)
    mov_affine = nib.load(str(mov)).affine
    if mask is not None:
        m = nib.load(str(mask))
        ijk = np.argwhere(np.asarray(m.dataobj) > 0)[::27]
        fijk = (np.linalg.inv(f.affine) @ m.affine @ np.c_[ijk, np.ones(len(ijk))].T)[:3]
        data = data / np.median(map_coordinates(data, fijk, order=1, mode="nearest"))
    return data, np.linalg.inv(f.affine) @ mov_affine @ tkr2vox


def load_registered(mov: Path, lta: Path, orig: Path) -> tuple[np.ndarray, np.ndarray]:
    """A volume in fMRIPrep T1w space, with the tkr-RAS(orig) -> mov-voxel matrix.

    The LTA maps the fMRIPrep T1w (src) onto ``orig.mgz`` (dst).  Any fMRIPrep image in the
    same T1w space (the T2w is resampled there) shares the voxel grid of the T1w, which is
    checked against the LTA's src geometry.
    """
    img = nib.load(str(mov))
    oimg = nib.load(str(orig))
    ltype, m = read_lta(lta)
    tkr2orig = np.linalg.inv(oimg.header.get_vox2ras_tkr())
    if ltype == 0:        # mov vox -> orig vox
        tkr2vox = np.linalg.inv(m) @ tkr2orig
    elif ltype == 1:      # mov scanner RAS -> orig scanner RAS
        tkr2vox = np.linalg.inv(img.affine) @ np.linalg.inv(m) @ oimg.affine @ tkr2orig
    else:
        raise ValueError(f"{lta}: unsupported LTA type {ltype}")
    src = next(l for l in open(lta) if l.startswith("volume ="))
    if [int(x) for x in src.split("=")[1].split()] != list(img.shape[:3]):
        raise ValueError(f"{mov.name}: grid {img.shape[:3]} differs from LTA src {src.strip()}")
    return np.asarray(img.dataobj, dtype=np.float32), tkr2vox


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
    floors = {}                       # fMRIPrep source -> out-of-field-of-view threshold
    for name in measures:
        m = MEASURES[name]
        if m.kind == "surface":
            p = sdir / m.source.format(hemi=hemi)
            out[name] = (np.asarray(nib.load(str(p)).dataobj, dtype=np.float64).reshape(-1)
                         if p.exists() else None)
            if out[name] is not None and len(out[name]) != len(white):
                raise ValueError(f"{p.name}: {len(out[name])} values for {len(white)} vertices")
        else:
            vols = [volumes.get(s) for s in m.volume_sources]
            if any(v is None for v in vols):
                out[name] = None
                continue
            xyz = depth_coords(m.depth, white, pial, normals)
            vals = [sample_source(v, xyz) for v in vols]
            for i, src in enumerate(m.volume_sources):
                if src.startswith("fmriprep:"):
                    if src not in floors:
                        gm = sample_source(volumes[src], depth_coords("gm50", white, pial, normals))
                        floors[src] = FOV_FLOOR * np.nanmedian(gm)
                    vals[i] = np.where(vals[i] >= floors[src], vals[i], np.nan)
            if m.kind == "ratio":
                with np.errstate(invalid="ignore", divide="ignore"):
                    out[name] = vals[0] / vals[1]
            else:
                out[name] = vals[0]
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
                    measures: tuple[str, ...], dk_codes: dict[str, str],
                    fmriprep_root: Path = DEFAULT_FMRIPREP_ROOT, reg_root: Path = DEFAULT_REG_ROOT
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
    for src in sorted({s for m in measures for s in MEASURES[m].volume_sources}):
        users = [m for m in measures if src in MEASURES[m].volume_sources]
        try:
            if src.startswith("fmriprep:"):
                p = fmriprep_image(fmriprep_root, sess.subject, sess.session, src.split(":")[1])
                lta = reg_root / sess.subject / sess.session / "bbr.lta"
                field = lta.with_name(BIASFIELD.get(src.split(":")[1], "-"))
                cost_file = lta.with_name("bbr.dat.mincost")
                cost = (float(cost_file.read_text().split()[0]) if cost_file.exists() else np.nan)
                for m in users:
                    qc[m]["bbr_mincost"] = cost
                missing = (f"fmriprep {src.split(':')[1]}" if p is None
                           else f"registration {lta.name}" if not lta.exists()
                           else f"registration failed (bbr cost {cost:.3f} > {MAX_BBR_COST})"
                           if cost > MAX_BBR_COST
                           else f"bias field {field.name}" if src.split(":")[1] in BIASFIELD
                           and not field.exists() else None)
                if missing is None:
                    volumes[src] = load_registered(p, lta, sdir / "mri" / "orig.mgz")
                    if src.split(":")[1] in BIASFIELD:
                        mask = fmriprep_image(fmriprep_root, sess.subject, sess.session, "")
                        volumes[src] = (*volumes[src], *load_biasfield(field, p, mask, volumes[src][1]))
            else:
                p = sdir / src
                missing = None if p.exists() else src
                if missing is None:
                    volumes[src] = load_volume(p)
        except Exception as exc:              # corrupt file: account for it, carry on
            missing = f"{src} unreadable: {exc}"
        if missing:
            for m in users:
                qc[m]["missing_input"] = qc[m]["missing_input"] or missing

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
    sess, fs_root, backfill_root, measures, dk_codes, fmriprep_root, reg_root = args
    try:
        return process_session(sess, Path(fs_root), Path(backfill_root), measures, dk_codes,
                               Path(fmriprep_root), Path(reg_root))
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
              fmriprep_root: Path = DEFAULT_FMRIPREP_ROOT, reg_root: Path = DEFAULT_REG_ROOT,
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
    tasks = [(s, str(fs_root), str(backfill_root), measures, codes, str(fmriprep_root),
              str(reg_root)) for s in sessions]
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
    ap.add_argument("--fmriprep-root", type=Path, default=DEFAULT_FMRIPREP_ROOT)
    ap.add_argument("--reg-root", type=Path, default=DEFAULT_REG_ROOT,
                    help="per-session bbr.lta (fMRIPrep T1w -> orig.mgz), from tools/t2reg.sbatch")
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
    ap.add_argument("--part-prefix", default="part",
                    help="chunk part files are <prefix>_<chunk>; use e.g. part_t2 for a second pass")
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
        part = a.part_name or f"{a.part_prefix}_{a.chunk:04d}"
    else:
        ap.error("give --chunk or --sessions")
    if not sessions:
        print("no sessions in this chunk")
        return 0
    run_chunk(sessions, part, out_dir, a.fs_root, a.backfill_root, measures, a.workers,
              a.fmriprep_root, a.reg_root,
              log=lambda m: print(m, flush=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
