"""
D1 step 0 on CSD3: inventory of the inputs that surf_metrics needs, per session.

For every session found by the ``abcd.hcp_stats`` discovery (primary tree +
backfill), record which FreeSurfer files exist, whether any T2 is present in
the FreeSurfer session, whether ``recon-all`` used ``-T2pial``, where the
subject-space HCP annotation is, and which ABCC derivatives (fMRIPrep T2w,
qsirecon tensor maps) exist for the same session.  Also picks the 20 pilot
sessions (5 per visit, distinct sites, all three manufacturers).

    sbatch --time=00:40:00 --cpus-per-task=16 tools/run_python.sbatch directions/d1_multimodal/csd3_discovery.py

Writes ``<release>/processed/hcp/surf_metrics_inventory.tsv`` and
``surf_metrics_pilot_sessions.txt``; prints the counts the discovery report needs.
"""

from __future__ import annotations

import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from abcd import paths
from abcd.hcp_stats import ATLAS, DEFAULT_FS_ROOT, DEFAULT_PARC_ROOT, discover_sessions

BACKFILL = Path("/rds/user/rajd2/hpc-work/abcd_development/legacy/hpc/work/parcellations_backfill/T1")
DERIV = Path("/rds/project/rds-CeXlNYOYMxw/derivatives")
FMRIPREP = DERIV / "abcc-fmriprep_v25.1.4"
QSIPREP = DERIV / "abcc-qsiprep_v0.21.4"
DKI = DERIV / "abcc-qsirecon-DIPYDKI"
GQI = DERIV / "abcc-qsirecon-DSIStudioGQI"

FS_FILES = {
    "surf": ["lh.white", "rh.white", "lh.pial", "rh.pial", "lh.area", "rh.area",
             "lh.w-g.pct.mgh", "rh.w-g.pct.mgh"],
    "stats": ["lh.w-g.pct.stats", "rh.w-g.pct.stats"],
    "mri": ["orig.mgz", "nu.mgz", "T1.mgz", "rawavg.mgz"],
    "label": ["lh.aparc.annot", "rh.aparc.annot",
              f"lh.{ATLAS}.annot", f"rh.{ATLAS}.annot"],
}
T2_NAME = re.compile(r"t2", re.IGNORECASE)


def _ls(p: Path) -> set[str]:
    try:
        return set(os.listdir(p))
    except OSError:
        return set()


def inventory(args) -> dict:
    sub, ses, source = args
    d = DEFAULT_FS_ROOT / sub / ses
    row = {"participant_id": sub, "session_id": ses, "parc_source": source}
    for sd, names in FS_FILES.items():
        have = _ls(d / sd)
        for n in names:
            row[f"fs_{sd}/{n}"] = n in have
    # any T2 anywhere in the session (names only; recon-all logs excluded)
    t2 = []
    for root, _dirs, files in os.walk(d):
        if "/scripts" in root:
            continue
        t2 += [os.path.relpath(os.path.join(root, f), d) for f in files if T2_NAME.search(f)]
    row["fs_t2_files"] = ";".join(sorted(t2))
    try:
        cmd = next((l for l in open(d / "scripts" / "recon-all.done") if l.startswith("CMDARGS")), "")
    except OSError:
        cmd = ""
    row["fs_T2pial"] = "-T2pial" in cmd
    row["fs_cmdargs_found"] = bool(cmd)
    bf = BACKFILL / sub / ses / ATLAS
    row["annot_backfill"] = all((bf / f"{h}.{ATLAS}.annot").exists() for h in ("lh", "rh"))
    # ABCC derivatives for the same session
    anat = _ls(FMRIPREP / sub / ses / "anat")
    row["fmriprep_T2w"] = any(f.endswith("desc-preproc_T2w.nii.gz") for f in anat)
    row["fmriprep_T1w"] = any(f.endswith("desc-preproc_T1w.nii.gz") for f in anat)
    row["fmriprep_to_fsnative_xfm"] = any("from-T1w_to-fsnative" in f for f in anat)
    row["qsiprep_dwi"] = any(f.endswith("space-T1w_desc-preproc_dwi.nii.gz")
                             for f in _ls(QSIPREP / sub / ses / "dwi"))
    dki = _ls(DKI / sub / ses / "dwi")
    row["dki_tensor_fa_T1w"] = any("space-T1w_model-tensor_param-fa_dwimap.nii.gz" in f for f in dki)
    row["dki_md_T1w"] = any("space-T1w_model-dki_param-md_dwimap.nii.gz" in f for f in dki)
    gqi = _ls(GQI / sub / ses / "dwi")
    row["gqi_tensor_fa_T1w"] = any("space-T1w_model-tensor_param-fa_dwimap.nii.gz" in f for f in gqi)
    row["gqi_tensor_md_T1w"] = any("space-T1w_model-tensor_param-md_dwimap.nii.gz" in f for f in gqi)
    return row


def pick_pilot(sessions: pd.DataFrame, rel: Path, n_per_visit: int = 5) -> pd.DataFrame:
    dyn = pd.read_csv(rel / "g" / "dyn" / "ab_g_dyn.tsv", sep="\t",
                      usecols=["participant_id", "session_id", "ab_g_dyn__design_site",
                               "ab_g_dyn__design_mr__manufact"])
    s = sessions.merge(dyn, on=["participant_id", "session_id"])
    out = []
    for ses, g in s.groupby("session_id"):
        g = g.sample(frac=1, random_state=20261027).drop_duplicates("ab_g_dyn__design_site")
        # round-robin over manufacturers so each visit covers all of them
        g = g.assign(k=g.groupby("ab_g_dyn__design_mr__manufact").cumcount()).sort_values("k")
        out.append(g.head(n_per_visit))
    return pd.concat(out).drop(columns="k")


def main() -> int:
    rel = paths.release_dir("7.0")
    out_dir = rel / "processed" / "hcp"
    sess = discover_sessions(DEFAULT_PARC_ROOT, {"backfill": BACKFILL})
    print(f"{len(sess)} sessions from hcp_stats discovery", flush=True)
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
    with ProcessPoolExecutor(workers) as ex:
        rows = list(ex.map(inventory, [(s.subject, s.session, s.source) for s in sess], chunksize=64))
    inv = pd.DataFrame(rows)
    inv.to_csv(out_dir / "surf_metrics_inventory.tsv", sep="\t", index=False)

    flags = [c for c in inv.columns if inv[c].dtype == bool]
    print("\n== counts (all sessions)")
    print(inv[flags].sum().to_string())
    print("\n== by visit")
    print(inv.groupby("session_id")[flags].sum().T.to_string())
    print("\n== annotation used by surf_metrics")
    prim = inv.parc_source == "primary"
    print(f"primary sessions {prim.sum()}, with FS label HCP annot both hemis "
          f"{(inv.loc[prim, f'fs_label/lh.{ATLAS}.annot'] & inv.loc[prim, f'fs_label/rh.{ATLAS}.annot']).sum()}")
    print(f"backfill sessions {(~prim).sum()}, with backfill annot {inv.loc[~prim, 'annot_backfill'].sum()}")
    print("\n== T2 file names in FreeSurfer sessions (top 20)")
    t2 = inv.fs_t2_files.str.split(";").explode()
    print(t2[t2 != ""].value_counts().head(20).to_string())
    print(f"sessions with any T2 file: {(inv.fs_t2_files != '').sum()}")

    pf = out_dir / "surf_metrics_pilot_sessions.tsv"
    if pf.exists():      # already picked (the pilot may have run before this inventory)
        pilot = pd.read_csv(pf, sep="\t")
    else:
        pilot = pick_pilot(inv[["participant_id", "session_id", "parc_source"]], rel)
        pilot.to_csv(pf, sep="\t", index=False)
    print("\n== pilot sessions")
    print(pilot.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
