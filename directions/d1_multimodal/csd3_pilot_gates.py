"""
D1 step 2 gates on the 20 pilot sessions (run by csd3_pilot.sbatch, after abcd.surf_metrics).

A.  own DK GWC (unweighted vertex mean) vs FreeSurfer stats/?h.w-g.pct.stats ``Mean``;
    A'. own HCP GWC vs rr480's ``?h_HCP.fsaverage.aparc.w-g.pct.stats`` where present.
B.  own HCP per-parcel vertex counts vs ``hcp_anatomical_stats_long.parquet`` (the .log tables).
C.  own per-vertex T1 samples vs FreeSurfer ``mri_vol2surf`` (lh, 3 sessions).
D.  own DK ``t1__gm50`` vs release ``mr_y_smri__t1__gm__dsk``, per region across sessions.
E.  seconds per session; projected core-hours for all sessions; T2-registration pilot timing.

Prints everything and writes ``<release>/processed/hcp/surf_metrics_pilot/gates.txt``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from abcd import paths
from abcd.dk_stats import _region_to_code
from abcd.hcp_stats import DEFAULT_FS_ROOT, DEFAULT_PARC_ROOT, region_from_structure
from abcd.surf_metrics import load_volume, sample_hemi

REL = paths.release_dir("7.0")
HCP = REL / "processed" / "hcp"
WORK = HCP / "surf_metrics_pilot"
N_ALL = 33795
OUT: list[str] = []


def say(s: str = "") -> None:
    print(s, flush=True)
    OUT.append(s)


def segstats_means(path: Path) -> pd.Series:
    """StructName -> Mean from an mri_segstats --sum table."""
    cols = None
    rows = {}
    for line in open(path):
        if line.startswith("# ColHeaders"):
            cols = line.split()[2:]
        elif not line.startswith("#") and line.strip():
            p = line.split()
            rows[p[cols.index("StructName")]] = float(p[cols.index("Mean")])
    return pd.Series(rows)


def corr_summary(x: pd.Series, y: pd.Series) -> str:
    d = (x - y).abs()
    return f"n={len(x)}, r={np.corrcoef(x, y)[0, 1]:.6f}, max|d|={d.max():.4f}, mean|d|={d.mean():.4f}"


def main() -> int:
    pilot = pd.read_csv(HCP / "surf_metrics_pilot_sessions.tsv", sep="\t")
    long = pd.read_parquet(HCP / "surf_metrics_parts" / "pilot.parquet")
    for c in ("participant_id", "session_id", "atlas", "measure", "hemi", "region"):
        long[c] = long[c].astype(str)
    qc = pd.read_csv(HCP / "surf_metrics_parts" / "pilot_qc.tsv", sep="\t")
    say(f"pilot: {pilot.shape[0]} sessions; surf_metrics status "
        f"{qc.status.value_counts().to_dict()}")
    say(pilot.groupby("session_id").size().to_string())
    codes = _region_to_code()

    # ---- A
    say("\n## Gate A: DK GWC vs FreeSurfer stats/?h.w-g.pct.stats (unweighted means)")
    own = long[(long.atlas == "dsk") & (long.measure == "gwc") & (long.region != "")]
    a, a2 = [], []
    for (sub, ses), g in own.groupby(["participant_id", "session_id"]):
        for hemi in ("lh", "rh"):
            fs = segstats_means(DEFAULT_FS_ROOT / sub / ses / "stats" / f"{hemi}.w-g.pct.stats")
            fs = fs[fs.index.isin(codes)].rename(index=codes)
            m = g[g.hemi == hemi].set_index("region").mean_u
            a.append(pd.DataFrame({"own": m, "fs": fs}).dropna())
            hp = (DEFAULT_PARC_ROOT / sub / ses / "HCP.fsaverage.aparc"
                  / f"{hemi}_HCP.fsaverage.aparc.w-g.pct.stats")
            if hp.exists():
                fh = segstats_means(hp)
                fh.index = [region_from_structure(i) for i in fh.index]
                fh = fh[fh.index.notna()]
                mh = long[(long.atlas == "hcp") & (long.measure == "gwc") & (long.region != "")
                          & (long.participant_id == sub) & (long.session_id == ses)
                          & (long.hemi == hemi)].set_index("region").mean_u
                a2.append(pd.DataFrame({"own": mh, "fs": fh}).dropna())
    a = pd.concat(a)
    passed = np.corrcoef(a.own, a.fs)[0, 1] > 0.999 and (a.own - a.fs).abs().max() < 0.5
    say(f"DK: {corr_summary(a.own, a.fs)}  -> {'PASS' if passed else 'FAIL'}")
    if a2:
        a2 = pd.concat(a2)
        say(f"A' HCP vs rr480 HCP w-g.pct.stats ({len(a2) // 180} hemispheres): "
            f"{corr_summary(a2.own, a2.fs)}")

    # ---- B
    say("\n## Gate B: HCP vertex counts vs .log tables (hcp_anatomical_stats_long.parquet)")
    keys = list(zip(pilot.participant_id, pilot.session_id))
    ref = pd.read_parquet(HCP / "hcp_anatomical_stats_long.parquet",
                          columns=["participant_id", "session_id", "hemi", "region", "nverts"],
                          filters=[("participant_id", "in", sorted(pilot.participant_id.unique()))])
    ref = ref[[k in set(keys) for k in zip(ref.participant_id, ref.session_id)]]
    mine = long[(long.atlas == "hcp") & (long.measure == "gwc") & (long.region != "")]
    b = mine.merge(ref, on=["participant_id", "session_id", "hemi", "region"],
                   how="outer", suffixes=("_own", "_log"), indicator=True)
    nmis = (b._merge != "both").sum()
    ndiff = (b.nverts_own != b.nverts_log).sum()
    say(f"{len(b)} parcel rows; unmatched {nmis}; vertex-count mismatches {ndiff} "
        f"-> {'PASS' if nmis == 0 and ndiff == 0 else 'FAIL'}")

    # ---- C
    say("\n## Gate C: per-vertex T1 (nu.mgz) vs FreeSurfer mri_vol2surf, lh")
    for f in sorted(WORK.glob("*.lh.nu.projfrac0.5.mgh")):
        s = f.name.split(".")[0]
        sub, ses = s.rsplit("_", 1)
        sdir = DEFAULT_FS_ROOT / sub / ses
        vals, _ = sample_hemi(sdir, "lh", ("t1__gm50", "t1__wm1"),
                              {"mri/nu.mgz": load_volume(sdir / "mri" / "nu.mgz")})
        for meas, fsf, how in (("t1__gm50", f, "projfrac 0.5 (along white normal; not the same point)"),
                               ("t1__gm50", WORK / f"{s}.lh.nu.mid.mgh", "--surf mid --projdist 0 (same definition)"),
                               ("t1__wm1", WORK / f"{s}.lh.nu.projdist-1.mgh", "projdist -1")):
            if not fsf.exists():
                continue
            fsv = np.asarray(nib.load(str(fsf)).dataobj).reshape(-1)
            ok = np.isfinite(vals[meas]) & np.isfinite(fsv) & (fsv > 0)
            r = np.corrcoef(vals[meas][ok], fsv[ok])[0, 1]
            say(f"{s} {meas}: n={ok.sum()}, r={r:.5f}, "
                f"median|d|={np.median(np.abs(vals[meas][ok] - fsv[ok])):.3f}  [{how}]")

    # ---- D
    say("\n## Gate D: own DK t1__gm50 vs release mr_y_smri__t1__gm__dsk")
    rt = pd.read_csv(REL / "y" / "mr_y_smri__t1__gm__dsk.tsv", sep="\t")
    rt = rt.merge(pilot[["participant_id", "session_id"]])
    say(f"{len(rt)} pilot sessions present in the release table")
    own = long[(long.atlas == "dsk") & (long.measure == "t1__gm50") & (long.region != "")]
    own = own.assign(col="mr_y_smri__t1__gm__dsk__" + own.region + "__" + own.hemi + "_mean")
    rl = rt.melt(id_vars=["participant_id", "session_id"], var_name="col", value_name="rel")
    d = own.merge(rl, on=["participant_id", "session_id", "col"]).dropna(subset=["rel"])
    per = d.groupby("col").apply(lambda g: np.corrcoef(g.mean_w, g.rel)[0, 1] if len(g) > 3 else np.nan,
                                 include_groups=False)
    say(f"per-region r across sessions (68 regions): median {per.median():.3f}, "
        f"IQR {per.quantile(.25):.3f}-{per.quantile(.75):.3f}, min {per.min():.3f}, "
        f"share > 0.8: {(per > 0.8).mean():.2f}")
    say(f"pooled region x session: {corr_summary(d.mean_w, d.rel)}")
    wc = long[(long.atlas == "dsk") & (long.measure == "t1__gm50") & (long.hemi == "all")]
    wc = wc.merge(rt[["participant_id", "session_id", "mr_y_smri__t1__gm__dsk_mean"]])
    say(f"whole-cortex mean across sessions: r={np.corrcoef(wc.mean_w, wc.mr_y_smri__t1__gm__dsk_mean)[0, 1]:.3f}"
        f" (own median {wc.mean_w.median():.1f}, release median {wc.mr_y_smri__t1__gm__dsk_mean.median():.1f})")
    say("no mr_y_dti__* / mr_y_rsi__* tables on CSD3: diffusion part of gate D not applicable")

    # ---- E
    say("\n## Gate E: timing")
    sec = qc.groupby(["participant_id", "session_id"]).seconds.first()
    say(f"surf_metrics (GWC + T1 gm50 + T1 wm1, HCP + DK, both hemis): {sec.mean():.1f} s/session "
        f"(median {sec.median():.1f}, max {sec.max():.1f}) on 4 workers sharing one node; "
        f"projected {sec.mean() * N_ALL / 3600:.0f} core-hours for {N_ALL} sessions")
    t2 = WORK / "t2reg_pilot.tsv"
    if t2.exists():
        t = pd.read_csv(t2, sep="\t")
        say("\n## T2 registration pilot (fMRIPrep T1w -> DAIRC FreeSurfer; mri_coreg + bbregister --t1)")
        say(t.status.value_counts().to_string())
        ok = t[t.status == "ok"]
        if len(ok):
            tot = ok.sec_coreg + ok.sec_bbregister
            say(f"seconds/session: coreg {ok.sec_coreg.median():.0f}, bbregister {ok.sec_bbregister.median():.0f}, "
                f"total median {tot.median():.0f} (max {tot.max():.0f})")
            say(f"bbregister min cost: median {ok.bbr_mincost.median():.3f}, max {ok.bbr_mincost.max():.3f}")
            say(f"projected {tot.mean() * N_ALL / 3600:.0f} core-hours for {N_ALL} sessions")
    (WORK / "gates.txt").write_text("\n".join(OUT) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
