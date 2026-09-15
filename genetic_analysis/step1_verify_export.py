#!/usr/bin/env python
"""Gate the CSD3-regenerated 7.0 phenotype export against the committed local run.

README_HPC.md Step 0 assumes the export is rsynced from the laptop.  It was
instead regenerated on CSD3 (01_export_pheno.sbatch) from the true 7.0 release
tables, so it has to be proved identical in sample construction to the local
re-run whose numbers are published in docs/vintage_comparison.csv and
README.md "Status".

Counts and the age centre are deterministic given the tables, the QC stack and
the sample rule, so they must match exactly.  The BLUPs come out of lme4 and
may differ in the last digits between machines, so global_slope is checked by
its correlation against the archived 6.0-vintage export on shared children --
the r = 0.9212 that docs/vintage_comparison.csv records.

    python genetic_analysis/step1_verify_export.py

Exit status 0 = every check passed and the export may be used downstream.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import argparse

REPO = Path(__file__).resolve().parent.parent
DK_RUN = REPO / "out" / "thickness_dsk_70_139406217085"
NEW = DK_RUN
OLD = REPO / "out" / "legacy_6.0_tabulated" / "thickness_dsk_70_139406217085"
VINTAGE = REPO / "docs" / "vintage_comparison.csv"

EXPECTED_PHENOTYPES = [
    "baseline_thickness", "global_slope", "slope_PC3", "slope_PC2", "slope_PC1",
]


def expected() -> dict[str, float]:
    v = pd.read_csv(VINTAGE)
    v = v[v.vintage == "7.0_tables"]
    out = dict(zip(v.quantity, v.value))
    both = pd.read_csv(VINTAGE)
    both = both[both.vintage == "both"]
    out.update(dict(zip(both.quantity, both.value)))
    return out


def main_hcp(run_dir: Path) -> int:
    """HCP-MMP gate: the DK export is the reference, not docs/vintage_comparison.

    Same children, same scans, same QC stack, same model -- only the atlas
    differs -- so sample construction must match the DK run almost exactly
    (the HCP table has 33,792 of the 33,794 DK sessions) and the whole-cortex
    phenotypes must correlate very highly with DK's on shared children.  The
    regions are 358 (360 minus the excluded hippocampal parcel H, both hemis).
    """
    checks = []
    def chk(name, got, want, ok): checks.append((name, got, want, bool(ok)))
    man = json.loads((run_dir / "manifest.json").read_text())
    dk = json.loads((DK_RUN / "manifest.json").read_text())
    fin, dfin = man["final"], dk["final"]
    chk("release", man["release"], "7.0", man["release"] == "7.0")
    chk("parcellation", man["config"]["parcellation"], "hcp", man["config"]["parcellation"] == "hcp")
    chk("exclude_regions", man["config"].get("exclude_regions"), ["H"], man["config"].get("exclude_regions") == ["H"])
    chk("n_regions", fin["n_regions"], 358, fin["n_regions"] == 358)
    for k in ("n_subjects", "n_scans", "n_families"):
        chk(f"{k} vs DK", fin[k], f"within 1% of {dfin[k]}", abs(fin[k] - dfin[k]) <= 0.01 * dfin[k])
    chk("age_centre vs DK", round(man["age_centre"], 3), round(dk["age_centre"], 3),
        abs(man["age_centre"] - dk["age_centre"]) < 0.02)
    gi = run_dir / "gcta_inputs" / "phenotypes_gcta.txt"
    if not gi.exists():
        chk("gcta_inputs/phenotypes_gcta.txt", "MISSING", "present", False)
    else:
        ph = pd.read_csv(gi, sep=r"\s+")
        cols = [c for c in ph.columns if c not in ("FID", "IID")]
        chk("export phenotypes", cols, EXPECTED_PHENOTYPES, cols == EXPECTED_PHENOTYPES)
        dkp = pd.read_csv(DK_RUN / "gcta_inputs" / "phenotypes_gcta.txt", sep=r"\s+")
        m = ph.merge(dkp, on="IID", suffixes=("_hcp", "_dk"))
        chk("subjects shared with DK", len(m), f">= 99% of {len(dkp)}", len(m) >= 0.99 * len(dkp))
        # Thresholds (set 2026-09-15 after inspecting the first HCP export):
        # the two whole-cortex means are different weightings of the same
        # cortex (68 large DK regions vs 358 small HCP parcels), so identity is
        # not expected.  Observed: baseline 0.978, global_slope 0.948, and
        # global_slope agreement RISES with visits (0.942 -> 0.953 for 2 -> 4),
        # i.e. the gap is slope measurement noise, not a systematic
        # difference.  0.90 is the floor below which a bug (wrong sessions,
        # wrong QC, a zero parcel leaking in) becomes the likelier reading; for
        # scale, the 6.0 -> 7.0 change of the SAME DK phenotype was r = 0.92.
        for c, lo in (("global_slope", 0.90), ("baseline_thickness", 0.95)):
            r = float(np.corrcoef(m[f"{c}_hcp"], m[f"{c}_dk"])[0, 1])
            chk(f"{c} r vs DK", round(r, 4), f"> {lo}", r > lo)
        # The slope PCs are atlas-specific decompositions and are NOT the same
        # phenotype across atlases (sign is arbitrary); report, do not gate.
        for c in ("slope_PC1", "slope_PC2", "slope_PC3"):
            r = float(np.corrcoef(m[f"{c}_hcp"], m[f"{c}_dk"])[0, 1])
            chk(f"{c} |r| vs DK (info)", round(abs(r), 4), "not gated", True)
    width = max(len(c[0]) for c in checks); failed = 0
    for name, got, want, ok in checks:
        failed += (not ok)
        print(f"{'PASS' if ok else 'FAIL'}  {name:<{width}}  got={got}  want={want}")
    print()
    if failed:
        print(f"{failed} of {len(checks)} HCP checks FAILED -- inspect before using this export."); return 1
    print(f"all {len(checks)} HCP checks passed -- the HCP-MMP export is consistent with the DK run."); return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parc", default="dsk", choices=["dsk", "hcp"],
                    help="dsk: gate against docs/vintage_comparison.csv; hcp: gate against the DK run")
    ap.add_argument("--run-dir", default=None, help="override the run directory")
    a = ap.parse_args()
    if a.parc == "hcp":
        run_dir = Path(a.run_dir) if a.run_dir else None
        if run_dir is None:
            import subprocess, os
            env = dict(os.environ, ABCD_CONFIG="ct_70_hcp_noglobal_mv2", PYTHONPATH=str(REPO / "src"))
            run_dir = Path(subprocess.check_output([sys.executable, "-m", "abcd.run_dir"], env=env, text=True).strip())
        return main_hcp(run_dir)
    if not NEW.exists():
        print(f"FAIL  export run dir missing: {NEW}")
        return 2
    exp = expected()
    man = json.loads((NEW / "manifest.json").read_text())
    fin = man["final"]
    visit_counts = fin.get("visit_counts", {})          # scans per visit
    n_by_visits = {int(k): int(vv) for k, vv
                   in (fin.get("n_visits_distribution") or {}).items()}

    checks: list[tuple[str, object, object, bool]] = []

    def chk(name, got, want, ok):
        checks.append((name, got, want, bool(ok)))

    chk("release", man["release"], "7.0", man["release"] == "7.0")
    chk("subjects_ge2_visits", fin["n_subjects"], int(exp["subjects_ge2_visits"]),
        fin["n_subjects"] == int(exp["subjects_ge2_visits"]))
    chk("scans_in_model", fin["n_scans"], int(exp["scans_in_model"]),
        fin["n_scans"] == int(exp["scans_in_model"]))
    chk("families", fin["n_families"], int(exp["families"]),
        fin["n_families"] == int(exp["families"]))
    chk("age_centre_years", round(man["age_centre"], 4), exp["age_centre_years"],
        abs(man["age_centre"] - exp["age_centre_years"]) < 5e-4)
    for k in (2, 3, 4):
        key = f"subjects_{k}_visits"
        if k in n_by_visits and key in exp:
            chk(key, n_by_visits[k], int(exp[key]), n_by_visits[k] == int(exp[key]))
    for v in ("v0", "v2", "v4", "v6"):
        key = f"scans_{v}"
        if v in visit_counts and key in exp:
            chk(key, visit_counts[v], int(exp[key]), visit_counts[v] == int(exp[key]))

    # --- the GCTA export itself -------------------------------------------
    gi = NEW / "gcta_inputs" / "phenotypes_gcta.txt"
    if not gi.exists():
        chk("gcta_inputs/phenotypes_gcta.txt", "MISSING", "present", False)
    else:
        ph = pd.read_csv(gi, sep=r"\s+")
        cols = [c for c in ph.columns if c not in ("FID", "IID")]
        chk("export phenotypes", cols, EXPECTED_PHENOTYPES, cols == EXPECTED_PHENOTYPES)
        chk("export rows", len(ph), int(exp["subjects_ge2_visits"]),
            len(ph) == int(exp["subjects_ge2_visits"]))
        chk("export IID spelling", str(ph.IID.iloc[0]), "NDAR_INV...",
            str(ph.IID.iloc[0]).startswith("NDAR_INV"))

        # --- against the archived 6.0-vintage export ------------------------
        old_gi = OLD / "gcta_inputs_v3" / "phenotypes_gcta.txt"
        if old_gi.exists():
            old = pd.read_csv(old_gi, sep=r"\s+")
            m = ph.merge(old, on="IID", suffixes=("_new", "_old"))
            chk("subjects_shared", len(m), int(exp["subjects_shared"]),
                len(m) == int(exp["subjects_shared"]))
            r = float(np.corrcoef(m.global_slope_new, m.global_slope_old)[0, 1])
            chk("global_slope r vs 6.0", round(r, 4),
                exp["global_slope_pearson_shared_subjects"],
                abs(r - exp["global_slope_pearson_shared_subjects"]) < 0.01)
            rb = float(np.corrcoef(m.baseline_thickness_new,
                                   m.baseline_thickness_old)[0, 1])
            chk("baseline_thickness r vs 6.0", round(rb, 4), "> global_slope r", rb > r)
        else:
            chk("archived 6.0 export", "MISSING", "present", False)

    width = max(len(c[0]) for c in checks)
    failed = 0
    for name, got, want, ok in checks:
        if not ok:
            failed += 1
        print(f"{'PASS' if ok else 'FAIL'}  {name:<{width}}  got={got}  want={want}")
    print()
    if failed:
        print(f"{failed} of {len(checks)} checks FAILED -- do not use this export "
              f"downstream; rsync the laptop's instead (README_HPC.md Step 0).")
        return 1
    print(f"all {len(checks)} checks passed -- export matches the committed local "
          f"re-run and may be used for steps 2-8.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
