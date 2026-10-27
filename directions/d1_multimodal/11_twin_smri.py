"""D1 step 7c: twin h2 of each new structural slope and its genetic correlation with the thinning slope, under the
Figure-1 model (M0) and the scanner-epoch model (M2) of 10_scanner_terms.R.

    python directions/d1_multimodal/11_twin_smri.py [n_boot]     (env abcd-spatial; default 300)

Slopes are residualised on male, age_first, age_span, n_visits + site, as D3 (03_build_traits_pairs.py). Thinning is
taken from the same model as the measure (M0 with M0, M2 with M2), oriented so that higher = faster thinning.
AE models (D3 ace.py); pair bootstrap CIs. Output: tables/twin_smri.tsv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
D3 = HERE.parent / "d3_twin_family"
sys.path.insert(0, str(D3))
import ace  # noqa: E402
from common import residualise  # noqa: E402

NB = int(sys.argv[1]) if len(sys.argv) > 1 else 300
B = pd.read_csv(HERE / "work" / "blups_scanner_terms.csv.gz")
T = pd.read_parquet(D3 / "work" / "traits.parquet")
pairs = pd.read_csv(D3 / "work" / "pairs.csv")
BR = ["male", "age_first", "age_span", "n_visits"]
MET = ["smri__area_hcp", "smri__vol_hcp", "smri__gwc_hcp", "smri__t1gm50_hcp", "smri__t1wm1_hcp", "smri__t1t2_dsk",
       "dti__is__md__gwc", "dti__is__md__gm"]
q = lambda s: np.nanpercentile(s, [2.5, 97.5])
rows = []
for mod in ("M0", "M2"):
    th = -B[(B.metric == "smri__thk_hcp") & (B.model == mod)].set_index("subject").slope
    T["thin"] = th.reindex(T.index); T["thin_r"] = residualise(T, "thin", BR, ["site"])
    for mt in MET:
        T["x"] = B[(B.metric == mt) & (B.model == mod)].set_index("subject").slope.reindex(T.index)
        T["x_r"] = residualise(T, "x", BR, ["site"])
        mz = ace.pair_arrays(T, pairs, ["x_r"], "MZ"); dz = ace.pair_arrays(T, pairs, ["x_r"], "DZ")
        u = ace.fit(mz, dz, "AE"); ub = ace.bootstrap(mz, dz, "AE", NB, seed=31, start=u["theta"], keys=["h2_1"])
        mz2 = ace.pair_arrays(T, pairs, ["thin_r", "x_r"], "MZ"); dz2 = ace.pair_arrays(T, pairs, ["thin_r", "x_r"], "DZ")
        f = ace.fit(mz2, dz2, "AE")
        bs = ace.bootstrap(mz2, dz2, "AE", NB, seed=32, start=f["theta"], keys=["h2_1", "h2_2", "rA", "rE", "rP"])
        r = dict(metric=mt, model=mod, h2=u["h2_1"], h2_lo=q(ub.h2_1)[0], h2_hi=q(ub.h2_1)[1], h2_thin=f["h2_1"],
                 rA=f["rA"], rA_lo=q(bs.rA)[0], rA_hi=q(bs.rA)[1], rE=f["rE"], rP=f["rP"], n_mz=f["n_mz"], n_dz=f["n_dz"])
        rows.append(r)
        print(f"{mod} {mt:18s} h2 {r['h2']:.2f} [{r['h2_lo']:.2f}, {r['h2_hi']:.2f}]  thin h2 {r['h2_thin']:.2f}  "
              f"rA {r['rA']:+.2f} [{r['rA_lo']:+.2f}, {r['rA_hi']:+.2f}]  rE {r['rE']:+.2f}  rP {r['rP']:+.2f}", flush=True)
pd.DataFrame(rows).to_csv(HERE / "tables" / "twin_smri.tsv", sep="\t", index=False)
