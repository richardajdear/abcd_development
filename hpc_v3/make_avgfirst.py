"""Does averaging regions BEFORE the LMM beat averaging the per-region BLUPs?

Motivation (user question, 2026-09-09): the four Experiment A phenotypes are
means of per-region slope BLUPs; if averaging over FEWER model fits (16
regions vs 68) adds noise, that could depress their heritability, and
averaging the thickness series first — one LMM on the region-mean CT —
might recover it.

Run as:  python hpc_v3/make_avgfirst.py          (repo root; needs R + lme4)

Pipeline: for each of {global, topDelta, topC3} x {both hemis, lh, rh}, build
the region-mean thickness per subject-visit from the settled run's
model_table.parquet, fit
    lmer(value ~ age_c + sex + (1 + age_c | subject) + (1 | site))
(the pipeline's own specification), and take fixef + BLUP as the avg-first
subject slope.  Committed outputs (summary-level only — per-subject scores
stay under hpc_v3/work/, which is gitignored):

  avgfirst_comparison.csv   r between constructions, lh/rh Spearman-Brown
                            consistency for both, r-with-global for both
  avgfirst_density.csv      60x60 binned counts of avg-first vs BLUP-mean
                            (what figure_avgfirst.py draws as the scatter)

ANSWER (n = 8,192): averaging first does NOT help.  r(methods) = 0.887 /
0.936 / 0.924 (global / topDelta / topC3) and the avg-first construction is
slightly LESS internally consistent (SB 0.860 vs 0.872, 0.796 vs 0.816,
0.773 vs 0.787).  The balanced design forces the fixed-effect parts to be
algebraically identical (see mem: every subject-visit has all 68 regions), so
the two constructions differ only in where BLUP shrinkage is applied — and
empirically that is a wash, with the settled mean-of-BLUPs marginally ahead.
The subset phenotypes' low h2 is therefore not a construction artefact.
Global avg-first reproduces the 2026-09-08 session's r = 0.887 exactly.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
WORK = HERE / "work"
RUN = REPO / "out/thickness_dsk_70_139406217085"


def export_series() -> None:
    m = pd.read_parquet(RUN / "model_table.parquet",
                        columns=["subject", "visit", "label", "value",
                                 "age_c", "sex", "site"])
    r = pd.read_csv(HERE / "phenotypes_v3_regions.csv")
    sets = {
        "global": [f"{h}_{b}" for b in r.base for h in ("lh", "rh")],
        "topDelta": [f"{h}_{b}" for b in r[r.in_topDelta].base
                     for h in ("lh", "rh")],
        "topC3": [f"{h}_{b}" for b in r[r.in_topC3].base for h in ("lh", "rh")],
    }
    for nm in list(sets):
        sets[f"{nm}_lh"] = [x for x in sets[nm] if x.startswith("lh_")]
        sets[f"{nm}_rh"] = [x for x in sets[nm] if x.startswith("rh_")]
    out = []
    for nm, labs in sets.items():
        s = (m[m.label.isin(labs)]
             .groupby(["subject", "visit", "age_c", "sex", "site"],
                      observed=True)["value"]
             .agg(["mean", "count"]).reset_index())
        assert s["count"].eq(len(labs)).all(), (nm, s["count"].unique())
        out.append(s.drop(columns="count")
                   .rename(columns={"mean": "value"}).assign(series=nm))
    WORK.mkdir(exist_ok=True)
    pd.concat(out).to_csv(WORK / "avgfirst_series.csv", index=False)


R_FIT = r"""
suppressMessages({library(data.table); library(lme4)})
d <- fread("hpc_v3/work/avgfirst_series.csv")
res <- list()
for (nm in unique(d$series)) {
  s <- d[series == nm]
  m <- lmer(value ~ age_c + sex + (1 + age_c | subject) + (1 | site),
            data = s, REML = TRUE,
            control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  re <- ranef(m)$subject
  res[[nm]] <- data.table(subject = rownames(re), series = nm,
                          slope = fixef(m)[["age_c"]] + re[["age_c"]])
  cat(nm, "fixef age_c =", round(fixef(m)[["age_c"]], 5), "\n")
}
fwrite(rbindlist(res), "hpc_v3/work/avgfirst_slopes.csv")
"""


def compare() -> None:
    alt = (pd.read_csv(WORK / "avgfirst_slopes.csv")
           .pivot(index="subject", columns="series", values="slope"))
    ph = pd.read_parquet(RUN / "phenotypes/phenotypes.parquet")
    sl = ph[ph.phenotype == "slope"].pivot(index="subject", columns="label",
                                           values="value")
    r = pd.read_csv(HERE / "phenotypes_v3_regions.csv")
    sets = {"global": list(sl.columns),
            "topDelta": [f"{h}_{b}" for b in r[r.in_topDelta].base
                         for h in ("lh", "rh")],
            "topC3": [f"{h}_{b}" for b in r[r.in_topC3].base
                      for h in ("lh", "rh")]}
    blup = pd.DataFrame({nm: sl[labs].mean(axis=1)
                         for nm, labs in sets.items()})
    alt = alt.reindex(blup.index)
    assert alt.notna().all().all()

    def sb(x, y):
        rr = x.corr(y)
        return 2 * rr / (1 + rr)

    rows, dens = [], []
    for nm in sets:
        b_lh = sl[[x for x in sets[nm] if x.startswith("lh_")]].mean(axis=1)
        b_rh = sl[[x for x in sets[nm] if x.startswith("rh_")]].mean(axis=1)
        rows.append(dict(
            phenotype=nm,
            r_methods=blup[nm].corr(alt[nm]),
            rho_methods=blup[nm].corr(alt[nm], method="spearman"),
            sb_blupmean=sb(b_lh, b_rh),
            sb_avgfirst=sb(alt[f"{nm}_lh"], alt[f"{nm}_rh"]),
            r_avgfirst_vs_global_avgfirst=alt[nm].corr(alt["global"]),
            r_blupmean_vs_global_blupmean=blup[nm].corr(blup["global"])))
        H, xe, ye = np.histogram2d(blup[nm], alt[nm], bins=60)
        nz = np.nonzero(H)
        dens.append(pd.DataFrame(dict(
            phenotype=nm, x_lo=xe[nz[0]], x_hi=xe[nz[0] + 1],
            y_lo=ye[nz[1]], y_hi=ye[nz[1] + 1], count=H[nz].astype(int))))
    pd.DataFrame(rows).round(4).to_csv(HERE / "avgfirst_comparison.csv",
                                       index=False)
    pd.concat(dens).to_csv(HERE / "avgfirst_density.csv", index=False)
    print(pd.DataFrame(rows).round(4).to_string(index=False))


def main() -> int:
    export_series()
    (WORK / "fit_avgfirst.R").write_text(R_FIT)
    subprocess.run(["Rscript", str(WORK / "fit_avgfirst.R")],
                   cwd=REPO, check=True)
    compare()
    return 0


if __name__ == "__main__":
    sys.exit(main())
