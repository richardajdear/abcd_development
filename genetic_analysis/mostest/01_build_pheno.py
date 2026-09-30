#!/usr/bin/env python
"""Step 16.1: regional phenotype files for the REGENIE + MOSTest discovery arm.

Writes, under --out (default genetic_analysis/work/results_70tab/mostest/pheno,
gitignored, PER-SUBJECT -- rule 16, never commit):

  pheno_<family>.txt       FID IID <k measures>   real data
  pheno_<family>_perm.txt  same columns, rows permuted (see below)
  covar.txt / covar_perm.txt   FID IID sex site baseline_age n_visits PC1..PC10
  keep_imp.txt             FID IID of the analysis set, in the imputed .fam's spelling
  keep_array.txt           the same children in the array .fam's spelling
  array_idmap.txt          old FID, old IID, new FID, new IID  (plink --update-ids)
  build_summary.json       counts and checks (summary level; copied to the
                           tracked table by 06_collect.py)

Families (DK, 68 regions each):
  slope     per-region LMM random slope BLUP   (re_slope, blups.parquet)
  ct        per-region LMM random intercept BLUP = regional thickness at the
            age centre (re_intercept) -- the POSITIVE CONTROL: cross-sectional,
            reliable, LDSC h2 z ~4.4 for its whole-cortex mean
  slopeols  per-child OLS slope of thickness on age, per region -- unshrunk
            sensitivity for the region-dependent BLUP shrinkage caveat
  global    global_slope_1lmm, baseline_thickness_1lmm (engine check against the
            existing GENESIS scan; not a MOSTest family)

Permutation.  MOSTest permutes genotypes across individuals once and re-runs
every univariate GWAS.  Shuffling the phenotype AND covariate rows jointly
against the genotype IDs is the same operation.  One permutation (seed 16) is
shared by all families.  The real-data and permuted files have identical ID
columns in identical order, so REGENIE sees the same sample either way.

IDs (rules 1-2).  Children are joined on the 8-character NDAR token.  REGENIE
needs FID/IID exactly as in the genotype .fam; the imputed PRS fileset and the
array fileset are matched separately, and the array IDs are rewritten to the
imputed spelling so step 1 and step 2 agree.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

TOK = re.compile(r"([A-Z0-9]{8})$")


def token(s: pd.Series) -> pd.Series:
    return s.astype(str).str.extract(TOK)[0]


def read_ws(path: Path, **kw) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", dtype={"FID": str, "IID": str}, **kw)


def read_fam(prefix: str) -> pd.DataFrame:
    f = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, usecols=[0, 1],
                    names=["FID", "IID"], dtype=str)
    f["_t"] = token(f.IID)
    if f._t.isna().any():
        raise SystemExit(f"FATAL: {prefix}.fam has IIDs without an 8-char token")
    if f._t.duplicated().any():
        raise SystemExit(f"FATAL: {prefix}.fam has duplicated tokens")
    return f


def wide(b: pd.DataFrame, col: str, prefix: str) -> pd.DataFrame:
    w = b.pivot(index="_t", columns="label", values=col)
    w.columns = [f"{prefix}_{c}" for c in w.columns]
    return w


def ols_slopes(mt: pd.DataFrame) -> pd.DataFrame:
    """Per child x region OLS slope of thickness on age (mm / year)."""
    mt = mt[["_t", "label", "age", "value"]].dropna()
    g = mt.groupby(["_t", "label"])
    n = g.age.transform("size")
    am = g.age.transform("mean"); vm = g.value.transform("mean")
    mt = mt.assign(xy=(mt.age - am) * (mt.value - vm), xx=(mt.age - am) ** 2, n=n)
    s = mt.groupby(["_t", "label"])[["xy", "xx"]].sum()
    s["slope"] = s.xy / s.xx
    w = s.slope.unstack("label")
    w.columns = [f"slopeols_{c}" for c in w.columns]
    return w


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="out/thickness_dsk_70_139406217085",
                    help="DK per-region LMM run dir (fits/blups.parquet, model_table.parquet)")
    ap.add_argument("--ref", default="genetic_analysis/work/results_70tab/prs_final_1lmm/pheno",
                    help="single-LMM export of the 8,596-child analysis set "
                         "(phenotypes_gcta.txt with FID = family id, covar_quant.txt, covar_categorical.txt)")
    ap.add_argument("--geno-imp", default="genetic_analysis/work/inputs/geno/abcd_imp_prs")
    ap.add_argument("--geno-array", required=True, help="array PLINK prefix (GENO_ARRAY)")
    ap.add_argument("--out", default="genetic_analysis/work/results_70tab/mostest/pheno")
    ap.add_argument("--seed", type=int, default=16)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    ref = Path(a.ref)

    # ---- analysis set and covariates --------------------------------------
    ph = read_ws(ref / "phenotypes_gcta.txt"); ph["_t"] = token(ph.IID)
    cq = read_ws(ref / "covar_quant.txt"); cq["_t"] = token(cq.IID)
    cc = read_ws(ref / "covar_categorical.txt"); cc["_t"] = token(cc.IID)
    need_g = {"global_slope_1lmm", "baseline_thickness_1lmm"}
    if not need_g <= set(ph.columns):
        raise SystemExit(f"FATAL: {ref}/phenotypes_gcta.txt lacks {need_g - set(ph.columns)}")
    base = ph[["_t", *sorted(need_g)]].merge(cq.drop(columns=["FID", "IID"]), on="_t") \
                                      .merge(cc.drop(columns=["FID", "IID"]), on="_t")
    if len(base) != len(ph):
        raise SystemExit(f"FATAL: covariates cover {len(base)} of {len(ph)} children")
    base["sex"] = base.sex.map({"F": 0, "M": 1}).astype("Int64")
    base["site"] = pd.factorize(base.site.astype(str), sort=True)[0]
    if base.sex.isna().any():
        raise SystemExit("FATAL: sex not coded F/M")

    # ---- regional phenotypes ------------------------------------------------
    run = Path(a.run)
    b = pd.read_parquet(run / "fits" / "blups.parquet")
    b["_t"] = token(b.subject)
    labels = sorted(b.label.unique())
    if len(labels) != 68:
        raise SystemExit(f"FATAL: expected 68 DK labels, found {len(labels)}")
    fam = {"slope": wide(b, "re_slope", "slope"), "ct": wide(b, "re_intercept", "ct")}
    mt = pd.read_parquet(run / "model_table.parquet", columns=["subject", "label", "age", "value"])
    mt["_t"] = token(mt.subject)
    fam["slopeols"] = ols_slopes(mt)
    fam["global"] = base.set_index("_t")[sorted(need_g)]

    # ---- genotype ID spellings --------------------------------------------
    imp = read_fam(a.geno_imp); arr = read_fam(a.geno_array)
    ids = base[["_t"]].merge(imp.rename(columns={"FID": "FID", "IID": "IID"}), on="_t", how="inner")
    ids = ids.merge(arr.rename(columns={"FID": "FID_arr", "IID": "IID_arr"}), on="_t", how="inner")
    for w in fam.values():
        ids = ids[ids._t.isin(w.dropna().index)]
    ids = ids.sort_values("_t").reset_index(drop=True)
    n = len(ids)
    if n < 8500:
        raise SystemExit(f"FATAL: only {n} children have all phenotypes and both genotype sets "
                         f"(analysis set {len(base)}, imputed {base._t.isin(imp._t).sum()}, "
                         f"array {base._t.isin(arr._t).sum()})")

    ids[["FID", "IID"]].to_csv(out / "keep_imp.txt", sep=" ", index=False, header=False)
    ids[["FID_arr", "IID_arr"]].to_csv(out / "keep_array.txt", sep=" ", index=False, header=False)
    ids[["FID_arr", "IID_arr", "FID", "IID"]].to_csv(out / "array_idmap.txt", sep=" ",
                                                     index=False, header=False)

    rng = np.random.default_rng(a.seed)
    perm = rng.permutation(n)
    covcols = ["sex", "site", "baseline_age", "n_visits", *[f"PC{i}" for i in range(1, 11)]]
    cov = base.set_index("_t").loc[ids._t, covcols].reset_index(drop=True)
    head = ids[["FID", "IID"]]
    pd.concat([head, cov], axis=1).to_csv(out / "covar.txt", sep=" ", index=False)
    pd.concat([head, cov.iloc[perm].reset_index(drop=True)], axis=1) \
      .to_csv(out / "covar_perm.txt", sep=" ", index=False)

    summ = dict(n=n, n_analysis_set=len(base), seed=a.seed, run=str(run), families={})
    for name, w in fam.items():
        v = w.loc[ids._t].reset_index(drop=True)
        if v.isna().any().any():
            raise SystemExit(f"FATAL: {name} has missing values after subsetting")
        pd.concat([head, v], axis=1).to_csv(out / f"pheno_{name}.txt", sep=" ", index=False,
                                            float_format="%.6g")
        pd.concat([head, v.iloc[perm].reset_index(drop=True)], axis=1) \
          .to_csv(out / f"pheno_{name}_perm.txt", sep=" ", index=False, float_format="%.6g")
        (out / f"traits_{name}.txt").write_text("\n".join(v.columns) + "\n")
        R = np.corrcoef(v.rank().values, rowvar=False) if v.shape[1] > 1 else np.ones((1, 1))
        iu = np.triu_indices(v.shape[1], 1)
        summ["families"][name] = dict(k=int(v.shape[1]),
                                      median_pheno_r=float(np.median(R[iu])) if iu[0].size else None)
    # consistency: per-region BLUP-mean vs single-LMM whole-cortex trait
    g = fam["global"].loc[ids._t]
    summ["r_mean_slope_blup_vs_global_slope_1lmm"] = float(
        np.corrcoef(fam["slope"].loc[ids._t].mean(axis=1), g.global_slope_1lmm)[0, 1])
    summ["r_mean_ct_blup_vs_baseline_thickness_1lmm"] = float(
        np.corrcoef(fam["ct"].loc[ids._t].mean(axis=1), g.baseline_thickness_1lmm)[0, 1])
    (out / "build_summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ, indent=2))
    if summ["r_mean_slope_blup_vs_global_slope_1lmm"] < 0.7:
        raise SystemExit("FATAL: regional slope BLUPs do not track the whole-cortex slope "
                         "(expect r ~0.8, README_HPC sec 1) -- wrong run dir or ID join")


if __name__ == "__main__":
    main()
