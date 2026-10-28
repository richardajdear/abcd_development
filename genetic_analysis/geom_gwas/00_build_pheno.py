"""D1 step 9: GENESIS/GCTA phenotype files for the geometry-slope GWAS.

Modelled on orderops/build_1lmm_pheno.py.  Reads the laptop export
(directions/d1_multimodal/work/gwas_pheno_geom.tsv, raw BLUPs), matches it on
the 8-char NDAR token to the HCP single-LMM reference set (FID/IID and row
order of prs_final_1lmm/pheno), z-scores each column over the matched set and
writes pheno/ (FID = family_id, for GENESIS) and pheno_fidiid/ (FID = IID, for
GCTA against the .grm.id; rule 2).  Covariate files are symlinked unchanged.

  python 00_build_pheno.py <geom.tsv> <prs_final_1lmm dir> <scan_geom dir>
"""
import os, sys
import pandas as pd

geom_f, ref_root, W = sys.argv[1:4]
TAGS = ["thk", "area", "vol", "thk_area", "thk_vol_area"]
COLS = [f"slope_{t}" for t in TAGS] + [f"int_{t}" for t in TAGS]

g = pd.read_csv(geom_f, sep="\t")
g["_t"] = g.subject.str.extract(r"([A-Z0-9]{8})$")[0]
assert g._t.notna().all() and g._t.is_unique

# Laptop child-level r between the slopes (brief), checked on the raw export.
LAPTOP = {("area", "thk"): .18, ("vol", "thk"): .73, ("vol", "area"): .65,
          ("thk_area", "thk"): .93, ("thk_area", "area"): .51, ("thk_area", "vol"): .82,
          ("thk_vol_area", "thk"): .90, ("thk_vol_area", "area"): .58,
          ("thk_vol_area", "vol"): .92, ("thk_vol_area", "thk_area"): .98}
r_raw = g[[f"slope_{t}" for t in TAGS]].corr()
bad = []
for (a, b), want in LAPTOP.items():
    got = r_raw.loc[f"slope_{a}", f"slope_{b}"]
    print(f"r(slope_{a}, slope_{b}) = {got:.3f}  (laptop {want:.2f})")
    if round(got, 2) != want: bad.append((a, b, got, want))

out = {}
for arm, sub in (("pheno", "pheno"), ("pheno_fidiid", "pheno_fidiid")):
    ref = pd.read_csv(f"{ref_root}/{sub}/phenotypes_gcta.txt", sep=r"\s+", dtype={"FID": str, "IID": str})
    ref["_t"] = ref.IID.str.extract(r"([A-Z0-9]{8})$")[0]
    m = ref.merge(g[["_t"] + COLS], on="_t", how="left")
    assert len(m) == len(ref) and (m.IID.values == ref.IID.values).all()
    for c in COLS: m[c] = (m[c] - m[c].mean()) / m[c].std()
    d = f"{W}/{arm}"; os.makedirs(d, exist_ok=True)
    m[["FID", "IID"] + COLS].to_csv(f"{d}/phenotypes_gcta.txt", sep=" ", index=False, na_rep="NA")
    pd.DataFrame([dict(name=c, mpheno=i + 1, priority=int(c.startswith("slope_")),
                       role=("reproduction control" if c == "slope_thk" else
                             "decomposition" if c in ("slope_area", "slope_vol") else
                             "secondary composite" if c.startswith("slope_") else "positive control (intercept)"),
                       n_nonmissing=int(m[c].notna().sum())) for i, c in enumerate(COLS)]
                 ).to_csv(f"{d}/phenotype_manifest.tsv", sep="\t", index=False)
    for f in ("covar_quant.txt", "covar_categorical.txt"):
        dst = f"{d}/{f}"
        if os.path.lexists(dst): os.remove(dst)
        os.symlink(os.path.realpath(f"{ref_root}/{sub}/{f}"), dst)
    out[arm] = m

m = out["pheno"]
print(f"\nreference set n = {len(m)}; matched = {m.slope_thk.notna().sum()}; "
      f"export rows = {len(g)}, export not in reference = {(~g._t.isin(m._t)).sum()}")
print("n missing per column:", {c: int(m[c].isna().sum()) for c in COLS})
r_s = m.slope_thk.corr(m.global_slope_1lmm); r_i = m.int_thk.corr(m.baseline_thickness_1lmm)
print(f"r(slope_thk, global_slope_1lmm) = {r_s:.4f}\nr(int_thk, baseline_thickness_1lmm) = {r_i:.4f}")
print("matched-set r between slopes:\n" + m[[f"slope_{t}" for t in TAGS]].corr().round(3).to_string())
if bad: print("WARNING: laptop r not reproduced to 2 dp:", bad)
print(f"GATE A: {'PASS' if r_s >= 0.99 else 'FAIL'} (r = {r_s:.4f}, threshold 0.99)")
sys.exit(0 if r_s >= 0.99 else 1)
