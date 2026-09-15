import pandas as pd, sys
run, ph, out = sys.argv[1:4]
b = pd.read_parquet(f"{run}/fits/blups.parquet"); assert b.label.nunique() == 1
b["_t"] = b.subject.str.extract(r"([A-Z0-9]{8})$")[0]
ref = pd.read_csv(f"{ph}/phenotypes_gcta.txt", sep=r"\s+", dtype={"FID": str, "IID": str})
ref["_t"] = ref.IID.str.extract(r"([A-Z0-9]{8})$")[0]
m = ref[["FID", "IID", "_t", "global_slope", "baseline_thickness"]].merge(
    b[["_t", "re_slope", "re_intercept"]], on="_t", how="inner")
assert len(m) == len(ref), (len(m), len(ref))
for c in ("re_slope", "re_intercept"): m[c] = (m[c] - m[c].mean()) / m[c].std()
print(f"corr(global_slope per-region-avg, single-LMM slope) = {m.global_slope.corr(m.re_slope):.4f}; "
      f"baseline: {m.baseline_thickness.corr(m.re_intercept):.4f}; n = {len(m)}")
out_df = m[["FID", "IID"]].assign(global_slope_1lmm=m.re_slope, baseline_thickness_1lmm=m.re_intercept)
out_df.to_csv(f"{out}/phenotypes_gcta.txt", sep=" ", index=False)
pd.DataFrame([dict(name="global_slope_1lmm", mpheno=1, priority=1, role="primary (single LMM)", n_nonmissing=len(m)),
              dict(name="baseline_thickness_1lmm", mpheno=2, priority=0, role="positive control (single LMM)", n_nonmissing=len(m))]
            ).to_csv(f"{out}/phenotype_manifest.tsv", sep="\t", index=False)
