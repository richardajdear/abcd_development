"""Align the 7.0 GCTA export to the cluster's ID conventions.

Adapted 2026-09-14 from the hpc_v3 version (README_HPC.md 2.1).  The original
aligned a v3 export against v2's published $PHENO and **copied v2's covariate
files**, which was right only while the ancestry PCs were identical between
vintages.  The 7.0 tables recomputed ab_g_stc__gen_pc__01..32 for all 11,670
genotyped children -- every value differs from 6.0 -- so the export's own
covariates are now the correct ones and the copy is dropped.  What is kept from
the original is its three-part logic, each part fixing a silent failure:

1. **IID spelling.**  The export writes ``NDAR_INV005V6D2C``; every ``.fam``,
   ``.grm.id`` and ``.profile`` on CSD3 spells that child ``sub-005V6D2C``.
   GENESIS normalises to the 8-char token and would not care; GCTA and PLINK
   match FID+IID literally and return zero subjects without an error
   (README_HPC.md 4 rule 1).

2. **Sample.**  Subset to the phenotyped-and-genotyped children by joining the
   8-char token against the array ``.fam`` -- 8,596 of the 8,716 phenotyped,
   up from 8,082 on the 6.0 tables.

3. **Scale.**  z-score every phenotype column *after* subsetting, so a beta is
   per SD of the analysis sample (rule 3).  ``global_slope`` has native
   variance ~1e-6, so the native mean/SD are written to ``align_report.tsv``:
   multiply a standardised beta by the SD there to get mm/yr.

Two FID conventions are written side by side (rule 2): ``<out>`` has
FID = family id, for ``prs_assoc.R`` and the Fulker step; ``<out>_fidiid`` has
FID = IID, matching the genotype ``.fam``, for GCTA.

    python genetic_analysis/setup/align_export.py <src-gcta_inputs> <out-dir>
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
FAM = ("/rds/project/rds-CeXlNYOYMxw/Data_Genetics/genotype_microarray/"
       "smokescreen/merged_chroms.fam")
EUR = REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep"


def token(s: pd.Series, what: str) -> pd.Series:
    """The 8-char NDAR token every ABCD id spelling ends with."""
    t = s.astype(str).str.extract(r"([A-Z0-9]{8})$")[0]
    if t.isna().any():
        raise SystemExit(f"{what}: {int(t.isna().sum())} ids do not end in an "
                         f"8-char token, e.g. {s[t.isna()].iloc[0]!r}")
    return t


def read_space(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep=r"\s+")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("src", help="the run's gcta_inputs/ directory")
    ap.add_argument("out", help="aligned output directory (pheno_70tab)")
    ap.add_argument("--fam", default=FAM, help="array .fam giving the cluster spelling")
    ap.add_argument("--eur", default=str(EUR), help="EUR-arm keep list, for the count")
    a = ap.parse_args(argv)
    src, out = Path(a.src), Path(a.out)
    out_fid = out.parent / f"{out.name}_fidiid"
    for d in (out, out_fid):
        d.mkdir(parents=True, exist_ok=True)

    phen = read_space(src / "phenotypes_gcta.txt")
    cols = [c for c in phen.columns if c not in ("FID", "IID")]
    n_export = len(phen)

    # 1 + 2: the .fam spelling, and the phenotyped-and-genotyped subset.
    fam = pd.read_csv(a.fam, sep=r"\s+", header=None,
                      usecols=[0, 1], names=["fam_FID", "fam_IID"], dtype=str)
    fam["_t"] = token(fam.fam_IID, "genotype .fam")
    if fam._t.duplicated().any():
        raise SystemExit("the genotype .fam has duplicate 8-char tokens")
    phen["_t"] = token(phen.IID, "export")

    keep = phen.merge(fam[["_t", "fam_IID"]], on="_t", how="inner")
    keep = keep.sort_values("fam_IID", kind="stable").reset_index(drop=True)
    family_id = keep.FID.astype(str)
    keep["IID"] = keep.fam_IID                      # sub-xxxxxxxx
    keep["FID"] = family_id                         # the real family id

    # 3: z-score over the analysis sample, recording the native units.
    rows = []
    for c in cols:
        m, s = keep[c].mean(), keep[c].std()
        if not s > 0:
            raise SystemExit(f"{c} has zero variance over the analysis sample")
        rows.append({"phenotype": c, "n": int(keep[c].notna().sum()),
                     "native_mean": m, "native_sd": s})
        keep[c] = (keep[c] - m) / s

    # covariates: the export's own (7.0 PCs), respelled and subset identically.
    cov = {}
    for f in ("covar_quant.txt", "covar_categorical.txt"):
        c = read_space(src / f)
        c["_t"] = token(c.IID, f)
        c = keep[["FID", "IID", "_t"]].merge(
            c.drop(columns=["FID", "IID"]), on="_t", how="left")
        missing = c.drop(columns=["FID", "IID", "_t"]).isna().all(axis=1).sum()
        if missing:
            raise SystemExit(f"{f}: {missing} analysis subjects have no covariates")
        cov[f] = c.drop(columns="_t")

    keep = keep[["FID", "IID"] + cols]

    # both FID conventions (rule 2)
    for d, fid in ((out, keep.FID), (out_fid, keep.IID)):
        k = keep.copy(); k["FID"] = fid.values
        k.to_csv(d / "phenotypes_gcta.txt", sep=" ", index=False, na_rep="NA")
        for f, frame in cov.items():
            c = frame.copy(); c["FID"] = fid.values
            c.to_csv(d / f, sep=" ", index=False, na_rep="NA")
        man = pd.read_csv(src / "phenotype_manifest.tsv", sep="\t")
        if "n_nonmissing" in man.columns:
            man["n_nonmissing"] = man.name.map(
                {r["phenotype"]: r["n"] for r in rows}).fillna(man.n_nonmissing)
        man.to_csv(d / "phenotype_manifest.tsv", sep="\t", index=False)

    # family map for setup/make_grm_famid.sh: IID (.fam spelling) -> family id.
    # It must come from THIS export -- gn_y_genrel's family/birth ids were
    # re-coded in 7.0, so v1's pheno_allanc/family_map.tsv is a wrong input.
    fmap = pd.DataFrame({"IID": keep.IID.values, "family_id": family_id.values})
    fmap.to_csv(out / "family_map.tsv", sep="\t", index=False)
    shutil.copy2(out / "family_map.tsv", out_fid / "family_map.tsv")

    # the record README_HPC.md Step 1 asks for
    rep = pd.DataFrame(rows)
    fam_sizes = family_id.value_counts()
    eur = pd.read_csv(a.eur, sep=r"\s+", header=None, usecols=[1],
                      names=["IID"], dtype=str)
    n_eur = keep.IID.isin(set(eur.IID)).sum()
    meta = pd.DataFrame([
        {"phenotype": "_n_phenotyped", "n": n_export, "native_mean": "", "native_sd": ""},
        {"phenotype": "_n_analysis", "n": len(keep), "native_mean": "", "native_sd": ""},
        {"phenotype": "_n_families", "n": int(fam_sizes.size), "native_mean": "", "native_sd": ""},
        {"phenotype": "_n_families_multi", "n": int((fam_sizes > 1).sum()), "native_mean": "", "native_sd": ""},
        {"phenotype": "_n_eur_arm", "n": int(n_eur), "native_mean": "", "native_sd": ""},
    ])
    rep = pd.concat([rep, meta], ignore_index=True)
    rep.to_csv(out / "align_report.tsv", sep="\t", index=False)
    shutil.copy2(out / "align_report.tsv", out_fid / "align_report.tsv")

    print(rep.to_string(index=False))
    print(f"\nphenotyped {n_export} -> phenotyped-and-genotyped {len(keep)} "
          f"({int(fam_sizes.size)} families, {int((fam_sizes > 1).sum())} "
          f"multi-member); EUR arm {int(n_eur)}")
    print(f"wrote {out}  (FID = family id)")
    print(f"wrote {out_fid}  (FID = IID, for GCTA)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
