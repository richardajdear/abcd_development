"""Synthetic fixture for testing the hpc_v2/ pipeline locally.

    python tools/make_test_genotypes.py          # build the v1 fixture first
    python tools/make_test_genotypes_v2.py       # then derive the v2 fixture

Everything lands in ``scratch/hpc_v2_test/`` (gitignored), shaped exactly as
``hpc_v2/config.local.sh.example`` expects:

    genotype/synthetic.{bed,bim,fam}    array fileset  = the v1 fixture, copied
    genotype_imp/synthetic_chr{1..22}.* "imputed" filesets = per-chromosome
                                        splits of the same genotypes (plink)
    grm/synthetic_full.grm.{bin,N.bin,id}  dense GRM over ALL subjects
                                        (plink --make-grm-bin; relatives IN)
    prs_profiles/score_SCZ_*.profile    synthetic PRS with a KNOWN within-family
                                        effect (see below)
    keep/eur.keep                       an arbitrary ~70% subject subset
    pheno/                              v1 export copied verbatim (real FIDs)

The PRS profiles are the one genuinely new piece.  They are built from the
fixture genotypes so that sibling scores are correlated (as real PRS are), and
the ``global_slope`` phenotype column is REWRITTEN as

    y = b_W * PRS_z + b_CONF * PRS_fammean_z + covariate noise

with b_W = -0.20 and b_CONF = -0.15: a large causal within-family effect plus
family-level confounding.  Large deliberately -- the fixture has ~700 pairs, so
an ABCD-sized effect (-0.04) would give the smoke test a ~10% chance of seeing
its own positive control.  At -0.20 the within-family z is ~5, so:

    - 06_prs_family must find beta_within ~ -0.20, significant;
    - beta_between must exceed it in magnitude (it absorbs the confounder);
    - beta_diff must be significantly nonzero.

That validates the decomposition end to end: a script that mixes up b_B/b_W,
loses the family IDs, or scales in the wrong place cannot pass.

Requires: plink (for the split + GRM), numpy, pandas.  No GCTA, no R.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
SEED = 20260823

B_WITHIN = -0.20   # causal within-family PRS effect (SD per SD)
B_CONF = -0.15     # extra family-mean effect = confounding the test must expose
PRS_THRESHOLDS = ["0.05", "0.5", "1"]  # mirrors the v1 local-test thresholds


def sh(*cmd: str) -> None:
    print("+", " ".join(map(str, cmd)))
    r = subprocess.run(list(map(str, cmd)), capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"FAILED ({r.returncode}):\n{r.stdout}\n{r.stderr}")


def read_bed_dosages(prefix: Path) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    """Minimal PLINK1 .bed reader (variant-major) -> dosage matrix n x m."""
    fam = pd.read_csv(prefix.with_suffix(".fam"), sep=r"\s+", header=None,
                      names=["FID", "IID", "PID", "MID", "SEX", "PHE"], dtype=str)
    bim = pd.read_csv(prefix.with_suffix(".bim"), sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"], dtype=str)
    n, m = len(fam), len(bim)
    raw = np.fromfile(prefix.with_suffix(".bed"), dtype=np.uint8)
    assert raw[0] == 0x6C and raw[1] == 0x1B and raw[2] == 0x01, "not variant-major PLINK1"
    body = raw[3:]
    bpv = (n + 3) // 4
    assert body.size == bpv * m, "bed size mismatch (run check_bfile_integrity)"
    # Decode 2-bit genotypes: 00=hom A1(2 doses), 10=het(1), 11=hom A2(0), 01=missing
    blocks = body.reshape(m, bpv)
    codes = np.zeros((m, bpv * 4), dtype=np.uint8)
    for shift, col in zip((0, 2, 4, 6), range(4)):
        codes[:, col::4] = (blocks >> shift) & 0b11
    codes = codes[:, :n]
    dos = np.empty((m, n), dtype=np.float32)
    dos[codes == 0] = 2.0
    dos[codes == 2] = 1.0
    dos[codes == 3] = 0.0
    dos[codes == 1] = np.nan
    return dos.T.copy(), fam, bim  # n x m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--v1-fixture", default=REPO / "scratch/hpc_test", type=Path)
    ap.add_argument("--out", default=REPO / "scratch/hpc_v2_test", type=Path)
    ap.add_argument("--plink", default="plink")
    args = ap.parse_args()

    v1 = args.v1_fixture
    out = args.out
    src = v1 / "genotype/synthetic"
    if not src.with_suffix(".bed").exists():
        sys.exit(f"v1 fixture not found at {src}.bed -- run tools/make_test_genotypes.py first")

    rng = np.random.default_rng(SEED)

    # --- 1. array fileset + phenotype export, copied from v1 -------------------
    (out / "genotype").mkdir(parents=True, exist_ok=True)
    for ext in (".bed", ".bim", ".fam"):
        shutil.copy2(src.with_suffix(ext), out / "genotype" / f"synthetic{ext}")
    shutil.copytree(v1 / "pheno", out / "pheno", dirs_exist_ok=True)

    # --- 2. per-chromosome "imputed" filesets ----------------------------------
    imp = out / "genotype_imp"
    imp.mkdir(exist_ok=True)
    bim = pd.read_csv(src.with_suffix(".bim"), sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"], dtype=str)
    for chrom in sorted(bim.CHR.astype(int).unique()):
        sh(args.plink, "--bfile", src, "--chr", str(chrom),
           "--make-bed", "--out", imp / f"synthetic_chr{chrom}")

    # --- 3. dense full-sample GRM (GCTA binary format, via plink) --------------
    (out / "grm").mkdir(exist_ok=True)
    sh(args.plink, "--bfile", src, "--make-grm-bin", "--out", out / "grm/synthetic_full")

    # --- 4. synthetic PRS profiles with a known within-family effect ------------
    dos, fam, bim = read_bed_dosages(out / "genotype" / "synthetic")
    m = dos.shape[1]
    w = rng.normal(0, 1, m) * (rng.random(m) < 0.10)     # 10% of SNPs weighted
    dosf = np.where(np.isnan(dos), np.nanmean(dos, axis=0, keepdims=True), dos)
    prs_raw = dosf @ w
    prs_z = (prs_raw - prs_raw.mean()) / prs_raw.std()

    d = fam[["FID", "IID"]].copy()
    d["PRS"] = prs_z
    fam_mean = d.groupby("FID")["PRS"].transform("mean")

    (out / "prs_profiles").mkdir(exist_ok=True)
    for disorder in ("SCZ", "MDD"):
        for thr in PRS_THRESHOLDS:
            # Different thresholds = same score + a little noise, as real
            # nested P+T scores are; MDD = an independent null score.
            if disorder == "SCZ":
                score = prs_z + rng.normal(0, 0.05, len(prs_z))
            else:
                score = rng.normal(0, 1, len(prs_z))
            prof = pd.DataFrame({
                "FID": d.FID, "IID": d.IID, "PHENO": -9,
                "CNT": m, "CNT2": int(m * 0.1),
                "SCORESUM": score,
            })
            p = out / "prs_profiles" / f"score_{disorder}_{thr}.profile"
            prof.to_csv(p, sep=" ", index=False)
    print(f"wrote {2 * len(PRS_THRESHOLDS)} PRS profiles")

    # --- 5. rewrite global_slope with the known effect --------------------------
    ph_path = out / "pheno/phenotypes_gcta.txt"
    ph = pd.read_csv(ph_path, sep=r"\s+", dtype={"FID": str, "IID": str})
    key = ph.IID.map(dict(zip(d.IID, range(len(d)))))
    assert key.notna().all(), "phenotype IDs not all present in .fam"
    idx = key.astype(int).to_numpy()
    y = (B_WITHIN * (prs_z[idx] - fam_mean.to_numpy()[idx])
         + (B_WITHIN + B_CONF) * fam_mean.to_numpy()[idx]
         + rng.normal(0, np.sqrt(1 - B_WITHIN**2 - (B_WITHIN + B_CONF)**2 / 2),
                      len(ph)))
    ph["global_slope"] = y
    ph.to_csv(ph_path, sep=" ", index=False)
    print(f"rewrote global_slope in {ph_path} with b_W={B_WITHIN}, "
          f"b_B={B_WITHIN + B_CONF} (family-mean coefficient)")

    # --- 6. EUR keep list --------------------------------------------------------
    (out / "keep").mkdir(exist_ok=True)
    fids = d.FID.unique()
    eur_fids = set(rng.choice(fids, size=int(len(fids) * 0.7), replace=False))
    eur = d[d.FID.isin(eur_fids)]
    eur[["FID", "IID"]].to_csv(out / "keep/eur.keep", sep=" ",
                               header=False, index=False)
    print(f"EUR keep list: {len(eur)} of {len(d)} subjects "
          f"({len(eur_fids)} whole families)")

    # --- 7. summary ---------------------------------------------------------------
    n_multi = (d.groupby("FID").size() >= 2).sum()
    print(f"\nfixture ready at {out}")
    print(f"  subjects: {len(d)}, multi-member families: {n_multi}")
    print("  point hpc_v2/config.local.sh at it "
          "(copy hpc_v2/config.local.sh.example)")


if __name__ == "__main__":
    main()
