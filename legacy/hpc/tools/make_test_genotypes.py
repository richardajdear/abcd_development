"""Synthetic genotypes for testing the hpc/ pipeline locally.

    python tools/make_test_genotypes.py [--n-subjects 1500] [--n-snps 5000]

Writes a PLINK1 binary fileset (``.bed``/``.bim``/``.fam``) plus the exported
phenotypes and covariates for the same subjects, into ``scratch/hpc_test/``
(gitignored).  With ``hpc/config.local.sh`` pointed at that directory, the whole
pipeline runs end to end on a laptop in a couple of minutes.

**What this does and does not test.**  The genotypes are simulated, so nothing
here says anything biological -- an h2 from this fixture is a property of the
simulation.  What it does test is every way the pipeline can be wired up wrong:
that exported IDs intersect the ``.fam``, that ``--mpheno`` selects the column
the manifest claims, that GCTA accepts the covariate files, that the summary
collectors parse real GCTA output rather than assumed formats, and that the
scripts run under ``run_all.sh`` in local mode.  Those are the failures that
would otherwise be discovered hours into a cluster job.

Two properties are deliberately inherited from the real data rather than
invented:

*Family structure.*  ``FID`` comes from ``model_table.parquet``, so the fixture
contains genuine twin and sibling groupings.  A GRM built on it therefore has a
realistic block of high off-diagonal values, ``--grm-cutoff`` actually removes
somebody, and the difference between the REML N and the fastGWA N shows up.
Relatives are given correlated genotypes (transmission from simulated parental
haplotypes), because unrelated genotypes under real family labels would let a
GRM-cutoff bug pass unnoticed.

*Subject IDs.*  Taken verbatim from the export, in the ``NDAR_INV`` form, which
is what makes the ID-intersection check meaningful.  Writing invented IDs would
turn the one check most likely to catch a real problem into a tautology.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import gcta_export, paths  # noqa: E402

SEED = 20260810

#: Fraction of causal SNPs, and the heritability they jointly explain in the
#: simulated phenotype.  Used only for the seeded-effect check: the pipeline
#: should recover *something* at these SNPs, confirming the phenotype column
#: actually reached the association test.
N_CAUSAL = 20
H2_SIM = 0.5


def _pick_subjects(pheno: pd.DataFrame, n: int, rng) -> pd.DataFrame:
    """Subsample whole families, never splitting one, multi-member first.

    Two requirements, in tension.  Families must not be split, or the inherited
    family structure is destroyed.  But most ABCD families contribute a single
    phenotyped subject, so sampling families uniformly yields a fixture that is
    almost entirely singletons -- and a fixture with few relatives cannot
    exercise ``--grm-cutoff`` or show the REML/fastGWA N difference, which is
    most of what the fixture exists to test.

    So: take every multi-member family first, then fill with singletons.  The
    result is deliberately *enriched* for relatedness relative to the real
    cohort.  That is the right bias for a plumbing test and the wrong bias for
    anything quantitative -- another reason no number from this fixture means
    anything biologically.
    """
    sizes = pheno.groupby("FID", sort=False).size()
    multi = sizes[sizes > 1].index.to_numpy()
    single = sizes[sizes == 1].index.to_numpy()
    rng.shuffle(multi)
    rng.shuffle(single)

    out, total = [], 0
    for f in np.concatenate([multi, single]):
        block = pheno[pheno.FID == f]
        if total + len(block) > n:
            if total >= n:
                break
            continue          # skip this one, a smaller family may still fit
        out.append(block)
        total += len(block)
    return pd.concat(out, ignore_index=True)


def simulate(pheno: pd.DataFrame, n_snps: int, rng) -> tuple[np.ndarray, np.ndarray]:
    """Genotype dosages (0/1/2), with real sibling relatedness within family.

    Explicit haplotype transmission: each family gets two parents, each parent
    two haplotypes drawn at population frequency, and each child inherits one
    haplotype from each parent independently per SNP.  Siblings therefore share
    ~50% of alleles by descent and the GRM's within-family off-diagonal lands
    near 0.5, as in the real data.

    The first attempt here instead perturbed the allele frequency per family and
    drew members independently around it.  That is much simpler and produces a
    within-family GRM off-diagonal of 0.01 -- i.e. no relatedness at all, since
    a shared frequency shift of that size induces almost no allele sharing.  A
    fixture built that way would let a broken ``--grm-cutoff`` pass, which is
    precisely the check the fixture exists to exercise.  Transmission is the
    cheapest construction that actually produces relatives.
    """
    maf = rng.uniform(0.05, 0.5, n_snps)

    n = len(pheno)
    G = np.empty((n, n_snps), dtype=np.int8)

    for _fid, block in pheno.groupby("FID", sort=False):
        idx = block.index.to_numpy()
        if len(idx) == 1:
            # Unrelated singleton: straight from population frequency.
            G[idx[0]] = rng.binomial(2, maf)
            continue
        # Four parental haplotypes (2 parents x 2 each).
        hap = rng.random((4, n_snps)) < maf          # bool: carries A1
        for i in idx:
            # One haplotype from each parent, chosen per SNP.
            from_mum = np.where(rng.random(n_snps) < 0.5, hap[0], hap[1])
            from_dad = np.where(rng.random(n_snps) < 0.5, hap[2], hap[3])
            G[i] = from_mum.astype(np.int8) + from_dad.astype(np.int8)
    return G, maf


def simulate_phenotype(G: np.ndarray, rng) -> tuple[np.ndarray, np.ndarray]:
    """A phenotype with a known genetic component, for the recovery check."""
    n, m = G.shape
    causal = rng.choice(m, size=min(N_CAUSAL, m), replace=False)
    Z = (G[:, causal] - G[:, causal].mean(0)) / (G[:, causal].std(0) + 1e-9)
    beta = rng.normal(0, 1, len(causal))
    g = Z @ beta
    g = (g - g.mean()) / g.std()
    e = rng.normal(0, 1, n)
    y = np.sqrt(H2_SIM) * g + np.sqrt(1 - H2_SIM) * e
    return y, causal


def write_plink(out_prefix: Path, pheno: pd.DataFrame, G: np.ndarray,
                maf: np.ndarray) -> None:
    """Write .bed/.bim/.fam by hand (no PLINK binary required).

    PLINK1 .bed packs four genotypes per byte, two bits each, in the order
    00=hom A1, 10=het, 11=hom A2, 01=missing -- note that 01 is *missing*, not
    a genotype, which is the detail that makes a naive 0/1/2 -> 0b00/0b01/0b10
    mapping produce a file that loads without error and is wrong.
    """
    n, m = G.shape
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # .fam: FID IID PAT MAT SEX PHENO  (-9 = missing phenotype)
    #
    # SEX must be joined from the categorical covariates, not read off the
    # phenotype frame -- that frame is FID/IID/phenotypes only, so a `.get("sex")`
    # on it silently falls back and writes every subject as female (GCTA then
    # reports "0 males", which is how this was caught).  Autosomal GRMs do not
    # use SEX, but writing a .fam that contradicts the covariate file is exactly
    # the kind of inconsistency a fixture should not normalise away.
    if "sex" not in pheno.columns:
        raise ValueError("write_plink needs a 'sex' column; merge it before calling")
    # astype(object) first: `sex` is a pandas Categorical, and .map() on one
    # returns a Categorical whose categories are {1,2}, so .fillna(0) raises
    # rather than filling.  PLINK's code for unknown sex is 0.
    sex_code = (pheno.sex.astype(object)
                .map({"M": 1, "F": 2})
                .fillna(0).astype(int))
    fam = pd.DataFrame({
        "FID": pheno.FID, "IID": pheno.IID,
        "PAT": 0, "MAT": 0,
        "SEX": sex_code,
        "PHENO": -9,
    })
    fam.to_csv(out_prefix.with_suffix(".fam"), sep=" ", header=False, index=False)

    # .bim: CHR SNP CM POS A1 A2, spread over 22 autosomes so --autosome keeps
    # them and MAGMA's positional annotation has something plausible to use.
    #
    # Sorted by chromosome, then position.  Interleaving them (chrom = i % 22)
    # is the obvious way to spread SNPs evenly and produces a file PLINK
    # rejects outright -- "split chromosome" -- because a real .bim is always in
    # genomic order.  GCTA is more tolerant, so this only shows up with PLINK,
    # which is reason enough to keep the fixture in the stricter format.
    chrom = np.repeat(np.arange(1, 23), int(np.ceil(m / 22)))[:m]
    pos = np.concatenate([
        100_000 + np.arange((chrom == c).sum()) * 50_000
        for c in range(1, 23) if (chrom == c).any()
    ])
    bim = pd.DataFrame({
        "CHR": chrom,
        "SNP": [f"rs{i:07d}" for i in range(m)],
        "CM": 0,
        "POS": pos,
        "A1": "A",
        "A2": "G",
    })
    bim.to_csv(out_prefix.with_suffix(".bim"), sep="\t", header=False, index=False)

    # .bed, SNP-major.
    code = np.array([0b00, 0b10, 0b11], dtype=np.uint8)   # 0,1,2 copies of A1
    with open(out_prefix.with_suffix(".bed"), "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))               # magic + SNP-major
        pad = (-n) % 4
        for j in range(m):
            col = code[G[:, j]]
            if pad:
                col = np.concatenate([col, np.zeros(pad, dtype=np.uint8)])
            packed = (col[0::4] | (col[1::4] << 2) | (col[2::4] << 4) | (col[3::4] << 6))
            fh.write(packed.astype(np.uint8).tobytes())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", default=None,
                    help="fitted run to take subjects/phenotypes from "
                         "(default: the run for $ABCD_CONFIG)")
    ap.add_argument("--n-subjects", type=int, default=1500)
    ap.add_argument("--n-snps", type=int, default=5000)
    ap.add_argument("--out", default=None,
                    help="default: <repo>/scratch/hpc_test")
    a = ap.parse_args(argv)

    from abcd.config import active_run_dir
    run_dir = Path(a.run_dir) if a.run_dir else active_run_dir()
    out_root = Path(a.out) if a.out else paths.REPO_ROOT / "scratch" / "hpc_test"

    rng = np.random.default_rng(SEED)

    # Real exported phenotypes/covariates: the same code path the cluster uses,
    # so a change to the export is reflected in the fixture automatically.
    built = gcta_export.build(run_dir)
    phen_all = built["phenotypes_gcta"]

    phen = _pick_subjects(phen_all.reset_index(drop=True), a.n_subjects, rng)
    phen = phen.reset_index(drop=True)

    # Sex lives in the categorical covariates; the .fam needs it (see write_plink).
    phen = phen.merge(built["covar_categorical"][["IID", "sex"]], on="IID", how="left")

    G, maf = simulate(phen, a.n_snps, rng)
    y_sim, causal = simulate_phenotype(G, rng)

    geno_prefix = out_root / "genotype" / "synthetic"
    write_plink(geno_prefix, phen, G, maf)

    # Phenotype file: the five real phenotypes, plus one simulated phenotype
    # with a known genetic basis.  The real ones test the plumbing; the
    # simulated one is the only column on which a *positive* result is
    # meaningful, so it is what confirms the association step works at all.
    pheno_dir = out_root / "pheno"
    pheno_dir.mkdir(parents=True, exist_ok=True)
    # Drop the covariate columns merged in for the .fam: they must not become
    # phenotype columns, or the manifest would list 'sex' as a phenotype and
    # --mpheno indices would shift.
    phen_out = phen.drop(columns=["sex"]).copy()
    phen_out["sim_h2_50"] = y_sim

    ids = set(phen.IID)
    for name, frame in built.items():
        if name == "phenotype_manifest":
            continue
        sub = frame[frame.IID.isin(ids)] if "IID" in frame else frame
        if name == "phenotypes_gcta":
            sub = phen_out
        sub.to_csv(pheno_dir / f"{name}.txt", sep=" ", index=False, na_rep="NA")

    # Manifest: rebuild it for the fixture's phenotype file, including the
    # simulated column, so mpheno indices match this file rather than the
    # cluster one.
    man = pd.DataFrame([
        {"name": c,
         "mpheno": phen_out.columns.get_loc(c) - 1,
         "priority": next((p["priority"] for p in gcta_export.PHENOTYPES
                           if p["name"] == c), 99),
         "role": next((p["role"] for p in gcta_export.PHENOTYPES
                       if p["name"] == c), "simulated"),
         "n_nonmissing": int(phen_out[c].notna().sum())}
        for c in phen_out.columns if c not in ("FID", "IID")
    ]).sort_values("priority", ignore_index=True)
    man.to_csv(pheno_dir / "phenotype_manifest.tsv", sep="\t", index=False)

    truth = pd.DataFrame({"snp": [f"rs{i:07d}" for i in causal]})
    truth.to_csv(out_root / "causal_snps.txt", index=False, header=False)

    n_fam = phen.FID.nunique()
    print(f"run_dir          {run_dir.name}")
    print(f"subjects         {len(phen)} in {n_fam} families "
          f"({len(phen) - n_fam} share a family with someone)")
    print(f"snps             {a.n_snps} over 22 autosomes")
    print(f"genotypes        {geno_prefix}.{{bed,bim,fam}}")
    print(f"phenotypes       {pheno_dir}/phenotypes_gcta.txt "
          f"({len(man)} columns incl. sim_h2_50)")
    print(f"causal snps      {out_root}/causal_snps.txt "
          f"({len(causal)} SNPs, simulated h2={H2_SIM})")
    print()
    print("Point the pipeline at it with hpc/config.local.sh:")
    print(f'  ABCD_HPC_ROOT="{out_root}"')
    print(f'  GENO="{geno_prefix}"')
    return 0


if __name__ == "__main__":
    sys.exit(main())
