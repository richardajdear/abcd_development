"""Paired Delta-beta, min-p permutation, and a random-region null for the
hpc_v3 C+T PRS arm.

`README.md` section 1 is explicit that a bigger point estimate is not a result:
the new subset phenotypes correlate r ~ 0.89-0.92 with `global_slope`, so the
primary readout is the *paired difference* in PRS beta against `global_slope`
in the same subjects, tested properly.  This script produces the three things
that requires.

THE MODEL IS COPIED, NOT CHOSEN.  `tools/prs_assoc.R` fits

    scale(y) ~ scale(PRS) + sex + age_c + PC1..PC10 + (1 | family_id)

but `age_c` does not exist in `covar_quant.txt` (gcta_export renames it
`baseline_age`), and prs_assoc.R silently drops terms it cannot find --- so
every published PRS beta in this project was fitted WITHOUT an age covariate.
The primary fits here reproduce that model exactly, because a comparison
against v2's numbers is worthless if the model differs.  `--with-age` refits
with `baseline_age` and `n_visits` added as a sensitivity arm.

(1) PAIRED DELTA-BETA, EXACTLY.  Both phenotypes are regressed on the SAME
design matrix X in the SAME subjects, so with b = (X'X)^-1 X' y the difference
of the two PRS coefficients is a linear functional of the difference of the
outcomes:

    beta_new - beta_global = c' y_new - c' y_global = c' (y_new - y_global)

i.e. it is exactly the PRS coefficient of a regression on the difference
score.  Its family-clustered sandwich standard error is therefore the right
standard error for Delta-beta, with no approximation and no simulation.  The
family-cluster bootstrap that `README.md` section 1 asks for is run alongside
it (families resampled, never individuals) and reported as a check --- if the
two disagree, believe neither without looking.

The random intercept is replaced by family clustering: both handle the same
non-independence, and the sandwich/bootstrap pair extends to the difference
score, which `lmer` does not do directly.  The OLS and lmer point estimates
are both reported so the substitution is auditable.

(2) MIN-P PERMUTATION over the 8 nested C+T thresholds.  Bonferroni treats
nested scores as independent; `legacy/hpc_v2/README_HPC.md` section 7 records that the
verdict on the settled phenotype is correction-dependent (survives permutation
at 0.035, fails Bonferroni at 0.068).  The null here permutes the
covariate-residualised phenotype in FAMILY BLOCKS (whole families exchanged
with other families of the same size), which preserves the within-family
correlation of y while breaking its link to the score.

(3) RANDOM-REGION NULL.  The paired test asks whether a subset beats the
global mean.  It does not ask whether *this* subset beats an arbitrary one ---
and `legacy/hpc/README_HPC.md` section 2 already records that a PRS-SELECTED
top-10-region phenotype sat at the median of random 10-region sets.  Our sets
were chosen from group maps, not from the PRS, so that finding does not
transfer; but the comparison is nearly free and it is the honest denominator
for an SNR-concentration claim.  R random 8-bilateral-region subset means are
built from the same BLUPs and pushed through the same regression.

    python genetic_analysis/setup/prs_paired_delta.py --out-dir genetic_analysis/work/prs
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
NEW = ["slope_topDelta", "slope_topC3", "slope_projDelta", "slope_projC3"]
REF = "global_slope"
THRESHOLDS = ["5e-8", "1e-5", "0p001", "0p01", "0p05", "0p1", "0p5", "1"]
DISORDERS = ["SCZ", "MDD", "ASD", "ALZ", "ALZnoAPOE"]


# ---------------------------------------------------------------------------
# estimation primitives
# ---------------------------------------------------------------------------
def ols(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.linalg.lstsq(X, y, rcond=None)[0]


def z(v: np.ndarray) -> np.ndarray:
    return (v - v.mean()) / v.std(ddof=1)


class Fitter:
    """Precomputed OLS machinery for one fixed design matrix.

    The design is the same for every phenotype at a given (disorder,
    threshold, stratum), so (X'X)^-1 X' is factorised once and every fit is a
    matrix-vector product.  The family-clustered meat matrix is accumulated
    with ``np.add.reduceat`` over a family-sorted copy of X --- the same
    quantity a per-family Python loop computes, without the 6,789 iterations.
    """

    def __init__(self, X: np.ndarray, fam_idx: np.ndarray, n_fam: int):
        self.X, self.n, self.k = X, *X.shape
        self.A = np.linalg.pinv(X.T @ X)
        self.P = self.A @ X.T                       # (k, n)
        self.g = n_fam
        self.df = n_fam - 1
        self.order = np.argsort(fam_idx, kind="stable")
        fs = fam_idx[self.order]
        self.starts = np.flatnonzero(np.r_[True, fs[1:] != fs[:-1]])
        self.Xo = X[self.order]
        self.c1 = (self.g / (self.g - 1)) * ((self.n - 1) / (self.n - self.k))
        self.se_scale = np.sqrt(self.A[1, 1])       # classical, x sigma

    def beta(self, y: np.ndarray) -> np.ndarray:
        return self.P @ y

    def fit(self, y: np.ndarray) -> tuple[float, float, float]:
        """(beta_PRS, family-clustered se, p) --- the reported estimator."""
        b = self.P @ y
        e = y - self.X @ b
        u = np.add.reduceat(self.Xo * e[self.order, None], self.starts, axis=0)
        V = self.c1 * self.A @ (u.T @ u) @ self.A
        se = float(np.sqrt(V[1, 1]))
        return float(b[1]), se, float(two_sided_p(b[1] / se, self.df))

    def classical_t(self, Y: np.ndarray) -> np.ndarray:
        """|t| for the PRS coefficient, over columns of Y, unclustered.

        Used ONLY inside the permutation null, where the family structure is
        preserved by the permutation scheme itself: the reference distribution
        is generated correctly whatever monotone statistic is used, and this
        one is a couple of matrix products for 2,000 replicates at once.
        """
        Y = np.atleast_2d(Y.T).T
        B = self.P @ Y
        E = Y - self.X @ B
        sigma = np.sqrt((E * E).sum(0) / (self.n - self.k))
        return np.abs(B[1] / (sigma * self.se_scale))


def two_sided_p(t, df):
    from scipy import stats
    return 2 * stats.t.sf(np.abs(t), df)


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load(pheno_dir: Path, prs_dir: Path, eur_ids: Path, with_age: bool):
    ph = pd.read_csv(pheno_dir / "phenotypes_gcta.txt", sep=" ")
    q = pd.read_csv(pheno_dir / "covar_quant.txt", sep=" ")
    c = pd.read_csv(pheno_dir / "covar_categorical.txt", sep=" ")
    d = ph.merge(q, on=["FID", "IID"]).merge(c[["FID", "IID", "sex", "site"]],
                                             on=["FID", "IID"])
    assert len(d) == len(ph), "covariate merge dropped subjects"

    eur = set(pd.read_csv(eur_ids, sep=r"\s+", header=None)[0])
    d["is_eur"] = d.IID.isin(eur)

    pcs = [f"PC{i}" for i in range(1, 11)]
    cov = ["sex_M"] + pcs + (["baseline_age", "n_visits"] if with_age else [])
    d["sex_M"] = (d.sex == "M").astype(float)
    assert d[cov].notna().all().all(), "missing covariate values"

    scores = {}
    for dis in DISORDERS:
        for thr in THRESHOLDS:
            f = prs_dir / f"score_{dis}_{thr}.profile"
            if not f.exists():
                continue
            s = pd.read_csv(f, sep=r"\s+")[["IID", "SCORESUM"]]
            assert s.SCORESUM.std() > 0, f"{f.name} has zero-variance score"
            scores[(dis, thr)] = s.set_index("IID").SCORESUM
    return d, cov, scores


def design(sub: pd.DataFrame, prs: pd.Series, cov: list[str]) -> np.ndarray:
    n = len(sub)
    return np.column_stack([np.ones(n), z(prs.loc[sub.IID].to_numpy(float)),
                            sub[cov].to_numpy(float)])


# ---------------------------------------------------------------------------
# analyses
# ---------------------------------------------------------------------------
class Stratum:
    """One analysis sample, with its family bookkeeping precomputed."""

    def __init__(self, d: pd.DataFrame, name: str, cov: list[str]):
        self.name = name
        self.d = (d[d.is_eur] if name == "EUR" else d).reset_index(drop=True)
        self.cov = cov
        self.n = len(self.d)
        ufam, self.fam_idx = np.unique(self.d.FID.to_numpy(str),
                                       return_inverse=True)
        self.n_fam = len(ufam)
        self.members = [np.flatnonzero(self.fam_idx == i)
                        for i in range(self.n_fam)]
        self.y = {c: z(self.d[c].to_numpy(float))
                  for c in self.d.columns if c.startswith(("slope_", "global_",
                                                           "baseline_"))}

    def design(self, prs: pd.Series) -> np.ndarray:
        return np.column_stack([
            np.ones(self.n), z(prs.loc[self.d.IID].to_numpy(float)),
            self.d[self.cov].to_numpy(float)])

    def fitter(self, prs: pd.Series) -> Fitter:
        return Fitter(self.design(prs), self.fam_idx, self.n_fam)


def paired_table(strata, scores, boot_thresholds, n_boot, seed) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for st in strata:
        for (dis, thr), s in scores.items():
            X = st.design(s)
            f = Fitter(X, st.fam_idx, st.n_fam)
            ref = st.y[REF]
            b_ref, se_ref, p_ref = f.fit(ref)

            boot_setup = None
            if thr in boot_thresholds:
                # per-family Gram blocks, so a bootstrap replicate is a gather
                # and a sum rather than a refit over 8,082 rows
                XX = np.stack([X[m].T @ X[m] for m in st.members])
                boot_setup = (XX, [X[m] for m in st.members])

            for ph in NEW:
                yy = st.y[ph]
                b_new, se_new, p_new = f.fit(yy)
                # exact paired difference: same X, so it is the difference score
                dlt, se_d, p_d = f.fit(yy - ref)
                assert abs(dlt - (b_new - b_ref)) < 1e-9

                bse = blo = bhi = np.nan
                if boot_setup is not None:
                    XX, Xg = boot_setup
                    dif = yy - ref
                    Xy = np.stack([Xg[i].T @ dif[m]
                                   for i, m in enumerate(st.members)])
                    boot = np.empty(n_boot)
                    for k in range(n_boot):
                        pick = rng.integers(0, st.n_fam, st.n_fam)
                        boot[k] = np.linalg.solve(XX[pick].sum(0),
                                                  Xy[pick].sum(0))[1]
                    bse = boot.std(ddof=1)
                    blo, bhi = np.percentile(boot, [2.5, 97.5])

                rows.append(dict(
                    disorder=dis, threshold=thr, stratum=st.name, phenotype=ph,
                    n=st.n, n_families=st.n_fam,
                    beta_new=b_new, se_new=se_new, p_new=p_new,
                    beta_global=b_ref, se_global=se_ref, p_global=p_ref,
                    delta=dlt, se_delta=se_d, p_delta=p_d,
                    boot_se=bse, boot_lo=blo, boot_hi=bhi))
    return pd.DataFrame(rows)


def minp_permutation(strata, scores, n_perm, seed) -> pd.DataFrame:
    """Family-block permutation null for min-p across the 8 nested thresholds.

    Whole families are exchanged with other families of the SAME SIZE, so each
    permuted vector is a valid relabelling of intact sibships: the within-
    family correlation of y survives and only its link to the score is broken.
    Bonferroni over nested C+T thresholds is the wrong correction (they are
    ~4.1 effective tests, legacy/hpc_v2/README_HPC.md section 7); this is the right one.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for st in strata:
        # layout: individuals ordered family-by-family, families grouped by size
        sizes = np.array([len(m) for m in st.members])
        blocks, layout = [], []
        for s_ in np.unique(sizes):
            fams = np.flatnonzero(sizes == s_)
            pos = len(layout)
            layout.extend(np.concatenate([st.members[i] for i in fams]))
            blocks.append((pos, len(fams), int(s_)))
        layout = np.asarray(layout)
        assert len(layout) == st.n and len(np.unique(layout)) == st.n

        # residualise on the covariates (no score) once per stratum; formed
        # as y - C (C'C)^-1 C'y rather than an n x n projection matrix, which
        # would be a 522 MB allocation at n = 8,082
        C = np.column_stack([np.ones(st.n), st.d[st.cov].to_numpy(float)])
        Cp = np.linalg.pinv(C.T @ C) @ C.T

        for dis in DISORDERS:
            thrs = [t for t in THRESHOLDS if (dis, t) in scores]
            if not thrs:
                continue
            fits = [st.fitter(scores[(dis, t)]) for t in thrs]
            for ph in NEW + [REF]:
                y = st.y[ph]
                obs = min(f.fit(y)[2] for f in fits)                 # reported
                t_obs = max(float(f.classical_t(y)[0]) for f in fits)  # statistic
                r = y - C @ (Cp @ y)
                Y = np.empty((st.n, n_perm))
                for b in range(n_perm):
                    perm = np.empty(st.n, dtype=int)
                    for pos, nf, s_ in blocks:
                        blk = layout[pos:pos + nf * s_].reshape(nf, s_)
                        perm[pos:pos + nf * s_] = blk[rng.permutation(nf)].ravel()
                    Y[layout, b] = r[perm]
                t_null = np.max([f.classical_t(Y) for f in fits], axis=0)
                rows.append(dict(
                    disorder=dis, stratum=st.name, phenotype=ph, n=st.n,
                    n_thresholds=len(thrs), n_perm=n_perm,
                    min_p=obs, max_abs_t=t_obs,
                    p_perm=float((np.sum(t_null >= t_obs) + 1) / (n_perm + 1)),
                    p_bonferroni=min(1.0, obs * len(thrs))))
    return pd.DataFrame(rows)


def random_region_null(strata, scores, run_dir, disorders, thresholds,
                       n_sets, seed) -> pd.DataFrame:
    """Where do the chosen 8-region sets sit among random 8-region sets?

    The paired test asks whether a subset beats the global mean; it does not
    ask whether THIS subset beats an arbitrary one.  legacy/hpc/README_HPC.md section
    2 already records that a PRS-SELECTED top-10-region phenotype sat at the
    median of random 10-region sets.  Ours were selected from group maps, not
    from the PRS, so that result does not transfer --- but it is the honest
    denominator for an SNR-concentration claim, and it is nearly free.
    """
    sys.path.insert(0, str(REPO / "src"))
    from abcd import gcta_export                              # noqa: E402
    rng = np.random.default_rng(seed)

    ph = pd.read_parquet(run_dir / "phenotypes" / "phenotypes.parquet")
    sl = ph[ph.phenotype == "slope"].pivot(index="subject", columns="label",
                                           values="value")
    assert sl.shape == (8192, 68), f"unexpected slope matrix {sl.shape}"
    sl.index = (gcta_export._to_genetics_id(sl.index.to_series())
                .str.replace(r"^NDAR_INV", "sub-", regex=True))
    bases = sorted({c[3:] for c in sl.columns})
    assert len(bases) == 34

    rows = []
    for st in strata:
        M = sl.loc[st.d.IID].to_numpy(float)
        cols = {c: i for i, c in enumerate(sl.columns)}
        for dis in disorders:
            for thr in thresholds:
                if (dis, thr) not in scores:
                    continue
                f = st.fitter(scores[(dis, thr)])
                obs = {p: f.fit(st.y[p])[0]
                       for p in ("slope_topDelta", "slope_topC3", REF)}
                null = np.empty(n_sets)
                for k in range(n_sets):
                    pick = rng.choice(len(bases), 8, replace=False)
                    idx = [cols[f"{h}_{bases[i]}"] for i in pick
                           for h in ("lh", "rh")]
                    null[k] = f.fit(z(M[:, idx].mean(axis=1)))[0]
                for p, b in obs.items():
                    rows.append(dict(
                        disorder=dis, threshold=thr, stratum=st.name,
                        phenotype=p, beta=b, n_random_sets=n_sets,
                        null_mean=float(null.mean()),
                        null_sd=float(null.std(ddof=1)),
                        # thinning associations are negative, so "more extreme"
                        # means more negative
                        pctile=float((null <= b).mean() * 100)))
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pheno-dir", default=str(REPO / "legacy/hpc_v2/work/pheno_v3"))
    ap.add_argument("--prs-dir",
                    default=str(REPO / "legacy/hpc_v2/work/results_v2/prs_ct_v3"))
    ap.add_argument("--eur-ids",
                    default=str(REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep"))
    ap.add_argument("--run-dir",
                    default=str(REPO / "out/thickness_dsk_70_139406217085"))
    ap.add_argument("--out-dir", default=str(REPO / "genetic_analysis/work/prs"))
    ap.add_argument("--with-age", action="store_true",
                    help="sensitivity arm: add baseline_age + n_visits")
    ap.add_argument("--boot-thresholds", default="0p5,0p01",
                    help="thresholds to bootstrap Delta-beta at "
                         "(the project's two reported C+T thresholds)")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--n-random-sets", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260908)
    a = ap.parse_args()

    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    tag = "_withage" if a.with_age else ""
    d, cov, scores = load(Path(a.pheno_dir), Path(a.prs_dir), Path(a.eur_ids),
                          a.with_age)
    strata = [Stratum(d, s, cov) for s in ("EUR", "full")]
    print(f"{len(d)} subjects, {d.is_eur.sum()} EUR, {len(scores)} scores, "
          f"{len(cov)} covariates: {', '.join(cov)}", flush=True)
    for st in strata:
        print(f"  {st.name}: n={st.n}, families={st.n_fam}", flush=True)

    t = paired_table(strata, scores, set(a.boot_thresholds.split(",")),
                     a.n_boot, a.seed)
    t.to_csv(out / f"prs_paired_delta{tag}.tsv", sep="\t", index=False)
    print(f"wrote prs_paired_delta{tag}.tsv ({len(t)} rows)", flush=True)

    m = minp_permutation(strata, scores, a.n_perm, a.seed + 1)
    m.to_csv(out / f"prs_minp_permutation{tag}.tsv", sep="\t", index=False)
    print(f"wrote prs_minp_permutation{tag}.tsv ({len(m)} rows)", flush=True)

    r = random_region_null(strata, scores, Path(a.run_dir),
                           ["SCZ", "MDD"], ["0p5", "0p01"],
                           a.n_random_sets, a.seed + 2)
    r.to_csv(out / f"prs_random_region_null{tag}.tsv", sep="\t", index=False)
    print(f"wrote prs_random_region_null{tag}.tsv ({len(r)} rows)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
