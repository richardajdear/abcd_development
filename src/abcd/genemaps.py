"""
Linking cortical maps to gene expression, and gene lists to GWAS.

The chain this module implements
-------------------------------
1. A regional map of developmental change (from ``R/fit_lmm.R``).
2. Its spatial correlation with reference expression components -- AHBA
   C1-C3 and the snRNA-seq PC1 -- tested against a spin null
   (:mod:`abcd.spatial`).
3. A ranked gene list obtained by correlating each gene's regional
   expression with the developmental map (PLS/correlation weights).
4. Enrichment of that ranked list in GWAS signal for SCZ and MDD.

Step 3 is where most of the risk lies, so two things are done explicitly.
First, the gene-to-map correlation is computed on the SAME regions as the
map, with the region set recorded, because AHBA coverage is partial (137 of
180 HCP-MMP parcels have expression data) and quietly using a different
region set for different genes would create spurious structure.  Second, the
null for step 4 is a *competitive* one: enrichment is tested against other
genes, not against a random-permutation baseline that ignores the fact that
gene sets differ systematically in length, expression level, and LD.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from . import paths, spatial


# --------------------------------------------------------------------------
# Reference maps and gene weights
# --------------------------------------------------------------------------

def ahba_components(parcellation: str = "hcp") -> pd.DataFrame:
    """AHBA differential-expression component scores per region.

    Returns a frame indexed by region (no hemisphere -- AHBA donors are
    overwhelmingly left-hemisphere, so the components are bilateral by
    construction) with columns ``C1, C2, C3``.
    """
    if parcellation != "hcp":
        raise NotImplementedError(
            "AHBA components are supplied in HCP-MMP space. Fitting the "
            "imaging model in 'hcp' avoids cross-parcellation averaging; "
            "see Release51Adapter._imaging_hcp."
        )
    p = paths.DATA_DIR / "ahba_dme_hcp_top8kgenes_scores.csv"
    df = pd.read_csv(p).set_index("label")
    return df[["C1", "C2", "C3"]]


def ahba_gene_weights() -> pd.DataFrame:
    """Gene loadings on the AHBA components (``weights.csv``, 7,973 genes)."""
    return pd.read_csv(paths.DATA_DIR / "weights.csv", index_col=0)


def snrnaseq_pc1(which: str = "PC1_herringV3") -> pd.Series:
    """snRNA-seq maturation PC1 gene loadings.

    Two independent datasets are supplied (``PC1_herringV3``, ``PC1_U01V2``);
    they correlate at rho = 0.54, so results should be reported for both
    rather than for a single arbitrary choice.
    """
    df = pd.read_csv(paths.DATA_DIR / "velmeshev_PC1_gene_loadings.csv")
    if which not in df.columns:
        raise KeyError(f"{which!r} not in {[c for c in df.columns if 'PC' in c]}")
    return df.set_index("feature_name")[which].dropna()


def ahba_expression(parcellation: str = "hcp",
                    path: str | Path | None = None) -> pd.DataFrame:
    """Full AHBA gene-by-region expression matrix, regions x genes.

    This is the matrix behind the published components, processed with the
    donor-normalisation / probe-selection pipeline of Dear et al. (2024).
    Prefer it over :func:`ahba_components` for gene-level work: the three
    components span only a 3-dimensional subspace of the 7,973-gene space,
    so any map projected through them yields a gene ranking that is a fixed
    linear combination of three vectors, no matter what the map looks like.

    The source files index regions by an integer code with no names attached.
    The row order is asserted to match :func:`ahba_components` by
    reconstructing the published component scores from expression x loadings;
    a mismatch raises rather than silently mislabelling every region.
    """
    if parcellation != "hcp":
        raise NotImplementedError(
            f"expression matrix wired for HCP-MMP only, got {parcellation!r}"
        )
    p = Path(path) if path is not None else paths.ahba_expression_path(parcellation)
    if not p.exists():
        raise FileNotFoundError(
            f"AHBA expression matrix not found at {p}. Set ABCD_AHBA_DIR to the "
            "directory holding {dk,hcp}_3d_ds5.csv."
        )
    expr = pd.read_csv(p, index_col=0)
    scores = ahba_components(parcellation)
    if len(expr) != len(scores):
        raise ValueError(
            f"expression has {len(expr)} regions, components have {len(scores)}"
        )
    expr.index = scores.index
    _assert_row_order(expr, scores)
    return expr


def _assert_row_order(expr: pd.DataFrame, scores: pd.DataFrame,
                      min_rho: float = 0.85) -> None:
    """Verify expression rows align with component scores, or raise.

    Reconstructs each component as z(expression) @ loadings and requires it to
    track the published scores.  Row misalignment is the failure mode that
    would invalidate every downstream result while producing no error, so it
    is checked rather than assumed.
    """
    w = ahba_gene_weights()
    common = expr.columns.intersection(w.index)
    if len(common) < 0.9 * len(w):
        raise ValueError(
            f"only {len(common)} of {len(w)} weighted genes present in matrix"
        )
    z = (expr[common] - expr[common].mean()) / expr[common].std()
    recon = z.to_numpy() @ w.loc[common].to_numpy()
    bad = {}
    for j, comp in enumerate(w.columns):
        if comp not in scores:
            continue
        rho = stats.spearmanr(recon[:, j], scores[comp].to_numpy()).statistic
        if abs(rho) < min_rho:
            bad[comp] = round(float(rho), 4)
    if bad:
        raise ValueError(
            "AHBA expression rows do not align with component scores "
            f"(reconstruction rho {bad}, need |rho| >= {min_rho}). The region "
            "order in the expression file does not match the scores file."
        )


def bilateral(map_: pd.Series) -> pd.Series:
    """Average ``lh_X`` and ``rh_X`` into a hemisphere-free region index.

    Reference expression maps are bilateral, so the imaging map must be too.
    Regions present in only one hemisphere are kept with that value.
    """
    s = map_.copy()
    stem = s.index.str.replace(r"^(lh|rh)_", "", regex=True)
    return s.groupby(stem).mean()


# --------------------------------------------------------------------------
# Map-level association
# --------------------------------------------------------------------------

def map_vs_components(dev_map: pd.Series, geom: spatial.ParcelGeometry,
                      components: pd.DataFrame | None = None,
                      n_perm: int = 1000, seed: int = 0) -> pd.DataFrame:
    """Spin-tested correlation of a developmental map with AHBA components.

    ``dev_map`` is indexed by hemisphere-prefixed label; it is averaged to
    bilateral regions internally.  The spin is performed in the *bilateral*
    space by restricting the geometry to left-hemisphere centroids, which
    keeps rotation and map in the same index.
    """
    comps = ahba_components() if components is None else components
    dev_b = bilateral(dev_map)

    lh = [l for l in geom.labels if l.startswith("lh_")]
    stem = [l[3:] for l in lh]
    idx = [geom.labels.index(l) for l in lh]
    geom_b = spatial.ParcelGeometry(
        labels=tuple(stem), coords=geom.coords[idx], space=geom.space,
        hemi=tuple(["lh"] * len(lh)),
    )
    spins = spatial.spin_indices(geom_b, n_perm=n_perm, seed=seed,
                                hemi_paired=False)

    rows = []
    for comp in comps.columns:
        common = dev_b.index.intersection(comps.index)
        x = dev_b.loc[common]
        y = comps[comp].loc[common]
        pos = {l: i for i, l in enumerate(geom_b.labels)}
        keep = [l for l in common if l in pos]
        sub = np.array([pos[l] for l in keep])
        sp_sub = np.clip(np.searchsorted(sub, spins[:, sub]), 0, len(sub) - 1)
        t = spatial.spatial_corr(x.loc[keep], y.loc[keep], null="spin",
                                 spins=sp_sub, n_perm=n_perm)
        tn = spatial.spatial_corr(x.loc[keep], y.loc[keep], null="naive",
                                  n_perm=n_perm, seed=seed)
        rows.append(dict(component=comp, n_regions=len(keep), rho=t.r,
                         p_spin=t.p, p_naive=tn.p,
                         null_sd_spin=t.null_sd, null_sd_naive=tn.null_sd))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Gene-level ranking
# --------------------------------------------------------------------------

@dataclass
class GeneRanking:
    """A ranked gene list with the provenance needed to interpret it."""

    scores: pd.Series
    n_regions: int
    regions: tuple[str, ...]
    source: str

    def top(self, n: int = 100, tail: str = "both") -> pd.Index:
        s = self.scores.sort_values()
        if tail == "positive":
            return s.index[-n:]
        if tail == "negative":
            return s.index[:n]
        return s.index[:n // 2].append(s.index[-(n // 2):])


def rank_genes_directly(map_: pd.Series,
                        geom: spatial.ParcelGeometry,
                        expr: pd.DataFrame | None = None,
                        n_perm: int = 1000,
                        seed: int = 0,
                        covariate: pd.Series | None = None,
                        ) -> tuple[pd.DataFrame, dict]:
    """Correlate a cortical map with every gene, with a rotation-based omnibus.

    Returns ``(per_gene, omnibus)``.  ``per_gene`` has one row per gene with
    the Spearman correlation and a spin p-value.  ``omnibus`` tests whether
    the *set* of gene correlations is stronger than chance.

    Two things make the omnibus necessary rather than optional.  First, the
    per-gene p-values are not independent -- genes share expression gradients,
    so a map that happens to align with the dominant gradient produces
    thousands of "significant" genes from one underlying alignment.  Counting
    genes below 0.05 therefore overstates evidence badly.  Second, the null
    must preserve the map's spatial autocorrelation, or any smooth map will
    correlate with any smooth gene.  Both are handled by rotating the map on
    the sphere and recomputing the entire gene vector per rotation: the
    resulting null is over *maps*, so it absorbs the gene-gene dependence.

    ``covariate`` residualises both the map and every gene on a regional
    nuisance variable (parcel size, residual noise, baseline thickness) before
    correlating, which is how a confounded signal is distinguished from a real
    one.  Correlations are computed on ranks, so residualisation is on ranks
    too -- a rank-partial correlation, not a Spearman of residuals.
    """
    expr = ahba_expression() if expr is None else expr
    m = bilateral(map_)
    pos = {lab: i for i, lab in enumerate(geom.labels)}
    keep = [lab for lab in m.index.intersection(expr.index) if lab in pos]
    if len(keep) < 30:
        raise ValueError(
            f"only {len(keep)} regions shared between map, expression and "
            "geometry; the spin null is meaningless at this size"
        )
    idx = np.array([pos[lab] for lab in keep])
    spins = spatial.spin_indices(geom, n_perm=n_perm, seed=seed, hemi_paired=False)
    spins = np.clip(np.searchsorted(idx, spins[:, idx]), 0, len(idx) - 1)

    X = expr.loc[keep].to_numpy(float)
    Xr = np.apply_along_axis(stats.rankdata, 0, X)
    Xr = (Xr - Xr.mean(0)) / Xr.std(0)
    cov = None
    if covariate is not None:
        cov = stats.rankdata(bilateral(covariate).loc[keep].to_numpy(float))
        cov = (cov - cov.mean()) / cov.std()
        Xr = Xr - np.outer(cov, (Xr * cov[:, None]).mean(0))
        Xr = Xr / Xr.std(0)

    def gene_vector(values: np.ndarray) -> np.ndarray:
        v = stats.rankdata(values)
        v = (v - v.mean()) / v.std()
        if cov is not None:
            v = v - cov * float((v * cov).mean())
            sd = v.std()
            # A covariate collinear with the map leaves nothing to correlate.
            # Return exact zeros rather than 0/0 -> nan, so the caller sees
            # "no residual signal" instead of a silently poisoned null.
            if sd < 1e-12:
                return np.zeros(Xr.shape[1])
            v = v / sd
        return (Xr * v[:, None]).mean(0)

    x = m.loc[keep].to_numpy(float)
    r_obs = gene_vector(x)
    null = np.array([gene_vector(x[s]) for s in spins])

    p_gene = (1 + (np.abs(null) >= np.abs(r_obs)).sum(0)) / (1 + n_perm)
    per_gene = pd.DataFrame({"gene": expr.columns, "rho": r_obs, "p_spin": p_gene})

    # Omnibus 1: mean |rho| over all genes -- threshold-free, the primary test.
    mean_obs = float(np.abs(r_obs).mean())
    mean_null = np.abs(null).mean(1)
    p_mean = (1 + int((mean_null >= mean_obs).sum())) / (1 + n_perm)
    # Omnibus 2: how many genes exceed the null's own 95th percentile.
    cut = float(np.percentile(np.abs(null), 95))
    n_obs = int((np.abs(r_obs) >= cut).sum())
    n_null = (np.abs(null) >= cut).sum(1)
    p_count = (1 + int((n_null >= n_obs).sum())) / (1 + n_perm)
    omnibus = dict(
        n_regions=len(keep), n_genes=len(r_obs), n_perm=n_perm,
        mean_abs_rho=mean_obs, null_mean_abs_rho=float(mean_null.mean()),
        p_mean=p_mean, cut_95=cut, n_genes_exceed=n_obs,
        null_median_exceed=float(np.median(n_null)), p_count=p_count,
        covariate=None if covariate is None else (covariate.name or "unnamed"),
    )
    return per_gene, omnibus


def rank_genes_by_component(dev_map: pd.Series,
                            components: pd.DataFrame | None = None,
                            weights: pd.DataFrame | None = None,
                            component: str = "C3") -> GeneRanking:
    """Project a developmental map onto genes through a component's loadings.

    The map is correlated with the component's regional scores, then each
    gene inherits ``sign(rho) * loading`` on that component.  This is the
    cheap route: it assumes the map's gene relevance is fully captured by one
    pre-computed component.  :func:`rank_genes_directly` does not assume that
    and is preferred when a gene-by-region expression matrix is available.
    """
    comps = ahba_components() if components is None else components
    w = ahba_gene_weights() if weights is None else weights
    dev_b = bilateral(dev_map)
    common = dev_b.index.intersection(comps.index)
    rho = stats.spearmanr(dev_b.loc[common], comps.loc[common, component]).statistic
    return GeneRanking(
        scores=np.sign(rho) * w[component],
        n_regions=len(common), regions=tuple(common),
        source=f"AHBA {component} loadings, signed by map rho={rho:+.3f}",
    )


def concordance(a: pd.Series, b: pd.Series, n_perm: int = 10000,
                seed: int = 0) -> dict:
    """Rank concordance of two gene scorings on their common genes.

    Used to ask whether the AHBA-derived and snRNA-seq-derived gene rankings
    agree -- the project's hypothesis (ii).  The null permutes gene labels,
    which is adequate here because genes (unlike cortical regions) are not
    spatially autocorrelated in this index.
    """
    common = a.index.intersection(b.index)
    x, y = a.loc[common].to_numpy(float), b.loc[common].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    r = float(stats.spearmanr(x, y).statistic)
    rng = np.random.default_rng(seed)
    null = np.array([stats.spearmanr(rng.permutation(x), y).statistic
                     for _ in range(n_perm)])
    return dict(rho=r, n_genes=int(ok.sum()),
                p=(1 + int((np.abs(null) >= abs(r)).sum())) / (1 + n_perm))
