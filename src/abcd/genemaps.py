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
