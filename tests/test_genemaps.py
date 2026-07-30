import pytest



def test_expression_row_order_is_validated():
    import pytest
    """The row-order guard must reject a shuffled expression matrix."""
    import numpy as np
    from abcd import genemaps as gm
    E = gm.ahba_expression()
    scores = gm.ahba_components()
    shuffled = E.iloc[np.random.default_rng(0).permutation(len(E))].set_axis(E.index)
    with pytest.raises(ValueError, match="do not align"):
        gm._assert_row_order(shuffled, scores)


def test_rank_genes_directly_omnibus(hcp_geom_lh, dev_map_hcp):
    from abcd import genemaps as gm
    E = gm.ahba_expression()
    per_gene, om = gm.rank_genes_directly(dev_map_hcp, hcp_geom_lh, expr=E, n_perm=50)
    assert len(per_gene) == E.shape[1]
    assert {"gene", "rho", "p_spin"} <= set(per_gene.columns)
    assert per_gene.rho.abs().max() <= 1.0
    assert 0 < om["p_mean"] <= 1 and om["n_regions"] == 137
    # a covariate identical to the map must annihilate the signal
    _, om_self = gm.rank_genes_directly(dev_map_hcp, hcp_geom_lh, expr=E,
                                        n_perm=50, covariate=dev_map_hcp)
    assert om_self["mean_abs_rho"] < 1e-6


def test_ahba_components_dsk_is_bilateral_and_dk_sized():
    from abcd import genemaps as gm
    """The DK score file ships with an ``lh_`` prefix but bilateral content.

    Callers compare it against a :func:`bilateral`-reduced imaging map, whose
    index has no hemisphere prefix, so the loader must strip it -- otherwise
    the index intersection is empty and every correlation is silently NaN.
    """
    dsk = gm.ahba_components("dsk")
    assert list(dsk.columns) == ["C1", "C2", "C3"]
    assert len(dsk) == 34, "DK cortex has 34 bilateral regions"
    assert not dsk.index.str.match(r"^(lh|rh)_").any()
    assert "bankssts" in dsk.index


def test_ahba_components_dsk_aligns_with_bilateral_imaging_map():
    import pandas as pd
    from abcd import genemaps as gm
    labels = [f"{h}_{r}" for r in gm.ahba_components("dsk").index for h in ("lh", "rh")]
    fake = pd.Series(range(len(labels)), index=labels, dtype=float)
    assert len(gm.bilateral(fake).index.intersection(
        gm.ahba_components("dsk").index)) == 34


def test_ahba_components_rejects_unknown_parcellation():
    from abcd import genemaps as gm
    with pytest.raises(ValueError, match="unknown parcellation"):
        gm.ahba_components("fsaverage5")
