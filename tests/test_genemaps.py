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
