"""
07_permutation_enrichment.py -- H2, test A: are SCZ / MDD GWAS gene sets enriched
in the transcriptomic thinning signature?  (Whitaker & Vertes 2016 style.)

For each gene-weight vector w (thinning-oriented: positive = higher expression
where adolescent thinning is faster) and each GWAS gene set S:
  * universe  = brain-expressed background for the vector's DS level
                (AHBA genes with a MAGMA Z in both SCZ and MDD) ∩ genes with a weight
  * statistic = mean w over S ∩ universe
  * null      = 10,000 random sets of the same size from the universe
                (a) unmatched, (b) matched on gene length decile (MAGMA START-STOP)
  * p         = two-sided empirical, plus one-sided "enriched in thinning direction"
  * effect    = z vs null; AUC (P[w_set > w_other], Mann-Whitney)
  * tails     = hypergeometric enrichment of S in the top and bottom deciles of w

Weight vectors
  lead ds0/ds25/ds50 : option 2 (dCT+CT) PLS2, -Z         [the signature]
  opt1/opt3/opt4 dCT components at ds25 (-Z)              [same signature, other Y]
  opt5 PLS2 (slopePC3-driven), oriented to correlate + with the lead
  opt2 PLS1 ds25 (-Z): the C1-like STATIC component        [within-analysis control]
  C1, C2, C3 (Dear 2024), NSPN PLS1_z, PLS2_z (Whitaker 2016) [benchmarks, native sign]

Output: results/permutation_enrichment.tsv (long), results/permutation_enrichment_wide.tsv
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES, REF = ROOT / "results", ROOT / "data" / "reference"
GS = REF / "gene_sets"
N_PERM, SEED = 10000, 0
rng = np.random.default_rng(SEED)

def load_Z(opt, ds, comp, sign=-1):
    t = pd.read_csv(RES / "pls_weights" / f"{opt}_{ds}.tsv", sep="\t", index_col=0)
    return sign * t[f"{comp}_Z"].dropna()

lead25 = load_Z("opt2_dCT_CT", "ds25", "PLS2")
opt5 = load_Z("opt5_slopePCs", "ds25", "PLS2", sign=1)
sh = opt5.index.intersection(lead25.index)
if stats.spearmanr(opt5.loc[sh], lead25.loc[sh]).statistic < 0:
    opt5 = -opt5
nspn = pd.read_csv(REF / "nspn_pls_gene_weights.csv").set_index("gene")
c123 = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)

vectors = {
    ("lead_opt2_PLS2", "ds0"): load_Z("opt2_dCT_CT", "ds0", "PLS2"),
    ("lead_opt2_PLS2", "ds25"): lead25,
    ("lead_opt2_PLS2", "ds50"): load_Z("opt2_dCT_CT", "ds50", "PLS2"),
    ("opt1_dCT_PLS1", "ds25"): load_Z("opt1_dCT", "ds25", "PLS1"),
    ("opt3_dCT_dT1T2_PLS2", "ds25"): load_Z("opt3_dCT_dT1T2", "ds25", "PLS2"),
    ("opt4_full4_PLS2", "ds25"): load_Z("opt4_full4", "ds25", "PLS2"),
    ("opt5_slopePCs_PLS2", "ds25"): opt5,
    ("control_opt2_PLS1_static", "ds25"): load_Z("opt2_dCT_CT", "ds25", "PLS1"),
    ("AHBA_C1", "ds0"): c123["C1"], ("AHBA_C2", "ds0"): c123["C2"], ("AHBA_C3", "ds0"): c123["C3"],
    ("NSPN_PLS1_z", "ds0"): nspn["PLS1_z"].dropna(), ("NSPN_PLS2_z", "ds0"): nspn["PLS2_z"].dropna(),
}

SETS = ["SCZ_prioritised", "SCZ_prio_finemap", "SCZ_locus_pool", "SCZ_pool_not_prio",
        "SCZ_trubetskoy_smr101", "MDD_highconf", "MDD_hc_finemap", "MDD_pool", "MDD_pool_not_hc",
        "MDD_howard2019_32", "MDD_howard2019_magma"]
sets = {s: set(l.strip() for l in open(GS / f"{s}.txt") if l.strip()) for s in SETS}
bg = {ds: set(l.strip() for l in open(GS / f"brain_background_{ds}.txt") if l.strip())
      for ds in ("ds0", "ds25", "ds50")}
# gene length from the SCZ MAGMA table (build 37)
mg = pd.read_csv(GS / "magma_SCZ_genes.tsv", sep="\t").dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
glen = (mg["STOP"] - mg["START"]).astype(float)

rows = []
for (vname, ds), w in vectors.items():
    univ = sorted(set(w.index) & bg[ds])
    wv = w.loc[univ]
    W = wv.to_numpy(); n = len(W)
    ranks = stats.rankdata(W)
    lbin = pd.qcut(glen.reindex(univ).fillna(glen.median()), 10, labels=False, duplicates="drop").to_numpy()
    bin_idx = {b: np.flatnonzero(lbin == b) for b in np.unique(lbin)}
    top = set(wv.nlargest(n // 10).index); bot = set(wv.nsmallest(n // 10).index)
    for sname, S in sets.items():
        Sg = [g for g in univ if g in S]
        k = len(Sg)
        if k < 5:
            continue
        idx = wv.index.get_indexer(Sg)
        obs = W[idx].mean()
        # unmatched null
        null = np.array([W[rng.choice(n, k, replace=False)].mean() for _ in range(N_PERM)])
        # length-matched null
        counts = pd.Series(lbin[idx]).value_counts()
        nullm = np.empty(N_PERM)
        for i in range(N_PERM):
            pick = np.concatenate([rng.choice(bin_idx[b], c, replace=False) for b, c in counts.items()])
            nullm[i] = W[pick].mean()
        def ptwo(nl): return (1 + (np.abs(nl - nl.mean()) >= abs(obs - nl.mean())).sum()) / (1 + N_PERM)
        def pone(nl): return (1 + (nl >= obs).sum()) / (1 + N_PERM)
        auc = (ranks[idx].sum() - k * (k + 1) / 2) / (k * (n - k))
        mw_p = stats.mannwhitneyu(W[idx], np.delete(W, idx), alternative="two-sided").pvalue
        nt, nb = len(top & set(Sg)), len(bot & set(Sg))
        exp = k * (n // 10) / n
        rows.append(dict(vector=vname, ds=ds, gene_set=sname, n_universe=n, n_set=k,
                         mean_Z_set=obs, mean_Z_universe=W.mean(),
                         z_unmatched=(obs - null.mean()) / null.std(), p2_unmatched=ptwo(null), p1_thinning_unmatched=pone(null),
                         z_lengthmatched=(obs - nullm.mean()) / nullm.std(), p2_lengthmatched=ptwo(nullm), p1_thinning_lengthmatched=pone(nullm),
                         auc=auc, p_mannwhitney=mw_p,
                         n_top_decile=nt, n_bottom_decile=nb, expected_decile=exp,
                         p_top_decile=stats.hypergeom.sf(nt - 1, n, n // 10, k),
                         p_bottom_decile=stats.hypergeom.sf(nb - 1, n, n // 10, k)))
    print(vname, ds, file=sys.stderr)

R = pd.DataFrame(rows)
R.to_csv(RES / "permutation_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
wide = R.pivot_table(index=["vector", "ds"], columns="gene_set", values="z_lengthmatched").round(2)
wide = wide[[s for s in SETS if s in wide.columns]]
wide.to_csv(RES / "permutation_enrichment_wide.tsv", sep="\t")
pd.set_option("display.width", 250)
print("z (length-matched null), positive = set genes have higher thinning-oriented weight")
print(wide.to_string())
print("\nlead ds25 detail:")
print(R[(R.vector == "lead_opt2_PLS2") & (R.ds == "ds25")][
    ["gene_set", "n_set", "z_lengthmatched", "p2_lengthmatched", "auc", "n_top_decile", "n_bottom_decile", "expected_decile", "p_top_decile"]]
    .round(4).to_string(index=False))
