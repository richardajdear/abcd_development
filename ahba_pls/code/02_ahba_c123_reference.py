"""02_ahba_c123_reference.py — AHBA C1-C3 (Dear et al. 2024) gene weights and DK scores, plus recomputed scores.

Inputs: ../data/weights.csv (gene x C1,C2,C3; 7,973 genes), ../data/ahba_dme_dsk_scores.csv (34 lh_ DK regions x C1..C3),
        AHBA_updated ahba_dk_lh_native_ds{0,25,50}.csv (33 LH regions x genes).
Recompute: z-score each gene across the 33 regions, restrict to genes shared with weights, score = Z @ w / ||w|| (per component).
Outputs: ahba_c123_gene_weights.csv, ahba_c123_dk_scores.csv, ahba_c123_scores_recomputed_ds{0,25,50}.csv, ahba_c123_recompute_check.json
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr

D = Path('/Users/richard/Git/abcd_development/data')
AHBA_UPD = Path('/Users/richard/Git/AHBA_updated/outputs/expression_levels')
OUT = Path('/Users/richard/Git/abcd_development/ahba_pls/data/reference')

w = pd.read_csv(D / 'weights.csv', index_col=0); w.index.name = 'gene'
assert w.shape == (7973, 3) and w.index.is_unique
w.to_csv(OUT / 'ahba_c123_gene_weights.csv')
sc = pd.read_csv(D / 'ahba_dme_dsk_scores.csv')
assert sc.columns[0] == 'label' and sc.label.str.startswith('lh_').all(), sc.columns.tolist()
sc = sc.set_index('label').sort_index()
sc.to_csv(OUT / 'ahba_c123_dk_scores.csv')

res = {'weights_shape': list(w.shape), 'scores_shape': list(sc.shape),
       'scores_regions_not_in_ds0': None, 'ds': {}}
for ds in (0, 25, 50):
    X = pd.read_csv(AHBA_UPD / f'ahba_dk_lh_native_ds{ds}.csv', index_col='label')
    assert X.shape[0] == 33
    Z = (X - X.mean()) / X.std(ddof=1)
    shared = Z.columns.intersection(w.index)
    W = w.loc[shared]; W = W / np.linalg.norm(W, axis=0)
    S = pd.DataFrame(Z[shared].values @ W.values, index=Z.index, columns=w.columns)
    S.index.name = 'label'
    S.to_csv(OUT / f'ahba_c123_scores_recomputed_ds{ds}.csv')
    common = S.index.intersection(sc.index)
    if res['scores_regions_not_in_ds0'] is None:
        res['scores_regions_not_in_ds0'] = sorted(set(sc.index) - set(S.index))
    r = {'n_shared_genes': int(len(shared)), 'n_regions': int(len(common))}
    for c in w.columns:
        r[f'{c}_spearman'] = round(float(spearmanr(S.loc[common, c], sc.loc[common, c])[0]), 3)
        r[f'{c}_pearson'] = round(float(pearsonr(S.loc[common, c], sc.loc[common, c])[0]), 3)
    res['ds'][ds] = r
json.dump(res, open(OUT / 'ahba_c123_recompute_check.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
