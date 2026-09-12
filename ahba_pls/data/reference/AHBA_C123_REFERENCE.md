# AHBA C1–C3 (Dear et al. 2024, Nat Neurosci) reference signatures

## NOTES
- Built by `code/02_ahba_c123_reference.py` (re-runnable).
- `ahba_c123_gene_weights.csv` = verbatim copy of `abcd_development/data/weights.csv` (7,973 genes × C1,C2,C3; index renamed `gene`).
- `ahba_c123_dk_scores.csv` = verbatim copy of `abcd_development/data/ahba_dme_dsk_scores.csv` (34 `lh_` DK regions × C1..C3,
  label column already ggseg-style). It includes `lh_frontalpole`, which the DK expression matrices lack (33 regions).
- Recomputed scores (`ahba_c123_scores_recomputed_ds{0,25,50}.csv`, 33 × 3): each gene z-scored across the 33 regions,
  restricted to genes shared with the weights, score = Z_shared @ w_shared with each weight column L2-normalised.
  The shipped scores were derived on a different pipeline (HCP-parcellated / earlier abagen build projected onto DK), so they are
  not expected to be identical; concordance on the 33 shared regions:

| DS level | n shared genes | C1 ρ | C2 ρ | C3 ρ | Pearson r C1 / C2 / C3 |
|---|---|---|---|---|---|
| ds0 | 7,971 | 0.991 | 0.944 | 0.906 | 0.991 / 0.969 / 0.891 |
| ds25 | 7,873 | 0.991 | 0.944 | 0.906 | 0.991 / 0.969 / 0.894 |
| ds50 | 6,969 | 0.99 | 0.952 | 0.907 | 0.99 / 0.971 / 0.898 |

- C3 (the target component) reproduces at Spearman ρ ≈ 0.91 at every DS level; DS filtering barely changes the recomputed scores
  because dropped genes carry small weights.
- **Use the recomputed scores** for H1 DK-score concordance with ABCD-PLS (identical regions/genes/z-scoring); the shipped scores
  are kept as the published reference.
