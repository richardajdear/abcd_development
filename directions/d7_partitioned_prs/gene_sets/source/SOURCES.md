# Gene-set sources

- `syngo1.3_genes.tsv`: SynGO release 1.3 (Koopmans et al. 2019, *Neuron* 103:217, doi:10.1016/j.neuron.2019.05.002), `genes.xlsx` of https://syngoportal.org/data/syngo1.3_complete_data.zip, downloaded 2026-09-30; columns hgnc_id, hgnc_symbol, ensembl_id, entrez_id (1,789 genes).
- Seidlitz et al. 2020 cell-class markers: `~/Git/AHBA/data/seidlitz_cell_genes.csv` (the file `ahba_pls/code/13_celltype_compare.py` uses).
- Everything else is already in the repo: `data/weights.csv` (AHBA C1-C3), `data/velmeshev_PC1_gene_loadings.csv`, `ahba_pls/results/hcp_summary_base_gene_weights.tsv` (ABCD PLS2, HCP base), `ahba_pls/data/reference/gene_sets/SCZ_locus_pool.txt`.
