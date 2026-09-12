# data/reference — reference signatures and gene sets (Phase 1, "Reference signatures" track)

## NOTES
- Three re-runnable scripts in `code/`: `01_nspn_reference.py` + `01b_nspn_symbol_harmonise.py` (NSPN), `02_ahba_c123_reference.py`
  (AHBA C1–C3), `03_gwas_gene_sets.py` (SCZ/MDD sets, MAGMA Z, backgrounds). mygene.info was reachable; responses are cached as JSON
  so re-runs are offline-deterministic.
- Detailed notes: `NSPN_REFERENCE.md`, `AHBA_C123_REFERENCE.md`, `gene_sets/GENE_SETS.md`.
- The AHBA DK matrices cover 33 of 34 LH regions; the missing region is `lh_frontalpole` (also in `../missing_region.txt`, written by the inputs track).
- Region labels are ggseg-style (`lh_bankssts`, …); all label-indexed files use the column name `label`.
