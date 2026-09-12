# Follow-up on the cluster: H3 and H4

Specification for the agent running on CSD3. Nothing here has been run. Phases 1–3
(the PLS derivation and the SCZ/MDD enrichment) are complete and local — see
`README.md`. This document defines the two arms that need ABCD genetics.

Read `../hpc_v2/README_HPC.md` §4 (setup), §7 (gotchas) and `../hpc_v3/README_HPC.md`
§1 (the export and its three silent mismatches) before submitting anything. The
conventions below are theirs, not new ones.

## What exists to work from

| file | what |
|:--|:--|
| `hpc/lead_pls2_dk_scores_68.csv` | the lead component's regional scores (gene-side, thinning-oriented), 68 DK labels; `lh_frontalpole`/`rh_frontalpole` are **NA** (no AHBA coverage) |
| `hpc/lead_pls2_gene_covar_entrez.txt` | MAGMA `--gene-covar` file, GENE = Entrez (NCBI37.3, same ids as `hpc/work/results/magma/*.genes.raw`): `thinning_Z_ds0/ds25/ds50` + `AHBA_C3` |
| `hpc/lead_pls2_gene_weights_symbol.tsv` | the same weights keyed by gene symbol, with rank and decile |
| `results/lead_signature_scores.csv` | 33-region scores with the reference maps alongside |

Orientation is fixed throughout: **positive weight / positive score = higher
expression where adolescent thinning is faster.** Do not re-sign anything.

## H3 — are the signature genes the genes driving global thinning?

**Test.** MAGMA gene-property regression in the *reverse* direction used by
`hpc_v3/README_HPC.md` §10: the ABCD phenotype's own gene-level Z regressed on the
signature weights.

1. Take the phenotype `.genes.raw` from the EUR arm (`hpc_v2/work/results_v2/magma_eur_v3/`,
   or regenerate for a phenotype that has no gene analysis yet — the gene analysis is the
   expensive step, so reuse where possible).
2. `magma --gene-results <pheno>.genes.raw --gene-covar ahba_pls/hpc/lead_pls2_gene_covar_entrez.txt --out <...>`
   Run the three DS columns as separate marginal models **and** one model with
   `--model condition=AHBA_C3`, which is the test that matters: does the ABCD-derived
   signature add anything over C3 in the *ABCD's own* genetics?
3. Phenotypes: `global_slope` (primary), `baseline_thickness`, `slope_PC1–3`,
   and the v3 set (`slope_topDelta`, `slope_topC3`, `slope_projDelta`, `slope_projC3`).

**Prior, stated so it cannot be forgotten.** Every forward and reverse MAGMA test in
`hpc_v3` §10 was null, including the phenotype built from the highest-C3 regions
(C3+ on `slope_topC3`: p = 0.905). The expectation here is a null result; the value
of running it is that the signature is a *continuous, data-driven* weighting derived
from this cohort's own thinning map, which none of the earlier covariates were.
Report it as null if it is null — do not hunt for the one cell below 0.05 among ~30 tests.

## H4 — a PLS-derived phenotype for GWAS

**Construction** (follow `hpc_v3/make_phenotypes_v3.py`, which already implements the
projection idiom and its characterisation):

```
slope_projPLS2[i] = z(slope_BLUP[i, regions]) · centre(score[regions])
```

with `score` from `hpc/lead_pls2_dk_scores_68.csv`, the 66 covered labels only, and
the same two variants as v3: raw, and residualised on `global_slope`. Selection uses
group-level maps only (the PLS was fitted on group-mean maps, no subject-level
variance and no genotypes), so there is no selection circularity — **but** note that
the PLS *was* fitted on the same subjects' group means, which `slope_projC3` was not.
Run the split-half optimism check in `docs/topH2_selection_optimism.csv` style: refit
the PLS on half the subjects' maps, project the other half, compare h².

**Characterise before submitting** (v3 §1 conventions): correlation with `global_slope`,
`slope_PC1–3`, `slope_projC3`, `slope_projDelta`; lh/rh Spearman–Brown internal
consistency; and the phenotype's distribution. Add to `gcta_inputs_v*/phenotypes_gcta.txt`
plus `phenotype_manifest.tsv`.

**Then**, in the pipeline's own order: GCTA REML (Zaitlen two-GRM, PC-Relate second
component — §11.9 of `hpc_v2/README_HPC.md`), GENESIS association scan, LDSC, and MAGMA
in both directions.

**Benchmarks to beat.** `global_slope` h² = 0.115 ± 0.056; `slope_projC3` h² = 0.121 ± 0.057
(v3 REML, Zaitlen/PC-Relate arm). A projection phenotype whose h² lands inside that
interval is not an improvement, and `hpc_v3` §7 found the v3 projection phenotypes were
*less* heritable by LDSC. State the comparison with its standard errors; a point estimate
alone is not an answer.

## Reporting back

Append a dated section to `README.md` §Results index here, and cross-reference from
`hpc_v*/README_HPC.md` as those files do for each other. Record the job ids and the
measured wall time (`hpc_v2` §11.11 keeps a ledger for sizing).
