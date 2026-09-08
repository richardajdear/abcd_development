# hpc_v3 — regional-subset phenotypes, and the order to run the experiments

**Written 2026-09-08.** This directory now holds two things:

1. **A new phenotype arm (Experiment A, this README):** genetics of the mean
   thinning rate over (i) the fastest-thinning regions and (ii) the highest
   AHBA-C3 regions, run through the existing v2 GENESIS pipeline.
2. **[`SETUP_CONTEXT.md`](SETUP_CONTEXT.md):** the method context for the
   power-borrowing route (MTAG / genomic SEM; MOSTest and JAGWAS ruled out).
   Read it before touching Experiments B–D below. Nothing here supersedes its
   findings; §5's "required first steps" are *reordered* by §4 of this file,
   with reasons.

Files:

| file | what |
|:---|:---|
| `make_phenotypes_v3.py` | builds the phenotypes + GCTA export (run locally) |
| `phenotypes_v3_regions.csv` | the two region sets with ΔCT and C3 values |
| `phenotypes_v3_characterization.csv` | 9×9 phenotype correlation matrix |
| `phenotypes_v3_consistency.csv` | lh/rh internal consistency + k-sensitivity |
| `h2_vs_scz_relevance.csv`, `phenotypic_corr_real.csv`, `regional_slope_eigen.csv` | evidence tables cited by SETUP_CONTEXT |

---

## 1. The new phenotypes

Both are built from the settled genetic run
(`ct_70_noglobal_mv2_genetic` → `out/thickness_dsk_70_139406217085`), from the
same per-subject regional `slope` values (BLUP + fixed) that `global_slope`
averages. Selection uses **group-level maps only** — no subject-level variance,
no genotypes — so there is no selection circularity for h², PRS, or rg.

| name | construction | role |
|:---|:---|:---|
| `slope_topDelta` | mean slope over the **8 bilateral regions** (16 labels) with the most negative group-mean thinning rate (`slope_total` in `docs/developmental_maps_noglobal.csv`, hemispheres averaged) | target |
| `slope_topC3` | mean slope over the 8 bilateral regions with the highest AHBA C3 score (`data/ahba_dme_dsk_scores.csv`) | target |
| `slope_projDelta` | z-scored slope map · mean-centred group thinning map | exploratory |
| `slope_projC3` | z-scored slope map · mean-centred C3 map | exploratory |

**Sign convention that matters:** high C3 = association cortex = *faster*
thinning — ρ(group slope, C3) = −0.546 across the 34 bilateral regions. So the
two sets overlap: **3 of 8 regions shared** (parsorbitalis, inferiorparietal,
rostralmiddlefrontal). `slope_topDelta` adds frontal pole, precuneus,
paracentral, medial orbitofrontal, superior parietal; `slope_topC3` adds
bankssts, supramarginal, pars triangularis, lateral orbitofrontal, inferior
temporal. Full table: `phenotypes_v3_regions.csv`.

### What the characterization says — read before interpreting any result

From `phenotypes_v3_characterization.csv` and `phenotypes_v3_consistency.csv`
(8,192 subjects):

- **The subset means are largely the same trait as `global_slope`:**
  r = 0.89 (`slope_topDelta`) and 0.92 (`slope_topC3`); the two subset means
  correlate 0.83 with each other. Experiment A is therefore an
  **SNR-concentration test** — does restricting to high-signal cortex sharpen
  the disorder association? — not a test of a new biological axis. Any
  write-up that treats `slope_topC3` and `global_slope` as independent
  findings is double-counting, exactly as SETUP_CONTEXT §4(b) warns for
  `slope_PC1`.
- **k-sensitivity is monotone and offers no free lunch:** at k = 4 bilateral
  regions r with global falls to 0.80/0.85 but internal consistency drops too
  (Spearman-Brown 0.73/0.75); at k = 17 consistency matches global (0.83) but
  r = 0.96 — indistinguishable from `global_slope`. k = 8 is the compromise:
  SB 0.82/0.79 vs 0.87 for global.
- **The projections are the genuinely different traits** (r with global 0.28
  and 0.41) and load on `slope_PC3` (0.52, 0.63) — the phenotype with the
  strongest AHBA coupling (C2: +0.854) but *negative* LDSC h² (SETUP_CONTEXT
  §7). They are also the noisiest (SB 0.63/0.65). They ride along at near-zero
  marginal cost; treat them as exploratory.
- **Success criterion is the disorder association, not h²** — the
  h²-vs-relevance anti-ranking (SETUP_CONTEXT §4c) is the reason this project
  does not optimise heritability. Concretely: the primary readout is the
  SCZ/MDD PRS association for the new phenotypes **compared against
  `global_slope` in the same subjects**. Because the estimates are correlated
  (r ≈ 0.9), test the *difference* properly: bootstrap subjects (resample
  families, not individuals) for the paired Δβ, or a seemingly-unrelated
  regression. A bigger point estimate alone is not a result.

### Rebuilding the export

```bash
export ABCD_CONFIG=ct_70_noglobal_mv2_genetic
PYTHONPATH=src python hpc_v3/make_phenotypes_v3.py
```

Writes `out/<run>/gcta_inputs_v3/`: `phenotypes_gcta.txt` (FID IID + 9
phenotype columns: the settled five, then the four above),
`covar_quant.txt` (baseline_age, n_visits, 10 ancestry PCs),
`covar_categorical.txt` (sex, site), `phenotype_manifest.tsv`, and
`phenotype_manifest_new_only.tsv` (just the four new rows — the v2 sbatch
scripts size their arrays from the manifest, so pointing `MANIFEST` at this
file runs only the new phenotypes while `--mpheno` still indexes the 9-column
file). Per-subject files stay under gitignored `out/` — the data-use rule in
SETUP_CONTEXT §8 applies.

---

## 2. Experiment A — the HPC run

Everything reuses `hpc_v2` (GENESIS, relatedness-aware) plus the v1 scoring
steps that consume v2 sumstats (`hpc/04_magma.sbatch`, `hpc/05_ldsc_rg.sbatch`,
`hpc/06_prs.sbatch`). **No new pipeline code should be needed**; if a step
needs modification, document it here per hpc_v2/README_HPC.md §10.

On CSD3, with `hpc_v2/config.local.sh` sourced:

| step | script | array | notes |
|:---|:---|:---|:---|
| 0 | copy `gcta_inputs_v3/` → `$ABCD_HPC_ROOT/pheno_v3/` | — | then `export PHENO=$ABCD_HPC_ROOT/pheno_v3/phenotypes_gcta.txt MANIFEST=$ABCD_HPC_ROOT/pheno_v3/phenotype_manifest_new_only.tsv COVAR_QUANT=... COVAR_CAT=...` |
| 1 | *(skip)* `01_gds`, `02_kinship` | — | **reuse** `work/results_v2/kinship/` (kinship_sparse.rds, pcair_pcs.tsv, pcair_unrelated.txt) — do not recompute (SETUP_CONTEXT §6) |
| 2 | `hpc_v2/03_null_model.sbatch` | 1–4 | one GENESIS null model per new phenotype |
| 3 | `hpc_v2/04_assoc.sbatch` | 1–88 | 4 phenotypes × 22 chr; sumstats out as `assoc/<pheno>.sumstats.tsv.gz` |
| 4 | `hpc_v2/05_reml_zaitlen.sbatch` | 1–4 | Zaitlen two-GRM h²; **report both kinship sparsity thresholds** — global_slope swung 0.183→0.115 between them (SETUP_CONTEXT §7) |
| 5 | LDSC: h² per new phenotype + rg vs SCZ, MDD | — | `hpc/05_ldsc_rg` pattern. Expectation from the global_slope precedent: h² z will likely be < 4 and rg **uninformative rather than null** — say so if it is; that is why PRS is the primary readout |
| 6 | PRS association: SCZ, MDD, **+ ASD and Alzheimer's controls** | — | `hpc/06_prs` pattern (population, EUR and all-ancestry strata) **and** `hpc_v2/06_prs_family.sbatch` (within-family Fulker decomposition). The AD control is mandatory: a non-specific confound is present in the population arm (SETUP_CONTEXT §7). Multiplicity across C+T thresholds: family-level permutation of min-p, not Bonferroni — or a continuous-shrinkage score |
| 7 | MAGMA: gene-based + gene-property, both directions | — | `hpc/04_magma.sbatch` pattern, including the `_vs_` gene-covar reverse regression (disorder gene Z on phenotype gene Z) |

Report λ_GC for every new GWAS — a hit count without it is not interpretable
(v1's 20 "hits" at λ = 1.107 were stratification).

### What to report (definition of done for A)

1. Per new phenotype: λ_GC, GENESIS h² (both GRM thresholds), LDSC h² (+z).
2. SCZ and MDD PRS β with ASD/AD controls, population + within-family, with
   the **paired Δβ vs `global_slope`** and a permutation-corrected p.
3. LDSC rg vs SCZ/MDD if and only if h² z supports it; otherwise the z and the
   statement that rg is uninformative at this N.
4. MAGMA gene-property both directions for SCZ/MDD.
5. Four new rows appended to the h²-vs-relevance table
   (`h2_vs_scz_relevance.csv`) — the anti-ranking pattern gains its first
   out-of-sample test: the subset phenotypes were *chosen* for biological
   relevance, so if the pattern holds they should sit low-h²/high-association.
6. An explicit statement of what the result licenses: these phenotypes share
   ~80–84 % of variance with `global_slope`; only the paired difference is a
   finding.

---

## 3. Experiments B–D — the power-borrowing route (from SETUP_CONTEXT)

- **B — LDSC phenotype × phenotype rg matrix + pilot regional GWAS
  (the go/no-go).** Needs sumstats for the 9 phenotypes (5 exist from v2;
  4 arrive with A) plus a pilot ~10 regional slopes (10 × 22 assoc tasks).
  Gate: if per-trait h² z is too low for stable rg, **MTAG and genomic SEM are
  blocked** — write that down and stop rather than proceeding on unstable
  covariances (SETUP_CONTEXT §5.1).
- **C — 68 per-region GWAS** (68 × 22 tasks, the big compute). Only if B
  passes. Reuse v2 kinship; report per-region λ_GC.
- **D — MTAG pilot → full; genomic SEM after.** Pilot: `global_slope` primary
  + 10 regional slopes; check `maxFDR` and whether SE(global_slope) actually
  falls. Scale to 68 only if it does. Genomic SEM last, and only if MTAG works
  and the phenotype needs reshaping. All four caveats in SETUP_CONTEXT §3
  apply, especially: MTAG's effective N is not a sample size, and nothing
  disorder-derived may enter the boost.

## 4. Recommended order: A → B → (gate) → C → D

**Run the phenotype tests (A) before any power-borrowing experiment.** Four
reasons:

1. **Cost asymmetry.** A adds 4 phenotypes to a pipeline that already runs —
   roughly the size of v2's original 5-phenotype pass, reusing kinship and
   null-model machinery. C alone is ~14× v2's association volume, and D sits
   behind a gate that may close.
2. **Gating.** The MTAG route's precondition (per-trait LDSC h² z, guidance
   ~z > 4) is *already doubtful*: `global_slope` sits at z = 1.16 and single
   regions will be noisier than their own average. A is not gated by anything.
3. **A's results redesign D.** If `slope_topC3` (or a projection) carries a
   stronger SCZ association than `global_slope`, *it* becomes the MTAG target
   — boosting the wrong primary trait wastes C's compute. And B needs A's
   sumstats anyway for the full 9×9 rg matrix.
4. **Hypothesis-directness.** The top-C3 phenotype is the direct test of the
   project's transcriptional hypothesis; MTAG/SEM are precision engineering.
   If A shows the C3-selected cortex carries no disorder signal beyond the
   global mean, that result reshapes what is worth boosting at all — and it
   is a publishable negative on its own terms.

One deliberate redundancy: B's rg matrix includes the four new phenotypes, so
A also buys down B's cost. The only scenario where A-first loses time is if
the subset phenotypes were pure noise — ruled out locally: internal
consistency 0.79–0.82, and r ≈ 0.9 with a phenotype whose PRS association is
already established.

## 5. Rules carried over (do not re-derive)

- Reuse v2 kinship/PC-AiR; never the v1 pooled GRM (SETUP_CONTEXT §6).
- Join IDs on the normalised NDAR token; family IDs from the phenotype
  export, never the `.fam` (§7).
- `scratch/hpc_test*` are synthetic fixtures with planted effects — never
  data (§7). `make_phenotypes_v3.py` asserts the 8,192 × 68 shape for this
  reason.
- Per-subject genotype-derived files are never committed (§8).
