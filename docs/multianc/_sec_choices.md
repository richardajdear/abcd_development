### Is the answer sensitive to *k*?

**For the conclusions that drive decisions: no.**

| conclusion | k=3 | k=4 | k=5 | sensitive? |
|---|---|---|---|---|
| `--grm-cutoff` keeps ~75 % of EUR-like | 75.8 % | 75.6 % | 75.4 % | no |
| …and ≤33 % of every other stratum | 2.7, 16.7 | 2.6, 22.2, 1.7 | 0.0, 21.3, 33.4, 1.8 | no |
| unrelated set is overwhelmingly European | 93.6 % | 93.7 % | 90.8 % | no |
| PC-AiR prunes far more evenly | ✔ | ✔ | ✔ (one exception) | no |
| within − between off-diagonal gap | 0.077 | 0.075 | 0.070 | no |
| clusters recover self-reported groups | partly | ✔ | ✔ | **yes** |

The only conclusion that moves with *k* is the *interpretation* of the clusters,
not the methodological finding. That is the expected behaviour: *k* controls how
finely PC space is partitioned, while the `--grm-cutoff` failure is a property of
the GRM's allele frequencies, which do not depend on any partition at all. The
strata are a **lens for seeing** the failure, not its cause.

### Why we used k=4

- **k=3 is too coarse.** It merges the Hispanic/Latino and Asian-ancestry groups
  into one cluster (70 % Hispanic *and* 14 % Asian), so a "stratum" would contain
  two populations with different allele frequencies — the exact confound the
  stratification exists to remove.
- **k=4 gives one cluster per recognisable group**, each validated against
  self-report, and all 5,656 anchor-European children land in one cluster with
  0.0 % anchor leakage into any other — at every k from 3 to 5.
- **k=5 is arguably better on purity** and we say so: it isolates an admixed
  cluster and raises `cluster1` from 79.6 % to 94.1 % self-reported Black. We did
  not adopt it because the extra cluster (n=773) has only 221 unrelated children,
  below the `STRAT_MIN_N = 300` floor for stratified heritability, so k=5 buys
  purity we cannot spend. **If a future analysis needs cleaner group definitions
  rather than more strata to fit, k=5 is the better choice** and
  `strata_k5.tsv` is already written.
- **k=2 and k=6 were also generated** (`strata_k2.tsv`, `strata_k6.tsv`). k=2
  collapses everything non-European into one bin (58 % anchor-EUR in the European
  cluster — badly impure). k=6 begins splitting the European cluster itself
  (27 anchor-European children leak into `cluster2`), which is over-fitting.

### Choices where we picked a number, and why

| decision | value | basis | what would change it |
|---|---|---|---|
| PCs entering k-means | 4 | PC1–4 carry 6.04/1.45/0.41/0.13 % of variance; PC5 onward are <0.1 % and are sampling noise | a sample with finer structure |
| PC scaling | ×√eigenvalue | unscaled eigenvectors weight a 0.13 % PC like a 6.04 % one; produced a demonstrably wrong split (§3a) | nothing — this is a correctness fix |
| k | 4 | one cluster per self-reported group; k=3 merges two, k=5 adds an unusable stratum | needing purity over usable strata → k=5 |
| relatedness pruning | PC-AiR | §5; `--grm-cutoff` is not ancestry-neutral | — |
| `STRAT_MIN_N` | 300 unrelated | below this, REML SEs span the parameter space and a 0.000 invites being read as evidence of absence | — |
| "substantially uncertain" | mean dosage ambiguity ≥ 0.05 | **our definition, not a standard** — see §6 | it is a comparison device; the ranking is robust, the cut is not load-bearing |
| imputation quality | `R² ≥ 0.8` | conventional; and it costs only ~10 % of stratum-common variants, so it is not the binding filter | — |
| fileset frequency floor | `MAC ≥ 10` | the union rule expressed in one pass (§6) | — |
| GRM frequency floor | `--maf 0.01`, plus a 0.001 arm | 0.01 is conventional; 0.001 is required because the published comparison used it | — |

### The compromise we did NOT resolve

**The pooled GRM still uses pooled allele frequencies.** A single GRM cannot
simultaneously have ancestry-specific allele frequencies *and* cross-ancestry
relatedness entries — the stratified GRMs get the former by having no
cross-ancestry pairs at all. Methods that do both exist (PC-Relate computes
individual-specific allele frequencies from PCs; `FastSparseGRM` builds a
block-diagonal ancestry-adjusted GRM), and **standard GCTA does not implement
them**. Our mitigations are in-sample PCs as covariates and the stratified arm as
a cross-check; neither repairs the pooled GRM itself. This is the principal
known weakness of the pooled heritability estimate and should be stated wherever
that estimate is reported.
