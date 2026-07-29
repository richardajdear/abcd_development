# Adolescent cortical development in ABCD 7.0: findings so far

Status: modelling framework validated on release 7.0; genetic analyses ready to
dispatch. Supersedes `REPORT_5.1_legacy.md`, which documents the 5.1 draft.

## 1. The framework reproduces across releases

Group-level developmental maps (regional mean subject slope, Desikan-Killiany,
68 regions) correlate ρ = 0.977 between 5.1 and 7.0, and ρ = 0.999 between the
≥2- and ≥3-visit filters within 7.0. The release upgrade and the adapter rewrite
did not disturb the cortical pattern, so downstream differences are attributable
to sample size and follow-up length rather than to pipeline changes.

Mean whole-cortex thinning is −0.0016 mm/yr on 7.0 (−0.0019 on 5.1).

## 2. `min_visits: 3` costs power; use ≥2

Slope reliability per subject rises with follow-up length, but retaining fewer
subjects costs more than the reliability gain returns. Effective N — the
quantity GWAS power scales with — peaks at the permissive filter:

| run | subjects | mean visits | mean slope reliability | effective N | vs 5.1 |
|---|---|---|---|---|---|
| 5.1, ≥2 visits | 6,937 | 2.30 | 0.165 | 1,145 | 1.00 |
| **7.0, ≥2 visits** | **8,192** | **2.86** | **0.179** | **1,463** | **1.28** |
| 7.0, ≥3 visits | 5,195 | 3.35 | 0.206 | 1,068 | 0.93 |
| 7.0, 4 visits | 1,830 | 4.00 | 0.249 | 456 | 0.40 |

Effective N = subjects × mean reliability. The ≥3-visit filter gives *less*
power than 5.1 did (0.93×) despite the larger release; ≥2 visits on 7.0 buys a
28% gain. Note that 7.0's headline advantage is longer follow-up (2.86 vs 2.30
mean visits) as much as more subjects.

`configs/ct_70_2visit.yaml` is therefore the default for phenotype work.
Figure: `release_power_tradeoff.png`.

## 3. The family random effect destroys the genetic signal — the key finding

**Any run intended for heritability or GWAS must set `family_effect: false`.**

Fitting `(1 | family_id)` partitions the between-family variance into its own
random effect. The subject-level BLUPs that remain are *within-family
deviations*, and the between-family component — exactly what genetic
relatedness explains — has been removed. With two phenotyped members per
family, the family effect centres them at zero, so their deviations must sum to
zero and anti-correlate by construction.

Measured on 7.0, whole-cortex thickness, 624 same-sex twin pairs and 578
sibling pairs:

| phenotype | model | r (same-sex twins) | r (siblings) | Falconer h² |
|---|---|---|---|---|
| baseline thickness *(positive control)* | with `(1\|family_id)` | +0.125 | **−0.117** | 0.485 |
| baseline thickness | without | +0.544 | +0.340 | 0.409 |
| developmental slope *(target)* | with `(1\|family_id)` | +0.130 | −0.009 | 0.278 |
| developmental slope | without | +0.304 | +0.162 | 0.283 |

Cortical thickness is strongly familial, so the negative sibling correlation is
not a weak result but a structurally impossible one — which is what exposed the
problem. Note the trap: Falconer's estimator differences two correlations and so
partly cancels the artefact, returning a *plausible-looking* h² = 0.485 from
invalid inputs. The h² value does not reveal the bug; only the pair correlations
do. `gcta_export.py` therefore guards on the model specification
(`FamilyEffectConflict`), not on the h² value, with regression tests in
`tests/test_gcta_export.py`.

Relatedness must instead be handled where it belongs: in the GRM for GCTA, or
via a mixed-model GWAS.

Figure: `family_effect_heritability.png`. Config: `configs/ct_70_genetic.yaml`.

### Headline result

**The developmental slope phenotype is heritable: h² ≈ 0.28** (twin/sibling
design, whole-cortex thickness). This is the number that makes the imaging
genetics of longitudinal change viable, and it is roughly two-thirds of the
baseline-thickness estimate (0.41) from the same subjects and estimator.

Caveat: zygosity is design-based, not genotype-confirmed. `twin_same_sex` mixes
MZ and DZ pairs, which *deflates* r_twin and so makes h² a conservative
underestimate rather than an inflated one. Genotype-based zygosity and a GCTA
GRM estimate are the next step; both need the cluster.

## 4. Developmental maps do not resemble AHBA C1–C3

Testing the group developmental map against the AHBA transcriptional components
in DK space, using the 5,000 spin permutations from the AHBA repo (rotating the
component, holding the imaging map fixed), bilateral averaging to match the
left-hemisphere component scores:

| run | C1 | C2 | C3 |
|---|---|---|---|
| 5.1, ≥2 visits | ρ=+0.024, p=0.65 | ρ=+0.011, p=0.68 | ρ=−0.038, p=0.62 |
| 7.0, ≥2 visits | ρ=+0.126, p=0.47 | ρ=−0.063, p=0.61 | ρ=−0.143, p=0.44 |
| 7.0, ≥3 visits | ρ=+0.132, p=0.46 | ρ=−0.055, p=0.62 | ρ=−0.140, p=0.44 |

All p > 0.44. This replicates the negative result the 5.1 draft recorded, now on
better data, so it is unlikely to be a power problem at the map level.

**Interpretation.** The version of the hypothesis that predicts the *group mean
developmental map* resembles C3 is not supported. This does not test the
hypothesis that actually motivates the project. C3 is a map of *where* genes are
expressed; the genetic hypothesis is about *which* genes carry variance in
developmental change between individuals. Those come apart: a gene set can drive
individual differences in maturation without its spatial expression gradient
matching the group mean rate of change. The subject-level heritable phenotype in
§3 is the right substrate, and the informative test is whether GWAS hits for the
slope phenotype are enriched in the C3 gene loadings and in snRNA-seq PC1 — not
whether the two cortical maps align.

I would treat §4 as a documented negative control on the map-level approach and
put the weight on §3 → GWAS → gene-level enrichment.

## 5. Next steps

1. **On the cluster** (see `hpc/README.md`): build the GRM, run GCTA-REML on
   `intercept` and `slope` for a genotype-based h², then fastGWA. Ancestry PCs
   are *not* in `covar_quant.txt` — add them first; ABCD is multi-ancestry.
2. Genotype-based zygosity to replace the design-based twin classes and
   sharpen the Falconer estimate.
3. Gene-level enrichment of slope GWAS against C3 loadings, snRNA-seq PC1, and
   SCZ/MDD GWAS (deferred by agreement until the GWAS exists).
4. Regional GWAS (68 phenotypes) once the whole-cortex run validates.
5. Additional imaging metrics — the assembly layer is metric-agnostic, so
   surface area, and functional measures, need only a config.

## Reproducing

```bash
export ABCD_ROOT=/path/to/ABCD
python -m abcd.assemble configs/ct_70_genetic.yaml     # -> out/<run_id>/
Rscript R/fit_lmm.R --run-dir out/<run_id> --cores 8
python -m abcd.phenotype out/<run_id>
python -m abcd.gcta_export out/<run_id>                # guarded
```
