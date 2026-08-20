### 6.1 "Substantially uncertain" — a definition we invented, and why

Imputation replaces a measured genotype with a **dosage**: a number in [0, 2]
giving the expected count of the alternate allele. A confident call sits near 0,
1 or 2; an uncertain one sits between.

The natural quality metric is Minimac's **R²** — its own estimate of how well
each variant was imputed. **But R² cannot answer our question.** Minimac writes
*one R² per variant, estimated across the whole cohort*. There is no per-ancestry
value, and it cannot be recomputed after the fact because it is an internal
quantity of the imputation, not a function of its output. **The number that would
reveal ancestry bias in an R² filter is the one number R² does not provide.**

So we measure what R² summarises. For each variant, within each stratum:

```
ambiguity = mean over children in that stratum of  min(|DS−0|, |DS−1|, |DS−2|)
```

- **0** — every child called confidently
- **0.5** — maximally uncertain (all dosages at 0.5 or 1.5)

and we call a variant **"substantially uncertain" in a stratum when its ambiguity
there is ≥ 0.05** — on average each child's dosage sits ≥0.05 from a whole number.

**This threshold is ours, not a standard.** It was chosen as "visibly off a clean
call" and is a *comparison device*: what matters is the contrast between strata,
not the absolute rate. The ranking is robust to the cut; the specific 0.05 is
not load-bearing, and no downstream filter uses it. It connects to something real
though — PLINK's hard-call threshold is 0.1, so dosages beyond that become
**missing** rather than rounded, and a stratum with higher ambiguity loses more
genotypes at conversion.

### 6.2 The measurement that had to be redone

Our first pass averaged ambiguity over **all** sampled variants and returned
0.0003–0.0007 for every stratum — apparently perfect, apparently identical.
**That conclusion was an artefact of the statistic.** 97.5 % of TOPMed variants
are near-monomorphic in this cohort, so nearly every child has dosage ≈0 and
ambiguity ≈0 *by construction*. The average was measuring how rare the variants
are, not how well they were imputed; it would have returned ≈0.0004 even if every
common variant were imputed badly.

The corrected version restricts to **variants common in the stratum being
assessed** — the only ones whose quality can affect that stratum's GRM entries.
This is recorded because the flawed version looked like a clean pass, and a clean
pass is exactly what nobody re-examines.

### 6.3 Results

RESULTS_PLACEHOLDER
