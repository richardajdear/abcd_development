# Retired analyses

Kept for provenance, not part of the reported results.

## `prscs/` — PRS-CS with the 1000 Genomes EUR LD reference
Superseded.  Its LD reference is 503 individuals; the UKB panel used by the
reported PRS-CS run is ~50,000.  This is not a detail: the 1000G run was the
single outlier that made the SCZ polygenic association look absent
(README_HPC §13.5), and §13.8 records that SBayesR on the same sumstats gave
p = 7.3e-04.  The reported four methods all use UKB LD.

## `prscsx/` — PRS-CSx, multi-ancestry
Dropped because it could not do the job it was brought in for.  PRS-CSx needs
per-ancestry sumstats, and the PGC SCZ `afram` and `latino` files ship 11
columns with no NCAS/NCON/NEFF and no daner equivalent, so the two ancestries
our target actually needs could not be included.  Running it on EUR+EAS only
addresses the wrong problem: the inflation in the pooled arm is target-side
ancestry confounding, which a discovery-side multi-ancestry prior does not
touch.  The within-ancestry standardisation is the correction that does.

## `prs_ct_v3/`, `assoc/`, `within_ancestry/`
Earlier C+T and within-ancestry runs.  Superseded by `../prs_final/`, which
recomputes every cell of the 4-method x 6-trait-arm grid from one set of
normalised sumstats, with the discovery GWAS matched to the target arm.
