# Directions memo — what mechanism are we proposing, and what would test it?

*2026-09-28. A planning document, not a results document. Every number quoted
here is from `README.md` / `genetic_analysis/current_results.tsv` as of this
date; nothing below has been run.*

## 1. The problem with the current story

The draft arc is: (fig 1) individual thinning rates track SCZ/MDD polygenic
risk and later symptoms; (fig 2) the group thinning map follows AHBA C3, which
is enriched for SCZ/MDD; (fig 3) the snRNA-seq maturation axis aligns with C3.

Read as a mechanism this reduces to three facts that do not chain together:

1. disorder risk genes are expressed in neurons (known);
2. polygenic risk shifts the global thinning rate by ~0.03 SD/SD (our result);
3. neuron-rich (C3-high) cortex thins fastest in the group mean (our result).

Our own tests say the chain breaks at the individual level. The PRS effect is
a uniform whole-cortex shift (regional PRS β maps do not resemble C3 or PLS2;
the row-centred C3 axis carries no SCZ/MDD signal). MAGMA with C3 loadings as a
gene property is null on the slope GWAS. Symptom-linked thinning follows the
normative map, not C3. So fig 2 currently explains *where normative thinning
happens*, not *why risk changes its rate*. The paper should say this plainly.

## 2. What the results do support: a developmental-tempo claim

Take the specificity panel as a whole. SCZ, MDD and ALZ-with-APOE push thinning
faster; **educational attainment pushes it slower, at the same magnitude**;
none touches baseline thickness. That is the signature of polygenic influence
on the *pace* of adolescent cortical maturation rather than on cortical
integrity. Related literature: protracted cortical development with higher
cognitive outcome (Shaw et al. 2006); accelerated maturation under adversity
(Tooley, Bassett & Mackey 2021); adolescence as a critical period whose closure
can be mistimed (Larsen & Luna 2018).

Under this reading, C3 and snRNA-seq PC1 are the transcriptional definition of
the *programme whose timing is shifted*: they mark which cortex is still
maturing in the 9–16 window and which genes constitute the maturation
programme. It is then coherent that they predict the group map and carry GWAS
signal without predicting individual differences in rate. The proposed
mechanism is: **adolescent-onset disorders are disorders of maturational
tempo, and the transcriptional programme being paced is the one that carries
their risk.**

### 2.1 What makes this falsifiable

"Tempo" is a claim about a *shared* latent process, so it makes predictions
outside the thickness measure that generated it. Each is a test that can fail:

| # | Prediction | Falsified if |
|:--|:--|:--|
| T1 | The PRS effect on the thinning slope depends on age/puberty (a phase advance shifts a curved trajectory, so the effect at 10 differs from the effect at 15; a pure rate scaling gives the same multiple at every age) | PRS × age and PRS × puberty are flat **and** T2 fails |
| T2 | The same scores shift other maturational slopes from the same scans in the "more mature" direction: T1w/T2w rising faster, WM FA rising faster, surface area falling faster, puberty earlier | the effect is thickness-specific |
| T3 | EA is opposite in sign on every index in T2, as it is for thickness | EA is opposite for thickness only |
| T4 | PRS-related extra thinning is proportional to the normative map (it happens where maturation is happening) | the PRS β map is a flat offset unrelated to dCT under a spin null. *Current evidence leans this way (uniform shift), but the per-parcel maps are noisy; test with the projected phenotype* |
| T5 | Faster thinners show smaller longitudinal cognitive gains, given baseline (premature closure of plasticity) | gains are unrelated or larger |
| T6 | Within MZ pairs, the faster-thinning twin has more later symptoms (slope → symptom is not explained by shared genes/environment) | the MZ-difference association is zero |

If T1–T3 fail the effect is thickness-specific and the tempo reading is wrong;
if T5–T6 fail the thinning is a marker without consequence. Either outcome is
publishable, which is the point.

### 2.2 Pruning or myelination — the two accounts and their imaging signatures

Adolescent "thinning" in T1 is at least two tissue processes: loss of
synaptic/dendritic neuropil, and intracortical myelination of deep layers that
moves the apparent GM/WM boundary outward (Natu et al. 2019; Whitaker et al.
2016). They predict different co-changes:

| measure (same scans) | myelination account | neuropil-loss account |
|:--|:--|:--|
| T1w/T2w ratio (or T1 GM intensity) in the thinning parcel | rises with thinning, within child | no co-change |
| grey/white contrast at the boundary | falls (GM becomes WM-like) | unchanged or rises |
| WM volume / adjacent WM FA | rises | no co-change |
| RSI restricted-directional (neurite density proxy) in cortex | unchanged | falls |
| regional pattern of the thinning | follows the myelin axis (AHBA C1 / T1w/T2w map; our slope PC2) | follows the neuronal axis (C3 / PLS2) |
| gene partition carrying the PRS effect | oligodendrocyte / myelin genes | synaptic (SynGO), C3-positive pole, complement |

Our between-child slope decomposition already separates the two candidate
axes: PC2 is the C1/myelin axis, and the C3 axis appears only after row
centring. The question becomes: on which axis does the PRS act, and which
co-changes accompany the risk-related thinning? C4A imputation being null is a
mild strike against complement-mediated pruning as the risk-carrying route;
it says nothing about myelination.

## 3. Directions, with what each buys, what it needs, and where it runs

### D1. Multimodal slopes from the same children (pace, and myelin vs neuropil)

*Buys:* T2, T3, T4 and the §2.2 discrimination. The single analysis that turns
"correlates with C3" into "we can name the tissue process the risk acts on".

*Needs:* per-scan, per-parcel measures beyond thickness.

- **Available now (DK, vendored):** T1 and T2 grey-matter intensity
  (`mr_y_smri__t1__gm__dsk`, `t2__gm__dsk`; the `t1t2_70_*` config already
  builds the ratio). So DK-level dCT vs dT1w/T2w, and PRS → T1w/T2w slope, can
  start immediately on the laptop.
- **In the release, not vendored:** DK area and volume (`mr_y_smri__area__dsk`,
  `vol__dsk`), DTI FA/MD (`mr_y_dti__*`), RSI (`mr_y_rsi__*`). Fetch from the
  release; same pipeline (`assemble` → `fit_lmm` → `phenotype`).
- **HCP-MMP:** thickness is vendored; area and volume were extracted by
  `hcp_stats.py` on CSD3 (`mr_y_smri__{thk,area,vol}__hcp.tsv`) but only `thk`
  is here — copy the other two. T1w/T2w in HCP-MMP is **not** in the release
  and is not produced by `mris_anatomical_stats`; it needs the T2 registered to
  the T1, intensity sampled at mid-depth (`mri_vol2surf`), and parcellated —
  a new CSD3 job over the FreeSurfer tree. Grey/white contrast
  (`?h.w-g.pct.mgh`) is a FreeSurfer output and can be parcellated with
  `mri_segstats` in the same job. Do the DK version first; only build HCP-MMP
  if the DK result needs the resolution.

*Analysis:* (i) fit slopes per modality with the settled specification;
(ii) within-child, across parcels, partial correlation of dCT with dT1w/T2w
and dGWC; (iii) PRS panel (SCZ, MDD, EA, ALZ±APOE) → each modality slope;
(iv) common-factor model across modality slopes, PRS → factor vs residual.

### D2. Timing versus rate: one-stage PRS × age, PRS × puberty

*Buys:* T1. *Needs:* nothing new — scan-level table, PRS export (8,596),
`ph_y_pds`. `tools/age_prs_interaction.R` is the template (README_HPC §4.2).
Add puberty tempo (PDS slope, or age at PDS stage) as an outcome. **Runs
locally now.**

### D3. Twin and sibling designs

*Buys:* T6 and the genetic-correlation questions LDSC cannot answer.
*Needs:* `gn_y_genrel` (vendored; zygosity codes 1/2/3, ~840 / ~1,290 / ~1,550
children — confirm code meanings against the data dictionary), the slope
phenotypes, and the outcomes. **Runs locally** (OpenMx or umx in R).

Why this is different from what has been done: GREML/LDSC/GWAS all estimate
*SNP-based* quantities and are power-limited by h²_SNP ≈ 0.18 at n ≈ 4,300
(LDSC h² z ≈ 0.4, so rg with anything is uninformative). The twin design does
not depend on SNP tagging and uses the full sample regardless of ancestry.
Three things become possible:

1. **Bivariate ACE** gives the *genetic correlation* between the thinning
   slope and another trait — puberty tempo, the T1w/T2w slope, CBCL change,
   cognitive gain — from cross-twin cross-trait covariances. This is the
   "shared aetiology" question we currently cannot ask at all.
2. **Within-family PRS** (regress sibling differences in slope on sibling
   differences in PRS) removes population stratification, assortative mating
   and genetic nurture. It matters most for the EA result, which in other
   cohorts shrinks by about half within family (Selzam 2019; Howe 2022). The
   informative quantity is the within/between ratio, not the within-family
   p-value — at β ≈ 0.03 the within-family test alone is underpowered.
3. **MZ-difference design** (T6): does the faster-thinning twin of an MZ pair
   have more symptoms at 15–17? All genetic and shared-environmental
   confounding is removed. No PRS design can do this; it is the closest thing
   to a causal test of slope → symptom that observational data allow.

### D4. Thinning rate → longitudinal cognitive gain

*Buys:* T5. *Needs:* NIH Toolbox tables (`nc_y_nihtb`; in the release, not
vendored; full battery at baseline and year 2, reduced at year 4). Model: gain
~ slope + baseline score + age + sex + site + SES. **Runs locally once the
table is fetched.** Cheap, and it is the claim a reader will remember.

### D5. Gene × environment and direction of effect

*Buys:* whether PRS and adversity act on the same pace variable (additive or
interactive), and whether the slope → symptom path is directional.
*Needs:* release tables for adversity/SES (demographics, ADI, family conflict,
life events — `abcd_p_demo`, `led_l_adi`, `fes`, `mh_y_le` families; not
vendored), plus CBCL at every wave (vendored). Cross-lagged model: baseline
symptoms → slope vs slope → later symptoms. Also the sensible covariate check
for D2/D3: the EA effect with SES in the model. **Local once fetched.**

### D6. Make fig 3 temporal: which cell-type programmes move between ~8 and ~18?

*Buys:* a cell-type-resolved prediction for §2.2 (oligodendrocyte-dominated
adolescent change → myelination; synaptic/excitatory-dominated → neuropil),
and disorder enrichment in the *adolescent-window* genes rather than the whole
axis. Lives in `transcriptional_maturation`, not here.

*The identifiability problem:* one time point per donor, few donors, so
age-related change is confounded with between-donor variation. The way round
it is that **every donor contributes every cell type**. Donor-level nuisance
(batch, PMI, sex, ancestry, agonal state) is shared across a donor's cell
types, so a cell-type × age effect is identifiable by contrasting cell types
*within* donor even when the marginal age effect is not: pseudobulk per donor
× cell type, fit `expression ~ cell_type × s(age) + (1 | donor)`, and read the
cell-type-specific age derivative in the 8–18 window. Equivalently, score each
donor's cells on the existing maturation PC1 per cell type and ask, per cell
type, where the score's age curve is still steep versus saturated. Replicate
the window across the two datasets (Herring, Velmeshev/U01), and borrow the
window definition from bulk time-courses with far more donors (BrainSpan,
PsychENCODE developmental) before trusting it in the snRNA-seq.

### D7. Partitioned polygenic scores

Already planned (README_HPC §4.4, §8.3 C3-D). Run as a *directional* check on
D1 and D6 — synaptic vs oligodendrocyte vs C3-pole partitions against random
partitions of equal SNP count — not as the route to the mechanism. With
β_total ≈ −0.03 each partition carries ~−0.01; only the contrast is testable.
HPC.

### Not recommended

More region-selection phenotypes (closed, README_HPC §3); more GWAS-side work
at this n; further attempts to make the regional PRS map resemble C3.

## 4. Suggested order

1. D2 (days, local, nothing to fetch) — establishes whether we are looking at
   timing or rate before anything else is built.
2. D1 on DK with the vendored T1/T2 tables (days, local) — the first
   myelination-vs-neuropil readout; decide from it whether HCP-MMP T1w/T2w
   and RSI are worth a CSD3 job.
3. Fetch NIH Toolbox, area/volume, DTI/RSI, adversity tables in one pass
   (needs release access) → D4, D5, remainder of D1.
4. D3 twin models on whatever slopes exist by then (local, R).
5. D6 in the snRNA-seq repo, in parallel.
6. D7 on HPC when the C3-axis pipeline (§8) is otherwise finished.

## 5. Revised arc if the tempo reading holds

- **Fig 1** (unchanged): individual thinning rate tracks polygenic risk in
  both directions (SCZ/MDD faster, EA slower), not baseline thickness, and
  tracks later symptoms.
- **Fig 2 — tempo:** PRS × age/puberty, multimodal slopes, twin genetic
  correlations; the process is mostly myelin / mostly neuropil / both.
- **Fig 3 — programme:** the maturation programme across adolescence
  (snRNA-seq, cell-type resolved), its adult spatial footprint (C3), its match
  to the normative thinning map, its disorder enrichment — stated as
  group-level, with the null individual-level C3 tests reported.
- **Fig 4 — consequence:** faster thinning → smaller cognitive gain, more
  depressive symptoms; MZ-difference test.
