# Directions memo — what mechanism are we proposing, and what would test it?

*2026-09-28, revised 2026-09-30 (§§2.3–2.4, T7–T8, D7–D8, §6). A planning document, not a results document. Every number quoted
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
| T1 | The PRS effect on the thinning slope depends on age/puberty (a phase advance shifts a curved trajectory, so the effect at 10 differs from the effect at 15; a pure rate scaling gives the same multiple at every age) | PRS × age and PRS × puberty are flat **and** T2 fails. **Tested 2026-09-30, D2: flat.** The normative rate peaks at 12–14, yet the PRS rate effect is same-signed before and after (change after 13: SCZ +0.16 µm/yr/SD, p 0.57; MDD +0.05, p 0.86), there is no level effect at 12.8 (bounds a phase advance to < ~1 month per SD), and PRS × pubertal stage is null given PRS × age. Rate scaling, not a timing shift — see [`directions/d2_prs_age_puberty/`](../directions/d2_prs_age_puberty/README.md). |
| T2 | The same scores shift other maturational slopes from the same scans in the "more mature" direction: T1w/T2w rising faster, WM FA rising faster, surface area falling faster, puberty earlier | the effect is thickness-specific |
| T3 | EA is opposite in sign on every index in T2, as it is for thickness | EA is opposite for thickness only |
| T4 | PRS-related extra thinning is proportional to the normative map (it happens where maturation is happening) | the PRS β map is a flat offset unrelated to dCT under a spin null. *Current evidence leans this way (uniform shift), but the per-parcel maps are noisy; test with the projected phenotype* |
| T5 | Faster thinners show smaller longitudinal cognitive gains, given baseline (premature closure of plasticity) | gains are unrelated or larger |
| T6 | Within MZ pairs, the faster-thinning twin has more later symptoms (slope → symptom is not explained by shared genes/environment) | the MZ-difference association is zero |
| T7 | The SCZ and MDD effects are carried by their *shared* genetic component (transdiagnostic liability): both attenuate in a joint model, the shared/subtracted scores carry the effect, and the EA effect is carried by its cognitive (Cog) rather than non-cognitive (NonCog) part | SCZ-not-MDD and MDD-not-SCZ scores each carry an independent effect, or EA acts through NonCog — then the slope indexes something other than one cognitive-neurodevelopmental axis, and the tempo claim must be stated per disorder |
| T8 | Faster thinning predicts the age-appropriate SCZ-spectrum phenotype, psychotic-like experiences (PQ-BC), as it predicts depressive symptoms | the slope predicts depressive symptoms only — then the SCZ PRS result is a genetic correlate with no SCZ-relevant phenotypic readout in this window |

T1 has failed on its own: within 9–17 the polygenic effect is a uniform scaling of the rate, so the *phase-advance* version of tempo is out and what remains to test is the *amplitude* version (more cortex lost per year across the whole window) via T2–T3. If T2–T3 also fail the effect is thickness-specific and the tempo reading is wrong;
if T5–T6 fail the thinning is a marker without consequence; T7 decides whether
the claim is written as transdiagnostic or disorder-specific. Either outcome of
each is publishable, which is the point.

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

### 2.3 The critical-period account of schizophrenia: how established is each link?

Larsen & Luna (2018, *Neurosci Biobehav Rev*, §"Schizophrenia") synthesise the
pruning and thinning literatures as: persistently immature cortical molecular
state → reduced perineuronal nets (PNNs) → prolonged plasticity → excess
pruning → excess thinning → connectivity and cognitive deficits. The links are
not equally supported, and the paper should lean only on the ones that are.

| link | evidence | status |
|:--|:--|:--|
| "Persistently immature" molecular state (D1/D2, NR1/NR3A, PV, GABA-A α1/α2 ratios at pre-adolescent levels; Catts et al. 2013) | adult, chronic, medicated post-mortem tissue; the "immature" reading is one interpretation of marker ratios that also have activity-dependent explanations (e.g. PV downregulation) | hypothesis |
| PNN reduction in SCZ (Mauney 2013; Enwright 2016; Berretta 2015) | small post-mortem series (~10–20 per group; amygdala, PFC, entorhinal) from two labs | replicated, low n; cause vs consequence unknown |
| PNN → plasticity (Pizzorusso 2002; Carulli 2010) | causal rodent work: chondroitinase reopens the visual-cortex critical period | established in rodent sensory cortex; PFC/human by analogy |
| prolonged plasticity → *more* pruning | largely conceptual (Feinberg 1982 lineage). Plasticity is not pruning: adolescent pruning is activity-dependent and microglia/complement-mediated. Genetic support for a pruning route is C4 (Sekar 2016) and SCZ-iPSC microglia eliminating more synapses (Sellgren 2019), neither of which involves PNNs | weakest link |
| fewer synapses in SCZ | post-mortem L3 DLPFC spine loss replicated (Glantz & Lewis 2000; Konopaske 2014); in vivo SV2A PET reduced in frontal/ACC with d ≈ 0.8–0.9 in chronic patients (Onwordi 2020) and in antipsychotic-naïve first-episode patients (Onwordi 2023) | established for adult patients; whether via excess *adolescent* pruning has never been measured — no human developmental synaptic-density series exists |
| synapse loss → MRI thinning | the gap Howes & Onwordi (2023) name. Mouse: spine density tracked the VBM signal but cortical thickness did not change (Keifer 2015); cellular parameters explain ~36 % of GMV variance and spine plasticity did not thicken cortex (Asan 2021); spines are a fraction of a percent of GM volume. Human: developmental thinning tracks myelination (Natu 2019) and the expression profiles of CA1-pyramidal, astrocyte and microglia marker genes, with the sign reversed in ageing (Vidal-Piñeiro 2020; Parker 2020) | not established; the imaging–histology bridge is the weak point |
| excess thinning at psychosis onset (Cannon 2015, NAPLS; ENIGMA case-control maps) | steeper PFC thinning in converters, correlated with pro-inflammatory cytokines | established as MRI phenomenon; cellular basis inferred |

The solid parts are the two ends — adult patients have fewer synaptic markers;
PNNs gate plasticity in rodents — and the developmental middle is inference.
That is where a longitudinal multimodal design (D1) can contribute and thickness
alone cannot. Two further consequences: (i) the PNN/GABA-immaturity part of the
account is SCZ-specific and our data do not speak to it; the synaptic-pace part
is shared with depression (SV2A is reduced in MDD too, Holmes 2019; stress-
induced dendritic atrophy is a core MDD model) and our data do; (ii) the
imaging–histology bridge is not made by swapping thickness for another
morphometric — see §2.4.

### 2.4 Should we measure grey-matter "density" rather than thickness?

Howes & Onwordi note that Keifer et al. saw spine-density change in VBM signal
but not in thickness, and suggest density is the closer proxy. Two cautions.
VBM "density" is a modulated tissue-probability map, not a histological
density: it mixes thickness, folding and partial-volume/intensity effects at the
GM/WM boundary, and Keifer's VBM was post-mortem mouse tissue at ~1 mm cortical
thickness. Whatever it indexed there is plausibly the *intensity* component —
which T1w/T2w and grey/white contrast measure more directly. And it inherits
the same interpretive ambiguity as thickness (Asan 2021).

FreeSurfer does not produce VBM. What it produces, and what is nearly free:

- **volume** per parcel (= thickness × area): DK volume is a release table;
  HCP-MMP volume is already extracted on CSD3 by `hcp_stats.py`, not yet
  copied. Volume and area slopes decompose thinning into a boundary shift vs
  a surface change;
- **grey/white contrast** (`?h.w-g.pct.mgh`): a FreeSurfer output, parcellated
  with `mri_segstats` in one CSD3 pass over the existing tree;
- **T1 and T2 GM intensity** per DK parcel: vendored now.

True VBM (CAT12/SPM on ~34k raw T1s) is feasible on CSD3 but is a new pipeline
with its own scanner-harmonisation problems and answers a less specific
question than the intensity measures. Decision: volume + GWC + T1w/T2w slopes
first (D1); VBM only if those are suggestive or a reviewer asks for the
Keifer-style measure.

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

*Analysis:* (i) fit slopes per modality with the settled specification,
including volume and area so thinning decomposes into boundary shift vs
surface change (§2.4); (ii) within-child, across parcels, partial correlation
of dCT with dT1w/T2w and dGWC; (iii) PRS panel (SCZ, MDD, EA, ALZ±APOE) → each
modality slope; (iv) common-factor model across modality slopes, PRS → factor
vs residual. The factor score is also the **tempo composite**: its reliability
should exceed the thickness slope's (~0.2), which raises every downstream β
(D2–D5, D7) more than any change on the genetic side can.

### D2. Timing versus rate: one-stage PRS × age, PRS × puberty

*Buys:* T1. *Needs:* nothing new — scan-level table, PRS export (8,596),
`ph_y_pds`. `tools/age_prs_interaction.R` is the template (README_HPC §4.2).
Add puberty tempo (PDS slope, or age at PDS stage) as an outcome. **Runs
locally now.**

**Done 2026-09-30 — [`directions/d2_prs_age_puberty/`](../directions/d2_prs_age_puberty/README.md).**
One-stage LMM on 26,597 scans / 8,596 children: SCZ PRS −0.28 µm/yr per SD
(p 7 × 10⁻⁴), MDD −0.17 (p 0.04), no level effect at 12.8; no age dependence
(hinge-at-13 and quadratic terms null; segment effects same-signed at <12 and
≥14); PRS × pubertal stage null given PRS × age. MDD PRS predicts earlier
puberty (+0.055 SD stage at 12.8, p 8 × 10⁻⁷) and earlier puberty predicts faster
thinning (−0.068 SD/SD), but adjusting for puberty leaves the PRS → slope effect
unchanged. EA/ALZ/ASD scores were not local; the README gives the rsync.

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

*Status 2026-09-30 — done, [`directions/d4_cognitive_gain/`](../directions/d4_cognitive_gain/README.md).*
Sign consistent with T5, but only for crystallised measures (SES-adjusted β +0.024 SD per SD
slower thinning, driven by picture vocabulary); fluid composite and all executive tasks null.
Robust to SES, image quality, scoring and interval. Baseline thickness's association with gain,
unlike the slope's, is mostly SES. Leading alternative reading: a shared EA-type influence on
both thinning rate and vocabulary growth. Needs the EA score and the D3 MZ-difference test.

### D5. Gene × environment and direction of effect

*Buys:* whether PRS and adversity act on the same pace variable (additive or
interactive), and whether the slope → symptom path is directional.
*Needs:* release tables for adversity/SES (demographics, ADI, family conflict,
life events — `abcd_p_demo`, `led_l_adi`, `fes`, `mh_y_le` families; not
vendored), plus CBCL at every wave (vendored). Cross-lagged model: baseline
symptoms → slope vs slope → later symptoms. Also the sensible covariate check
for D2/D3: the EA effect with SES in the model. **Local once fetched.**

*Status 2026-09-30 — run except the EA step:*
[`directions/d5_adversity_direction/`](../directions/d5_adversity_direction/README.md).
SES/ADI relate to the rate only between sites, and within site SES predicts
thickness level, not rate. Parent-reported life events predict faster thinning
(−0.034, robust to site and scan quality). SES + adversity explain 4 % of the SCZ
polygenic effect, and there is no PRS × environment term. Between children,
thinning predicts later symptoms and baseline symptoms do not predict thinning.
Within-child cross-lags are small in both directions and need an LCM-SR check. The
EA-with-SES test awaits the CSD3 EA score file.

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

### D7. Partitioned polygenic scores from the single-cell programme

*Implemented as a CSD3 pipeline in [`directions/d7_partitioned_prs/`](../directions/d7_partitioned_prs/README.md) (SBayesRC posterior partitioned by gene; pre-registered sets and readout there). Not yet run.*

*Buys:* the one test that closes the single-cell → genetics → imaging triangle
at the individual level. The informative result is a partition whose
association with the slope (or symptoms) is *stronger* than the full score's,
relative to size-matched random partitions — that would say we have found the
molecular pathways that are actually predictive, not just enriched.

*Power arithmetic, so the design is honest about it.* The SE of a PRS → slope
β depends on n and the outcome, not on the score: ≈ 0.011 in the pooled arm
for any partition. If the thinning-relevant genetic signal is a component of
disorder liability, and a partition captures a share *s* of that signal while
carrying a fraction *f* of the disorder's SNP variance, then
β_partition ≈ β_full × s/√f:

| s (share of thinning signal) | f (share of disorder h²) | β_partition / β_full | detectable vs full score? |
|:--|:--|:--|:--|
| 1.0 | 0.25 | 2.0 (0.03 → 0.06) | yes |
| 0.75 | 0.25 | 1.5 (0.03 → 0.045) | marginal (~1.5 SE; the scores are correlated) |
| 0.5 | 0.25 | 1.0 | no gain |

So only **coarse partitions with strong enrichment** are testable: most of the
thinning signal has to sit in a minority of the genome. Fine-grained pathway
scans will return noise and must not be run as discovery.

*Design:*

1. Pre-specify three or four contrasts: brain-expressed vs not (positive
   control; must enrich); neuronal vs glial; the adolescent-window maturation
   genes from D6 vs the rest of the axis; synaptic (SynGO) vs
   oligodendrocyte/myelin (the §2.2 discrimination). The ST12 locus genes and
   the C3 poles (README_HPC §4.4, §8.3 C3-D) stay as secondary contrasts.
2. Partition scores with PRS-CS or SBayesRC weights restricted to each
   partition and its complement (PRSet, or SBayesRC annotation-partitioned
   weights).
3. Null: random partitions matched on SNP count **and** on share of disorder
   h². Report the enrichment ratio β_partition / β_expected with a CI — one
   interpretable number per contrast, not a p-value sweep.
4. Run the same partitions against **symptoms** in the full ~11k (no imaging
   needed; the MDD-score → depressive-symptom effect is somewhat larger than
   the slope effect), and against the D1 tempo composite once it exists.

*Where the power actually is.* The single-cell → genetics link is powered by
the disorder GWAS (n ≈ 10⁵): are adolescent-window, cell-type-resolved
programme genes enriched for SCZ/MDD heritability (S-LDSC / MAGMA)? That is
already partly done with C3 weights (`ahba_pls/`, disorder panel) and should be
redone with the D6 gene sets. The defensible architecture is then: single-cell
names the programme; the disorder GWAS shows the programme carries risk; ABCD
shows polygenic risk shifts the tempo of the process the programme executes;
the partitioned PRS is the one individual-level check, pre-registered, and
reported as underpowered if that is what it is. HPC.

### D8. Specificity: schizophrenia, depression, and the direction of EA

*Buys:* T7 and T8 — whether the paper is written as transdiagnostic or
disorder-specific. SCZ–MDD rg ≈ 0.3–0.35 and individual PRS correlation
≈ 0.1–0.15, so equal βs are consistent with a shared component *or* two
separate routes. *Needs:* existing scores; public summary statistics for the
derived GWAS; PQ-BC (in the release, not vendored). **Local**, apart from
scoring the new summary statistics, which follows the existing PRS pipeline.

1. **Joint model** (README_HPC §4.1): SCZ + MDD + EA + ALZ-noAPOE + APOE ε4 in
   one regression on the slope. Mutual attenuation → shared; independence →
   separate.
2. **Subtracted and shared scores.** GWAS-by-subtraction (Demange 2021) or
   mtCOJO for SCZ-conditional-on-MDD and MDD-conditional-on-SCZ; a shared
   score from a cross-disorder GWAS (PGC-CDG2 2019) or a genomic-SEM factor
   GWAS (Grotzinger 2022: psychotic vs internalising factors). Score each and
   test on the slope and the tempo composite. This directly asks "is it the
   shared part".
3. **EA split into Cog and NonCog** (Demange 2021; summary statistics public).
   Cog is negatively genetically correlated with SCZ, NonCog positively. If
   slower thinning tracks Cog, the triad SCZ/MDD/EA collapses onto one
   cognitive-neurodevelopmental axis; if NonCog, it is something else
   (personality/SES-like), and the SES analysis in D5 becomes the priority.
4. **Psychotic-like experiences** (PQ-BC, yearly) as the age-appropriate
   SCZ-spectrum outcome: slope → PLEs at 15–17 given baseline, alongside the
   depressive-symptom result (T8).

*Framing.* SCZ+, MDD+, ALZ/APOE+, EA− with nothing on baseline thickness is the
pattern of a **transdiagnostic liability**, not a disorder mechanism. It agrees
with ENIGMA's shared pyramidal-cell/dendritic thinning signature across six
disorders (Patel et al. 2021), and synaptic loss is itself not SCZ-specific
(§2.3). Write the paper that way, and let items 1–3 say how much, if anything,
is disorder-specific.

### Not recommended

More region-selection phenotypes (closed, README_HPC §3); more GWAS-side work
at this n; further attempts to make the regional PRS map resemble C3.

## 4. Suggested order

1. D2 (days, local, nothing to fetch) — establishes whether we are looking at
   timing or rate before anything else is built. D8 item 1 (joint model) runs
   in the same session on the same export.
2. D1 on DK with the vendored T1/T2 tables (days, local) — the first
   myelination-vs-neuropil readout; decide from it whether HCP-MMP T1w/T2w,
   GWC and RSI are worth a CSD3 job.
3. Fetch NIH Toolbox, area/volume, DTI/RSI, adversity, PQ-BC tables in one
   pass (needs release access) → D4, D5, D8 item 4, remainder of D1 including
   the tempo composite.
4. D8 items 2–3: obtain Cog/NonCog, subtracted and cross-disorder summary
   statistics, score with the existing PRS pipeline (HPC), test locally.
5. D3 twin models on whatever slopes exist by then (local, R).
6. D6 in the snRNA-seq repo, in parallel; its gene sets feed D7.
7. D7 on HPC, pre-registered contrasts only, once D6 gene sets and the tempo
   composite exist.

## 5. Revised arc if the tempo reading holds

- **Fig 1** (unchanged): individual thinning rate tracks polygenic risk in
  both directions (SCZ/MDD faster, EA slower), not baseline thickness, and
  tracks later symptoms.
- **Fig 2 — tempo:** PRS × age/puberty, multimodal slopes and the tempo
  composite, twin genetic correlations; the process is mostly myelin / mostly
  neuropil / both.
- **Fig 3 — programme:** the maturation programme across adolescence
  (snRNA-seq, cell-type resolved), its adult spatial footprint (C3), its match
  to the normative thinning map, its disorder enrichment in the disorder
  GWAS, and the pre-registered partitioned-PRS check — stated as group-level
  where it is, with the null individual-level C3 tests reported.
- **Fig 4 — consequence and specificity:** faster thinning → smaller
  cognitive gain, more depressive symptoms and PLEs; MZ-difference test;
  shared vs disorder-specific scores and Cog/NonCog.

## 6. References cited in this memo

Resolved against CrossRef on 2026-09-30 (title match confirmed for every entry); the DOI is the identifier to use.

- Asan 2021 — *Scientific Reports* 11 (2021). doi:10.1038/s41598-021-83491-8
- Berretta 2015 — *Schizophrenia Research* 167:18 (2015). doi:10.1016/j.schres.2014.12.040
- Cannon 2015 — *Biological Psychiatry* 77:147 (2015). doi:10.1016/j.biopsych.2014.05.023
- Carulli 2010 — *Brain* 133:2331 (2010). doi:10.1093/brain/awq145
- Catts 2013 — *Frontiers in Cellular Neuroscience* 7 (2013). doi:10.3389/fncel.2013.00060
- Demange 2021 — *Nature Genetics* 53:35 (2021). doi:10.1038/s41588-020-00754-2
- Enwright 2016 — *Neuropsychopharmacology* 41:2206 (2016). doi:10.1038/npp.2016.24
- Feinberg 1982 — *Journal of Psychiatric Research* 17:319 (1982). doi:10.1016/0022-3956(82)90038-3
- Glantz & Lewis 2000 — *Archives of General Psychiatry* 57:65 (2000). doi:10.1001/archpsyc.57.1.65
- Grotzinger 2022 — *Nature Genetics* 54:548 (2022). doi:10.1038/s41588-022-01057-4
- Holmes 2019 — *Nature Communications* 10 (2019). doi:10.1038/s41467-019-09562-7
- Howes & Onwordi 2023 — *Molecular Psychiatry* 28:1843 (2023). doi:10.1038/s41380-023-02043-w
- Keifer 2015 — *Nature Communications* 6 (2015). doi:10.1038/ncomms8582
- Konopaske 2014 — *JAMA Psychiatry* 71:1323 (2014). doi:10.1001/jamapsychiatry.2014.1582
- Larsen & Luna 2018 — *Neuroscience & Biobehavioral Reviews* 94:179 (2018). doi:10.1016/j.neubiorev.2018.09.005
- Mauney 2013 — *Biological Psychiatry* 74:427 (2013). doi:10.1016/j.biopsych.2013.05.007
- Natu 2019 — *PNAS* 116:20750 (2019). doi:10.1073/pnas.1904931116
- Onwordi 2020 — *Nature Communications* 11 (2020). doi:10.1038/s41467-019-14122-0
- Onwordi 2023 — *Biological Psychiatry* 95:639 (2024). doi:10.1016/j.biopsych.2023.05.022
- Parker 2020 — *JAMA Psychiatry* 77:1127 (2020). doi:10.1001/jamapsychiatry.2020.1495
- Patel 2021 — *JAMA Psychiatry* 78:47 (2021). doi:10.1001/jamapsychiatry.2020.2694
- PGC-CDG2 2019 — *Cell* 179:1469 (2019). doi:10.1016/j.cell.2019.11.020
- Pizzorusso 2002 — *Science* 298:1248 (2002). doi:10.1126/science.1072699
- Sekar 2016 — *Nature* 530:177 (2016). doi:10.1038/nature16549
- Sellgren 2019 — *Nature Neuroscience* 22:374 (2019). doi:10.1038/s41593-018-0334-7
- Shaw 2006 — *Nature* 440:676 (2006). doi:10.1038/nature04513
- Tooley 2021 — *Nature Reviews Neuroscience* 22:372 (2021). doi:10.1038/s41583-021-00457-5
- Vidal-Piñeiro 2020 — *Scientific Reports* 10 (2020). doi:10.1038/s41598-020-78471-3
- Whitaker 2016 — *PNAS* 113:9105 (2016). doi:10.1073/pnas.1601745113
