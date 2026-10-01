# Directions memo — what mechanism are we proposing, and what would test it?

*2026-09-28; revised 2026-09-30 (§§2.3–2.4, T7–T8, D7–D8, §6) and 2026-10-01 (§3b scorecard after D2–D5, §§4–5). A planning document, not a results document. Every number quoted
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

*Status 2026-09-30 — done, [`directions/d2_prs_age_puberty/`](../directions/d2_prs_age_puberty/README.md).*
**T1 fails.** One-stage LMM on 26,597 scans: SCZ PRS × age −0.28 µm/yr per SD
(p 7e-4), MDD −0.17 (p 0.04); no level effect at 12.8 (bounds any phase advance
to < 1.1 months per SD); the rate effect is the same before and after the 12–14
peak; PRS × puberty is null once age is in the model. MDD PRS predicts earlier
puberty (+0.055 SD), and earlier puberty faster thinning (−0.068), but adjusting
for puberty leaves PRS → slope unchanged. Reading: uniform rate scaling across
9–17, not a timing shift. EA/ALZ/ASD scores (T3) not on the laptop.

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

*Status 2026-09-30 — done, [`directions/d3_twin_family/`](../directions/d3_twin_family/README.md)
(271 MZ / 432 DZ / 579 sibling pairs imaged; FIML twin models in `ace.py`, OpenMx
check pending on CSD3).* ΔCT twin A = 0.46 [0.37, 0.54], C = 0 (SNP h² 0.18).
Bivariate AE: rA with puberty timing −0.24 [−0.33, −0.15] (first support for
T2); with the T1w/T2w slope 0.08 [−0.12, 0.25]; with CBCL change 0.04–0.11, none
excluding 0, while the ΔCT–symptom covariance sits in E (internalising rE −0.16
[−0.27, −0.03]). **T6 holds:** within 250 MZ pairs the slope → symptom βs keep
the population sign and are 2.4–5.5× larger (2 of 4 CIs exclude 0). Within-family
PRS is not assessable (MDE 2.9–4.2× the population effect).

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

**Status (2026-09-30): run.** Code and tables in
[`directions/d3_twin_family/`](../directions/d3_twin_family/README.md). In brief: twin
h² of the thinning rate 0.46 [0.37, 0.54] with C = 0 (vs SNP h² 0.18). Its genetic
correlation with puberty timing is rA = −0.24 [−0.33, −0.15] (earlier puberty, faster
thinning; first support for T2), and it shares no detectable genetics with symptom change.
The slope–symptom association does not shrink within 250 MZ pairs (T6 not falsified,
wide CIs). Within-family PRS is not assessable at 1,087 families (MDE ≈ 3–4× the
population effect).

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

## 3b. Scorecard and synthesis after D2–D5 (2026-10-01)

D2, D3, D4 and D5 are run (laptop, ABCD 7.0, HCP-MMP single-LMM trait); D7 and
the MOSTest discovery arm are in progress on CSD3; D1, D6 and D8 are not
started. Every number here is a row of a `directions/*/results|tables/*.tsv`.

### 3b.1 Predictions

| # | prediction | result | verdict |
|:--|:--|:--|:--|
| T1 | PRS effect on the slope depends on age/puberty (phase advance) | rate effect identical across the 12–14 rate peak; no level effect at 12.8 (phase advance < 1.1 months/SD); PRS × PDS null given age (D2) | **fails** — uniform rate scaling |
| T2 | the same liability shifts other maturational clocks | genetic correlation ΔCT–puberty timing rA −0.24 [−0.33, −0.15]; MDD PRS → earlier puberty +0.055 SD (p 8e-7), SCZ PRS weaker (+0.015, n.s.); but ΔCT–T1w/T2w-slope rA 0.08 [−0.12, 0.25] (D3, D2) | **partly holds** — shared with the pubertal clock, not with the myelin readout |
| T3 | EA opposite in sign on every readout | not run: EA/ALZ/ASD score profiles are on CSD3 only | **blocked** (D2, D3, D4, D5 all name it) |
| T4 | PRS-related thinning follows the normative map | not re-tested; prior evidence is a uniform shift | open |
| T5 | faster thinners gain less cognitively | sign as predicted but small (crystallised +0.024 SD per SD, vocabulary +0.031, both FDR < 0.05); fluid/executive null; robust to SES, quality, interval (D4) | **weakly holds, wrong domain** — looks like a shared EA-type influence, not plasticity closure |
| T6 | slope → symptom survives within MZ pairs | within-MZ β same sign, 2.4–5.5× population, 2/4 CIs exclude 0; ACE: rA(ΔCT, symptom change) ≈ 0, covariance in E (D3); between children thinning precedes symptoms and baseline symptoms do not predict thinning (D5) | **holds** |
| T7 | SCZ/MDD effects are one shared component; EA acts via Cog | not run (D8) | open |
| T8 | slope predicts psychotic-like experiences | not run (D8; `mh_y_pps` fetched) | open |

### 3b.2 What the four studies say together

The results separate the thinning rate into **two components with different
causes and different consequences**.

**A heritable rate component, scaled by polygenic risk, that does not reach
symptoms.** Twin h² of the rate is 0.46 with no shared-environment term (D3);
SCZ and MDD scores scale it uniformly across 9–17 with no timing shift (D2);
the effect is independent of SES and adversity (4 % attenuation for SCZ, no
PRS × environment term; D5). Its genetics overlap with pubertal timing (rA
−0.24) and not detectably with the T1w/T2w slope (rA 0.08) or with symptom
change (rA 0.04–0.11, all CIs spanning 0). Puberty is a parallel correlate,
not the mediator: adjusting for pubertal stage leaves PRS → slope unchanged
(D2).

**A non-shared-environmental component that does reach symptoms.** The
slope–symptom covariance sits in E in the bivariate model and the association
does not shrink within MZ pairs (D3); between children it runs forward
(thinning → later symptoms; D5); parent-reported negative life events predict
faster thinning (−0.034, robust to site and scan quality; D5). Within children,
a small symptoms → later thinning lag also appears (RI-CLPM, −0.03 to −0.04),
to be re-tested with LCM-SR before it is called a finding.

Consequences for the memo's claims:

1. **The tempo hypothesis survives only in its amplitude form.** Risk carriers
   lose more cortex per year across the whole window (≈ 2 µm per SD of score
   over eight years); they are not further along the same trajectory. The
   "critical-period closure brought forward" reading of §2 is not what the
   data show within 9–17. Whether the extra loss is thinning that would
   otherwise occur after 17 is a question for the 8-year follow-up.
2. **The genetic–symptom chain is broken for a structural reason, not for lack
   of power.** The part of the rate that polygenic risk scales and the part
   that tracks symptoms are different parts of its variance. So the PRS →
   thinning → symptom mediation that fig 1 implies should not be claimed;
   what can be claimed is that genetic liability and adverse experience both
   accelerate thinning, additively, and that the symptom-linked acceleration is
   the experience-dependent one. This matches adversity-accelerated maturation
   (Tooley 2021) and the stress-driven synapse-loss model of depression
   (§2.3) better than the pruning-genetics model of schizophrenia.
3. **Pruning versus myelination: a first, weak lean toward neuropil.** The
   rate's genetics are not the T1w/T2w slope's genetics (rA 0.08). If the
   T1w/T2w slope is an adequate myelination readout (DK only, no bias-field
   correction, reliability 0.27–0.37), the heritable thinning component is not
   myelination. D1's intensity and diffusion slopes decide this; the D3
   pipeline gives each new modality its twin h² and rA with ΔCT for free.
4. **SES acts on level, polygenic risk on rate.** SES predicts thickness at
   12.8 (+0.062) and not the rate within site (D5); C for the rate is 0 while C
   for thickness is 0.08 (D3). The level and the rate have different inputs,
   which is why fig 1's "rate not level" result for the scores is informative
   rather than incidental.
5. **MDD and SCZ already look different (relevant to D8).** The MDD score
   carries gene–environment correlation with adversity (+0.093) and lower SES,
   predicts earlier puberty strongly, and attenuates 18 % with the environment
   in the model; the SCZ score does none of these. MDD's thinning effect may
   run partly through correlated environment and the pubertal clock; SCZ's
   looks like a direct polygenic effect on rate.
6. **The cognitive consequence is the EA axis, not plasticity.** T5's effect
   is on vocabulary, the most education-loaded measure and the one EA predicts
   best, and slower thinners already score higher at 10. With EA predicting
   slower thinning, the parsimonious reading is a shared EA-type influence on
   both. The EA score (T3) and an MZ-difference test of vocabulary gain decide
   it; the bivariate twin model cannot (cognitive gain is too noisy to
   decompose: rMZ < rDZ).

Caveats carried forward: the pubertal-timing rA could be shared body-size/BMI
genetics rather than a maturational clock (test with BMI as a third trait);
the twin sample is small and 86 % same-sex DZ; D4 is concurrent development,
not prediction; the life-events result is one informant and one wave; D5's
170 coefficients are uncorrected.

### 3b.3 What this changes in the plan

- **One rsync unblocks four directions.** The EA, ALZ ± APOE and ASD
  `score_*.profile` files from CSD3 (`results_70tab_hcp/prs_final_1lmm/` and
  `legacy/hpc_v2/.../PRSCS/EA/`; commands in the D2 and D5 READMEs) give T3 in
  D2, the EA triangle in D4, the EA-with-SES step in D5 and within-family EA in
  D3. This is the highest-value step and costs no compute.
- **D1 rises in priority.** With T1 failed, the amplitude claim needs a tissue
  process, and D3 has already produced a prior (rA with T1w/T2w ≈ 0). Run the
  diffusion/intensity slopes (fetch list in `abcd70_fetch_list.csv`), and put
  every new slope through the D3 twin pipeline for h² and rA with ΔCT.
- **New test T9 (links D3 to D5):** within MZ pairs, does the twin with more
  negative life events thin faster, and does that difference carry the
  symptom difference? This is the E-component test the co-twin result asks
  for, it is cheap (same pairs, `mh_p_ple` by child), and it is the only design
  here that can show adversity → thinning without genetic confounding.
- **D8 gains a concrete hypothesis** from point 5: MDD-specific and
  SCZ-specific (subtracted) scores should differ in their relation to puberty,
  adversity and SES, not only in β on the slope.
- **Add BMI** (`ph_y_anthr`, vendored) as a third trait to the D3 bivariate
  models to test whether the ΔCT–puberty rA is body-size genetics.
- **Write-up:** fig 2 becomes a two-component figure (heritable, PRS-scaled,
  puberty-linked rate vs environment-linked, symptom-predictive rate), and the
  fig 1 → fig 2 text should stop implying a PRS → thinning → symptom chain.

## 4. Suggested order (revised 2026-10-01; original order in git history)

Done: D2, D3, D4, D5 (laptop). Running on CSD3: D7, MOSTest discovery arm.

1. **Pull the EA / ALZ ± APOE / ASD score profiles from CSD3** (rsync in the D2
   and D5 READMEs) and re-run D2 `01–03`, D4 with EA as covariate and PRS
   triangle, D5's EA step, D3 `05` with EA. Closes T3 and the D4 reading in a
   day, no compute.
2. **D1 on the fetched tables** — diffusion (`is` DTI, RSI in cortical GM and
   superficial WM), volume, area, GWC, T1w/T2w slopes with the settled
   specification; each through the D3 twin pipeline (h², rA with ΔCT) and the
   PRS panel; then the cross-modal tempo composite.
3. **T9 and BMI in D3**: MZ-discordant life events → ΔCT → symptoms; BMI as a
   third trait against the ΔCT–puberty rA.
4. **D8**: joint model now (scores exist); Cog/NonCog, subtracted and
   cross-disorder scores via the existing PRS pipeline on CSD3; PQ-BC outcome
   (`mh_y_pps`) locally. Test the point-5 hypothesis (MDD vs SCZ routes).
5. **LCM-SR** replacement for the D5 RI-CLPM before any within-child claim.
6. D6 in the snRNA-seq repo; its gene sets become D7's K4 slot.
7. D7 and MOSTest results when CSD3 returns them; read D7 against the
   pre-registration only.

## 5. Revised arc (2026-10-01)

- **Fig 1** (unchanged): individual thinning rate tracks polygenic risk in
  both directions (SCZ/MDD faster, EA slower), not baseline thickness, and
  tracks later symptoms.
- **Fig 2 — two components of the rate.** (a) heritable (twin 0.46, C 0),
  scaled uniformly by polygenic risk across 9–17 with no timing shift, shares
  genes with pubertal timing, independent of SES/adversity; (b) non-shared
  environmental, tracks life events and later symptoms, survives the MZ
  co-twin control. The PRS → thinning → symptom chain is not claimed.
- **Fig 3 — tissue and programme.** D1 multimodal slopes and their twin
  genetics (myelin vs neuropil), the snRNA-seq maturation programme and its
  adult footprint (C3) against the normative map, disorder-GWAS enrichment of
  the D6 gene sets, the pre-registered D7 check — group-level where it is,
  individual-level nulls reported.
- **Fig 4 — consequences and specificity.** Cognitive gain (EA axis vs
  plasticity), PLEs and depressive symptoms, shared vs disorder-specific
  scores, MDD's puberty/adversity route vs SCZ's direct one.

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
