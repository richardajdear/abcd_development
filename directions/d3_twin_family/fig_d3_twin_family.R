#!/usr/bin/env Rscript
# fig_d3_twin_family.R -- one figure for DIRECTIONS.md D3 (twin and family designs).
# Reads only results/*.tsv written by 04-06; every printed statistic is read from them.
# Run from the repo root: Rscript directions/d3_twin_family/fig_d3_twin_family.R   (env r)
suppressMessages({library(ggplot2); library(patchwork); library(data.table)})
args <- commandArgs(FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
RES <- file.path(HERE, "results"); dir.create(file.path(HERE, "figures"), showWarnings = FALSE)
BASE <- 7
CR <- fread(file.path(HERE, "..", "..", "genetic_analysis", "current_results.tsv"))
snp <- unique(CR[readout == "GREML h2" & atlas == "hcp" & construction == "1lmm", .(trait, estimate)])
snp_slope <- snp[trait == "global_slope", estimate]; snp_ct <- snp[trait == "baseline_thickness", estimate]
PCN <- fread(file.path(RES, "pair_counts.tsv"))
dz_opp <- PCN[grepl("opposite", stage), DZ] / PCN[grepl("one pair", stage), DZ]
th <- theme_bw(base_size = BASE) + theme(
  panel.grid.minor = element_blank(), plot.title = element_text(size = BASE, face = "bold"),
  plot.subtitle = element_text(size = BASE - 1, colour = "grey25"), axis.text = element_text(size = BASE - 1),
  axis.title = element_text(size = BASE - 0.5), legend.text = element_text(size = BASE - 1),
  legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title.position = "plot")

LAB <- c(dCT_r = "Thinning rate (\u0394CT)", dCT_rq = "\u0394CT, quality-adjusted", CT0_r = "Cortical thickness",
         dT1T2_r = "T1w/T2w slope", pub_timing_r = "Puberty timing", pub_tempo_r = "Puberty tempo",
         depress_chg_r = "CBCL depressive \u0394", internal_chg_r = "CBCL internalising \u0394",
         external_chg_r = "CBCL externalising \u0394", pfactor_chg_r = "CBCL p-factor \u0394",
         fluid_gain_r = "Fluid cognition gain", cryst_gain_r = "Crystallised cognition gain",
         total_gain_r = "Total cognition gain")
ORD <- names(LAB)

# ---- a: twin correlations ------------------------------------------------------------
tc <- fread(file.path(RES, "twin_correlations.tsv"))
tcl <- melt(tc[, .(trait, rMZ, rDZ, rSIB)], id.vars = "trait", variable.name = "pair", value.name = "r")
tcl[, pair := factor(sub("^r", "", pair), levels = c("MZ", "DZ", "SIB"))]
tcl[, trait := factor(trait, levels = rev(ORD), labels = rev(LAB[ORD]))]
nmz <- tc[trait == "dCT_r", nMZ]; ndz <- tc[trait == "dCT_r", nDZ]; nsib <- tc[trait == "dCT_r", nSIB]
PCOL <- c(MZ = "#B2182B", DZ = "#2171b5", SIB = "grey55")
pa <- ggplot(tcl, aes(r, trait, colour = pair)) + geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(size = 1.4) + scale_colour_manual(values = PCOL) + th +
  labs(x = "co-twin / co-sibling correlation", y = NULL,
       title = sprintf("a   Twin similarity (%d MZ, %d DZ, %d sibling pairs)", nmz, ndz, nsib),
       subtitle = sprintf("\u0394CT: rMZ %.2f vs rDZ %.2f; siblings (scanned apart) %.2f",
                          tc[trait == "dCT_r", rMZ], tc[trait == "dCT_r", rDZ], tc[trait == "dCT_r", rSIB])) +
  theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))

# ---- b: ACE decomposition -------------------------------------------------------------
u <- fread(file.path(RES, "ace_univariate.tsv"))[model == "ACE"]
ul <- melt(u[, .(trait, A = h2, C = c2, E = e2)], id.vars = "trait", variable.name = "comp", value.name = "v")
ul[, comp := factor(comp, levels = c("E", "C", "A"))]
ul[, trait := factor(trait, levels = rev(ORD), labels = rev(LAB[ORD]))]
ue <- u[, .(trait = factor(trait, levels = rev(ORD), labels = rev(LAB[ORD])), h2, h2_lo, h2_hi)]
d0 <- u[trait == "dCT_r"]; c0 <- u[trait == "CT0_r"]
pb <- ggplot(ul, aes(v, trait, fill = comp)) + geom_col(width = 0.72) +
  geom_errorbar(data = ue, aes(xmin = h2_lo, xmax = h2_hi, y = trait), inherit.aes = FALSE,
                width = 0.3, linewidth = 0.3, orientation = "y") +
  scale_fill_manual(values = c(A = "#B2182B", C = "#fdae61", E = "grey85"), breaks = c("A", "C", "E")) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) + th +
  labs(x = "share of variance (ACE; bar = 95% CI of A)", y = NULL,
       title = sprintf("b   \u0394CT is heritable: A = %.2f [%.2f, %.2f], C = %.2f",
                       d0$h2, d0$h2_lo, d0$h2_hi, d0$c2),
       subtitle = sprintf("vs SNP h\u00b2 %.2f (GREML); thickness: twin A = %.2f vs SNP h\u00b2 %.2f",
                          snp_slope, c0$h2, snp_ct)) +
  theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0), axis.text.y = element_blank())

# ---- c: bivariate genetic / environmental correlations (AE) ---------------------------
bv <- fread(file.path(RES, "ace_bivariate.tsv"))[model == "AE"]
bl <- rbind(bv[, .(trait1, trait2, comp = "rA (genetic)", est = rA, lo = rA_lo, hi = rA_hi)],
            bv[, .(trait1, trait2, comp = "rE (unique env.)", est = rE, lo = rE_lo, hi = rE_hi)])
bl[, trait2 := factor(trait2, levels = rev(ORD), labels = rev(LAB[ORD]))]
bl[, which := ifelse(trait1 == "dCT_r", "\u0394CT", "\u0394CT, quality-adjusted")]
pt <- bv[trait1 == "dCT_r" & trait2 == "pub_timing_r"]
ci <- bv[trait1 == "dCT_r" & trait2 == "internal_chg_r"]
nb_bi <- unique(bv$n_boot)
nsig <- bv[trait1 == "dCT_r" & grepl("_chg_r$", trait2), sum(rA_lo > 0 | rA_hi < 0)]
pc <- ggplot(bl, aes(est, trait2, colour = comp, shape = which)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), size = 0.18, linewidth = 0.3) +
  scale_colour_manual(values = c("rA (genetic)" = "#B2182B", "rE (unique env.)" = "grey45")) +
  scale_shape_manual(values = c(16, 1)) + th +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(x = "correlation with \u0394CT (bivariate AE; higher \u0394CT = slower thinning)", y = NULL,
       title = sprintf("c   Genes for earlier puberty are genes for faster thinning (rA = %.2f)", pt$rA),
       subtitle = sprintf("95%% CI [%.2f, %.2f]; symptom change: %d of 4 rA exclude 0; internalising rE %.2f",
                          pt$rA_lo, pt$rA_hi, nsig, ci$rE)) +
  theme(legend.position = "top", legend.box = "horizontal", legend.margin = margin(0, 8, -4, 8))

# ---- d: within-family PRS -------------------------------------------------------------
wf <- fread(file.path(RES, "within_family_prs.tsv"))[outcome == "dCT"]
ws <- fread(file.path(RES, "within_family_sample.tsv"))
nf <- ws[grepl("families", set), n_families]
wl <- rbind(wf[, .(score, est = pop_full, se = pop_full_se, kind = "population (all)")],
            wf[, .(score, est = pop_fam, se = pop_fam_se, kind = "population (families)")],
            wf[, .(score, est = within, se = within_se, kind = "within family")],
            wf[, .(score, est = between, se = between_se, kind = "between family")])
wl[, kind := factor(kind, levels = rev(c("population (all)", "population (families)", "within family", "between family")))]
SL <- c(SCZ25_prscs = "SCZ 2025, PRS-CS", SCZ25_sbrc = "SCZ 2025, SBayesRC", MDD_prscs = "MDD, PRS-CS",
        MDD_sbrc = "MDD, SBayesRC")
wl[, score := factor(score, levels = rev(names(SL)), labels = rev(SL))]
mde <- range(wf$mde80_over_pop)
pd_ <- ggplot(wl, aes(est, score, colour = kind)) + geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_pointrange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), position = position_dodge(width = 0.7),
                  size = 0.18, linewidth = 0.3) +
  scale_colour_manual(values = c("population (all)" = "black", "population (families)" = "grey55",
                                 "within family" = "#B2182B", "between family" = "#2171b5"),
                      breaks = c("population (all)", "population (families)", "within family", "between family")) +
  th + guides(colour = guide_legend(nrow = 2)) +
  labs(x = "\u03b2 of score on \u0394CT (SD/SD; negative = faster thinning)", y = NULL,
       title = sprintf("d   Within-family PRS effects: not assessable at %d families", nf),
       subtitle = sprintf("detectable within-family effect (80%% power) = %.1f\u2013%.1f\u00d7 the population effect",
                          mde[1], mde[2])) +
  theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))

# ---- e: co-twin control ---------------------------------------------------------------
cc <- fread(file.path(RES, "cotwin_control.tsv"))[adjust == "none"]
cc[, design := factor(design, levels = rev(c("population", "within MZ", "within DZ", "within SIB")))]
OL <- c(depress = "Depressive", internal = "Internalising", external = "Externalising", pfactor = "p-factor")
cc[, outcome := factor(outcome, levels = rev(names(OL)), labels = rev(OL))]
npz <- fread(file.path(RES, "cotwin_control.tsv"))[design == "within MZ" & adjust == "none" & outcome == "depress", n_pairs]
cq <- fread(file.path(RES, "cotwin_control.tsv"))
ratio_mz <- cq[adjust == "none" & design == "within MZ", beta] / cq[adjust == "none" & design == "population", beta]
n_mz_sig <- cq[adjust == "none" & design == "within MZ", sum(hi < 0 | lo > 0)]
qchg <- cq[adjust == "quality" & design == "population" & outcome == "depress", beta]
pe <- ggplot(cc, aes(beta, outcome, colour = design)) + geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.72), size = 0.18, linewidth = 0.3) +
  scale_colour_manual(values = c("population" = "black", "within MZ" = "#B2182B", "within DZ" = "#2171b5",
                                 "within SIB" = "grey55"),
                      breaks = c("population", "within MZ", "within DZ", "within SIB")) + th +
  labs(x = "\u03b2 of \u0394CT on symptoms at 15\u201317 | baseline (per SD; negative = faster thinning, more symptoms)",
       y = NULL,
       title = sprintf("e   Co-twin control: the association does not shrink within %d MZ pairs", npz),
       subtitle = sprintf("within-MZ \u03b2 = %.1f\u2013%.1f\u00d7 population; CI excludes 0 for %d of 4; quality-adjusted depressive %.3f",
                          min(ratio_mz), max(ratio_mz), n_mz_sig, qchg)) +
  theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))

nb_uni <- 500  # 04_twin_ace.py default; bivariate count read from the table (n_boot)
ratio_h2 <- d0$h2 / snp_slope
cap <- paste(
  "\u2022 ABCD 7.0; 8,716 children with 2\u20134 scans. \u0394CT = single-LMM random slope of HCP-MMP mean thickness (the Figure-1 trait).",
  sprintf("\u2022 Pairs: genetically inferred zygosity (gn_y_genrel), both imaged; one pair per birth event / family. Only %.0f%% of DZ pairs are opposite-sex.", 100 * dz_opp),
  "\u2022 Traits residualised on sex, site, age and visit schedule (brain: + age span, n visits); CBCL / cognition as change given baseline.",
  sprintf("\u2022 a\u2013c: full-information ML twin models (ace.py); 95%% CIs from %d (univariate) / %d (bivariate) pair bootstraps.", nb_uni, nb_bi),
  "\u2022 Twins are scanned the same day: scan-level residuals correlate equally in MZ and DZ pairs (shared session error, loads on C).",
  "\u2022 d: pooled-arm (_zanc) scores; Mundlak between/within split, families with \u2265 2 genotyped children, MZ co-twins dropped; cluster SEs.",
  "\u2022 e: CBCL_late ~ \u0394CT + CBCL_base + age (+ sex + site); within-pair = pair fixed effects. Quality = mean and slope of log topological defects.",
  sep = "\n")
fig <- ((pa | pb) + plot_layout(widths = c(1.15, 0.85))) / pc / (pd_ | pe) +
  plot_layout(heights = c(1.15, 1, 0.8)) +
  plot_annotation(
    title = sprintf("Twin and family designs: thinning rate is %.1f\u00d7 as heritable in twins as in SNPs, shares genes with puberty timing, and shows no sign of familial confounding of its symptom link",
                    ratio_h2),
    subtitle = paste(sprintf("\u0394CT has twin h\u00b2 %.2f with no shared environment (b); its genetics overlap with earlier puberty (c) but not detectably with symptom change,", d0$h2),
                     "whose correlation with \u0394CT sits in the unique environment (c) and does not shrink within MZ pairs (e). Within-family PRS tests are underpowered (d).",
                     sep = "\n"),
    caption = cap,
    theme = theme(plot.title = element_text(size = BASE + 1, face = "bold"),
                  plot.subtitle = element_text(size = BASE, colour = "grey20"),
                  plot.caption = element_text(size = BASE - 1.5, hjust = 0, lineheight = 1.4, colour = "grey30")))
ggsave(file.path(HERE, "figures", "fig_d3_twin_family.png"), fig, width = 10, height = 10.4, dpi = 220)
cat("wrote figures/fig_d3_twin_family.png\n")
