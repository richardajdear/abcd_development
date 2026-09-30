#!/usr/bin/env Rscript
# D2 figure -- reads only tables/*.tsv written by 02_fit_models.R; recomputes nothing.
# Run (env r, from the repo root):  Rscript directions/d2_prs_age_puberty/03_figure.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
D <- "directions/d2_prs_age_puberty"
terms <- fread(file.path(D, "tables/onestage_terms.tsv"))
seg   <- fread(file.path(D, "tables/prs_rate_by_age.tsv"))
pt    <- fread(file.path(D, "tables/puberty_tempo.tsv"))
vc    <- fread(file.path(D, "tables/variance_components.tsv"))
COL <- c(SCZ = "#762a83", MDD = "#1b7837")
BASE <- 7.6
th <- theme_bw(base_size = BASE) + theme(panel.grid.minor = element_blank(), legend.position = "none",
                                        plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = BASE - 0.8))
um <- function(x) x * 1000    # mm/yr -> um/yr

# ---- a. normative thinning rate by age segment (SCZ pooled sbrc seg3 fixed effects) ------------
s3 <- terms[disorder == "SCZ" & method == "sbrc" & arm == "pooled" & model == "seg3"]
b <- setNames(s3$est, s3$term)
norm <- data.table(segment = factor(c("<12", "12-14", ">=14"), levels = c("<12", "12-14", ">=14")),
                   rate = um(c(b["age_c"], b["age_c"] + b["h12"], b["age_c"] + b["h12"] + b["h14"])))
lin <- terms[disorder == "SCZ" & method == "sbrc" & arm == "pooled" & model == "linear" & term == "age_c", est]
pa <- ggplot(norm, aes(segment, rate)) + geom_col(fill = "grey55", width = 0.6) +
  geom_hline(yintercept = um(lin), linetype = 2, linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", rate)), vjust = 1.4, colour = "white", size = 2.4) +
  labs(title = "a   Normative thinning rate peaks at 12-14",
       subtitle = sprintf("mean HCP-MMP thickness, um/yr; dashed = linear fit %.1f", um(lin)),
       x = "age segment (years)", y = "thinning rate (um/yr)") + th

# ---- b. PRS effect on the rate, overall and by segment -----------------------------------------
ov <- terms[method == "sbrc" & arm == "pooled" & model == "linear" & term == "prs:age_c",
            .(disorder, segment = "9-17 (linear)", est, se, p)]
sg <- seg[method == "sbrc" & arm == "pooled" & model == "seg3", .(disorder, segment, est, se, p)]
fb <- rbind(ov, sg)
fb[, segment := factor(segment, levels = rev(c("9-17 (linear)", "<12", "12-14", ">=14")))]
fb[, disorder := factor(disorder, levels = c("SCZ", "MDD"))]
h13 <- terms[method == "sbrc" & arm == "pooled" & model == "hinge13" & term == "prs:h13"]
pb <- ggplot(fb, aes(y = segment, x = um(est), xmin = um(est - 1.96 * se), xmax = um(est + 1.96 * se), colour = disorder)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.25, linewidth = 0.4) +
  scale_colour_manual(values = COL) +
  labs(title = "b   PRS scales the rate; no age dependence",
       subtitle = sprintf("PRS x rate change after 13: SCZ %+.2f (p %.2f), MDD %+.2f (p %.2f)",
                          um(h13[disorder == "SCZ", est]), h13[disorder == "SCZ", p],
                          um(h13[disorder == "MDD", est]), h13[disorder == "MDD", p]),
       x = "PRS effect on thinning rate (um/yr per SD of score; 95% CI)", y = NULL) +
  th + theme(legend.position = c(0.85, 0.2), legend.title = element_blank(), legend.background = element_blank())

# ---- c. level effect at 12.8 (bounds any phase shift) -----------------------------------------
lv <- terms[method == "sbrc" & arm == "pooled" & model == "linear" & term == "prs", .(disorder, est, se, p)]
lv[, disorder := factor(disorder, levels = c("SCZ", "MDD"))]
# a phase advance of D years on a linear trajectory thins the cortex by |rate| * D at every age
dmax <- (abs(lv$est) + 1.96 * lv$se) / abs(lin)
pc <- ggplot(lv, aes(y = disorder, x = um(est), xmin = um(est - 1.96 * se), xmax = um(est + 1.96 * se), colour = disorder)) +
  geom_vline(xintercept = 0, linewidth = 0.3) + geom_pointrange(size = 0.25, linewidth = 0.4) +
  scale_colour_manual(values = COL) +
  labs(title = "c   No thickness difference at 12.8",
       subtitle = sprintf("bounds any phase advance to < %.1f months/SD (SCZ), < %.1f (MDD)", 12 * dmax[1], 12 * dmax[2]),
       x = "PRS effect on thickness at age 12.8 (um per SD; 95% CI)", y = NULL) + th

# ---- d. puberty: PRS -> pubertal timing; timing -> thinning; PRS -> thinning given puberty ------
pp <- pt[method == "sbrc" & arm == "pooled"]
d1 <- pp[outcome == "pds_stage_z" & term == "prs", .(disorder, what = "PRS -> pubertal stage at 12.8", est, se, p)]
d2 <- pp[outcome == "slope_z" & term == "prs", .(disorder, what = "PRS -> thinning slope", est, se, p)]
d3 <- pp[outcome == "slope_z|puberty" & term == "prs", .(disorder, what = "PRS -> thinning slope | puberty", est, se, p)]
d4 <- pt[disorder == "none" & outcome == "slope_z" & term == "pds_stage_z", .(disorder = "none", what = "pubertal stage -> thinning slope", est, se, p)]
fd <- rbind(d1, d2, d3, d4)
fd[, what := factor(what, levels = rev(c("PRS -> pubertal stage at 12.8", "pubertal stage -> thinning slope",
                                         "PRS -> thinning slope", "PRS -> thinning slope | puberty")))]
fd[, disorder := factor(disorder, levels = c("SCZ", "MDD", "none"))]
pd_ <- ggplot(fd, aes(y = what, x = est, xmin = est - 1.96 * se, xmax = est + 1.96 * se, colour = disorder)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.25, linewidth = 0.4) +
  scale_colour_manual(values = c(COL, none = "grey30")) +
  labs(title = "d   Puberty tracks both, not the PRS effect",
       subtitle = sprintf("PRS -> slope given puberty: SCZ %.3f -> %.3f; MDD %.3f -> %.3f",
                          d2[disorder == "SCZ", est], d3[disorder == "SCZ", est], d2[disorder == "MDD", est], d3[disorder == "MDD", est]),
       x = "standardised effect (SD per SD; 95% CI)", y = NULL) + th

n_obs <- terms[method == "sbrc" & arm == "pooled" & model == "linear" & disorder == "SCZ", n_obs][1]
n_ch  <- terms[method == "sbrc" & arm == "pooled" & model == "linear" & disorder == "SCZ", n_children][1]
tau1  <- um(vc[disorder == "SCZ" & method == "sbrc" & arm == "pooled" & model == "linear" & grp == "subject" & var1 == "age_c" & (is.na(var2) | var2 == ""), sdcor])
cap <- paste(
  sprintf("\u2022 %d scans of %d genotyped children (ABCD 7.0, HCP-MMP, pooled arm); scores: SCZ 2025 and MDD, SBayesRC, z within ancestry cluster", n_obs, n_ch),
  "\u2022 one-stage LMM: mean thickness ~ PRS x age basis + sex + PC1-10 + (1 + age | child) + (1 | site) + (1 | family); age centred at 12.797",
  sprintf("\u2022 a: piecewise-linear age basis with knots at 12 and 14; b: PRS x each segment; the SD of the child slope is %.1f um/yr, so -0.28 um/yr = -0.07 SD", tau1),
  "\u2022 c: PRS main effect = thickness difference at the age centre; a phase advance of D years on a trajectory thinning at 20 um/yr would show as 20 D um",
  "\u2022 d: pubertal stage = child intercept of PDS (youth report) ~ age x sex + (1 + age | child) + (1 | site) over all yearly waves, i.e. PDS at 12.8; child-level models carry sex, PC1-10, (1 | family)",
  "\u2022 EUR-arm and PRS-CS versions in tables/; both agree in sign and magnitude", sep = "\n")
fig <- (pa | pb) / (pc | pd_) + plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Polygenic risk for SCZ and MDD scales the rate of adolescent thinning uniformly from 9 to 17; it does not shift its timing",
                  subtitle = paste("The normative rate peaks at 12-14 (a). A phase advance would flip the sign of the PRS effect across that peak and thin the cortex at 12.8;",
                                   "neither happens (b, c). Pubertal timing relates to both PRS and thinning but explains none of the PRS effect (d).", sep = "\n"),
                  caption = cap,
                  theme = theme(plot.title = element_text(face = "bold", size = BASE + 2), plot.subtitle = element_text(size = BASE),
                                plot.caption = element_text(hjust = 0, lineheight = 1.45, size = BASE - 1.2)))
ggsave(file.path(D, "figures/fig_d2_prs_age_puberty.png"), fig, width = 8.2, height = 6.4, dpi = 200, bg = "white")
cat("written\n")
