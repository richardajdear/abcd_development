#!/usr/bin/env Rscript
# D4 figure: thinning rate and baseline thickness vs cognitive gain (baseline -> year 6).
# Reads only results/d4_assoc.tsv and results/d4_sample.tsv (03_assoc.py); recomputes nothing.
# Run from repo root:  Rscript directions/d4_cognitive_gain/code/fig_d4_cognitive_gain.R  (env ahba-pls-r)
suppressMessages({library(ggplot2); library(dplyr); library(patchwork); library(readr)})
D4 <- "directions/d4_cognitive_gain"
A <- read_tsv(file.path(D4, "results/d4_assoc.tsv"), show_col_types = FALSE)
S <- read_tsv(file.path(D4, "results/d4_sample.tsv"), show_col_types = FALSE)
sv <- function(k) S$value[S[[1]] == k]
RED <- "#B2182B"; BLUE <- "#2171b5"
lab <- c(fluid = "Fluid composite", cryst = "Crystallised composite", total = "Total composite",
         flanker = "Flanker", pattern = "Pattern comparison", picseq = "Picture sequence",
         picvocab = "Picture vocabulary", reading = "Oral reading", cardsort = "Card sort",
         listsort = "List sorting")
ord <- rev(names(lab))
th <- theme_bw(base_size = 7.6) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 8.2), plot.subtitle = element_text(size = 7),
        legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(3, "mm"))
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

# ---- a: slope vs baseline thickness, M2 (SES-adjusted) --------------------------------
a <- A |> filter(model == "M2", term %in% c("global_slope", "baseline_thickness")) |>
  mutate(outcome = factor(outcome, ord),
         term = factor(term, c("baseline_thickness", "global_slope"),
                       c("CT (baseline thickness)", "\u0394CT (thinning rate)")))
sl <- a |> filter(term == "\u0394CT (thinning rate)")
n_sig_sl <- sum(sl$p_fdr < 0.05); n_sig_ct <- sum(a$p_fdr[a$term != "\u0394CT (thinning rate)"] < 0.05)
ct1 <- A |> filter(model == "M1", term == "baseline_thickness"); n_sig_ct1 <- sum(ct1$p_fdr < 0.05)
ct_cr1 <- ct1$est[ct1$outcome == "cryst"]; ct_cr2 <- a$est[a$outcome == "cryst" & a$term != "\u0394CT (thinning rate)"]
sl_cr1 <- A$est[A$model == "M1" & A$term == "global_slope" & A$outcome == "cryst"]
cr <- sl |> filter(outcome == "cryst"); pv <- sl |> filter(outcome == "picvocab")
fl <- sl |> filter(outcome == "fluid")
pa <- ggplot(a, aes(est, outcome, colour = term)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi, shape = p_fdr < 0.05), size = 0.45, linewidth = 0.4,
                  position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  scale_colour_manual(values = setNames(c(BLUE, RED), levels(a$term))) +
  scale_y_discrete(labels = lab) +
  labs(x = "\u03b2 (SD of year-6 score per SD of brain phenotype)", y = NULL,
       title = sprintf("a   Slower thinning goes with larger crystallised gain (\u03b2 = %.3f)", cr$est),
       subtitle = sprintf(paste0("SES-adjusted. \u0394CT: %d of 10 FDR < 0.05 (vocabulary \u03b2 = %.3f), fluid null (p = %.2f).\n",
                                 "CT: %d of 10 without SES, %d with it (crystallised %.3f \u2192 %.3f); \u0394CT %.3f \u2192 %.3f."),
                          n_sig_sl, pv$est, fl$p, n_sig_ct1, n_sig_ct, ct_cr1, ct_cr2, sl_cr1, cr$est)) + th

# ---- b: robustness of the thinning-rate effect, composites + vocabulary ---------------
mods <- c(M1 = "ANCOVA", M2 = "+ SES", M3 = "+ CT, image quality", DIFF = "difference score",
          M1_AGECORR = "age-corrected scores", M1_3SCANS = "\u2265 3 scans", M1_Y4 = "year 4 (not year 6)")
b <- A |> filter(term == "global_slope", model %in% names(mods),
                 outcome %in% c("fluid", "cryst", "picvocab")) |>
  mutate(model = factor(model, rev(names(mods)), rev(mods)),
         outcome = factor(outcome, c("fluid", "cryst", "picvocab"), lab[c("fluid", "cryst", "picvocab")]))
cr_all <- b |> filter(outcome == lab["cryst"])
pb <- ggplot(b, aes(est, model)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi, shape = p < 0.05), colour = RED, size = 0.45, linewidth = 0.4) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  facet_wrap(~outcome, nrow = 1) +
  labs(x = "\u03b2 for \u0394CT (SD per SD)", y = NULL,
       title = sprintf("b   The crystallised effect survives every specification (\u03b2 %.3f to %.3f)",
                       min(cr_all$est), max(cr_all$est)),
       subtitle = sprintf("Crystallised p %.3f\u2013%.3f across %d specifications; fluid never p < 0.05.",
                          min(cr_all$p), max(cr_all$p), nrow(cr_all))) + th

fig <- (pa | pb) + plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    title = "Slower individual thinning goes with larger crystallised, not fluid, cognitive gain from age 10 to 16 \u2014 a small, SES-independent effect",
    subtitle = paste0("The sign matches the plasticity prediction (faster thinning, smaller gain), but the effect is small and confined to vocabulary-type measures (a).\n",
                     "Unlike baseline thickness, whose association with gain is largely SES, it is unchanged by SES, image quality, scoring or interval (b)."),
    caption = paste(
      sprintf("\u2022 ABCD 7.0 imaging children (n = %s) with NIH Toolbox at baseline and year 6 (n = %s crystallised, %s fluid); mean age %.1f \u2192 %.1f y.",
              format(sv("n_imaging"), big.mark = ","), format(sv("n_cryst_both"), big.mark = ","),
              format(sv("n_fluid_both"), big.mark = ","), sv("age_base_mean"), sv("age_late_mean")),
      "\u2022 Brain: single LMM on per-scan HCP-MMP cortical mean (Figure 1 phenotypes), random slope (\u0394CT) and intercept (CT), z-scored.",
      "\u2022 Sign: \u0394CT negative = faster thinning, so \u03b2 > 0 means slower thinning goes with larger gain; the plasticity prediction (T5) is \u03b2 > 0.",
      "\u2022 Model: y_year6 ~ brain + y_baseline + age_baseline + age_year6 + sex + site; M2 adds caregiver education, household income, ADI.",
      "\u2022 NIH Toolbox uncorrected standard scores; \u03b2 in SD of the year-6 score; OLS with family-clustered SE; filled = FDR < 0.05 (a) or p < 0.05 (b), BH within model over 10 outcomes.",
      "\u2022 Fluid and total composites exist only at baseline and year 6; year-4 row uses the tasks given at year 4.",
      sep = "\n"),
    theme = theme(plot.title = element_text(face = "bold", size = 9), plot.subtitle = element_text(size = 7.4),
                  plot.caption = element_text(hjust = 0, size = 6.4, lineheight = 1.45)))
ggsave(file.path(D4, "figures/fig_d4_cognitive_gain.png"), fig, width = 10.5, height = 5.2, dpi = 220)
cat("wrote", file.path(D4, "figures/fig_d4_cognitive_gain.png"), "\n")
