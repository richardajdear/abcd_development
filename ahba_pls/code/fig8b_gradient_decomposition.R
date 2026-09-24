#!/usr/bin/env Rscript
# fig8b_gradient_decomposition.R -- where the flatter gradient comes from: later
# symptoms vs each child's thinning rate in the normatively fast / mid / slow
# thirds of cortex. Reads results/gradient_scores_decomposition.tsv (26). Group-level.
# Writes figures/fig_gradient_decomposition.png
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
DC <- read.delim(file.path(ROOT, "results", "gradient_scores_decomposition.tsv"))
TXT <- 7; FADE <- 0.3
OUT <- c(totprob = "Total problems", pfactor = "p-factor", internal = "Internalising", anxdisord = "Anxiety (DSM)",
         rulebreak = "Rule-breaking", external = "Externalising")
TIER <- c(thin_fast = "fast-thinning third", thin_mid = "middle third", thin_slow = "slow-thinning third")
TCOL <- c("fast-thinning third" = "#b2182b", "middle third" = "grey45", "slow-thinning third" = "#2166ac")
d <- DC |> filter(outcome %in% names(OUT), term %in% names(TIER)) |>
  mutate(o = factor(OUT[outcome], rev(OUT)), tier = factor(TIER[term], TIER),
         lo = y_sd_beta - 1.96 * se * y_sd_beta / beta, hi = y_sd_beta + 1.96 * se * y_sd_beta / beta)
panel <- function(specs, title, sub) {
  ggplot(filter(d, spec %in% specs), aes(y_sd_beta, o, colour = tier, alpha = p < 0.05)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.45, position = position_dodge(0.6)) +
    geom_point(size = 1.5, position = position_dodge(0.6)) +
    scale_colour_manual(values = TCOL, name = NULL, limits = unname(TIER), drop = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
    labs(x = "\u03b2 per SD of thinning rate (SD of symptoms at 15\u201317)", y = NULL, title = title, subtitle = sub) +
    theme_classic(base_size = TXT) +
    theme(plot.title = element_text(size = TXT + 0.5, face = "bold"), plot.subtitle = element_text(size = TXT - 0.5, colour = "grey30"),
          legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank())
}
sf <- DC |> filter(outcome == "totprob", spec == "fast + slow | CT, QC", term == "slow - fast")
pa <- panel(c("fast alone", "mid alone", "slow alone"), "a  Each third on its own",
            "Which regions do children with rising symptoms thin more in?")
pb <- panel("fast + slow | CT, QC", "b  Fast and slow thirds together",
            sprintf("| baseline CT, image quality. Total problems: slow \u2212 fast = %.3f, p = %.1g", sf$y_sd_beta, sf$p))
fig <- (pa | (pb + theme(legend.position = "none")))   # one legend: b uses the same colours
fig <- fig + plot_annotation(
  title = "Children with rising general symptoms thin more in normally slow-thinning cortex; fast-thinning cortex thins at a typical rate",
  subtitle = paste0("Parcels split into thirds by normative thinning rate (fast 26, middle 20, slow 13 \u00b5m/yr). Later CBCL ~ tier rate(s) + baseline CBCL + sex + site + ages + scans; no global-thinning covariate.\n",
                    "Jointly (b), faster thinning in slow cortex goes with more problems and, at the same slow-cortex rate, faster fast-cortex thinning with fewer \u2014 hence the flatter gradient.\nRule-breaking rises with thinning in all three thirds alike, so its link to an atypical pattern (gradient fit) is not a flattening."),
  theme = theme(plot.title = element_text(size = TXT + 2.5, face = "bold"), plot.subtitle = element_text(size = TXT, colour = "grey25", lineheight = 1.15)))
ggsave(file.path(ROOT, "figures", "fig_gradient_decomposition.png"), fig, width = 9.5, height = 4.3, dpi = 300, bg = "white")
cat("wrote fig_gradient_decomposition.png\n")
