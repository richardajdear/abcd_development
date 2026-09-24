#!/usr/bin/env Rscript
# fig_cbcl_candidate.R -- DRAFT of a symptoms panel for the empty slot of the
# HCP-MMP summary slide. Standalone, not yet placed in fig5_hcp_summary.R.
# Reads results/cbcl_explore_{assoc,spin}.tsv (7.0; coefficients only).
# Writes figures/fig_cbcl_candidate.png.
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
A <- read.delim(file.path(ROOT, "results", "cbcl_explore_assoc.tsv"))
S <- read.delim(file.path(ROOT, "results", "cbcl_explore_spin.tsv"))
A$adj_global <- as.logical(A$adj_global)          # pandas writes True/False
stopifnot(all(A$source == "7.0 p/mh_p_cbcl.tsv"))
FADE <- 0.3; TXT <- 7
BR <- c(global_thin = "global thinning", pls1_thick_proj = "PLS1 (baseline CT)", pls2_proj = "PLS2 (thinning | global)")
COL <- c("global thinning" = "grey35", "PLS1 (baseline CT)" = "#2166ac", "PLS2 (thinning | global)" = "#b2182b")
pick <- function(d) d |> filter((brain == "global_thin" & !adj_global) | (brain == "pls1_thick_proj" & !adj_global) |
                                (brain == "pls2_proj" & adj_global)) |>
  mutate(br = factor(BR[brain], rev(BR)))
CB <- c(pfactor = "p-factor", internal = "Internalising", external = "Externalising",
        depress = "Depressive (DSM)", thought = "Thought problems")
DX <- c(mdd_youth_DX = "MDD (youth report)", mdd_parent_DX = "MDD (parent report)",
        psychosis_parent_DX = "Psychosis spectrum (parent)")
a <- pick(A |> filter(model == "change", outcome %in% names(CB))) |>
  mutate(o = factor(CB[outcome], rev(CB)), lo = y_sd_beta - 1.96 * se * y_sd_beta / beta,
         hi = y_sd_beta + 1.96 * se * y_sd_beta / beta, est = y_sd_beta)
k <- pick(A |> filter(model == "ksads_logit", outcome %in% names(DX))) |>
  mutate(o = factor(sprintf("%s  (%d cases)", DX[outcome], n_cases), rev(sprintf("%s  (%d cases)", DX[names(DX)],
                    sapply(names(DX), \(x) n_cases[outcome == x][1])))),
         est = exp(beta), lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se))
stopifnot(nrow(a) == 15, nrow(k) == 9)
forest <- function(d, ref, xl, title) {
  ggplot(d, aes(est, o, colour = br, alpha = p < 0.05)) +
    geom_vline(xintercept = ref, colour = "grey55", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.45,
                  position = position_dodge(width = 0.6)) +
    geom_point(size = 1.4, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = COL, name = NULL, breaks = unname(BR)) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
    labs(x = xl, y = NULL, title = title) +
    theme_classic(base_size = TXT) +
    theme(plot.title = element_text(size = TXT, face = "bold", hjust = 0), legend.position = "bottom",
          axis.line.y = element_blank(), axis.ticks.y = element_blank())
}
spin <- S |> filter(brain == "pls2_proj", model == "change")
pa <- forest(a, 0, "\u03b2 per SD (symptom change, ages ~10 \u2192 15\u201317)", "CBCL symptom change")
pk <- forest(k, 1, "odds ratio per SD", "KSADS lifetime diagnosis")
fig <- (pa | pk) + plot_layout(guides = "collect") +
  plot_annotation(title = "f  Symptoms (exploratory)",
                  caption = sprintf("PLS2 projection vs spin-rotated maps (CBCL change): min p_spin = %.2f over %d outcomes. n \u2248 8,200.",
                                    min(spin$p_spin), nrow(spin)),
                  theme = theme(plot.title = element_text(size = TXT + 0.8, face = "bold"),
                                plot.caption = element_text(size = TXT - 1, colour = "grey35"))) &
  theme(legend.position = "bottom")
ggsave(file.path(ROOT, "figures", "fig_cbcl_candidate.png"), fig, width = 7.6, height = 2.9, dpi = 300, bg = "white")
cat("wrote fig_cbcl_candidate.png\n")
