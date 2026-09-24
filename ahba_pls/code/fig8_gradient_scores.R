#!/usr/bin/env Rscript
# fig8_gradient_scores.R -- EXPLORATORY: per-child "how closely does this child's
# thinning follow the normative gradient" scores vs later CBCL symptoms.
# Reads results/gradient_scores_{assoc,bins}.tsv (26_gradient_scores.py),
# results/hcp_summary_maps.csv, data/hcp_polygons.csv. Group-level only.
# Writes figures/fig_gradient_scores.png
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"; RES <- file.path(ROOT, "results")
A <- read.delim(file.path(RES, "gradient_scores_assoc.tsv"))
B <- read.delim(file.path(RES, "gradient_scores_bins.tsv"))
G <- read.csv(file.path(RES, "hcp_summary_maps.csv"))
poly <- read.csv(file.path(ROOT, "data", "hcp_polygons.csv")) |> filter(view %in% c("lateral", "medial"))
TXT <- 7; FADE <- 0.3
th <- theme_classic(base_size = TXT) +
  theme(plot.title = element_text(size = TXT + 0.5, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = TXT - 0.5, colour = "grey30"), strip.background = element_blank())
OUT <- c(rulebreak = "Rule-breaking", external = "Externalising", totprob = "Total problems", pfactor = "p-factor",
         depress = "Depressive (DSM)", withdep = "Withdrawn/depressed", internal = "Internalising",
         anxdep = "Anxious/depressed", anxdisord = "Anxiety (DSM)")
SC <- c(grad_fit = "gradient fit (pattern r)", grad_slope = "gradient slope (steepness)", pls2_beyond_grad = "PLS2 beyond gradient")
SCOL <- c("gradient fit (pattern r)" = "#b2182b", "gradient slope (steepness)" = "#2166ac", "PLS2 beyond gradient" = "grey45")

# a: the normative gradient itself
gm <- G |> mutate(label = tolower(label), g = -dCT * 1000)
pd <- poly |> mutate(label = tolower(label)) |> inner_join(gm, by = "label")
pa <- ggplot(pd, aes(x, y, group = interaction(view, label, group, subgroup), fill = g)) +
  geom_polygon(colour = "grey40", linewidth = 0.05) + coord_fixed(expand = FALSE) +
  scale_fill_gradientn(colours = c("white", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"), name = "\u00b5m/yr") +
  labs(title = "a  The normative thinning gradient",
       subtitle = "group mean thinning rate; each child's pattern is compared with this map") +
  theme_void(base_size = TXT) +
  theme(plot.title = element_text(size = TXT + 0.5, face = "bold"), plot.subtitle = element_text(size = TXT - 0.5, colour = "grey30"),
        legend.key.height = unit(10, "pt"), legend.key.width = unit(5, "pt"))

# b: forest, fully adjusted, all children
fb <- A |> filter(sample == "all", adjust == "+global+CT+QC", score %in% names(SC)) |>
  mutate(sc = factor(SC[score], SC), o = factor(OUT[outcome], rev(OUT)),
         lo = y_sd_beta - 1.96 * se * y_sd_beta / beta, hi = y_sd_beta + 1.96 * se * y_sd_beta / beta)
stopifnot(nrow(fb) == 27)
pb <- ggplot(fb, aes(y_sd_beta, o, colour = sc, alpha = p < 0.05)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.45, position = position_dodge(0.65)) +
  geom_point(size = 1.4, position = position_dodge(0.65)) +
  scale_colour_manual(values = SCOL, name = NULL) + scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  labs(x = "\u03b2 per SD of score (SD of symptoms at 15\u201317)", y = NULL, title = "b  Scores vs symptoms at 15\u201317",
       subtitle = "| baseline symptoms, global thinning, baseline CT & its gradient, image quality; n \u2248 8,200") +
  th + theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(colour = guide_legend(nrow = 1))

dec <- function(sc, outs, title, sub) {
  d <- B |> filter(score == sc, outcome %in% outs) |> mutate(o = factor(OUT[outcome], OUT[outs]))
  ggplot(d, aes(decile, mean_resid_sd, colour = o)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_resid_sd - 1.96 * se, ymax = mean_resid_sd + 1.96 * se), width = 0, linewidth = 0.35,
                  position = position_dodge(0.5)) +
    geom_point(size = 1.1, position = position_dodge(0.5)) + geom_line(linewidth = 0.3, position = position_dodge(0.5)) +
    scale_x_continuous(breaks = 1:10) + scale_colour_manual(values = c("#b2182b", "#ef8a62", "grey35"), name = NULL) +
    labs(x = sub, y = "adjusted symptoms at 15\u201317 (SD)", title = title) + th + theme(legend.position = "bottom")
}
pc <- dec("grad_fit", c("rulebreak", "external", "pfactor"), "c  By decile of gradient fit", "decile (1 = least typical pattern)")
pdd <- dec("grad_slope", c("totprob", "pfactor", "anxdisord"), "d  By decile of gradient slope", "decile (1 = flattest gradient)")

fig <- ((pa / pb) + plot_layout(heights = c(0.55, 1)) | (pc / pdd)) + plot_layout(widths = c(1.1, 1)) +
  plot_annotation(
    title = "Children whose thinning departs from the normative pattern have more symptoms later \u2014 mostly externalising",
    subtitle = paste0("Per child, across 179 parcels: gradient fit = correlation of their thinning with the group map; gradient slope = regression slope on it (1 = typical, < 1 = flatter).\n",
                      "Models adjust for the child's overall thinning rate (the scores correlate with it at 0.17\u20130.29). PLS2 beyond gradient = loading on PLS2 after the gradient (137 parcels)."),
    theme = theme(plot.title = element_text(size = TXT + 2.5, face = "bold"),
                  plot.subtitle = element_text(size = TXT, colour = "grey25", lineheight = 1.15)))
ggsave(file.path(ROOT, "figures", "fig_gradient_scores.png"), fig, width = 10.5, height = 6.8, dpi = 300, bg = "white")
cat("wrote fig_gradient_scores.png\n")
