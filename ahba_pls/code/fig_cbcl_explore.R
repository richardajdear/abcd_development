#!/usr/bin/env Rscript
# fig_cbcl_explore.R [70|51] -- EXPLORATORY heatmap of every brain-phenotype x
# CBCL association from 23_cbcl_explore.py (OLS models only; clinical and KSADS
# logistic fits are in the table). Reads results/cbcl_explore[_51]_assoc.tsv
# (coefficients only). Writes figures/fig_cbcl_explore[_51].png.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
SRC <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(SRC)) SRC <- "70"
TAG <- if (SRC == "70") "cbcl_explore" else "cbcl_explore_51"
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
A <- read.delim(file.path(ROOT, "results", paste0(TAG, "_assoc.tsv")))
A$adj_global <- as.logical(A$adj_global)          # pandas writes True/False
A$q <- p.adjust(A$p, "BH")
BR <- c(global_thin = "global thinning", pls2_proj = "PLS2 projection", pls2_top = "PLS2 top-decile mean",
        pls2_topbottom = "PLS2 top \u2212 bottom", c3_proj = "C3 projection",
        pls1_proj = "PLS1 projection (thinning)", pls1_thick_proj = "PLS1 projection (baseline CT)")
OUTS <- c("totprob", "internal", "external", "pfactor", "anxdep", "withdep", "somatic", "social", "thought",
          "attention", "rulebreak", "aggressive", "depress", "anxdisord", "somaticpr", "adhd", "opposit", "conduct")
MOD <- if (SRC == "70") {
  c(base = "baseline (age ~10)", y3 = "year 3 (age ~13)", late = "years 5\u20137 (age ~15\u201317)",
    change = "years 5\u20137 | baseline", trajectory = "symptom slope, all 8 waves")
} else {
  c(base = "baseline (age ~10)", y3 = "year 3 (age ~13)", change = "year 3 | baseline", y4 = "year 4 (40% of sample)")
}
d <- A |> filter(model %in% names(MOD), outcome %in% OUTS) |>
  mutate(br = paste0(BR[brain], ifelse(adj_global, "  | global", "")),
         br = factor(br, rev(unique(c(rbind(BR, paste0(BR, "  | global")))))),
         outcome = factor(outcome, OUTS), model = factor(MOD[model], MOD),
         star = ifelse(q < 0.05, "**", ifelse(p < 0.05, "*", ""))) |>
  filter(!is.na(br))
lim <- max(abs(d$y_sd_beta), na.rm = TRUE) * c(-1, 1)
p <- ggplot(d, aes(outcome, br, fill = y_sd_beta)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = star), size = 2.4, vjust = 0.75) +
  scale_fill_distiller(palette = "RdBu", limits = lim, name = "\u03b2 (SD/SD)") +
  facet_wrap(~model, nrow = 1) +
  labs(x = NULL, y = NULL,
       title = sprintf("Exploratory: child-level thinning phenotypes vs CBCL (%s CBCL, 7.0 imaging; n \u2248 8,000)",
                       if (SRC == "70") "7.0" else "5.1"),
       subtitle = sprintf("OLS, site FE + sex + age, family-clustered SE; raw scores log1p. * p < 0.05, ** BH q < 0.05 over all %s fits (incl. logistic). Positive = more thinning in the pattern, more symptoms.",
                          format(nrow(A), big.mark = ","))) +
  theme_minimal(base_size = 7) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 8.5), plot.subtitle = element_text(size = 6.5, colour = "grey30"))
out <- file.path(ROOT, "figures", paste0(sub("cbcl_explore", "fig_cbcl_explore", TAG), ".png"))
ggsave(out, p, width = if (SRC == "70") 14 else 12, height = 4.2, dpi = 250, bg = "white")
cat("wrote", out, "\n")
