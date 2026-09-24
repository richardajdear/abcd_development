#!/usr/bin/env Rscript
# fig_cbcl_explore.R -- EXPLORATORY heatmap of every brain-phenotype x CBCL
# association from 23_cbcl_explore.py. Reads results/cbcl_explore_assoc.tsv
# (coefficients only). Writes figures/fig_cbcl_explore.png.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
A <- read.delim(file.path(ROOT, "results", "cbcl_explore_assoc.tsv"))
A$q <- p.adjust(A$p, "BH")
BR <- c(global_thin = "global thinning", pls2_proj = "PLS2 projection", pls2_top = "PLS2 top-decile mean",
        pls2_topbottom = "PLS2 top \u2212 bottom", c3_proj = "C3 projection",
        pls1_proj = "PLS1 projection (thinning)", pls1_thick_proj = "PLS1 projection (baseline CT)")
OUTS <- c("totprob", "internal", "external", "pfactor", "anxdep", "withdep", "somatic", "social", "thought",
          "attention", "rulebreak", "aggressive", "depress", "anxdisord", "somaticpr", "adhd", "opposit", "conduct")
MOD <- c(base = "baseline (age ~10)", fu = "year 3 (age ~13)", change = "year 3 | baseline", fu_year4 = "year 4 (40% of sample)")
d <- A |> filter(model %in% names(MOD), outcome %in% OUTS) |>
  mutate(br = paste0(BR[brain], ifelse(adj_global, "  | global", "")),
         br = factor(br, rev(unique(c(rbind(BR, paste0(BR, "  | global")))))),
         outcome = factor(outcome, OUTS), model = factor(MOD[model], MOD),
         star = ifelse(q < 0.1, "**", ifelse(p < 0.05, "*", "")))
d <- d |> filter(!is.na(br))
lim <- max(abs(d$y_sd_beta)) * c(-1, 1)
p <- ggplot(d, aes(outcome, br, fill = y_sd_beta)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = star), size = 2.4, vjust = 0.75) +
  scale_fill_distiller(palette = "RdBu", limits = lim, name = "\u03b2 (SD/SD)") +
  facet_wrap(~model, nrow = 1) +
  labs(x = NULL, y = NULL,
       title = "Exploratory: child-level thinning phenotypes vs CBCL (5.1 CBCL, 7.0 imaging; n \u2248 8,000)",
       subtitle = "OLS, site FE + sex + age, family-clustered SE; raw scores log1p. * p < 0.05, ** BH q < 0.1 over all 988 fits. Positive = more thinning in the pattern, more symptoms.") +
  theme_minimal(base_size = 7) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 8.5), plot.subtitle = element_text(size = 6.5, colour = "grey30"))
ggsave(file.path(ROOT, "figures", "fig_cbcl_explore.png"), p, width = 12, height = 4.2, dpi = 250, bg = "white")
cat("wrote fig_cbcl_explore.png\n")
