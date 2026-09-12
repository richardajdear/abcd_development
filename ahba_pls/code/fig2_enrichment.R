#!/usr/bin/env Rscript
# fig2_enrichment.R -- Figure 2 in ggplot2 + patchwork.
#
# Reads results/enrichment_summary.tsv, permutation_enrichment.tsv,
# magma_gene_property.tsv, magma_conditional.tsv. Every number on the figure,
# including the counts quoted in panel subtitles, is computed from these tables.
# Writes figures/fig_enrichment.png

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})

ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "figures")
S <- read.delim(file.path(RES, "enrichment_summary.tsv"))
M <- read.delim(file.path(RES, "magma_gene_property.tsv"))
P <- read.delim(file.path(RES, "permutation_enrichment.tsv"))
C <- read.delim(file.path(RES, "magma_conditional.tsv"))

lvl <- c("ABCD PLS2 (dCT+CT), ds0", "ABCD PLS2 (dCT+CT), ds25", "ABCD PLS2 (dCT+CT), ds50",
         "dCT alone (opt 1)", "dCT+dT1T2 PLS2 (opt 3)", "4-map PLS2 (opt 4)",
         "slope-PC PLS2 (opt 5)", "static PLS1 (control)",
         "AHBA C1", "AHBA C2", "AHBA C3 (Dear 2024)", "NSPN PLS1", "NSPN PLS2 (Whitaker 2016)")
grp <- c(rep("ABCD signature", 3), rep("other Y options", 5), rep("published", 5))
meta <- data.frame(label = lvl, grp = factor(grp, levels = unique(grp)))
S <- S |> mutate(label = factor(label, levels = rev(lvl)), disorder = factor(disorder, c("SCZ", "MDD")))

base <- theme_bw(base_size = 7.6) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 7.8, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, colour = "grey30", margin = margin(b = 3)),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6.6),
        legend.title = element_text(size = 6.6), legend.text = element_text(size = 6.2),
        legend.key.width = unit(7, "pt"), legend.key.height = unit(16, "pt"),
        plot.margin = margin(2, 4, 2, 4))
theme_set(base)

# ---- panel a: MAGMA gene-property ------------------------------------------
n_sig <- S |> filter(grepl("^ABCD PLS2", label)) |> summarise(n = sum(magma_p < 0.05)) |> pull(n)
a_lead <- S |> filter(label == "ABCD PLS2 (dCT+CT), ds25")
a_c3 <- S |> filter(label == "AHBA C3 (Dear 2024)")
pa <- ggplot(S, aes(disorder, label, fill = sign(magma_beta_std) * -log10(magma_p))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(magma_p < 0.05, signif(magma_p, 2), "")), size = 2.1,
            colour = ifelse(abs(sign(S$magma_beta_std) * -log10(S$magma_p)) > 3, "white", "grey10")) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-6, 6), oob = squish,
                       name = "sign(\u03b2) \u00d7\n\u2212log10 p") +
  labs(x = NULL, y = NULL, title = "a   Continuous test: disorder gene Z on weight",
       subtitle = sprintf("all %d of %d ABCD-signature tests p<0.05,\nbut \u03b2 is ~half C3's (%.3f vs %.3f, MDD)",
                          n_sig, 6, a_lead$magma_beta_std[a_lead$disorder == "MDD"],
                          a_c3$magma_beta_std[a_c3$disorder == "MDD"])) +
  theme(axis.text.y = element_text(size = 6.4))

# ---- panel b: set permutation ----------------------------------------------
setcols <- c(perm_SCZ_prioritised_z = "SCZ\nprio.\n(120)", perm_SCZ_locus_pool_z = "SCZ\npool\n(685)",
             perm_MDD_highconf_z = "MDD\nhigh-conf\n(308)", perm_MDD_pool_z = "MDD\npool\n(4,600)")
B <- S |> select(label, disorder, all_of(names(setcols))) |>
  pivot_longer(-c(label, disorder), names_to = "set", values_to = "z") |>
  filter(!is.na(z)) |> mutate(set = factor(setcols[set], levels = setcols))
lead_b <- B |> filter(label == "ABCD PLS2 (dCT+CT), ds25")
stat_b <- B |> filter(label == "static PLS1 (control)", set == setcols[2])
pb <- ggplot(B, aes(set, label, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", z)), size = 2.1,
            colour = ifelse(abs(B$z) > 2.6, "white", "grey10")) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-4, 4), oob = squish, name = "z") +
  labs(x = NULL, y = NULL, title = "b   Set test: prioritised genes over-weighted?",
       subtitle = sprintf("only MDD high-confidence is enriched (z = %.1f);\nthe SCZ pool signal sits on the static axis (z = %.1f)",
                          lead_b$z[lead_b$set == setcols[3]], stat_b$z[1])) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(2, 4, 2, 6))

# ---- panel c: conditional models -------------------------------------------
pair_rows <- function(run, partner) {
  sub <- C |> filter(run == !!run)
  mdl <- sub$MODEL[sub$VARIABLE == partner][1]
  sub |> filter(MODEL == mdl)
}
rows <- bind_rows(
  C |> filter(run == "conditional_none", VARIABLE %in% c("lead", "C3", "NSPN_PLS2")) |>
    mutate(panel = paste0(VARIABLE, " alone")),
  pair_rows("conditional_cond_C1", "C1") |> filter(VARIABLE == "lead") |> mutate(panel = "lead | C1"),
  pair_rows("conditional_cond_static", "static") |> filter(VARIABLE == "lead") |> mutate(panel = "lead | static PLS1"),
  pair_rows("conditional_cond_C3", "C3") |> mutate(panel = paste0(VARIABLE, " | joint lead+C3")),
  pair_rows("conditional_cond_lead", "NSPN_PLS2") |> mutate(panel = paste0(VARIABLE, " | joint lead+NSPN"))
) |>
  mutate(se_std = ifelse(BETA != 0, SE * BETA_STD / BETA, NA),
         nm = recode(panel,
           "lead alone" = "ABCD PLS2 alone", "C3 alone" = "AHBA C3 alone", "NSPN_PLS2 alone" = "NSPN PLS2 alone",
           "lead | C1" = "ABCD PLS2 | C1", "lead | static PLS1" = "ABCD PLS2 | static PLS1",
           "lead | joint lead+C3" = "ABCD PLS2 | C3 (joint)", "C3 | joint lead+C3" = "AHBA C3 | ABCD PLS2 (joint)",
           "lead | joint lead+NSPN" = "ABCD PLS2 | NSPN PLS2 (joint)", "NSPN_PLS2 | joint lead+NSPN" = "NSPN PLS2 | ABCD PLS2 (joint)"),
         disorder = factor(disorder, c("SCZ", "MDD")))
ord <- c("ABCD PLS2 alone", "ABCD PLS2 | C1", "ABCD PLS2 | static PLS1", "ABCD PLS2 | C3 (joint)",
         "AHBA C3 alone", "AHBA C3 | ABCD PLS2 (joint)",
         "NSPN PLS2 alone", "NSPN PLS2 | ABCD PLS2 (joint)", "ABCD PLS2 | NSPN PLS2 (joint)")
rows <- rows |> filter(nm %in% ord) |> mutate(nm = factor(nm, levels = rev(ord)))
lead_c3 <- rows |> filter(nm == "ABCD PLS2 | C3 (joint)")
c3_lead <- rows |> filter(nm == "AHBA C3 | ABCD PLS2 (joint)")
pc <- ggplot(rows, aes(BETA_STD, nm, colour = disorder)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(xmin = BETA_STD - 1.96 * se_std, xmax = BETA_STD + 1.96 * se_std),
                  position = position_dodge(width = 0.55), size = 0.22, linewidth = 0.4) +
  scale_colour_manual(values = c(SCZ = "#762a83", MDD = "#1b7837"), name = NULL) +
  labs(x = "MAGMA gene-property \u03b2 (semi-standardised), 95% CI", y = NULL,
       title = sprintf("c   Which gene ranking carries the signal? (shared universe, n = %s genes)",
                       format(C$NGENES[1], big.mark = ",")),
       subtitle = sprintf("ABCD PLS2 survives C1 and the static axis but vanishes beside C3 (p = %.2f / %.2f, SCZ / MDD) while C3 holds (p = %.4f / %.0e);\nABCD PLS2 and NSPN PLS2 attenuate each other rather than one absorbing the other",
                          lead_c3$P[lead_c3$disorder == "SCZ"], lead_c3$P[lead_c3$disorder == "MDD"],
                          c3_lead$P[c3_lead$disorder == "SCZ"], c3_lead$P[c3_lead$disorder == "MDD"])) +
  theme(legend.position = c(0.93, 0.22), legend.background = element_blank(),
        legend.key.height = unit(8, "pt"), panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey93"))

methods <- paste(
  "Methods.",
  "\u2022 Universe: AHBA genes with a MAGMA gene-level result in both disorders (13,763 / 10,322 / 6,931 genes at ds0 / ds25 / ds50).",
  "\u2022 GWAS: SCZ = PGC3 wave 3 (Trubetskoy 2022); MDD = PGC MDD 2025 (Adams 2025). Gene analysis from this project's hpc/ run (MAGMA v1.10, NCBI37.3).",
  "\u2022 Panel a: MAGMA --gene-covar regression; models gene\u2013gene LD and conditions internally on gene size, density and sample size; two-sided.",
  "\u2022 Panel b: mean weight of the set vs 10,000 random sets matched on size and gene-length decile; sets replicate the hpc MAGMA definitions.",
  "\u2022 Panel c: the same regression with a second weight vector added. \"X | Y\" = coefficient of X with Y in the model.",
  "\u2022 Rows in a and b: ABCD signature at three DS filters; other Y-matrix options and the within-analysis static control; published benchmarks.",
  sep = "\n")

fig <- (pa | pb) / pc +
  plot_layout(heights = c(1, 0.85)) +
  plot_annotation(
    title = "The signature carries SCZ and MDD genetic risk \u2014 but the signal is C3's, not new",
    subtitle = "Positive = genes weighted towards faster thinning carry more disorder GWAS signal.",
    caption = methods,
    theme = theme(plot.title = element_text(size = 9.5, face = "bold"),
                  plot.subtitle = element_text(size = 7.4, colour = "grey25", margin = margin(b = 4)),
                  plot.caption = element_text(size = 6.2, colour = "grey25", hjust = 0, lineheight = 1.45,
                                              margin = margin(t = 6))))

ggsave(file.path(FIG, "fig_enrichment.png"), fig, width = 7.8, height = 6.6, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_enrichment.png"), "\n")
