#!/usr/bin/env Rscript
# fig2_enrichment_combined.R -- Figure 2: which gene rankings carry SCZ/MDD
# genetic risk, and which survive with another in the model.  MAGMA only.
#
# Reads results/magma_all_marginal.tsv and magma_all_joint.tsv (one shared
# universe, written by code/14_magma_all_options.py) and fits nothing.
# Writes figures/fig_enrichment_combined.png

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})

ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "figures")
M <- read.delim(file.path(RES, "magma_all_marginal.tsv"))
J <- read.delim(file.path(RES, "magma_all_joint.tsv"))
# the universe is shared across vectors WITHIN a disorder; the two GWAS cover
# slightly different gene sets, so assert one N per disorder and quote both.
ng <- M |> group_by(disorder) |> summarise(n = n_distinct(NGENES), N = max(NGENES), .groups = "drop")
ngj <- J |> group_by(disorder) |> summarise(n = n_distinct(NGENES), N = max(NGENES), .groups = "drop")
stopifnot(all(ng$n == 1), all(ngj$n == 1), all(ng$N == ngj$N[match(ng$disorder, ngj$disorder)]))
NG <- paste(sprintf("%s (%s)", format(ng$N, big.mark = ","), ng$disorder), collapse = " / ")

pal <- c(SCZ = "#762a83", MDD = "#1b7837")
base <- theme_bw(base_size = 7.6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey93"),
        plot.title = element_text(size = 7.9, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, colour = "grey30", margin = margin(b = 3)),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6.7),
        legend.text = element_text(size = 6.5), legend.title = element_blank(),
        plot.margin = margin(2, 5, 2, 5))
theme_set(base)

lab <- c(ABCD_PLS2_HCP = "ABCD PLS2, HCP-MMP (137 parcels)",
         ABCD_PLS2_DK = "ABCD PLS2, DK (33 regions)",
         ABCD_dCT_HCP = "ABCD dCT alone, HCP-MMP",
         ABCD_dCT_DK = "ABCD dCT alone, DK",
         ABCD_dCT_dT1T2_DK = "ABCD dCT+dT1T2 PLS2, DK",
         AHBA_C3 = "AHBA C3 (Dear 2024)",
         NSPN_PLS2 = "NSPN PLS2 (Whitaker 2016)",
         AHBA_C1 = "AHBA C1 (static gradient)")
grp <- c(rep("ABCD-derived", 5), rep("published", 3))

A <- M |> filter(VARIABLE %in% names(lab)) |>
  mutate(nm = factor(lab[VARIABLE], levels = rev(lab)),
         disorder = factor(disorder, c("SCZ", "MDD")),
         sig = ifelse(P < 0.05, "p < 0.05", "n.s."))
# the claim is computed, not asserted: which vectors reach p<0.05 for BOTH
# disorders, and which do not (NSPN PLS2 is SCZ-only; C1 fails both).
both <- A |> group_by(VARIABLE) |> summarise(n_sig = sum(P < 0.05), .groups = "drop")
abcd_v <- names(lab)[grepl("^ABCD", names(lab))]
n_abcd_both <- sum(both$n_sig[both$VARIABLE %in% abcd_v] == 2)
part <- both |> filter(n_sig < 2) |> pull(VARIABLE)
part_txt <- paste(sapply(part, function(v) {
  r <- A |> filter(VARIABLE == v)
  hit <- as.character(r$disorder[r$P < 0.05])
  sprintf("%s (%s)", sub(" \\(.*", "", lab[v]),
          if (length(hit)) paste(hit, "only") else "neither")
}), collapse = ", ")
subtitle_a <- sprintf("All %d ABCD-derived rankings are enriched for both disorders and AHBA C3 is the strongest (MDD \u03b2 = %.3f).\nNot significant for both: %s",
                      n_abcd_both, A$BETA_STD[A$VARIABLE == "AHBA_C3" & A$disorder == "MDD"], part_txt)

pa <- ggplot(A, aes(BETA_STD, nm, colour = disorder)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_hline(yintercept = 3.5, colour = "grey80", linewidth = 0.3, linetype = "22") +
  geom_pointrange(aes(xmin = BETA_STD - 1.96 * se_std, xmax = BETA_STD + 1.96 * se_std,
                      alpha = sig),
                  position = position_dodge(width = 0.6), size = 0.22, linewidth = 0.42) +
  scale_colour_manual(values = pal) +
  scale_alpha_manual(values = c("p < 0.05" = 1, "n.s." = 0.35), guide = "none") +
  labs(x = "MAGMA gene-property \u03b2 (semi-standardised), 95% CI", y = NULL,
       title = "a   Does the gene ranking carry disorder risk? (each vector alone)",
       subtitle = subtitle_a) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(t = -4, b = -2), legend.key.height = unit(7, "pt"))

# ---- panel b: joint models, both coefficients side by side -----------------
short <- c(ABCD_PLS2_HCP = "ABCD PLS2 HCP", ABCD_PLS2_DK = "ABCD PLS2 DK",
           AHBA_C3 = "AHBA C3", AHBA_C1 = "AHBA C1", NSPN_PLS2 = "NSPN PLS2")
pair_lab <- c("ABCD_PLS2_HCP + ABCD_PLS2_DK" = "HCP vs DK\n(does the finer atlas add?)",
              "ABCD_PLS2_HCP + AHBA_C1" = "HCP vs C1\n(is it the static gradient?)",
              "ABCD_PLS2_HCP + NSPN_PLS2" = "HCP vs NSPN PLS2",
              "ABCD_PLS2_DK + NSPN_PLS2" = "DK vs NSPN PLS2",
              "ABCD_PLS2_HCP + AHBA_C3" = "HCP vs C3",
              "ABCD_PLS2_DK + AHBA_C3" = "DK vs C3",
              "NSPN_PLS2 + AHBA_C3" = "NSPN PLS2 vs C3")
B <- J |> filter(pair %in% names(pair_lab)) |>
  mutate(pr = factor(pair_lab[pair], levels = rev(pair_lab)),
         vec = short[VARIABLE],
         disorder = factor(disorder, c("SCZ", "MDD")),
         sig = ifelse(P < 0.05, "p < 0.05", "n.s."))
pb <- ggplot(B, aes(BETA_STD, pr, colour = disorder)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_pointrange(aes(xmin = BETA_STD - 1.96 * se_std, xmax = BETA_STD + 1.96 * se_std,
                      shape = vec, alpha = sig),
                  position = position_dodge(width = 0.72), size = 0.2, linewidth = 0.4) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_alpha_manual(values = c("p < 0.05" = 1, "n.s." = 0.35), guide = "none") +
  scale_shape_manual(values = c(16, 21, 17, 2, 15), name = NULL) +
  facet_wrap(~ disorder, nrow = 1) +
  labs(x = "MAGMA gene-property \u03b2 in the JOINT model, 95% CI", y = NULL,
       title = "b   Which of a pair keeps its signal when both are in the model? (filled = p < 0.05)",
       subtitle = "For MDD the HCP ranking beats DK and NSPN PLS2 outright; for SCZ the pairs merely attenuate each other. Both are still beaten by C3.") +
  theme(strip.background = element_rect(fill = "grey96", colour = NA),
        strip.text = element_text(size = 6.9, face = "bold"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(t = -4, b = -2), legend.key.height = unit(7, "pt"),
        panel.spacing = unit(6, "pt"))

hcp_dk <- B |> filter(pair == "ABCD_PLS2_HCP + ABCD_PLS2_DK", disorder == "MDD")
methods <- paste(
  "Methods.",
  sprintf("\u2022 One MAGMA gene-property run per panel row, all on the SAME universe (n = %s genes with a gene-level result for both disorders and a weight in every vector).", NG),
  "\u2022 GWAS: SCZ = PGC3 wave 3 (Trubetskoy 2022); MDD = PGC MDD 2025 (Adams 2025). Gene analysis from this project's hpc/ run (MAGMA v1.10, NCBI37.3).",
  "\u2022 \u03b2 is the change in disorder gene Z per SD of gene weight; MAGMA models gene-gene LD and conditions internally on gene size, density and sample size; two-sided.",
  "\u2022 Panel b fits both covariates in one model, so each pair shows two coefficients: the vector that keeps a non-zero \u03b2 carries the signal the other only shares.",
  "\u2022 Sign convention: positive = genes expressed more where adolescent thinning is faster (C1/C3 and NSPN PLS2 keep their published sign).",
  sep = "\n")

fig <- (pa / pb) + plot_layout(heights = c(1, 1.25)) +
  plot_annotation(
    title = "SCZ and MDD enrichment of every candidate gene ranking, and which survives a head-to-head",
    subtitle = sprintf("The HCP-MMP thinning signature is the strongest ABCD-derived ranking for MDD (joint \u03b2 = %.3f vs %.3f for DK, p = %.3f vs %.2f); AHBA C3 still leads overall.",
                       hcp_dk$BETA_STD[hcp_dk$VARIABLE == "ABCD_PLS2_HCP"],
                       hcp_dk$BETA_STD[hcp_dk$VARIABLE == "ABCD_PLS2_DK"],
                       hcp_dk$P[hcp_dk$VARIABLE == "ABCD_PLS2_HCP"],
                       hcp_dk$P[hcp_dk$VARIABLE == "ABCD_PLS2_DK"]),
    caption = methods,
    theme = theme(plot.title = element_text(size = 9.4, face = "bold"),
                  plot.subtitle = element_text(size = 7.3, colour = "grey25", margin = margin(b = 4)),
                  plot.caption = element_text(size = 6.2, colour = "grey25", hjust = 0,
                                              lineheight = 1.45, margin = margin(t = 6))))

ggsave(file.path(FIG, "fig_enrichment_combined.png"), fig, width = 7.8, height = 6.4, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_enrichment_combined.png"), "\n")
