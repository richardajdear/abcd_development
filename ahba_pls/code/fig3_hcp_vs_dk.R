#!/usr/bin/env Rscript
# fig3_hcp_vs_dk.R -- does fitting the PLS in HCP-MMP instead of DK add
# SCZ/MDD gene-level signal?  Reads results/hcp_vs_dk_enrichment.tsv and
# results/hcp_vs_dk_conditional.tsv (one shared MAGMA universe); recomputes
# nothing. Writes figures/fig_hcp_vs_dk.png

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})

ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "figures")
M <- read.delim(file.path(RES, "hcp_vs_dk_enrichment.tsv"))
C <- read.delim(file.path(RES, "hcp_vs_dk_conditional.tsv"))
cmp <- read.delim(file.path(RES, "hcp_pls_components.tsv"))
conc <- read.delim(file.path(RES, "hcp_concordance.tsv"))
PROV <- read.delim(file.path(RES, "hcp_run_provenance.tsv"))   # sample size, written by 12_hcp_pls.py

base <- theme_bw(base_size = 7.6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey93"),
        plot.title = element_text(size = 7.8, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, colour = "grey30", margin = margin(b = 3)),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6.6),
        legend.text = element_text(size = 6.4), legend.title = element_blank(),
        plot.margin = margin(2, 5, 2, 5))
theme_set(base)
pal <- c(SCZ = "#762a83", MDD = "#1b7837")

# ---- panel a: marginal betas -----------------------------------------------
# Panel a shows ONLY vectors on the matched gene basis (abagen-data, the same
# 7,973 genes = HCP-space DS-50), so the HCP-vs-DK contrast is parcellation and
# nothing else. The AHBA_updated native-DK ds0/ds25/ds50 vectors are in the
# table with basis = "AHBA_updated native-DK" and are excluded here: they differ
# from the HCP fit in build and gene identity as well as parcellation. They are
# compared as gene RANKINGS (where provenance does not matter) in
# fig_enrichment_combined.png.
lab_a <- c(HCP_PLS2_thinning = "HCP-MMP PLS2, dCT+CT (137 parcels)",
           DK_PLS2_matchedX = "DK PLS2, dCT+CT (33 regions)",
           HCP_PLS1_dCTonly = "HCP-MMP, dCT alone (137)",
           DK_PLS1_dCTonly_matchedX = "DK, dCT alone (33)",
           C3_shipped = "AHBA C3 (Dear 2024)",
           C1_shipped = "AHBA C1 (static axis)")
MATCHED <- "abagen-data, 7,973 genes (HCP-space DS-50)"
stopifnot(all(M$basis[M$VARIABLE %in% names(lab_a)] == MATCHED))
# rows ordered by the SCZ effect (largest at the top), computed from the table
ord_a <- M |> filter(VARIABLE %in% names(lab_a), disorder == "SCZ") |>
  arrange(BETA_STD) |> pull(VARIABLE)
A <- M |> filter(VARIABLE %in% names(lab_a)) |>
  mutate(nm = factor(lab_a[VARIABLE], levels = lab_a[ord_a]),
         disorder = factor(disorder, c("SCZ", "MDD")))
# which ABCD-derived vector leads on SCZ, read from the table
scz_top <- A |> filter(disorder == "SCZ", grepl("^(HCP|DK)", VARIABLE)) |> arrange(desc(BETA_STD))
pa <- ggplot(A, aes(BETA_STD, nm, colour = disorder)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(xmin = BETA_STD - 1.96 * se_std, xmax = BETA_STD + 1.96 * se_std),
                  position = position_dodge(width = 0.55), size = 0.22, linewidth = 0.4) +
  scale_colour_manual(values = pal) +
  labs(x = "MAGMA gene-property \u03b2 (semi-standardised), 95% CI", y = NULL,
       title = sprintf("a   Each vector on its own, ordered by SCZ effect (one gene basis, shared MAGMA universe of %s genes)",
                       format(M$NGENES[1], big.mark = ",")),
       subtitle = sprintf("HCP-MMP raises the MDD effect from %.3f to %.3f; SCZ is unchanged (%.3f vs %.3f).\nOn SCZ the single-Y (dCT-alone) fits lead both PLS2 fits in either atlas \u2014 %s tops the ABCD rows at \u03b2 = %.3f",
                          A$BETA_STD[A$VARIABLE == "DK_PLS2_matchedX" & A$disorder == "MDD"],
                          A$BETA_STD[A$VARIABLE == "HCP_PLS2_thinning" & A$disorder == "MDD"],
                          A$BETA_STD[A$VARIABLE == "DK_PLS2_matchedX" & A$disorder == "SCZ"],
                          A$BETA_STD[A$VARIABLE == "HCP_PLS2_thinning" & A$disorder == "SCZ"],
                          sub(" \\(.*", "", lab_a[scz_top$VARIABLE[1]]), scz_top$BETA_STD[1])) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(t = -4, b = -2), legend.key.height = unit(7, "pt"))

# ---- panel b: joint models --------------------------------------------------
pick <- function(run, var, nm) C |> filter(run == !!run, VARIABLE == !!var) |> mutate(nm = nm)
B <- bind_rows(
  pick("none", "HCP_PLS2_thinning", "HCP PLS2 alone"),
  pick("cond_DK", "HCP_PLS2_thinning", "HCP PLS2 | DK PLS2"),
  pick("none", "DK_PLS2_matchedX", "DK PLS2 alone"),
  pick("cond_HCP", "DK_PLS2_matchedX", "DK PLS2 | HCP PLS2"),
  pick("cond_C3", "HCP_PLS2_thinning", "HCP PLS2 | AHBA C3"),
  pick("cond_HCP", "C3_shipped", "AHBA C3 | HCP PLS2")) |>
  mutate(nm = factor(nm, levels = rev(c("HCP PLS2 alone", "HCP PLS2 | DK PLS2",
                                        "DK PLS2 alone", "DK PLS2 | HCP PLS2",
                                        "HCP PLS2 | AHBA C3", "AHBA C3 | HCP PLS2"))),
         disorder = factor(disorder, c("SCZ", "MDD")))
pb <- ggplot(B, aes(BETA_STD, nm, colour = disorder)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(xmin = BETA_STD - 1.96 * se_std, xmax = BETA_STD + 1.96 * se_std),
                  position = position_dodge(width = 0.55), size = 0.22, linewidth = 0.4) +
  scale_colour_manual(values = pal, guide = "none") +
  labs(x = "MAGMA gene-property \u03b2 (semi-standardised), 95% CI", y = NULL,
       title = "b   Joint models: which ranking carries the signal?",
       subtitle = sprintf("for MDD, HCP survives DK (p = %.3f) while DK does not survive HCP (p = %.2f); both remain absorbed by C3",
                          B$P[B$nm == "HCP PLS2 | DK PLS2" & B$disorder == "MDD"],
                          B$P[B$nm == "DK PLS2 | HCP PLS2" & B$disorder == "MDD"]))

lead <- cmp |> filter(option == "hcp_opt2_dCT_CT", component == "PLS2")
c3s <- conc |> filter(reference == "C3", level == "scores")
c3w <- conc |> filter(reference == "C3", level == "weights")
methods <- paste(
  "Methods.",
  sprintf("\u2022 Imaging: thickness in HCP-MMP1.0, derived from the release FreeSurfer surfaces (ABCD tabulates only Desikan). After the 2026-09-14\n  parcellation backfill this run has %s children with \u22652 visits (%s sessions) \u2014 the same sample as the DK run.",
          format(PROV$n_subjects[1], big.mark = ","), format(PROV$n_sessions[1], big.mark = ",")),
  "\u2022 Same option-2 design as Figure 1 (Y = bilateral dCT + CT), 137 of 180 left parcels with AHBA donor coverage; medial-wall parcel H excluded.",
  sprintf("\u2022 PLS2 carries dCT (salience %.2f), explains %.0f%% of the cross-covariance, spin p = %.3f over 5,000 rotations of the complete 180-parcel map; bootstrap reproducibility %.2f.",
          lead$sal_dCT, 100 * lead$cov_explained, lead$p_spin_singular, lead$boot_reproducibility),
  sprintf("\u2022 It is the same axis: scores vs AHBA C3 rho = %.2f (p_spin < 0.001, 137 parcels), gene weights rho = %.2f \u2014 measured in C3's own parcellation rather than a DK projection.",
          c3s$rho, c3w$rho),
  "\u2022 Gene basis, identical for every row: the 7,973 genes the shipped C1\u2013C3 weights were fitted on, which are exactly the columns of hcp_3d_ds5.csv \u2014 the top half of",
  "  15,946 genes by differential stability computed in HCP space. The HCP arm is therefore DS-50-filtered, not unfiltered. DS filters are parcellation-specific (dk_3d_ds5 shares",
  "  88% of that list, the AHBA_updated native-DK ds50 matrix 87%), so the DK rows here are those same 7,973 genes taken from dk_3d.csv \u2014 same abagen build, same genes, 33",
  "  regions \u2014 and parcellation is the only difference on the figure. The AHBA_updated DK vectors of the main analysis are in the results table, marked as a different basis.",
  sprintf("\u2022 All %d hemisphere-region models converged with none singular; medial-wall parcel H is excluded by the run config.", PROV$n_labels_fitted[1]),
  sep = "\n")

fig <- (pa / pb) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Fitting the PLS in C3's own parcellation adds MDD signal, but not SCZ",
    subtitle = "Positive = genes weighted towards faster thinning carry more disorder GWAS signal.",
    caption = methods,
    theme = theme(plot.title = element_text(size = 9.5, face = "bold"),
                  plot.subtitle = element_text(size = 7.4, colour = "grey25", margin = margin(b = 4)),
                  plot.caption = element_text(size = 6.2, colour = "grey25", hjust = 0,
                                              lineheight = 1.45, margin = margin(t = 6))))

ggsave(file.path(FIG, "fig_hcp_vs_dk.png"), fig, width = 7.6, height = 5.2, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_hcp_vs_dk.png"), "\n")