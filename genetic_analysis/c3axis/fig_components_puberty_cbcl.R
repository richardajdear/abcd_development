# Run from the repo root: Rscript genetic_analysis/c3axis/fig_components_puberty_cbcl.R   (reads the tables written by 01/04)
suppressMessages({library(ggplot2); library(dplyr); library(patchwork); library(scales)})
AS <- read.delim("genetic_analysis/c3axis/table_components_puberty_cbcl.tsv"); SH <- read.delim("genetic_analysis/c3axis/table_component_split_half.tsv")
TXT <- 8
lvl <- unique(AS$label)
OL <- c(pds_p_timing_mean="timing (parent, mean)", pds_p_timing_base="timing (parent, baseline)", pds_p_tempo="tempo (parent)",
        pds_y_timing_mean="timing (youth, mean)", pds_y_tempo="tempo (youth)")
AS <- AS |> mutate(label = factor(label, rev(lvl)),
  z = ifelse(family == "cbcl_onset", log(pmax(effect, 1e-9)), effect),
  out = ifelse(outcome %in% names(OL), OL[outcome], sub("^cbcl_base_|^onset_", "", outcome)),
  star = ifelse(q_fdr_family < 0.05, "*", ifelse(p < 0.05, "\u00b7", "")))
tile <- function(fm, ttl, sub, leg) {
  d <- AS |> filter(family == fm); lim <- max(abs(d$z)) * c(-1, 1)
  ggplot(d, aes(out, label, fill = z)) + geom_tile(colour = "white") +
    geom_text(aes(label = paste0(sprintf("%.2f", effect), star)), size = 2.1) +
    scale_fill_distiller(palette = "RdBu", limits = lim, name = leg) +
    labs(x = NULL, y = NULL, title = ttl, subtitle = sub) + theme_minimal(base_size = TXT) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(), plot.title = element_text(face = "bold"),
          legend.key.width = unit(5, "pt"))
}
cnt <- function(fm) with(AS[AS$family == fm, ], sprintf("%d of %d tests p < 0.05, %d survive FDR", sum(p < 0.05), length(p), sum(q_fdr_family < 0.05)))
pa <- tile("puberty", "a   Puberty: small effects on C1 and the global factor; no FDR-surviving effect on the C3 axis (rc5)", cnt("puberty"), "\u03b2 (SD/SD)")
pb <- tile("cbcl_baseline", "b   CBCL at baseline: nothing beyond chance", cnt("cbcl_baseline"), "\u03b2 (SD/SD)")
pc <- tile("cbcl_onset", "c   CBCL onset (ages ~15-17): global factor leans internalising; nothing survives FDR", cnt("cbcl_onset"), "log OR")
SH <- SH |> mutate(decomposition = factor(decomposition, unique(decomposition)))
pd <- ggplot(SH, aes(factor(comp), same_index_median, colour = decomposition, group = decomposition)) +
  geom_linerange(aes(ymin = same_index_p05, ymax = same_index_median), position = position_dodge(0.5), alpha = 0.5) +
  geom_point(position = position_dodge(0.5), size = 1.6) +
  geom_text(aes(y = 0.02, label = sprintf("%.1f", var_explained_pct)), position = position_dodge(0.5), size = 2, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 1)) + scale_colour_manual(values = c("#636363", "#4292c6", "#d95f02"), name = NULL) +
  labs(x = "component", y = "split-half |r| of loadings (median; line to 5th pct)",
       title = "d   Split-half reliability of the component maps (100 random halves); numbers = % variance explained",
       subtitle = with(SH[SH$decomposition == "LH only, row-centred PCA (= DME)" & SH$comp %in% 4:5, ],
                       sprintf("Row-centred 4 and 5 swap order in %.0f%% of splits (%.2f%% vs %.2f%% variance); their shared plane and its C3-aligned direction are stable",
                               100 * (1 - index_stable_frac[1]), var_explained_pct[1], var_explained_pct[2]))) +
  theme_bw(base_size = TXT) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), legend.position = "bottom")
cap <- paste(
  "\u2022 Components: plain PCA of standardised HCP-MMP slope BLUPs (358 parcels, PC1-5) and PCA of LEFT-hemisphere row-centred slopes (rc1-7, identical to DME); child scores standardised. n = 8,716.",
  "\u2022 Signs: PC1 high = faster global thinning; PC2/rc1 +C1, PC3/rc2 +C2, rc5 +C3; rc4/rc7 high = amplified normative dCT contrast; others A-P axis.",
  "\u2022 Puberty: parent/youth Pubertal Development Scale (ph_p_pds, ph_y_pds); timing = within-sex, within-wave residual on cubic age; tempo = per-child PDS slope (>= 3 waves), adjusted for baseline timing.",
  "\u2022 CBCL (mh_p_cbcl): baseline log1p raw scores, p-factor = PC1 of 8 syndromes; onset = below threshold at ses-00A and at/above at any of 05A-07A vs below at every wave (T >= 60 broadband, >= 65 syndrome/DSM; p-factor top decile).",
  "\u2022 Models: outcome ~ component + sex + site + age at first scan + follow-up span + n visits; OLS or logistic, family-clustered SE. * FDR q < 0.05 within panel; \u00b7 p < 0.05.", sep = "\n")
fig <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1, 1)) + plot_annotation(
  title = "Slope components vs puberty and CBCL: puberty tracks the C1 axis and the global factor; the C3 axis is unrelated to either",
  subtitle = "Effects are small (|\u03b2| \u2264 0.10). The C3 component (rc5) shows no FDR-surviving association with puberty or symptoms; the global factor (PC1) leans towards internalising onset.",
  caption = cap, theme = theme(plot.title = element_text(size = TXT + 3, face = "bold"), plot.caption = element_text(hjust = 0, lineheight = 1.4, size = TXT - 0.8)))
ggsave("genetic_analysis/c3axis/fig_components_puberty_cbcl.png", fig, width = 15, height = 11, dpi = 160, bg = "white")
