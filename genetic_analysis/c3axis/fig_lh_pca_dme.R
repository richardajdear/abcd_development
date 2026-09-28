# Run from the repo root: Rscript genetic_analysis/c3axis/fig_lh_pca_dme.R   (reads the tables written by 01/04)
# LH-only PCA vs DME of between-child slope covariance (HCP-MMP, 179 LH parcels).
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})
WD <- "genetic_analysis/c3axis"; AP <- "ahba_pls/data"
mp <- read.csv(file.path(WD, "lh_component_maps.csv"))
LS <- read.delim(file.path(WD, "lh_components_vs_refs.tsv"))
LM <- read.delim(file.path(WD, "lh_method_comparison.tsv"))
CG <- read.delim(file.path(WD, "lh_child_global_coupling.tsv"))
poly <- read.csv(file.path(AP, "hcp_polygons.csv")) |> filter(view %in% c("lateral", "medial"))
TXT <- 8
g <- function(m, c, t) { r <- LS[LS$method == m & LS$comp == c & LS$ref == t, ]; stopifnot(nrow(r) == 1); r }
pf <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
DM <- "DME alpha=1"; RC <- "PCA, row-centred (constant removed)"

brain <- function(v, title) {
  vals <- setNames(mp[[v]], mp$label); l <- max(abs(vals), na.rm = TRUE) * c(-1, 1)
  d <- poly |> mutate(val = unname(vals[label]))
  ggplot(d, aes(x, y, group = interaction(view, label, group, subgroup), fill = val)) +
    geom_polygon(colour = "grey35", linewidth = 0.04) + coord_fixed(expand = FALSE) +
    scale_fill_distiller(palette = "RdBu", limits = l, na.value = "grey85", breaks = pretty_breaks(3), name = NULL) +
    labs(title = title) + theme_void(base_size = TXT) +
    theme(plot.title = element_text(size = TXT - 0.3, face = "bold", hjust = 0.5), legend.position = "bottom",
          legend.key.height = unit(3, "pt"), legend.key.width = unit(18, "pt"), plot.margin = margin(1, 2, 1, 2))
}
row <- wrap_plots(list(
  brain("PCA_c1", sprintf("PCA 1: global (C1 %.2f, C3 %.2f)", g("PCA", "c1", "C1")$rho, g("PCA", "c1", "C3")$rho)),
  brain("DME_c1", sprintf("DME 1 (C1 %.2f)", g(DM, "c1", "C1")$rho)),
  brain("DME_c2", sprintf("DME 2 (C2 %.2f)", g(DM, "c2", "C2")$rho)),
  brain("DME_c5", sprintf("DME 5 (C3 %.2f)", g(DM, "c5", "C3")$rho)),
  brain("DME_c7", sprintf("DME 7 (dCT %.2f)", g(DM, "c7", "dCT_lh")$rho)),
  brain("C3", "AHBA C3"), brain("dCT_lh", "dCT, LH (mm/yr)")), nrow = 1)
pa <- wrap_elements(full = row + plot_annotation(title = "a   Left hemisphere only: selected components and reference maps",
                                                 theme = theme(plot.title = element_text(size = TXT + 1, face = "bold"))))

ML <- c(PCA = "PCA", `PCA, row-centred (constant removed)` = "PCA on row-centred slopes", `DME alpha=1` = "DME (alpha = 1)")
tl <- LS |> mutate(method = factor(ML[method], ML), comp = factor(sub("c", "", comp), 8:1),
                   ref = factor(recode(ref, CT_lh = "CT", dCT_lh = "dCT"), c("CT", "dCT", "PLS1", "PLS2", "C1", "C2", "C3")),
                   sig = p_spin < 0.05)
c5 <- g(DM, "c5", "C3"); p8 <- g("PCA", "c8", "C3")
pb <- ggplot(tl, aes(ref, comp, fill = rho)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho), fontface = ifelse(sig, "bold", "plain"), alpha = sig), size = 2.2) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), name = "Spearman rho") + facet_wrap(~method) +
  labs(x = NULL, y = "component",
       title = sprintf("b   Row-centring makes PCA reproduce DME: C3 emerges at component 5 (DME rho = %.2f, %s) instead of PCA's PC8 (%.2f)",
                       c5$rho, pf(c5$p_spin), p8$rho),
       subtitle = "In plain PCA the global factor (PC1) absorbs part of C1 and C3; once each child's mean slope is removed, C1, C2 and C3 separate") +
  theme_minimal(base_size = TXT) + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
                                         legend.key.width = unit(5, "pt"))

MC <- c(PCA = "PCA", `PCA, global regressed` = "PCA, global mean regressed (per-parcel beta)",
        `PCA, row-centred (constant removed)` = "PCA, global mean subtracted (row-centred)",
        `DME alpha=0` = "DME alpha = 0", `DME alpha=0.5` = "DME alpha = 0.5", `DME alpha=1` = "DME alpha = 1",
        `DME alpha=1, global regressed` = "DME alpha = 1, global mean regressed")
lc <- LM |> filter(method %in% names(MC)) |> mutate(method = factor(MC[method], MC),
        fam = ifelse(grepl("^DME", method), "DME", "PCA"), fam = ifelse(grepl("regressed", method), "global regressed", fam))
rmin <- min(CG$r_DME_vs_rowcentredPCA[1:7])
pc <- ggplot(lc, aes(k, abs(C3), colour = method)) + geom_line() + geom_point(size = 1.1) +
  scale_x_continuous(breaks = 1:10) +
  scale_colour_manual(values = c("#000000", "#999999", "#e41a1c", "#9ecae1", "#4292c6", "#08519c", "#fdae6b"), name = NULL) +
  labs(x = "component index", y = "|rho| with AHBA C3",
       title = sprintf("c   What DME removes is the constant, not the regionally-varying global factor (DME vs row-centred PCA: |r| \u2265 %.3f for components 1-7)", rmin),
       subtitle = "Subtracting each child's mean (red) matches DME at any alpha (blues); regressing the mean out with per-parcel betas (grey, orange) weakens C3 and pushes it to component 7") +
  theme_bw(base_size = TXT) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
                                    legend.position = "right")

cap <- paste(
  "\u2022 Input: standardised child \u00d7 parcel slope BLUPs for the 179 left-hemisphere HCP-MMP parcels (thickness_hcp_70_aa6e91efba82, n = 8,716). AHBA C1-C3/PLS maps cover 137 parcels; CT/dCT are the LH fixed effects.",
  "\u2022 DME: brainspace, normalized-angle kernel, sparsity 0 (Dear et al. 2024 settings). Row-centred: each child's cortex-wide mean of the standardised slopes subtracted uniformly from all parcels.",
  "\u2022 Global mean regressed: each parcel's slope regressed on the child's mean slope with its own beta (removes the whole PC1, including its regional variation).",
  "\u2022 Bold: spin p < 0.05 (5,000 rotations). Signs of components are arbitrary (oriented to the anterior-posterior axis; PCA 1 signed so that high = strong global loading).", sep = "\n")
fig <- pa / pb / pc + plot_layout(heights = c(0.42, 1.15, 0.7)) + plot_annotation(
  title = "Left hemisphere only: DME is PCA on each child's slopes minus their cortex-wide mean, and on that basis C3 is component 5",
  subtitle = "Removing the global mode uniformly releases C1, C2 and C3 as separate components (b). Regressing it out parcel by parcel also removes the C3-weighted part of the global factor (c).",
  caption = cap, theme = theme(plot.title = element_text(size = TXT + 3, face = "bold"), plot.subtitle = element_text(size = TXT + 1),
                               plot.caption = element_text(hjust = 0, lineheight = 1.4, size = TXT - 0.5)))
ggsave(file.path(WD, "fig_lh_pca_dme.png"), fig, width = 13.5, height = 12, dpi = 170, bg = "white")
