#!/usr/bin/env Rscript
# fig1_signature_both.R -- Figure 1: the thinning signature derived in BOTH
# parcellations (DK, 33 regions; HCP-MMP, 137 parcels), its agreement with the
# two published signatures, and the cell-class profile of three rankings.
#
# Reads only saved tables and fits nothing; every printed statistic comes from
# the table plotted.
#   data/y_maps_bilateral_34.csv, data/dk_polygons.csv, data/hcp_polygons.csv
#   results/lead_signature_scores.csv, lead_signature_weights.tsv
#   results/hcp_pls_scores.csv, hcp_pls_weights.tsv, hcp_pls_components.tsv
#   results/concordance_scores.tsv, concordance_weights.tsv, hcp_concordance.tsv
#   results/pls_components.tsv, celltype_compare.tsv
#   data/reference/{nspn_dk_maps_bilateral_34,ahba_c123_*}
# Writes figures/fig_signature_both.png
#
# HCP polygons come from ggsegGlasser (installed in ahba_pls/.Rlib; see README).

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})

ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); DATA <- file.path(ROOT, "data")
REF <- file.path(DATA, "reference"); FIG <- file.path(ROOT, "figures")
LEAD_OPT <- "opt2_dCT_CT"; LEAD_COMP <- "PLS2"; LEAD_DS <- "ds25"

y34    <- read.csv(file.path(DATA, "y_maps_bilateral_34.csv"), row.names = 1)
dkS    <- read.csv(file.path(RES, "lead_signature_scores.csv"))
dkW    <- read.delim(file.path(RES, "lead_signature_weights.tsv"))
hcpS   <- read.csv(file.path(RES, "hcp_pls_scores.csv"), row.names = 1)
hcpW   <- read.delim(file.path(RES, "hcp_pls_weights.tsv"))
conc   <- read.delim(file.path(RES, "concordance_scores.tsv"))
wconc  <- read.delim(file.path(RES, "concordance_weights.tsv"))
hconc  <- read.delim(file.path(RES, "hcp_concordance.tsv"))
comps  <- read.delim(file.path(RES, "pls_components.tsv"))
hcomps <- read.delim(file.path(RES, "hcp_pls_components.tsv"))
cts    <- read.delim(file.path(RES, "celltype_compare.tsv"))
c3w    <- read.csv(file.path(REF, "ahba_c123_gene_weights.csv"))
nspn34 <- read.csv(file.path(REF, "nspn_dk_maps_bilateral_34.csv"), row.names = 1)
c3dk   <- read.csv(file.path(REF, sprintf("ahba_c123_scores_recomputed_%s.csv", LEAD_DS)), row.names = 1)
dkpoly <- read.csv(file.path(DATA, "dk_polygons.csv"))
hcppoly <- read.csv(file.path(DATA, "hcp_polygons.csv"))

# concordance rho for the DK fit is stored in pls.py's raw orientation (dCT
# salience positive); the figure is in the thinning orientation, so flip it.
cc <- function(ref) {
  r <- conc |> filter(option == LEAD_OPT, ds == LEAD_DS, component == LEAD_COMP, reference == ref)
  list(rho = -r$rho[1], p = r$p_spin[1])
}
hc <- function(ref, lvl) {
  r <- hconc |> filter(reference == ref, level == lvl); list(rho = r$rho[1], p = r$p_spin[1])
}
wcc <- function(ref) {
  r <- wconc |> filter(option == LEAD_OPT, ds == LEAD_DS, component == LEAD_COMP, reference == ref)
  list(rho = -r$rho[1], n = r$n_genes[1])
}
pf <- function(p) if (is.na(p)) "" else if (p < 1e-3) "p_spin<0.001" else sprintf("p_spin=%.3f", p)

# ------------------------------------------------------------------ brains ----
dk_lh <- dkpoly |> filter(hemi == "left", view %in% c("lateral", "medial")) |>
  mutate(grp = interaction(label, view, group, subgroup, drop = TRUE),
         region = sub("^lh_", "", label))
hcp_lh <- hcppoly |> mutate(grp = interaction(label, view, group, subgroup, drop = TRUE))

brain <- function(poly, values, title, diverging = TRUE) {
  d <- data.frame(region = sub("^lh_", "", names(values)), value = as.numeric(values))
  pd <- poly |> left_join(d, by = "region")
  sc <- if (diverging) {
    lim <- max(abs(pd$value), na.rm = TRUE)
    scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-lim, lim), na.value = "grey78")
  } else scale_fill_viridis_c(na.value = "grey78")
  ggplot(pd, aes(x, y, group = grp, fill = value)) +
    geom_polygon(colour = "grey40", linewidth = 0.05) + sc +
    coord_fixed(expand = FALSE) + labs(title = title) + theme_void(base_size = 7) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6.9, hjust = 0, margin = margin(b = 0.5, l = 1)),
          plot.margin = margin(0.5, 1, 0.5, 1))
}
blank_note <- function(txt) ggplot() + annotate("text", 0, 0, label = txt, size = 2.15,
                                                colour = "grey45", lineheight = 1.15) +
  theme_void() + theme(plot.margin = margin(0.5, 1, 0.5, 1))

dk_scores <- setNames(dkS$thinning_scores_gene_side, dkS$label)
hcp_scores <- setNames(hcpS$thinning_score, rownames(hcpS))
hcpY <- read.csv(file.path(RES, "hcp_y_maps_180.csv"), row.names = 1)
hcp_dCT <- setNames(hcpY$dCT, rownames(hcpY))   # all parcels, not just AHBA-covered
hcp_c3 <- setNames(hcpS$C3, rownames(hcpS))

b1 <- brain(dk_lh, setNames(y34$dCT, rownames(y34)), "thinning rate (dCT)", diverging = FALSE)
b2 <- brain(dk_lh, dk_scores, "ABCD PLS2 scores")
b3 <- brain(dk_lh, setNames(nspn34$PLS2, rownames(nspn34)), "NSPN PLS2 (Whitaker 2016)")
b4 <- brain(dk_lh, setNames(c3dk$C3, rownames(c3dk)), "AHBA C3 (Dear 2024)")
b5 <- brain(hcp_lh, hcp_dCT, "thinning rate (dCT)", diverging = FALSE)
b6 <- brain(hcp_lh, hcp_scores, "ABCD PLS2 scores")
b7 <- blank_note("no HCP-MMP map:\nNSPN PLS2 exists only\nin 308-region /DK space")
b8 <- brain(hcp_lh, hcp_c3, "AHBA C3 (native space)")

row_lab <- function(txt, sub) ggplot() +
  annotate("text", 0, 0.20, label = txt, size = 2.35, fontface = "bold", hjust = 0.5) +
  annotate("text", 0, -0.32, label = sub, size = 1.9, colour = "grey40", hjust = 0.5) +
  xlim(-1, 1) + ylim(-1, 1) + coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 1, 0, 1))

# ---------------------------------------------------------------- scatters ----
sm <- theme_bw(base_size = 6.6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.15, colour = "grey94"),
        plot.title = element_text(size = 6.6, face = "bold", margin = margin(b = 0.5)),
        plot.subtitle = element_text(size = 6.1, colour = "grey30", margin = margin(b = 1.5)),
        axis.title = element_text(size = 6.2), axis.text = element_text(size = 5.8),
        plot.margin = margin(1, 3, 1, 3))

scat <- function(x, y, xl, yl, title, sub, col = "grey25") {
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_point(size = 0.55, colour = col, alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.3, colour = "#b2182b", formula = y ~ x) +
    labs(x = xl, y = yl, title = title, subtitle = sub) + sm
}
hexp <- function(x, y, xl, yl, title, sub) {
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_hex(bins = 34, linewidth = 0) +
    scale_fill_gradient(low = "grey88", high = "grey15", guide = "none") +
    labs(x = xl, y = yl, title = title, subtitle = sub) + sm
}

dk_c3 <- cc(sprintf("C3_recomputed_%s", LEAD_DS)); dk_nspn <- cc("NSPN_PLS2")
hcp_c3s <- hc("C3", "scores"); hcp_c3w <- hc("C3", "weights")
dk_c3w <- wcc("C3"); dk_nspnw <- wcc("NSPN_PLS2_z"); hcp_dkw <- hc("DK_PLS2_ds25", "weights")

s1 <- scat(dkS$C3_recomputed, dkS$thinning_scores_gene_side, "AHBA C3 score", "PLS2 score",
           "e   DK scores vs C3", sprintf("rho=%.2f, %s", dk_c3$rho, pf(dk_c3$p)))
s2 <- scat(dkS$NSPN_PLS2, dkS$thinning_scores_gene_side, "NSPN PLS2 score", "PLS2 score",
           "f   DK scores vs NSPN", sprintf("rho=%.2f, %s", dk_nspn$rho, pf(dk_nspn$p)))
s3 <- scat(hcpS$C3, hcpS$thinning_score, "AHBA C3 score", "PLS2 score",
           "g   HCP scores vs C3", sprintf("rho=%.2f, %s", hcp_c3s$rho, pf(hcp_c3s$p)),
           col = "#2166ac")

# each matrix covers a different gene set, so drop missing weights before
# intersecting -- the resulting n must match the concordance tables (asserted).
dkz <- with(dkW[!is.na(dkW$thinning_Z_ds25), ], setNames(thinning_Z_ds25, gene))
hcpz <- with(hcpW[!is.na(hcpW$hcp_opt2_dCT_CT_PLS2_Z), ],
             setNames(-hcp_opt2_dCT_CT_PLS2_Z, gene))
c3v <- with(c3w[!is.na(c3w$C3), ], setNames(C3, c3w[[1]][!is.na(c3w$C3)]))
g1 <- intersect(names(dkz), names(c3v)); g2 <- intersect(names(hcpz), names(c3v))
g3 <- intersect(names(dkz), names(hcpz))
stopifnot(length(g1) == dk_c3w$n, length(g2) == hcp_c3w$n, length(g3) == hcp_dkw$n)
h1 <- hexp(c3v[g1], dkz[g1], "C3 gene weight", "DK PLS2 gene Z", "h   DK weights vs C3",
           sprintf("rho=%.2f, n=%.1fk genes", dk_c3w$rho, length(g1) / 1000))
h2 <- hexp(c3v[g2], hcpz[g2], "C3 gene weight", "HCP PLS2 gene Z", "i   HCP weights vs C3",
           sprintf("rho=%.2f, n=%.1fk genes", hcp_c3w$rho, length(g2) / 1000))
h3 <- hexp(dkz[g3], hcpz[g3], "DK PLS2 gene Z", "HCP PLS2 gene Z", "j   HCP vs DK weights",
           sprintf("rho=%.2f, n=%.1fk genes", hcp_dkw$rho, length(g3) / 1000))

# --------------------------------------------------------------- cell classes -
ord <- cts |> filter(vector == "AHBA C3") |> arrange(z) |> pull(cell_class)
CT <- cts |> mutate(cell_class = factor(cell_class, levels = ord),
                    vector = factor(vector, c("ABCD PLS2, DK", "ABCD PLS2, HCP-MMP", "AHBA C3")),
                    sig = ifelse(p_perm < 0.05, "p < 0.05", "n.s."))
astro <- CT |> filter(cell_class == "Astro")
pk <- ggplot(CT, aes(z, cell_class, colour = vector, shape = vector)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_point(aes(alpha = sig), position = position_dodge(width = 0.62), size = 1.25) +
  scale_colour_manual(values = c("ABCD PLS2, DK" = "#1f78b4", "ABCD PLS2, HCP-MMP" = "#33a02c",
                                 "AHBA C3" = "grey20")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_alpha_manual(values = c("p < 0.05" = 1, "n.s." = 0.3), guide = "none") +
  labs(x = "marker enrichment z (vs 20,000 size-matched random gene sets)", y = NULL,
       title = "k   Cell-class profile of the three gene rankings",
       subtitle = sprintf("Neuronal up, glial down in all three \u2014 except astrocytes, which flip sign in HCP-MMP (z = %+.1f) against DK (%+.1f) and C3 (%+.1f)",
                          astro$z[astro$vector == "ABCD PLS2, HCP-MMP"],
                          astro$z[astro$vector == "ABCD PLS2, DK"],
                          astro$z[astro$vector == "AHBA C3"])) +
  theme_bw(base_size = 7.2) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey94"),
        panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey94"),
        plot.title = element_text(size = 7.6, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.4, colour = "grey30", margin = margin(b = 2)),
        axis.title = element_text(size = 6.8), axis.text = element_text(size = 6.4),
        legend.position = "right", legend.title = element_blank(),
        legend.text = element_text(size = 6.3), legend.key.height = unit(9, "pt"),
        plot.margin = margin(2, 4, 2, 4))

# ------------------------------------------------------------------ assemble --
dkl <- comps |> filter(option == LEAD_OPT, ds == LEAD_DS, component == LEAD_COMP)
hl <- hcomps |> filter(option == "hcp_opt2_dCT_CT", component == "PLS2")
# each brain row carries its OWN plot_layout(widths=); a single widths= on the
# outer "/" of two rows would be applied to the 2-slot outer layout and error.
WID <- c(0.62, 1, 1, 1, 1)
r1 <- (row_lab("Desikan\u2013\nKilliany", "33 of 34\nparcels") | b1 | b2 | b3 | b4) +
  plot_layout(widths = WID)
r2 <- (row_lab("HCP-MMP\n(Glasser)", "137 of 180\nparcels") | b5 | b6 | b7 | b8) +
  plot_layout(widths = WID)
scatters <- (s1 | s2 | s3 | h1 | h2 | h3) + plot_layout(nrow = 1)

methods <- paste(
  "Methods.",
  "\u2022 Y = bilateral thinning rate (dCT) + baseline thickness (CT); PLS-SVD of gene expression against it. PLS1 takes the static thickness gradient, so the thinning signature is PLS2.",
  sprintf("\u2022 DK fit: AHBA_updated native-DK matrix at ds25 (33 regions, 12,007 genes), spin p = %.3f. HCP fit: hcp_3d.csv (137 parcels, 7,973 genes), spin p = %.3f \u2014 5,000 rotations of the complete map in each case.",
          dkl$p_spin_singular[1], hl$p_spin_singular[1]),
  sprintf("\u2022 Gene weights are bootstrap Z over 1,000 region resamples (reproducibility %.2f DK, %.2f HCP). Panels e\u2013j are Spearman; gene-level n differs because the two expression matrices differ.",
          dkl$boot_reproducibility[1], hl$boot_reproducibility[1]),
  "\u2022 Imaging: 8,192 subjects (DK) and 5,947 (HCP-MMP, derived from the FreeSurfer surfaces, coverage still incomplete), \u22652 visits, release 7.0 tables.",
  "\u2022 Panel k: Seidlitz et al. 2020 marker compilation; each vector against its own universe. Marker genes are co-expressed, so the independent-gene null is anti-conservative \u2014 read the pattern, not the absolute z.",
  "\u2022 Everything is in the thinning orientation: positive = expressed more where adolescent thinning is faster. Grey parcels have no AHBA donor coverage.",
  sep = "\n")

# "/" flattens nested stacks, so the outer layout sees FOUR rows (the two brain
# rows, the scatter strip, the cell-class panel) and needs four heights.
fig <- (r1 / r2 / scatters / pk) +
  plot_layout(heights = c(0.62, 0.62, 0.95, 1.05)) +
  plot_annotation(
    title = "One transcriptomic signature of adolescent thinning, derived twice from the ABCD",
    subtitle = "The same axis appears in both parcellations and matches both published signatures; the HCP-MMP version differs from the DK one mainly in its astrocyte loading.",
    caption = methods,
    theme = theme(plot.title = element_text(size = 9.6, face = "bold"),
                  plot.subtitle = element_text(size = 7.4, colour = "grey25", margin = margin(b = 4)),
                  plot.caption = element_text(size = 6.1, colour = "grey25", hjust = 0,
                                              lineheight = 1.42, margin = margin(t = 6))))

ggsave(file.path(FIG, "fig_signature_both.png"), fig, width = 8.2, height = 7.4, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_signature_both.png"), "\n")
