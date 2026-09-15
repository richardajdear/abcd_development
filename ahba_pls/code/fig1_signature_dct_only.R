#!/usr/bin/env Rscript
# fig1_signature_dct_only.R -- the single-Y variant of fig1_signature_both:
# Y = thinning rate (dCT) ALONE, in both parcellations.  Same panel structure
# as the main figure so the two can be read side by side.
#
# Reads only saved tables and fits nothing; every printed statistic comes from
# the table plotted.
#   results/dct_only_{scores_dk,scores_hcp,weights,components,concordance}  (16_dct_only.py)
#   results/hcp_y_maps_180.csv, celltype_all_options.tsv, hcp_run_provenance.tsv
#   data/{y_maps_bilateral_34,dk_polygons,hcp_polygons}.csv
#   data/reference/{nspn_dk_maps_bilateral_34,ahba_c123_scores_recomputed_ds25,
#                   ahba_c123_gene_weights}.csv
# Writes figures/fig_signature_dct_only.png

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})

ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); DATA <- file.path(ROOT, "data")
REF <- file.path(DATA, "reference"); FIG <- file.path(ROOT, "figures")
DS <- "ds25"

y34    <- read.csv(file.path(DATA, "y_maps_bilateral_34.csv"), row.names = 1)
dkS    <- read.csv(file.path(RES, "dct_only_scores_dk.csv"), row.names = 1)
hcpS   <- read.csv(file.path(RES, "dct_only_scores_hcp.csv"), row.names = 1)
W      <- read.delim(file.path(RES, "dct_only_weights.tsv"))
CO     <- read.delim(file.path(RES, "dct_only_components.tsv"))
CC     <- read.delim(file.path(RES, "dct_only_concordance.tsv"))
cts    <- read.delim(file.path(RES, "celltype_all_options.tsv"))
c3w    <- read.csv(file.path(REF, "ahba_c123_gene_weights.csv"))
nspn34 <- read.csv(file.path(REF, "nspn_dk_maps_bilateral_34.csv"), row.names = 1)
c3dk   <- read.csv(file.path(REF, sprintf("ahba_c123_scores_recomputed_%s.csv", DS)), row.names = 1)
dkpoly <- read.csv(file.path(DATA, "dk_polygons.csv"))
hcppoly <- read.csv(file.path(DATA, "hcp_polygons.csv"))
hcpY   <- read.csv(file.path(RES, "hcp_y_maps_180.csv"), row.names = 1)
PROV   <- read.delim(file.path(RES, "hcp_run_provenance.tsv"))

# every statistic is looked up from the concordance table by (parcellation, reference, level)
g <- function(parc, ref, lvl) {
  r <- CC |> filter(parcellation == parc, reference == ref, level == lvl)
  stopifnot(nrow(r) == 1)
  list(rho = r$rho[1], p = r$p_spin[1], n = r$n[1])
}
cmp <- function(parc) { r <- CO |> filter(parcellation == parc); stopifnot(nrow(r) == 1); r }
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

b1 <- brain(dk_lh, setNames(y34$dCT, rownames(y34)), "thinning rate (dCT)", diverging = FALSE)
b2 <- brain(dk_lh, setNames(dkS$thinning_score, rownames(dkS)), "dCT-only component scores")
b3 <- brain(dk_lh, setNames(nspn34$PLS2, rownames(nspn34)), "NSPN PLS2 (Whitaker 2016)")
b4 <- brain(dk_lh, setNames(c3dk$C3, rownames(c3dk)), "AHBA C3 (Dear 2024)")
b5 <- brain(hcp_lh, setNames(hcpY$dCT, rownames(hcpY)), "thinning rate (dCT)", diverging = FALSE)
b6 <- brain(hcp_lh, setNames(hcpS$thinning_score, rownames(hcpS)), "dCT-only component scores")
b7 <- blank_note("no HCP-MMP map:\nNSPN PLS2 exists only\nin 308-region /DK space")
b8 <- brain(hcp_lh, setNames(hcpS$C3, rownames(hcpS)), "AHBA C3 (native space)")

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

dk_c3 <- g("DK", "C3", "scores"); dk_nspn <- g("DK", "NSPN_PLS2", "scores")
hcp_c3s <- g("HCP", "C3", "scores")
dk_c3w <- g("DK", "C3", "weights"); dk_c1w <- g("DK", "C1", "weights")
hcp_c3w <- g("HCP", "C3", "weights"); hcp_c1w <- g("HCP", "C1", "weights")
hcp_dkw <- g("HCP", "DK_dCT_PLS1", "weights")
dk_leadw <- g("DK", "DK_opt2_PLS2", "weights"); hcp_leadw <- g("HCP", "HCP_opt2_PLS2", "weights")

s1 <- scat(dkS$C3_recomputed, dkS$thinning_score, "AHBA C3 score", "dCT-only score",
           "e   DK scores vs C3", sprintf("rho=%.2f, %s", dk_c3$rho, pf(dk_c3$p)))
s2 <- scat(dkS$NSPN_PLS2, dkS$thinning_score, "NSPN PLS2 score", "dCT-only score",
           "f   DK scores vs NSPN", sprintf("rho=%.2f, %s", dk_nspn$rho, pf(dk_nspn$p)))
s3 <- scat(hcpS$C3, hcpS$thinning_score, "AHBA C3 score", "dCT-only score",
           "g   HCP scores vs C3", sprintf("rho=%.2f, %s", hcp_c3s$rho, pf(hcp_c3s$p)),
           col = "#2166ac")

dkz <- with(W[!is.na(W$dk_thinning_Z), ], setNames(dk_thinning_Z, gene))
hcpz <- with(W[!is.na(W$hcp_thinning_Z), ], setNames(hcp_thinning_Z, gene))
c3v <- with(c3w[!is.na(c3w$C3), ], setNames(C3, c3w[[1]][!is.na(c3w$C3)]))
g1 <- intersect(names(dkz), names(c3v)); g2 <- intersect(names(hcpz), names(c3v))
g3 <- intersect(names(dkz), names(hcpz))
stopifnot(length(g1) == dk_c3w$n, length(g2) == hcp_c3w$n, length(g3) == hcp_dkw$n)

h1 <- hexp(c3v[g1], dkz[g1], "C3 gene weight", "DK dCT-only gene Z", "h   DK weights vs C3",
           sprintf("vs C3 %.2f; vs C1 %.2f", dk_c3w$rho, dk_c1w$rho))
h2 <- hexp(c3v[g2], hcpz[g2], "C3 gene weight", "HCP dCT-only gene Z", "i   HCP weights vs C3",
           sprintf("vs C3 %.2f; vs C1 %.2f", hcp_c3w$rho, hcp_c1w$rho))
h3 <- hexp(dkz[g3], hcpz[g3], "DK dCT-only gene Z", "HCP dCT-only gene Z",
           "j   HCP vs DK weights",
           sprintf("rho=%.2f, %.1fk genes", hcp_dkw$rho, length(g3) / 1000))

# --------------------------------------------------------------- cell classes -
VEC <- c("ABCD dCT alone, DK", "ABCD dCT alone, HCP-MMP", "AHBA C3")
cs <- cts |> filter(universe == "shared", vector %in% VEC)
stopifnot(length(unique(cs$n_universe)) == 1, setequal(unique(cs$vector), VEC))
NU <- cs$n_universe[1]
ord <- cs |> filter(vector == "AHBA C3") |> arrange(z) |> pull(cell_class)
CTd <- cs |> mutate(cell_class = factor(cell_class, levels = ord),
                    vector = factor(vector, VEC),
                    sig = ifelse(p_perm < 0.05, "p < 0.05", "n.s."))
astro <- CTd |> filter(cell_class == "Astro")
pk <- ggplot(CTd, aes(z, cell_class, colour = vector, shape = vector)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_point(aes(alpha = sig), position = position_dodge(width = 0.62), size = 1.25) +
  scale_colour_manual(values = setNames(c("#1f78b4", "#33a02c", "grey20"), VEC)) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_alpha_manual(values = c("p < 0.05" = 1, "n.s." = 0.3), guide = "none") +
  labs(x = sprintf("marker enrichment z (vs 20,000 random gene sets, shared universe of %s genes)",
                   format(NU, big.mark = ",")), y = NULL,
       title = "k   Cell-class profile of the dCT-only rankings",
       subtitle = sprintf("Neuronal up, glial down in all three; astrocytes again split by atlas (HCP-MMP z = %+.1f, DK %+.1f, C3 %+.1f)",
                          astro$z[astro$vector == "ABCD dCT alone, HCP-MMP"],
                          astro$z[astro$vector == "ABCD dCT alone, DK"],
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
cd <- cmp("DK"); ch <- cmp("HCP")
WID <- c(0.62, 1, 1, 1, 1)
r1 <- (row_lab("Desikan\u2013\nKilliany", sprintf("%d of 34\nparcels", cd$n_regions)) | b1 | b2 | b3 | b4) +
  plot_layout(widths = WID)
r2 <- (row_lab("HCP-MMP\n(Glasser)", sprintf("%d of 180\nparcels", ch$n_regions)) | b5 | b6 | b7 | b8) +
  plot_layout(widths = WID)
scatters <- (s1 | s2 | s3 | h1 | h2 | h3) + plot_layout(nrow = 1)

methods <- paste(
  "Methods.",
  "\u2022 Y = bilateral thinning rate (dCT) ALONE. With one Y column the PLS has a single component, so PLS1 is the vector of gene\u2013map correlations \u2014 no static",
  "  component is estimated, and nothing absorbs the baseline-thickness gradient (compare fig_signature_both.png, where Y = dCT + CT and the signature is PLS2).",
  sprintf("\u2022 DK fit: AHBA_updated native-DK matrix at %s (%d regions, %s genes), spin p = %.3f, bootstrap reproducibility %.2f. HCP fit: hcp_3d.csv (%d parcels, %s genes), spin p = %.3f, reproducibility %.2f.",
          DS, cd$n_regions, format(cd$n_genes, big.mark = ","), cd$p_spin_singular, cd$boot_reproducibility,
          ch$n_regions, format(ch$n_genes, big.mark = ","), ch$p_spin_singular, ch$boot_reproducibility),
  "  Spin p is the component's own singular value against 5,000 rotations of the complete map; NEITHER single-Y component reaches 0.05, so the map-level covariance is not spatially specific.",
  sprintf("\u2022 Agreement with the two-Y signature of the main figure: gene weights rho = %.2f in DK and %.2f in HCP-MMP.",
          dk_leadw$rho, hcp_leadw$rho),
  sprintf("\u2022 Imaging: %s children with \u22652 QC-passing visits, true 7.0 tabulated tables; both parcellations run on the same sample.",
          format(PROV$n_subjects[1], big.mark = ",")),
  "\u2022 Panel k: Seidlitz et al. 2020 markers on one shared universe. Marker genes are co-expressed, so the independent-gene null is anti-conservative \u2014 read the pattern, not the absolute z.",
  "\u2022 Everything is in the thinning orientation: positive = expressed more where adolescent thinning is faster. Grey parcels have no AHBA donor coverage.",
  sep = "\n")

fig <- (r1 / r2 / scatters / pk) +
  plot_layout(heights = c(0.62, 0.62, 0.95, 1.05)) +
  plot_annotation(
    title = "Thinning rate alone: the same axis at 137 parcels, a mixture at 33",
    subtitle = sprintf("With no static map in Y, the DK component confounds the axis with the static gradient (weights rho %.2f vs C3, %.2f vs C1; scores %s).\nAt 137 parcels the confound largely disappears (%.2f vs C3, %.2f vs C1) and the vector is rho %.2f with the two-Y signature \u2014 but neither single-Y component is spin-significant.",
                       dk_c3w$rho, dk_c1w$rho, pf(dk_c3$p), hcp_c3w$rho, hcp_c1w$rho, hcp_leadw$rho),
    caption = methods,
    theme = theme(plot.title = element_text(size = 9.6, face = "bold"),
                  plot.subtitle = element_text(size = 7.4, colour = "grey25", margin = margin(b = 4)),
                  plot.caption = element_text(size = 6.1, colour = "grey25", hjust = 0,
                                              lineheight = 1.42, margin = margin(t = 6))))

ggsave(file.path(FIG, "fig_signature_dct_only.png"), fig, width = 8.2, height = 7.4, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_signature_dct_only.png"), "\n")
