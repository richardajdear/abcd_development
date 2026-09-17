#!/usr/bin/env Rscript
# fig1_signature_designs.R -- the thinning signature in both parcellations, with
# EVERY pairwise comparison shown as two square pair matrices: regional maps and
# gene vectors.  Each matrix carries Desikan-Killiany BELOW the diagonal and
# HCP-MMP ABOVE it, which works because NSPN PLS2 now exists in HCP-MMP space
# too (18_nspn_to_hcp.py).
#
# Reads only saved tables and fits nothing; every statistic printed in a panel
# comes from the pairs table for that panel.
#   results/design_grid_points_scores.csv    regional values, long
#   results/design_grid_points_weights.tsv   gene vectors, wide
#   results/design_grid_pairs_scores.tsv     rho + spin p for every map pair
#   results/design_grid_pairs_weights.tsv    rho for every gene-vector pair
#   results/design_grid_components.tsv, hcp_y_maps_180.csv, hcp_run_provenance.tsv
#   data/{y_maps_bilateral_34,dk_polygons,hcp_polygons}.csv
#   data/reference/{nspn_dk_maps_bilateral_34,nspn_hcp_maps,
#                   ahba_c123_scores_recomputed_ds25}.csv
# Writes figures/fig_signature_designs.png
#
# The four-feature design (CT + dCT + T1w/T2w + dT1w/T2w) is computed by
# 17_design_grid.py but not plotted: it exists in DK only and its component is
# not spin-significant.  Its numbers are in design_grid_concordance.tsv.

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
                  library(scales); library(grid)})

ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); DATA <- file.path(ROOT, "data")
REF <- file.path(DATA, "reference"); FIG <- file.path(ROOT, "figures")

PS   <- read.delim(file.path(RES, "design_grid_pairs_scores.tsv"))
PW   <- read.delim(file.path(RES, "design_grid_pairs_weights.tsv"))
SP   <- read.csv(file.path(RES, "design_grid_points_scores.csv"), check.names = FALSE)
WP   <- read.delim(file.path(RES, "design_grid_points_weights.tsv"), check.names = FALSE)
CMP  <- read.delim(file.path(RES, "design_grid_components.tsv"))
PROV <- read.delim(file.path(RES, "hcp_run_provenance.tsv"))
y34  <- read.csv(file.path(DATA, "y_maps_bilateral_34.csv"), row.names = 1)
hcpY <- read.csv(file.path(RES, "hcp_y_maps_180.csv"), row.names = 1)
# the cached DK polygon table holds both hemispheres and four views; the HCP one
# is already left-lateral/medial. Match them, or the DK row renders four tiny
# views per map instead of two.
dkp  <- read.csv(file.path(DATA, "dk_polygons.csv")) |>
  filter(hemi == "left", view %in% c("lateral", "medial"))
hcpp <- read.csv(file.path(DATA, "hcp_polygons.csv")) |>
  filter(view %in% c("lateral", "medial"))
stopifnot(length(unique(dkp$view)) == 2, length(unique(hcpp$view)) == 2)
nspn34  <- read.csv(file.path(REF, "nspn_dk_maps_bilateral_34.csv"), row.names = 1)
nspnH   <- read.csv(file.path(REF, "nspn_hcp_maps.csv"), row.names = 1)
c3dk <- read.csv(file.path(REF, "ahba_c123_scores_recomputed_ds25.csv"), row.names = 1)

S_VARS <- c("dCT rate", "dCT + CT", "dCT alone", "NSPN PLS2", "AHBA C3")
# AHBA C1 and the second snRNAseq dataset (PC1_U01V2) are computed in
# design_grid_pairs_weights.tsv but are not columns here: C1 is the
# static-gradient control rather than an axis under comparison, and the U01
# version is a replication of the plotted Herring one (they agree at rho 0.54).
W_VARS <- c("dCT + CT", "dCT alone", "NSPN PLS2", "AHBA C3", "snRNAseq PC1")
ODD <- c(a = "dCT rate", b = "snRNAseq PC1")   # shaded rows/columns, see below
PAL <- c(DK = "#d6604d", HCP = "#2166ac")   # light red / blue, one per parcellation
stopifnot(setequal(unique(c(PS$var_x, PS$var_y)), S_VARS),
          all(W_VARS %in% c(PW$var_x, PW$var_y)))

pf <- function(p) if (is.na(p)) "" else if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)

# ------------------------------------------------------------------- brains ----
brain <- function(poly, vals, title = NULL, diverging = TRUE) {
  d <- poly |> mutate(v = unname(vals[label]))
  lim <- if (diverging) max(abs(d$v), na.rm = TRUE) * c(-1, 1) else range(d$v, na.rm = TRUE)
  ggplot(d, aes(x, y, group = interaction(view, label, group, subgroup), fill = v)) +
    geom_polygon(colour = "grey35", linewidth = 0.07) +
    coord_fixed(expand = FALSE) +
    (if (diverging) scale_fill_distiller(palette = "RdBu", limits = lim, na.value = "grey78")
     else scale_fill_viridis_c(na.value = "grey78")) +
    # No facet by view: in the ggseg polygon space the lateral and medial views
    # already occupy adjacent, non-overlapping x ranges (DK 747-1379 and
    # 1411-2043), so one coord_fixed panel draws them side by side at full size.
    # Faceting instead gives each view the UNION x range, which draws every
    # surface at 49% of its panel width -- and coord_fixed then halves the height
    # too (free scales are not an option: coord_fixed rejects them).
    theme_void() +
    theme(legend.position = "none", strip.text = element_blank(),
          plot.title = element_text(size = 7.0, hjust = 0.5, margin = margin(b = 1)),
          plot.margin = margin(0, 1, 0, 1)) +
    labs(title = title)
}
row_lab <- function(txt, sub) ggplot() +
  annotate("text", 0, 0.45, label = txt, size = 2.1, fontface = "bold", hjust = 0.5) +
  annotate("text", 0, -0.48, label = sub, size = 1.7, colour = "grey40", hjust = 0.5) +
  xlim(-1, 1) + ylim(-1, 1) + coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 1, 0, 1))

sc <- function(parc, v) {
  x <- SP |> filter(parcellation == parc, variable == v)
  setNames(x$value, x$label)
}
dk_maps <- list(
  brain(dkp, setNames(y34$dCT, rownames(y34)), "dCT rate (mm/yr)", diverging = FALSE),
  brain(dkp, sc("DK", "dCT + CT"), "dCT + CT"), brain(dkp, sc("DK", "dCT alone"), "dCT alone"),
  brain(dkp, setNames(nspn34$PLS2, rownames(nspn34)), "NSPN PLS2"),
  brain(dkp, setNames(c3dk$C3, rownames(c3dk)), "AHBA C3"))
hcp_maps <- list(
  brain(hcpp, setNames(hcpY$dCT, rownames(hcpY)), diverging = FALSE),
  brain(hcpp, sc("HCP", "dCT + CT")), brain(hcpp, sc("HCP", "dCT alone")),
  brain(hcpp, setNames(nspnH$PLS2, rownames(nspnH))), brain(hcpp, sc("HCP", "AHBA C3")))

WID <- c(0.62, rep(1, length(dk_maps)))
brow <- function(lab, sub, panels)
  Reduce(`|`, panels, init = row_lab(lab, sub)) + plot_layout(widths = WID)
r1 <- brow("Desikan\u2013\nKilliany", "33 of 34\nparcels", dk_maps)
r2 <- brow("HCP-MMP\n(Glasser)", "137 of 180\nparcels", hcp_maps)

# -------------------------------------------------------------- pair matrix ----
# One cell per ordered (row, column) pair of variables: below the diagonal the
# DK version of the pair, above it the HCP-MMP version, on the diagonal the
# variable name. Both triangles carry real, different data for regional maps;
# for gene vectors the reference-vs-reference cells (C3/C1/NSPN) are
# parcellation-independent and so are identical in the two triangles.
shade <- function(vars_, odd) {
  # every cell in the odd variable's row or column, including its diagonal cell,
  # so the band is continuous
  expand.grid(row = vars_, col = vars_, stringsAsFactors = FALSE) |>
    filter(row == odd | col == odd) |>
    mutate(row = factor(row, levels = vars_), col = factor(col, levels = vars_))
}

cells <- function(vars_) {
  expand.grid(row = vars_, col = vars_, stringsAsFactors = FALSE) |>
    mutate(i = match(row, vars_), j = match(col, vars_),
           parcellation = ifelse(i > j, "DK", ifelse(i < j, "HCP", NA)),
           row = factor(row, levels = vars_), col = factor(col, levels = vars_))
}

stat_of <- function(P, r, c, parc) {
  h <- P |> filter(parcellation == parc,
                   (var_x == r & var_y == c) | (var_x == c & var_y == r))
  stopifnot(nrow(h) == 1)
  h
}

# --- regional maps
SC <- cells(S_VARS) |> filter(!is.na(parcellation))
spts <- do.call(rbind, lapply(seq_len(nrow(SC)), function(k) {
  z <- SC[k, ]; x <- sc(z$parcellation, as.character(z$col)); y <- sc(z$parcellation, as.character(z$row))
  lab <- intersect(names(x), names(y))
  data.frame(row = z$row, col = z$col, parcellation = z$parcellation, x = x[lab], y = y[lab])
}))
slab <- SC |> rowwise() |>
  mutate(st = list(stat_of(PS, as.character(row), as.character(col), parcellation))) |>
  mutate(txt = sprintf("%.2f", st$rho), sub = pf(st$p_spin)) |> ungroup()

mat_theme <- theme_bw(base_size = 6.4) +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill = "grey96", colour = NA),
        strip.text = element_text(size = 5.9, margin = margin(1.2, 1.2, 1.2, 1.2)),
        strip.text.y.right = element_text(angle = 90),
        panel.spacing = unit(1.2, "pt"), plot.margin = margin(1, 1, 1, 1),
        plot.title = element_text(size = 7.6, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.2, colour = "grey30", margin = margin(b = 2)))

pa <- ggplot(spts, aes(x, y)) +
  geom_rect(data = shade(S_VARS, ODD[["a"]]), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey93", inherit.aes = FALSE) +
  geom_point(aes(colour = parcellation), size = 0.28, alpha = 0.75) +
  geom_smooth(aes(colour = parcellation), method = "lm", se = FALSE,
              linewidth = 0.25, formula = y ~ x) +
  geom_text(data = slab, aes(x = -Inf, y = Inf, label = txt), hjust = -0.25, vjust = 1.35,
            size = 2.05, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = slab, aes(x = -Inf, y = Inf, label = sub), hjust = -0.3, vjust = 3.1,
            size = 1.6, colour = "grey35", inherit.aes = FALSE) +
  scale_colour_manual(values = PAL, guide = "none") +
  facet_grid(row ~ col, scales = "free", switch = NULL) + mat_theme +
  labs(title = "a   Regional maps \u2014 every pair",
       subtitle = "lower triangle DK (33 regions, red), upper HCP-MMP (137, blue)  \u00b7  bold rho, spin p below")

# --- gene vectors
WC <- cells(W_VARS) |> filter(!is.na(parcellation))
wcol <- function(v, parc) if (v %in% c("dCT + CT", "dCT alone")) paste0(v, "|", parc) else v
wpts <- do.call(rbind, lapply(seq_len(nrow(WC)), function(k) {
  z <- WC[k, ]
  x <- WP[[wcol(as.character(z$col), z$parcellation)]]
  y <- WP[[wcol(as.character(z$row), z$parcellation)]]
  ok <- is.finite(x) & is.finite(y)
  data.frame(row = z$row, col = z$col, parcellation = z$parcellation, x = x[ok], y = y[ok])
}))
wlab <- WC |> rowwise() |>
  mutate(st = list(stat_of(PW, as.character(row), as.character(col), parcellation))) |>
  mutate(txt = sprintf("%.2f", st$rho), sub = sprintf("n=%.1fk", st$n / 1000)) |> ungroup()
# the three reference-vs-reference cells are the same data in both triangles
REFS <- c("NSPN PLS2", "AHBA C3", "snRNAseq PC1")
ref_ref <- WC |> filter(row %in% REFS, col %in% REFS, parcellation == "HCP") |>
  mutate(mark = "=")

pb <- ggplot(wpts, aes(x, y)) +
  geom_rect(data = shade(W_VARS, ODD[["b"]]), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey93", inherit.aes = FALSE) +
  geom_hex(aes(fill = parcellation, alpha = after_stat(count)), bins = 22, linewidth = 0) +
  geom_text(data = wlab, aes(x = -Inf, y = Inf, label = txt), hjust = -0.25, vjust = 1.35,
            size = 2.05, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = wlab, aes(x = -Inf, y = Inf, label = sub), hjust = -0.3, vjust = 3.1,
            size = 1.6, colour = "grey35", inherit.aes = FALSE) +
  geom_text(data = ref_ref, aes(x = Inf, y = Inf, label = mark), hjust = 1.4, vjust = 1.4,
            size = 2.2, colour = "grey55", inherit.aes = FALSE) +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_alpha_continuous(range = c(0.25, 1), guide = "none") +
  facet_grid(row ~ col, scales = "free") + mat_theme +
  labs(title = "b   Gene weights \u2014 every pair",
       subtitle = "same triangles  \u00b7  bold rho, shared-gene n below  \u00b7  \u201c=\u201d: parcellation-independent")

# ------------------------------------------------------------------ caption ---
d1 <- CMP |> filter(design == "dCT + CT", parcellation == "DK")
h1 <- CMP |> filter(design == "dCT + CT", parcellation == "HCP")
d2 <- CMP |> filter(design == "dCT alone", parcellation == "DK")
h2 <- CMP |> filter(design == "dCT alone", parcellation == "HCP")
d4 <- CMP |> filter(grepl("T1T2", design), parcellation == "DK")
nspn_hcp_n <- sum(!is.na(nspnH$PLS2))
pc1_pair <- PW |> filter(parcellation == "DK", grepl("snRNAseq", var_x), grepl("snRNAseq", var_y))
stopifnot(nrow(pc1_pair) == 1)
pc1_agree <- pc1_pair$rho[1]
pc1_n <- max(PW$n[grepl("snRNAseq", PW$var_x) | grepl("snRNAseq", PW$var_y)])
methods <- paste(
  "Methods.",
  "\u2022 Y designs: dCT + CT (PLS2, the lead signature) and dCT alone (a single Y column, so PLS1 is simply the vector of gene\u2013map correlations and nothing absorbs the",
  "  baseline-thickness gradient). A four-feature design (CT + dCT + T1w/T2w + dT1w/T2w) is computed by 17_design_grid.py but not shown: DK-only, and its C3-aligned",
  sprintf("  component is not spin-significant (p = %.3f). Its numbers are in design_grid_concordance.tsv.", d4$p_spin),
  sprintf("\u2022 Spin p is each component's own singular value against 5,000 rotations: %.3f (dCT+CT, DK), %.3f (dCT+CT, HCP), %.3f (dCT alone, DK), %.3f (dCT alone, HCP). Only the",
          d1$p_spin, h1$p_spin, d2$p_spin, h2$p_spin),
  "  dCT+CT fits clear 0.05, so the single-Y fits are usable gene rankings but carry no independent spatial claim.",
  sprintf("\u2022 NSPN PLS2 in HCP-MMP space (%d of 180 parcels) is resampled from the published 308-region map through fsaverage vertices, the route in AHBA/notebooks/", nspn_hcp_n),
  "  MT_whitakervertes.ipynb (code/18_nspn_to_hcp.py; validation against the published DK table r = 0.985, see data/reference/NSPN_HCP.md). The resampling is not",
  "  lossy in the way the weak HCP-NSPN cells might suggest: the target grid is FINER than the source (180 vs 152 parcels per hemisphere), each HCP parcel averages a",
  "  median of 4 source parcels, and an HCP-native map survives a round trip through the 308 grid at r = 0.91 (code/19_resample_check.py). What differs is the ABCD",
  "  component itself: its DK and HCP score maps correlate only rho = 0.74 across parcellations, and NSPN PLS2 tracks the DK one.",
  sprintf("\u2022 snRNAseq PC1 = the snRNA-seq maturation axis from the transcriptional_maturation project (data/velmeshev_PC1_gene_loadings.csv, %s genes). Two independent",
          format(pc1_n, big.mark = ",")),
  sprintf("  datasets exist and agree at rho = %.2f; the Herring one is plotted and the U01 one is in design_grid_pairs_weights.tsv. It is the only axis here not derived from the AHBA.",
          pc1_agree),
  sprintf("\u2022 Gene weights are bootstrap Z over 1,000 resamples; n = %s genes (DK) and %s (HCP-MMP), and each panel's n is the overlap of its two vectors.",
          format(d1$n_genes, big.mark = ","), format(h1$n_genes, big.mark = ",")),
  sprintf("\u2022 Imaging: %s children with \u22652 QC-passing visits, true 7.0 tabulated tables; both parcellations run on the same sample.",
          format(PROV$n_subjects[1], big.mark = ",")),
  sprintf("\u2022 The shaded row and column of each matrix is the odd one out: %s in a (an imaging map, not a gene-derived axis) and %s in b (the only axis not derived from the AHBA).",
          ODD[["a"]], ODD[["b"]]),
  "\u2022 Components are in the thinning orientation (multiplied by \u2212sign of their dCT salience): positive = expressed more where thinning is faster. dCT rate itself is",
  "  mm/yr, so it runs the other way \u2014 hence the negative correlations in its row. Grey parcels have no AHBA donor coverage.",
  sep = "\n")
methods_panel <- ggplot() +
  annotate("text", x = 0, y = 1, label = methods, hjust = 0, vjust = 1,
           size = 2.05, lineheight = 1.42, colour = "grey25") +
  xlim(0, 1) + ylim(0, 1) + coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(4, 2, 0, 2))

fig <- (r1 / r2 / (pa | pb) / methods_panel) +
  # aspect.ratio = 1 pins the matrices' panel shape, so these weights are set to
  # the physical height each row actually needs (inches): brains ~0.7, the 5x5
  # matrices ~4.35 at this width, methods ~1.3.
  # a brain cell is ~1.65 in wide and the two views together span 1295 x 425
  # polygon units, so each brain row needs ~0.55 in of height.
  plot_layout(heights = c(0.58, 0.58, 4.35, 2.9)) +
  plot_annotation(
    title = "One transcriptomic axis of adolescent thinning, and every pairwise comparison behind it",
    subtitle = sprintf("Adding baseline CT to Y recovers AHBA C3 in both parcellations (DK rho %.2f, HCP %.2f); thinning rate alone recovers it only at 137 parcels (%.2f vs %.2f),\nbecause at 33 regions its gene weights load on the static gradient C1 as heavily as on C3 (%.2f vs %.2f). NSPN PLS2 matches the DK fit (%.2f) but not the\nHCP one (%.2f) \u2014 not a resampling artefact (round trip r = 0.91) but the two ABCD fits' score maps differing (rho 0.74 across parcellations).",
                       stat_of(PS, "dCT + CT", "AHBA C3", "DK")$rho, stat_of(PS, "dCT + CT", "AHBA C3", "HCP")$rho,
                       stat_of(PS, "dCT alone", "AHBA C3", "HCP")$rho, stat_of(PS, "dCT alone", "AHBA C3", "DK")$rho,
                       stat_of(PW, "dCT alone", "AHBA C1", "DK")$rho, stat_of(PW, "dCT alone", "AHBA C3", "DK")$rho,
                       stat_of(PS, "dCT + CT", "NSPN PLS2", "DK")$rho, stat_of(PS, "dCT + CT", "NSPN PLS2", "HCP")$rho),
    theme = theme(plot.title = element_text(size = 9.6, face = "bold"),
                  plot.subtitle = element_text(size = 7.3, colour = "grey25",
                                               margin = margin(b = 4))))

ggsave(file.path(FIG, "fig_signature_designs.png"), fig, width = 8.6, height = 9.5, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_signature_designs.png"), "\n")
