#!/usr/bin/env Rscript
# fig1_signature_designs.R -- the signature under THREE Y-matrix designs and in
# both parcellations, with all the concordance scatters in one figure.
#
# Replaces the cell-class panel of fig_signature_both.png (that comparison now
# has its own figure, fig_celltypes.png) with the dCT-only and four-feature
# scatter sets, so the three designs can be read against each other directly.
#
# Reads only saved tables and fits nothing; every printed statistic comes from
# the table plotted.
#   results/design_grid_{scores,weights,components,concordance}   (17_design_grid.py)
#   results/hcp_y_maps_180.csv, hcp_run_provenance.tsv
#   data/{y_maps_bilateral_34,dk_polygons,hcp_polygons}.csv
#   data/reference/{nspn_dk_maps_bilateral_34,ahba_c123_scores_recomputed_ds25,
#                   ahba_c123_gene_weights}.csv
# Writes figures/fig_signature_designs.png

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})

ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); DATA <- file.path(ROOT, "data")
REF <- file.path(DATA, "reference"); FIG <- file.path(ROOT, "figures")
DS <- "ds25"

GS   <- read.csv(file.path(RES, "design_grid_scores.csv"))
GW   <- read.delim(file.path(RES, "design_grid_weights.tsv"))
GC   <- read.delim(file.path(RES, "design_grid_components.tsv"))
CC   <- read.delim(file.path(RES, "design_grid_concordance.tsv"))
PROV <- read.delim(file.path(RES, "hcp_run_provenance.tsv"))
y34  <- read.csv(file.path(DATA, "y_maps_bilateral_34.csv"), row.names = 1)
hcpY <- read.csv(file.path(RES, "hcp_y_maps_180.csv"), row.names = 1)
nspn34 <- read.csv(file.path(REF, "nspn_dk_maps_bilateral_34.csv"), row.names = 1)
c3dk <- read.csv(file.path(REF, sprintf("ahba_c123_scores_recomputed_%s.csv", DS)), row.names = 1)
c3w  <- read.csv(file.path(REF, "ahba_c123_gene_weights.csv"))
hcpS1 <- read.csv(file.path(RES, "dct_only_scores_hcp.csv"), row.names = 1)   # carries C1-C3 in HCP space
dkpoly <- read.csv(file.path(DATA, "dk_polygons.csv"))
hcppoly <- read.csv(file.path(DATA, "hcp_polygons.csv"))

# SHOW_FOUR toggles the four-feature (CT + dCT + T1w/T2w + dT1w/T2w) row. It is
# off by default: that design's component is not spin-significant (p = 0.098) and
# it exists in DK only, so the figure reads more cleanly without it. The row is
# still computed by 17_design_grid.py and its numbers are in the README table and
# in design_grid_*.tsv -- set this to TRUE to put it back on the figure.
SHOW_FOUR <- as.logical(Sys.getenv("SHOW_FOUR", "FALSE"))
ALL_DESIGNS <- c("dCT + CT", "dCT alone", "CT + dCT + T1T2 + dT1T2")
stopifnot(setequal(unique(GS$design), ALL_DESIGNS))
DESIGNS <- if (SHOW_FOUR) ALL_DESIGNS else ALL_DESIGNS[1:2]
have <- GC |> transmute(key = paste(design, parcellation)) |> pull(key)

# statistic lookups -- nothing on this figure is typed in
g <- function(des, parc, lvl, ref) {
  r <- CC |> filter(design == des, parcellation == parc, level == lvl, reference == ref)
  stopifnot(nrow(r) == 1); list(rho = r$rho[1], p = r$p_spin[1], n = r$n[1])
}
cmp <- function(des, parc) { r <- GC |> filter(design == des, parcellation == parc)
  stopifnot(nrow(r) == 1); r }
pf <- function(p) if (is.na(p)) "" else if (p < 1e-3) "p_spin<0.001" else sprintf("p_spin=%.3f", p)
scores_of <- function(des, parc) { d <- GS |> filter(design == des, parcellation == parc)
  setNames(d$score, d$label) }

# ------------------------------------------------------------------ brains ----
dk_lh <- dkpoly |> filter(hemi == "left", view %in% c("lateral", "medial")) |>
  mutate(grp = interaction(label, view, group, subgroup, drop = TRUE),
         region = sub("^lh_", "", label))
hcp_lh <- hcppoly |> mutate(grp = interaction(label, view, group, subgroup, drop = TRUE))

brain <- function(poly, values, title = NULL, diverging = TRUE) {
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
          plot.title = element_text(size = 6.3, hjust = 0.5, margin = margin(b = 0.5)),
          plot.margin = margin(0.5, 1, 0.5, 1))
}
note <- function(txt, size = 1.95) ggplot() +
  annotate("text", 0, 0, label = txt, size = size, colour = "grey45", lineheight = 1.15) +
  theme_void() + theme(plot.margin = margin(0.5, 1, 0.5, 1))
row_lab <- function(txt, sub) ggplot() +
  annotate("text", 0, 0.34, label = txt, size = 2.25, fontface = "bold", hjust = 0.5) +
  annotate("text", 0, -0.40, label = sub, size = 1.8, colour = "grey40", hjust = 0.5) +
  xlim(-1, 1) + ylim(-1, 1) + coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 1, 0, 1))

# columns: thinning rate, one per design shown, then the two published maps
des_hdr <- c("dCT + CT" = "dCT + CT", "dCT alone" = "dCT alone",
             "CT + dCT + T1T2 + dT1T2" = "four features")[DESIGNS]
dk_brains <- c(
  list(brain(dk_lh, setNames(y34$dCT, rownames(y34)), "thinning rate (dCT)", diverging = FALSE)),
  lapply(DESIGNS, function(d) brain(dk_lh, scores_of(d, "DK"), des_hdr[[d]])),
  list(brain(dk_lh, setNames(nspn34$PLS2, rownames(nspn34)), "NSPN PLS2"),
       brain(dk_lh, setNames(c3dk$C3, rownames(c3dk)), "AHBA C3")))
hcp_brains <- c(
  list(brain(hcp_lh, setNames(hcpY$dCT, rownames(hcpY)), diverging = FALSE)),
  lapply(DESIGNS, function(d)
    if (paste(d, "HCP") %in% have) brain(hcp_lh, scores_of(d, "HCP"))
    else note("no T1w/T2w in\nHCP-MMP\n(see caption)")),
  list(note("NSPN PLS2 exists\nonly in DK /\n308-region space"),
       brain(hcp_lh, setNames(hcpS1$C3, rownames(hcpS1)))))

WID <- c(0.55, rep(1, length(dk_brains)))
brow <- function(lab, sub, panels)
  Reduce(`|`, panels, init = row_lab(lab, sub)) + plot_layout(widths = WID)
r1 <- brow("Desikan\u2013\nKilliany", "33 of 34\nparcels", dk_brains)
r2 <- brow("HCP-MMP\n(Glasser)", "137 of 180\nparcels", hcp_brains)

# ---------------------------------------------------------------- scatters ----
sm <- theme_bw(base_size = 6.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.15, colour = "grey94"),
        plot.title = element_text(size = 6.4, face = "bold", margin = margin(b = 0.5)),
        plot.subtitle = element_text(size = 6.0, colour = "grey30", margin = margin(b = 1.5)),
        axis.title = element_text(size = 6.0), axis.text = element_text(size = 5.6),
        plot.margin = margin(1, 3, 1, 3))
PAL <- c(DK = "grey25", HCP = "#2166ac")

scat <- function(x, y, xl, yl, title, sub, parc) {
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_point(size = 0.5, colour = PAL[[parc]], alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.3, colour = "#b2182b", formula = y ~ x) +
    labs(x = xl, y = yl, title = title, subtitle = sub) + sm
}
hexp <- function(x, y, xl, yl, title, sub) {
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_hex(bins = 32, linewidth = 0) +
    scale_fill_gradient(low = "grey88", high = "grey15", guide = "none") +
    labs(x = xl, y = yl, title = title, subtitle = sub) + sm
}

c3vec <- with(c3w[!is.na(c3w$C3), ], setNames(C3, c3w[[1]][!is.na(c3w$C3)]))
wvec <- function(des, parc) {
  col <- make.names(paste0(des, "|", parc))
  stopifnot(col %in% names(GW))
  v <- GW[[col]]; z <- setNames(v, GW$gene); z[!is.na(z)]
}

letters_i <- 0
nxt <- function() { letters_i <<- letters_i + 1; letters[letters_i] }
panels <- list()
for (des in DESIGNS) {
  short <- c("dCT + CT" = "dCT+CT", "dCT alone" = "dCT alone",
             "CT + dCT + T1T2 + dT1T2" = "4 features")[[des]]
  # scores vs C3, DK then HCP
  sdk <- g(des, "DK", "scores", "C3")
  panels[[length(panels) + 1]] <- scat(
    c3dk[names(scores_of(des, "DK")), "C3"], scores_of(des, "DK"),
    "AHBA C3 score", sprintf("%s score", short),
    sprintf("%s   DK scores vs C3", nxt()), sprintf("rho=%.2f, %s", sdk$rho, pf(sdk$p)), "DK")
  if (paste(des, "HCP") %in% have) {
    shc <- g(des, "HCP", "scores", "C3")
    panels[[length(panels) + 1]] <- scat(
      hcpS1[names(scores_of(des, "HCP")), "C3"], scores_of(des, "HCP"),
      "AHBA C3 score", sprintf("%s score", short),
      sprintf("%s   HCP scores vs C3", nxt()), sprintf("rho=%.2f, %s", shc$rho, pf(shc$p)), "HCP")
  } else {
    panels[[length(panels) + 1]] <- note("HCP-MMP: not fitted\n(no T1w/T2w parcellation)", 2.0)
  }
  # gene weights vs C3, DK then HCP
  wdk <- g(des, "DK", "weights", "C3"); wdk_c1 <- g(des, "DK", "weights", "C1")
  zdk <- wvec(des, "DK"); sh <- intersect(names(zdk), names(c3vec))
  stopifnot(length(sh) == wdk$n)
  panels[[length(panels) + 1]] <- hexp(
    c3vec[sh], zdk[sh], "C3 gene weight", sprintf("%s gene Z", short),
    sprintf("%s   DK weights vs C3", nxt()),
    sprintf("vs C3 %.2f; vs C1 %.2f", wdk$rho, wdk_c1$rho))
  if (paste(des, "HCP") %in% have) {
    whc <- g(des, "HCP", "weights", "C3"); whc_c1 <- g(des, "HCP", "weights", "C1")
    zhc <- wvec(des, "HCP"); sh2 <- intersect(names(zhc), names(c3vec))
    stopifnot(length(sh2) == whc$n)
    panels[[length(panels) + 1]] <- hexp(
      c3vec[sh2], zhc[sh2], "C3 gene weight", sprintf("%s gene Z", short),
      sprintf("%s   HCP weights vs C3", nxt()),
      sprintf("vs C3 %.2f; vs C1 %.2f", whc$rho, whc_c1$rho))
  } else {
    panels[[length(panels) + 1]] <- note("HCP-MMP: not fitted\n(no T1w/T2w parcellation)", 2.0)
  }
}

des_lab <- function(des) {
  d <- cmp(des, "DK")
  # one line for the design name (the column is wide enough), two for the stats,
  # with enough separation that the blocks cannot collide at any row height
  txt <- sub("CT \\+ dCT \\+ T1T2 \\+ dT1T2", "four features", des)
  ggplot() +
    annotate("text", 0, 0.48, label = txt, size = 2.2, fontface = "bold", hjust = 0.5) +
    annotate("text", 0, -0.42, label = sprintf("%s\nspin p = %.3f", d$component, d$p_spin),
             size = 1.75, colour = "grey40", hjust = 0.5) +
    xlim(-1, 1) + ylim(-1, 1) + coord_cartesian(clip = "off") +
    theme_void() + theme(plot.margin = margin(0, 1, 0, 1))
}
srow <- function(i) (des_lab(DESIGNS[i]) | panels[[4*i-3]] | panels[[4*i-2]] |
                       panels[[4*i-1]] | panels[[4*i]]) + plot_layout(widths = c(0.5, 1, 1, 1, 1))

# ------------------------------------------------------------------ caption ---
d1 <- cmp(DESIGNS[1], "DK"); h1 <- cmp(DESIGNS[1], "HCP"); d2 <- cmp(DESIGNS[2], "DK")
d3 <- cmp(ALL_DESIGNS[3], "DK")   # computed even when not plotted
four_bullets <- if (SHOW_FOUR) c(
  sprintf("\u2022 Third design: CT + dCT + T1w/T2w + dT1w/T2w (four columns, four components; the component shown is chosen by %s).", d3$selection),
  "\u2022 The four-feature design is DK-only: ABCD tabulates T1w/T2w in Desikan space, and the locally-derived HCP-MMP parcellation covers thickness only, so a 137-parcel\n  T1w/T2w map would need a new surface run (tools/hcp_backfill.sbatch does thickness).") else
  sprintf("\u2022 A third design (CT + dCT + T1w/T2w + dT1w/T2w, the closest analogue of the NSPN PNAS Y matrix) is computed by 17_design_grid.py but not shown: its C3-aligned\n  component explains %.0f%% of the cross-covariance and is not spin-significant (p = %.3f), and it is DK-only because there is no HCP-MMP T1w/T2w parcellation. Its\n  numbers are in design_grid_*.tsv and the README table; SHOW_FOUR=TRUE puts the row back on the figure.",
          100 * d3$cov_explained, d3$p_spin)
methods <- paste(
  "Methods.",
  "\u2022 Y designs: dCT + CT (PLS2, the lead signature) and dCT alone (a single Y column, so PLS1 is simply the vector of gene\u2013map correlations and nothing absorbs the\n  baseline-thickness gradient).",
  paste(four_bullets, collapse = "\n"),
  sprintf("\u2022 Spin p is each component's own singular value against 5,000 rotations of the complete map: %.3f (dCT+CT, DK), %.3f (dCT+CT, HCP), %.3f (dCT alone, DK), %.3f (dCT alone, HCP).\n  Only the dCT+CT fits clear 0.05 \u2014 the single-Y fits are usable gene rankings but carry no independent spatial claim.",
          d1$p_spin, h1$p_spin, d2$p_spin, cmp(DESIGNS[2], "HCP")$p_spin),
  sprintf("\u2022 Gene weights are bootstrap Z over %s resamples; n = %s genes (DK, AHBA_updated ds25 matrix) and %s (HCP-MMP, abagen-data at the shipped C1\u2013C3 gene list).",
          "1,000", format(d1$n_genes, big.mark = ","), format(h1$n_genes, big.mark = ",")),
  sprintf("\u2022 Imaging: %s children with \u22652 QC-passing visits, true 7.0 tabulated tables; both parcellations run on the same sample.",
          format(PROV$n_subjects[1], big.mark = ",")),
  "\u2022 Everything is in the thinning orientation (each component multiplied by \u2212sign of its dCT salience): positive = expressed more where adolescent thinning is faster.",
  "  Grey parcels have no AHBA donor coverage. Cell-class profiles of these rankings are in fig_celltypes.png.",
  sep = "\n")

# one scatter row per design shown; `/` flattens, so heights must count every row
srows <- Reduce(`/`, lapply(seq_along(DESIGNS), srow))
methods_panel <- ggplot() +
  annotate("text", x = 0, y = 1, label = methods, hjust = 0, vjust = 1,
           size = 2.05, lineheight = 1.42, colour = "grey25") +
  xlim(0, 1) + ylim(0, 1) + coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(4, 2, 0, 2))

fig <- (r1 / r2 / srows / methods_panel) +
  # The brain panels use coord_fixed, so their aspect caps how tall their row can
  # draw; if their relative weight is small the cap shrinks the WHOLE grid
  # proportionally and the surplus shows up as dead space above the caption.
  # Weighting the brain rows near 1 makes the scatter rows the binding constraint.
  # Brains (coord_fixed, two views wide) need ~0.4x the height of a scatter row;
  # giving them more only pads their row.
  plot_layout(heights = c(0.42, 0.42, rep(1, length(DESIGNS)), 0.70)) +
  plot_annotation(
    title = if (SHOW_FOUR)
      "One axis, three ways of asking for it \u2014 and the static map is what makes the difference"
      else "One axis, two ways of asking for it \u2014 and the static map is what makes the difference",
    subtitle = sprintf("Adding baseline CT to the Y matrix (row 1) recovers AHBA C3 at DK resolution (scores rho %.2f, p_spin \u2264 0.001); thinning rate alone (row 2) does not\n(%.2f, %s), because its gene weights load on the static gradient C1 as heavily as on C3 (%.2f vs %.2f, against %.2f vs %.2f for dCT + CT).\nAt 137 parcels that confound largely disappears (panels f and h): rho %.2f vs C3 and only %.2f vs C1, without any static map in Y.",
                       g(DESIGNS[1], "DK", "scores", "C3")$rho,
                       g(DESIGNS[2], "DK", "scores", "C3")$rho, pf(g(DESIGNS[2], "DK", "scores", "C3")$p),
                       g(DESIGNS[2], "DK", "weights", "C1")$rho, g(DESIGNS[2], "DK", "weights", "C3")$rho,
                       g(DESIGNS[1], "DK", "weights", "C1")$rho, g(DESIGNS[1], "DK", "weights", "C3")$rho,
                       g(DESIGNS[2], "HCP", "weights", "C3")$rho, g(DESIGNS[2], "HCP", "weights", "C1")$rho),
    theme = theme(plot.title = element_text(size = 9.6, face = "bold"),
                  plot.subtitle = element_text(size = 7.3, colour = "grey25",
                                               margin = margin(b = 4))))

# title+caption take a fixed ~2.3 in; the panel stack scales with its own height weights
# measured fixed rows (title, subtitle, caption, margins) take ~3.0 in on this
# device; the rest is shared by the flexible panel rows at ~1.45 in per unit.
H <- 0.9 + 1.45 * (2 * 0.42 + length(DESIGNS) + 0.70)
ggsave(file.path(FIG, "fig_signature_designs.png"), fig, width = 8.6, height = H, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_signature_designs.png"), "\n")
