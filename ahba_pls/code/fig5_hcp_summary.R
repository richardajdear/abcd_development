#!/usr/bin/env Rscript
# fig5_hcp_summary.R -- one-slide summary of the HCP-MMP arm (Y = dCT + CT,
# 137 parcels): maps, what each PLS component tracks, agreement with AHBA C1/C3,
# SCZ/MDD gene-property enrichment, and cell-class / layer enrichment.
#
# Reads only saved tables and fits nothing:
#   results/hcp_summary_{maps.csv,map_pairs,gene_pairs,gene_weights,sets}.tsv  (22_hcp_summary.py)
#   results/magma_disorder_panel.tsv                                       (21_magma_disorder_panel.py)
#   results/hcp_pls_components.tsv, hcp_run_provenance.tsv, data/hcp_polygons.csv
# Writes figures/fig_hcp_summary.png (16:9 slide)
#
# Conventions: PLS1 is oriented to align with AHBA C1 (positive = thinner
# baseline cortex); PLS2 to faster thinning. Non-significant panels are drawn at
# alpha 0.3 (b, c: spin p >= 0.05; d: p >= 0.05; e: BH q >= 0.05).

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})

# Usage: Rscript code/fig5_hcp_summary.R [3d_ds5|ds5|base]   (AHBA matrix, see 22_hcp_summary.py)
VARIANT <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(VARIANT)) VARIANT <- "3d_ds5"
VLAB <- c(`3d_ds5` = "hcp_3d_ds5: \u22653-donor region filter + DS5 gene filter",
          ds5 = "hcp_ds5: no region filter, DS5 gene filter",
          base = "hcp_base: no region or gene filter")
stopifnot(VARIANT %in% names(VLAB))
P <- if (VARIANT == "3d_ds5") "hcp_summary" else paste0("hcp_summary_", VARIANT)
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "figures")
rd <- function(x) { f <- file.path(RES, paste0(P, "_", x)); if (grepl("csv$", x)) read.csv(f) else read.delim(f) }
maps <- rd("maps.csv"); MP <- rd("map_pairs.tsv"); GP <- rd("gene_pairs.tsv"); GW <- rd("gene_weights.tsv")
ST <- rd("sets.tsv"); MG <- rd("magma.tsv"); CMP <- rd("components.tsv")
poly <- read.csv(file.path(ROOT, "data", "hcp_polygons.csv")) |> filter(view %in% c("lateral", "medial"))

cm <- function(k, col) CMP[CMP$component == k, col]
stopifnot(cm("PLS1", "sal_CT") > 0.9, cm("PLS2", "sal_dCT") > 0.9, all(CMP$variant == VARIANT))
pf <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
mpair <- function(x, y) { r <- MP[MP$x == x & MP$y == y, ]; stopifnot(nrow(r) == 1); r }
gpair <- function(x, y) { r <- GP[GP$x == x & GP$y == y, ]; stopifnot(nrow(r) == 1); r }
stopifnot(mpair("C1", "PLS1")$rho > 0.9)                     # orientation check
FADE <- 0.3
# (parcel systems stay in the maps table -- hcp_parcel_systems.csv -- but are not plotted)

TXT <- 6.4
base_theme <- theme_classic(base_size = TXT) +
  theme(plot.title = element_text(size = TXT + 0.8, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = TXT - 0.3, colour = "grey30", hjust = 0, margin = margin(b = 2)),
        strip.background = element_blank(), strip.text = element_text(size = TXT, face = "bold"),
        axis.text = element_text(size = TXT - 0.8), legend.text = element_text(size = TXT - 0.6),
        legend.title = element_text(size = TXT - 0.4), plot.margin = margin(2, 4, 2, 4))

# ------------------------------------------------------------------ a: maps ----
brain <- function(v, title, fill_scale) {
  vals <- setNames(maps[[v]], maps$label)
  d <- poly |> mutate(val = unname(vals[label]))
  ggplot(d, aes(x, y, group = interaction(view, label, group, subgroup), fill = val)) +
    geom_polygon(colour = "grey35", linewidth = 0.05) + coord_fixed(expand = FALSE) + fill_scale +
    labs(title = title) + theme_void(base_size = TXT) +
    theme(plot.title = element_text(size = TXT + 0.6, face = "bold", hjust = 0.5),
          legend.position = "bottom", legend.key.height = unit(3, "pt"), legend.key.width = unit(22, "pt"),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-2, 0, 0, 0),
          plot.margin = margin(1, 3, 2, 3))
}
div <- function(v) { l <- max(abs(maps[[v]]), na.rm = TRUE) * c(-1, 1)
  scale_fill_distiller(palette = "RdBu", limits = l, na.value = "grey82", breaks = pretty_breaks(3), name = NULL) }
# white-anchored: CT white (thinnest) -> blue (thickest); dCT white (0 mm/yr, no
# thinning) -> red (fastest thinning), i.e. warm = high thinning rate
# (ColorBrewer Blues / Reds stops; a two-stop gradient washes out the mid-range)
ct_scale  <- scale_fill_gradientn(colours = c("white", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
                                  na.value = "grey82", breaks = pretty_breaks(3), name = NULL)
dct_scale <- scale_fill_gradientn(colours = c("#67000d", "#cb181d", "#fb6a4a", "#fcbba1", "white"),
                                  limits = c(min(maps$dCT), 0), na.value = "grey82",
                                  breaks = pretty_breaks(3), name = NULL)
pa <- (brain("CT", "Baseline thickness (CT, mm)", ct_scale) |
       brain("PLS1", sprintf("PLS1 (%.0f%% of covariance)", 100 * cm("PLS1", "cov_explained")), div("PLS1")) |
       brain("C1", "AHBA C1", div("C1"))) /
      (brain("dCT", "Thinning rate (dCT, mm/yr)", dct_scale) |
       brain("PLS2", sprintf("PLS2 (%.0f%% of covariance)", 100 * cm("PLS2", "cov_explained")), div("PLS2")) |
       brain("C3", "AHBA C3", div("C3")))
pa <- wrap_elements(full = pa + plot_annotation(title = "a  Maps: PLS1 resembles C1, PLS2 resembles C3",
                                                  theme = theme(plot.title = element_text(size = TXT + 0.8, face = "bold"))))

# --------------------------------------------------- b: CT vs dCT -------------
ctd <- mpair("CT", "dCT"); ctd_sig <- ctd$p_spin < 0.05
pctd <- ggplot(maps, aes(CT, dCT)) +
  geom_point(size = 1.05, stroke = 0, colour = "grey20", alpha = if (ctd_sig) 1 else FADE) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4,
              colour = if (ctd_sig) "grey10" else alpha("grey10", FADE)) +
  annotate("text", -Inf, Inf, label = sprintf("rho = %.2f, %s", ctd$rho, sub("p ", "p_spin ", pf(ctd$p_spin))),
           hjust = -0.06, vjust = 1.2, size = 1.9, fontface = "bold", alpha = if (ctd_sig) 1 else FADE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = "CT (mm)", y = "dCT (mm/yr)", title = "b  Thickness vs thinning") +
  base_theme

# ------------------------------------------ c, d: region scatters by system ----
long_maps <- function(xs) maps |> filter(!is.na(PLS1)) |>
  pivot_longer(all_of(xs), names_to = "xv", values_to = "x") |>
  pivot_longer(c(PLS1, PLS2), names_to = "yv", values_to = "y") |>
  rowwise() |> mutate(sig = mpair(xv, yv)$p_spin < 0.05) |> ungroup() |>
  mutate(xv = factor(xv, xs), yv = factor(yv, c("PLS1", "PLS2")))
lab_maps <- function(xs) expand.grid(xv = xs, yv = c("PLS1", "PLS2"), stringsAsFactors = FALSE) |>
  rowwise() |> mutate(r = mpair(xv, yv)$rho, p = mpair(xv, yv)$p_spin,
                      txt = sprintf("rho = %.2f, %s", r, sub("p ", "p_spin ", pf(p))), sig = p < 0.05) |>
  ungroup() |> mutate(xv = factor(xv, xs), yv = factor(yv, c("PLS1", "PLS2")))
scat <- function(xs, xl) {
  d <- long_maps(xs); lab <- lab_maps(xs)
  ggplot(d, aes(x, y)) +
    geom_point(aes(alpha = sig), size = 1.05, stroke = 0, colour = "grey20") +
    # alpha on geom_smooth reaches only the ribbon, so the line fades through its colour
    geom_smooth(data = \(z) filter(z, sig), method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey10", linewidth = 0.4) +
    geom_smooth(data = \(z) filter(z, !sig), method = "lm", formula = y ~ x, se = FALSE,
                colour = alpha("grey10", FADE), linewidth = 0.4) +
    geom_text(data = lab, aes(-Inf, Inf, label = txt, alpha = sig), hjust = -0.06, vjust = 1.2, size = 1.9,
              fontface = "bold", inherit.aes = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    facet_grid(yv ~ xv, scales = "free", switch = "both") + labs(x = xl, y = "PLS score") + base_theme +
    theme(strip.placement = "outside", panel.spacing = unit(6, "pt"))
}
pb <- scat(c("CT", "dCT"), "imaging map") + theme(legend.position = "none") +
  labs(title = "b  What each component tracks")
pc1 <- scat(c("C1", "C3"), "AHBA score") + theme(legend.position = "none") +
  labs(title = "c  Regions vs AHBA")

gl <- GW |> pivot_longer(c(C1, C3), names_to = "xv", values_to = "x") |>
  pivot_longer(c(PLS1, PLS2), names_to = "yv", values_to = "y") |>
  rowwise() |> mutate(sig = mpair(xv, yv)$p_spin < 0.05) |> ungroup() |>
  mutate(xv = factor(xv, c("C1", "C3")), yv = factor(yv, c("PLS1", "PLS2")))
glab <- expand.grid(xv = c("C1", "C3"), yv = c("PLS1", "PLS2"), stringsAsFactors = FALSE) |>
  rowwise() |> mutate(txt = sprintf("rho = %.2f", gpair(xv, yv)$rho), sig = mpair(xv, yv)$p_spin < 0.05) |>
  ungroup() |> mutate(xv = factor(xv, c("C1", "C3")), yv = factor(yv, c("PLS1", "PLS2")))
pc2 <- ggplot(gl, aes(x, y)) + geom_hex(aes(alpha = sig), bins = 36, linewidth = 0) +
  scale_fill_gradient(low = "grey85", high = "grey5", trans = "log10", guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  geom_text(data = glab, aes(-Inf, Inf, label = txt, alpha = sig), hjust = -0.06, vjust = 1.2, size = 1.9,
            fontface = "bold", inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
  facet_grid(yv ~ xv, scales = "free", switch = "both") +
  labs(title = "   \u2026and genes vs AHBA",
       x = "AHBA gene weight", y = "PLS gene Z") + base_theme +
  theme(strip.placement = "outside", panel.spacing = unit(6, "pt"))

# ------------------------------------------------------------ d: MAGMA ------
VEC <- c(ABCD_PLS1_HCP = "PLS1", AHBA_C1 = "AHBA C1", ABCD_PLS2_HCP = "PLS2", AHBA_C3 = "AHBA C3")
GA  <- c(SCZ25_META = "SCZ (2025, multi-ancestry)", MDD_div = "MDD (2025, multi-ancestry)",
         ASD = "ASD (SPARK + iPSYCH + PGC)")
GA  <- GA[names(GA) %in% MG$gene_analysis]           # ASD appears once its gene analysis is local
d <- MG |> filter(VARIABLE %in% names(VEC), gene_analysis %in% names(GA)) |>
  mutate(v = factor(VEC[VARIABLE], rev(VEC)), g = factor(GA[gene_analysis], GA),
         sig = P < 0.05, lab = ifelse(P < 0.001, sprintf("%.0e", P), sprintf("%.3f", P)))
stopifnot(nrow(d) == 4 * length(GA))
slope <- MG |> filter(VARIABLE %in% names(VEC), grepl("^global_slope_", gene_analysis))
stopifnot(nrow(slope) == 16)
cat(sprintf("[%s] global slope: min p = %.3f over %d tests (not plotted)\n", VARIANT, min(slope$P), nrow(slope)))
pd <- ggplot(d, aes(BETA_STD, v, alpha = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = BETA_STD - 1.96 * se_std, xmax = BETA_STD + 1.96 * se_std), width = 0,
                orientation = "y", linewidth = 0.45, colour = "#b2182b") +
  geom_point(size = 1.4, colour = "#b2182b") +
  geom_text(aes(label = lab), vjust = -0.75, size = 1.85, colour = "#b2182b") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  facet_wrap(~g, nrow = 1) + labs(x = "MAGMA gene-property \u03b2 (standardised, 95% CI)", y = NULL,
    title = "e  Genetic risk") +
  base_theme + theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
                     panel.spacing.x = unit(10, "pt"))

# ---------------------------------------------- e: cell classes and layers ----
SET_ORD <- c("Neuro-Ex", "Neuro-In", "Astro", "Micro", "Oligo", "OPC", "Endo", "Per",
             "L1", "L2", "L3", "L4", "L5", "L6", "WM")
ZMAX <- 20
e <- ST |> filter(set != "Neuro") |>
  mutate(vector = factor(recode(vector, C1 = "AHBA C1", C3 = "AHBA C3"), c("PLS1", "AHBA C1", "PLS2", "AHBA C3")),
         set = factor(set, SET_ORD), kind = factor(kind, c("cell class", "layer"), c("cell classes", "cortical layers")),
         zc = pmax(pmin(z, ZMAX), -ZMAX), az = abs(zc), sig = q_bh < 0.05)
stopifnot(!any(is.na(e$set)))
pe <- ggplot(e, aes(set, vector)) +
  geom_point(aes(size = az, fill = zc, alpha = sig), shape = 22, colour = "grey30", stroke = 0.15) +
  scale_size_area(max_size = 9.5, limits = c(0, ZMAX), breaks = c(5, 10, 20), labels = c("5", "10", "\u226520"), name = "|z|") +
  scale_fill_distiller(palette = "RdBu", limits = c(-ZMAX, ZMAX), breaks = c(-ZMAX, 0, ZMAX),
                       labels = c("\u2212", "0", "+"), name = "direction") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), labels = c(`TRUE` = "q < 0.05", `FALSE` = "n.s."),
                     name = NULL) +
  scale_y_discrete(limits = rev) + facet_grid(~kind, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "d  Cell classes and cortical layers") +
  guides(size = guide_legend(order = 1, override.aes = list(fill = "grey60")),
         fill = guide_colourbar(order = 2, barwidth = unit(4, "pt"), barheight = unit(22, "pt")),
         alpha = guide_legend(order = 3, override.aes = list(size = 3, fill = "#b2182b"))) +
  base_theme + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     panel.grid.major = element_line(colour = "grey94", linewidth = 0.25),
                     legend.key.size = unit(8, "pt"), legend.spacing.y = unit(2, "pt"))

# ------------------------------------------------------------- assemble -------
# methods live in the README (HCP-MMP summary slide section), not on the slide
# pctd (CT vs dCT) is still built above but not placed; the CBCL panel will
# take the slot beside the MAGMA panel once its analysis is settled.
row1 <- (pa | pb) + plot_layout(widths = c(1.55, 1))
row2 <- (pc1 | pc2 | pe) + plot_layout(widths = c(1, 1, 1.9))
row3 <- (pd | plot_spacer()) + plot_layout(widths = c(1, 1.55))
fig <- (row1 / row2 / row3) + plot_layout(heights = c(1.15, 0.95, 0.72))
fig <- fig + plot_annotation(
    title = "HCP-MMP: the ABCD thinning component is the AHBA C3 axis, and it carries SCZ and MDD risk",
    subtitle = sprintf(paste0("PLS of AHBA expression against adolescent thickness and thinning separates a static component (PLS1 = AHBA C1, thin cortex) from a thinning component (PLS2 = AHBA C3).\n",
                              "PLS2 is enriched for SCZ and MDD genetic risk and for neuronal and upper-layer genes, but C3 remains the stronger disorder ranking.   ",
                              "AHBA matrix: %s (%d parcels \u00d7 %s genes)."),
                       VLAB[[VARIANT]], cm("PLS1", "n_parcels"), format(cm("PLS1", "n_genes"), big.mark = ",")),
    theme = theme(plot.title = element_text(size = 10.5, face = "bold"),
                  plot.subtitle = element_text(size = 7.3, colour = "grey25", lineheight = 1.15, margin = margin(b = 4))))
OUT <- file.path(FIG, if (VARIANT == "3d_ds5") "fig_hcp_summary.png" else paste0("fig_hcp_summary_", VARIANT, ".png"))
ggsave(OUT, fig, width = 13.33, height = 7.9, dpi = 300, bg = "white")
cat("wrote", OUT, "\n")
