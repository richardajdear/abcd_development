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

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})

ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "figures")
maps <- read.csv(file.path(RES, "hcp_summary_maps.csv"))
MP   <- read.delim(file.path(RES, "hcp_summary_map_pairs.tsv"))
GP   <- read.delim(file.path(RES, "hcp_summary_gene_pairs.tsv"))
GW   <- read.delim(file.path(RES, "hcp_summary_gene_weights.tsv"))
ST   <- read.delim(file.path(RES, "hcp_summary_sets.tsv"))
MG   <- read.delim(file.path(RES, "magma_disorder_panel.tsv"))
CMP  <- read.delim(file.path(RES, "hcp_pls_components.tsv")) |> filter(option == "hcp_opt2_dCT_CT")
PRV  <- read.delim(file.path(RES, "hcp_run_provenance.tsv"))
poly <- read.csv(file.path(ROOT, "data", "hcp_polygons.csv")) |> filter(view %in% c("lateral", "medial"))

cm <- function(k, col) CMP[CMP$component == k, col]
stopifnot(cm("PLS1", "sal_CT") > 0.9, cm("PLS2", "sal_dCT") > 0.9)
pf <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
mpair <- function(x, y) { r <- MP[MP$x == x & MP$y == y, ]; stopifnot(nrow(r) == 1); r }
gpair <- function(x, y) { r <- GP[GP$x == x & GP$y == y, ]; stopifnot(nrow(r) == 1); r }
TXT <- 6.4
base_theme <- theme_classic(base_size = TXT) +
  theme(plot.title = element_text(size = TXT + 0.8, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = TXT - 0.3, colour = "grey30", hjust = 0, margin = margin(b = 2)),
        strip.background = element_blank(), strip.text = element_text(size = TXT, face = "bold"),
        axis.text = element_text(size = TXT - 0.8), legend.text = element_text(size = TXT - 0.8),
        legend.title = element_text(size = TXT - 0.5), plot.margin = margin(2, 4, 2, 4))

# ------------------------------------------------------------------ a: maps ----
brain <- function(v, title, sub, diverging = TRUE) {
  vals <- setNames(maps[[v]], maps$label)
  d <- poly |> mutate(val = unname(vals[label]))
  lim <- if (diverging) max(abs(d$val), na.rm = TRUE) * c(-1, 1) else range(d$val, na.rm = TRUE)
  ggplot(d, aes(x, y, group = interaction(view, label, group, subgroup), fill = val)) +
    geom_polygon(colour = "grey35", linewidth = 0.06) + coord_fixed(expand = FALSE) +
    (if (diverging) scale_fill_distiller(palette = "RdBu", limits = lim, na.value = "grey80",
                                         breaks = pretty_breaks(3), name = NULL)
     else scale_fill_viridis_c(option = "magma", na.value = "grey80", breaks = pretty_breaks(3), name = NULL)) +
    labs(title = title, subtitle = sub) + theme_void(base_size = TXT) +
    theme(plot.title = element_text(size = TXT + 0.8, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = TXT - 0.4, colour = "grey30", hjust = 0.5),
          legend.position = "bottom", legend.key.height = unit(3, "pt"), legend.key.width = unit(20, "pt"),
          legend.margin = margin(0, 0, 0, 0), plot.margin = margin(1, 3, 1, 3))
}
pa <- brain("CT", "Baseline thickness (CT)", "mm, age 10", FALSE) |
  brain("dCT", "Thinning rate (dCT)", sprintf("mm/yr; %d of %d parcels thin", sum(maps$dCT < 0), nrow(maps)), FALSE) |
  brain("PLS1", "PLS1 score", sprintf("%.0f%% of covariance; tracks CT", 100 * cm("PLS1", "cov_explained"))) |
  brain("PLS2", "PLS2 score", sprintf("%.0f%% of covariance; tracks thinning", 100 * cm("PLS2", "cov_explained")))

# -------------------------------------------- scatter helpers (b and c) -------
scat <- function(d, xl, yl, lab) {
  ggplot(d, aes(x, y)) + geom_point(size = 0.55, colour = "grey25", alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#b2182b", linewidth = 0.4) +
    geom_text(data = lab, aes(-Inf, Inf, label = txt), hjust = -0.06, vjust = 1.2, size = 1.95,
              fontface = "bold", inherit.aes = FALSE) +
    # 30% headroom so the annotation never sits on the points
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    facet_grid(yv ~ xv, scales = "free", switch = "both") + labs(x = xl, y = yl) + base_theme +
    theme(strip.placement = "outside", panel.spacing = unit(6, "pt"))
}
long_maps <- function(xs) maps |> filter(!is.na(PLS1)) |>
  pivot_longer(all_of(xs), names_to = "xv", values_to = "x") |>
  pivot_longer(c(PLS1, PLS2), names_to = "yv", values_to = "y") |>
  mutate(xv = factor(xv, xs), yv = factor(yv, c("PLS1", "PLS2")))
lab_maps <- function(xs) expand.grid(xv = xs, yv = c("PLS1", "PLS2"), stringsAsFactors = FALSE) |>
  rowwise() |> mutate(txt = sprintf("rho = %.2f, %s", mpair(xv, yv)$rho, sub("p ", "p_spin ", pf(mpair(xv, yv)$p_spin)))) |>
  ungroup() |> mutate(xv = factor(xv, xs), yv = factor(yv, c("PLS1", "PLS2")))

pb <- scat(long_maps(c("CT", "dCT")), "imaging map", "PLS score", lab_maps(c("CT", "dCT"))) +
  labs(title = "b  What each component tracks",
       subtitle = "PLS1 = where cortex is thick; PLS2 = where it thins fastest")

pc1 <- scat(long_maps(c("C1", "C3")), "AHBA score", "PLS score", lab_maps(c("C1", "C3"))) +
  labs(title = "c  Regions vs AHBA", subtitle = "PLS1 is C1 (sign-flipped); PLS2 is C3")

gl <- GW |> pivot_longer(c(C1, C3), names_to = "xv", values_to = "x") |>
  pivot_longer(c(PLS1, PLS2), names_to = "yv", values_to = "y") |>
  mutate(xv = factor(xv, c("C1", "C3")), yv = factor(yv, c("PLS1", "PLS2")))
glab <- expand.grid(xv = c("C1", "C3"), yv = c("PLS1", "PLS2"), stringsAsFactors = FALSE) |>
  rowwise() |> mutate(txt = sprintf("rho = %.2f", gpair(xv, yv)$rho)) |> ungroup() |>
  mutate(xv = factor(xv, c("C1", "C3")), yv = factor(yv, c("PLS1", "PLS2")))
pc2 <- ggplot(gl, aes(x, y)) + geom_hex(bins = 38, linewidth = 0) +
  scale_fill_gradient(low = "grey88", high = "grey10", trans = "log10", guide = "none") +
  geom_text(data = glab, aes(-Inf, Inf, label = txt), hjust = -0.06, vjust = 1.2, size = 1.95,
            fontface = "bold", inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
  facet_grid(yv ~ xv, scales = "free", switch = "both") +
  labs(title = "   \u2026and genes vs AHBA", subtitle = sprintf("same pattern in %s gene weights (bootstrap Z)",
                                                     format(GP$n[1], big.mark = ",")),
       x = "AHBA gene weight", y = "PLS gene Z") + base_theme +
  theme(strip.placement = "outside", panel.spacing = unit(6, "pt"))

# ------------------------------------------------------------ d: MAGMA ------
VEC <- c(ABCD_PLS1_HCP = "PLS1", ABCD_PLS2_HCP = "PLS2", AHBA_C1 = "AHBA C1", AHBA_C3 = "AHBA C3")
GA  <- c(SCZ25_META = "SCZ (2025, multi-ancestry)", MDD_EUR = "MDD (2025, EUR)")
d <- MG |> filter(VARIABLE %in% names(VEC), gene_analysis %in% names(GA)) |>
  mutate(v = factor(VEC[VARIABLE], rev(VEC)), g = factor(GA[gene_analysis], GA),
         sig = P < 0.05, lab = ifelse(P < 0.001, sprintf("%.0e", P), sprintf("%.3f", P)))
stopifnot(nrow(d) == 8)
slope <- MG |> filter(VARIABLE %in% names(VEC), grepl("^slope_", gene_analysis))
stopifnot(nrow(slope) == 16)
pd <- ggplot(d, aes(BETA_STD, v, colour = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = BETA_STD - 1.96 * se_std, xmax = BETA_STD + 1.96 * se_std), width = 0, orientation = "y", linewidth = 0.45) +
  geom_point(size = 1.4) + geom_text(aes(label = lab), vjust = -0.75, size = 1.85, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey55"), guide = "none") +
  facet_wrap(~g, nrow = 1) + labs(x = "MAGMA gene-property \u03b2 (standardised, 95% CI)", y = NULL,
    title = "d  Genetic risk", subtitle = "PLS2 and C3 carry SCZ and MDD risk; C1 carries SCZ only") +
  base_theme + theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25))

# ---------------------------------------------- e: cell classes and layers ----
SET_ORD <- c("Neuro-Ex", "Neuro-In", "Astro", "Micro", "Oligo", "OPC", "Endo", "Per",
             "L1", "L2", "L3", "L4", "L5", "L6", "WM")
e <- ST |> filter(set != "Neuro") |>
  mutate(vector = factor(recode(vector, C1 = "AHBA C1", C3 = "AHBA C3"), c("PLS1", "PLS2", "AHBA C1", "AHBA C3")),
         set = factor(set, SET_ORD), kind = factor(kind, c("cell class", "layer"), c("cell classes", "cortical layers")),
         zc = pmax(pmin(z, 15), -15), lab = sprintf("%.0f%s", z, ifelse(q_bh < 0.05, "", "\u00b7")))
stopifnot(!any(is.na(e$set)))
pe <- ggplot(e, aes(set, vector, fill = zc)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = paste0(sprintf("%.0f", z), ifelse(q_bh < 0.05, "*", "")), colour = abs(zc) > 9),
            size = 1.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_distiller(palette = "RdBu", limits = c(-15, 15), name = "z", breaks = c(-15, 0, 15),
                       labels = c("\u2264\u221215", "0", "\u226515")) +
  scale_y_discrete(limits = rev) + facet_grid(~kind, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "e  Cell classes and cortical layers",
       subtitle = "PLS2 follows C3 (neurons and L2\u20133 up; oligodendrocytes, microglia and WM down) except astrocytes; PLS1 \u2248 \u2212C1") +
  base_theme + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.key.width = unit(4, "pt"), legend.key.height = unit(14, "pt"))

# ------------------------------------------------------------- assemble -------
smin <- min(slope$P)
methods <- paste(
  sprintf("\u2022 ABCD 7.0: %s children, %s sessions; thickness intercept (CT) and slope (dCT) per parcel from one mixed model each; %d bilateral HCP-MMP parcels, %d with AHBA coverage.",
          format(PRV$n_subjects, big.mark = ","), format(PRV$n_sessions, big.mark = ","), PRV$n_parcels_pls, cm("PLS1", "n_regions")),
  sprintf("\u2022 PLS: AHBA expression (7,973 genes, C1\u2013C3 gene set) against Y = [dCT, CT]. Spin p for the singular values: PLS1 %s, PLS2 %s; bootstrap reproducibility %.2f / %.2f.",
          pf(cm("PLS1", "p_spin_singular")), pf(cm("PLS2", "p_spin_singular")), cm("PLS1", "boot_reproducibility"), cm("PLS2", "boot_reproducibility")),
  "\u2022 Orientation: PLS1 positive = thicker; PLS2 positive = faster thinning. Map rho are Spearman with spin p (5,000 rotations); gene rho over all genes.",
  sprintf("\u2022 d: MAGMA gene-property, two-sided; SCZ gene Z per ancestry on matched 1000G panels then --meta. Global-slope GWAS (GENESIS, EUR): no vector associated (min p = %.2f), not shown.", smin),
  "\u2022 e: mean gene Z of each marker set vs 20,000 random sets; * = BH q < 0.05 within vector. Cell classes: Seidlitz 2020; layers: Maynard 2021 (FDR < 0.05, t > 0).",
  sep = "\n")
methods_panel <- ggplot() + annotate("text", 0, 1, label = methods, hjust = 0, vjust = 1, size = 1.95,
                                     lineheight = 1.35, colour = "grey25") +
  xlim(0, 1) + ylim(0, 1) + coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(2, 4, 0, 4))

row1 <- wrap_elements(full = pa + plot_annotation(title = "a  Maps", theme = theme(plot.title = element_text(size = TXT + 0.8, face = "bold"))))
row2 <- pb | pc1 | pc2
row3 <- (pd | pe) + plot_layout(widths = c(1, 1.45))
fig <- (row1 / row2 / row3 / methods_panel) + plot_layout(heights = c(0.78, 1.0, 0.74, 0.34)) +
  plot_annotation(
    title = "HCP-MMP: the ABCD thinning component is the AHBA C3 axis, and it carries SCZ and MDD risk",
    subtitle = "PLS of AHBA expression against adolescent thickness and thinning separates a static component (PLS1 = AHBA C1, thick cortex) from a thinning component (PLS2 = AHBA C3). PLS2 is enriched for\nSCZ and MDD genetic risk and for neuronal and upper-layer genes \u2014 but C3 remains the stronger disorder ranking.",
    theme = theme(plot.title = element_text(size = 10.5, face = "bold"),
                  plot.subtitle = element_text(size = 7.3, colour = "grey25", lineheight = 1.15, margin = margin(b = 4))))
ggsave(file.path(FIG, "fig_hcp_summary.png"), fig, width = 13.33, height = 7.7, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_hcp_summary.png"), "\n")
