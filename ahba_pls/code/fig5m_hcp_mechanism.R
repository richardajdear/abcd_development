#!/usr/bin/env Rscript
# fig5m_hcp_mechanism.R -- the base-AHBA summary slide with a MECHANISM row:
# where SCZ/MDD gene-level risk sits on the transcriptomic axis (PLS2 / C3) vs
# where symptom-linked extra thinning happens (normally slow-thinning cortex).
# Derived from fig5_hcp_summary.R (variant base); panel c shows only the matched
# pairs (PLS1-C1, PLS2-C3); text sizes follow genetic_analysis/fig1.R (BASE = 7).
#
# Reads only saved tables and fits nothing:
#   results/hcp_summary_base_*                            (22_hcp_summary.py base)
#   results/tier_profile.tsv                              (27_tier_profile.py)
#   results/gradient_scores_decomposition.tsv, gradient_tiers.csv (26_gradient_scores.py)
#   data/hcp_polygons.csv
# Writes figures/fig_hcp_mechanism.png
#
# Conventions: PLS1 aligned with AHBA C1 (positive = thinner baseline cortex);
# PLS2 to faster thinning. Faded (alpha 0.3) = not significant (b, c: spin p;
# d: BH q; e: p; f: spin p of the feature map vs thinning rate; g: p).

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})

# Usage: Rscript code/fig5_hcp_summary.R [3d_ds5|ds5|base]   (AHBA matrix, see 22_hcp_summary.py)
VARIANT <- "base"
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

# text sizes as genetic_analysis/fig1.R: BASE 7; titles 7 bold, axis text 6, axis
# titles 6.5, legend text 6, in-panel text 5.5 pt
BASE <- 7; TXT <- BASE; GT <- (BASE - 1.5) / .pt
base_theme <- theme_classic(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0, margin = margin(b = 3)),
        strip.background = element_blank(), strip.text = element_text(size = BASE - 0.5, face = "bold"),
        axis.text = element_text(size = BASE - 1, colour = "grey20"), axis.title = element_text(size = BASE - 0.5),
        legend.text = element_text(size = BASE - 1), legend.title = element_text(size = BASE - 1),
        plot.margin = margin(2, 4, 2, 4))

# ------------------------------------------------------------------ a: maps ----
brain <- function(v, title, fill_scale) {
  vals <- setNames(maps[[v]], maps$label)
  d <- poly |> mutate(val = unname(vals[label]))
  ggplot(d, aes(x, y, group = interaction(view, label, group, subgroup), fill = val)) +
    geom_polygon(colour = "grey35", linewidth = 0.05) + coord_fixed(expand = FALSE) + fill_scale +
    labs(title = title) + theme_void(base_size = TXT) +
    theme(plot.title = element_text(size = BASE - 0.5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = BASE - 1.5),
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
                                                  theme = theme(plot.title = element_text(size = BASE, face = "bold"))))

# --------------------------------------------------- b: CT vs dCT -------------
ctd <- mpair("CT", "dCT"); ctd_sig <- ctd$p_spin < 0.05
pctd <- ggplot(maps, aes(CT, dCT)) +
  geom_point(size = 1.05, stroke = 0, colour = "grey20", alpha = if (ctd_sig) 1 else FADE) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4,
              colour = if (ctd_sig) "grey10" else alpha("grey10", FADE)) +
  annotate("text", -Inf, Inf, label = sprintf("rho = %.2f, %s", ctd$rho, sub("p ", "p_spin ", pf(ctd$p_spin))),
           hjust = -0.06, vjust = 1.2, size = GT, fontface = "bold", alpha = if (ctd_sig) 1 else FADE) +
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
    geom_text(data = lab, aes(-Inf, Inf, label = txt, alpha = sig), hjust = -0.06, vjust = 1.2, size = GT,
              fontface = "bold", inherit.aes = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    facet_grid(yv ~ xv, scales = "free", switch = "both") + labs(x = xl, y = "PLS score") + base_theme +
    theme(strip.placement = "outside", panel.spacing = unit(6, "pt"))
}
pb <- scat(c("CT", "dCT"), "imaging map") + theme(legend.position = "none") +
  labs(title = "b  What each component tracks")
# c: only the matched pairs (PLS1-C1, PLS2-C3); the off-diagonal pairs stay in
# hcp_summary_base_map_pairs.tsv / gene_pairs.tsv
PAIRS <- tibble(xv = c("C1", "C3"), yv = c("PLS1", "PLS2"), pair = c("PLS1 vs AHBA C1", "PLS2 vs AHBA C3"))
mlong <- bind_rows(maps |> filter(!is.na(PLS1)) |> transmute(pair = PAIRS$pair[1], x = C1, y = PLS1),
                   maps |> filter(!is.na(PLS2)) |> transmute(pair = PAIRS$pair[2], x = C3, y = PLS2)) |>
  filter(!is.na(x)) |> left_join(PAIRS, by = "pair") |> rowwise() |> mutate(sig = mpair(xv, yv)$p_spin < 0.05) |> ungroup()
mlab <- PAIRS |> rowwise() |> mutate(txt = sprintf("rho = %.2f\n%s", mpair(xv, yv)$rho, sub("p ", "p_spin ", pf(mpair(xv, yv)$p_spin))),
                                     sig = mpair(xv, yv)$p_spin < 0.05) |> ungroup()
pc1 <- ggplot(mlong, aes(x, y)) +
  geom_point(aes(alpha = sig), size = 0.9, stroke = 0, colour = "grey20") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey10", linewidth = 0.4) +
  geom_text(data = mlab, aes(-Inf, Inf, label = txt, alpha = sig), hjust = -0.06, vjust = 1.15, size = GT,
            fontface = "bold", lineheight = 0.9, inherit.aes = FALSE) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.42))) +
  facet_wrap(~pair, nrow = 1, scales = "free") + labs(x = "AHBA score", y = "PLS score", title = "c  Regions vs AHBA") +
  base_theme + theme(panel.spacing = unit(6, "pt"))

gl <- bind_rows(GW |> transmute(pair = "PLS1 vs AHBA C1", x = C1, y = PLS1),
                GW |> transmute(pair = "PLS2 vs AHBA C3", x = C3, y = PLS2)) |> filter(!is.na(x), !is.na(y))
glab <- PAIRS |> rowwise() |> mutate(txt = sprintf("rho = %.2f", gpair(xv, yv)$rho)) |> ungroup()
pc2 <- ggplot(gl, aes(x, y)) + geom_hex(bins = 34, linewidth = 0) +
  scale_fill_gradient(low = "grey85", high = "grey5", trans = "log10", guide = "none") +
  geom_text(data = glab, aes(-Inf, Inf, label = txt), hjust = -0.06, vjust = 1.3, size = GT,
            fontface = "bold", inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  facet_wrap(~pair, nrow = 1, scales = "free") +
  labs(title = "   \u2026and genes", x = "AHBA gene weight", y = "PLS gene Z") + base_theme +
  theme(panel.spacing = unit(6, "pt"))

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
  geom_text(aes(label = lab), vjust = -0.75, size = GT, colour = "#b2182b") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  facet_wrap(~g, nrow = 1) + labs(x = "MAGMA gene-property \u03b2 (standardised, 95% CI)", y = NULL,
    title = "e  SCZ / MDD genetic risk (MAGMA)") +
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
  scale_size_area(max_size = 8, limits = c(0, ZMAX), breaks = c(5, 10, 20), labels = c("5", "10", "\u226520"), name = "|z|") +
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

# ------------------------------------ f: what each third of cortex is made of ----
# thirds of the NORMATIVE thinning map (26_gradient_scores.py); features are
# parcel maps of mean z-expression (27_tier_profile.py). Tier means are
# descriptive; fading follows the spin p of each feature map vs thinning rate.
TIERCOL <- c(slow = "#2166ac", mid = "grey55", fast = "#b2182b")
TIERLAB <- c(slow = "slow-thinning third", mid = "middle third", fast = "fast-thinning third")
TP  <- read.delim(file.path(RES, "tier_profile.tsv"))
GT3 <- read.csv(file.path(RES, "gradient_tiers.csv"))
DC  <- read.delim(file.path(RES, "gradient_scores_decomposition.tsv"))
um  <- GT3 |> group_by(tier) |> summarise(r = mean(normative_thinning_um_yr), .groups = "drop")
tl  <- setNames(sprintf("%s\n(%.0f \u00b5m/yr)", TIERLAB[um$tier], um$r), um$tier)
tb  <- poly |> mutate(key = tolower(label)) |>
  left_join(GT3 |> transmute(key = tolower(label), tier), by = "key")
stopifnot(sum(!is.na(unique(tb[c("key", "tier")])$tier)) == 179)
pf1 <- ggplot(tb, aes(x, y, group = interaction(view, label, group, subgroup), fill = tier)) +
  geom_polygon(colour = "grey35", linewidth = 0.05) + coord_fixed(expand = FALSE) +
  scale_fill_manual(values = TIERCOL, breaks = names(TIERCOL), labels = tl[names(TIERCOL)], na.value = "grey82", name = NULL) +
  labs(title = "f  What each third of cortex is made of") + theme_void(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold", margin = margin(b = 3)),
        legend.position = "right", legend.text = element_text(size = BASE - 1.5, lineheight = 0.9),
        legend.key.size = unit(7, "pt"), plot.margin = margin(2, 4, 0, 4))

FEAT <- c(PLS2 = "PLS2 score", C3 = "AHBA C3 score",
          `Neuro-Ex` = "excitatory neurons", `Neuro-In` = "inhibitory neurons", Oligo = "oligodendrocytes",
          L2 = "layer 2", L3 = "layer 3", WM = "white matter",
          `SCZ top 500` = "SCZ top-500 genes", `MDD top 500` = "MDD top-500 genes")
FGRP <- c(PLS2 = "axis", C3 = "axis", `Neuro-Ex` = "cells", `Neuro-In` = "cells", Oligo = "cells",
          L2 = "layers", L3 = "layers", WM = "layers", `SCZ top 500` = "risk genes", `MDD top 500` = "risk genes")
h <- TP |> filter(feature %in% names(FEAT)) |>
  pivot_longer(c(mean_slow, mean_mid, mean_fast), names_to = "tier", values_to = "m") |>
  mutate(tier = factor(sub("mean_", "", tier), c("slow", "mid", "fast")), xi = as.integer(tier),
         f = factor(FEAT[feature], rev(FEAT)), grp = factor(FGRP[feature], c("axis", "cells", "layers", "risk genes")),
         sig = p_spin < 0.05)
stopifnot(nrow(h) == 3 * length(FEAT))
hr <- h |> distinct(f, grp, rho_with_thinning, p_spin, sig) |>
  mutate(txt = sprintf("%.2f  %s", rho_with_thinning, ifelse(p_spin < 0.001, "<.001", sub("^0", "", sprintf("%.3f", p_spin)))))
ML <- 0.8
pf2 <- ggplot(h, aes(xi, f)) +
  geom_tile(aes(fill = pmax(pmin(m, ML), -ML), alpha = sig), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.2f", m), alpha = sig), size = GT) +
  geom_text(data = hr, aes(x = 4.25, label = txt, alpha = sig), hjust = 0, size = GT) +
  # column header in the top facet only (annotate() would repeat it in every facet)
  geom_text(data = data.frame(grp = factor("axis", levels(h$grp))), aes(x = 4.25, y = Inf, label = "\u03c1 with rate, p_spin"),
            inherit.aes = FALSE, hjust = 0, vjust = -0.6, size = GT, fontface = "bold") +
  scale_x_continuous(breaks = 1:3, labels = c("slow", "middle", "fast"), limits = c(0.5, 5.9), expand = c(0, 0), position = "top") +
  scale_fill_distiller(palette = "RdBu", limits = c(-ML, ML), breaks = c(-ML, 0, ML),
                       labels = c("\u2212", "0", "+"), name = "mean z") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  coord_cartesian(clip = "off") + labs(x = "third of normative thinning", y = NULL) +
  base_theme + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                     strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "plain", size = BASE - 1.5),
                     panel.spacing.y = unit(3, "pt"), legend.key.height = unit(14, "pt"), legend.key.width = unit(5, "pt"),
                     plot.margin = margin(12, 4, 2, 4))
pF <- (pf1 / pf2) + plot_layout(heights = c(0.36, 1))

# ------------------------------------------- g: symptoms by third of cortex ----
OUTG <- c(totprob = "Total problems", pfactor = "p-factor", internal = "Internalising", anxdisord = "Anxiety (DSM)")
g <- DC |> filter(outcome %in% names(OUTG), spec %in% c("fast alone", "mid alone", "slow alone")) |>
  mutate(tier = factor(sub("thin_", "", term), c("fast", "mid", "slow")), o = factor(OUTG[outcome], rev(OUTG)),
         s = se * y_sd_beta / beta, lo = y_sd_beta - 1.96 * s, hi = y_sd_beta + 1.96 * s, sig = p < 0.05)
stopifnot(nrow(g) == 3 * length(OUTG))
jt <- DC |> filter(outcome %in% names(OUTG), spec == "fast + slow | CT, QC", term == "slow - fast") |>
  mutate(o = factor(OUTG[outcome], rev(OUTG)), sig = p < 0.05,
         txt = sprintf("slow \u2212 fast (joint)\n%+.3f, %s", y_sd_beta, ifelse(p < 0.001, "p < .001", sub("p = 0", "p = ", pf(p)))))
XR <- max(g$hi) * 1.08
pg <- ggplot(g, aes(y_sd_beta, o, colour = tier, alpha = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.5, position = position_dodge(0.62)) +
  geom_point(size = 1.5, position = position_dodge(0.62)) +
  geom_text(data = jt, aes(x = XR, y = o, label = txt, alpha = sig), inherit.aes = FALSE, hjust = 0, size = GT, lineheight = 0.9) +
  scale_colour_manual(values = TIERCOL, breaks = c("slow", "mid", "fast"), labels = TIERLAB[c("slow", "mid", "fast")], name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.02, 0.04)) +
  coord_cartesian(xlim = c(min(g$lo), XR * 1.9), clip = "off") +
  labs(x = "\u03b2 per SD of thinning rate in that third\n(SD of CBCL at 15\u201317 | baseline CBCL)", y = NULL,
       title = "g  Rising symptoms go with faster thinning in slow cortex") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  base_theme + theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                     legend.margin = margin(0, 0, 0, 0))

# ------------------------------------------------------------- assemble -------
row1 <- (pa | pb) + plot_layout(widths = c(1.55, 1))
row2 <- (pc1 | pc2 | pe) + plot_layout(widths = c(0.78, 0.78, 1.55))
row3 <- (pd | pF | pg) + plot_layout(widths = c(1, 1.1, 1))
fig <- (row1 / row2 / row3) + plot_layout(heights = c(1.1, 0.72, 1.12))
fig <- fig + plot_annotation(
    title = "The thinning axis (PLS2 = AHBA C3) carries SCZ/MDD gene-level risk; children whose symptoms rise thin faster in normally slow, PLS2-low cortex",
    subtitle = sprintf(paste0("a\u2013e: PLS of AHBA expression (%s: %d parcels \u00d7 %s genes) against ABCD baseline thickness and thinning; PLS2 = AHBA C3 is neuronal and upper-layer (d) and enriched for SCZ and MDD risk genes (e).\n",
                              "f\u2013g: the slowest-thinning third of cortex sits at the PLS2-low pole \u2014 fewer neuronal and L3 genes, more oligodendrocyte \u2014 and is where children whose symptoms rise over adolescence thin faster (g).\n",
                              "The MAGMA enrichment is a genome-wide shift across many genes: the top risk genes themselves are not concentrated in fast-thinning cortex (f, faded)."),
                       "hcp_base", cm("PLS1", "n_parcels"), format(cm("PLS1", "n_genes"), big.mark = ",")),
    theme = theme(plot.title = element_text(size = BASE + 2, face = "bold"),
                  plot.subtitle = element_text(size = BASE - 0.5, colour = "grey25", lineheight = 1.15, margin = margin(b = 4))))
OUT <- file.path(FIG, "fig_hcp_mechanism.png")
ggsave(OUT, fig, width = 12, height = 9.4, dpi = 300, bg = "white")
cat("wrote", OUT, "\n")
