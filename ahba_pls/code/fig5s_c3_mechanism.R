#!/usr/bin/env Rscript
# fig5s_c3_mechanism.R -- simplified mechanism slide, hcp_3d_ds5 fit only:
#   a  CT / PLS1 / C1 and dCT / PLS2 / C3 maps
#   b  the two strong pairs (PLS1 vs CT, PLS2 vs dCT)
#   c  Spearman matrix of the six region maps, spin p   (29_c3_mechanism_inputs.py)
#   d  MAGMA gene-property for C1 and C3 (SCZ, MDD, ASD, ALZ, EA, intelligence)
#   e  the bottom / top 20% of the AHBA C3 region score (c3_q20.csv): schematic
#      of their marker-gene expression (c3_q20_profile.tsv), the thinning rate of
#      every parcel by group (spin-tested difference), and later CBCL vs each
#      child's thinning in each group (gradient_scores_c3q20.tsv, 26)
# Reads only saved tables and fits nothing. Faded (alpha 0.3) = not significant
# (c: spin p; d, e: p). Writes figures/fig_hcp_c3_mechanism_3d_ds5.png
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})

VARIANT <- "3d_ds5"
VFILE <- c(`3d_ds5` = "hcp_3d_ds5", ds5 = "hcp_ds5", `3d` = "hcp_3d", base = "hcp_base")
VLAB <- c(`3d_ds5` = "hcp_3d_ds5: \u22653-donor region filter + DS5 gene filter",
          ds5 = "hcp_ds5: no region filter, DS5 gene filter",
          `3d` = "hcp_3d: \u22653-donor region filter, no gene filter",
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
BASE <- 9; TXT <- BASE; GT <- (BASE - 1.5) / .pt
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


# --------------------------------------------- b: the two strong pairs only -------
B2 <- tibble(xv = c("CT", "dCT"), yv = c("PLS1", "PLS2"), pair = c("PLS1 vs baseline CT", "PLS2 vs thinning rate"))
bl <- bind_rows(maps |> filter(!is.na(PLS1)) |> transmute(pair = B2$pair[1], x = CT, y = PLS1),
                maps |> filter(!is.na(PLS2)) |> transmute(pair = B2$pair[2], x = dCT, y = PLS2)) |>
  mutate(pair = factor(pair, B2$pair))
blab <- B2 |> rowwise() |> mutate(txt = sprintf("rho = %.2f\n%s", mpair(xv, yv)$rho, sub("p ", "p_spin ", pf(mpair(xv, yv)$p_spin)))) |>
  ungroup() |> mutate(pair = factor(pair, B2$pair))
stopifnot(all(sapply(seq_len(2), function(i) mpair(B2$xv[i], B2$yv[i])$p_spin < 0.05)))
pb <- ggplot(bl, aes(x, y)) +
  geom_point(size = 1.0, stroke = 0, colour = "grey20") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey10", linewidth = 0.4) +
  geom_text(data = blab, aes(-Inf, Inf, label = txt), hjust = -0.06, vjust = 1.15, size = GT, fontface = "bold",
            lineheight = 0.9, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  facet_wrap(~pair, ncol = 1, scales = "free", strip.position = "top") +
  labs(x = "imaging map (CT mm; dCT mm/yr)", y = "PLS score", title = "b  What each component tracks") +
  base_theme + theme(panel.spacing = unit(8, "pt"))

# ------------------------------------------- c: region-score correlation matrix ----
SM <- read.delim(file.path(RES, "c3_score_matrix.tsv"))
# ordered so the two clusters sit in diagonal blocks: static (CT, PLS1, C1) and thinning (dCT, PLS2, C3)
MV <- c(CT = "CT", PLS1 = "PLS1", C1 = "AHBA C1", dCT = "dCT", PLS2 = "PLS2", C3 = "AHBA C3")
stopifnot(nrow(SM) == choose(length(MV), 2), all(c(SM$x, SM$y) %in% names(MV)))
cmx <- SM |> mutate(ix = match(x, names(MV)), iy = match(y, names(MV))) |>
  transmute(i = pmax(ix, iy), j = pmin(ix, iy), rho, p_spin) |>          # lower triangle: row i, column j
  mutate(X = j, Y = length(MV) + 1 - i, sig = p_spin < 0.05, txt = sprintf("%.2f", rho))
dg <- data.frame(X = seq_along(MV), Y = rev(seq_along(MV)), lab = sub("AHBA ", "", unname(MV)))   # C1 / C3 = AHBA components
K <- length(MV)
pc <- ggplot(cmx, aes(X, Y)) +
  geom_tile(aes(fill = rho, alpha = sig), colour = "white", linewidth = 0.8) +
  geom_tile(data = dg, fill = "grey93", colour = "white", linewidth = 0.8) +
  geom_text(data = dg, aes(label = lab), size = GT + 0.2, fontface = "bold", colour = "grey15") +
  geom_text(aes(label = txt, alpha = sig, fontface = ifelse(sig, "bold", "plain")), size = GT + 0.4) +
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = K - 2.5, ymax = K + 0.5, fill = NA, colour = "grey15", linewidth = 0.6, linetype = "22") +
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 0.5, ymax = 3.5, fill = NA, colour = "grey15", linewidth = 0.6, linetype = "22") +
  annotate("text", x = 3.4, y = K + 0.25, label = "static", hjust = 1, vjust = 1, size = GT + 0.2, fontface = "italic", colour = "grey30") +
  annotate("text", x = 6.4, y = 3.25, label = "thinning", hjust = 1, vjust = 1, size = GT + 0.2, fontface = "italic", colour = "grey30") +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), breaks = c(-1, 0, 1), name = "Spearman \u03c1") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  coord_fixed(expand = FALSE) +
  labs(x = NULL, y = NULL, title = "c  Regional maps: two clusters",
       subtitle = "C1, C3 = AHBA; bold = p_spin < 0.05") +
  theme_void(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold", margin = margin(b = 2)),
        plot.subtitle = element_text(size = BASE - 1.5, colour = "grey35", margin = margin(b = 4)),
        legend.position = "bottom", legend.title = element_text(size = BASE - 1, vjust = 0.8),
        legend.text = element_text(size = BASE - 1.5),
        legend.key.width = unit(26, "pt"), legend.key.height = unit(5, "pt"), plot.margin = margin(4, 8, 4, 4))

# ---------------------------------------------- d: MAGMA, C1 / C3 only ----------
VEC <- c(AHBA_C1 = "AHBA C1", AHBA_C3 = "AHBA C3")
GA <- c(SCZ25_META = "Schizophrenia", MDD_div = "Depression", ASD = "Autism",
        ALZ = "Alzheimer's disease", EA = "Educational attainment", INT = "Intelligence")
GGRP <- c(SCZ25_META = "psychiatric", MDD_div = "psychiatric", ASD = "psychiatric",
          ALZ = "neuro-\ndegen.", EA = "cognitive", INT = "cognitive")
stopifnot(all(names(GA) %in% MG$gene_analysis))
dm <- MG |> filter(VARIABLE %in% names(VEC), gene_analysis %in% names(GA)) |>
  mutate(v = factor(VEC[VARIABLE], VEC), g = factor(GA[gene_analysis], rev(GA)),
         grp = factor(GGRP[gene_analysis], unique(GGRP)), sig = P < 0.05,
         lab = ifelse(P < 0.001, sprintf("%.0e", P), formatC(signif(P, 2), format = "fg", digits = 2)))
stopifnot(nrow(dm) == 2 * length(GA))
pd <- ggplot(dm, aes(BETA_STD, g, alpha = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = BETA_STD - 1.96 * se_std, xmax = BETA_STD + 1.96 * se_std), width = 0,
                orientation = "y", linewidth = 0.5, colour = "#b2182b") +
  geom_point(size = 1.6, colour = "#b2182b") +
  geom_text(aes(label = lab), vjust = -0.7, size = GT - 0.2, colour = "#b2182b") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.05, 0.1)) +
  facet_grid(grp ~ v, scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "MAGMA gene-property \u03b2 (std., 95% CI)", y = NULL, title = "d  Gene-level genetic risk (MAGMA)") +
  base_theme + theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
                     strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = BASE - 1),
                     panel.spacing.x = unit(8, "pt"), panel.spacing.y = unit(3, "pt"))

# ------------------------- e: the two ends of AHBA C3 (bottom / top 20%) ----------
# groups: bottom 20% / middle 60% / top 20% of the C3 region score (c3_q20.csv,
# 26_gradient_scores.py); box values = group means (c3_q20_profile.tsv, 29_...)
GRPCOL <- c(`C3-low` = "#2166ac", `C3-mid` = "grey60", `C3-high` = "#b2182b")
GRPLAB <- c(`C3-low` = "bottom 20% of C3", `C3-mid` = "middle 60%", `C3-high` = "top 20% of C3")
Q <- read.csv(file.path(RES, "c3_q20.csv"))
PQ <- read.delim(file.path(RES, "c3_q20_profile.tsv"))
stopifnot(nrow(Q) == 137, all(c("C3-low", "C3-mid", "C3-high") %in% Q$tier))
bx <- range(poly$x); by <- range(poly$y); sc <- 4.6 / diff(bx)
fb <- poly |> mutate(key = tolower(label)) |>
  left_join(Q |> transmute(key = tolower(label), tier), by = "key") |>
  mutate(X = 4.95 + (x - bx[1]) * sc, Y = 1.25 + (y - by[1]) * sc)
stopifnot(length(unique(fb$key[!is.na(fb$tier)])) == 137)
cen <- fb |> filter(view == "lateral", tier %in% c("C3-low", "C3-high")) |> group_by(tier) |>
  summarise(X = mean(X), Y = mean(Y), .groups = "drop")
lo <- setNames(PQ$C3.low, PQ$feature); hi <- setNames(PQ$C3.high, PQ$feature)
thin_lo <- lo[["normative_thinning_um_yr"]]; thin_hi <- hi[["normative_thinning_um_yr"]]
p_thin <- PQ$p_spin_high_minus_low[PQ$feature == "normative_thinning_um_yr"]
FEATS <- c(`Neuro-Ex` = "excitatory neurons", `Neuro-In` = "inhibitory neurons", L2 = "layer 2", L3 = "layer 3",
           L5 = "layer 5", Oligo = "oligodendrocytes", Astro = "astrocytes", Micro = "microglia",
           L1 = "layer 1", WM = "white matter")
stopifnot(all(names(FEATS) %in% PQ$feature))
box <- function(x0, tier, vals, head) {
  yy <- 3.0 - 0.29 * (seq_along(FEATS) - 1); z0 <- x0 + 3.0; k <- 0.7
  v <- unname(vals[names(FEATS)])
  list(rect = data.frame(x0 = x0, x1 = x0 + 4.65, y0 = 0.1, y1 = 3.9, tier = tier),
       head = data.frame(x = x0 + 0.12, y = c(3.66, 3.36), tier = tier, face = c("bold", "plain"),
                         lab = c(head, sprintf("thinning %.1f \u00b5m/yr \u00b7 CT %.2f mm",
                                               vals[["normative_thinning_um_yr"]], vals[["CT"]]))),
       lab  = data.frame(x = x0 + 0.12, y = yy, lab = unname(FEATS)),
       bar  = data.frame(xmin = pmin(z0, z0 + k * v), xmax = pmax(z0, z0 + k * v), y = yy, v = v, tx = x0 + 4.55),
       axis = data.frame(x = z0, y0 = 0.2, y1 = 3.15))
}
nlo <- sum(Q$tier == "C3-low"); nhi <- sum(Q$tier == "C3-high")
BL <- box(0.0, "C3-low", lo, sprintf("bottom 20%% of C3 (%d parcels)", nlo))
BR <- box(9.9, "C3-high", hi, sprintf("top 20%% of C3 (%d parcels)", nhi))
BB <- lapply(names(BL), function(k) bind_rows(BL[[k]], BR[[k]])); names(BB) <- names(BL)
arr <- data.frame(x = c(4.65, 9.9), y = c(2.3, 2.3), tier = c("C3-low", "C3-high"),
                  xend = cen$X[match(c("C3-low", "C3-high"), cen$tier)], yend = cen$Y[match(c("C3-low", "C3-high"), cen$tier)])
pe1 <- ggplot() +
  geom_polygon(data = fb, aes(X, Y, group = interaction(view, label, group, subgroup), fill = tier),
               colour = "grey35", linewidth = 0.05) +
  scale_fill_manual(values = GRPCOL, na.value = "grey88", guide = "none") +
  geom_rect(data = BB$rect, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, colour = tier), fill = NA, linewidth = 0.6) +
  geom_curve(data = arr, aes(x, y, xend = xend, yend = yend, colour = tier), curvature = -0.2,
             arrow = arrow(length = unit(5, "pt"), type = "closed"), linewidth = 0.6) +
  scale_colour_manual(values = GRPCOL, guide = "none") +
  geom_text(data = BB$head, aes(x, y, label = lab, fontface = face, colour = tier), hjust = 0, size = GT + 0.4) +
  geom_text(data = BB$lab, aes(x, y, label = lab), hjust = 0, size = GT, colour = "grey15") +
  geom_segment(data = BB$axis, aes(x = x, xend = x, y = y0, yend = y1), colour = "grey55", linewidth = 0.3) +
  geom_rect(data = BB$bar, aes(xmin = xmin, xmax = xmax, ymin = y - 0.1, ymax = y + 0.1),
            fill = ifelse(BB$bar$v > 0, "#d6604d", "#4393c3")) +
  geom_text(data = BB$bar, aes(tx, y, label = sprintf("%+.2f", v)), hjust = 1, size = GT - 0.3, colour = "grey25") +
  annotate("text", x = 7.25, y = 0.55, size = GT, colour = "grey30", lineheight = 0.95,
           label = "bars: mean z-expression of marker genes\n(descriptive: C3 is built from the same expression)") +
  coord_fixed(xlim = c(-0.05, 14.6), ylim = c(0.05, 3.95), expand = FALSE, clip = "off") +
  theme_void(base_size = BASE) + theme(plot.margin = margin(2, 2, 6, 2))

# e2: normative thinning of every parcel, by C3 group
Qp <- Q |> mutate(tier = factor(tier, names(GRPLAB)))
qm <- Qp |> group_by(tier) |> summarise(m = mean(normative_thinning_um_yr), s = sd(normative_thinning_um_yr) / sqrt(n()),
                                        n = n(), .groups = "drop")
stopifnot(abs(qm$m[qm$tier == "C3-low"] - thin_lo) < 0.01, abs(qm$m[qm$tier == "C3-high"] - thin_hi) < 0.01)   # table is 4 s.f.
set.seed(1)
pe2 <- ggplot(Qp, aes(tier, normative_thinning_um_yr, colour = tier)) +
  geom_jitter(width = 0.18, height = 0, size = 1.1, alpha = 0.6, stroke = 0) +
  geom_errorbar(data = qm, aes(x = tier, ymin = m - 1.96 * s, ymax = m + 1.96 * s), inherit.aes = FALSE,
                width = 0.25, linewidth = 0.7, colour = "grey10") +
  geom_point(data = qm, aes(tier, m), inherit.aes = FALSE, size = 2.2, colour = "grey10") +
  annotate("text", x = 2, y = max(Qp$normative_thinning_um_yr) * 1.06, size = GT + 0.2, fontface = "bold",
           label = sprintf("top \u2212 bottom = %+.1f \u00b5m/yr (%.0f%% slower at the bottom)\np_spin = %.3f",
                           thin_hi - thin_lo, 100 * (1 - thin_lo / thin_hi), p_thin)) +
  scale_colour_manual(values = GRPCOL, guide = "none") +
  scale_x_discrete(labels = sub(" of C3", "\nof C3", GRPLAB)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = NULL, y = "normative thinning (\u00b5m/yr)", subtitle = "Thinning rate of each parcel") +
  base_theme + theme(plot.subtitle = element_text(size = BASE, face = "bold"))

# e3: later symptoms vs each child's thinning in each C3 group
CG <- read.delim(file.path(RES, "gradient_scores_c3q20.tsv"))
OUTG <- c(totprob = "Total\nproblems", pfactor = "p-factor")
g3 <- CG |> filter(outcome %in% names(OUTG), spec %in% c("low alone", "mid alone", "high alone")) |>
  mutate(tier = factor(term, c("C3-high", "C3-mid", "C3-low")), o = factor(OUTG[outcome], rev(OUTG)),
         s = se * y_sd_beta / beta, lo = y_sd_beta - 1.96 * s, hi = y_sd_beta + 1.96 * s, sig = p < 0.05)
stopifnot(nrow(g3) == 3 * length(OUTG), !any(is.na(g3$tier)))
jt <- CG |> filter(outcome %in% names(OUTG), spec == "low + high | CT, QC", term == "low - high") |>
  mutate(o = factor(OUTG[outcome], rev(OUTG)), sig = p < 0.05,
         txt = sprintf("bottom \u2212 top\n(joint model)\n%+.3f, %s", y_sd_beta, sub("p = 0", "p = ", pf(p))))
stopifnot(nrow(jt) == length(OUTG))
XR <- max(g3$hi) * 1.1
pe3 <- ggplot(g3, aes(y_sd_beta, o, colour = tier, alpha = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.65, position = position_dodge(0.6)) +
  geom_point(size = 2.1, position = position_dodge(0.6)) +
  geom_text(data = jt, aes(x = XR, y = o, label = txt, alpha = sig), inherit.aes = FALSE, hjust = 0, size = GT, lineheight = 0.9) +
  scale_colour_manual(values = GRPCOL, breaks = names(GRPLAB), labels = unname(GRPLAB), name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.02, 0.04)) +
  coord_cartesian(xlim = c(min(g3$lo), XR * 1.75), clip = "off") +
  labs(x = "\u03b2 per SD of the child's thinning rate in that group\n(SD of CBCL at 15\u201317 | baseline CBCL)", y = NULL,
       subtitle = "Later symptoms vs thinning in each group") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  base_theme + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.subtitle = element_text(size = BASE, face = "bold"),
                     legend.position = "bottom", legend.key.size = unit(9, "pt"), legend.margin = margin(0, 0, 0, 0))

spec_ns <- all(jt$p >= 0.05)
pE <- (pe1 / (pe2 | pe3)) + plot_layout(heights = c(0.85, 1)) +
  plot_annotation(title = sprintf("e  The two ends of AHBA C3: the bottom 20%% thins %.0f%% slower, but symptom-linked thinning %s",
                                  100 * (1 - thin_lo / thin_hi),
                                  if (spec_ns) "is not specific to either end" else "differs between the ends"),
                  theme = theme(plot.title = element_text(size = BASE, face = "bold")))
pE <- wrap_elements(full = pE)

# -------------------------------------------- claims for the figure subtitle ----
e3 <- dm |> filter(v == "AHBA C3", sig, BETA_STD > 0) |> pull(gene_analysis)
SH <- c(SCZ25_META = "SCZ", MDD_div = "MDD", ASD = "ASD", ALZ = "Alzheimer's", EA = "EA", INT = "intelligence")
sub_txt <- sprintf(paste0(
  "a\u2013c: PLS of AHBA expression (hcp_3d_ds5: %d parcels \u00d7 %s genes) against ABCD thickness and thinning:\nPLS1 = AHBA C1 tracks baseline thickness, PLS2 = AHBA C3 tracks thinning rate (a\u2013c). C3 genes carry %s signal (d).\n",
  "e: the bottom 20%% of C3 is glia- and white-matter-rich and thins %.0f%% slower than the neuron- and upper-layer-rich top 20%% (p_spin = %.3f).\n%s"),
  cm("PLS1", "n_parcels"), format(cm("PLS1", "n_genes"), big.mark = ","), paste(SH[e3], collapse = ", "),
  100 * (1 - thin_lo / thin_hi), p_thin,
  if (spec_ns) "Children whose symptoms rise thin faster at both ends alike: the symptom link is not specific to either end of C3."
  else "The symptom-linked thinning differs between the two ends of C3.")

# ------------------------------------------------------------- assemble -------
row1 <- (pa | pb) + plot_layout(widths = c(2.55, 1))
left <- (pc / pd) + plot_layout(heights = c(1.15, 1))
row2 <- (wrap_elements(full = left) | pE) + plot_layout(widths = c(1.05, 1.7))
fig <- (wrap_elements(full = row1) / row2) + plot_layout(heights = c(0.62, 1.2))
fig <- fig + plot_annotation(
    title = "AHBA C3 separates fast-thinning, neuron-rich cortex from slow-thinning, glia-rich cortex, and carries psychiatric gene-level risk",
    subtitle = sub_txt,
    theme = theme(plot.title = element_text(size = BASE + 3, face = "bold"),
                  plot.subtitle = element_text(size = BASE - 0.5, colour = "grey25", lineheight = 1.15, margin = margin(b = 4))))
OUT <- file.path(FIG, "fig_hcp_c3_mechanism_3d_ds5.png")
ggsave(OUT, fig, width = 12.5, height = 12.5, dpi = 300, bg = "white")
cat("wrote", OUT, "\n")
