#!/usr/bin/env Rscript
# fig5s_c3_mechanism.R -- simplified mechanism slide, hcp_3d_ds5 fit only:
#   a  CT / PLS1 / C1 and dCT / PLS2 / C3 maps
#   b  the two strong pairs (PLS1 vs CT, PLS2 vs dCT)
#   c  Spearman matrix of the six region maps, spin p   (29_c3_mechanism_inputs.py)
#   d  cell-class / layer enrichment of AHBA C1 and C3 gene weights
#   e  MAGMA gene-property for C1 and C3 (SCZ, MDD, ASD, ALZ, EA, intelligence)
#   f  schematic: thirds of the AHBA C3 region score, their normative thinning
#      and marker-gene expression (c3_tiers.csv, c3_tier_profile.tsv)
#   g  later CBCL vs thinning in each C3 third (gradient_scores_c3tiers.tsv, 26)
# Reads only saved tables and fits nothing. Faded (alpha 0.3) = not significant
# (c: spin p; d: BH q; e, g: p). Writes figures/fig_hcp_c3_mechanism_3d_ds5.png
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
MV <- c(CT = "CT", dCT = "dCT", PLS1 = "PLS1", PLS2 = "PLS2", C1 = "AHBA C1", C3 = "AHBA C3")
stopifnot(nrow(SM) == choose(length(MV), 2), all(c(SM$x, SM$y) %in% names(MV)))
cmx <- SM |> mutate(ix = match(x, names(MV)), iy = match(y, names(MV))) |>
  transmute(r = MV[pmax(ix, iy)], c = MV[pmin(ix, iy)], rho, p_spin) |>      # lower triangle
  mutate(r = factor(r, rev(MV)), c = factor(c, MV), sig = p_spin < 0.05, txt = sprintf("%.2f", rho))
pc <- ggplot(cmx, aes(c, r)) +
  geom_tile(aes(fill = rho, alpha = sig), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt, alpha = sig, fontface = ifelse(sig, "bold", "plain")), size = GT + 0.3) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), breaks = c(-1, 0, 1), name = "Spearman\n\u03c1") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) + coord_fixed() +
  labs(x = NULL, y = NULL, title = "c  Regional maps: correlation matrix",
       subtitle = "bold = p_spin < 0.05; faded = n.s.") +
  base_theme + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                     axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(),
                     plot.subtitle = element_text(size = BASE - 1.5, colour = "grey35"),
                     legend.key.height = unit(14, "pt"), legend.key.width = unit(5, "pt"))

# --------------------------------------- d: cell classes and layers, C1 / C3 only ----
SET_ORD <- c("Neuro-Ex", "Neuro-In", "Astro", "Micro", "Oligo", "OPC", "Endo", "Per",
             "L1", "L2", "L3", "L4", "L5", "L6", "WM")
VORD <- c("AHBA C1", "AHBA C3")
ZMAX <- 20
e <- ST |> filter(set != "Neuro", vector %in% c("C1", "C3")) |>
  mutate(vector = factor(recode(vector, C1 = "AHBA C1", C3 = "AHBA C3"), VORD),
         set = factor(set, rev(SET_ORD)), kind = factor(kind, c("cell class", "layer"), c("cell classes", "layers")),
         zc = pmax(pmin(z, ZMAX), -ZMAX), az = abs(zc), sig = q_bh < 0.05)
stopifnot(!any(is.na(e$set)), nrow(e) == length(SET_ORD) * length(VORD))
pd <- ggplot(e, aes(vector, set)) +
  geom_point(aes(size = az, fill = zc, alpha = sig), shape = 22, colour = "grey30", stroke = 0.15) +
  scale_size_area(max_size = 6.5, limits = c(0, ZMAX), breaks = c(5, 10, 20), labels = c("5", "10", "\u226520"), name = "|z|") +
  scale_fill_distiller(palette = "RdBu", limits = c(-ZMAX, ZMAX), breaks = c(-ZMAX, 0, ZMAX),
                       labels = c("\u2212", "0", "+"), name = "direction") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), labels = c(`TRUE` = "q < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_x_discrete(position = "top") +
  facet_grid(kind ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, title = "d  Cell classes and layers") +
  guides(size = guide_legend(order = 1, override.aes = list(fill = "grey60")),
         fill = guide_colourbar(order = 2, barwidth = unit(4, "pt"), barheight = unit(20, "pt")),
         alpha = guide_legend(order = 3, override.aes = list(size = 3, fill = "#b2182b"))) +
  base_theme + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                     axis.text.x.top = element_text(size = BASE - 1),
                     strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = BASE - 1),
                     panel.grid.major = element_line(colour = "grey94", linewidth = 0.25),
                     legend.key.size = unit(7, "pt"), legend.spacing.y = unit(1, "pt"))

# ---------------------------------------------- e: MAGMA, C1 / C3 only ----------
VEC <- c(AHBA_C1 = "AHBA C1", AHBA_C3 = "AHBA C3")
GA <- c(SCZ25_META = "Schizophrenia", MDD_div = "Depression", ASD = "Autism",
        ALZ = "Alzheimer's disease", EA = "Educational attainment", INT = "Intelligence")
GGRP <- c(SCZ25_META = "psychiatric", MDD_div = "psychiatric", ASD = "psychiatric",
          ALZ = "neuro-\ndegenerative", EA = "cognitive", INT = "cognitive")
stopifnot(all(names(GA) %in% MG$gene_analysis))
dm <- MG |> filter(VARIABLE %in% names(VEC), gene_analysis %in% names(GA)) |>
  mutate(v = factor(VEC[VARIABLE], VEC), g = factor(GA[gene_analysis], rev(GA)),
         grp = factor(GGRP[gene_analysis], unique(GGRP)), sig = P < 0.05,
         lab = ifelse(P < 0.001, sprintf("%.0e", P), formatC(signif(P, 2), format = "fg", digits = 2)))
stopifnot(nrow(dm) == 2 * length(GA))
pe <- ggplot(dm, aes(BETA_STD, g, alpha = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = BETA_STD - 1.96 * se_std, xmax = BETA_STD + 1.96 * se_std), width = 0,
                orientation = "y", linewidth = 0.45, colour = "#b2182b") +
  geom_point(size = 1.3, colour = "#b2182b") +
  geom_text(aes(label = lab), vjust = -0.65, size = GT - 0.15, colour = "#b2182b") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.05, 0.1)) +
  facet_grid(grp ~ v, scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "MAGMA gene-property \u03b2 (std., 95% CI)", y = NULL, title = "e  Gene-level genetic risk (MAGMA)") +
  base_theme + theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
                     strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = BASE - 1),
                     panel.spacing.x = unit(6, "pt"), panel.spacing.y = unit(3, "pt"))

# ---------------------------------- f: schematic -- what the C3 thirds are ---------
# thirds of the AHBA C3 region score (c3_tiers.csv, 26_gradient_scores.py); the
# box values are tier means (c3_tier_profile.tsv, 29_c3_mechanism_inputs.py)
TIERCOL <- c(`C3-low` = "#2166ac", `C3-mid` = "grey60", `C3-high` = "#b2182b")
C3T <- read.csv(file.path(RES, "c3_tiers.csv"))
PR  <- read.delim(file.path(RES, "c3_tier_profile.tsv"))
stopifnot(nrow(C3T) == 137)
bx <- range(poly$x); by <- range(poly$y); sc <- 4.6 / diff(bx)
fb <- poly |> mutate(key = tolower(label)) |>
  left_join(C3T |> transmute(key = tolower(label), tier), by = "key") |>
  mutate(X = 4.9 + (x - bx[1]) * sc, Y = 1.35 + (y - by[1]) * sc)
stopifnot(length(unique(fb$key[!is.na(fb$tier)])) == 137)
cen <- fb |> filter(view == "lateral", !is.na(tier)) |> group_by(tier) |>
  summarise(X = mean(X), Y = mean(Y), .groups = "drop")
th <- setNames(PR$C3.low, PR$feature); th2 <- setNames(PR$C3.high, PR$feature)
thin_lo <- th[["normative_thinning_um_yr"]]; thin_hi <- th2[["normative_thinning_um_yr"]]
ct_lo <- th[["CT"]]; ct_hi <- th2[["CT"]]
c3d <- SM |> filter(x == "dCT", y == "C3"); stopifnot(nrow(c3d) == 1)
FEATS <- c(`Neuro-Ex` = "excitatory neurons", `Neuro-In` = "inhibitory neurons", L2 = "layer 2", L3 = "layer 3",
           L5 = "layer 5", Oligo = "oligodendrocytes", Astro = "astrocytes", Micro = "microglia",
           L1 = "layer 1", WM = "white matter")
stopifnot(all(names(FEATS) %in% PR$feature))
box <- function(x0, tier, vals, thin, ct, head) {
  yy <- 3.05 - 0.3 * (seq_along(FEATS) - 1); z0 <- x0 + 2.95; k <- 0.8     # 0.8 units per SD
  v <- unname(vals[names(FEATS)])
  list(rect = data.frame(x0 = x0, x1 = x0 + 4.6, y0 = 0.02, y1 = 3.95, tier = tier),
       head = data.frame(x = x0 + 0.12, y = c(3.72, 3.45),
                         lab = c(head, sprintf("thinning %.1f \u00b5m/yr \u00b7 CT %.2f mm", thin, ct)),
                         face = c("bold", "plain"), tier = tier),
       lab  = data.frame(x = x0 + 0.12, y = yy, lab = unname(FEATS)),
       bar  = data.frame(xmin = pmin(z0, z0 + k * v), xmax = pmax(z0, z0 + k * v), y = yy, v = v,
                         tx = x0 + 4.5, hj = 1),
       axis = data.frame(x = z0, y0 = 0.12, y1 = 3.2))
}
nlo <- sum(C3T$tier == "C3-low"); nhi <- sum(C3T$tier == "C3-high")
L <- box(0.0, "C3-low", th, thin_lo, ct_lo, sprintf("bottom third of C3 (%d parcels)", nlo))
R <- box(9.9, "C3-high", th2, thin_hi, ct_hi, sprintf("top third of C3 (%d parcels)", nhi))
BB <- lapply(names(L), function(k) bind_rows(L[[k]], R[[k]])); names(BB) <- names(L)
arr <- data.frame(x = c(4.6, 9.9), y = c(2.3, 2.3), xend = cen$X[match(c("C3-low", "C3-high"), cen$tier)],
                  yend = cen$Y[match(c("C3-low", "C3-high"), cen$tier)], tier = c("C3-low", "C3-high"))
slower <- 100 * (1 - thin_lo / thin_hi)
pF <- ggplot() +
  geom_polygon(data = fb, aes(X, Y, group = interaction(view, label, group, subgroup), fill = tier),
               colour = "grey35", linewidth = 0.05) +
  scale_fill_manual(values = TIERCOL, na.value = "grey88", guide = "none") +
  geom_rect(data = BB$rect, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, colour = tier), fill = NA, linewidth = 0.5) +
  geom_curve(data = arr, aes(x, y, xend = xend, yend = yend, colour = tier), curvature = -0.2,
             arrow = arrow(length = unit(4, "pt"), type = "closed"), linewidth = 0.5) +
  scale_colour_manual(values = TIERCOL, guide = "none") +
  geom_text(data = BB$head, aes(x, y, label = lab, fontface = face, colour = tier), hjust = 0, size = GT + 0.3) +
  geom_text(data = BB$lab, aes(x, y, label = lab), hjust = 0, size = GT, colour = "grey15") +
  geom_segment(data = BB$axis, aes(x = x, xend = x, y = y0, yend = y1), colour = "grey55", linewidth = 0.3) +
  geom_rect(data = BB$bar, aes(xmin = xmin, xmax = xmax, ymin = y - 0.11, ymax = y + 0.11),
            fill = ifelse(BB$bar$v > 0, "#d6604d", "#4393c3")) +
  geom_text(data = BB$bar, aes(tx, y, label = sprintf("%+.2f", v), hjust = hj), size = GT - 0.3, colour = "grey25") +
  annotate("text", x = 7.25, y = 0.75, size = GT + 0.3, lineheight = 0.95,
           label = sprintf("bottom third thins %.0f%% slower than the top third\n(C3 vs dCT: \u03c1 = %.2f, p_spin = %.3f)\nat a similar baseline thickness\n\nbars: mean z-expression of marker genes",
                           slower, c3d$rho, c3d$p_spin)) +
  coord_fixed(xlim = c(-0.05, 14.55), ylim = c(0, 4.0), expand = FALSE, clip = "off") +
  labs(title = "f  What AHBA C3 separates: slow-thinning, glia-rich cortex vs fast-thinning, neuron-rich cortex") +
  theme_void(base_size = BASE) + theme(plot.title = element_text(size = BASE, face = "bold", margin = margin(b = 2)),
                                       plot.margin = margin(2, 4, 2, 4))

# ----------------------------- g: symptoms vs thinning in each C3 third ----------
CG <- read.delim(file.path(RES, "gradient_scores_c3tiers.tsv"))
OUTG <- c(totprob = "Total problems", pfactor = "p-factor")
g3 <- CG |> filter(outcome %in% names(OUTG), spec %in% c("low alone", "mid alone", "high alone")) |>
  mutate(tier = factor(term, c("C3-high", "C3-mid", "C3-low")), o = factor(OUTG[outcome], rev(OUTG)),
         s = se * y_sd_beta / beta, lo = y_sd_beta - 1.96 * s, hi = y_sd_beta + 1.96 * s, sig = p < 0.05)
stopifnot(nrow(g3) == 3 * length(OUTG), !any(is.na(g3$tier)))
jt <- CG |> filter(outcome %in% names(OUTG), spec == "low + high | CT, QC", term == "low - high") |>
  mutate(o = factor(OUTG[outcome], rev(OUTG)), sig = p < 0.05,
         txt = sprintf("low \u2212 high (joint)\n%+.3f, %s", y_sd_beta, sub("p = 0", "p = ", pf(p))))
stopifnot(nrow(jt) == length(OUTG))
XR <- max(g3$hi) * 1.08
pG <- ggplot(g3, aes(y_sd_beta, o, colour = tier, alpha = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.55, position = position_dodge(0.6)) +
  geom_point(size = 1.7, position = position_dodge(0.6)) +
  geom_text(data = jt, aes(x = XR, y = o, label = txt, alpha = sig), inherit.aes = FALSE, hjust = 0, size = GT, lineheight = 0.9) +
  scale_colour_manual(values = TIERCOL, breaks = c("C3-low", "C3-mid", "C3-high"),
                      labels = c("bottom third of C3", "middle third", "top third of C3"), name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  coord_cartesian(xlim = c(min(g3$lo), XR * 1.55), clip = "off") +
  labs(x = "\u03b2 per SD of thinning rate in that third\n(SD of CBCL at 15\u201317 | baseline CBCL)", y = NULL,
       title = "g  Is symptom-linked thinning specific to one pole of C3?") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(alpha = 1))) +
  base_theme + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                     legend.position = "bottom", legend.key.size = unit(8, "pt"), legend.margin = margin(0, 0, 0, 0))

# -------------------------------------------- claims for the figure subtitle ----
e3 <- dm |> filter(v == "AHBA C3", sig, BETA_STD > 0) |> pull(gene_analysis)
spec_ns <- all(jt$p >= 0.05)
SH <- c(SCZ25_META = "SCZ", MDD_div = "MDD", ASD = "ASD", ALZ = "Alzheimer's", EA = "EA", INT = "intelligence")
sub_txt <- sprintf(paste0(
  "a\u2013c: PLS of AHBA expression (hcp_3d_ds5: %d parcels \u00d7 %s genes) against ABCD thickness and thinning: PLS1 = AHBA C1 tracks baseline thickness, PLS2 = AHBA C3 tracks thinning rate.\n",
  "d\u2013f: C3-high cortex is neuron- and upper-layer-rich and thins fastest; C3-low cortex is glia- and white-matter-rich and thins %.0f%% slower. C3 genes carry %s signal (e).\n",
  "g: %s"),
  cm("PLS1", "n_parcels"), format(cm("PLS1", "n_genes"), big.mark = ","), slower, paste(SH[e3], collapse = ", "),
  if (spec_ns) "children whose symptoms rise thin faster in both poles of C3 alike; the low \u2212 high contrast is null for total problems and p-factor, so the symptom link is not specific to either pole."
  else "the symptom-linked thinning differs between the C3 poles (low \u2212 high contrast p < 0.05).")

# ------------------------------------------------------------- assemble -------
row1 <- (pa | pb) + plot_layout(widths = c(2.6, 1))
row2 <- (pc | pd | pe) + plot_layout(widths = c(1.15, 0.75, 1.2))
row3 <- (pF | pG) + plot_layout(widths = c(2.1, 1))
fig <- (wrap_elements(full = row1) / wrap_elements(full = row2) / wrap_elements(full = row3)) +
  plot_layout(heights = c(1.0, 1.05, 1.15))
fig <- fig + plot_annotation(
    title = "AHBA C3 separates fast-thinning, neuron-rich cortex from slow-thinning, glia-rich cortex, and carries psychiatric gene-level risk",
    subtitle = sub_txt,
    theme = theme(plot.title = element_text(size = BASE + 2, face = "bold"),
                  plot.subtitle = element_text(size = BASE - 0.5, colour = "grey25", lineheight = 1.15, margin = margin(b = 4))))
OUT <- file.path(FIG, "fig_hcp_c3_mechanism_3d_ds5.png")
ggsave(OUT, fig, width = 12, height = 9.6, dpi = 300, bg = "white")
cat("wrote", OUT, "\n")
