#!/usr/bin/env Rscript
# fig1.R -- Figure 1 DRAFT: adolescent cortical thinning and the polygenic
# scores that predict it.  HCP-MMP1.0 throughout, ABCD release 7.0.
#
# Usage (repo root):  LC_ALL=en_US.UTF-8 Rscript genetic_analysis/fig1.R
#   (the UTF-8 locale is required: µ, →, β and · are in the source strings)
#
# Reads only committed tables in genetic_analysis/fig1_inputs/ (plus the lh
# Glasser polygons in ahba_pls/data/hcp_polygons.csv) and fits nothing; every
# printed statistic is read from the table plotted, via sprintf.
#   hcp70_group_maps.csv              a   per-parcel CT at the age centre and ΔCT
#   hcp70_age_by_visit.csv,
#   hcp70_scans_per_child.csv         b   design
#   hcp70_scans.csv (gitignored, individual-level; falls back to
#   hcp70_scans_hist2d.csv), hcp70_scan_summary.csv       c   trajectories
#   hcp70_regional_slope_reliability.csv,
#   hcp70_global_slope_splithalf.csv  d   reliability
#   hcp70_prs_key_arms.tsv            e,f PRS (single-LMM columns)
#   hcp70_magma_locus_sets.tsv (+ work/results_70tab_hcp/magma_pooled/
#   (MAGMA gene-set panel removed; gene-set results go to SI)
# Writes docs/figures/fig1.png and docs/figures/fig1_caption.md.

stopifnot("run with LC_ALL=en_US.UTF-8" = l10n_info()$`UTF-8`)
suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
HERE <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
REPO <- dirname(HERE)
IN <- file.path(HERE, "fig1_inputs")
OUT <- file.path(REPO, "docs/figures/fig1.png")
CAP <- file.path(REPO, "docs/figures/fig1_caption.md")
POLY <- file.path(REPO, "ahba_pls/data/hcp_polygons.csv")

# ---------------------------------------------------------------- style ----
BASE <- 7
SLOPE_C <- "#B2182B"; BASE_C <- "#2171b5"   # ΔCT red, CT blue (as panels a, e)
METHODS <- c(CT = "C+T", PRSCS = "PRS-CS", SBayesRC = "SBayesRC")  # SBayesR in SI
METHOD_C <- c("C+T" = "grey25", "PRS-CS" = "#3A9AD9", "SBayesRC" = "#C0569E")
VISIT_C <- c(v0 = "grey82", v2 = "grey62", v4 = "grey42", v6 = "grey20")
VISIT_L <- c(v0 = "baseline", v2 = "2-year", v4 = "4-year", v6 = "6-year")

th <- theme_bw(base_size = BASE) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        axis.text = element_text(size = BASE - 1, colour = "grey20"),
        axis.title = element_text(size = BASE - 0.5),
        plot.title = element_text(size = BASE, face = "bold", hjust = 0,
                                  margin = margin(b = 1)),
        plot.subtitle = element_text(size = BASE - 1, colour = "grey30",
                                     hjust = 0, margin = margin(b = 3)),
        plot.title.position = "plot",
        legend.text = element_text(size = BASE - 1),
        legend.title = element_blank(),
        legend.key.size = unit(7, "pt"),
        legend.background = element_blank(),
        plot.margin = margin(3, 5, 3, 3))

stars <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
fmt_p <- function(p) ifelse(p < 1e-3, sprintf("%.0e", p),
                     ifelse(p < 0.01, sprintf("%.3f", p), sprintf("%.2f", p)))

# ------------------------------------------------------------ a: maps ------
maps <- read.csv(file.path(IN, "hcp70_group_maps.csv"))
poly <- read.csv(POLY)
# lateral and medial side by side, aspect preserved (one scale for both views)
sx <- max(tapply(poly$x, poly$view, function(v) diff(range(v))))
poly <- poly |>
  group_by(view) |>
  mutate(x = (x - min(x)) / sx + ifelse(view == "medial", 1.06, 0),
         y = (y - min(y)) / sx) |>
  ungroup() |>
  mutate(grp = interaction(label, view, group, subgroup, drop = TRUE))
n_missing <- length(setdiff(unique(poly$label), maps$label))

# both maps and their colourbars in ONE fixed-aspect ggplot, colours
# precomputed with the fig5_hcp_summary palettes (white->blue CT, red->white dCT)
CT_COLS  <- c("white", "#c6dbef", "#6baed6", "#2171b5", "#08306b")
DCT_COLS <- c("#67000d", "#cb181d", "#fb6a4a", "#fcbba1", "white")
lo <- unname(floor(quantile(maps$ct_at_centre, 0.02) * 10) / 10)
hi <- unname(ceiling(quantile(maps$ct_at_centre, 0.98) * 10) / 10)
um <- maps$slope_mm_per_yr * 1000
vmin <- unname(floor(quantile(um, 0.02) / 5) * 5)
pal <- function(cols, v, a, b) {
  f <- scales::gradient_n_pal(cols); f(scales::rescale(pmin(pmax(v, a), b), from = c(a, b)))
}
H <- max(poly$y); W <- max(poly$x); GAP <- 0.30          # map height/width, data units
row_off <- c(base = 0, thin = -(H + GAP))
map_df <- bind_rows(
  poly |> left_join(transmute(maps, label, col = pal(CT_COLS, ct_at_centre, lo, hi)), by = "label") |>
    mutate(y = y + row_off["base"], grp = paste("b", grp)),
  poly |> left_join(transmute(maps, label, col = pal(DCT_COLS, slope_mm_per_yr * 1000, vmin, 0)),
                    by = "label") |>
    mutate(y = y + row_off["thin"], grp = paste("t", grp))) |>
  mutate(col = ifelse(is.na(col), "grey85", col))
# vertical colourbars right of each map
bar <- function(cols, a, b, off, lab, fmt) {
  n <- 60; v <- seq(a, b, length.out = n); hgt <- H * 0.8; y0 <- off + H * 0.1
  list(tiles = data.frame(x = W + 0.10, y = y0 + (seq_len(n) - 0.5) * hgt / n,
                          fill = pal(cols, v, a, b), h = hgt / n),
       text = data.frame(x = W + 0.17, y = c(y0, y0 + hgt), lab = sprintf(fmt, c(a, b))),
       title = data.frame(x = W + 0.07, y = y0 + hgt + 0.08, lab = lab))
}
b1 <- bar(CT_COLS, lo, hi, row_off["base"], "mm", "%.1f")
b2 <- bar(DCT_COLS, vmin, 0, row_off["thin"], "µm/yr", "%.0f")
tiles <- rbind(b1$tiles, b2$tiles); btxt <- rbind(b1$text, b2$text); bttl <- rbind(b1$title, b2$title)
PA_ASPECT <- (W + 0.42) / (H + 0.02 - row_off["thin"])    # width / height of the map block
ROW1_H_IN <- 1.75   # effective map-block height, calibrated on the render
# the "a   " tag sits inside panel a's title; shift the second title by its width
tag_in <- grid::convertWidth(grid::grobWidth(grid::textGrob("a   ",
            gp = grid::gpar(fontsize = BASE, fontface = "bold"))), "inches", valueOnly = TRUE)
TAG_DX <- 1.16 * tag_in * (H + 0.02 - row_off["thin"]) / ROW1_H_IN   # 1.16: calibrated on the render
pa <- ggplot() +
  geom_polygon(data = map_df, aes(x, y, group = grp, fill = col), colour = "white",
               linewidth = 0.06) +
  geom_tile(data = tiles, aes(x, y, fill = fill, height = h), width = 0.06) +
  geom_text(data = btxt, aes(x, y, label = lab), hjust = 0, size = (BASE - 1.5) / .pt) +
  geom_text(data = bttl, aes(x, y, label = lab), hjust = 0, vjust = 0, size = (BASE - 1.5) / .pt) +
  annotate("text", x = TAG_DX, y = row_off["thin"] + H + 0.06, label = "Thinning rate (ΔCT)",
           hjust = 0, vjust = 0, fontface = "bold", size = BASE / .pt) +
  scale_fill_identity() +
  coord_fixed(xlim = c(0, W + 0.42), ylim = c(row_off["thin"], H + 0.02), expand = FALSE,
              clip = "off") +
  guides(x = "none", y = "none") +
  labs(title = "a   Cortical thickness (CT)") +
  theme_void(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0, margin = margin(b = 4)),
        plot.title.position = "plot", plot.margin = margin(3, 3, 3, 3))

# ------------------------------------------------------------ b: design ----
age <- read.csv(file.path(IN, "hcp70_age_by_visit.csv"))
names(age)[1] <- "visit"
per <- read.csv(file.path(IN, "hcp70_scans_per_child.csv"))
n_kids <- sum(per$n_children); n_scans <- sum(age$count)
share <- setNames(per$n_children / n_kids, per$n_scans)
scans_f <- file.path(IN, "hcp70_scans.csv")
if (file.exists(scans_f)) {
  sc <- read.csv(scans_f)
  pb <- ggplot(sc, aes(age, fill = visit)) +
    geom_histogram(binwidth = 0.25, boundary = 8, colour = NA)
  ymax <- max(table(cut(sc$age, seq(8, 19, 0.25)))) * 1.0
} else {
  xs <- seq(8, 18.5, 0.05)
  dd <- do.call(rbind, lapply(seq_len(nrow(age)), function(i) data.frame(
    visit = age$visit[i], age = xs,
    n = age$count[i] * 0.25 * dnorm(xs, age$mean[i], age$std[i]))))
  pb <- ggplot(dd, aes(age, n, fill = visit)) + geom_area(position = "identity")
  ymax <- max(dd$n)
}
pb <- pb +
  annotate("text", x = age$mean, y = ymax * c(1.12, 1.05, 1.12, 1.05),
           label = VISIT_L[age$visit], colour = c("grey45", VISIT_C[-1]),
           size = (BASE - 1.5) / .pt) +
  scale_fill_manual(values = VISIT_C, guide = "none") +
  scale_x_continuous(breaks = seq(8, 18, 2), limits = c(8, 18.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), labels = comma) +
  coord_cartesian(ylim = c(0, ymax * 1.2)) +
  labs(x = "age at scan (years)", y = "scans",
       title = sprintf("b   %s children, %s scans", comma(n_kids), comma(n_scans))) +
  th


# per-child phenotype distributions (fig1_prep_1lmm.R): the two traits of e-h
cd <- read.csv(file.path(IN, "hcp70_child_density.csv"))
cs <- read.csv(file.path(IN, "hcp70_child_summary.csv"))
csv_ <- function(v, f) cs[[f]][cs$var == v]
dens_panel <- function(v, colour, xlab, title, fmt) {
  d <- cd[cd$var == v, ]
  m <- csv_(v, "mean"); sdv <- csv_(v, "sd")
  ggplot(d, aes(x, density)) +
    geom_area(fill = colour, alpha = 0.35, colour = colour, linewidth = 0.4) +
    annotate("segment", x = m, xend = m, y = 0, yend = approx(d$x, d$density, m)$y,
             linewidth = 0.3, colour = colour) +
    annotate("text", x = max(d$x), y = max(d$density) * 1.62, label = sprintf(fmt, m, sdv),
             size = (BASE - 1.5) / .pt, colour = "grey20", hjust = 1, vjust = 1,
             lineheight = 0.9) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.65))) +
    labs(x = xlab, y = NULL, title = title) +
    th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
               axis.line.y = element_blank(), panel.grid.major.y = element_blank())
}
pb1 <- dens_panel("thickness_mm", BASE_C,
                  sprintf("CT at age %.1f (mm)", csv_("thickness_mm", "age_centre")),
                  "e   Whole-cortex traits", "mean %.2f\nSD %.3f mm")
pb2 <- dens_panel("slope_um_per_yr", SLOPE_C, "ΔCT (µm / year)", NULL,
                  "mean %.1f\nSD %.1f µm/yr")

pe <- (pb1 / pb2)

# ------------------------------------------------------ c: trajectories ----
summ <- read.csv(file.path(IN, "hcp70_scan_summary.csv"))
sv <- setNames(summ$value, summ$metric)
YL <- c(2.3, 3.0)
trend <- data.frame(age = c(8.3, 18.2)) |>
  mutate(ct = sv["ols_intercept_mm"] + sv["ols_slope_mm_per_yr"] * age)
HI_C <- c("#E08214", "#762A83", "#1B7837")               # highlighted children: not yellow, CT blue or ΔCT red
if (file.exists(scans_f)) {
  set.seed(7)
  kids <- sample(unique(sc$sid), 250)
  set.seed(25)                                   # highlighted-child draw (grey sample unchanged)
  # highlighted children, chosen on the fitted traits so they differ visibly in
  # BOTH CT and ΔCT: a thick fast thinner and a thin slow thinner (4 scans each)
  # and a near-average child (3 scans); smooth trajectories only
  tr <- read.csv(file.path(IN, "hcp70_child_traits.csv"))
  cen <- csv_("thickness_mm", "age_centre")
  zz <- function(v) (v - mean(v)) / sd(v)
  cand <- sc |> group_by(sid, n_visits) |> summarise(.groups = "drop") |>
    inner_join(mutate(tr, zct = zz(ct_mm), zdct = zz(dct_um_per_yr)), by = "sid") |>
    inner_join(sc |> inner_join(tr, by = "sid") |>
                 mutate(res = mean_ct - (ct_mm + dct_um_per_yr / 1000 * (age - cen))) |>
                 group_by(sid) |> summarise(rsd = sd(res)), by = "sid") |>
    filter(rsd < 0.03)
  pick <- function(n, f) {                 # x[sample.int()] -- sample(x, 1) on a length-1 x samples 1:x
    x <- filter(cand, n_visits == n, !!f)$sid
    stopifnot(length(x) > 0); x[sample.int(length(x), 1)]
  }
  hi3 <- c(pick(4, quote(zct > 1.0 & zct < 1.8 & zdct < -1.5 & zdct > -2.5)),
           pick(4, quote(zct < -1.0 & zct > -1.8 & zdct > 1.5 & zdct < 2.5)),
           pick(3, quote(abs(zct) < 0.3 & abs(zdct) < 0.3)))
  fit_hi <- tr |> filter(sid %in% hi3) |> mutate(child = factor(match(sid, hi3))) |>
    inner_join(sc |> group_by(sid) |> summarise(a0 = min(age), a1 = max(age)), by = "sid") |>
    tidyr::pivot_longer(c(a0, a1), values_to = "age") |>
    mutate(ct = ct_mm + dct_um_per_yr / 1000 * (age - cen))
  lines_bg <- sc |> filter(sid %in% kids)
  lines_hi <- sc |> filter(sid %in% hi3) |> mutate(child = factor(match(sid, hi3)))
  pc <- ggplot() +
    geom_line(data = lines_bg, aes(age, mean_ct, group = sid, colour = "one child (250 shown)"),
              linewidth = 0.2, alpha = 0.7)
} else {
  h2 <- read.csv(file.path(IN, "hcp70_scans_hist2d.csv"))
  pc <- ggplot() + geom_tile(data = h2, aes(age_bin + 0.125, ct_bin + 0.01, alpha = n),
                             fill = "grey50") + scale_alpha(guide = "none")
}
pc <- pc +
  {if (file.exists(scans_f)) list(
    geom_line(data = fit_hi, aes(age, ct, group = sid), colour = HI_C[fit_hi$child],
              linewidth = 0.45),
    geom_point(data = lines_hi, aes(age, mean_ct), colour = HI_C[lines_hi$child], size = 1.1))} +
  scale_colour_manual(values = c("one child (250 shown)" = "grey72"),
                      breaks = "one child (250 shown)") +
  scale_x_continuous(breaks = seq(8, 18, 2)) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  coord_cartesian(xlim = c(8.5, 18)) +
  labs(x = "age (years)", y = "cortex-wide CT (mm)",
       title = "c   CT ~ age + sex + (age|child) + (1|site)") +
  th + theme(plot.title = element_text(size = BASE - 0.5, face = "bold"),legend.position = c(0.02, 0.01), legend.justification = c(0, 0),
             legend.key.width = unit(8, "pt"), legend.key.height = unit(6, "pt"),
             legend.text = element_text(size = BASE - 1.5),
             legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

# ------------------------------------------------------- d: reliability ----
rel <- read.csv(file.path(IN, "hcp70_regional_slope_reliability.csv"))
cw <- read.csv(file.path(IN, "hcp70_global_slope_reliability_1lmm.csv"))  # fig1_prep_1lmm.R
med <- rel |> group_by(n_visits) |> summarise(m = median(reliability, na.rm = TRUE))
pd_ <- ggplot() +
  geom_boxplot(data = rel, aes(factor(n_visits), reliability),
               outlier.shape = NA, width = 0.55, linewidth = 0.3,
               colour = "grey35", fill = "white", fatten = 0) +
  geom_segment(data = med, aes(x = as.numeric(factor(n_visits)) - 0.27,
                               xend = as.numeric(factor(n_visits)) + 0.27,
                               y = m, yend = m), colour = SLOPE_C, linewidth = 0.6) +
  geom_line(data = cw, aes(factor(n_visits), mean, group = 1), linewidth = 0.4) +
  geom_point(data = cw, aes(factor(n_visits), mean), shape = 15, size = 1.3) +
  annotate("text", x = 0.55, y = 0.66, label = "cortex-wide ΔCT\n(single LMM, mean)",
           hjust = 0, vjust = 1, size = (BASE - 1.5) / .pt, lineheight = 0.9) +
  annotate("text", x = 0.55, y = 0.50, label = "one parcel\n(358 parcels)",
           hjust = 0, vjust = 1, size = (BASE - 1.5) / .pt, colour = SLOPE_C, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.6, 0.2)) +
  labs(x = "scans per child", y = expression(Delta * "CT reliability  1 - v/" * tau^2),
       title = "d   ΔCT reliability") +
  th

# ------------------------------------------------------------ e, f: PRS ----
ROWS <- tribble(
  ~trait_arm,   ~label,                     ~arm,
  "SCZ25_EUR",  "schizophrenia · EUR",      "EUR",
  "SCZ25_META", "schizophrenia · pooled",   "pooled",
  "MDD_eur",    "depression · EUR",         "EUR",
  "MDD_pooled", "depression · pooled",      "pooled",
  "ALZ",        "Alzheimer's",              "EUR",
  "ALZ_noAPOE", "Alzheimer's, no APOE",     "EUR",
  "EA",         "education",                "EUR",
  "ASD",        "autism",                   "EUR") |>
  mutate(y = rev(seq_len(n())))
prs <- read.delim(file.path(IN, "hcp70_prs_key_arms.tsv")) |>
  filter(method %in% names(METHODS)) |>
  inner_join(ROWS, by = "trait_arm") |>
  mutate(method = factor(METHODS[method], levels = METHODS),
         yy = y + c(0.25, 0, -0.25)[as.integer(method)],
         arm = factor(arm, levels = c("EUR", "pooled")),
         sig = p_adj_1lmm < 0.05)
count_sig <- function(ph, arm_) {
  d <- prs |> filter(phenotype == ph, trait_arm == arm_)
  sprintf("%d/%d", sum(d$sig), nrow(d))
}
bands <- ROWS |> filter(y %% 2 == 0)
prs1 <- prs |> filter(method == METHODS[["PRSCS"]]) |>
  mutate(trait = ifelse(phenotype == "global_slope", "ΔCT", "CT"),
         yy = y + ifelse(trait == "ΔCT", 0.17, -0.17))
pf <- ggplot(prs1, aes(beta_1lmm, yy, colour = trait)) +
  geom_rect(data = bands, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5), fill = "grey95") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_errorbarh(aes(xmin = beta_1lmm - 1.96 * se_1lmm, xmax = beta_1lmm + 1.96 * se_1lmm),
                 height = 0, linewidth = 0.35) +
  geom_point(aes(shape = arm), fill = "white", size = 1.25, stroke = 0.45) +
  geom_point(data = filter(prs1, sig), aes(shape = arm, fill = trait), size = 1.25,
             stroke = 0.45, show.legend = FALSE) +
  geom_text(data = filter(prs1, sig), aes(x = beta_1lmm + 1.96 * se_1lmm, label = stars(p_adj_1lmm)),
            hjust = -0.15, vjust = 0.75, size = (BASE - 0.5) / .pt, show.legend = FALSE) +
  scale_colour_manual(values = c("ΔCT" = SLOPE_C, "CT" = BASE_C), guide = "none") +
  scale_fill_manual(values = c("ΔCT" = SLOPE_C, "CT" = BASE_C), guide = "none") +
  scale_shape_manual(values = c(EUR = 21, pooled = 23),
                     labels = c(EUR = "EUR arm (n = 4,308)",
                                pooled = "pooled, within-ancestry z (n = 8,596)"), name = NULL) +
  scale_y_continuous(breaks = ROWS$y, labels = ROWS$label, expand = expansion(add = 0.5)) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.1))) +
  labs(x = "PRS-CS score, β per SD (95% CI)", y = NULL, title = "f   Polygenic scores") +
  th + theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(),
             axis.ticks.y = element_blank(),
             axis.text.y = element_text(size = BASE - 1, colour = "grey10"),
             axis.title.x = element_text(size = BASE - 1.5))

# ---------------------------------------------------------- g: symptoms ----
# fig1_prep_cbcl.py: same single-LMM phenotypes as e/f, each fitted alone;
# outcome definitions from ahba_pls/code/23_cbcl_explore.py
cb <- read.delim(file.path(IN, "hcp70_cbcl_assoc.tsv")) |> filter(fit == "alone", outcome != "thought")
CBL <- c(pfactor = "p-factor", internal = "internalising", external = "externalising",
         depress = "depressive (DSM)")
DXL <- c(mdd_youth_DX = "MDD, youth", mdd_parent_DX = "MDD, parent",
         psychosis_parent_DX = "psychosis, parent")
ncase <- cb |> filter(model == "ksads_logit") |> group_by(outcome) |> summarise(k = first(n_cases))
DXL2 <- setNames(sprintf("%s (%s)", DXL, comma(ncase$k[match(names(DXL), ncase$outcome)])), names(DXL))
cb <- cb |>
  mutate(block = ifelse(model == "change", "CBCL change (β per SD)", "KSADS dx (OR per SD)"),
         lab = ifelse(model == "change", CBL[outcome], DXL2[outcome]),
         lab = factor(lab, rev(c(CBL, DXL2))),
         trait = ifelse(brain == "global_slope", "ΔCT", "CT"),
         dy = ifelse(brain == "global_slope", 0.17, -0.17),
         yy = as.numeric(lab) + dy,
         ref = ifelse(model == "change", 0, 1))
hp <- function(o, b) cb$p[cb$outcome == o & cb$brain == b]
he <- function(o, b) cb$est[cb$outcome == o & cb$brain == b]
h_block <- function(d, title, xlab_, key = FALSE) {
  lv <- levels(droplevels(d$lab))
  d <- d |> mutate(yy = match(as.character(lab), lv) + dy)
  ggplot(d, aes(est, yy, colour = trait)) +
    geom_vline(aes(xintercept = ref), linewidth = 0.3, colour = "grey30") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.35) +
    geom_point(shape = 21, fill = "white", size = 1.25, stroke = 0.45) +
    geom_point(data = filter(d, p < 0.05), aes(fill = trait), shape = 21, size = 1.25,
               stroke = 0.45) +
    geom_text(data = filter(d, p < 0.05), aes(x = hi, label = stars(p)), hjust = -0.15,
              vjust = 0.75, size = (BASE - 0.5) / .pt, show.legend = FALSE) +
    scale_colour_manual(values = c("ΔCT" = SLOPE_C, "CT" = BASE_C), breaks = c("ΔCT", "CT"),
                        name = NULL, guide = if (key) "legend" else "none") +
    scale_fill_manual(values = c("ΔCT" = SLOPE_C, "CT" = BASE_C), guide = "none") +
    scale_y_continuous(breaks = seq_along(lv), labels = lv, limits = c(0.5, length(lv) + 0.5),
                       expand = expansion(0)) +
    labs(x = xlab_, y = NULL, title = title) +
    th + theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.y = element_text(size = BASE - 1, colour = "grey10"),
               axis.title.x = element_text(size = BASE - 1.5))
}
n_cb <- n_distinct(cb$outcome[cb$model == "change"]); n_dx <- n_distinct(cb$outcome[cb$model != "change"])
ph1 <- h_block(filter(cb, model == "change"), "g   Symptoms", "CBCL change, β per SD", key = TRUE)
ph2 <- h_block(filter(cb, model != "change"), NULL, "KSADS diagnosis, OR per SD")
ph <- (ph1 / ph2) + plot_layout(heights = c(n_cb, n_dx))

# -------------------------------------------------------- methods panel ----
bullets <- c(
  sprintf("\u2022  ABCD release 7.0, HCP-MMP1.0 (358 parcels; %d lh polygon%s without a value drawn grey); %s children with >= 2 QC-passing scans.",
          n_missing, ifelse(n_missing == 1, "", "s"), comma(n_kids)),
  "\u2022  Slope model: thickness ~ age + sex + (1 + age | child) + (1 | site), per parcel (a, d) and on the per-scan cortical mean (c, e–g); the trait is the child's age slope.",
  "\u2022  PRS: PRS-CS shown (C+T and SBayesRC in SI); EUR arm scored with European discovery GWAS, pooled arm with multi-ancestry GWAS, z-scored within ancestry cluster; model score + age + sex + 10 PCs + (1 | family); filled = p < 0.05.")
SI_NOTE <- "- MAGMA gene-set tests (SCZ/MDD common-variant and exome sets, three LD references) are in Supplementary Information: no set is enriched for ΔCT once ABCD's own EUR children are the LD reference.\n- Sensitivity analyses (DK parcellation, mean of per-parcel slopes, PGC3 2022 GWAS, SBayesR, per-ancestry strata) in Supplementary Information."
ptxt <- ggplot() +
  annotate("text", x = 0, y = 1,
           label = paste(vapply(bullets, function(b) paste(strwrap(b, 150, exdent = 3),
                                collapse = "\n"), ""), collapse = "\n"),
           hjust = 0, vjust = 1, size = (BASE - 1.5) / .pt, lineheight = 1.3,
           colour = "grey25") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(2, 3, 0, 3))

# ----------------------------------------------------------- assemble ------
# row-1 panel height ~ 1.95 in; other columns 1 : 1.5 : 0.78 share the rest of
# the 7.2 in width; choose panel-a's relative width so its cell is as wide as
# the map block needs at that height (+ the colourbar labels)
pa_in <- ROW1_H_IN * PA_ASPECT
PA_W <- pa_in / ((7.2 - pa_in) / (0.9 + 1.3 + 0.73))
row1 <- (pa | pb | pc | pd_) + plot_layout(widths = c(PA_W, 0.9, 1.3, 0.73))
row2 <- (pe | pf | ph) + plot_layout(widths = c(PA_W, (1 + 1.15 + 0.78) / 2, (1 + 1.15 + 0.78) / 2), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0))
fig <- (wrap_elements(full = row1) / wrap_elements(full = row2)) +
  plot_layout(heights = c(0.8, 1.0)) +
  plot_annotation(
    title = "Individual rates of adolescent cortical thinning track genetic risk and symptoms of adolescent-onset mental illness",
    subtitle = paste0("Every region thins across adolescence (a). From 2–4 scans per child (b) a mixed model separates each child's cortical thickness (CT) from their rate of thinning (ΔCT; c–e).\n",
                      "Faster thinning goes with higher polygenic risk for schizophrenia and depression (both ancestry arms) and for Alzheimer's disease (f),\n",
                      "and with rising depressive and externalising symptoms and parent-reported depression (g); CT shows a different, smaller set of associations."),
    theme = theme(plot.title = element_text(size = BASE + 1.5, face = "bold"),
                  plot.subtitle = element_text(size = BASE - 0.5, colour = "grey30",
                                               margin = margin(b = 4))))
ggsave(OUT, fig, width = 7.2, height = 5.0, dpi = 300, bg = "white")
cat(OUT, "\n")

# ------------------------------------------------------------ caption ------
cwm <- setNames(cw$mean, cw$n_visits)
k <- function(arm_, ph = "global_slope") count_sig(ph, arm_)
sp <- function(arm_, ph) fmt_p(prs1$p_adj_1lmm[prs1$trait_arm == arm_ & prs1$phenotype == ph])
cap <- c(
  "**Figure 1 | Individual rates of adolescent cortical thinning track genetic risk and symptoms of adolescent-onset mental illness.**",
  sprintf("CT and ΔCT are the intercept and slope of one linear mixed model per region (a) or on the per-scan cortex-wide mean (b–h): CT is thickness at the sample-mean scan age (%.1f years), not at the first scan, and ΔCT is the per-year change (negative = thinning). Both are fitted from all of a child's scans.", csv_("thickness_mm", "age_centre")),
  "",
  sprintf("- **a** Group maps per parcel (bilateral mean, left hemisphere shown): mean per-child CT (%.1f–%.1f mm) and ΔCT (%.0f to %.1f µm / year).", min(maps$ct_at_centre), max(maps$ct_at_centre), min(um), max(um)),
  sprintf("- **b** Age at scan by visit; %s children, %s scans; 2 / 3 / 4 scans per child for %s / %s / %s children.",
          comma(n_kids), comma(n_scans), comma(per$n_children[per$n_scans == 2]),
          comma(per$n_children[per$n_scans == 3]), comma(per$n_children[per$n_scans == 4])),
  "- **c** Per-scan cortical mean thickness against age for 250 random children (grey) and three highlighted children (points, colour): a child with thick cortex and fast thinning, one with thin cortex and slow thinning (four scans each) and a near-average child (three scans). Each coloured line is that child's fitted trajectory from the mixed model in the title, fixed effects plus the child's random intercept and slope (site effect omitted), so its height at 12.8 years is the child's CT and its slope is their ΔCT. Lines are shrunk toward the population trajectory, most for children with fewer or closely spaced scans.",
  sprintf("- **d** ΔCT reliability 1 − v/τ², where v is a child's conditional (posterior) variance of the slope random effect and τ² the between-child slope variance: boxes, per-parcel LMMs (358 parcels; each value the mean over children; medians %.2f / %.2f / %.2f for 2 / 3 / 4 scans); squares, the single LMM on the cortical mean (mean over children %.2f / %.2f / %.2f).",
          med$m[1], med$m[2], med$m[3], cwm["2"], cwm["3"], cwm["4"]),
  sprintf("- **e** Distributions of whole-cortex CT and ΔCT across the %s children (single LMM on the per-scan cortical mean; CT at the sample-mean scan age, %.1f years). Per-child estimates are shrunk toward the mean (ΔCT SD %.1f µm/yr against a model between-child SD of %.1f µm/yr); every child's estimated ΔCT is negative, and CT and ΔCT are nearly independent across children (r = %.2f).",
          comma(n_kids), csv_("thickness_mm", "age_centre"), csv_("slope_um_per_yr", "sd"),
          csv_("slope_um_per_yr", "model_sd"),
          with(read.csv(file.path(IN, "hcp70_child_traits.csv")), cor(ct_mm, dct_um_per_yr))),
  sprintf("- **f** PRS-CS polygenic score association with ΔCT (red) and CT (blue) on each row; β per SD of score with 95%% CI; filled = p < 0.05. Schizophrenia (2025 GWAS): ΔCT p = %s EUR / %s pooled, CT p = %s / %s. C+T and SBayesRC give the same ΔCT pattern (SI). * p < 0.05, ** p < 0.01, *** p < 0.001 (uncorrected).",
          sp("SCZ25_EUR", "global_slope"), sp("SCZ25_META", "global_slope"),
          sp("SCZ25_EUR", "baseline_thickness"), sp("SCZ25_META", "baseline_thickness")),
  sprintf("- **g** Symptoms against the same single-LMM phenotypes (ΔCT red, CT blue; each fitted alone; filled = p < 0.05; * p < 0.05, ** p < 0.01, *** p < 0.001). CBCL (parent report, log1p raw scale sums): symptoms at ages ~15–17 (mean of the available waves 5–7) adjusted for the same score at baseline (ANCOVA; not a difference score and not a per-child symptom slope), y_late ~ phenotype + y_baseline + age_late + sex + site, β in SD of y_late. KSADS: diagnosis (present or past) at any administered session from baseline to year 6, i.e. lifetime and including baseline cases, logistic with age at the last CBCL, odds ratio per SD. Both phenotypes are signed as in e–f (ΔCT negative = faster thinning), so β < 0 or OR < 1 means more symptoms with faster thinning or with thinner cortex. Estimates are unchanged when both phenotypes enter one model (CT–ΔCT r = 0.06). Family-clustered SEs; n = %s (CBCL change) and %s (KSADS). Faster thinning goes with rising depressive (p = %s), externalising (p = %s) and general (p-factor, p = %s) symptoms and parent-reported MDD (OR %.2f, p = %s); thinner cortex (CT) goes with youth-reported MDD (OR %.2f, p = %s) and parent-reported psychosis spectrum (OR %.2f, p = %s). Exploratory, uncorrected.",
          comma(max(cb$n[cb$model == "change"])), comma(max(cb$n[cb$model == "ksads_logit"])),
          fmt_p(hp("depress", "global_slope")), fmt_p(hp("external", "global_slope")), fmt_p(hp("pfactor", "global_slope")), he("mdd_parent_DX", "global_slope"), fmt_p(hp("mdd_parent_DX", "global_slope")),
          he("mdd_youth_DX", "baseline_thickness"), fmt_p(hp("mdd_youth_DX", "baseline_thickness")),
          he("psychosis_parent_DX", "baseline_thickness"), fmt_p(hp("psychosis_parent_DX", "baseline_thickness"))),
  "", "Methods", sub("^\u2022  ", "- ", bullets), SI_NOTE)
writeLines(cap, CAP)
cat(CAP, "\n")
