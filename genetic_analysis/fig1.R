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
#   hcp70_group_maps.csv              a   per-parcel baseline CT and thinning
#   hcp70_age_by_visit.csv,
#   hcp70_scans_per_child.csv         b   design
#   hcp70_scans.csv (gitignored, individual-level; falls back to
#   hcp70_scans_hist2d.csv), hcp70_scan_summary.csv       c   trajectories
#   hcp70_regional_slope_reliability.csv,
#   hcp70_global_slope_splithalf.csv  d   reliability
#   hcp70_prs_key_arms.tsv            e,f PRS (single-LMM columns)
#   hcp70_magma_locus_sets.tsv (+ work/results_70tab_hcp/magma_pooled/
#   table_magma_pooled.tsv when committed)                g   MAGMA
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
POOLED_MAGMA <- file.path(HERE, "work/results_70tab_hcp/magma_pooled/table_magma_pooled.tsv")

# ---------------------------------------------------------------- style ----
BASE <- 7
SLOPE_C <- "#B2182B"; BASE_C <- "grey40"
METHODS <- c(CT = "C+T", PRSCS = "PRS-CS", SBayesRC = "SBayesRC")  # SBayesR in SI
METHOD_C <- c("C+T" = "grey25", "PRS-CS" = "#3A9AD9", "SBayesRC" = "#C0569E")
VISIT_C <- c(v0 = "#cfe0f2", v2 = "#8fb8de", v4 = "#4f86c6", v6 = "#1f4e8c")
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

brain <- function(values, title, subtitle, scale) {
  d <- poly |> left_join(values, by = "label")
  ggplot(d, aes(x, y, group = grp, fill = value)) +
    geom_polygon(colour = "white", linewidth = 0.06) + scale +
    coord_fixed(expand = FALSE) +
    guides(fill = guide_colourbar(barwidth = unit(60, "pt"), barheight = unit(3, "pt"),
                                  title.position = "left", title.vjust = 1)) +
    labs(title = title) +
    theme_void(base_size = BASE) +
    theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0,
                                    margin = margin(b = 1)),
          plot.subtitle = element_text(size = BASE - 1, colour = "grey30",
                                       hjust = 0, margin = margin(b = 2)),
          plot.title.position = "plot",
          legend.position = "bottom", legend.title = element_text(size = BASE - 1.5),
          legend.text = element_text(size = BASE - 1.5),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "pt"),
          plot.margin = margin(3, 3, 2, 3))
}
lo <- unname(floor(quantile(maps$baseline_ct, 0.02) * 10) / 10)
hi <- unname(ceiling(quantile(maps$baseline_ct, 0.98) * 10) / 10)
um <- maps$slope_mm_per_yr * 1000
vmin <- unname(floor(quantile(um, 0.02) / 5) * 5)
pa1 <- brain(
  data.frame(label = maps$label, value = maps$baseline_ct),
  "a   Baseline thickness",
  sprintf("at age ~10; %.1f–%.1f mm across parcels", min(maps$baseline_ct), max(maps$baseline_ct)),
  scale_fill_gradient(low = "white", high = "#08306B", limits = c(lo, hi),
                      oob = squish, na.value = "grey85", name = "mm",
                      breaks = c(lo, hi)))
pa2 <- brain(
  data.frame(label = maps$label, value = um),
  "     Thinning rate",
  sprintf("every parcel thins: %.0f to %.1f µm/yr", min(um), max(um)),
  scale_fill_distiller(palette = "Reds", direction = -1, limits = c(vmin, 0),
                       oob = squish, na.value = "grey85", name = "µm / yr",
                       breaks = c(vmin, 0)))

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

# ------------------------------------------------------ c: trajectories ----
summ <- read.csv(file.path(IN, "hcp70_scan_summary.csv"))
sv <- setNames(summ$value, summ$metric)
trend <- data.frame(age = c(8.3, 18.2)) |>
  mutate(ct = sv["ols_intercept_mm"] + sv["ols_slope_mm_per_yr"] * age,
         kind = "population trend (OLS)")
KIND_C <- c("one child" = "grey72", "child with 4 scans" = "#1f4e8c",
            "population trend (OLS)" = "black")
if (file.exists(scans_f)) {
  set.seed(7)
  kids <- sample(unique(sc$sid), 250)
  hi3 <- sample(unique(sc$sid[sc$n_visits == 4]), 3)
  lines_bg <- sc |> filter(sid %in% kids) |> mutate(kind = "one child")
  lines_hi <- sc |> filter(sid %in% hi3) |> mutate(kind = "child with 4 scans")
  pc <- ggplot() +
    geom_line(data = lines_bg, aes(age, mean_ct, group = sid, colour = kind),
              linewidth = 0.2, alpha = 0.7) +
    geom_line(data = lines_hi, aes(age, mean_ct, group = sid, colour = kind),
              linewidth = 0.45) +
    geom_point(data = lines_hi, aes(age, mean_ct, colour = kind), size = 0.7)
} else {
  h2 <- read.csv(file.path(IN, "hcp70_scans_hist2d.csv"))
  pc <- ggplot() + geom_tile(data = h2, aes(age_bin + 0.125, ct_bin + 0.01,
                                            alpha = n), fill = "grey50") +
    scale_alpha(guide = "none")
}
pc <- pc +
  geom_line(data = trend, aes(age, ct, colour = kind), linewidth = 0.8) +
  scale_colour_manual(values = KIND_C,
                      labels = c("one child" = "one child (250 shown)",
                                 "child with 4 scans" = "child with 4 scans (3)",
                                 "population trend (OLS)" = "population trend (OLS)")) +
  scale_x_continuous(breaks = seq(8, 18, 2)) +
  coord_cartesian(xlim = c(8.5, 18), ylim = c(2.3, 3.0)) +
  labs(x = "age (years)", y = "mean cortical thickness (mm)",
       title = "c   Each child's slope is the trait") +
  th + theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
             legend.key.width = unit(10, "pt"))

# ------------------------------------------------------- d: reliability ----
rel <- read.csv(file.path(IN, "hcp70_regional_slope_reliability.csv"))
sh <- read.csv(file.path(IN, "hcp70_global_slope_splithalf.csv")) |> filter(n_visits > 0)
med <- rel |> group_by(n_visits) |> summarise(m = median(reliability))
pd_ <- ggplot() +
  geom_boxplot(data = rel, aes(factor(n_visits), reliability),
               outlier.shape = NA, width = 0.55, linewidth = 0.3,
               colour = "grey35", fill = "white", fatten = 0) +
  geom_segment(data = med, aes(x = as.numeric(factor(n_visits)) - 0.27,
                               xend = as.numeric(factor(n_visits)) + 0.27,
                               y = m, yend = m), colour = SLOPE_C, linewidth = 0.6) +
  geom_line(data = sh, aes(factor(n_visits), spearman_brown, group = 1),
            linewidth = 0.4) +
  geom_point(data = sh, aes(factor(n_visits), spearman_brown), shape = 15,
             size = 1.3) +
  annotate("text", x = 3.35, y = sh$spearman_brown[sh$n_visits == 4] - 0.07,
           label = "cortex-wide\n(split-half)", hjust = 1, vjust = 1,
           size = (BASE - 1.5) / .pt, lineheight = 0.9) +
  annotate("text", x = 3.35, y = max(med$m) + 0.22, label = "one parcel\n(358)",
           hjust = 1, size = (BASE - 1.5) / .pt, colour = SLOPE_C, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "scans per child", y = "slope reliability",
       title = "d   Slope reliability") +
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
forest <- function(ph, title, subtitle, colour, show_y) {
  d <- prs |> filter(phenotype == ph)
  ggplot(d, aes(beta_1lmm, yy, colour = method)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
              fill = "grey95") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
    geom_errorbarh(aes(xmin = beta_1lmm - 1.96 * se_1lmm,
                       xmax = beta_1lmm + 1.96 * se_1lmm), height = 0,
                   linewidth = 0.35) +
    geom_point(aes(shape = arm), fill = "white", size = 1.25, stroke = 0.45) +
    geom_point(data = filter(d, sig), aes(shape = arm, fill = method),
               size = 1.25, stroke = 0.45, show.legend = FALSE) +
    scale_colour_manual(values = METHOD_C) +
    scale_fill_manual(values = METHOD_C, guide = "none") +
    scale_shape_manual(values = c(EUR = 21, pooled = 23),
                       labels = c(EUR = "EUR arm (n = 4,308)",
                                  pooled = "pooled, within-ancestry z (n = 8,596)")) +
    scale_y_continuous(breaks = ROWS$y, labels = if (show_y) ROWS$label else NULL,
                       expand = expansion(add = 0.5)) +
    scale_x_continuous(limits = c(-0.09, 0.09), breaks = c(-0.05, 0, 0.05)) +
    labs(x = "β per SD of score (95% CI)", y = NULL, title = title) +
    th + theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.y = element_text(size = BASE - 1, colour = "grey10"),
               plot.title = element_text(face = "bold", size = BASE))
}
pe <- forest("global_slope", "e   PRS → thinning rate",
             sprintf("schizophrenia %s methods EUR, %s pooled",
                     count_sig("global_slope", "SCZ25_EUR"),
                     count_sig("global_slope", "SCZ25_META")), SLOPE_C, TRUE)
pf <- forest("baseline_thickness", "f   PRS → baseline thickness",
             sprintf("schizophrenia %s EUR, %s pooled",
                     count_sig("baseline_thickness", "SCZ25_EUR"),
                     count_sig("baseline_thickness", "SCZ25_META")), BASE_C, FALSE)

# ----------------------------------------------------------- g: MAGMA ------
SETS <- c(SCZ_locus_pool = "SCZ loci", MDD_highconf = "MDD high-confidence")
mg <- read.delim(file.path(IN, "hcp70_magma_locus_sets.tsv"))
if (file.exists(POOLED_MAGMA)) {
  pl <- read.delim(POOLED_MAGMA) |> filter(kind == "gene-set", variable %in% names(SETS))
  if ("construction" %in% names(pl)) pl <- filter(pl, construction == "1lmm")
  pl$arm <- "pooled"
  mg <- bind_rows(mg, pl[, intersect(names(mg), names(pl))])
}
slots <- expand.grid(variable = names(SETS),
                     phenotype = c("global_slope", "baseline_thickness"),
                     arm = c("EUR", "pooled"), stringsAsFactors = FALSE) |>
  mutate(off = case_when(phenotype == "global_slope" & arm == "EUR" ~ 0.30,
                         phenotype == "global_slope" ~ 0.10,
                         arm == "EUR" ~ -0.10, TRUE ~ -0.30),
         y = ifelse(variable == "SCZ_locus_pool", 2, 1) + off) |>
  left_join(mg, by = c("variable", "phenotype", "arm")) |>
  mutate(trait = ifelse(phenotype == "global_slope", "thinning rate", "baseline thickness"),
         arm = factor(arm, levels = c("EUR", "pooled")))
ngenes <- mg |> group_by(variable) |> summarise(n = first(n_genes))
glab <- sprintf("%s\n(%s genes)", SETS, comma(ngenes$n[match(names(SETS), ngenes$variable)]))
mp <- function(v, ph, a = "EUR") {
  x <- mg$p[mg$variable == v & mg$phenotype == ph & mg$arm == a]
  if (length(x)) fmt_p(x) else "n.d."
}
have <- filter(slots, !is.na(beta))
pg <- ggplot(have, aes(beta, y, colour = trait)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5, fill = "grey95") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_errorbarh(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se),
                 height = 0, linewidth = 0.35) +
  geom_point(aes(shape = arm), fill = "white", size = 1.25, stroke = 0.45,
             show.legend = FALSE) +
  geom_point(data = filter(have, p < 0.05), aes(shape = arm, fill = trait),
             size = 1.25, stroke = 0.45, show.legend = FALSE) +
  geom_text(data = filter(have, p < 0.05),
            aes(x = beta + 1.96 * se + 0.012, y = y, label = sprintf("p = %s", fmt_p(p))),
            hjust = 0, size = (BASE - 1.5) / .pt, show.legend = FALSE) +
  geom_text(data = filter(slots, is.na(beta)), aes(0.005, y, label = "n.d."),
            inherit.aes = FALSE, hjust = 0, size = (BASE - 2) / .pt, colour = "grey55") +
  scale_colour_manual(values = c("thinning rate" = SLOPE_C, "baseline thickness" = BASE_C)) +
  scale_fill_manual(values = c("thinning rate" = SLOPE_C, "baseline thickness" = BASE_C),
                    guide = "none") +
  scale_shape_manual(values = c(EUR = 21, pooled = 23), guide = "none") +
  scale_y_continuous(breaks = c(2, 1), labels = glab, limits = c(0.5, 2.5),
                     expand = expansion(0)) +
  scale_x_continuous(limits = c(-0.16, 0.44), breaks = c(0, 0.2, 0.4)) +
  labs(x = "enrichment β (95% CI)", y = NULL, title = "g   MAGMA gene sets") +
  th + theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(),
             axis.ticks.y = element_blank(),
             axis.text.y = element_text(size = BASE - 1, colour = "grey10"),
             legend.position = "none")

# -------------------------------------------------------- methods panel ----
pooled_done <- any(mg$arm == "pooled")
bullets <- c(
  sprintf("\u2022  ABCD release 7.0, HCP-MMP1.0 (358 parcels; %d lh polygon%s without a value drawn grey); %s children with >= 2 QC-passing scans.",
          n_missing, ifelse(n_missing == 1, "", "s"), comma(n_kids)),
  "\u2022  Slope model: thickness ~ age + sex + (1 + age | child) + (1 | site), per parcel (a, d) and on the per-scan cortical mean (c, e–g); the trait is the child's age slope.",
  "\u2022  PRS: C+T, PRS-CS, SBayesRC; EUR arm scored with European discovery GWAS, pooled arm with multi-ancestry GWAS, z-scored within ancestry cluster; model score + age + sex + 10 PCs + (1 | family); filled = p < 0.05 (C+T corrected over its thresholds).",
  sprintf("\u2022  MAGMA competitive gene-set test; EUR arm on 1000 Genomes EUR LD; pooled arm on the ABCD analysis sample as LD reference%s.",
          ifelse(pooled_done, "", " (pending, n.d.)")))
SI_NOTE <- "- Sensitivity analyses (DK parcellation, mean of per-parcel slopes, PGC3 2022 GWAS, SBayesR, per-ancestry strata) in Supplementary Information."
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
row1 <- ((pa1 / pa2) | pb | pc | pd_) + plot_layout(widths = c(1.35, 1, 1.15, 0.75))
row2 <- (pe | pf | pg) + plot_layout(widths = c(1, 0.72, 0.8), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0))
fig <- (row1 / row2) +
  plot_layout(heights = c(1, 1.0)) +
  plot_annotation(
    title = "Polygenic risk for schizophrenia predicts the rate, not the baseline level, of adolescent cortical thinning",
    subtitle = paste0("Every parcel thins (a) and each child's cortex-wide slope is measured reliably (b–d). Schizophrenia PRS predicts faster thinning in both ancestry arms (e)\n",
                      "but not baseline thickness (f); gene-level MAGMA evidence is weaker and arm-dependent (g)."),
    theme = theme(plot.title = element_text(size = BASE + 1.5, face = "bold"),
                  plot.subtitle = element_text(size = BASE - 0.5, colour = "grey30",
                                               margin = margin(b = 4))))
ggsave(OUT, fig, width = 7.2, height = 5.1, dpi = 300, bg = "white")
cat(OUT, "\n")

# ------------------------------------------------------------ caption ------
sb <- setNames(sh$spearman_brown, sh$n_visits)
k <- function(arm_, ph = "global_slope") count_sig(ph, arm_)
cap <- c(
  "**Figure 1 | Polygenic risk for schizophrenia predicts the rate, not the baseline level, of adolescent cortical thinning.**",
  "",
  sprintf("- **a** Group maps per parcel (bilateral mean, left hemisphere shown): baseline thickness and mean per-child thinning rate (%.0f to %.1f µm / year).", min(um), max(um)),
  sprintf("- **b** Age at scan by visit; %s children, %s scans; 2 / 3 / 4 scans per child for %s / %s / %s children.",
          comma(n_kids), comma(n_scans), comma(per$n_children[per$n_scans == 2]),
          comma(per$n_children[per$n_scans == 3]), comma(per$n_children[per$n_scans == 4])),
  sprintf("- **c** Per-scan cortical mean vs age for 250 random children (grey), three with four scans (blue), population OLS trend %.0f µm / year (black).", sv["ols_slope_mm_per_yr"] * 1000),
  sprintf("- **d** Slope reliability: single-parcel model-based reliability across 358 parcels (boxes, medians %.2f / %.2f / %.2f for 2 / 3 / 4 scans) vs split-half consistency of the cortex-wide slope (squares; %.2f / %.2f / %.2f).",
          med$m[1], med$m[2], med$m[3], sb["2"], sb["3"], sb["4"]),
  sprintf("- **e, f** PRS association with thinning rate (e) and baseline thickness (f); β per SD with 95%% CI. Schizophrenia (2025 GWAS): thinning %s EUR / %s pooled, baseline %s / %s. Alzheimer's %s with APOE, %s without; depression %s / %s; education %s (opposite sign); autism %s.",
          k("SCZ25_EUR"), k("SCZ25_META"), k("SCZ25_EUR", "baseline_thickness"),
          k("SCZ25_META", "baseline_thickness"), k("ALZ"), k("ALZ_noAPOE"),
          k("MDD_eur"), k("MDD_pooled"), k("EA"), k("ASD")),
  sprintf("- **g** MAGMA enrichment in SCZ curated loci (Trubetskoy 2022 ST12, %s genes) and MDD high-confidence genes (MDD2025 prioritised; %s genes), p thinning / baseline; p printed on points with p < 0.05. SCZ: EUR %s / %s, pooled %s / %s. MDD: EUR %s / %s, pooled %s / %s. EUR arm on 1000 Genomes EUR LD; pooled arm on the ABCD analysis sample as its own LD reference.",
          comma(ngenes$n[ngenes$variable == "SCZ_locus_pool"][1]),
          comma(ngenes$n[ngenes$variable == "MDD_highconf"][1]),
          mp("SCZ_locus_pool", "global_slope"), mp("SCZ_locus_pool", "baseline_thickness"),
          mp("SCZ_locus_pool", "global_slope", "pooled"), mp("SCZ_locus_pool", "baseline_thickness", "pooled"),
          mp("MDD_highconf", "global_slope"), mp("MDD_highconf", "baseline_thickness"),
          mp("MDD_highconf", "global_slope", "pooled"), mp("MDD_highconf", "baseline_thickness", "pooled")),
  "", "Methods", sub("^\u2022  ", "- ", bullets), SI_NOTE)
writeLines(cap, CAP)
cat(CAP, "\n")
