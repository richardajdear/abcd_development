#!/usr/bin/env Rscript
# fig1.R -- Figure 1 DRAFT: adolescent cortical thinning and the polygenic
# scores that predict it.  HCP-MMP1.0 throughout, ABCD release 7.0.
#
# Usage (repo root):  LC_ALL=en_US.UTF-8 Rscript genetic_analysis/fig1.R
#   (the UTF-8 locale is required: µ, →, β and · are in the source strings)
#
# Reads only committed tables in genetic_analysis/fig1_inputs/ and fits nothing;
# every printed statistic is read from the table plotted, via sprintf.
#   ahba_pls/data/hcp_polygons.csv (lh lateral outline)    a   schematic brain silhouettes
#   hcp70_age_by_visit.csv,
#   hcp70_scans_per_child.csv         b   design
#   hcp70_scans_per_child_all.csv     c   scans per child, incl. 1-scan children (from the run manifest)
#   hcp70_scans.csv (gitignored, individual-level; falls back to
#   hcp70_scans_hist2d.csv), hcp70_scan_summary.csv       d   trajectories
#   hcp70_regional_slope_reliability.csv,
#   hcp70_global_slope_reliability_1lmm.csv   SI  ΔCT reliability (was panel d)
#   hcp70_child_summary.csv               caption: age centre, ΔCT SD (whole-cortex trait panel removed)
#   hcp70_prs_key_arms.tsv            e   PRS (single-LMM columns; fig1_prep_prs.py)
#   hcp70_cbcl_assoc.tsv              f   symptoms at 15–17 y (fig1_prep_cbcl.py)
#   hcp70_dep_bins.csv                g   depressive symptoms vs ΔCT (fig1_prep_cbcl.py)
#   hcp70_prs_mediation.tsv           h, i  mediation through ΔCT (figS_prep_prs_symptoms.py)
#   (MAGMA gene-set panel removed; gene-set results go to SI)
# Writes docs/figures/fig1.png, docs/figures/fig1_caption.md and the former
# panel d as docs/figures/figS_dct_reliability.png.

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
SI_D <- file.path(REPO, "docs/figures/figS_dct_reliability.png")   # former panel d
POLY <- file.path(REPO, "ahba_pls/data/hcp_polygons.csv")   # lh lateral outlines for the a schematic

# ---------------------------------------------------------------- style ----
BASE <- 10
SLOPE_C <- "#440154"; BASE_C <- "#22A884"   # viridis stops: ΔCT dark purple, CT teal-green (a, e-h)
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
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        plot.margin = margin(3, 5, 3, 3))

stars <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
fmt_p <- function(p) ifelse(p < 1e-3, sprintf("%.0e", p),
                     ifelse(p < 0.01, sprintf("%.3f", p), sprintf("%.2f", p)))

# ------------------------------------------------------- a: schematic ------
# graphical abstract, one row per data type (imaging, genetics, symptoms):
# data glyph -> schematic of the processing step -> the per-child measure;
# right-angle links on the right join ΔCT to polygenic scores and to symptoms
n_kids_a <- sum(read.csv(file.path(IN, "hcp70_scans_per_child.csv"))$n_children)
n_scans_a <- sum(read.csv(file.path(IN, "hcp70_age_by_visit.csv"))$count)
FS <- (BASE - 1.5) / .pt
lat <- read.csv(POLY) |> filter(view == "lateral")         # lh lateral: frontal left, occipital right
BW <- 1.0
sc_ <- BW / diff(range(lat$x))
lat <- lat |> mutate(x = (x - min(x)) * sc_, y = (y - min(y)) * sc_,
                     grp = interaction(label, group, subgroup, drop = TRUE))
# cerebellum (with folia) and brainstem, tucked under the occipital / temporal lobes
ell <- function(cx, cy, rx, ry, n = 60) { t <- seq(0, 2 * pi, length.out = n); data.frame(x = cx + rx * cos(t), y = cy + ry * sin(t)) }
cereb <- ell(0.79, 0.0, 0.19, 0.115)
folia <- bind_rows(lapply(1:4, function(k) { t <- seq(-0.9, 0.9, length.out = 25)
  data.frame(k = k, x = 0.79 + (0.17 - 0.03 * k) * sin(t), y = 0.0 + 0.06 - 0.04 * k + 0.02 * cos(t)) }))
stem <- data.frame(x = c(0.56, 0.66, 0.65, 0.58), y = c(0.12, 0.12, -0.08, -0.08))
BR_FILL <- scales::gradient_n_pal(c(BASE_C, "#8FD3BE"))(seq(0, 1, length.out = 4))
BR_Y0 <- 4.18; BR_DX <- 0.5
br_layers <- unlist(lapply(1:4, function(i) {
  ox <- 0.02 + (i - 1) * BR_DX; oy <- BR_Y0
  cb <- mutate(cereb, x = x + ox, y = y + oy); st <- mutate(stem, x = x + ox, y = y + oy)
  fo <- mutate(folia, x = x + ox, y = y + oy); cx <- mutate(lat, x = x + ox, y = y + oy, grp = paste(i, grp))
  sub_c <- scales::gradient_n_pal(c(BR_FILL[i], "grey45"))(0.35)          # cerebellum / stem: greyer than cortex
  list(geom_polygon(data = st, aes(x, y), fill = sub_c, colour = "white", linewidth = 0.6),
       geom_polygon(data = cb, aes(x, y), fill = sub_c, colour = "white", linewidth = 0.6),
       geom_path(data = fo, aes(x, y, group = k), colour = alpha("white", 0.8), linewidth = 0.2),
       geom_polygon(data = cx, aes(x, y, group = grp), fill = "white", colour = "white", linewidth = 0.9),
       geom_polygon(data = cx, aes(x, y, group = grp), fill = BR_FILL[i], colour = alpha("white", 0.35),
                    linewidth = 0.05))
}), recursive = FALSE)
# DNA double helix
HY <- 2.45
hx <- seq(0, 2.2, length.out = 200)
dna <- bind_rows(data.frame(x = 0.1 + hx, y = HY + 0.3 * sin(hx * 2 * pi / 0.85), s = "a"),
                 data.frame(x = 0.1 + hx, y = HY - 0.3 * sin(hx * 2 * pi / 0.85), s = "b"))
rx <- seq(0.06, 2.15, by = 0.0708)
rungs <- data.frame(x = 0.1 + rx, y0 = HY + 0.3 * sin(rx * 2 * pi / 0.85), y1 = HY - 0.3 * sin(rx * 2 * pi / 0.85),
                    col = rep(c("#E08214", "#2166AC", "#B2182B", "#8C8C8C"), length.out = length(rx)))
# questionnaire icon
QX <- 0.8; QY <- 0.42
q_lines <- data.frame(y = QY + c(0.72, 0.5, 0.28))
# processing-step schematics, x 3.15-5.25
PX0 <- 3.15; PX1 <- 5.25
mini <- data.frame(k = rep(1:3, each = 2), x = rep(c(3.35, 5.05), 3)) |>
  mutate(y = 4.42 + c(0.14, 0, -0.14)[k] + c(-0.05, -0.18, -0.32)[k] * (x - 4.2))
mpts <- data.frame(k = rep(1:3, each = 3), x = rep(c(3.5, 4.2, 4.9), 3)) |>
  mutate(y = 4.42 + c(0.14, 0, -0.14)[k] + c(-0.05, -0.18, -0.32)[k] * (x - 4.2) +
           c(0.03, -0.04, 0.02, -0.02, 0.04, -0.02, 0.02, -0.03, 0.02))
set.seed(3)                                                # Manhattan plot (illustrative)
nm_ <- 140; mh <- data.frame(i = seq_len(nm_)) |>
  mutate(x = PX0 + 0.15 + (i - 1) / (nm_ - 1) * (PX1 - PX0 - 0.25), chr = (i - 1) %/% 14,
         h = rexp(nm_, 1 / 0.12))
for (pk in c(22, 61, 97, 124)) mh$h[pk + (-2:2)] <- c(0.25, 0.45, 0.62, 0.42, 0.22) + c(0, 0, runif(1, 0, 0.1), 0, 0)
mh <- mh |> mutate(h = pmin(h, 0.68), y = HY - 0.36 + h, sig = h > 0.5)
MH_SIG <- HY - 0.36 + 0.5
set.seed(5)                                                # late vs baseline symptoms (illustrative)
bb <- data.frame(x = runif(26, 0.1, 0.9)) |> mutate(y = 0.15 + 0.6 * x + rnorm(26, 0, 0.09))
bb <- bind_rows(bb, data.frame(x = 0.42, y = 0.15 + 0.6 * 0.42 + 0.27))       # the highlighted child
SX <- function(v) PX0 + 0.25 + v * (PX1 - PX0 - 0.4); SY <- function(v) QY - 0.02 + v * 0.95
bb <- mutate(bb, X = SX(x), Y = SY(y), hi = row_number() == n())
chip <- function(x0, x1, y0, y1, fill, col, lab, tcol = col)
  list(annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill, colour = col, linewidth = 0.35),
       annotate("text", x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = lab, colour = tcol,
                size = FS + 0.1, fontface = "bold", lineheight = 0.85))
arr <- function(x, y, xend, yend, col = "grey35")
  annotate("segment", x = x, y = y, xend = xend, yend = yend, colour = col, linewidth = 0.4,
           arrow = arrow(length = unit(3, "pt"), type = "closed"))
sub_ <- function(x, y, lab) annotate("text", x = x, y = y, label = lab, size = FS - 0.4, colour = "grey35",
                                     lineheight = 0.9)
ttl_ <- function(y, lab) annotate("text", x = (PX0 + PX1) / 2 + 0.05, y = y, label = lab, size = FS,
                                  fontface = "bold", colour = "grey20")
axes_ <- function(y0, y1) list(
  annotate("segment", x = PX0 + 0.1, xend = PX0 + 0.1, y = y0, yend = y1, colour = "grey50", linewidth = 0.3),
  annotate("segment", x = PX0 + 0.1, xend = PX1, y = y0, yend = y0, colour = "grey50", linewidth = 0.3))
CX0 <- 5.8; CX1 <- 8.2                                     # per-child measure chips
DCT_Y <- 3.98; SYM_Y <- QY + 0.48
L1 <- CX1 + 0.4; L2 <- CX1 + 0.8                           # link tracks (non-crossing)
link <- function(pts) annotate("path", x = pts$x, y = pts$y, colour = SLOPE_C, linewidth = 0.45,
                               linejoin = "mitre", arrow = arrow(length = unit(3, "pt"), type = "closed"))
pa <- ggplot() +
  # ---- imaging
  br_layers +
  sub_(1.27, BR_Y0 - 0.4, sprintf("MRI, 9–17 y\n%s children · %s scans", comma(n_kids_a), comma(n_scans_a))) +
  arr(2.68, 4.42, PX0 - 0.05, 4.42) +
  axes_(3.85, 5.0) +
  geom_line(data = mini, aes(x, y, group = k), colour = SLOPE_C, linewidth = 0.45) +
  geom_point(data = mpts, aes(x, y), colour = "grey30", size = 0.45) +
  ttl_(5.22, "mixed model") +
  arr(PX1 + 0.05, 4.65, CX0 - 0.05, 4.86, col = BASE_C) + arr(PX1 + 0.05, 4.2, CX0 - 0.05, DCT_Y, col = SLOPE_C) +
  chip(CX0, CX1, 4.62, 5.22, alpha(BASE_C, 0.10), BASE_C, "CT intercept") +
  chip(CX0, CX1, DCT_Y - 0.3, DCT_Y + 0.3, alpha(SLOPE_C, 0.12), SLOPE_C, "ΔCT slope") +
  # ---- genetics
  geom_path(data = dna, aes(x, y, group = s), colour = "grey25", linewidth = 0.5) +
  geom_segment(data = rungs, aes(x = x, xend = x, y = y0, yend = y1), colour = rungs$col, linewidth = 0.4) +
  sub_(1.2, HY - 0.55, sprintf("genotypes\n%s children", comma(8596))) +
  arr(2.45, HY, PX0 - 0.05, HY) +
  axes_(HY - 0.36, HY + 0.4) +
  geom_point(data = mh, aes(x, y), colour = ifelse(mh$sig, "grey10", ifelse(mh$chr %% 2 == 0, "grey55", "grey78")),
             size = 0.35, stroke = 0) +
  annotate("segment", x = PX0 + 0.1, xend = PX1, y = MH_SIG, yend = MH_SIG, colour = "grey30", linewidth = 0.25,
           linetype = "22") +
  ttl_(HY + 0.62, "GWAS") +
  arr(PX1 + 0.05, HY, CX0 - 0.05, HY) +
  chip(CX0, CX1, HY - 0.36, HY + 0.36, "grey96", "grey40", "polygenic\nscores", tcol = "grey15") +
  # ---- symptoms
  annotate("rect", xmin = QX, xmax = QX + 0.8, ymin = QY, ymax = QY + 0.95, fill = "white", colour = "grey35", linewidth = 0.35) +
  annotate("rect", xmin = QX + 0.25, xmax = QX + 0.55, ymin = QY + 0.88, ymax = QY + 1.02, fill = "grey35", colour = NA) +
  annotate("point", x = QX + 0.16, y = q_lines$y, shape = 22, size = 0.9, stroke = 0.3, colour = "grey35", fill = "white") +
  annotate("segment", x = QX + 0.3, xend = QX + 0.68, y = q_lines$y, yend = q_lines$y, colour = "grey55", linewidth = 0.35) +
  sub_(1.2, QY - 0.32, "CBCL, 9–10 y\nand 15–17 y") +
  arr(1.75, SYM_Y, PX0 - 0.05, SYM_Y) +
  axes_(QY - 0.05, QY + 1.0) +
  annotate("segment", x = SX(0.05), xend = SX(0.95), y = SY(0.15 + 0.6 * 0.05), yend = SY(0.15 + 0.6 * 0.95),
           colour = "grey30", linewidth = 0.35) +
  annotate("segment", x = SX(0.42), xend = SX(0.42), y = SY(0.15 + 0.6 * 0.42), yend = SY(0.15 + 0.6 * 0.42 + 0.27),
           colour = SLOPE_C, linewidth = 0.35, linetype = "22") +
  geom_point(data = filter(bb, !hi), aes(X, Y), colour = "grey60", size = 0.45) +
  geom_point(data = filter(bb, hi), aes(X, Y), colour = SLOPE_C, size = 0.9) +
  annotate("text", x = PX1 - 0.02, y = QY - 0.14, label = "9–10 y", hjust = 1, size = FS - 0.6, colour = "grey40") +
  annotate("text", x = PX0 + 0.18, y = QY + 1.0, label = "15–17 y", hjust = 0, vjust = 1, size = FS - 0.6,
           colour = "grey40") +
  ttl_(QY + 1.2, "beyond baseline") +
  arr(PX1 + 0.05, SYM_Y, CX0 - 0.05, SYM_Y) +
  chip(CX0, CX1, SYM_Y - 0.36, SYM_Y + 0.36, "grey96", "grey40", "symptoms\n15–17 y", tcol = "grey15") +
  # ---- associations: polygenic scores -> ΔCT (inner track), ΔCT -> symptoms (outer track)
  link(data.frame(x = c(CX1 + 0.04, L1, L1, CX1 + 0.06), y = c(HY, HY, DCT_Y - 0.12, DCT_Y - 0.12))) +
  link(data.frame(x = c(CX1 + 0.04, L2, L2, CX1 + 0.06), y = c(DCT_Y + 0.12, DCT_Y + 0.12, SYM_Y, SYM_Y))) +
  coord_fixed(xlim = c(0, L2 + 0.15), ylim = c(-0.15, 5.4), expand = FALSE, clip = "off") +
  labs(title = "a   Study design") +
  theme_void(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0, margin = margin(b = 3)),
        plot.title.position = "plot", plot.margin = margin(3, 6, 3, 3))

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
  annotate("text", x = age$mean, y = ymax * c(1.42, 1.2, 1.42, 1.2),
           label = VISIT_L[age$visit], colour = c("grey45", VISIT_C[-1]),
           size = (BASE - 1.5) / .pt) +
  scale_fill_manual(values = VISIT_C, guide = "none") +
  scale_x_continuous(breaks = seq(8, 18, 2), limits = c(8, 18.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), labels = comma) +
  coord_cartesian(ylim = c(0, ymax * 1.5)) +
  labs(x = "age at scan (years)", y = "scans",
       title = sprintf("b   %s children, %s scans", comma(n_kids), comma(n_scans))) +
  th


# per-child phenotype distributions (fig1_prep_1lmm.R): the two traits of e-h
cd <- read.csv(file.path(IN, "hcp70_child_density.csv"))
cs <- read.csv(file.path(IN, "hcp70_child_summary.csv"))
csv_ <- function(v, f) cs[[f]][cs$var == v]
# (whole-cortex trait distributions, former panel e, removed; cd / cs kept for the caption)

# ------------------------------------------------------- c: scans per child ----
# hcp70_scans_per_child_all.csv: children with 4 / 3 / 2 / 1 QC-passing scans,
# from out/thickness_hcp_70_aa6e91efba82/manifest.json (n_visits_distribution;
# 1-scan children = n_subjects_before - n_subjects_after of the min_visits rule).
# Only children with >= 2 scans enter the slope model.
spc <- read.csv(file.path(IN, "hcp70_scans_per_child_all.csv")) |>
  mutate(nlab = factor(n_scans, levels = c(4, 3, 2, 1)), grp = ifelse(in_model == 1, "in model", "excluded"))
p_spc <- ggplot(spc, aes(nlab, n_children, fill = grp)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = comma(n_children)), vjust = -0.3, size = (BASE - 1.5) / .pt, colour = "grey20") +
  annotate("text", x = 4, y = spc$n_children[spc$n_scans == 1] * 2.05, label = "not\nmodelled", size = (BASE - 2) / .pt,
           colour = "grey45", lineheight = 0.85) +
  scale_fill_manual(values = c("in model" = "grey35", "excluded" = "grey80"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22)), labels = comma) +
  labs(x = "QC-passing scans per child", y = "children", title = "c   Scans per child") +
  th + theme(panel.grid.major.x = element_blank())

# ------------------------------------------------------ d: trajectories ----
summ <- read.csv(file.path(IN, "hcp70_scan_summary.csv"))
sv <- setNames(summ$value, summ$metric)
YL <- c(2.3, 3.0)
trend <- data.frame(age = c(8.3, 18.2)) |>
  mutate(ct = sv["ols_intercept_mm"] + sv["ols_slope_mm_per_yr"] * age)
HI_C <- c("#d95f02", "#1b9e77", "#7570b3")               # highlighted children (fast, slow, average): Figure 2 panel c colours (Dark2)
if (file.exists(scans_f)) {
  set.seed(7)
  kids <- sample(unique(sc$sid), 250)
  set.seed(25)                                   # highlighted-child draw (grey sample unchanged)
  # highlighted children, chosen on the fitted ΔCT so their slopes differ as much
  # as the (shrunk) distribution allows: a fast thinner (z < -2), a near-average
  # child and a slow thinner (z > 2); 4 scans spanning >= 5.5 y, smooth
  # trajectories only. CT offsets (fast thinner thinner, slow thinner thicker)
  # keep the three fitted lines apart so they fan out with age.
  tr <- read.csv(file.path(IN, "hcp70_child_traits.csv"))
  cen <- csv_("thickness_mm", "age_centre")
  zz <- function(v) (v - mean(v)) / sd(v)
  cand <- sc |> group_by(sid, n_visits) |> summarise(.groups = "drop") |>
    inner_join(mutate(tr, zct = zz(ct_mm), zdct = zz(dct_um_per_yr)), by = "sid") |>
    inner_join(sc |> inner_join(tr, by = "sid") |>
                 mutate(res = mean_ct - (ct_mm + dct_um_per_yr / 1000 * (age - cen))) |>
                 group_by(sid) |> summarise(rsd = sd(res), span = max(age) - min(age)),
               by = "sid") |>
    filter(rsd < 0.03, span >= 5.5)
  pick <- function(n, f) {                 # x[sample.int()] -- sample(x, 1) on a length-1 x samples 1:x
    x <- filter(cand, n_visits == n, !!f)$sid
    stopifnot(length(x) > 0); x[sample.int(length(x), 1)]
  }
  hi3 <- c(pick(4, quote(zct < -0.5 & zct > -1.5 & zdct < -2.0)),
           pick(4, quote(zct > 0.5 & zct < 1.5 & zdct > 2.0)),
           pick(4, quote(abs(zct) < 0.3 & abs(zdct) < 0.2)))
  fit_hi <- tr |> filter(sid %in% hi3) |> mutate(child = factor(match(sid, hi3))) |>
    inner_join(sc |> group_by(sid) |> summarise(a0 = min(age), a1 = max(age)), by = "sid") |>
    tidyr::pivot_longer(c(a0, a1), values_to = "age") |>
    mutate(ct = ct_mm + dct_um_per_yr / 1000 * (age - cen))
  hi_dct <- tr$dct_um_per_yr[match(hi3, tr$sid)]
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
              linewidth = 0.75),
    geom_point(data = lines_hi, aes(age, mean_ct), colour = HI_C[lines_hi$child], size = 1.7))} +
  scale_colour_manual(values = c("one child (250 shown)" = "grey72"),
                      breaks = "one child (250 shown)") +
  scale_x_continuous(breaks = seq(8, 18, 2)) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  coord_cartesian(xlim = c(8.5, 18)) +
  labs(x = "age (years)", y = "cortex-wide CT (mm)",
       title = "d   CT ~ age + sex + (age|child) + (1|site)") +
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
  "ASD",        "autism",                   "EUR")
prs_raw <- read.delim(file.path(IN, "hcp70_prs_key_arms.tsv"))  # fig1_prep_prs.py
# bipolar / ADHD are in the PRS table but shown in SI only
ROWS <- ROWS |> mutate(y = rev(seq_len(n())))
prs <- prs_raw |>
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
                 height = 0, linewidth = 0.5) +
  geom_point(aes(shape = arm), fill = "white", size = 1.9, stroke = 0.6) +
  geom_point(data = filter(prs1, sig), aes(shape = arm, fill = trait), size = 1.9,
             stroke = 0.45, show.legend = FALSE) +
  geom_text(data = filter(prs1, sig), aes(x = beta_1lmm + 1.96 * se_1lmm, label = stars(p_adj_1lmm)),
            hjust = -0.15, vjust = 0.75, size = (BASE - 0.5) / .pt, show.legend = FALSE) +
  scale_colour_manual(values = c("ΔCT" = SLOPE_C, "CT" = BASE_C), breaks = c("ΔCT", "CT"), name = NULL,
                      guide = guide_legend(order = 1, override.aes = list(shape = 21, fill = "white"))) +
  scale_fill_manual(values = c("ΔCT" = SLOPE_C, "CT" = BASE_C), guide = "none") +
  scale_shape_manual(values = c(EUR = 21, pooled = 23),
                     labels = c(EUR = "EUR arm (n = 4,308)",
                                pooled = "pooled, within-ancestry z (n = 8,596)"), name = NULL,
                     guide = guide_legend(order = 2)) +
  scale_y_continuous(breaks = ROWS$y, labels = ROWS$label, expand = expansion(add = 0.5)) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.1))) +
  labs(x = "PRS-CS score, β per SD (95% CI)", y = NULL, title = "e   Polygenic scores") +
  th + theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(),
             axis.ticks.y = element_blank(),
             axis.text.y = element_text(size = BASE - 1, colour = "grey10"),
             axis.title.x = element_text(size = BASE - 1.5))

# ---------------------------------------------------------- g: symptoms ----
# fig1_prep_cbcl.py: same single-LMM phenotypes as e/f, each fitted alone;
# outcome definitions from ahba_pls/code/23_cbcl_explore.py
cb_all <- read.delim(file.path(IN, "hcp70_cbcl_assoc.tsv")) |> filter(fit == "alone")
CBL <- c(pfactor = "p-factor", internal = "internalising", external = "externalising",
         depress = "depressive (DSM)")   # thought problems and PQ-BC rows stay in the table (SI)
cb <- cb_all |>
  filter(model == "change", outcome %in% names(CBL)) |>
  mutate(lab = factor(CBL[outcome], rev(CBL)),
         trait = ifelse(brain == "global_slope", "ΔCT", "CT"),
         dy = ifelse(brain == "global_slope", 0.17, -0.17),
         yy = as.numeric(lab) + dy)
hp <- function(o, b) cb_all$p[cb_all$outcome == o & cb_all$brain == b & cb_all$model == "change"]
he <- function(o, b) cb_all$est[cb_all$outcome == o & cb_all$brain == b & cb_all$model == "change"]
ph <- ggplot(cb, aes(est, yy, colour = trait)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5) +
  geom_point(shape = 21, fill = "white", size = 1.9, stroke = 0.6) +
  geom_point(data = filter(cb, p < 0.05), aes(fill = trait), shape = 21, size = 1.9, stroke = 0.6) +
  geom_text(data = filter(cb, p < 0.05), aes(x = hi, label = stars(p)), hjust = -0.15,
            vjust = 0.75, size = (BASE - 0.5) / .pt, show.legend = FALSE) +
  scale_colour_manual(values = c("ΔCT" = SLOPE_C, "CT" = BASE_C), breaks = c("ΔCT", "CT"), name = NULL,
                      guide = "none") +
  scale_fill_manual(values = c("ΔCT" = SLOPE_C, "CT" = BASE_C), guide = "none") +
  scale_y_continuous(breaks = seq_along(levels(cb$lab)), labels = levels(cb$lab),
                     limits = c(0.5, length(CBL) + 0.5), expand = expansion(0)) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.08))) +
  labs(x = "symptoms at 15–17 y, adjusted for baseline\nβ per SD of trait (95% CI)", y = NULL,
       title = "f   Symptoms") +
  th + theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(),
             axis.ticks.y = element_blank(),
             axis.text.y = element_text(size = BASE - 1, colour = "grey10"),
             axis.title.x = element_text(size = BASE - 1.5))

# --------------------------------------------- h: depressive symptoms ----
# fig1_prep_cbcl.py: y = DSM-depressive score at 15-17 y residualised on its
# baseline score + age + sex + site (SD units); x = ΔCT. Binned scatter: per-child
# residual SD ~1 against a ~0.03 SD/SD slope, so individual points would sit
# off the visible axis; means of 20 equal-count ΔCT bins with 95% CIs.
hb <- read.csv(file.path(IN, "hcp70_dep_bins.csv"))
xr <- c(-25, -15)
fit_h <- data.frame(x = xr) |> mutate(y = hb$fit_intercept[1] + hb$fit_slope_per_um[1] * x)
p_h <- ggplot()
p_h <- p_h +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_line(data = fit_h, aes(x, y), colour = SLOPE_C, linewidth = 0.75) +
  geom_errorbar(data = hb, aes(x, ymin = y - 1.96 * se, ymax = y + 1.96 * se), width = 0,
                colour = "grey20", linewidth = 0.45) +
  geom_point(data = hb, aes(x, y), colour = "grey20", size = 1.4) +
  annotate("text", x = xr[2], y = 0.3, hjust = 1, vjust = 1, size = (BASE - 1.5) / .pt,
           label = sprintf("%s children, p = %s", comma(hb$n_total[1]), fmt_p(hb$fit_p[1])),
           colour = "grey25") +
  annotate("text", x = xr[1], y = -0.27, hjust = 0, vjust = 0, size = (BASE - 1.5) / .pt,
           label = "\u2190 faster thinning", colour = SLOPE_C) +
  scale_x_continuous(limits = xr, breaks = seq(-24, -16, 2)) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  labs(x = "ΔCT (µm / year)", y = "depressive symptoms\nbeyond baseline (SD)",
       title = "g   Depressive symptoms at 15–17 y") +
  th + theme(axis.title = element_text(size = BASE - 1.5), plot.margin = margin(3, 5, 14, 3))

# ------------------------------------------- g, h: mediation through ΔCT ----
# figS_prep_prs_symptoms.py: product-of-coefficients mediation, pooled arm,
# beyond-baseline covariates; 95% CIs from family-resampling bootstrap. For each
# score the outcome shown is its strongest total (score -> symptom) effect among
# the four symptom rows of e; every outcome and both arms are in SI.
MED <- read.delim(file.path(IN, "hcp70_prs_mediation.tsv"))
DIS_C <- c(SCZ = "#E08214", MDD = "#C0569E")              # disorder colours, as Figure 2
DIS_L <- c(SCZ = "SCZ\nPRS", MDD = "MDD\nPRS")
OUT_S <- c(pfactor = "p-factor", internal = "internalising", external = "externalising", depress = "depressive")
pick_med <- function(tr) MED |> filter(trait == tr, arm == "pooled", outcome %in% names(CBL)) |>
  slice_max(c_total, n = 1)
med_g <- pick_med("SCZ"); med_h <- pick_med("MDD")
# reliability of the per-child ΔCT estimate: var(BLUP) / model between-child
# variance (= mean per-child reliability 1 - v/tau^2, 0.30). The BLUP regressed
# on the score is shrunk by ~R, so a·b / R approximates the indirect effect on the
# true thinning rate (CI scaled by 1/R; uncertainty in R ignored)
REL <- (csv_("slope_um_per_yr", "sd") / csv_("slope_um_per_yr", "model_sd"))^2
path_panel <- function(r, title) {
  st <- function(est, p) sprintf("%+.3f%s", est, ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
  FSp <- (BASE - 2) / .pt
  box <- function(x0, x1, y0, y1, f, cl) annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = f,
                                                 colour = cl, linewidth = 0.4)
  ar <- function(x, y, xe, ye, cl) annotate("segment", x = x, y = y, xend = xe, yend = ye, colour = cl,
                                           linewidth = 0.5, arrow = arrow(length = unit(3.5, "pt"), type = "closed"))
  dc <- DIS_C[r$trait]
  pct <- if (r$c_p < 0.05) sprintf("%.1f%% of the total", 100 * r$prop) else "total not significant"
  ggplot() +
    box(0, 2.4, 0.55, 1.65, alpha(dc, 0.12), dc) +
    annotate("text", x = 1.2, y = 1.1, label = DIS_L[r$trait], size = FSp, fontface = "bold", colour = dc,
             lineheight = 0.85) +
    box(3.7, 6.3, 2.55, 3.45, alpha(SLOPE_C, 0.12), SLOPE_C) +
    annotate("text", x = 5, y = 3.0, label = "ΔCT", size = FSp + 0.3, fontface = "bold", colour = SLOPE_C) +
    box(6.9, 10, 0.55, 1.65, "grey96", "grey40") +
    annotate("text", x = 8.45, y = 1.1, label = paste0(OUT_S[r$outcome], "\n15–17 y"), size = FSp, fontface = "bold",
             colour = "grey15", lineheight = 0.85) +
    ar(1.7, 1.7, 3.65, 2.7, SLOPE_C) + ar(6.35, 2.7, 8.2, 1.7, SLOPE_C) +
    ar(2.45, 1.1, 6.85, 1.1, "grey35") +
    annotate("text", x = 2.5, y = 2.4, label = paste("a", st(r$a, r$a_p)), hjust = 1, size = FSp, colour = SLOPE_C) +
    annotate("text", x = 7.5, y = 2.4, label = paste("b", st(r$b, r$b_p)), hjust = 0, size = FSp, colour = SLOPE_C) +
    annotate("text", x = 4.65, y = 1.25, label = sprintf("c′ %+.3f", r$c_direct), size = FSp - 0.15, vjust = 0,
             colour = "grey25") +
    annotate("text", x = 4.65, y = 0.95, label = sprintf("total %s", st(r$c_total, r$c_p)), size = FSp - 0.15,
             vjust = 1, colour = "grey40") +
    annotate("text", x = 5, y = 0.2, vjust = 1, size = FSp, colour = SLOPE_C, lineheight = 0.95,
             label = sprintf("via ΔCT: a·b = %+.4f (%+.4f to %+.4f)\n%s", r$indirect, r$ind_lo, r$ind_hi, pct)) +
    annotate("text", x = 5, y = -0.95, vjust = 1, size = FSp - 0.15, colour = "grey40", lineheight = 0.95,
             label = sprintf("reliability-corrected (R = %.2f): %+.4f%s", REL, r$indirect / REL,
                             if (r$c_p < 0.05) sprintf(", %.0f%% of the total", 100 * r$prop / REL) else "")) +
    coord_cartesian(xlim = c(0, 10), ylim = c(-1.35, 3.55), expand = FALSE, clip = "off") +
    labs(title = title) + theme_void(base_size = BASE) +
    theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0, margin = margin(b = 4)),
          plot.title.position = "plot", plot.margin = margin(3, 8, 6, 3))
}
p_g <- path_panel(med_g, "h   Schizophrenia PRS, via ΔCT")
p_hm <- path_panel(med_h, "i   Depression PRS, via ΔCT")

# -------------------------------------------------------- methods panel ----
bullets <- c(
  sprintf("\u2022  ABCD release 7.0, HCP-MMP1.0 (358 parcels, cortex-wide mean thickness per scan); %s children with >= 2 QC-passing scans.", comma(n_kids)),
  "\u2022  Slope model: thickness ~ age + sex + (1 + age | child) + (1 | site), on the per-scan cortical mean (d–i); the trait is the child's age slope.",
  "\u2022  PRS: PRS-CS shown (C+T and SBayesRC in SI); EUR arm scored with European discovery GWAS, pooled arm with multi-ancestry GWAS, z-scored within ancestry cluster; model score + age + sex + 10 PCs + (1 | family); filled = p < 0.05.")
SI_NOTE <- "- ΔCT reliability by scans per child is in figS_dct_reliability.png.\n- Score → symptom associations for all outcomes and both arms, and mediation through ΔCT for every outcome, are in figS_prs_symptoms.png.\n- MAGMA gene-set tests (SCZ/MDD common-variant and exome sets, three LD references) are in Supplementary Information: no set is enriched for ΔCT once ABCD's own EUR children are the LD reference.\n- Sensitivity analyses (DK parcellation, mean of per-parcel slopes, PGC3 2022 GWAS, SBayesR, per-ancestry strata) in Supplementary Information."
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
# 10.8 x 7.5 in (Figure 2 width, ~13% shorter). Row 1: a | b over c | d; row 2: e | f over g | h over i (mediation)
pbc <- (pb / p_spc) + plot_layout(heights = c(1, 1))
row1 <- (pa | pbc | pc) + plot_layout(widths = c(2.25, 0.85, 1.0))
# f, g, h, i in ONE 2 x 2 grid so the mediation diagrams share rows with f and g
# (titles line up exactly); g and i are both wrapped so their titles sit at the
# top of the same row (g keeps its y-title next to its own axis)
pfghi <- ph + wrap_elements(full = p_g) + wrap_elements(full = p_h) + wrap_elements(full = p_hm) +
  plot_layout(ncol = 2, widths = c(1.1, 1.15), heights = c(1, 1.15))
row2 <- (pf | pfghi) + plot_layout(widths = c(1.05, 2.25), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.justification = "left",
        legend.box.just = "left", legend.margin = margin(6, 0, 0, 0))
fig <- (wrap_elements(full = row1) / wrap_elements(full = row2)) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Individual rates of adolescent cortical thinning track genetic risk and symptoms of adolescent-onset mental illness",
    theme = theme(plot.title = element_text(size = BASE + 2, face = "bold", margin = margin(b = 4))))
ggsave(OUT, fig, width = 10.8, height = 7.5, dpi = 300, bg = "white")
ggsave(SI_D, pd_ + labs(title = "ΔCT reliability"), width = 2.4, height = 2.2, dpi = 300, bg = "white")
cat(OUT, "\n")

# ------------------------------------------------------------ caption ------
cwm <- setNames(cw$mean, cw$n_visits)
k <- function(arm_, ph = "global_slope") count_sig(ph, arm_)
sp <- function(arm_, ph) fmt_p(prs1$p_adj_1lmm[prs1$trait_arm == arm_ & prs1$phenotype == ph])
cap <- c(
  "**Figure 1 | Individual rates of adolescent cortical thinning track genetic risk and symptoms of adolescent-onset mental illness.**",
  sprintf("CT and ΔCT are the intercept and slope of one linear mixed model on the per-scan cortex-wide mean: CT is thickness at the sample-mean scan age (%.1f years), not at the first scan, and ΔCT is the per-year change (negative = thinning). Both are fitted from all of a child's scans. Per-child estimates are shrunk toward the mean (ΔCT SD %.1f µm/yr against a model between-child SD of %.1f µm/yr), and CT and ΔCT are nearly independent across children (r = %.2f).", csv_("thickness_mm", "age_centre"),
          csv_("slope_um_per_yr", "sd"), csv_("slope_um_per_yr", "model_sd"),
          with(read.csv(file.path(IN, "hcp70_child_traits.csv")), cor(ct_mm, dct_um_per_yr))),
  "",
  sprintf("- **a** Study design. Cortex-wide mean thickness from %s scans of %s children (brain silhouettes are schematic, not data) enters one linear mixed model with a random intercept and age slope per child, giving each child's CT intercept and ΔCT slope. Polygenic scores built from published GWAS are tested against each trait (e), each trait against symptoms at 15–17 years adjusted for baseline symptoms (f, g), and the score → ΔCT → symptom path by mediation (h, i).", comma(n_scans), comma(n_kids)),
  sprintf("- **b** Age at scan by visit; %s children, %s scans.", comma(n_kids), comma(n_scans)),
  sprintf("- **c** Children by number of QC-passing scans: %s / %s / %s with 4 / 3 / 2 scans enter the slope model; %s children with a single scan pass QC but are not modelled (a slope needs at least two scans).",
          comma(spc$n_children[spc$n_scans == 4]), comma(spc$n_children[spc$n_scans == 3]),
          comma(spc$n_children[spc$n_scans == 2]), comma(spc$n_children[spc$n_scans == 1])),
  sprintf("- **d** Per-scan cortical mean thickness against age for 250 random children (grey) and three highlighted children (points, colour), chosen for well-separated slopes: a fast thinner (ΔCT %.1f µm/yr), a near-average child (%.1f) and a slow thinner (%.1f), four scans each spanning at least 5.5 years. Each coloured line is that child's fitted trajectory from the mixed model in the title, fixed effects plus the child's random intercept and slope (site effect omitted), so its height at 12.8 years is the child's CT and its slope is their ΔCT. Lines are shrunk toward the population trajectory, most for children with fewer or closely spaced scans.", hi_dct[1], hi_dct[3], hi_dct[2]),
  sprintf("- **e** PRS-CS polygenic score association with ΔCT (purple) and CT (green) on each row; β per SD of score with 95%% CI; filled = p < 0.05. Schizophrenia (2025 GWAS): ΔCT p = %s EUR / %s pooled, CT p = %s / %s. C+T and SBayesRC give the same ΔCT pattern; bipolar and ADHD scores are in SI. * p < 0.05, ** p < 0.01, *** p < 0.001 (uncorrected).",
          sp("SCZ25_EUR", "global_slope"), sp("SCZ25_META", "global_slope"),
          sp("SCZ25_EUR", "baseline_thickness"), sp("SCZ25_META", "baseline_thickness")),
  sprintf("- **f** Symptoms at ages ~15–17 (mean of the available waves 5–7) against the same single-LMM traits (ΔCT purple, CT green; each fitted alone; filled = p < 0.05; * p < 0.05, ** p < 0.01, *** p < 0.001). Model: y_late ~ trait + y_baseline + age_late + sex + site (OLS, family-clustered SEs), β in SD of y_late per SD of trait, so each estimate is the association with later symptoms beyond the child's baseline level (not a difference score or a symptom slope). CBCL parent report, log1p raw scale sums; p-factor is the first principal component of the eight syndrome scales. Traits are signed as in e (ΔCT negative = faster thinning), so β < 0 means more symptoms with faster thinning or with thinner cortex. n = %s. Faster thinning goes with more depressive (β %.3f, p = %s), externalising (p = %s) and general (p-factor, p = %s) symptoms; CT shows none of these. CBCL thought problems and youth-reported psychotic-like experiences (PQ-BC) are in SI. Exploratory, uncorrected.",
          comma(max(cb$n)), he("depress", "global_slope"),
          fmt_p(hp("depress", "global_slope")), fmt_p(hp("external", "global_slope")), fmt_p(hp("pfactor", "global_slope"))),
  sprintf("- **g** The strongest association in f at the child level: DSM-oriented depressive score at 15–17 years, residualised on the baseline score, age, sex and site (SD of the late score), against ΔCT (%s children; points, means and 95%% CIs in equal-width ΔCT bins of %.0f µm/yr (%s–%s children per bin; the %s children beyond the outer bin edges are omitted from the bins but included in the fit); purple, OLS fit on all children, p = %s with family-clustered SEs). The slope here (%.3f SD per SD of ΔCT) is slightly smaller than in f (%.3f) because ΔCT itself is not residualised on the covariates.",
          comma(hb$n_total[1]), hb$width[1], comma(min(hb$n)), comma(max(hb$n)), comma(hb$n_outside[1]), fmt_p(hb$fit_p[1]),
          hb$fit_slope_per_um[1] * csv_("slope_um_per_yr", "sd"), he("depress", "global_slope")),
  sprintf("- **h, i** Mediation of each polygenic score's association with later symptoms through ΔCT (pooled arm, %s and %s children), each for the outcome in f on which that score has its largest total effect (%s symptoms for schizophrenia, %s symptoms for depression). Product of coefficients with the covariates of f plus 10 genetic PCs: a, score → ΔCT (ΔCT in SD; negative = faster thinning); b, ΔCT → symptoms with the score in the model; c′, the direct effect; c, the total effect. The indirect effect a·b, with a 95%% CI from %d bootstrap replicates resampling families, is %+.4f for schizophrenia (%.1f%% of the total) and %+.4f for depression (%.1f%%). Because ΔCT is a shrunk per-child estimate, a·b is attenuated by its reliability R (mean %.2f; 0.17 / 0.32 / 0.38 for children with 2 / 3 / 4 scans). Dividing by R gives the reliability-corrected values shown in grey, %+.4f for schizophrenia (%.0f%% of the total) and %+.4f for depression (%.0f%%); this approximate correction treats R as known and the same for every child. a·b > 0 with both a and b negative: higher scores, faster thinning, more symptoms. * p < 0.05, ** p < 0.01, *** p < 0.001. All outcomes and both arms are in SI.",
          comma(med_g$n), comma(med_h$n), OUT_S[med_g$outcome], OUT_S[med_h$outcome], med_g$n_boot,
          med_g$indirect, 100 * med_g$prop, med_h$indirect, 100 * med_h$prop,
          REL, med_g$indirect / REL, 100 * med_g$prop / REL, med_h$indirect / REL, 100 * med_h$prop / REL),
  "", "Methods", sub("^\u2022  ", "- ", bullets), SI_NOTE)
writeLines(cap, CAP)
cat(CAP, "\n")
