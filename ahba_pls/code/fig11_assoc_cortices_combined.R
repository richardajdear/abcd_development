#!/usr/bin/env Rscript
# fig11_assoc_cortices_combined.R -- EXPLORATORY: regional maps of thinning vs
# later symptoms (p-factor, total problems; 30_cbcl_assoc_cortices.py) and vs
# SCZ 2025 / MDD polygenic scores (pooled arms, SBayesRC; 31_prs_assoc_cortices.py)
# on the 22 Glasser cortices, each tested against normative dCT, PLS2 and AHBA C3,
# with the parcel-vs-cortex comparison over the four model variants.
# Reads results/cbcl_assoc_cortex_{maps,corr}.tsv, prs_assoc_cortex_{maps,corr}.tsv, hcp_summary_maps.csv,
# data/hcp_polygons.csv, data/reference/hcp_cortices/hcp_parcel_systems.csv.
# Group-level only. Writes figures/fig_assoc_cortices_combined.png
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"; RES <- file.path(ROOT, "results")
BASE <- 9; GT <- (BASE - 1.5) / .pt; FADE <- 0.3
MAPV <- "relative_ct"                                           # the variant shown in a and b
MLAB <- c(absolute = "unadjusted", relative = "| global thinning", absolute_ct = "| baseline CT",
          relative_ct = "| global thinning + baseline CT")
OUT <- c(pfactor = "p-factor", totprob = "Total problems",
         SBayesRC_SCZ25_META = "SCZ 2025 PRS", SBayesRC_MDD_pooled = "MDD PRS")
KIND <- c(pfactor = "symptoms at 15\u201317 | baseline", totprob = "symptoms at 15\u201317 | baseline",
          SBayesRC_SCZ25_META = "SBayesRC, pooled arm", SBayesRC_MDD_pooled = "SBayesRC, pooled arm")
REFS <- c(dCT = "normative dCT\n(higher = slower thinning)", PLS2 = "PLS2", C3 = "AHBA C3")

# one long table per kind: `outcome` is a CBCL scale or a PRS arm
PM <- read.delim(file.path(RES, "prs_assoc_cortex_maps.tsv"))
M  <- bind_rows(read.delim(file.path(RES, "cbcl_assoc_cortex_maps.tsv")) |> select(outcome, map, cortex_id, r),
                PM |> filter(level == "cortex") |> transmute(outcome = arm, map, cortex_id = as.integer(region), r)) |>
  filter(outcome %in% names(OUT))
CR <- bind_rows(read.delim(file.path(RES, "cbcl_assoc_cortex_corr.tsv")) |> select(outcome, map, level, reference, rho, p_spin, p_perm),
                read.delim(file.path(RES, "prs_assoc_cortex_corr.tsv")) |>
                  transmute(outcome = arm, map, level, reference, rho, p_spin, p_perm)) |>
  filter(outcome %in% names(OUT))
G  <- read.csv(file.path(RES, "hcp_summary_maps.csv")) |> mutate(key = tolower(label))
SY <- read.csv(file.path(ROOT, "data", "reference", "hcp_cortices", "hcp_parcel_systems.csv")) |> mutate(key = tolower(label))
poly <- read.csv(file.path(ROOT, "data", "hcp_polygons.csv")) |> filter(view %in% c("lateral", "medial")) |> mutate(key = tolower(label))
stopifnot(n_distinct(M$cortex_id) == 22, setequal(unique(M$outcome), names(OUT)),
          nrow(CR) == length(OUT) * 4 * 2 * 3)

# reference maps per cortex = mean over the cortex's parcels (PLS2, C3: AHBA-covered only), as in script 30
RC <- G |> inner_join(SY |> select(key, cortex_id), by = "key") |>
  group_by(cortex_id) |> summarise(dCT = mean(dCT, na.rm = TRUE), PLS2 = mean(PLS2, na.rm = TRUE),
                                   C3 = mean(C3, na.rm = TRUE), .groups = "drop")
cid <- SY |> select(key, cortex_id)

base_theme <- theme_classic(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0),
        strip.background = element_blank(), strip.text = element_text(size = BASE - 0.5, face = "bold"),
        axis.text = element_text(size = BASE - 1.5, colour = "grey20"), axis.title = element_text(size = BASE - 1),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3))

# ------------------------------------------------------------ a: cortex maps ----
brain <- function(vals, title, fill) {
  d <- poly |> left_join(cid, by = "key") |> left_join(vals, by = "cortex_id")
  ggplot(d, aes(x, y, group = interaction(view, label, group, subgroup), fill = v)) +
    geom_polygon(colour = "grey35", linewidth = 0.05) + coord_fixed(expand = FALSE) + fill +
    labs(title = title) + theme_void(base_size = BASE) +
    theme(plot.title = element_text(size = BASE - 0.5, face = "bold", hjust = 0.5),
          legend.position = "bottom", legend.text = element_text(size = BASE - 2),
          legend.key.height = unit(4, "pt"), legend.key.width = unit(24, "pt"),
          legend.margin = margin(0, 0, 0, 0), plot.margin = margin(1, 3, 2, 3))
}
div <- function(v) { m <- max(abs(v), na.rm = TRUE); b <- signif(0.75 * m, 1)
  scale_fill_distiller(palette = "RdBu", limits = c(-m, m), na.value = "grey82", breaks = c(-b, 0, b), name = NULL) }
am <- lapply(names(OUT), function(o) {
  v <- M |> filter(outcome == o, map == MAPV) |> transmute(cortex_id, v = r)
  brain(v, sprintf("%s\n(partial r)", OUT[[o]]), div(v$v))
})
dct_scale <- scale_fill_gradientn(colours = c("#67000d", "#cb181d", "#fb6a4a", "#fcbba1", "white"),
                                  limits = c(min(RC$dCT), 0), na.value = "grey82", breaks = scales::pretty_breaks(3), name = NULL)
ar <- list(brain(RC |> transmute(cortex_id, v = dCT), "Normative dCT (mm/yr)", dct_scale),
           brain(RC |> transmute(cortex_id, v = PLS2), "PLS2", div(RC$PLS2)),
           brain(RC |> transmute(cortex_id, v = C3), "AHBA C3", div(RC$C3)))
SYM <- c("pfactor", "totprob"); PRS <- c("SBayesRC_SCZ25_META", "SBayesRC_MDD_pooled")
names(am) <- names(OUT)
hdr <- function(txt) ggplot() + annotate("text", 0, 0, label = txt, hjust = 0, size = (BASE + 0.5) / .pt, fontface = "bold") +
  xlim(0, 1) + theme_void() + theme(plot.margin = margin(2, 0, 0, 2))
brow <- function(ps, title) wrap_elements(full = wrap_plots(ps, nrow = 1) + plot_annotation(title = title,
  theme = theme(plot.title = element_text(size = BASE, face = "bold"))))

# ------------------------------------------------ scatters vs reference -------
scat <- function(os, title, ylab) {
  sb <- M |> filter(map == MAPV, outcome %in% os) |> select(outcome, cortex_id, r) |>
    inner_join(RC |> pivot_longer(-cortex_id, names_to = "reference", values_to = "x"), by = "cortex_id",
               relationship = "many-to-many") |>
    mutate(outcome = factor(OUT[outcome], OUT[os]), reference = factor(REFS[reference], REFS))
  lb <- CR |> filter(level == "cortex", map == MAPV, outcome %in% os) |>
    mutate(outcome = factor(OUT[outcome], OUT[os]), reference = factor(REFS[reference], REFS),
           sig = p_spin < 0.05 & p_perm < 0.05,
           txt = sprintf("rho = %.2f\np_spin %.3f, p_perm %.3f", rho, p_spin, p_perm))
  stopifnot(nrow(lb) == length(os) * length(REFS))
  sb <- sb |> left_join(lb |> select(outcome, reference, sig), by = c("outcome", "reference"))
  ggplot(sb, aes(x, r, alpha = sig)) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_point(size = 1.7, colour = "grey15", stroke = 0) +
    geom_smooth(data = \(z) filter(z, sig), method = "lm", formula = y ~ x, se = FALSE, colour = "grey10", linewidth = 0.5) +
    geom_smooth(data = \(z) filter(z, !sig), method = "lm", formula = y ~ x, se = FALSE, colour = alpha("grey10", FADE), linewidth = 0.5) +
    geom_text(data = lb, aes(-Inf, Inf, label = txt), hjust = -0.05, vjust = 1.15, size = GT - 0.2, fontface = "bold", lineheight = 0.9) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.42))) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    facet_grid(outcome ~ reference, scales = "free", switch = "y") +
    labs(x = "reference map, cortex mean", y = ylab, title = title) +
    base_theme + theme(strip.placement = "outside", panel.spacing = unit(7, "pt"), plot.title.position = "plot")
}

# ---------------------------------------- parcel vs cortex, all variants ------
dumb <- function(os, title) {
  cc <- CR |> filter(outcome %in% os) |>
    mutate(outcome = factor(OUT[outcome], OUT[os]), reference = factor(sub("\n.*", "", REFS[reference]), sub("\n.*", "", REFS)),
           map = factor(MLAB[map], rev(MLAB)), level = factor(level, c("parcel", "cortex"), c("179 parcels", "22 cortices")),
           sig = p_spin < 0.05 & p_perm < 0.05)
  ggplot(cc, aes(rho, map, colour = level, shape = sig)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_line(aes(group = map), colour = "grey75", linewidth = 0.4) +
    geom_point(size = 2.1, stroke = 0.8) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p_spin and p_perm < 0.05", `FALSE` = "n.s."), name = NULL, drop = FALSE) +
    scale_colour_manual(values = c("179 parcels" = "grey55", "22 cortices" = "#b2182b"), name = NULL) +
    scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5)) +
    facet_grid(outcome ~ reference) +
    labs(x = "Spearman \u03c1 with the reference map", y = NULL, title = title) +
    base_theme + theme(legend.position = "bottom", panel.spacing = unit(7, "pt"), plot.title.position = "plot",
                       panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.25))
}

pa  <- (plot_spacer() | brow(ar, "a  Reference maps, averaged into the 22 Glasser cortices") | plot_spacer()) + plot_layout(widths = c(0.5, 3, 0.5))
pbL <- brow(am[SYM], sprintf("b  Thinning vs later symptoms (%s)", MLAB[[MAPV]]))
pbR <- brow(am[PRS], sprintf("c  Thinning vs polygenic scores (%s)", MLAB[[MAPV]]))
pdL <- scat(SYM, "d  Symptom maps vs the reference maps (22 cortices)", "partial r: symptoms at 15\u201317 | baseline")
pdR <- scat(PRS, "e  PRS maps vs the reference maps (22 cortices)", "partial r: SBayesRC score, pooled arm")
pfL <- dumb(SYM, "f  Symptoms: 179 parcels vs 22 cortices, all four models")
pfR <- dumb(PRS, "g  PRS: 179 parcels vs 22 cortices, all four models")
colL <- (hdr("Symptoms at 15\u201317 given baseline (CBCL)") / pbL / pdL / pfL) + plot_layout(heights = c(0.06, 0.42, 1, 0.95))
colR <- (hdr("SCZ 2025 and MDD polygenic scores (SBayesRC, pooled arm)") / pbR / pdR / pfR) + plot_layout(heights = c(0.06, 0.42, 1, 0.95))

hit <- CR |> filter(p_spin < 0.05, p_perm < 0.05) |>
  mutate(txt = sprintf("%s %s vs %s at %s (\u03c1 = %.2f)", OUT[outcome], MLAB[map], reference,
                       ifelse(level == "cortex", "22 cortices", "179 parcels"), rho))
hit <- hit |> mutate(sym = outcome %in% c("pfactor", "totprob")); sh <- hit
stopifnot(all(hit$reference[hit$sym] == "dCT"), all(hit$rho[hit$sym] > 0))
sub_txt <- paste0(
  "Symptoms: per region, CBCL at 15\u201317 ~ thinning + baseline CBCL + sex + site + ages + scans (n \u2248 8,200). PRS: per region, thinning ~ SBayesRC score + sex + baseline age + PC1\u201310 (pooled arms, score standardised within ancestry, n = 8,596).\n",
  "Both: [+ global thinning] [+ that region's baseline CT]; 22 cortices = Glasser 2016, bilateral mean. Nulls: parcel-level spin re-averaged into regions, and permutation (outcome residuals / score).\n",
  sprintf("Tests passing both nulls: symptoms %d of %d, all vs dCT (\u03c1 %.2f to %.2f, positive = extra thinning where cortex normally thins least); ",
          sum(sh$sym), nrow(CR) / 2, min(hit$rho[hit$sym]), max(hit$rho[hit$sym])),
  sprintf("PRS %d of %d: %s.", sum(!hit$sym), nrow(CR) / 2, paste(hit$txt[!hit$sym], collapse = "; ")))

fig <- (wrap_elements(full = pa) / ((colL | colR) + plot_layout(guides = "collect"))) + plot_layout(heights = c(0.18, 1)) &
  theme(legend.position = "bottom")
fig <- fig + plot_annotation(title = "Regional thinning linked to later symptoms follows normative thinning; thinning linked to SCZ / MDD polygenic scores does not clearly follow dCT, PLS2 or AHBA C3",
                  subtitle = sub_txt,
                  theme = theme(plot.title = element_text(size = BASE + 2, face = "bold"),
                                plot.subtitle = element_text(size = BASE - 1, colour = "grey25", lineheight = 1.15)))
ggsave(file.path(ROOT, "figures", "fig_assoc_cortices_combined.png"), fig, width = 15, height = 14, dpi = 300, bg = "white")
cat("wrote fig_assoc_cortices_combined.png\n")
