#!/usr/bin/env Rscript
# fig9_cbcl_assoc_cortices.R -- EXPLORATORY: thinning-symptom association maps
# on the 22 Glasser cortices (30_cbcl_assoc_cortices.py) for p-factor and total
# problems, against PLS2, AHBA C3 and normative dCT, and the parcel-vs-cortex
# comparison. Reads results/cbcl_assoc_cortex_{maps,corr}.tsv, hcp_summary_maps.csv,
# data/hcp_polygons.csv, data/reference/hcp_cortices/hcp_parcel_systems.csv.
# Group-level only. Writes figures/fig_cbcl_assoc_cortices.png
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"; RES <- file.path(ROOT, "results")
BASE <- 9; GT <- (BASE - 1.5) / .pt; FADE <- 0.3
MAPV <- "relative_ct"                                           # the variant shown in a and b
MLAB <- c(absolute = "unadjusted", relative = "| global thinning", absolute_ct = "| baseline CT",
          relative_ct = "| global thinning + baseline CT")
OUT <- c(pfactor = "p-factor", totprob = "Total problems")
REFS <- c(dCT = "normative dCT\n(higher = slower thinning)", PLS2 = "PLS2", C3 = "AHBA C3")

M  <- read.delim(file.path(RES, "cbcl_assoc_cortex_maps.tsv"))
CR <- read.delim(file.path(RES, "cbcl_assoc_cortex_corr.tsv"))
G  <- read.csv(file.path(RES, "hcp_summary_maps.csv")) |> mutate(key = tolower(label))
SY <- read.csv(file.path(ROOT, "data", "reference", "hcp_cortices", "hcp_parcel_systems.csv")) |> mutate(key = tolower(label))
poly <- read.csv(file.path(ROOT, "data", "hcp_polygons.csv")) |> filter(view %in% c("lateral", "medial")) |> mutate(key = tolower(label))
stopifnot(n_distinct(M$cortex_id) == 22, all(names(OUT) %in% M$outcome))

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
  brain(v, sprintf("%s: partial r", OUT[[o]]), div(v$v))
})
dct_scale <- scale_fill_gradientn(colours = c("#67000d", "#cb181d", "#fb6a4a", "#fcbba1", "white"),
                                  limits = c(min(RC$dCT), 0), na.value = "grey82", breaks = scales::pretty_breaks(3), name = NULL)
ar <- list(brain(RC |> transmute(cortex_id, v = dCT), "Normative dCT (mm/yr)", dct_scale),
           brain(RC |> transmute(cortex_id, v = PLS2), "PLS2", div(RC$PLS2)),
           brain(RC |> transmute(cortex_id, v = C3), "AHBA C3", div(RC$C3)))
pa <- wrap_elements(full = wrap_plots(c(am, ar), nrow = 1) +
  plot_annotation(title = sprintf("a  Association maps on the 22 cortices (%s) and the reference maps averaged the same way", MLAB[[MAPV]]),
                  theme = theme(plot.title = element_text(size = BASE, face = "bold"))))

# ------------------------------------------------ b: cortex map vs reference ----
sb <- M |> filter(map == MAPV, outcome %in% names(OUT)) |> select(outcome, cortex_id, r) |>
  inner_join(RC |> pivot_longer(-cortex_id, names_to = "reference", values_to = "x"), by = "cortex_id",
             relationship = "many-to-many") |>
  mutate(outcome = factor(OUT[outcome], OUT), reference = factor(REFS[reference], REFS))
lb <- CR |> filter(level == "cortex", map == MAPV, outcome %in% names(OUT)) |>
  mutate(outcome = factor(OUT[outcome], OUT), reference = factor(REFS[reference], REFS),
         sig = p_spin < 0.05 & p_perm < 0.05,
         txt = sprintf("rho = %.2f\np_spin %.3f, p_perm %.3f", rho, p_spin, p_perm))
stopifnot(nrow(lb) == length(OUT) * length(REFS))
sb <- sb |> left_join(lb |> select(outcome, reference, sig), by = c("outcome", "reference"))
pb <- ggplot(sb, aes(x, r, alpha = sig)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(size = 1.9, colour = "grey15", stroke = 0) +
  geom_smooth(data = \(z) filter(z, sig), method = "lm", formula = y ~ x, se = FALSE, colour = "grey10", linewidth = 0.5) +
  geom_smooth(data = \(z) filter(z, !sig), method = "lm", formula = y ~ x, se = FALSE, colour = alpha("grey10", FADE), linewidth = 0.5) +
  geom_text(data = lb, aes(-Inf, Inf, label = txt), hjust = -0.05, vjust = 1.15, size = GT, fontface = "bold", lineheight = 0.9) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.38))) +
  facet_grid(outcome ~ reference, scales = "free", switch = "y") +
  labs(x = "reference map, cortex mean", y = "partial r of later symptoms with thinning",
       title = "b  22 cortices: each symptom map against the three reference maps (faded = not both p < 0.05)") +
  base_theme + theme(strip.placement = "outside", panel.spacing = unit(8, "pt"))

# ---------------------------------------- c: parcel vs cortex, all variants ----
cc <- CR |> filter(outcome %in% names(OUT)) |>
  mutate(outcome = factor(OUT[outcome], OUT), reference = factor(sub("\n.*", "", REFS[reference]), sub("\n.*", "", REFS)),
         map = factor(MLAB[map], rev(MLAB)), level = factor(level, c("parcel", "cortex"),
                                                             c("179 parcels", "22 cortices")),
         sig = p_spin < 0.05 & p_perm < 0.05)
pc <- ggplot(cc, aes(rho, map, colour = level, shape = sig)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(aes(group = map), colour = "grey75", linewidth = 0.4) +
  geom_point(size = 2.3, stroke = 0.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p_spin and p_perm < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_colour_manual(values = c("179 parcels" = "grey55", "22 cortices" = "#b2182b"), name = NULL) +
  facet_grid(outcome ~ reference) +
  labs(x = "Spearman \u03c1 of the association map with the reference map", y = NULL,
       title = "c  Does averaging into 22 cortices sharpen the pattern? (all four model variants)") +
  base_theme + theme(legend.position = "bottom", panel.spacing = unit(8, "pt"),
                     panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.25))

# -------------------------------------------------------------- assemble ------
g <- function(o, lev, ref, m = MAPV) CR[CR$outcome == o & CR$level == lev & CR$reference == ref & CR$map == m, ]
sub_txt <- sprintf(paste0(
  "Per region: CBCL at 15\u201317 ~ thinning + baseline CBCL + sex + site + ages + scans [+ global thinning] [+ baseline CT]; n \u2248 8,200. Cortices = Glasser 2016 (22 per hemisphere, bilateral mean).\n",
  "With global thinning and baseline CT adjusted, averaging into cortices raises the correlation with normative dCT (p-factor \u03c1 %.2f \u2192 %.2f, total problems %.2f \u2192 %.2f) but not the evidence\n",
  "(p_spin %.3f \u2192 %.3f and %.3f \u2192 %.3f: 22 regions leave far fewer degrees of freedom). PLS2 and AHBA C3 stay null at both resolutions."),
  g("pfactor", "parcel", "dCT")$rho, g("pfactor", "cortex", "dCT")$rho, g("totprob", "parcel", "dCT")$rho, g("totprob", "cortex", "dCT")$rho,
  g("pfactor", "parcel", "dCT")$p_spin, g("pfactor", "cortex", "dCT")$p_spin, g("totprob", "parcel", "dCT")$p_spin, g("totprob", "cortex", "dCT")$p_spin)
fig <- (pa / pb / pc) + plot_layout(heights = c(0.5, 1.05, 0.9)) +
  plot_annotation(title = "At the 22-cortex level, symptom-linked thinning tracks where cortex normally thins least \u2014 not PLS2 or AHBA C3",
                  subtitle = sub_txt,
                  theme = theme(plot.title = element_text(size = BASE + 2, face = "bold"),
                                plot.subtitle = element_text(size = BASE - 1, colour = "grey25", lineheight = 1.15)))
ggsave(file.path(ROOT, "figures", "fig_cbcl_assoc_cortices.png"), fig, width = 12, height = 11, dpi = 300, bg = "white")
cat("wrote fig_cbcl_assoc_cortices.png\n")
