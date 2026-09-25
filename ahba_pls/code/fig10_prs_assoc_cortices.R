#!/usr/bin/env Rscript
# fig10_prs_assoc_cortices.R -- EXPLORATORY: regional maps of SCZ / MDD polygenic
# score vs adolescent thinning (31_prs_assoc_cortices.py), PRS-CS and SBayesRC,
# and every map-vs-reference test (normative dCT, PLS2, AHBA C3) at 22 cortices
# and 179 parcels. Reads results/prs_assoc_cortex_{maps,corr}.tsv,
# data/hcp_polygons.csv, data/reference/hcp_cortices/hcp_parcel_systems.csv.
# Group-level only. Writes figures/fig_prs_assoc_cortices.png
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"; RES <- file.path(ROOT, "results")
BASE <- 9; GT <- (BASE - 2) / .pt; FADE <- 0.3
M  <- read.delim(file.path(RES, "prs_assoc_cortex_maps.tsv"))
CR <- read.delim(file.path(RES, "prs_assoc_cortex_corr.tsv"))
SY <- read.csv(file.path(ROOT, "data", "reference", "hcp_cortices", "hcp_parcel_systems.csv")) |> mutate(key = tolower(label))
poly <- read.csv(file.path(ROOT, "data", "hcp_polygons.csv")) |> filter(view %in% c("lateral", "medial")) |> mutate(key = tolower(label))

ARM <- c(SCZ25_META = "SCZ 2025, pooled", SCZ25_EUR = "SCZ 2025, EUR", MDD_pooled = "MDD, pooled", MDD_eur = "MDD, EUR")
METH <- c(PRSCS = "PRS-CS", SBayesRC = "SBayesRC")
MLAB <- c(absolute = "unadjusted", relative = "| global thinning", absolute_ct = "| baseline CT",
          relative_ct = "| global thinning\n+ baseline CT")
REF <- c(dCT = "dCT", PLS2 = "PLS2", C3 = "C3")
split_arm <- function(d) d |> mutate(method = sub("_.*", "", arm), trait = sub("^[^_]+_", "", arm)) |>
  filter(trait %in% names(ARM)) |>
  mutate(method = factor(METH[method], METH), trait = factor(ARM[trait], ARM))
M <- split_arm(M); CR <- split_arm(CR)
stopifnot(nlevels(droplevels(CR$trait)) == 4, nlevels(droplevels(CR$method)) == 2)

# ------------------------------------------------ a: cortex maps (unadjusted) ----
MAPV <- "absolute"
cm <- M |> filter(level == "cortex", map == MAPV) |> transmute(method, trait, cortex_id = as.integer(region), r)
lim <- max(abs(cm$r)) * c(-1, 1)
pd_ <- poly |> left_join(SY |> select(key, cortex_id), by = "key") |>
  inner_join(cm, by = "cortex_id", relationship = "many-to-many")
pa <- ggplot(pd_, aes(x, y, group = interaction(view, label, group, subgroup), fill = r)) +
  geom_polygon(colour = "grey35", linewidth = 0.04) + coord_fixed(expand = FALSE) +
  scale_fill_distiller(palette = "RdBu", limits = lim, breaks = c(-0.04, -0.02, 0, 0.02, 0.04), name = "partial r") +
  facet_grid(method ~ trait, switch = "y") +
  labs(title = "a  Score vs each child's thinning rate, averaged in each of the 22 Glasser cortices (unadjusted; red = higher score, faster thinning)") +
  theme_void(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold", margin = margin(b = 4)),
        strip.text = element_text(size = BASE - 0.5, face = "bold", margin = margin(2, 2, 2, 2)),
        strip.text.y.left = element_text(angle = 90),
        legend.position = "right", legend.key.width = unit(6, "pt"), legend.key.height = unit(26, "pt"),
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 1.5),
        panel.spacing = unit(6, "pt"))

# --------------------------------------------- b: every map-vs-reference test ----
tb <- CR |> mutate(ref = factor(REF[reference], REF), map = factor(MLAB[map], MLAB),
                   level = factor(level, c("cortex", "parcel"), c("22", "179")),
                   col = interaction(ref, level, sep = "\n", lex.order = TRUE),
                   row = factor(paste(trait, method, sep = " \u00b7 "),
                                rev(as.vector(outer(ARM, METH, paste, sep = " \u00b7 ")))),
                   sig = p_spin < 0.05 & p_perm < 0.05, txt = sprintf("%.2f", rho))
stopifnot(nrow(tb) == 8 * 4 * 3 * 2)
pb <- ggplot(tb, aes(col, row)) +
  geom_tile(aes(fill = rho, alpha = sig), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt, alpha = sig, fontface = ifelse(sig, "bold", "plain")), size = GT) +
  scale_fill_distiller(palette = "RdBu", limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5), name = "Spearman \u03c1") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
  facet_wrap(~map, nrow = 1) +
  labs(x = "reference map (normative dCT, PLS2, AHBA C3) \u00b7 resolution (22 cortices, 179 parcels)", y = NULL,
       title = "b  Does the PRS map follow normative dCT, PLS2 or AHBA C3? (bold = p_spin and p_perm < 0.05)") +
  theme_minimal(base_size = BASE) +
  theme(plot.title = element_text(size = BASE, face = "bold"), plot.title.position = "plot",
        panel.grid = element_blank(), strip.text = element_text(size = BASE - 0.5, face = "bold"),
        axis.text.x = element_text(size = BASE - 1.5, lineheight = 0.9), axis.title.x = element_text(size = BASE - 1, margin = margin(t = 4)), axis.text.y = element_text(size = BASE - 1),
        legend.position = "right", legend.key.width = unit(6, "pt"), legend.key.height = unit(26, "pt"),
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 1.5))

nb <- sum(tb$sig)
# a test "replicates" when it passes both nulls under BOTH scoring methods
repl <- tb |> group_by(trait, map, level, ref) |> filter(all(sig), n() == 2) |>
  summarise(rho = paste(sprintf("%.2f", rho[order(method)]), collapse = " / "), .groups = "drop")
rtxt <- if (nrow(repl)) paste(sprintf("%s %s vs %s, %s regions (%s)", repl$trait, gsub("\n", " ", repl$map), repl$ref,
                                      repl$level, repl$rho), collapse = "; ") else "none"
fig <- (pa / pb) + plot_layout(heights = c(0.42, 1)) +
  plot_annotation(
    title = "Polygenic-score maps of adolescent thinning: no robust alignment with PLS2 or AHBA C3; MDD (EUR) tracks normative thinning once baseline CT is adjusted",
    subtitle = sprintf(paste0(
      "Per region: thinning ~ score + sex + baseline age + PC1\u201310 [+ global thinning] [+ that region's baseline CT] (the fixed part of the Figure 1 PRS model; OLS). Pooled arms: n = 8,596, score standardised within ancestry;\n",
      "EUR arms: n = 4,308, raw score. Nulls: spin of the parcel-level reference re-averaged into cortices (5,000); permutation of the residualised score (1,000). %d of %d tests pass both.\n",
      "Passing under both PRS-CS and SBayesRC (\u03c1 PRS-CS / SBayesRC): %s.\nNegative \u03c1 with dCT = the score accelerates thinning where cortex normally thins fastest."),
      nb, nrow(tb), rtxt),
    theme = theme(plot.title = element_text(size = BASE + 2, face = "bold"),
                  plot.subtitle = element_text(size = BASE - 1, colour = "grey25", lineheight = 1.15)))
ggsave(file.path(ROOT, "figures", "fig_prs_assoc_cortices.png"), fig, width = 14, height = 9.8, dpi = 300, bg = "white")
cat("wrote fig_prs_assoc_cortices.png\n")
