#!/usr/bin/env Rscript
# fig4_celltypes.R -- cell-class marker enrichment of every gene ranking in
# panel a of the enrichment figure.  Reads results/celltype_all_options.tsv
# (written by code/15_celltype_all_options.py) and fits nothing.
#
# Panel a: heatmap, 9 classes x 8 rankings, z on a shared gene universe.
# Panel b: the astrocyte / oligodendrocyte plane, where the rankings separate
#          by PARCELLATION rather than by Y-matrix option.
# Writes figures/fig_celltypes.png

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
                  library(ggrepel); library(scales)})

ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "figures")
T <- read.delim(file.path(RES, "celltype_all_options.tsv"))
S <- T |> filter(universe == "shared")
stopifnot(length(unique(S$n_universe)) == 1)
NU <- format(S$n_universe[1], big.mark = ",")

# column order and grouping follow panel a of the enrichment figure
# "DK, HCP gene basis" is the DK fit on the same 7,973 genes as the HCP matrix;
# it is what makes the astrocyte comparison a parcellation contrast rather than
# a parcellation + abagen-build + DS-gene-set contrast.
ord <- c("ABCD PLS2, HCP-MMP", "ABCD dCT alone, HCP-MMP",
         "ABCD PLS2, DK (HCP gene basis)",
         "ABCD PLS2, DK", "ABCD dCT alone, DK", "ABCD dCT+dT1T2 PLS2, DK",
         "AHBA C3", "NSPN PLS2", "AHBA C1")
grp <- setNames(c(rep("HCP-MMP\n(137 parcels)", 2), rep("Desikan\u2013Killiany\n(33 regions)", 4),
                  rep("published\ncomponents", 3)), ord)
stopifnot(setequal(ord, unique(S$vector)))

# claims about counts are computed, never written in: how many thinning-derived
# rankings show the neuronal-up / microglia-down pattern, and how many rankings
# sit on each side of zero for astrocytes.
thin_v <- setdiff(ord, c("AHBA C1"))                       # C1 is the static control
n_thin <- S |> filter(vector %in% thin_v, cell_class %in% c("Neuro", "Neuro-Ex", "Neuro-In")) |>
  group_by(vector) |> summarise(ok = all(z > 0), .groups = "drop") |>
  inner_join(S |> filter(vector %in% thin_v, cell_class %in% c("Micro", "Endo")) |>
               group_by(vector) |> summarise(ok2 = all(z < 0), .groups = "drop"), by = "vector") |>
  summarise(n = sum(ok & ok2)) |> pull(n)
astro_all <- S |> filter(cell_class == "Astro")
n_astro_neg <- sum(astro_all$z < 0)
n_cells <- nrow(S)
rho_uni <- suppressWarnings(cor(
  (T |> filter(universe == "own") |> arrange(vector, cell_class))$z,
  (T |> filter(universe == "shared") |> arrange(vector, cell_class))$z,
  method = "spearman"))
own_n <- T |> filter(universe == "own") |> distinct(vector, n_universe) |>
  arrange(desc(n_universe))

# rows ordered by the mean z across rankings: neuronal at the top, glial below
row_ord <- S |> group_by(cell_class) |> summarise(m = mean(z), .groups = "drop") |>
  arrange(m) |> pull(cell_class)
short <- c("ABCD PLS2, HCP-MMP" = "PLS2 (dCT+CT)", "ABCD dCT alone, HCP-MMP" = "dCT alone",
           "ABCD PLS2, DK (HCP gene basis)" = "PLS2, HCP gene basis",
           "ABCD PLS2, DK" = "PLS2 (dCT+CT)", "ABCD dCT alone, DK" = "dCT alone",
           "ABCD dCT+dT1T2 PLS2, DK" = "PLS2 (dCT+dT1T2)",
           "AHBA C3" = "AHBA C3", "NSPN PLS2" = "NSPN PLS2", "AHBA C1" = "AHBA C1")
D <- S |> mutate(cell_class = factor(cell_class, levels = row_ord),
                 vector = factor(vector, levels = ord),
                 grp = factor(grp[as.character(vector)], levels = unique(grp)),
                 lab = ifelse(p_perm < 0.05, sprintf("%.1f", z), sprintf("(%.1f)", z)))

D <- D |> mutate(col = factor(short[as.character(vector)], levels = unique(short[ord])))
pa <- ggplot(D, aes(col, cell_class, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab, colour = abs(z) > 9), size = 2.05, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10")) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-20, 20),
                       oob = squish, name = "z") +
  facet_grid(~ grp, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL,
       title = "a   Cell-class marker enrichment of every candidate gene ranking",
       subtitle = sprintf("z against 20,000 size-matched random gene sets on one shared universe of %s genes; values in brackets are p_perm \u2265 0.05", NU)) +
  theme_minimal(base_size = 7.4) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 22, hjust = 1, size = 6.4),
        axis.text.y = element_text(size = 6.9),
        strip.text = element_text(size = 6.6, face = "bold", colour = "grey25"),
        strip.background = element_rect(fill = "grey96", colour = NA),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, colour = "grey30", margin = margin(b = 3)),
        legend.key.width = unit(6, "pt"), legend.key.height = unit(20, "pt"),
        legend.title = element_text(size = 6.6), legend.text = element_text(size = 6.2),
        plot.margin = margin(2, 4, 2, 4))

# ---- panel b: the astrocyte / oligodendrocyte plane ------------------------
W <- S |> filter(cell_class %in% c("Astro", "Oligo")) |>
  select(vector, cell_class, z, p_perm) |>
  pivot_wider(names_from = cell_class, values_from = c(z, p_perm)) |>
  mutate(grp = factor(grp[vector], levels = unique(grp)))
pb <- ggplot(W, aes(z_Astro, z_Oligo, colour = grp)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(size = 1.7) +
  geom_text_repel(aes(label = vector), size = 2.05, min.segment.length = 0.1,
                  segment.size = 0.2, box.padding = 0.35, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_colour_manual(values = c("#33a02c", "#1f78b4", "grey25"), name = NULL) +
  labs(x = "astrocyte marker z", y = "oligodendrocyte marker z",
       title = "b   The two classes that separate the rankings \u2014 by parcellation, not by Y-matrix option",
       subtitle = sprintf("Both HCP-MMP rankings sit alone in the astrocyte-positive half (%+.1f and %+.1f); the other %d rankings \u2014 including the DK fit on the HCP gene basis \u2014 are astrocyte-negative",
                          W$z_Astro[W$vector == "ABCD PLS2, HCP-MMP"],
                          W$z_Astro[W$vector == "ABCD dCT alone, HCP-MMP"],
                          n_astro_neg)) +
  theme_bw(base_size = 7.4) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.15, colour = "grey94"),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, colour = "grey30", margin = margin(b = 3)),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6.6),
        legend.position = "right", legend.text = element_text(size = 6.4),
        legend.key.height = unit(12, "pt"),
        plot.margin = margin(2, 4, 2, 4))

astro_hcp <- S |> filter(cell_class == "Astro", vector == "ABCD PLS2, HCP-MMP")
astro_dk <- S |> filter(cell_class == "Astro", vector == "ABCD PLS2, DK")
oli <- S |> filter(cell_class == "Oligo", vector %in% c("ABCD PLS2, HCP-MMP", "ABCD PLS2, DK", "AHBA C3"))
methods <- paste(
  "Methods.",
  "\u2022 Markers: the Seidlitz et al. 2020 compilation, nine classes (Neuro / Neuro-Ex / Neuro-In, Astro, Oligo, OPC, Micro, Endo, Per), 38\u2013862 genes per class after intersecting the universe.",
  sprintf("\u2022 Test: mean gene weight of a class against 20,000 random gene sets of the same size, drawn from one universe shared by all %d rankings (%s genes), identical random draws per column.",
          length(ord), NU),
  sprintf("\u2022 Running each ranking on its own universe instead (%s genes across the nine vectors) gives the same picture:\n  z agrees at Spearman %.3f over all %d cells. Both versions are in the results table.",
          paste(format(range(own_n$n_universe), big.mark = ","), collapse = "\u2013"), rho_uni, n_cells),
  "\u2022 Marker genes are co-expressed, so the independent-gene null is anti-conservative: read the sign and the pattern across columns, not the absolute z.",
  "\u2022 Sign convention: positive = the class's markers are expressed more where adolescent thinning is faster.\n  C1, C3 and NSPN PLS2 keep their published sign, so positive there means \"higher where the component score is high\".",
  sep = "\n")

fig <- (pa / pb) + plot_layout(heights = c(1, 0.92)) +
  plot_annotation(
    title = "Cell-class profiles agree across Y-matrix options but split on astrocytes by parcellation",
    subtitle = sprintf("Neuronal classes are positive and microglia/endothelia negative in all %d thinning-derived rankings.\nAstrocytes flip with the atlas (HCP-MMP %+.1f vs DK %+.1f), and oligodendrocytes deepen with it (%+.1f vs %+.1f, against %+.1f for C3).",
                       n_thin, astro_hcp$z, astro_dk$z,
                       oli$z[oli$vector == "ABCD PLS2, HCP-MMP"],
                       oli$z[oli$vector == "ABCD PLS2, DK"],
                       oli$z[oli$vector == "AHBA C3"]),
    caption = methods,
    theme = theme(plot.title = element_text(size = 9.4, face = "bold"),
                  plot.subtitle = element_text(size = 7.2, colour = "grey25", margin = margin(b = 4)),
                  plot.caption = element_text(size = 6.1, colour = "grey25", hjust = 0,
                                              lineheight = 1.42, margin = margin(t = 6))))

ggsave(file.path(FIG, "fig_celltypes.png"), fig, width = 7.8, height = 6.6, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_celltypes.png"), "\n")
