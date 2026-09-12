#!/usr/bin/env Rscript
# fig1_lead_signature.R -- Figure 1 in ggplot2 + patchwork + ggseg.
#
# Reads only saved tables (nothing is recomputed here); every statistic printed
# on the figure is read from the table at plot time.
#   data/y_maps_bilateral_34.csv              Y maps
#   data/dk_polygons.csv                      DK polygons (from ggseg, cached)
#   results/lead_signature_scores.csv         region scores + references
#   results/lead_signature_weights.tsv        gene weights (thinning-oriented Z)
#   results/lead_celltype_enrichment.tsv      cell-class marker enrichment
#   results/concordance_scores.tsv            spin rho / p for the score maps
#   results/ymaps_vs_nspn.tsv                 spin rho / p for the raw Y maps
#   results/pls_components.tsv                component stats
#   results/lead_overlap_with_references.tsv  decile overlaps
#   data/reference/*                          NSPN and C1-C3 references
# Writes figures/fig_lead_signature.png

suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(ggrepel); library(scales)
})

ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), ".."))
if (is.na(ROOT)) ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results"); DATA <- file.path(ROOT, "data")
REF <- file.path(DATA, "reference"); FIG <- file.path(ROOT, "figures")
LEAD_OPT <- "opt2_dCT_CT"; LEAD_COMP <- "PLS2"; LEAD_DS <- "ds25"

# ---------------------------------------------------------------- data --------
y34    <- read.csv(file.path(DATA, "y_maps_bilateral_34.csv"), row.names = 1)
scores <- read.csv(file.path(RES, "lead_signature_scores.csv"))
wts    <- read.delim(file.path(RES, "lead_signature_weights.tsv"))
ct     <- read.delim(file.path(RES, "lead_celltype_enrichment.tsv"))
conc   <- read.delim(file.path(RES, "concordance_scores.tsv"))
wconc  <- read.delim(file.path(RES, "concordance_weights.tsv"))
ymaps  <- read.delim(file.path(RES, "ymaps_vs_nspn.tsv"))
comps  <- read.delim(file.path(RES, "pls_components.tsv"))
ovl    <- read.delim(file.path(RES, "lead_overlap_with_references.tsv"))
c3w    <- read.csv(file.path(REF, "ahba_c123_gene_weights.csv"))
nspnw  <- read.csv(file.path(REF, "nspn_pls_gene_weights.csv"))
poly   <- read.csv(file.path(DATA, "dk_polygons.csv"))

lead <- comps |> filter(option == LEAD_OPT, ds == LEAD_DS, component == LEAD_COMP)
# concordance rho is stored in pls.py's raw orientation (dCT salience positive);
# the figure shows the thinning orientation, so flip the sign.
cc <- function(ref) {
  r <- conc |> filter(option == LEAD_OPT, ds == LEAD_DS, component == LEAD_COMP, reference == ref)
  list(rho = -r$rho[1], p = r$p_spin[1])
}
ym <- function(map, ref) {
  r <- ymaps |> filter(y_map == map, reference == ref); list(rho = r$rho[1], p = r$p_spin[1])
}
fold <- function(ref, tail_) {
  r <- ovl |> filter(reference == ref, tail == tail_); r$fold[1]
}
pfmt <- function(p) if (p < 1e-3) "p_spin < 0.001" else sprintf("p_spin = %.3f", p)
rfmt <- function(x) sprintf("rho = %.2f", x)

# ---------------------------------------------------------------- brains ------
poly_lh <- poly |> filter(hemi == "left", view %in% c("lateral", "medial")) |>
  mutate(grp = interaction(label, view, group, subgroup, drop = TRUE),
         region = sub("^lh_", "", label))
brain_panel <- function(values, title, subtitle, diverging = TRUE) {
  d <- data.frame(region = sub("^lh_", "", names(values)), value = as.numeric(values))
  pd <- poly_lh |> left_join(d, by = "region")
  sc <- if (diverging) {
    lim <- max(abs(pd$value), na.rm = TRUE)
    scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-lim, lim), na.value = "grey92")
  } else {
    scale_fill_viridis_c(na.value = "grey92")
  }
  ggplot(pd, aes(x, y, group = grp, fill = value)) +
    geom_polygon(colour = "grey35", linewidth = 0.08) +
    sc + coord_fixed(expand = FALSE) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 7) +
    theme(legend.position = "none",
          plot.title = element_text(size = 7.4, face = "bold", hjust = 0, margin = margin(b = 1)),
          plot.subtitle = element_text(size = 6.3, colour = "grey30", hjust = 0, margin = margin(b = 1, l = 2)),
          plot.margin = margin(1, 2, 1, 2))
}

s_lh <- setNames(scores$thinning_scores_gene_side, scores$label)
nspn34 <- read.csv(file.path(REF, "nspn_dk_maps_bilateral_34.csv"), row.names = 1)
c3r <- read.csv(file.path(REF, sprintf("ahba_c123_scores_recomputed_%s.csv", LEAD_DS)), row.names = 1)

pa <- brain_panel(setNames(y34$dCT, rownames(y34)), "a   ABCD thinning rate",
                  "fastest in association cortex", diverging = FALSE)
pb <- brain_panel(s_lh, "b   ABCD PLS2 gene scores", "red = signature genes over-expressed")
pc <- brain_panel(setNames(nspn34$PLS2, rownames(nspn34)), "c   NSPN PLS2 (Whitaker 2016)",
                  "same spatial pattern, n = 297")
pd_ <- brain_panel(setNames(c3r$C3, rownames(c3r)), "d   AHBA C3 (Dear 2024)",
                   "imaging-free version of the axis")

# ---------------------------------------------------------------- scatters ----
base_theme <- theme_bw(base_size = 7.6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
        plot.title = element_text(size = 7.6, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, colour = "grey30", margin = margin(b = 2)),
        axis.title = element_text(size = 7), plot.margin = margin(2, 4, 2, 2))
theme_set(base_theme)

lab_regions <- c("entorhinal", "precentral")
sc_panel <- function(xv, yv, xlab, ylab, title, subtitle) {
  d <- data.frame(x = xv, y = yv, region = sub("^lh_", "", scores$label)) |>
    mutate(lab = ifelse(region %in% lab_regions, region, NA))
  ggplot(d, aes(x, y)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey75", linewidth = 0.4) +
    geom_point(size = 0.9, colour = "grey15") +
    geom_text_repel(aes(label = lab), size = 1.9, colour = "grey35", na.rm = TRUE,
                    min.segment.length = 0, segment.size = 0.2, box.padding = 0.3) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle)
}

e <- ym("dCT", "NSPN_CT_delta"); f <- cc("C3_recomputed_ds25"); g <- cc("NSPN_PLS2")
pe <- sc_panel(scores$NSPN_CT_delta, scores$dCT, "NSPN thinning rate (a.u.)", "ABCD thinning rate (mm/yr)",
               sprintf("e   Thinning maps agree (%s)", rfmt(e$rho)),
               sprintf("independent cohorts, %s", pfmt(e$p)))
pf <- sc_panel(scores$C3_recomputed, scores$thinning_scores_gene_side, "AHBA C3 score", "ABCD PLS2 score",
               sprintf("f   Scores track AHBA C3 (%s)", rfmt(f$rho)),
               sprintf("strongest of any reference map, %s", pfmt(f$p)))
pg <- sc_panel(scores$NSPN_PLS2, scores$thinning_scores_gene_side, "NSPN PLS2 score", "ABCD PLS2 score",
               sprintf("g   ...and NSPN PLS2 (%s)", rfmt(g$rho)),
               sprintf("recovered without the PLS2 genes, %s", pfmt(g$p)))

# ---------------------------------------------------------------- genes -------
z25 <- setNames(wts$thinning_Z_ds25, wts$gene)
# rho and n come from concordance_weights.tsv (thinning orientation = -rho there);
# the panel only draws the points.
cw <- function(ref) {
  r <- wconc |> filter(option == LEAD_OPT, ds == LEAD_DS, component == LEAD_COMP, reference == ref)
  list(rho = -r$rho[1], n = r$n_genes[1])
}
hex_panel <- function(ref_vec, xlab, title, subtitle, ref_key) {
  sh <- intersect(names(z25), names(ref_vec)); sh <- sh[!is.na(ref_vec[sh]) & !is.na(z25[sh])]
  d <- data.frame(x = as.numeric(ref_vec[sh]), y = as.numeric(z25[sh]))
  st <- cw(ref_key)
  stopifnot(length(sh) == st$n)
  rho <- st$rho
  ggplot(d, aes(x, y)) + geom_hex(bins = 38) +
    scale_fill_gradient(low = "grey88", high = "grey12", trans = "log10", guide = "none") +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
    labs(x = xlab, y = "ABCD PLS2 gene Z (bootstrap)",
         title = sprintf(title, rho), subtitle = sprintf(subtitle, st$n / 1000))
}
ph <- hex_panel(setNames(c3w$C3, c3w[[1]]), "AHBA C3 gene weight",
                "h   Gene weights track C3 (rho = %.2f)",
                sprintf("%%.1fk shared genes; top deciles overlap %.1fx", fold("C3", "top")), "C3")
pi_ <- hex_panel(setNames(nspnw$PLS2_z, nspnw$gene), "NSPN PLS2 gene Z",
                 "i   ...and NSPN PLS2 (rho = %.2f)",
                 sprintf("%%.1fk shared genes; top deciles overlap %.1fx", fold("NSPN_PLS2_z", "top")), "NSPN_PLS2_z")

ctp <- ct |> mutate(cell_class = factor(cell_class, levels = cell_class[order(z)]),
                    dir = ifelse(z > 0, "neuronal", "glial / vascular"))
pj <- ggplot(ctp, aes(z, cell_class, fill = dir)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_text(aes(label = paste0("n=", n_markers), hjust = ifelse(z > 0, -0.15, 1.15)),
            size = 1.8, colour = "grey35") +
  scale_fill_manual(values = c("neuronal" = "#b2182b", "glial / vascular" = "#2166ac"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  labs(x = "marker enrichment z", y = NULL,
       title = "j   Signature is neuronal, not glial",
       subtitle = sprintf("%d of %d Seidlitz 2020 sets p_perm<0.001\n(%s weakest, p = %.2f)",
                          sum(ct$p_perm < 0.001), nrow(ct),
                          ct$cell_class[which.max(ct$p_perm)], max(ct$p_perm)))

# ---------------------------------------------------------------- assemble ----
methods <- paste(
  "Methods.",
  "\u2022 Imaging (Y): bilateral DK maps from the settled ABCD 7.0 mixed model \u2014 thinning rate dCT (age slope) and baseline thickness CT (intercept).",
  "\u2022 Expression (X): AHBA left-hemisphere DK expression, 12,007 genes at the 25% differential-stability filter, z-scored per gene; 33 of 34 regions (no frontal pole).",
  sprintf("\u2022 Model: PLS-SVD of X'Y. PLS1 takes the static thickness gradient; PLS2 carries dCT (salience %.2f) and explains %.0f%% of the cross-covariance.",
          lead$sal_dCT[1], 100 * lead$cov_explained[1]),
  sprintf("\u2022 Inference: 5,000 spin rotations of the complete 34-parcel map (PLS2 p = %.3f); gene weights = Z over 1,000 region bootstraps, Procrustes-aligned.",
          lead$p_spin_singular[1]),
  "\u2022 Orientation: positive score / weight = expressed where thinning is FASTER. References: Whitaker 2016 PLS2 (symbols harmonised) and Dear 2024 C1\u2013C3.",
  sep = "\n")

fig <- (pa | pb | pc | pd_) /
  (pe | pf | pg) /
  (ph | pi_ | pj) +
  plot_layout(heights = c(0.62, 1, 1)) +
  plot_annotation(
    title = "An ABCD-derived transcriptomic signature of adolescent cortical thinning",
    subtitle = paste0("The second PLS component of gene expression against ABCD thinning reproduces both prior signatures \u2014 NSPN-PLS2 and AHBA-C3 \u2014 in regional\n",
                      "scores (e\u2013g) and in gene weights (h, i), and is neuronal rather than glial (j). Its disorder enrichment is in Figure 2."),
    caption = methods,
    theme = theme(plot.title = element_text(size = 9.5, face = "bold"),
                  plot.subtitle = element_text(size = 7.4, colour = "grey25", margin = margin(b = 4)),
                  plot.caption = element_text(size = 6.2, colour = "grey25", hjust = 0, lineheight = 1.45,
                                              margin = margin(t = 6))))

ggsave(file.path(FIG, "fig_lead_signature.png"), fig, width = 7.8, height = 7.6, dpi = 300, bg = "white")
cat("wrote", file.path(FIG, "fig_lead_signature.png"), "\n")
