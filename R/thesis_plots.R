#' Plotting helpers for notebooks/03_thesis_figures.qmd
#'
#' These are a reworked subset of the plotting code from the original
#' `thesis_abcd.qmd` rough analysis (`ABCD/code/{plot_mixedeffects.r,
#' plot_pcs.r,brain_plots.r}` in the sibling ABCD repo). They have been
#' adapted, not copied, for three reasons:
#'
#'  1. Data now comes from this repo's model-fit pipeline
#'     (`out/<run_id>/{model_table,phenotypes,fits}`) instead of raw ABCD
#'     release files, so the functions take tidy data frames with fixed
#'     column names (`region`, `age`, `value`, `pred`, `slope`, `sex`,
#'     `subject`) rather than release-specific columns.
#'  2. HCP-MMP (Glasser) regions render via `ggsegGlasser`. That package's
#'     current GitHub `main` targets a newer sf-based ggseg ("cortical_atlas")
#'     that is incompatible with the `ggseg 1.6.6` installed here (classic
#'     `brain_atlas`); this repo installs the older, compatible tag
#'     `LCBC-UiO/ggsegGlasser@1.0.01` instead (see `data/README.md` / repo
#'     memory for the install command). Its glasser labels are
#'     `{lh,rh}_{L,R}_<region>` (e.g. `lh_L_V1`), one `L_`/`R_` letter off
#'     from the pipeline's own `{lh,rh}_<region>` labels, so `plot_map()`
#'     rewrites labels before handing off to ggseg. Any atlas string other
#'     than "dk"/"hcp" falls back to a centroid scatter (parcel centres on
#'     the sphere, from `data/<atlas>_centroids.csv`), so the function still
#'     degrades gracefully for an atlas without ggseg support.
#'  3. Correlation annotations are computed by hand (`add_corr_labels`)
#'     rather than via `ggpubr::stat_cor`, so this file has one fewer
#'     dependency than the original; `ggpubr` is still used elsewhere for
#'     convenience but is not required by this file.
#'
#' Callers must set a global `DATA_DIR` (path to this repo's `data/`) before
#' sourcing, e.g. `DATA_DIR <- file.path(ROOT, "data")`, so `plot_map()` can
#' find the centroid files for the scatter fallback.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(ggseg)
  library(ggsegGlasser)
  library(pals)
})

.rdbu <- function(n = 100) rev(pals::brewer.rdbu(n))
.spectral <- function(n = 100) rev(pals::brewer.spectral(n))

#' Drop an atlas's non-cortical placeholder rows (medial wall / corpus
#' callosum), which ggseg otherwise renders as its own polygon. Left in,
#' `ggseg()`'s join against our (necessarily incomplete) per-region data gives
#' those rows an NA value for every one of our columns -- including the
#' faceting variable -- which shows up as a spurious extra "NA" facet that
#' doesn't fit the rest of the layout.
.drop_medial_wall <- function(atlas) {
  atlas$data <- atlas$data %>% filter(!is.na(region))
  atlas
}

#' The full set of valid (left-hemisphere, non-medial-wall) region labels an
#' atlas can render, in this pipeline's own `{lh,rh}_<region>` naming.
#'
#' Used to fill in any region our own map data doesn't cover (e.g. AHBA only
#' has scores for 137 of 180 HCP-MMP regions) with an explicit NA *before*
#' handing off to `ggseg()`. Without this, `ggseg()`'s internal join adds
#' those atlas-only regions back in with no map assigned at all (there is no
#' existing per-map row to attach the geometry to), and ggplot renders that
#' as its own spurious "NA" facet alongside the real ones.
.atlas_labels_ours <- function(atlas) {
  if (atlas == "dk") {
    d <- .drop_medial_wall(dk)$data
    unique(d$label[d$hemi == "left"])
  } else if (atlas == "hcp") {
    d <- .drop_medial_wall(glasser)$data
    str_replace(unique(d$label[d$hemi == "left"]), "^lh_L_", "lh_")
  } else {
    NULL
  }
}

# --------------------------------------------------------------------------
# Brain maps
# --------------------------------------------------------------------------

#' Plot one or more region -> value maps.
#'
#' @param maps data.frame with a `label` column (or rownames) of parcel
#'   labels in this pipeline's own format (e.g. `lh_precentral`, `lh_V1`) and
#'   one column per map to draw.
#' @param atlas "dk" (ggseg DK surface), "hcp" (ggsegGlasser HCP-MMP surface),
#'   or anything else (centroid scatter fallback using
#'   `data/<atlas>_centroids.csv`).
plot_map <- function(maps, atlas = "dk", title = "", ncol = 3,
                     colorscale = c("fixed", "symmetric", "zero"),
                     limits = c(-3, 3), labels = c("-3sd", "+3sd"), name = "") {
  colorscale <- match.arg(colorscale)
  if ("label" %in% colnames(maps)) {
    maps <- maps %>% remove_rownames() %>% column_to_rownames("label")
  }
  df <- maps %>%
    rownames_to_column("label") %>%
    pivot_longer(-label, names_to = "map", values_to = "value") %>%
    mutate(map = factor(map, levels = unique(map)))

  atlas_labels <- .atlas_labels_ours(atlas)
  if (!is.null(atlas_labels)) {
    df <- df %>% complete(label = union(label, atlas_labels), map)
  }

  if (colorscale == "fixed") {
    m_min <- limits[1]; m_max <- limits[2]
  } else if (colorscale == "symmetric") {
    m_max <- max(abs(df$value), na.rm = TRUE); m_min <- -m_max
  } else {
    m_min <- 0; m_max <- max(df$value, na.rm = TRUE)
  }

  if (atlas == "dk") {
    p <- df %>%
      ggseg(atlas = .drop_medial_wall(dk), hemi = "left", mapping = aes(fill = value),
            colour = "grey", size = .1, show.legend = TRUE) +
      scale_fill_gradientn(colors = .rdbu(100), limits = c(m_min, m_max), oob = squish,
                           breaks = c(m_min, m_max), labels = labels, name = name,
                           guide = guide_colorbar(barheight = .3, barwidth = 3)) +
      theme_void()
  } else if (atlas == "hcp") {
    # ggsegGlasser labels are {lh,rh}_{L,R}_<region>; ours are {lh,rh}_<region>.
    df <- df %>%
      mutate(label = str_replace(label, "^lh_", "lh_L_"),
            label = str_replace(label, "^rh_", "rh_R_"))
    p <- df %>%
      ggseg(atlas = .drop_medial_wall(glasser), hemi = "left", mapping = aes(fill = value),
            colour = "grey", size = .05, show.legend = TRUE) +
      scale_fill_gradientn(colors = .rdbu(100), limits = c(m_min, m_max), oob = squish,
                           breaks = c(m_min, m_max), labels = labels, name = name,
                           guide = guide_colorbar(barheight = .3, barwidth = 3)) +
      theme_void()
  } else {
    centroids <- read_csv(file.path(DATA_DIR, paste0(atlas, "_centroids.csv")),
                          show_col_types = FALSE) %>%
      filter(hemi == "lh") %>%
      mutate(label = str_remove(label, "^lh_"))
    df <- df %>% mutate(label = str_remove(label, "^lh_")) %>% inner_join(centroids, by = "label")
    p <- df %>%
      ggplot(aes(x = x, y = y, color = value)) +
      geom_point(size = 2.4) +
      coord_fixed() +
      scale_color_gradientn(colors = .rdbu(100), limits = c(m_min, m_max), oob = squish,
                            breaks = c(m_min, m_max), labels = labels, name = name,
                            guide = guide_colorbar(barheight = .3, barwidth = 3)) +
      theme_void()
  }

  p + facet_wrap(~map, ncol = ncol, dir = "v") +
    theme(legend.position = "bottom", legend.title = element_text(vjust = 1),
         strip.text = element_text(size = 8), plot.title = element_text(hjust = .5),
         plot.tag = element_blank()) +   # brain maps never carry a patchwork panel tag
    ggtitle(title)
}

# --------------------------------------------------------------------------
# Mixed-effects trajectory panel
# --------------------------------------------------------------------------

plot_raw_points <- function(df, regions_order) {
  df %>%
    mutate(region = factor(region, levels = regions_order, ordered = TRUE)) %>%
    ggplot(aes(x = age, y = value, color = sex)) +
    facet_grid(. ~ region) +
    geom_point(size = .4, alpha = .2, stroke = NA) +
    geom_smooth(method = "lm", linewidth = .5, se = FALSE) +
    scale_x_continuous("Age", breaks = c(9, 11, 13, 15)) +
    scale_y_continuous("Thickness (mm)") +
    theme_classic() +
    theme(legend.title = element_blank(), text = element_text(size = 8),
         strip.background = element_blank())
}

plot_subject_sample <- function(df_sample, regions_order) {
  df_sample %>%
    mutate(region = factor(region, levels = regions_order, ordered = TRUE)) %>%
    ggplot(aes(x = age, y = value, color = slope, group = subject)) +
    facet_grid(. ~ region) +
    geom_point(size = .3) +
    geom_line(linewidth = .2) +
    scale_color_viridis_c("Slope\n(mm/yr)") +
    scale_x_continuous("Age", breaks = c(9, 11, 13, 15)) +
    scale_y_continuous("Thickness (mm)") +
    theme_classic() +
    theme(text = element_text(size = 8), strip.background = element_blank())
}

plot_subject_preds <- function(preds_sample, regions_order) {
  preds_sample %>%
    mutate(region = factor(region, levels = regions_order, ordered = TRUE)) %>%
    ggplot(aes(x = age, y = pred, color = slope, group = subject)) +
    facet_grid(. ~ region) +
    geom_point(size = .3) +
    geom_line(linewidth = .2) +
    scale_color_viridis_c("Slope\n(mm/yr)") +
    scale_x_continuous("Age", breaks = c(9, 11, 13, 15)) +
    scale_y_continuous("Model-implied thickness (mm)") +
    theme_classic() +
    theme(text = element_text(size = 8), strip.background = element_blank())
}

# --------------------------------------------------------------------------
# Structural covariance matrix + PCA
# --------------------------------------------------------------------------

plot_corr_matrix <- function(R, region_order, name = "r", lim = 1, legend_title = "R") {
  R %>%
    mutate(x = factor(x, ordered = TRUE, levels = region_order),
          y = factor(y, ordered = TRUE, levels = region_order)) %>%
    ggplot(aes(x = x, y = y, fill = r)) +
    geom_raster() +
    scale_fill_gradientn(name = legend_title, colors = .spectral(100),
                         limits = c(-lim, lim), oob = squish, breaks = c(-lim, 0, lim)) +
    guides(fill = guide_colorbar(direction = "horizontal", barheight = .5)) +
    theme_classic() +
    theme(aspect.ratio = 1, axis.line = element_blank(), axis.text = element_blank(),
         axis.ticks = element_blank(), axis.title = element_blank(),
         legend.position = c(.5, -.1), legend.title = element_text(vjust = 1, face = "italic")) +
    ggtitle(name)
}

#' Per-facet Pearson correlation, for annotating scatter plots without ggpubr.
add_corr_labels <- function(df, xcol, ycol, facet_cols) {
  df %>%
    group_by(across(all_of(facet_cols))) %>%
    summarise(
      r = suppressWarnings(cor(.data[[xcol]], .data[[ycol]], use = "pairwise.complete.obs")),
      label = sprintf("r=%.2f", r),
      .groups = "drop"
    )
}

plot_ahba_scatter <- function(df, highlight = NULL) {
  if (!("highlight" %in% colnames(df))) {
    df$highlight <- FALSE
  }
  corr_labels <- add_corr_labels(df, "AHBA_score", "ABCD_score", c("ABCD", "AHBA"))

  df %>%
    ggplot(aes(x = AHBA_score, y = ABCD_score)) +
    facet_grid(ABCD ~ AHBA, scales = "free_y") +
    geom_point(aes(fill = ABCD_score), shape = 21, stroke = .5, color = "darkgrey") +
    scale_fill_gradientn(colors = .rdbu(100), guide = "none") +
    geom_smooth(aes(color = highlight), method = "lm", se = FALSE) +
    scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "forestgreen"), guide = "none") +
    geom_text(data = corr_labels, aes(label = label), x = -Inf, y = Inf,
             hjust = -.2, vjust = 1.4, size = 3, color = "black", inherit.aes = FALSE) +
    xlab("AHBA z-score") + ylab("ABCD z-score") +
    theme_classic() +
    theme(panel.spacing = unit(1, "lines"), strip.background = element_blank())
}
