#!/usr/bin/env Rscript
# fig7_cbcl_assoc.R -- EXPLORATORY: regional maps of the linear association
# between thinning rate and CBCL symptoms at ages ~15-17 given baseline symptoms,
# and their agreement with PLS2 / C3 / the normative thinning map.
# Reads results/cbcl_assoc_{maps,map_corr,map_partial}.tsv (25_cbcl_assoc_maps.py),
# results/hcp_summary_maps.csv, data/hcp_polygons.csv. Group-level only.
# Writes figures/fig_casecontrol.png
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results")
M  <- read.delim(file.path(RES, "cbcl_assoc_maps.tsv")) |> mutate(d = r)
CR <- read.delim(file.path(RES, "cbcl_assoc_map_corr.tsv")) |> mutate(p_label = p_perm)
PT <- read.delim(file.path(RES, "cbcl_assoc_map_partial.tsv"))
MAPV <- "relative_ct"
G  <- read.csv(file.path(RES, "hcp_summary_maps.csv"))
poly <- read.csv(file.path(ROOT, "data", "hcp_polygons.csv")) |> filter(view %in% c("lateral", "medial"))
TXT <- 7; FADE <- 0.3
OUT <- c(internal = "Internalising", depress = "Depressive (DSM)", anxdep = "Anxious/depressed",
         anxdisord = "Anxiety (DSM)", withdep = "Withdrawn/depressed", pfactor = "p-factor",
         rulebreak = "Rule-breaking")
stopifnot(setequal(unique(M$outcome), names(OUT)))
nn <- M |> distinct(outcome, n)
olab <- setNames(sprintf("%s  (n = %s)", OUT[nn$outcome], format(nn$n, big.mark = ",")), nn$outcome)
base_theme <- theme_classic(base_size = TXT) +
  theme(plot.title = element_text(size = TXT + 0.5, face = "bold", hjust = 0), strip.background = element_blank(),
        strip.text = element_text(size = TXT - 0.5))

# ---- a: maps (absolute d, cases - controls) ---------------------------------
SHOW <- c("pfactor", "anxdisord", "depress", "rulebreak")
md <- M |> filter(map == MAPV, outcome %in% SHOW) |>
  mutate(label = tolower(label))
lim <- max(abs(md$d)) * c(-1, 1)
pd <- poly |> mutate(label = tolower(label)) |> inner_join(md, by = "label", relationship = "many-to-many") |>
  mutate(outcome = factor(olab[outcome], olab[SHOW]))
pa <- ggplot(pd, aes(x, y, group = interaction(view, label, group, subgroup), fill = d)) +
  geom_polygon(colour = "grey40", linewidth = 0.05) + coord_fixed(expand = FALSE) +
  scale_fill_distiller(palette = "RdBu", limits = lim, name = "partial r\n(faster thinning,\nmore symptoms > 0)") +
  facet_wrap(~outcome, ncol = 1) + labs(title = "a  Thinning vs symptoms at 15\u201317", subtitle = "| baseline symptoms, global thinning, baseline CT") +
  theme_void(base_size = TXT) +
  theme(plot.title = element_text(size = TXT + 0.5, face = "bold"),
        plot.subtitle = element_text(size = TXT - 0.5, colour = "grey30"),
        strip.text = element_text(size = TXT - 0.5, margin = margin(b = 4)),
        panel.spacing.y = unit(6, "pt"),
        legend.key.height = unit(14, "pt"), legend.key.width = unit(6, "pt"))

# ---- b: internalising relative map vs PLS2 -----------------------------------
bi <- M |> filter(map == MAPV, outcome == "anxdisord") |> mutate(label = tolower(label)) |>
  inner_join(G |> mutate(label = tolower(label)) |> select(label, dCT), by = "label")
rb <- CR |> filter(outcome == "anxdisord", map == MAPV, reference == "dCT")
pb <- ggplot(bi, aes(dCT, d)) + geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(size = 1, colour = "grey30", stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey10", linewidth = 0.4) +
  annotate("text", -Inf, Inf, hjust = -0.08, vjust = 1.3, size = 2.2, fontface = "bold",
           label = sprintf("rho = %.2f, p_spin = %.3f, p_perm = %.3f", rb$rho, rb$p_spin, rb$p_label)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = "normative dCT (mm/yr; higher = slower thinning)", y = "partial r, anxiety (DSM)",
       title = "b  Anxiety map vs the normative thinning map") + base_theme

# ---- c: every map vs the three references ------------------------------------
cc <- CR |> filter(reference %in% c("PLS2", "C3", "dCT")) |>
  mutate(sig = p_spin < 0.05 & p_label < 0.05,
         outcome = factor(OUT[outcome], rev(OUT)),
         reference = factor(recode(reference, dCT = "dCT (higher = slower thinning)"), c("PLS2", "C3", "dCT (higher = slower thinning)")),
         map = factor(recode(map, absolute = "unadjusted", relative = "| global thinning",
                             absolute_ct = "| baseline CT", relative_ct = "| global thinning + baseline CT"),
                      c("unadjusted", "| global thinning", "| baseline CT", "| global thinning + baseline CT")))
pc <- ggplot(cc, aes(rho, outcome, colour = map, alpha = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 1.5, position = position_dodge(width = 0.7)) +
  scale_colour_manual(values = c(unadjusted = "grey25", "| global thinning" = "#b2182b",
                                 "| baseline CT" = "#4393c3", "| global thinning + baseline CT" = "#762a83"), name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  facet_wrap(~reference, nrow = 1) + scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
  labs(x = "Spearman rho of the association map with the reference map", y = NULL,
       title = "c  Map agreement (opaque: p_spin and p_perm < 0.05)") +
  base_theme + theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(colour = guide_legend(nrow = 2))

fig <- (pa | (pb / pc + plot_layout(heights = c(1, 1.25)))) + plot_layout(widths = c(0.75, 1)) +
  plot_annotation(
    title = "Later symptoms track thinning that departs from the normative gradient, not the PLS2 pattern",
    subtitle = paste0("Per parcel: CBCL (log1p raw sum) at ages ~15\u201317 ~ thinning rate + baseline CBCL + sex + site + ages + scans [+ global thinning] [+ that parcel's baseline CT]; n \u2248 8,200.\n",
                      "Internalising and anxiety maps correlate negatively with PLS2, but that is the dCT pattern (\u03c1 with PLS2 given dCT \u2248 0); with baseline CT, they align with slower-thinning cortex (dCT \u03c1 0.29\u20130.41)."),
    theme = theme(plot.title = element_text(size = TXT + 2.5, face = "bold"),
                  plot.subtitle = element_text(size = TXT, colour = "grey25", lineheight = 1.15)))
ggsave(file.path(ROOT, "figures", "fig_cbcl_assoc.png"), fig, width = 10, height = 6.8, dpi = 300, bg = "white")
cat("wrote fig_cbcl_assoc.png\n")
