#!/usr/bin/env Rscript
# fig6_casecontrol.R -- EXPLORATORY: regional case-control maps of adolescent
# thinning (children who developed symptoms vs those who did not) and their
# agreement with PLS2 / C3 / the normative thinning map.
# Reads results/casecontrol_{maps,map_corr}.tsv (24_casecontrol_maps.py),
# results/hcp_summary_maps.csv, data/hcp_polygons.csv. Group-level only.
# Writes figures/fig_casecontrol.png
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})
ROOT <- "/Users/richard/Git/abcd_development/ahba_pls"
RES <- file.path(ROOT, "results")
M  <- read.delim(file.path(RES, "casecontrol_maps.tsv"))
CR <- read.delim(file.path(RES, "casecontrol_map_corr.tsv"))
G  <- read.csv(file.path(RES, "hcp_summary_maps.csv"))
poly <- read.csv(file.path(ROOT, "data", "hcp_polygons.csv")) |> filter(view %in% c("lateral", "medial"))
TXT <- 7; FADE <- 0.3
OUT <- c(internal = "Internalising", depress = "Depressive (DSM)", anxdep = "Anxious/depressed",
         anxdisord = "Anxiety (DSM)", withdep = "Withdrawn/depressed", pfactor = "p-factor",
         rulebreak = "Rule-breaking", mdd_parent = "MDD (KSADS, parent)")
stopifnot(setequal(unique(M$outcome), names(OUT)))
nn <- M |> distinct(outcome, n_case, n_ctrl)
olab <- setNames(sprintf("%s\n%s cases / %s controls", OUT[nn$outcome], format(nn$n_case, big.mark = ","),
                         format(nn$n_ctrl, big.mark = ",")), nn$outcome)
base_theme <- theme_classic(base_size = TXT) +
  theme(plot.title = element_text(size = TXT + 0.5, face = "bold", hjust = 0), strip.background = element_blank(),
        strip.text = element_text(size = TXT - 0.5))

# ---- a: maps (absolute d, cases - controls) ---------------------------------
SHOW <- c("internal", "depress", "pfactor", "mdd_parent")
md <- M |> filter(map == "absolute", outcome %in% SHOW) |>
  mutate(label = tolower(label))
lim <- max(abs(md$d)) * c(-1, 1)
pd <- poly |> mutate(label = tolower(label)) |> inner_join(md, by = "label", relationship = "many-to-many") |>
  mutate(outcome = factor(olab[outcome], olab[SHOW]))
pa <- ggplot(pd, aes(x, y, group = interaction(view, label, group, subgroup), fill = d)) +
  geom_polygon(colour = "grey40", linewidth = 0.05) + coord_fixed(expand = FALSE) +
  scale_fill_distiller(palette = "RdBu", limits = lim, name = "Cohen's d\n(cases thin\nfaster > 0)") +
  facet_wrap(~outcome, ncol = 1) + labs(title = "a  Case \u2212 control difference in thinning rate") +
  theme_void(base_size = TXT) +
  theme(plot.title = element_text(size = TXT + 0.5, face = "bold"), strip.text = element_text(size = TXT - 0.5),
        legend.key.height = unit(14, "pt"), legend.key.width = unit(6, "pt"))

# ---- b: internalising relative map vs PLS2 -----------------------------------
bi <- M |> filter(map == "relative", outcome == "internal") |> mutate(label = tolower(label)) |>
  inner_join(G |> mutate(label = tolower(label)) |> select(label, PLS2), by = "label") |> filter(!is.na(PLS2))
rb <- CR |> filter(outcome == "internal", map == "relative", reference == "PLS2")
pb <- ggplot(bi, aes(PLS2, d)) + geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(size = 1, colour = "grey30", stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey10", linewidth = 0.4) +
  annotate("text", -Inf, Inf, hjust = -0.08, vjust = 1.3, size = 2.2, fontface = "bold",
           label = sprintf("rho = %.2f, p_spin = %.3f, p_label = %.3f", rb$rho, rb$p_spin, rb$p_label)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = "PLS2 score", y = "Cohen's d, internalising | global thinning", title = "b  Internalising vs PLS2") + base_theme

# ---- c: every map vs the three references ------------------------------------
cc <- CR |> filter(reference %in% c("PLS2", "C3", "dCT")) |>
  mutate(sig = p_spin < 0.05 & p_label < 0.05,
         outcome = factor(OUT[outcome], rev(OUT)),
         reference = factor(recode(reference, dCT = "normative dCT"), c("PLS2", "C3", "normative dCT")),
         map = factor(recode(map, absolute = "absolute", relative = "| global thinning"), c("absolute", "| global thinning")))
pc <- ggplot(cc, aes(rho, outcome, colour = map, alpha = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 1.6, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(absolute = "grey25", "| global thinning" = "#b2182b"), name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = FADE), guide = "none") +
  facet_wrap(~reference, nrow = 1) + scale_x_continuous(breaks = c(-0.3, 0, 0.3)) +
  labs(x = "Spearman rho of the case-control map with the reference map", y = NULL,
       title = "c  Map agreement (opaque: p_spin and p_label < 0.05)") +
  base_theme + theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank())

fig <- (pa | (pb / pc + plot_layout(heights = c(1, 1.25)))) + plot_layout(widths = c(0.75, 1)) +
  plot_annotation(
    title = "Children who develop symptoms thin slightly faster overall, but not in the PLS2 pattern",
    subtitle = paste0("Case = below the CBCL borderline cut-off at baseline and at/above it at ages ~15\u201317 (KSADS: lifetime diagnosis); ",
                      "control = below it at every wave.\nOnly internalising reaches p < 0.05 against PLS2, in the opposite direction: ",
                      "cases' extra thinning sits where PLS2 is low."),
    theme = theme(plot.title = element_text(size = TXT + 2.5, face = "bold"),
                  plot.subtitle = element_text(size = TXT, colour = "grey25", lineheight = 1.15)))
ggsave(file.path(ROOT, "figures", "fig_casecontrol.png"), fig, width = 10, height = 6.4, dpi = 300, bg = "white")
cat("wrote fig_casecontrol.png\n")
