#!/usr/bin/env Rscript
# fig_d5.R -- D5 summary figure. Reads ONLY results/table_d5_*.tsv (written by
# 02_models.R); recomputes nothing.
#   Rscript directions/d5_adversity_direction/code/fig_d5.R      # repo root, env ahba-pls-r
suppressMessages({library(ggplot2); library(patchwork); library(dplyr); library(readr)})
D5  <- "directions/d5_adversity_direction"
RES <- file.path(D5, "results"); OUT <- file.path(D5, "figures/fig_d5_adversity_direction.png")
rd  <- function(f) read_tsv(file.path(RES, f), show_col_types = FALSE)
A0 <- rd("table_d5_A0_reproduction.tsv"); A1 <- rd("table_d5_A1_env_slope.tsv")
A2 <- rd("table_d5_A2_rGE.tsv");          A3 <- rd("table_d5_A3_attenuation.tsv")
B12 <- rd("table_d5_B12_direction.tsv");  B3 <- rd("table_d5_B3_riclpm.tsv")

BASE <- 7
th <- theme_bw(base_size = BASE) + theme(
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
  plot.title = element_text(size = BASE, face = "bold"), plot.title.position = "plot",
  plot.subtitle = element_text(size = BASE - 1, colour = "grey25"),
  axis.text = element_text(size = BASE - 1), axis.title = element_text(size = BASE - 0.5),
  legend.text = element_text(size = BASE - 1), legend.title = element_blank(),
  legend.position = "bottom", legend.key.size = unit(7, "pt"), legend.margin = margin(0, 0, 0, 0))
ci <- function(d) mutate(d, lo = beta - 1.96 * se, hi = beta + 1.96 * se)
star <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
dodge <- position_dodge(width = 0.6)
vline <- geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50")
fmtp <- function(p) ifelse(p < .001, "< .001", sub("^0", "", sprintf("= %.3f", p)))

ENVLAB <- c(income = "Household income", parent_edu = "Parental education", adi = "Area deprivation (ADI)",
            conflict_y = "Family conflict (youth)", conflict_p = "Family conflict (parent)",
            bad_events_y = "Negative life events (youth)", bad_events_p = "Negative life events (parent)",
            ses = "SES composite", adversity = "Adversity composite")

# ---- a: environment -> slope ------------------------------------------------------
a <- A1 |> filter(outcome == "global_slope_1lmm", variant %in% c("base", "+site", "+site+qc")) |> ci() |>
  mutate(env = factor(ENVLAB[term], rev(ENVLAB)),
         variant = factor(variant, c("base", "+site", "+site+qc"),
                          c("genetics-arm covariates", "+ site", "+ site + scan quality")))
g <- function(t, v) a |> filter(term == t, variant == v)
sub_a <- sprintf("SES and ADI effects are between-site (ADI %.3f -> %.3f with site); parent-reported events survive (%.3f, p %s)",
                 g("adi", "genetics-arm covariates")$beta, g("adi", "+ site")$beta,
                 g("bad_events_p", "+ site + scan quality")$beta, fmtp(g("bad_events_p", "+ site + scan quality")$p))
pa <- ggplot(a, aes(beta, env, colour = variant)) + vline +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = dodge, size = 0.12, linewidth = 0.35) +
  scale_colour_manual(values = c("grey20", "#2166ac", "#92c5de")) +
  labs(title = "a   Environment and the thinning rate", subtitle = sub_a,
       x = "beta (SD thinning slope per SD; < 0 = faster thinning)", y = NULL) + th

# ---- b: gene-environment correlation ------------------------------------------------
SC <- c(SCZ_pooled = "SCZ, pooled", MDD_pooled = "MDD, pooled", SCZ_eur = "SCZ, EUR", MDD_eur = "MDD, EUR",
        EA_eur = "EA, EUR")
SC <- SC[names(SC) %in% A0$score]
b <- A2 |> filter(env %in% c("ses", "adversity")) |> ci() |>
  mutate(score = factor(SC[score], rev(SC)), env = factor(ENVLAB[env], ENVLAB[c("ses", "adversity")]))
bm <- b |> filter(score == "MDD, pooled", env == "Adversity composite")
pb <- ggplot(b, aes(beta, score, colour = env)) + vline +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = dodge, size = 0.12, linewidth = 0.35) +
  scale_colour_manual(values = c("#1b7837", "#b2182b")) +
  labs(title = "b   Polygenic scores and the environment",
       subtitle = sprintf("MDD risk goes with more adversity (pooled %.3f, p %s) and lower SES",
                          bm$beta, fmtp(bm$p)),
       x = "beta (SD environment per SD score)", y = NULL) + th

# ---- c: PRS -> slope with environment in the model ----------------------------------
MOD <- c(`PRS alone` = "score alone", `+ses+adversity` = "+ SES + adversity",
         `+ses+adversity+site` = "+ SES + adversity + site")
cc <- A3 |> filter(model %in% names(MOD)) |> ci() |>
  mutate(score = factor(SC[score], rev(SC)), model = factor(MOD[model], MOD))
att <- A3 |> filter(model == "+ses+adversity") |> mutate(s = SC[score])
A4 <- rd("table_d5_A4_gxe.tsv") |> filter(role == "PRS x environment")
sub_c <- sprintf("Environment explains %.0f%% of the SCZ and %.0f%% of the MDD effect (pooled)\nPRS x environment interaction: %d/%d tests p < .05",
                 att$attenuation_pct[att$score == "SCZ_pooled"], att$attenuation_pct[att$score == "MDD_pooled"],
                 sum(A4$p < .05), nrow(A4))
pc <- ggplot(cc, aes(beta, score, colour = model)) + vline +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = dodge, size = 0.12, linewidth = 0.35) +
  geom_text(aes(x = hi, label = star(p)), position = dodge, hjust = -0.2, size = (BASE - 1.5) / .pt,
            show.legend = FALSE) +
  scale_colour_manual(values = c("grey20", "#e08214", "#8073ac")) +
  labs(title = "c   Polygenic effect on thinning, net of environment", subtitle = sub_c,
       x = "beta (SD thinning slope per SD score)", y = NULL) + th

# ---- d: between-child direction ------------------------------------------------------
OUTL <- c(pfactor = "p-factor", internal = "Internalising", external = "Externalising", depress = "Depressive (DSM)")
dd <- B12 |> filter(model == "base") |> ci() |>
  mutate(outcome = factor(OUTL[outcome], rev(OUTL)),
         direction = factor(direction, c("slope -> symptoms at 15-17 | baseline", "baseline symptoms -> subsequent slope"),
                            c("thinning rate -> symptoms at 15-17 | baseline", "baseline symptoms -> thinning rate")))
nf <- dd |> filter(direction == levels(direction)[1]) ; nr <- dd |> filter(direction == levels(direction)[2])
pd <- ggplot(dd, aes(beta, outcome, colour = direction)) + vline +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = dodge, size = 0.12, linewidth = 0.35) +
  geom_text(aes(x = lo, label = star(p)), position = dodge, hjust = 1.3, size = (BASE - 1.5) / .pt,
            show.legend = FALSE) +
  scale_colour_manual(values = c("#b2182b", "grey55")) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(title = "d   Between children: forward, not reverse",
       subtitle = sprintf("Thinning predicts later symptoms (%d/4 p < .05); baseline symptoms predict thinning in %d/4",
                          sum(nf$p < .05), sum(nr$p < .05)),
       x = "beta (SD per SD; < 0 = faster thinning with more symptoms)", y = NULL) + th

# ---- e: within-child cross-lags --------------------------------------------------------
e <- B3 |> filter(label %in% c("cxy", "cyx")) |>
  mutate(se_std = se * std_all / est, lo = std_all - 1.96 * se_std, hi = std_all + 1.96 * se_std,
         outcome = factor(OUTL[outcome], rev(OUTL[c("internal", "external", "depress")])),
         path = factor(path, c("CT_t -> symptoms_t+2y", "symptoms_t -> CT_t+2y"),
                       c("thickness -> symptoms 2 y later", "symptoms -> thickness 2 y later")),
         thickness = factor(thickness, c("age+sex+site", "age+sex+site+scan quality")))
ep <- e |> filter(thickness == "age+sex+site")
pe <- ggplot(e, aes(std_all, outcome, colour = path, linetype = thickness, group = interaction(path, thickness))) + vline +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.75), size = 0.12, linewidth = 0.35,
                  key_glyph = "path") +
  scale_colour_manual(values = c("#b2182b", "grey35")) +
  scale_linetype_manual(values = c("solid", "22"), labels = c("thickness | age, sex, site", "dashed: + scan quality")) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2, override.aes = list(colour = "grey20"))) +
  labs(title = "e   Within children: small lags in both directions",
       subtitle = sprintf("RI-CLPM: symptoms -> thinner cortex %d/3; thinner cortex -> symptoms %d/3 (p < .05)",
                          sum(ep$pvalue[ep$label == "cyx"] < .05), sum(ep$pvalue[ep$label == "cxy"] < .05)),
       x = "standardised cross-lagged path (< 0 = thinner / more symptoms)", y = NULL) + th

n_all <- max(A1$n); n_eur <- A0$n[A0$score == "SCZ_eur"]
cap <- paste(
  sprintf("\u2022 ABCD 7.0, HCP-MMP; %s genotyped children (EUR arm %s). Thinning rate = standardised random slope of the single LMM on per-scan cortical mean thickness (Figure 1).", format(n_all, big.mark = ","), format(n_eur, big.mark = ",")),
  "\u2022 Environment at baseline: income, parental education, ADI national percentile, Family Environment Scale conflict; negative life-event counts at year 1. Composites = mean of z scores (SES: income, education, -ADI).",
  "\u2022 a-c: lmer, outcome and predictors standardised, sex + baseline age + PC1-10 + (1 | family), as in tools/prs_assoc.R; PRS-CS scores in matched arms (pooled = within-ancestry z; EUR = raw).",
  sprintf("\u2022 Reproduction check: all four PRS -> thinning betas match the cluster table (max |diff| = %.1e). Scan quality = log1p FreeSurfer surface-defect count.", max(A0$abs_diff)),
  "\u2022 d: symptoms at 15-17 (mean of years 5-7) ~ slope + baseline symptoms + age + sex + site (OLS, family-clustered SE), as Figure 1g; reverse = slope ~ baseline (year 0) symptoms.",
  "\u2022 e: random-intercept cross-lagged panel model, whole-cortex thickness and CBCL at years 0/2/4/6, each residualised within wave; lags equal over time; MLR, FIML, family clusters.",
  "\u2022 CBCL = parent report, log1p raw sums. Error bars 95% CI; * p < .05, ** < .01, *** < .001, uncorrected.", sep = "\n")

fig <- (pa | pb) / (pc | pd | pe) + plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Adversity and polygenic risk act on the thinning rate largely separately; the thinning -> symptom link runs forward",
    subtitle = "Parent-reported life events predict faster thinning (a); SES effects are between-site. Environment explains little of the polygenic effect (b, c).\nBetween children, thinning predicts later symptoms but baseline symptoms do not predict thinning (d); within children, cross-lags are small in both directions (e).",
    caption = cap,
    theme = theme(plot.title = element_text(size = BASE + 1, face = "bold"),
                  plot.subtitle = element_text(size = BASE, colour = "grey20"),
                  plot.caption = element_text(size = BASE - 1.2, hjust = 0, lineheight = 1.4)))
ggsave(OUT, fig, width = 11, height = 7.6, dpi = 220)
cat(OUT, "\n")
