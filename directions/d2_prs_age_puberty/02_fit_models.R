#!/usr/bin/env Rscript
# D2 step 2 -- one-stage longitudinal models: does polygenic risk change the RATE of cortical
# thinning uniformly across 9-17 (rate scaling), or does its effect depend on age / pubertal
# stage (a shift in developmental timing)?  docs/DIRECTIONS.md D2, prediction T1.
#
# Trait: per-scan mean HCP-MMP cortical thickness (the Figure 1 single-LMM trait), age centred
# at 12.797.  All models carry sex + PC1-10 and (1 + age_c | subject) + (1 | site) + (1 | family_id).
#
#   linear    mean_ct ~ prs * age_c                       prs:age_c = PRS effect on the rate (mm/yr per SD)
#   hinge13   mean_ct ~ prs * (age_c + h13)               h13 = max(age - 13, 0): rate before 13 and its
#                                                         change after 13.  prs:h13 is the T1 test.
#   quad      mean_ct ~ prs * (age_c + age_c^2)           curvature of the PRS effect
#   seg3      mean_ct ~ prs * (age_c + h12 + h14)         rate effects in <12, 12-14, >14
#   puberty   mean_ct ~ prs * age_c + prs * pds_c         pds_c = PDS mean at the scan, centred; does the
#                                                         PRS x age term move once PRS x puberty is in?
#   puberty0  the linear model on the puberty rows        the matched baseline for that comparison
#   pubonly   mean_ct ~ age_c + prs * pds_c               PRS effect per unit pubertal stage at fixed age
#
# Puberty tempo: PDS mean (youth) ~ age_c * sex + (1 + age_c | subject) + (1 | site) over all yearly
# waves; per-child BLUPs = PDS at 12.8 (stage) and PDS slope (tempo).  Then PRS -> stage / tempo,
# and stage / tempo -> the Figure 1 thinning slope, each with (1 | family_id).
#
# Inputs: the gitignored tables from 01_build_tables.py.  Outputs (aggregate only, committed):
#   tables/onestage_terms.tsv        every fixed-effect term of every model
#   tables/prs_rate_by_age.tsv       PRS effect on the rate per age segment (hinge13, seg3), with SE
#   tables/variance_components.tsv   tau0, tau1, site, family, residual per model
#   tables/puberty_tempo.tsv         PRS -> PDS stage/tempo; PDS stage/tempo -> thinning slope
#   tables/trajectories.tsv          model-implied mean_ct at ages 9-17 for PRS = -1, 0, +1 SD
# Run (env r, from the repo root):  Rscript directions/d2_prs_age_puberty/02_fit_models.R
suppressMessages({library(data.table); library(lme4); library(lmerTest)})
WORK <- "genetic_analysis/work/results_70tab_hcp/d2_prs_age_puberty"
OUT  <- "directions/d2_prs_age_puberty/tables"
CEN  <- 12.797
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

sc <- fread(file.path(WORK, "scans.csv"))
sc[, `:=`(sex = factor(sex), site = factor(site), family_id = factor(family_id), subject = factor(subject))]
sc[, `:=`(h12 = pmax(age - 12, 0), h13 = pmax(age - 13, 0), h14 = pmax(age - 14, 0), age_c2 = age_c^2)]
sc[, pds_c := pds_mean_y - mean(pds_mean_y, na.rm = TRUE)]
sc[, pds_p_c := pds_mean_p - mean(pds_mean_p, na.rm = TRUE)]
pcs <- paste0("PC", 1:10)
covs <- paste(c("sex", pcs), collapse = " + ")
re   <- "(1 + age_c | subject) + (1 | site) + (1 | family_id)"
ctrl <- lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)

fit <- function(rhs, d) {
  f <- as.formula(paste("mean_ct ~", rhs, "+", covs, "+", re))
  lmer(f, data = d, REML = TRUE, control = ctrl)
}
terms_of <- function(m, key) {
  s <- as.data.table(coef(summary(m)), keep.rownames = "term")
  setnames(s, c("term", "est", "se", "df", "t", "p"))
  cbind(key, s[!term %in% c("(Intercept)", "sexM", "sexF", pcs)], n_obs = nobs(m), n_children = nlevels(droplevels(m@frame$subject)))
}
vc_of <- function(m, key) {
  v <- as.data.table(VarCorr(m))
  cbind(key, v[, .(grp, var1, var2, vcov, sdcor)])
}
# linear combination of fixed effects: est, se, z, p
lincomb <- function(m, w) {
  b <- fixef(m); V <- as.matrix(vcov(m)); w <- w[names(b)]; w[is.na(w)] <- 0
  est <- sum(w * b); se <- sqrt(as.numeric(t(w) %*% V %*% w))
  data.table(est = est, se = se, z = est / se, p = 2 * pnorm(-abs(est / se)))
}
seg_rows <- function(m, key, segs) {
  rbindlist(lapply(names(segs), function(nm) {
    w <- setNames(rep(0, length(fixef(m))), names(fixef(m))); w[segs[[nm]]] <- 1
    cbind(key, segment = nm, lincomb(m, w))
  }))
}

TERMS <- list(); VC <- list(); SEG <- list(); TRAJ <- list()
grid <- CJ(disorder = c("scz", "mdd"), method = c("sbrc", "prscs"), arm = c("pooled", "eur"))
for (i in seq_len(nrow(grid))) {
  g <- grid[i]; col <- sprintf("%s_%s_%s", g$disorder, g$method, g$arm)
  d <- copy(sc)[!is.na(get(col))]; d[, prs := get(col)]
  key <- data.table(disorder = toupper(g$disorder), method = g$method, arm = g$arm)
  cat(sprintf("== %s: %d scans, %d children\n", col, nrow(d), uniqueN(d$subject)))

  m1 <- fit("prs * age_c", d)
  TERMS[[length(TERMS) + 1]] <- terms_of(m1, cbind(key, model = "linear")); VC[[length(VC) + 1]] <- vc_of(m1, cbind(key, model = "linear"))
  if (g$method == "sbrc") {
    fx <- fixef(m1); psex <- mean(d$sex == levels(d$sex)[2]); sexnm <- grep("^sex", names(fx), value = TRUE)
    TRAJ[[length(TRAJ) + 1]] <- rbindlist(lapply(c(-1, 0, 1), function(z) {
      a <- seq(9, 17, 0.5); t <- a - CEN
      cbind(key, prs_sd = z, data.table(age = a, mean_ct = fx[["(Intercept)"]] + fx[["age_c"]] * t + fx[["prs"]] * z +
                                          fx[["prs:age_c"]] * z * t + fx[[sexnm]] * psex))
    }))
  }

  m2 <- fit("prs * (age_c + h13)", d)
  TERMS[[length(TERMS) + 1]] <- terms_of(m2, cbind(key, model = "hinge13")); VC[[length(VC) + 1]] <- vc_of(m2, cbind(key, model = "hinge13"))
  SEG[[length(SEG) + 1]] <- seg_rows(m2, cbind(key, model = "hinge13"),
                                     list("<13" = "prs:age_c", ">=13" = c("prs:age_c", "prs:h13")))

  m3 <- fit("prs * (age_c + age_c2)", d)
  TERMS[[length(TERMS) + 1]] <- terms_of(m3, cbind(key, model = "quad"))

  m4 <- fit("prs * (age_c + h12 + h14)", d)
  TERMS[[length(TERMS) + 1]] <- terms_of(m4, cbind(key, model = "seg3"))
  SEG[[length(SEG) + 1]] <- seg_rows(m4, cbind(key, model = "seg3"),
                                     list("<12" = "prs:age_c", "12-14" = c("prs:age_c", "prs:h12"),
                                          ">=14" = c("prs:age_c", "prs:h12", "prs:h14")))

  dp <- d[!is.na(pds_c)]
  m50 <- fit("prs * age_c", dp)
  TERMS[[length(TERMS) + 1]] <- terms_of(m50, cbind(key, model = "puberty0"))
  m5 <- fit("prs * age_c + prs * pds_c", dp)
  TERMS[[length(TERMS) + 1]] <- terms_of(m5, cbind(key, model = "puberty")); VC[[length(VC) + 1]] <- vc_of(m5, cbind(key, model = "puberty"))
  m5b <- fit("age_c + prs * pds_c", dp)
  TERMS[[length(TERMS) + 1]] <- terms_of(m5b, cbind(key, model = "pubonly"))
  if (g$method == "sbrc") {
    dq <- d[!is.na(pds_p_c)]
    m5p <- fit("prs * age_c + prs * pds_p_c", dq)
    TERMS[[length(TERMS) + 1]] <- terms_of(m5p, cbind(key, model = "puberty_parent"))
  }
}
terms <- rbindlist(TERMS); vc <- rbindlist(VC); seg <- rbindlist(SEG)
# express rate effects per SD of the child slope (sqrt tau1 of the linear model, same arm)
tau1 <- vc[model == "linear" & grp == "subject" & var1 == "age_c" & is.na(var2), .(disorder, method, arm, sd_slope = sdcor)]
terms <- tau1[terms, on = c("disorder", "method", "arm")]
terms[, est_per_sd_slope := ifelse(grepl("prs:", term) & !grepl("pds", term), est / sd_slope, NA_real_)]
seg <- tau1[seg, on = c("disorder", "method", "arm")]
seg[, `:=`(est_per_sd_slope = est / sd_slope, se_per_sd_slope = se / sd_slope)]
fwrite(terms, file.path(OUT, "onestage_terms.tsv"), sep = "\t")
fwrite(vc, file.path(OUT, "variance_components.tsv"), sep = "\t")
fwrite(seg, file.path(OUT, "prs_rate_by_age.tsv"), sep = "\t")
fwrite(rbindlist(TRAJ), file.path(OUT, "trajectories.tsv"), sep = "\t")

# ---- puberty tempo ---------------------------------------------------------------------------
pl <- fread(file.path(WORK, "pds_long.csv"))[!is.na(pds_mean_y) & !is.na(age_c)]
ch <- fread(file.path(WORK, "children.csv"))
pl <- ch[, .(subject, sex, site, family_id)][pl, on = "subject"]
pl[, `:=`(sex = factor(sex), site = factor(site), subject = factor(subject))]
mp <- lmer(pds_mean_y ~ age_c * sex + (1 + age_c | subject) + (1 | site), data = pl, REML = TRUE, control = ctrl)
vp <- as.data.table(VarCorr(mp)); print(vp)
rp <- ranef(mp)$subject
pds <- data.table(subject = rownames(rp), pds_stage = rp[["(Intercept)"]], pds_tempo = rp[["age_c"]])
pds[, `:=`(pds_stage_z = scale(pds_stage)[, 1], pds_tempo_z = scale(pds_tempo)[, 1])]
cat(sprintf("PDS tempo LMM: %d waves, %d children; r(stage, tempo) = %.2f\n", nobs(mp), nrow(pds), cor(pds$pds_stage, pds$pds_tempo)))
ch <- pds[ch, on = "subject"]
ch[, `:=`(sex = factor(sex), family_id = factor(family_id), site = factor(site),
          slope_z = scale(global_slope_1lmm)[, 1], ct_z = scale(baseline_thickness_1lmm)[, 1])]
PT <- list()
lm_row <- function(f, d, key) {
  m <- lmer(f, data = d, REML = TRUE, control = ctrl)
  s <- as.data.table(coef(summary(m)), keep.rownames = "term"); setnames(s, c("term", "est", "se", "df", "t", "p"))
  cbind(key, s[!term %in% c("(Intercept)", "sexM", "sexF", pcs)], n = nobs(m))
}
for (i in seq_len(nrow(grid))) {
  g <- grid[i]; col <- sprintf("%s_%s_%s", g$disorder, g$method, g$arm)
  d <- copy(ch)[!is.na(get(col)) & !is.na(pds_tempo_z)]; d[, prs := get(col)]
  key <- data.table(disorder = toupper(g$disorder), method = g$method, arm = g$arm)
  for (y in c("pds_stage_z", "pds_tempo_z"))
    PT[[length(PT) + 1]] <- lm_row(as.formula(paste(y, "~ prs +", covs, "+ (1 | family_id)")), d, cbind(key, outcome = y))
  # thinning slope on PRS with and without puberty (does puberty carry the PRS effect?)
  PT[[length(PT) + 1]] <- lm_row(as.formula(paste("slope_z ~ prs +", covs, "+ (1 | family_id)")), d, cbind(key, outcome = "slope_z"))
  PT[[length(PT) + 1]] <- lm_row(as.formula(paste("slope_z ~ prs + pds_stage_z + pds_tempo_z +", covs, "+ (1 | family_id)")), d,
                                 cbind(key, outcome = "slope_z|puberty"))
}
# puberty -> thinning, all genotyped children (no PRS), pooled arm rows
d <- ch[!is.na(pds_tempo_z)]
for (y in c("slope_z", "ct_z"))
  PT[[length(PT) + 1]] <- lm_row(as.formula(paste(y, "~ pds_stage_z + pds_tempo_z +", covs, "+ (1 | family_id)")), d,
                                 cbind(data.table(disorder = "none", method = "none", arm = "all"), outcome = y))
fwrite(rbindlist(PT), file.path(OUT, "puberty_tempo.tsv"), sep = "\t")
fwrite(rbind(vc, cbind(data.table(disorder = "none", method = "none", arm = "all", model = "pds_tempo_lmm"),
                       vp[, .(grp, var1, var2, vcov, sdcor)])),
       file.path(OUT, "variance_components.tsv"), sep = "\t")
cat("done\n")
