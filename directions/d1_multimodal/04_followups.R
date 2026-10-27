#!/usr/bin/env Rscript
# D1 step 4 -- two follow-ups on 03: (a) are the PRS effects on grey-matter diffusion slopes the thinning effect
# re-measured? Figure-1 model + site, with and without the child's thinning slope as a covariate.
# (b) how many independent axes do the 47 child slopes span? PCA of site-residualised, z-scored slopes.
# Outputs: tables/dmri_prs_adj_thin.tsv, tables/slope_pca.tsv. Run (env r, repo root).
suppressMessages({library(data.table); library(lme4); library(lmerTest)})
W <- "directions/d1_multimodal/work"; OUT <- "directions/d1_multimodal/tables"
B <- fread(file.path(W, "blups_dmri.csv.gz"))
thin <- B[metric == "smri__thk_hcp", .(subject, thin_fast = -slope)]
ch <- fread("genetic_analysis/work/results_70tab_hcp/d9_sex_moderation/children.csv")
ch[, `:=`(sex = factor(sex, levels = c("F", "M")), FID = factor(FID), site = factor(site), age_c = baseline_age)]
ch <- thin[ch, on = "subject"]
cov0 <- paste(c("sex", "age_c", paste0("PC", 1:10), "site"), collapse = " + ")
MET <- c("dti__is__md__gm", "dti__is__td__gm", "rsi__fni__gm", "dti__is__fa__gwc", "rsi__rnd__gwc", "rsi__hnd__gwc", "rsi__hni__gm")
R <- list()
for (mt in MET) for (dis in c("scz", "mdd")) for (arm in c("eur", "pooled")) {
  col <- sprintf("%s_prscs_%s", dis, arm)
  d <- B[metric == mt, .(subject, yv = slope)][ch, on = "subject", nomatch = 0][!is.na(get(col)) & !is.na(PC1) & !is.na(thin_fast)]
  if (arm == "eur") d <- d[eur == 1]
  d[, `:=`(prs = as.numeric(scale(get(col))), y = as.numeric(scale(yv)), tz = as.numeric(scale(thin_fast)))]
  for (adj in c("site", "site+thin")) {
    f <- paste("y ~ prs +", cov0, if (adj == "site+thin") "+ tz" else "", "+ (1 | FID)")
    s <- coef(summary(lmer(as.formula(f), d, control = lmerControl(calc.derivs = FALSE))))
    R[[length(R) + 1]] <- data.table(metric = mt, disorder = toupper(dis), arm, model = adj, est = s["prs", 1], se = s["prs", 2],
                                     p = s["prs", 5], b_thin = if (adj == "site+thin") s["tz", 1] else NA_real_, n = nrow(d))
  }
}
R <- rbindlist(R); fwrite(R, file.path(OUT, "dmri_prs_adj_thin.tsv"), sep = "\t")
print(dcast(R, metric + disorder + arm ~ model, value.var = c("est", "p")), digits = 2)
# (b) PCA
ws <- dcast(B, subject ~ metric, value.var = "slope"); site1 <- unique(B[, .(subject, site)], by = "subject")
ws <- site1[ws, on = "subject"]; ws <- ws[complete.cases(ws)]
mets <- setdiff(names(ws), c("subject", "site"))
X <- sapply(mets, function(m) scale(resid(lm(ws[[m]] ~ factor(ws$site))))[, 1])
pc <- prcomp(X); ve <- pc$sdev^2 / sum(pc$sdev^2)
cat(sprintf("\nPCA of %d site-residualised slopes on %d children: variance explained %s\n", length(mets), nrow(X),
            paste(sprintf("%.2f", ve[1:8]), collapse = " ")))
L <- data.table(metric = mets, pc$rotation[, 1:5])
fwrite(cbind(L, var_explained = NA), file.path(OUT, "slope_pca.tsv"), sep = "\t")
for (k in 1:4) { o <- L[order(-abs(get(paste0("PC", k))))][1:7]; cat(sprintf("PC%d top: %s\n", k, paste(sprintf("%s %+.2f", o$metric, o[[paste0("PC", k)]]), collapse = "; "))) }
cat(sprintf("thickness loading on PC1-5: %s\n", paste(sprintf("%+.2f", unlist(L[metric == "smri__thk_hcp", 2:6])), collapse = " ")))
