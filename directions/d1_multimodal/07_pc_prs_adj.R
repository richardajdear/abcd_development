#!/usr/bin/env Rscript
# D1 step 5c -- PRS on diffusion PC2 / PC4 (05_pls_composite.py scores), with and without the thinning slope as a
# covariate (Figure-1 model + site), plus the correlation of every diffusion PC with the thinning slope.
# Output: tables/pc_prs_adj_thin.tsv. Run (env r, repo root).
suppressMessages({library(data.table); library(lme4); library(lmerTest)})
W <- "directions/d1_multimodal/work"; OUT <- "directions/d1_multimodal/tables"
B <- fread(file.path(W, "blups_dmri.csv.gz")); thin <- B[metric == "smri__thk_hcp", .(subject, thin_fast = -slope)]
pc <- fread(file.path(W, "pls_child.csv.gz"))[, .(subject, dPC1 = PC1, dPC2 = PC2, dPC3 = PC3, dPC4 = PC4, dPC5 = PC5)]
x <- merge(pc, thin, by = "subject")
cat("r(dPC_k, faster thinning):", paste(sprintf("dPC%d %+.3f", 1:5, sapply(1:5, function(k) cor(x[[paste0("dPC", k)]], x$thin_fast, use = "complete.obs"))), collapse = "; "), "\n")
ch <- fread("genetic_analysis/work/results_70tab_hcp/d9_sex_moderation/children.csv")
ch[, `:=`(sex = factor(sex, levels = c("F", "M")), FID = factor(FID), site = factor(site), age_c = baseline_age)]
ch <- x[ch, on = "subject"]
cov0 <- paste(c("sex", "age_c", paste0("PC", 1:10), "site"), collapse = " + ")
R <- list()
for (yv in c("dPC2", "dPC4")) for (dis in c("scz", "mdd")) for (meth in c("prscs", "sbrc")) for (arm in c("eur", "pooled")) {
  col <- sprintf("%s_%s_%s", dis, meth, arm)
  d <- ch[!is.na(get(col)) & !is.na(get(yv)) & !is.na(PC1) & !is.na(thin_fast)]
  if (arm == "eur") d <- d[eur == 1]
  d[, `:=`(prs = as.numeric(scale(get(col))), y = as.numeric(scale(get(yv))), tz = as.numeric(scale(thin_fast)))]
  for (adj in c("site", "site+thin")) {
    f <- paste("y ~ prs +", cov0, if (adj == "site+thin") "+ tz" else "", "+ (1 | FID)")
    s <- coef(summary(lmer(as.formula(f), d, control = lmerControl(calc.derivs = FALSE))))
    R[[length(R) + 1]] <- data.table(outcome = yv, disorder = toupper(dis), method = meth, arm, model = adj,
                                     est = s["prs", 1], p = s["prs", 5], n = nrow(d))
  }
}
R <- rbindlist(R); fwrite(R, file.path(OUT, "pc_prs_adj_thin.tsv"), sep = "\t")
R[, cell := sprintf("%+.3f (p %.3f)", est, p)]
print(dcast(R, outcome + disorder + method + arm ~ model, value.var = "cell"))
