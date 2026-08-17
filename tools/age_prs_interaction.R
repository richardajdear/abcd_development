suppressMessages({library(data.table); library(lme4); library(lmerTest)})
pcs <- paste0("PC", 1:10)
out <- list()
for (dis in c("SCZ","MDD")) {
  d <- fread(sprintf("/tmp/agexprs_%s.csv", dis))
  d[, `:=`(sex=factor(sex), site=factor(site), family_id=factor(family_id),
           subject=factor(subject))]
  # age_c is centred at 12.44. The interaction prs_z:age_c IS the test:
  # main effect = PRS effect on thickness AT 12.44 (the intercept),
  # interaction = PRS effect on the RATE of change (the slope).
  f <- as.formula(paste("value ~ prs_z * age_c + sex +", paste(pcs, collapse=" + "),
                        "+ (1 + age_c | subject) + (1 | site) + (1 | family_id)"))
  m <- lmer(f, data=d, REML=TRUE,
            control=lmerControl(optimizer="bobyqa", calc.derivs=FALSE))
  s <- summary(m)$coefficients
  for (term in c("prs_z","age_c","prs_z:age_c")) {
    out[[length(out)+1]] <- data.table(disorder=dis, term=term,
      est=s[term,"Estimate"], se=s[term,"Std. Error"],
      t=s[term,"t value"], p=s[term,"Pr(>|t|)"])
  }
  cat(sprintf("%s done: n=%d obs, %d subjects\n", dis, nrow(d), uniqueN(d$subject)))
}
res <- rbindlist(out)
fwrite(res, "/tmp/agexprs_results.csv")
print(res)
