library(survival)

# Various edge cases for concordance
fit1 <- coxph(Surv(time, status) ~ age + offset(ph.ecog*0) +strata(sex), lung)
fit2 <- coxph(Surv(time, status) ~ age + ph.ecog +strata(sex), lung)

test <- concordance(fit1, fit2)
