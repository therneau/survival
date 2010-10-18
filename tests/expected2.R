library(survival)
#
# A Cox model with a factor, followed by survexp.  
#
pfit2 <- coxph(Surv(time, status > 0) ~ trt + log(bili) +
          log(protime) + age + platelet + sex, data = pbc)
esurv <- survexp(~ trt, ratetable = pfit2, data = pbc)
