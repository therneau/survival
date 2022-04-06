library(survival)
#
# Check that the survival curves from a Cox model with beta=0 
#  match ordinary survival
#
#  Aalen
surv1 <- survfit(Surv(time,status) ~ sex, data=lung, stype=2)
fit1 <- coxph(Surv(time, status) ~ age + strata(sex), data=lung, iter=0,
              ties='breslow')
fit1$var <- 0*fit1$var   #sneaky, causes the extra term in the Cox variance
                         # calculation to be zero
surv2 <- survfit(fit1, stype=2)
surv3 <- survfit(fit1)

arglist <- c('n', 'time', 'n.risk','n.event', 'n.censor', 'surv', 'strata',
             'std.err', 'upper', 'lower')
all.equal(unclass(surv1)[arglist], unclass(surv2)[arglist])
all.equal(unclass(surv1)[arglist], unclass(surv3)[arglist])


# Efron method
surv1 <- survfit(Surv(time,status) ~ sex, data=lung, stype=2, ctype=2)
surv2 <- survfit(fit1, ctype=2)
all.equal(unclass(surv1)[arglist], unclass(surv2)[arglist])

# Kaplan-Meier
surv1 <- survfit(Surv(time,status) ~ sex, data=lung)
surv2 <- survfit(fit1, stype=1)
all.equal(unclass(surv1)[arglist], unclass(surv2)[arglist])


# Now add some random weights
rwt <- runif(nrow(lung), .5, 3)
surv1 <- survfit(Surv(time,status) ~ sex, data=lung, stype=2, weight=rwt,
                 robust=FALSE)
fit1 <- coxph(Surv(time, status) ~ age + strata(sex), data=lung, iter=0,
              ties='breslow', weight=rwt, robust=FALSE)
fit1$var <- 0*fit1$var   #sneaky
surv2 <- survfit(fit1, stype=2, ctype=1)
surv3 <- survfit(fit1)

all.equal(unclass(surv1)[arglist], unclass(surv2)[arglist])
all.equal(unclass(surv1)[arglist], unclass(surv3)[arglist])


# Efron method
surv1 <- survfit(Surv(time,status) ~ sex, data=lung, stype=2, ctype=2,
                 weight=rwt, robust=FALSE)
surv2 <- survfit(fit1, ctype=2, robust=FALSE)
all.equal(unclass(surv1)[arglist], unclass(surv2)[arglist])

# Kaplan-Meier
surv1 <- survfit(Surv(time,status) ~ sex, data=lung, weight=rwt, robust=FALSE)
surv2 <- survfit(fit1, stype=1, robust=FALSE)
all.equal(unclass(surv1)[arglist], unclass(surv2)[arglist])

