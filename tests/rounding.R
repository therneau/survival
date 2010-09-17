library(survival)
#
# Survival curves could fail with data that was almost exact.
#  The calculations use both unique() and table(), which don't 
#  necessarily give the same number of values.
# Check that the routine handles this properly
#  

tdata <- data.frame(time=c(1,2, sqrt(2)^2, 2, sqrt(2)^2),
                    status=rep(1,5), 
                    group=c(1,1,1,2,2))
fit <- survfit(Surv(time, status) ~ group, data=tdata)

all.equal(sum(fit$strata), length(fit$time))
