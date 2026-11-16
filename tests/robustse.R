library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)
#
# Single-state survfit with a robust variance.  The infinitesimal jackknife
#  se of S is sqrt(colSums(U^2)), U = influence.surv.  Whatever scale the
#  std.err component uses (the logse flag), summary.survfit reports se(S)
#  and the log-scale confidence limit is S exp(-z se(S)/S).
#
zz <- qnorm(.975)
for (stype in 1:2) {
    for (ctype in 1:2) {
        fit <- survfit(Surv(time, status) ~ 1, data=aml, stype=stype,
                       ctype=ctype, id=1:nrow(aml), robust=TRUE, influence=TRUE)
        seS <- sqrt(colSums(fit$influence.surv^2))
        keep <- (fit$surv > 0)
        print(c(stype=stype, ctype=ctype, logse=fit$logse))
        print(aeq(fit$std.err[keep], 
                  (if (stype==1) seS else seS/fit$surv)[keep]))
        print(aeq(summary(fit)$std.err, seS[fit$n.event > 0]))
        print(aeq(fit$lower[keep], 
                  (fit$surv * exp(-zz * seS/fit$surv))[keep]))
    }
}
