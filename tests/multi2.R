library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

#
# Check that estimates agree with single state models
#  Use a simplified version of the myeloid data set
tdata <- tmerge(myeloid[,1:2], myeloid, id=id, death=event(futime,death),
                priortx = tdc(txtime), sct= event(txtime))
tdata$event <- factor(with(tdata, sct + 2*death), 0:2, 
                      c("censor", "sct", "death"))
fit <- coxph(Surv(tstart, tstop, event) ~ trt, tdata, id=id, 
             iter=3, x=TRUE)

fit12 <- coxph(Surv(tstart, tstop, event=='sct') ~ trt, tdata,
               subset=(priortx==0))
fit13 <- coxph(Surv(tstart, tstop, event=='death') ~ trt, tdata,
               subset=(priortx==0))
fit23 <- coxph(Surv(tstart, tstop, event=='death') ~ trt, tdata,
               subset=(priortx==1))
aeq(coef(fit), c(coef(fit12), coef(fit13), coef(fit23)))
aeq(fit$loglik, fit12$loglik + fit13$loglik + fit23$loglik)
aeq(fit$naive.var, diag(c(fit12$var, fit13$var, fit23$var)))

ii <- fit$strata==1
tfit <- coxph(fit$y[ii,] ~ fit$x[ii,])
aeq(tfit$loglik, fit12$loglik)   # check that x, y, strata are correct
ii <- fit$strata==2
tfit <- coxph(fit$y[ii,] ~ fit$x[ii,])
aeq(tfit$loglik, fit13$loglik)   # check that x, y, strata are correct
ii <- fit$strata==3
tfit <- coxph(fit$y[ii,] ~ fit$x[ii,])
aeq(tfit$loglik, fit23$loglik)   # check that x, y, strata are correct

# check out model.matrix and model.frame
fita <- coxph(Surv(tstart, tstop, event) ~ trt, tdata, id=id)
fitb <- coxph(Surv(tstart, tstop, event) ~ trt, tdata, id=id, model=TRUE)
all.equal(model.frame(fita), fitb$model)

fitc <- coxph(Surv(tstart, tstop, event) ~ trt, tdata, id=id, x=TRUE)

