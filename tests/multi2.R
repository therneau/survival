library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# Check that estimates from a multi-state model agree with single state models
#  Use a simplified version of the myeloid data set
tdata <- tmerge(myeloid[,1:3], myeloid, id=id, death=event(futime,death),
                priortx = tdc(txtime), sct= event(txtime))
tdata$event <- factor(with(tdata, sct + 2*death), 0:2,
                      c("censor", "sct", "death"))
fit <- coxph(Surv(tstart, tstop, event) ~ trt + sex, tdata, id=id,
             iter=4, x=TRUE, robust=FALSE)

fit12 <- coxph(Surv(tstart, tstop, event=='sct') ~ trt + sex, tdata,
               subset=(priortx==0), iter=4, x=TRUE)
fit13 <- coxph(Surv(tstart, tstop, event=='death') ~ trt + sex, tdata,
               subset=(priortx==0), iter=4, x=TRUE)
fit23 <- coxph(Surv(tstart, tstop, event=='death') ~ trt + sex, tdata,
               subset=(priortx==1), iter=4, x=TRUE)
aeq(coef(fit), c(coef(fit12), coef(fit13), coef(fit23))) 
aeq(fit$loglik, fit12$loglik + fit13$loglik + fit23$loglik)
temp <- matrix(0, 6,6)
temp[1:2, 1:2] <- fit12$var
temp[3:4, 3:4] <- fit13$var
temp[5:6, 5:6] <- fit23$var
aeq(fit$var, temp)

# check out model.frame
fita <- coxph(Surv(tstart, tstop, event) ~ trt, tdata, id=id)
fitb <- coxph(Surv(tstart, tstop, event) ~ trt, tdata, id=id, model=TRUE)
all.equal(model.frame(fita), fitb$model)
# model.frame fails due to an interal rule in R, factors vs characters
#  result when the xlev arg is in the call.  So model.frame(fita) has trt
#  as a factor, not character.

#check residuals
indx1 <- which(fit$rmap[,2] ==1)
indx2 <- which(fit$rmap[,2] ==2)
indx3 <- which(fit$rmap[,2] ==3)
aeq(residuals(fit), c(residuals(fit12), residuals(fit13), residuals(fit23)))
aeq(residuals(fit)[indx1], residuals(fit12))
aeq(residuals(fit)[indx2], residuals(fit13))
aeq(residuals(fit)[indx3], residuals(fit23))

# score residuals
temp <- residuals(fit, type='score')
aeq(temp[indx1, 1:2], residuals(fit12, type='score'))
aeq(temp[indx2, 3:4], residuals(fit13, type='score'))
aeq(temp[indx3, 5:6], residuals(fit23, type='score'))

all(temp[indx1, 3:6] ==0)
all(temp[indx2, c(1,2,5,6)] ==0)
all(temp[indx3, 1:4]==0)

temp <- residuals(fit, type="dfbeta")
all(temp[indx1, 3:6] ==0)
all(temp[indx2, c(1,2,5,6)] ==0)
all(temp[indx3, 1:4]==0)
aeq(temp[indx1, 1:2], residuals(fit12, type='dfbeta'))
aeq(temp[indx2, 3:4], residuals(fit13, type='dfbeta'))
aeq(temp[indx3, 5:6], residuals(fit23, type='dfbeta'))

temp <- residuals(fit, type="dfbetas")
all(temp[indx1, 3:6] ==0)
all(temp[indx2, c(1,2,5,6)] ==0)
all(temp[indx3, 1:4]==0)
aeq(temp[indx1, 1:2], residuals(fit12, type='dfbetas'))
aeq(temp[indx2, 3:4], residuals(fit13, type='dfbetas'))
aeq(temp[indx3, 5:6], residuals(fit23, type='dfbetas'))

# Schoenfeld and scaled shoenfeld have one row per event
sr1 <- residuals(fit12, type="schoenfeld")
sr2 <- residuals(fit13, type="schoenfeld")
sr3 <- residuals(fit23, type="schoenfeld")
end <- rep(1:3, c(nrow(sr1), nrow(sr2), nrow(sr3)))
temp <- residuals(fit, type="schoenfeld")
aeq(temp[end==1, 1:2], sr1)
aeq(temp[end==2, 3:4], sr2)
aeq(temp[end==3, 5:6], sr3)
all(temp[end==1, 3:6] ==0)
all(temp[end==2, c(1,2,5,6)] ==0)
all(temp[end==3, 1:4] ==0)


#The scaled Schoenfeld don't agree, due to the use of a robust
#  variance in fit, regular variance in fit12, fit13 and fit23
#Along with being scaled by different event counts
xfit <- fit
xfit$var <- xfit$naive.var
if (FALSE) {
    xfit <- fit
    xfit$var <- xfit$naive.var  # fixes the first issue
    temp <- residuals(xfit, type="scaledsch")
    aeq(d1* temp[sindx1, 1:2], residuals(fit12, type='scaledsch'))
    aeq(temp[sindx2, 3:4], residuals(fit13, type='scaledsch'))
    aeq(temp[sindx3, 5:6], residuals(fit23, type='scaledsch'))
}

if (FALSE) { # the predicted values are a work in progress
# predicted values differ because of different centering
c0 <-  sum(fit$mean * coef(fit))
c12 <- sum(fit12$mean * coef(fit12))
c13 <- sum(fit13$mean* coef(fit13))
c23 <- sum(fit23$mean * coef(fit23))

aeq(predict(fit)+c0, c(predict(fit12)+c12, predict(fit13)+c13, 
                       predict(fit23)+c23))
aeq(exp(predict(fit)), predict(fit, type='risk'))

# expected survival is independent of centering
aeq(predict(fit, type="expected"), c(predict(fit12, type="expected"),
                                     predict(fit13, type="expected"),
                                     predict(fit23, type="expected")))
}
# predict(type='terms') is a matrix, centering changes as well
if (FALSE) {
    temp <- predict(fit, type='terms')
    all(temp[indx1, 3:6] ==0)
    all(temp[indx2, c(1,2,5,6)] ==0)
    all(temp[indx3, 1:4]==0)
    aeq(temp[indx1, 1:2], predict(fit12, type='terms'))
    aeq(temp[indx2, 3:4], predict(fit13, type='terms'))
    aeq(temp[indx3, 5:6], predict(fit23, type='terms'))
} # end of prediction section

# The global and per strata zph tests will differ for the KM or rank
#  transform, because the overall and subset will have a different list
#  of event times, which changes the transformed value for all of them.
# But identity and log are testable.
test_a <- cox.zph(fit, transform="log",global=FALSE)
test_a12 <- cox.zph(fit12, transform="log",global=FALSE)
test_a13 <- cox.zph(fit13, transform="log", global=FALSE)
test_a23 <-  cox.zph(fit23, transform="log", global=FALSE)
aeq(test_a$y[test_a$strata==1, 1:2], test_a12$y)

aeq(test_a$table[1:2,], test_a12$table)
aeq(test_a$table[3:4,], test_a13$table)
aeq(test_a$table[5:6,], test_a23$table)

# check cox.zph fit - transform = 'identity'
test_b <- cox.zph(fit, transform="identity",global=FALSE)
test_b12 <- cox.zph(fit12, transform="identity",global=FALSE)
test_b13 <- cox.zph(fit13, transform="identity", global=FALSE)
test_b23 <-  cox.zph(fit23, transform="identity", global=FALSE)

aeq(test_b$table[1:2,], test_b12$table)
aeq(test_b$table[3:4,], test_b13$table)
aeq(test_b$table[5:6,], test_b23$table)

# check out subscripting of a multi-state zph
cname <- c("table", "x", "time", "y", "var")
sapply(cname, function(x) aeq(test_b[1:2]$x, test_b12$x))
sapply(cname, function(x) aeq(test_b[3:4]$x, test_b13$x))
sapply(cname, function(x) aeq(test_b[5:6]$x, test_b23$x))

# check model.matrix
mat1 <- model.matrix(fit)
all.equal(mat1, fit$x)

# Check that the internal matix agrees (uses stacker, which is not exported)
mat2 <- model.matrix(fit12)
mat3 <- model.matrix(fit13)
mat4 <- model.matrix(fit23)

# first reconstruct istate
tcheck <- survcheck(Surv(tstart, tstop, event) ~ 1, tdata, id=id)
temp <- survival:::stacker(fit$cmap, fit$smap, as.numeric(tcheck$istate), fit$x,
                          fit$y, NULL, fit$states)
aeq(temp$X[temp$transition==1, 1:2], mat2)
aeq(temp$X[temp$transition==2, 3:4], mat3)
aeq(temp$X[temp$transition==3, 5:6], mat4)


