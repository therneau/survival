#
# Test residuals for a coxphms model.
# This is done by explicitly fitting the three individual transitions and 
#  computing the residuals on those. 
# Use the same myeloid data set as multi2.R
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

tdata <- tmerge(myeloid[,1:4], myeloid, id=id, death=event(futime,death),
                priortx = tdc(txtime), sct= event(txtime))
tdata$event <- factor(with(tdata, sct + 2*death), 0:2,
                      c("censor", "sct", "death"))

# make the processing a little harder by making sex missing for 3 people,
#  but don't use it as a covariate in the entry to sct transition.
# Obs 428 is a 2:3 transition and sex is missing, so it will end up omitted
#  in fit, listed in na.action, ditto in fit23.  
# The other 3 NA are at risk for entry:death so appear as residuals in fit12,
#  but are part of na.omit in fit13.
tdata$sex[tdata$id %in% 273:275] <- NA   # obs 425 to 428

check <- survcheck(Surv(tstart, tstop, event) ~ 1, tdata, id=id)
fit <- coxph(list(Surv(tstart, tstop, event) ~ trt, 
                  1:3 + 2:3 ~ sex,
                  1:2 + 2:3 ~ flt3), tdata, id=id)
aeq(check$transitions, fit$transitions + c(0,0,0, 0,0,0, 0,1,0))

# this is a difference between coxph, which has the more sophisticated list
#  form of a formula, and survfit/survcheck which do not.
# survcheck with ~1 on the right uses all 1009 obs, 1 more obs than fit, but 
#  with ~sex on the right it will use 1004, 3 less than fit
aeq(fit$n, 1008)

# Multi-state coxph defaults to breslow rather than efron
fit12 <- coxph(Surv(tstart, tstop, event=='sct') ~ trt + flt3, tdata,
               subset=(priortx==0), iter=4, x=TRUE, method='breslow')
fit13 <- coxph(Surv(tstart, tstop, event=='death') ~ trt + sex, tdata,
               subset=(priortx==0), iter=4, x=TRUE, method= 'breslow')
fit23 <- coxph(Surv(tstart, tstop, event=='death') ~ trt + sex+ flt3, tdata,
               subset=(priortx==1), iter=4, x=TRUE, method="breslow")

# martingale residuals
rr1 <- resid(fit)
# One row per obs (except row 428), one col per transtition
# Obs 1-2 start in state (s0) so contribute to the 1:2 and 1:3 transitions,
#  columns 1 and 2 of rr1, while rr1[1:2, 3] =0.
# Obs 3 starts in the sct state and only contributes to the 2:3 transition
all(rr1[1:2,3] ==0)
 
temp1 <- 0* rr1  # on row per retained obs, one col per state
tdat2 <- tdata[-fit$na.action,] # version of data without obs 428
temp1[tdat2$priortx==0, 1] <- resid(fit12)
temp1[tdat2$priortx==0 & !is.na(tdat2$sex), 2] <- resid(fit13)
temp1[tdat2$priortx==1 & !is.na(tdat2$sex), 3] <- resid(fit23)
aeq(rr1, temp1)

# collapsed martingale
aeq(resid(fit, collapse=TRUE), rowsum(temp1, tdat2$id))

# Schoenfeld residuals, which have one row per observed event.  The rows
#  for transition 1:2, 1:3, and 2:3 follow one after the other in a 
#  single matrix
rr2 <- resid(fit, type="schoen")
transition <- attr(rr2, "transition")
temp2 <- 0 * rr2
temp2[transition=="1:2", c(1,3,4)] <- resid(fit12, type="schoen")
temp2[transition=="1:3", 1:2]      <- resid(fit13, type="schoen")
temp2[transition=="2:3", ]         <- resid(fit23, type="schoen")
aeq(rr2, temp2)

# dfbeta residuals
rr3 <- resid(fit, type='dfbeta')
temp3 <- 0*rr3
temp3[tdat2$priortx==0, c(1,3,4),1] <- resid(fit12, type='dfbeta')
temp3[tdat2$priortx==0 &!is.na(tdat2$sex), 1:2, 2] <- 
    resid(fit13, type='dfbeta')
temp3[tdat2$priortx==1 & !is.na(tdat2$sex), 1:4, 3] <- 
    resid(fit23, type='dfbeta')
aeq(rr3, temp3)

# rowsum requires a matrix, fold temp3 into one, then unfold it
temp3b <- rowsum(matrix(temp3, ncol=prod(dim(temp3)[-1])), tdat2$id, 
                 reorder=FALSE)
temp3b <- array(temp3b, dim=c(nrow(temp3b), dim(temp3)[-1]))
aeq(resid(fit, type='dfbeta', collapse=TRUE), temp3b)


# More complex formula
fit2 <- coxph(list(Surv(tstart, tstop, event) ~ trt,
                   1:3 + 2:3 ~ trt + strata(sex),
                   1:2 + 2:3 ~ flt3:sex), tdata, id=id)

