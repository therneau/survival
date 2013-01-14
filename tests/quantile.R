#
# Formal test of the quantile routine for survfit
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# There are 8 cases:  strata Y/N,  ncol(surv) >1,  conf.int = T/F
#  Subcase: the quantile exactly agrees with a horizontal segment of 
#    the curve or not.
# First do the 4 cases where fit$surv is a vector
#
test1 <- data.frame(time=  c(9, 3,1,1,6,6,8, 10),
                    status=c(1,NA,1,0,1,1,0,  0),
                    x=     c(0, 2,1,1,1,0,0,  0))

# True survival = (6/7) * (3/5) * (1/2) for overall
#   The q's are chosen to include a point < first jump, mid, after last jump,
#   and exact intersections with the "flats" of the curve.
#   
qq <- c(13/14, 6/7, 2/3, .5, 9/35, .1)
 
# Nothing on the right hand side, simple survival (no strata)
fit1 <- survfit(Surv(time, status) ~ 1, test1, conf.type='none')
aeq(quantile(fit1, 1-qq), c(1, 3.5, 6, 9, 9.5, NA))  #without conf.int

fit2 <-  survfit(Surv(time, status) ~ 1, test1)  #with conf.int
aeq(quantile(fit2, 1-qq), 
    list(quantile = c(1, 3.5, 6, 9, 9.5, NA),
         lower = c(1,1,1,6,6,9),
         upper = rep(as.numeric(NA), 6)), check.attributes=FALSE)
aeq(quantile(fit2, 1-qq, FALSE), c(1, 3.5, 6, 9, 9.5, NA)) 


# Now a variable on the right (strata in the result)
#  curve 0: (t=6, S=3/4),  (t=9, S=3/8)
#  curve 1: (t=1, S=2/3),  (t=6, S= 0)
fit1 <- survfit(Surv(time, status) ~ x, test1, conf.type='none') 
aeq(quantile(fit1, 1-qq),
    matrix(c(6,6,9,9,NA,NA,  1,1,3.5, 6,6,6), nrow=2, byrow=T))

fit2 <- survfit(Surv(time, status) ~ x, test1)
aeq(quantile(fit2, 1-qq, FALSE),
    matrix(c(6,6,9,9,NA,NA,  1,1,3.5, 6,6,6), nrow=2, byrow=T))

temp <- quantile(fit2, 1-qq)
aeq(temp$quantile, matrix(c(6,6,9,9,NA,NA,  1,1,3.5, 6,6,6), nrow=2, byrow=T))
aeq(temp$lower,    matrix(c(6,6,6,6,9,9,    1,1,1,1, NA,NA), nrow=2, byrow=T))
aeq(temp$upper,    rep(as.numeric(NA), 12))

# Second major case set -- a survfit object where fit$surv is a matrix
#  This arises from coxph models
#  There is only 1 subject with ph.ecog=3 which is a nice edge case
cfit <- coxph(Surv(time, status) ~ age + strata(ph.ecog), lung)
sfit <- survfit(cfit, newdata=data.frame(age=c(50, 70)))
qtot <- quantile(sfit, qq)
for (i in 1:4) {
    for (j in 1:2) {
        temp <- quantile(sfit[i,j], qq)
        print(c(aeq(qtot$quantile[i,j,], temp$quantile),
                aeq(qtot$upper[i,j,], temp$upper),
                aeq(qtot$lower[i,j,], temp$lower)))
    }
}
temp <- quantile(sfit, qq, conf.int=FALSE)
all.equal(qtot$quantile, temp)

#
# Third case -- a survfitms object, which results from cumulative
#  incidence curves.
#
tdata <- data.frame(time=c(1,2,2,3,3,3,5,6),
                    status = c(0,1,0,1,0,1,0,1),
                    event =  c(1,1,2,2,1,2,3,2),
                    grp = c(1,2,1,2,1,2,1,2))

fit1 <- survfit(Surv(time, status*event, type='mstate') ~1, tdata)
temp <- quantile(fit1, c(.1, .2, .5))
aeq(temp$quantile, matrix(c(2, NA, NA, 3,3,6), nrow=2, byrow=TRUE))
aeq(temp$lower   , matrix(c(2,  2, NA, 3,3,3), nrow=2, byrow=TRUE))
aeq(temp$upper   , c(NA,6, rep(NA,4)))

fit2 <- survfit(Surv(time, status*event, type='mstate') ~1, tdata, 
                conf.int=FALSE)
temp <- quantile(fit2, c(.1, .2, .5))
aeq(temp, matrix(c(2, NA, NA, 3,3,6), nrow=2, byrow=TRUE))

# Use a larger data set for the multi-group + multi-column case, the MGUS data
#  However, it has almost no censoring, so add a little to make the
#  quantiles not be exactly even percentiles
mdata <- data.frame(time=mgus1$stop,
                    status=mgus1$status,
                    event= mgus1$event,
                    sex=mgus1$sex,
                    stat2= factor(ifelse(mgus1$status==0, 0, 
                                         as.numeric(mgus1$event)),
                                  levels=0:2, 
                                  labels=c("censor", levels(mgus1$event)))
                    )[mgus1$start==0,]
mdata$stat2[seq(1, nrow(mdata), by=5)] <- "censor"

fit3 <- survfit(Surv(time, stat2) ~sex, mdata)
temp1 <- quantile(fit3, 0:10/20)
temp2 <- quantile(fit3, 0:10/20, conf.int=FALSE)
aeq(temp1$quantile, temp2)

for (i in 1:2) {
    for (j in 1:2){
        temp3 <- quantile(fit3[i,j], 0:10/20)
        print(c(aeq(temp1$quantile[i,j,], temp3$quantile),
                aeq(temp1$upper[i,j,], temp3$upper),
                aeq(temp1$lower[i,j,], temp3$lower)))
    }
}

# Do one set of quantiles by brute force        
zz <- 1:fit3$strata[1]
temp3 <- double(10)
tt <- fit3$time[zz]
for (i in 1:10) temp3[i] <- min(tt[fit3$prev[zz,2] > i/20])
aeq(temp3, temp2[1,2,2:11])



