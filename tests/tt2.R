# A reprise of tt.R, using (time1, time2) data.
library(survival)
library(splines)
aeq <- function(x, y) all.equal(as.vector(x), as.vector(y))

# A contrived example for the tt function
#
mkdata <- function(n, beta) {
    age <- round(runif(n, 20, 60))
    x <- rbinom(n, 1, .5)

    futime <- rep(40, n)   # everyone has 40 years of follow-up
    entry  <- pmax(0, seq(-10, 30, length=n))  # 1/4 enter at 0
    entry  <- round(entry)
    status <- rep(0, n)
    dtime <-  runif(n/2, 1, 40)  # 1/2 of them die
    dtime <- sort(dtime)

    # The risk is set to beta[1]*x + beta[2]* f(current_age)
    #   where f= 0 up to age 40, rises linear to age 70, flat after that
    for (i in 1:length(dtime)) {
        atrisk <- (futime >= dtime[i] & entry < dtime[i])
        c.age <- age + dtime
        age2 <- pmin(30, pmax(0, c.age-40))
        xbeta <- beta[1]*x + beta[2]*age2
        
        # Select a death according to risk
        risk <- ifelse(atrisk, exp(xbeta), 0)
        dead <- sample(1:n, 1, prob=risk/sum(risk))
        
        futime[dead] <- dtime[i]
        status[dead] <- 1
    }
    out <- data.frame(time1= entry, time2=round(futime,1), status=status, 
                      age=age, x=x, risk=risk,
               casewt = sample(1:5, n, replace=TRUE),
               grp = sample(1:15, n, replace=TRUE), id= 1:n)
    subset(out, time1 < time2)
}

set.seed(1953)  # a good year
# Make n larger for the (time1, time2) case; more stress.
tdata <- mkdata(250, c(log(1.5), 2/30))   # data set has many ties
#tdata <- mkdata(100, c(log(1.5), 2/30))   # data set has many ties
tdata$strat <- floor(tdata$grp/10)

dtime <- sort(unique(tdata$time2[tdata$status==1]))
data2 <- survSplit(Surv(time1, time2, status) ~., tdata, cut=dtime)
data2$c.age <- data2$age + data2$time2  # current age

# fit1 uses data at the event times, fit2$c.age might have a 
#  wider range due to censorings.  To make the two fits agree
#  fix the knots.  I know a priori that 20 to 101 will cover it.
ns2 <- function(x) ns(x, Boundary.knots=c(20, 101), knots=c(45, 60, 75))

fit1 <- coxph(Surv(time1, time2, status)~ x + tt(age), tdata,
              tt= function(x, t, ...) ns2(x+t))

fit2 <- coxph(Surv(time1, time2, status) ~ x + ns2(c.age), data2)

aeq(coef(fit1), coef(fit2))
aeq(vcov(fit1), vcov(fit2))

#
# Check that cluster, weight, and offset were correctly expanded
#
fit3a <- coxph(Surv(time1, time2, status)~ x + tt(age), tdata, weights=casewt, 
              tt= function(x, t, ...) ns2(x+t), x=TRUE)
fit3b <-  coxph(Surv(time1, time2, status) ~ x + ns2(c.age), data2,
                weights=casewt)
aeq(coef(fit3a), coef(fit3b))
aeq(vcov(fit3a), vcov(fit3b))

fit4a <- coxph(Surv(time1, time2, status)~ x + tt(age), tdata,
              tt= function(x, t, ...) ns2(x+t), cluster=grp)
fit4b <-  coxph(Surv(time1, time2, status) ~ x + ns2(c.age), data2,
                cluster=grp)
fit4c <- coxph(Surv(time1, time2, status) ~ x + ns2(c.age) + cluster(grp),
               data2)
aeq(coef(fit4a), coef(fit4b))
aeq(vcov(fit4a), vcov(fit4b))
aeq(coef(fit4a), coef(fit4c))
aeq(vcov(fit4a), vcov(fit4c))

fit5a <- coxph(Surv(time1, time2, status)~ x + tt(age) + offset(grp/10), tdata,
              tt= function(x, t, ...) ns2(x+t),)
fit5b <-  coxph(Surv(time1, time2, status) ~ x + ns2(c.age)+ offset(grp/10)
                , data=data2)
aeq(coef(fit5a), coef(fit5b))
aeq(vcov(fit5a), vcov(fit5b))

# Check that strata is correct
fit6a <- coxph(Surv(time1, time2, status) ~ x + tt(age) + strata(strat), tdata,
               tt = function(x, t, ...) (x+t)^2, x=TRUE)
fit6b <- coxph(Surv(time1, time2, status) ~ x + I(c.age^2) +strata(strat), data2)
aeq(coef(fit6a), coef(fit6b))
aeq(vcov(fit6a), vcov(fit6b))
