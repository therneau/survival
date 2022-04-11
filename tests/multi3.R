library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# Check that a multi-state model, correctly set up, gives the same
# solution as a time-dependent covariate.
# This is a stronger test than mstrata: there the covariate which was mapped
#  into a state was constant, here it is time-dependent.
#
# First build the TD data set from pbcseq, with a categorical bilirubin
pbc1 <- pbcseq
pbc1$bili4 <- cut(pbc1$bili, c(0,1, 2,4, 100), 
                  c("normal", "1-2x", "2-4x", ">4"))
ptemp <- subset(pbc1, !duplicated(id))  # first row of each

pbc2 <- tmerge(ptemp[, c("id", "age", "sex")], ptemp, id,
               death= event(futime, status==2))

pbc2 <- tmerge(pbc2, pbc1, id=id, bili = tdc(day, bili),
                 bili4 = tdc(day, bili4), bstat = event(day, as.numeric(bili4)))
btemp <- with(pbc2, ifelse(death, 5, bstat))

# a row with the same starting and ending bili4 level is not an event
b2 <- ifelse(((as.numeric(pbc2$bili4)) == btemp), 0, btemp)
pbc2$bstat <- factor(b2, 0:5,
                     c("censor", "normal", "1-2x", "2-4x", ">4", "death"))
check1 <- survcheck(Surv(tstart, tstop, bstat) ~ 1, istate= bili4,
                    id = id, data=pbc2)
check1$transitions
all.equal(as.character(pbc2$bili4), as.character(check1$istate))
# the above verifies that I created the data set correctly

# Standard coxph fit with a time dependent bili4 variable.
fit1 <- coxph(Surv(tstart, tstop, death) ~ age + bili4, pbc2)

# An additive multi-state fit, where bili4 is a state
#  The three forms below should all give identical models
fit2 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   c(1:4):5 ~ age / common + shared), id= id, istate=bili4,
              data=pbc2)
fit2b <- coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   1:5 + 2:5 + 3:5 + 4:5 ~ age / common + shared), 
              id= id, istate=bili4, data=pbc2)
fit2c <- coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   0:5 ~ age / common + shared), 
              id= id, istate=bili4, data=pbc2)

# Make sure the names are correct and the coefficients match
aeq(coef(fit1), coef(fit2))
aeq(names(coef(fit2)), c("age", "ph(2:5/1:5)", "ph(3:5/1:5)", "ph(4:5/1:5)"))
all.equal(coef(fit2), coef(fit2b))
all.equal(coef(fit2), coef(fit2c))

# Now a model with a separate age effect for each bilirubin group
fit3  <- coxph(Surv(tstart, tstop, death) ~ age*bili4, pbc2)
fit3b <- coxph(Surv(tstart, tstop, death) ~ bili4/age, pbc2)
fit4 <-  coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   c(1:4):5 ~ age / shared), id= id, istate=bili4,
              data=pbc2)
all.equal(fit3$loglik, fit3b$loglik)
all.equal(fit3$loglik, fit4$loglik)

# The coefficients are quite different due to different codings for dummy vars
# Unpack the interaction, first 4 coefs will be the age effect within each
#  bilirubin group
temp <- c(coef(fit3)[1] + c(0, coef(fit3)[5:7]), coef(fit3)[2:4])
names(temp)[1:4] <- c("age1", "age2", "age3", "age4")
aeq(temp, coef(fit3b)[c(4:7, 1:3)])
aeq(temp, coef(fit4))

# Third, a model with separate baseline hazards for each bili group
fit5 <- coxph(Surv(tstart, tstop, death) ~ strata(bili4)/age, pbc2,
              cluster=id)
fit6 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1, 0:5 ~ age),
                   id=id, istate=bili4, pbc2)
aeq(coef(fit5), coef(fit6))
aeq(fit5$var, fit6$var)
aeq(fit5$naive.var, fit6$naive.var)
