options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Fit models to the Colon cancer data used in Lin
#
fitc1 <- coxph(Surv(time, status) ~ rx + extent + node4 + cluster(id)
	        + strata(etype), colon)
fitc1

fitc2 <- coxph(Surv(time, status) ~ rx + extent + node4 + 
	       frailty(id, dist='gauss', trace=T)
	        + strata(etype), colon)
fitc2

fitc3 <- coxph(Surv(time, status) ~ rx + extent + node4 + frailty(id, trace=T)
	        + strata(etype), colon)
fitc3

fitc4 <- coxph(Surv(time, status) ~ rx + extent + node4 + frailty(id, df=30)
	        + strata(etype), colon)
fitc4

# Do a fit, removing the no-event people
temp <- tapply(colon$status, colon$id, sum)
keep <- !(is.na(match(colon$id, names(temp[temp>0])))) 
fitc5 <- coxph(Surv(time, status) ~ rx + extent + node4 +cluster(id)
	       + strata(etype), colon, subset=keep)
fitc5

#
# Do the factor fit, but first remove the no-event people
#
#  Ha!  This routine has a factor with 506 levels.  It uses all available
#    memory, and can't finish in my patience window.  Commented out.

#fitc6 <- coxph(Surv(time, status) ~ rx + extent + node4 + factor(id), colon,
#	       subset=keep)

# Try the AIC
fitc7 <- coxph(Surv(time, status) ~ rx + extent + node4 +
	          frailty(id, method='aic', caic=T), colon)
fitc7

if (FALSE) {
    # This fails to converge, with different bad paths in different
    #  releases (due to changes in eps parameters).  It's not a bug per
    #  se, just too unstable a solution for the code to find.  
    # Remove it from the tests 2/2/09
    fitc8 <- coxph(Surv(time, status) ~ rx + extent + node4 +
	          frailty(id, method='aic', dist='t', caic=T), colon)

    fitc8
    }



