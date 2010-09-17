# $Id: frailty.gammacon.S 11166 2008-11-24 22:10:34Z therneau $
# Correct the loglik for a gamma frailty
#  Term2 is the hard one, discussed in section 3.5 of the report
# The penalty function only adds \vu \sum(w_j) to the CoxPL, so this
#  does a bit more than equation 15.
#
frailty.gammacon <- function(d, nu) {
    nfrail <- length(d)
    maxd <- max(d)
    if (nu > 1e7*maxd) term1 <- sum(d*d)/nu  #second order Taylor series
    else               term1 <- sum(d + nu*log(nu/(nu+d)))  #easy part
   
    tbl <- table(factor(d[d>0], levels=1:maxd))
    ctbl<- rev(cumsum(rev(tbl)))   
    dlev<- 1:maxd
    term2.numerator <- nu + rep(dlev-1, ctbl)
    term2.denom     <- nu + rep(dlev, tbl*dlev)
    term2 <- sum(log(term2.numerator/term2.denom))

    term1 + term2
    }
   
