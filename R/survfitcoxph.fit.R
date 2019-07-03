# The survfitcoxph.fit function was replaced by coxsurv.fit, which has a few
#   more arguments and a more logical naming (it fits in with other calls).
# Enough people use rms, however, which calls survfitcoxph.fit, that we are
#  giving this as a temporary pass through.  

survfitcoxph.fit <- function(y, x, wt, x2, risk, newrisk, strata, se.fit,
                              survtype, vartype, varmat, id, y2, strata2,
                              unlist=TRUE) {
    .Depricated("coxsurv.fit", "survival")

    Call <- match.call()
    
    if (missing(survtype)) {
        stype <- 1
        ctype <- 1
    } else {
        stype <- c(1,2,2)[survtype]
        ctype <- c(1,1,2)[survtype]
    }

    indx <- match(c("y", "x", "wt", "x2", "y2", "risk", "strata", "strata2",
                    "se.fit", "varmat"), names(Call), nomatch=0)
    temp <-Call[c(1, "indx")]
    temp[[1]] <- as.name(survival:::coxsurv.fit)

    temp$ctype <- ctype
    temp$stype <- stype
    if (!missing(newrisk)) temp$risk2 <- newrisk
    if (!missing(id)) temp$id2 <- id
    temp$unlist <- unlist

    eval(temp, parent.frame())
}
