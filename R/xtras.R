vcov.coxph <- function (object, complete=TRUE, ...) {
    # conform to the standard vcov results
    vmat <- object$var
    vname <- names(object$coefficients)
    dimnames(vmat) <- list(vname, vname)
    if (!complete && any(is.na(coef(object)))) {
        keep <- !is.na(coef(object))
        vmat[keep, keep, drop=FALSE]
        }
    else vmat
}

vcov.survreg<-function (object, complete=TRUE, ...) {
    if (!complete && any(is.na(coef(object)))) {
        keep <- !is.na(coef(object))
        vv <- object$var[keep, keep, drop=FALSE]
    }
    else vv <- object$var
    vname <- names(coef(object))   # add dimnames
    extra <- ncol(vv) - length(vname)
    if (extra ==1) vname <- c(vname, "Log(scale)")
    else if(extra >1) 
        vname <- c(vname, paste("Log(scale[", names(object$scale), "])", sep=''))
    dimnames(vv) <- list(vname, vname)
    vv
}

# The extractAIC methods for coxph and survreg objects are defined
#  in the stats package.  Don't reprise them here.
extractAIC.coxph.penal<- function(fit,scale,k=2,...){
    edf<-sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

extractAIC.coxph.null <- function(fit, scale, k=2, ...) {
    c(0, -2*fit$loglik[1])
}

labels.survreg <- function(object, ...) attr(object,"term.labels")

rep.Surv <- function(x, ...) {
    indx <- rep(1:nrow(x), ...)
    x[indx,]
}

# This function is just like all.vars -- except that it does not recur
#  on the $ sign, it follows both arguments of +, *, - and : in order to
#  track formulas, all arguments of Surv, and only the first of things 
#  like ns().  And - it works only on formulas.
# This is used to generate a warning in coxph if the same variable is used
#  on both sides, so perfection is not required of the function.
# Changed from terms.inner to innerterms at CRAN's request, as the former
#  created a false positive as an undocumented method for the terms generic
innerterms <- function(x) {
    if (inherits(x, "formula")) {
        if (length(x) ==3) c(innerterms(x[[2]]), innerterms(x[[3]]))
        else innerterms(x[[2]])
    }
    else if (inherits(x, "call") && 
             (x[[1]] != as.name("$") && x[[1]] != as.name("["))) {
        if (x[[1]] == '+' || x[[1]]== '*' || x[[1]] == '-' || x[[1]] ==':') {
            # terms in a model equation, unary minus only has one argument
            if (length(x)==3) c(innerterms(x[[2]]), innerterms(x[[3]]))
            else innerterms(x[[2]])
        }
        else if (x[[1]] == as.name("Surv"))
                 unlist(lapply(x[-1], innerterms))
        else if (length(x) ==2) innerterms(x[[2]])
        else character(0)
    }
    else(deparse(x))
}
   
# If a subject had (start, stop) observations of (1,2) (2,10) (10,15) (20,25),
#  say, code often wants to distiguish intervals that are "real" censoring
#  from a simple split due to a time dependent covariate.
# This routine returns 1*(first in a sequence) + 2*(last in a sequence),
#  which for the above is 1,0,2,3.  This assumes no overlapping intervals

survflag <- function(y, id) {
    if (!inherits(y, "Surv")) stop("y must be a Surv object")
    if (nrow(y) != length(id)) stop("length mismatch")
    if (ncol(y) != 3) stop("y needs to be of (tstart, tstop) form")
  
    n <- nrow(y)
    indx <- order(id, y[,2])  # sort the data by time within id
    y2 <- y[indx,]
    id2 <- id[indx]

    newid <- (id2[-n] != id2[-1])
    gap <-  (y2[-n,2] < y2[-1,1]) 
   
    flag <- 1L*c(TRUE, newid | gap) + 2L*c(newid | gap, TRUE)
    flag[indx] <- flag   # return it to data order
    flag
}


# Dummy methods, to create an informative error message
coef.survfit <- function(object, ...) 
    stop("coef method not applicable for survfit objects")
vcov.survfit <- function(object, ...) 
    stop("vcov method not applicable for survfit objects")
confint.survfit <- function(object, ...)
    stop(paste("confint method not defined for survfit objects," ,
         "use quantile for confidence intervals of the median survival"))

