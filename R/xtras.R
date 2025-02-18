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

vcov.coxphms <- function (object, complete=TRUE, matrix=FALSE, ...) {
    # conform to the standard vcov results
    vmat <- object$var
    vname <- names(object$coefficients)
    dimnames(vmat) <- list(vname, vname)
    cmap <- object$cmap
    if (!complete && any(is.na(coef(object)))) {
        keep <- !is.na(coef(object))
        vmat <- vmat[keep, keep, drop=FALSE]
        cmap <- cmap[keep,]
    }
    
    if (!matrix) vmat
    else {
        v2 <- array(0, dim=c(dim(vout), ncol(cmap)),
                    dimnames= c(dimnames(vout), transition=list(colnames(cmap))))
        for (i in 1:ncol(cmap)) {
            j <- cmap[,i]
            v2[j>0, j>0, i] <- vmat[j,j]
        }
        v2
    }
}

vcov.survreg<-function (object, complete=TRUE, ...) {
    if (!complete && any(is.na(coef(object)))) {
        keep <- !is.na(coef(object))
        vv <- object$var[keep, keep, drop=FALSE]
        vname <- names(coef(object))[keep]
    }
    else {
        vv <- object$var
        vname <- names(coef(object))   # add dimnames
    }
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

labels.survreg <- function(object, ...) attr(object$terms, "term.labels")
labels.coxph <- function(object, ...) attr(object$terms, "term.labels")
labels.aareg <- function(object, ...) attr(object$terms, "term.labels")


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
# To support users who like the 'extended Kaplan-Meier' , i.e. a person can
#  switch curves midstream, we need to define someone as censored from the
#  first curve and an entry to the second whenever such a switch occurs.
# This routine returns 1*(first in a sequence) + 2*(last in a sequence),
#  which for the above is 1,0,2,3.  This assumes no overlapping intervals.

survflag <- function(y, id, group) {
    if (!inherits(y, "Surv")) stop("y must be a Surv object")
    if (nrow(y) != length(id)) stop("length mismatch")
    if (ncol(y) != 3) stop("y needs to be of (tstart, tstop) form")
  
    n <- nrow(y)
    if (missing(group))
        indx <- order(id, y[,2])  # sort the data by time within id
    else indx <- order(group, id, y[,2])
    y2 <- y[indx,]
    id2 <- id[indx]

    if (missing(group)) newid <- (id2[-n] != id2[-1])
    else {
        group2 <-as.numeric(group)[indx]  #normally group is a factor
        newid <- ((id2[-n] != id2[-1]) | (group2[-n] != group2[-1]))
    }       
    gap <-  (y2[-n,2] < y2[-1,1]) 

    flag <- unname(1L*c(TRUE, newid | gap) + 2L*c(newid | gap, TRUE))
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

# This is self defense for my functions agains the survival:: aficiandos.
# Replace survival::strata with strata, survival:cluster with cluster, etc.
# We don't need to replace survival::Surv, it is legal either way, and there
#  is a lot of usage of this.
#  Then update the envionment of the formula
removeDoubleColonSurv <- function(formula)
{
    sname <- c("Surv", "strata", "cluster", "pspline", "tt",
               "frailty", "ridge", "frailty", "frailty.gaussian",
               "frailty.gamma", "frailty.t")
    cname <- paste0("survival::", sname[-1])
    # three counts: survival::sname(), sname(), sname as variable
    found1 <- found2 <- found3 <- NULL 
    fix <- function(expr) {
        if (is.call(expr)) {
            if (!is.na(i <- match(deparse1(expr[[1]]), sname)))
                    found2 <<- c(found2, sname[i])
            else if (!is.na(i <- match(deparse1(expr[[1]]), cname))) {
                found1 <<- c(found1, sname[i+1])
                # remove the survival:: part
                expr[[1]] <- str2lang(paste0(sname[i+1], '()'))[[1]]
            }
            for(i in seq_along(expr)[-1]) {
                # the !is.null is to deal with the rare case of 
                # Surv(time, status) ~ NULL, which is legal (I was surprised:
                # and My.stepwise does it). Setting a term to NULL removes it,
                # leaving an invalid formula with no right hand side.
                if (!is.null(expr[[i]])) expr[[i]] <- fix(expr[[i]])
             }
        } else if (is.name(expr) && 
                   !is.na(i <- match(as.character(expr), sname))) 
            found3 <<- c(found3, sname[i])
        expr
    }
    newform <- fix(formula)

    # If something is used as a name and also a function, we can't add the
    #  function to our env, e.g.,  coxph(Surv(time,stat) ~ strata + strata(group)
    #  which is exactly what EPI:eff.match does. We will find the function as
    #  a match for both of them, since our env is searched first.
    # If some user has their own strata function and calls EPI:eff.match, they
    #  are SOL, we can't save them
    found <- unique(c(found1, found2))
    if (length(found3)) found <- found[!(found %in% found2)]

    # Return a list of 2 parts: the new formula (with updated environment),
    # and whether the "call" component should be rewritten. There are three
    # opinions wrt the second. a. If a survival:: was removed, then we want
    # that to also disappear from model printouts, a reinforcement for the
    # users that they shouldn't type that. (length of found1 > 0)
    # b. Option a is dishonest, the call should be what you typed  
    # c. Option a causes more downstream techncial troubles than it is worth.
    # We currenty opt for c.
    
    if (length(found) >0) { # most often true 
        list(formula = addSurvFun(newform, found), newcall=FALSE)
       #list(newform = addSurvFun(newform, found), newcall= !is.null(found1))
    } else NULL # don't return a new formula
}

# The second part of my defense. Because model.frame is not a part of the
#  survival package, invocations of Surv, strata, etc within a formula are
#  not guarranteed to use my version; if a user had their own local copy of Surv
#  it would use that!  To ensure I get the proper ones, we insert a new 
#  environment into the call chain, and attach it to the formula.
# Because this routine is a part of the surival package I can refer to
#  strata and etc below without resorting to the survival:: form.
#  (In fact, we found out that the :: form can fail, i.e., if another
#  package has Imports:survival in the DESCRIPTION file but does not
#  have import(survival) in the NAMESPACE.)
#  Thus get(i) rather than get(paste0("survival::", i))
#
addSurvFun <- function(formula, found) {
    myenv <- new.env(parent= environment(formula))
    tt <- function(x) x
    for (i in found) 
        assign(i, get(i), envir= myenv)
    environment(formula) <- myenv
    formula
}

# This is useful for a timeline data set; for a counting process one the
#  tdc function in tmerge already does this.
# Replace any NA with the most recent non-NA, for each id separately
#  Better known as "last value carried forward"
#
lvcf <- function(x, id, time) {
    if (!missing(time)) indx <- order(id, time)
    else indx <- order(id)   

    for (i in seq(along=x)) {
        j <- indx[i]
        if (!is.na(x[j]) || i==1 || id[j]!= id[jlag]) current <- x[j]
        else x[j] <- current
        jlag <- j
    }
    x
}
