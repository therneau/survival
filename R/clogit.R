## conditional logistic regression
##
## case ~ exposure + strata(matching)
##
clogit<-function(formula, data, weights, subset, na.action,
                 method=c("exact","approximate", "efron", "breslow"),
                 ... ) {
    
    Call <- match.call()  # how we were called

    # Create a call to model.frame() that contains the formula (required)
    #  and the data argument (if present).
    # It's only job is to find out the number of rows in the data
    #  before subset or na.action are applied.
    indx <- match(c("formula", "data"), names(Call), nomatch=0)
    if (indx[1]==0) stop("A formula argument is required")
    mf <- Call[c(1,indx)]
    mf[[1]] <- as.name("model.frame")
    mf$na.action <- "na.pass"
    nrows<-NROW(eval(mf, parent.frame()))
    method <- match.arg(method)

    # Catch the rare case of a person asking for robust variance, and give
    #  them a nicer warning than will occur if they fall through to the
    #  coxph call
    if (missing(data)) temp <- terms(formula, special='cluster')
    else temp <- terms(formula, special="cluster", data=data)
    if (!is.null(attr(temp, 'specials')$cluster) && method=="exact")
        stop("robust variance plus the exact method is not supported")

    # Now build a call to coxph with the formula fixed up to have
    #  our special left hand side.
    coxcall <- Call
    coxcall[[1]] <- as.name("coxph")
    newformula <- formula
    newformula[[2]] <- substitute(Surv(rep(1,nn),case),
                                list(case=formula[[2]],nn=nrows))
    environment(newformula) <- environment(formula)
    coxcall$formula<-newformula

    coxcall$method <- switch(method, exact="exact",
                             efron="efron",
                             "breslow")
    if (!is.null(coxcall$weights)) {
        coxcall$weights <- NULL
        warning("Weights are ignored in clogit")
    }
    coxcall<-eval(coxcall, sys.frame(sys.parent()))
    coxcall$userCall<-sys.call()
    
    class(coxcall)<-c("clogit","coxph")
    coxcall
}


print.clogit <- function(x,...){

    x$call<-x$userCall
    NextMethod()

}

survfit.clogit <- function(formula, ...)
    stop("predicted survival curves are not defined for a clogit model")

