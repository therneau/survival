# Automatically generated from the noweb directory
survexp <- function(formula, data,
        weights, subset, na.action, rmap, times,
        method=c("ederer", "hakulinen", "conditional", "individual.h", 
                 "individual.s"),
        cohort=TRUE,  conditional=FALSE,
        ratetable=survival::survexp.us, scale=1, se.fit,
        model=FALSE, x=FALSE, y=FALSE) {
    Call <- match.call()
        
    # keep the first element (the call), and the following selected arguments
    indx <- match(c('formula', 'data', 'weights', 'subset', 'na.action'),
                      names(Call), nomatch=0)
    if (indx[1] ==0) stop("A formula argument is required")
    tform <- Call[c(1,indx)]  # only keep the arguments we wanted
    tform[[1L]] <- quote(stats::model.frame)  # change the function called
        
    Terms <- if(missing(data)) terms(formula)
             else              terms(formula, data=data)
    if (!missing(rmap)) {
        rcall <- substitute(rmap)
        if (!is.call(rcall) || rcall[[1]] != as.name('list'))
            stop ("Invalid rcall argument")
        }
    else rcall <- NULL   # A ratetable, but no rcall argument

    # Check that there are no illegal names in rcall, then expand it
    #  to include all the names in the ratetable
    if (is.ratetable(ratetable))   {
        varlist <- names(dimnames(ratetable))
        if (is.null(varlist)) varlist <- attr(ratetable, "dimid") # older style
    }
    else if(inherits(ratetable, "coxph") && !inherits(ratetable, "coxphms")) {
        ## Remove "log" and such things, to get just the list of
        #   variable names
        varlist <- all.vars(delete.response(ratetable$terms))
        }
    else stop("Invalid rate table")

    temp <- match(names(rcall)[-1], varlist) # 2,3,... are the argument names
    if (any(is.na(temp)))
        stop("Variable not found in the ratetable:", (names(rcall))[is.na(temp)])
        
    if (any(!(varlist %in% names(rcall)))) {
        to.add <- varlist[!(varlist %in% names(rcall))]
        temp1 <- paste(text=paste(to.add, to.add, sep='='), collapse=',')
        if (is.null(rcall)) rcall <- parse(text=paste("list(", temp1, ")"))[[1]]
        else {
            temp2 <- deparse(rcall)
            rcall <- parse(text=paste("c(", temp2, ",list(", temp1, "))"))[[1]]
            }
        }
    # Create a temporary formula, used only in the call to model.frame
    newvar <- all.vars(rcall)
    if (length(newvar) > 0) {
        temp <- paste(paste(deparse(Terms), collapse=""),  
                       paste(newvar, collapse='+'), sep='+')
        tform$formula <- as.formula(temp, environment(Terms))
        }
    mf <- eval(tform, parent.frame())
    n <- nrow(mf)
    if (n==0) stop("Data set has 0 rows")
    if (!missing(se.fit) && se.fit)
        warning("se.fit value ignored")

    weights <- model.extract(mf, 'weights')
    if (length(weights) ==0) weights <- rep(1.0, n)
    if (inherits(ratetable, 'ratetable') && any(weights !=1))
        warning("weights ignored")

    if (any(attr(Terms, 'order') >1))
            stop("Survexp cannot have interaction terms")
    if (!missing(times)) {
        if (any(times<0)) stop("Invalid time point requested")
        if (length(times) >1 )
            if (any(diff(times)<0)) stop("Times must be in increasing order")
        }
    Y <- model.extract(mf, 'response')
    no.Y <- is.null(Y)
    if (no.Y) {
        if (missing(times)) {
            if (is.ratetable(ratetable)) 
                stop("either a times argument or a response is needed")
            }
        else newtime <- times
        }
    else {
        if (is.matrix(Y)) {
            if (is.Surv(Y) && attr(Y, 'type')=='right') Y <- Y[,1]
            else stop("Illegal response value")
            }
        if (any(Y<0)) stop ("Negative follow up time")
    #    if (missing(npoints)) temp <- unique(Y)
    #    else                  temp <- seq(min(Y), max(Y), length=npoints)
        temp <- unique(Y)
        if (missing(times)) newtime <- sort(temp)
        else  newtime <- sort(unique(c(times, temp[temp<max(times)])))
        }

    if (!missing(method)) method <- match.arg(method)
    else {
        # the historical defaults and older arguments
        if (!missing(conditional) && conditional) method <- "conditional"
        else {
            if (no.Y) method <- "ederer"
            else method <- "hakulinen"
            }
        if (!missing(cohort) && !cohort) method <- "individual.s"
        }
    if (no.Y && (method!="ederer")) 
        stop("a response is required in the formula unless method='ederer'")
    ovars <- attr(Terms, 'term.labels')
    # rdata contains the variables matching the ratetable
    rdata <- data.frame(eval(rcall, mf), stringsAsFactors=TRUE)  
    if (is.ratetable(ratetable)) {
        israte <- TRUE
        if (no.Y) {
            Y <- rep(max(times), n)
            }
        rtemp <- match.ratetable(rdata, ratetable)
        R <- rtemp$R
        }
    else if (inherits(ratetable, 'coxph')) {
        israte <- FALSE
        Terms <- ratetable$terms
        }
    else if (inherits(ratetable, "coxphms"))
        stop("survexp not defined for multi-state coxph models")
    else stop("Invalid ratetable")
    if (substring(method, 1, 10) == "individual") { #individual survival
        if (no.Y) stop("for individual survival an observation time must be given")
        if (israte)
             temp <- survexp.fit (1:n, R, Y, max(Y), TRUE, ratetable)
        else {
            rmatch <- match(names(data), names(rdata))
            if (any(is.na(rmatch))) rdata <- cbind(rdata, data[,is.na(rmatch)])
            temp <- survexp.cfit(1:n, rdata, Y, 'individual', ratetable)
        }
        if (method == "individual.s") xx <- temp$surv
        else xx <- -log(temp$surv)
        names(xx) <- row.names(mf)
        na.action <- attr(mf, "na.action")
        if (length(na.action)) return(naresid(na.action, xx))
        else return(xx)
        }
    if (length(ovars)==0)  X <- rep(1,n)  #no categories
    else {
        odim <- length(ovars)
        for (i in 1:odim) {
            temp <- mf[[ovars[i]]]
            ctemp <- class(temp)
            if (!is.null(ctemp) && ctemp=='tcut')
                stop("Can't use tcut variables in expected survival")
            }
        X <- strata(mf[ovars])
        }

    #do the work
    if (israte)
        temp <- survexp.fit(as.numeric(X), R, Y, newtime,
                           method=="conditional", ratetable)
    else {
        temp <- survexp.cfit(as.numeric(X), rdata, Y, method, ratetable, weights)
        newtime <- temp$time
        }
    if (missing(times)) {
        n.risk <- temp$n
        surv <- temp$surv
        }
    else {
        if (israte) keep <- match(times, newtime)
        else {
            # The result is from a Cox model, and it's list of
            #  times won't match the list requested in the user's call
            # Interpolate the step function, giving survival of 1
            #  for requested points that precede the Cox fit's
            #  first downward step.  The code is like summary.survfit.
            n <- length(temp$time)
            keep <- approx(temp$time, 1:n, xout=times, yleft=0,
                           method='constant', f=0, rule=2)$y
            }

        if (is.matrix(temp$surv)) {
            surv <- (rbind(1,temp$surv))[keep+1,,drop=FALSE]
            n.risk <- temp$n[pmax(1,keep),,drop=FALSE]
             }
        else {
            surv <- (c(1,temp$surv))[keep+1]
            n.risk <- temp$n[pmax(1,keep)]
            }
        newtime <- times
        }
    newtime <- newtime/scale
    if (is.matrix(surv)) {
        dimnames(surv) <- list(NULL, levels(X))
        out <- list(call=Call, surv= drop(surv), n.risk=drop(n.risk),
                        time=newtime)
        }
    else {
         out <- list(call=Call, surv=c(surv), n.risk=c(n.risk),
                       time=newtime)
         }
    if (model) out$model <- mf
    else {
        if (x) out$x <- X
        if (y) out$y <- Y
        }
    if (israte && !is.null(rtemp$summ)) out$summ <- rtemp$summ
    if (no.Y) out$method <- 'Ederer'
    else if (conditional) out$method <- 'conditional'
    else                  out$method <- 'cohort'
    class(out) <- c('survexp', 'survfit')
    out
}
