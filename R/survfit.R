# Automatically generated from the noweb directory
survfit <- function(formula, ...) {
    UseMethod("survfit")
}

survfit.formula <- function(formula, data, weights, subset, 
                            na.action, stype=1, ctype=1, 
                            id, cluster, robust, istate, 
                            timefix=TRUE, etype, error, ...) {

    Call <- match.call()
    Call[[1]] <- as.name('survfit')  #make nicer printout for the user
    # create a copy of the call that has only the arguments we want,
    #  and use it to call model.frame()
    indx <- match(c('formula', 'data', 'weights', 'subset','na.action',
                    'istate', 'id', 'cluster', "etype"), names(Call), nomatch=0)
    #It's very hard to get the next error message other than malice
    #  eg survfit(wt=Surv(time, status) ~1) 
    if (indx[1]==0) stop("a formula argument is required")
    temp <- Call[c(1, indx)]
    temp[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(temp)

    Terms <- terms(formula, c("strata", "cluster"))
    ord <- attr(Terms, 'order')
    if (length(ord) & any(ord !=1))
            stop("Interaction terms are not valid for this function")

    n <- nrow(mf)
    Y <- model.response(mf)
    if (inherits(Y, "Surv2")) {
        # this is Surv2 style data
        # if there are any obs removed due to missing, remake the model frame
        if (length(attr(mf, "na.action"))) {
            temp$na.action <- na.pass
            mf <- eval.parent(temp)
        }
        if (!is.null(attr(Terms, "specials")$cluster))
            stop("cluster() cannot appear in the model statement")
        new <- surv2data(mf)
        mf <- new$mf
        istate <- new$istate
        id <- new$id
        Y <- new$y
        if (anyNA(mf[-1])) { #ignore the response variable still found there
            if (missing(na.action)) temp <- get(getOption("na.action"))(mf[-1])
            else temp <- na.action(mf[-1])
            omit <- attr(temp, "na.action")
            mf <- mf[-omit,]
            Y <- Y[-omit]
            id <- id[-omit]
            istate <- istate[-omit]
        }                      
        n <- nrow(mf)
    }       
    else {
        if (!is.Surv(Y)) stop("Response must be a survival object")
        id <- model.extract(mf, "id")
        istate <- model.extract(mf, "istate")
    }
    if (n==0) stop("data set has no non-missing observations")

    casewt <- model.extract(mf, "weights")
    if (is.null(casewt)) casewt <- rep(1.0, n)
    else {
        if (!is.numeric(casewt)) stop("weights must be numeric")
        if (any(!is.finite(casewt))) stop("weights must be finite") 
        if (any(casewt <0)) stop("weights must be non-negative")
        casewt <- as.numeric(casewt)  # transform integer to numeric
    }

    if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

    cluster <- model.extract(mf, "cluster")
    temp <- untangle.specials(Terms, "cluster")
    if (length(temp$vars)>0) {
        if (length(cluster) >0) stop("cluster appears as both an argument and a model term")
        if (length(temp$vars) > 1) stop("can not have two cluster terms")
        cluster <- mf[[temp$vars]]
        Terms <- Terms[-temp$terms]
    }

    ll <- attr(Terms, 'term.labels')
    if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
    else X <- strata(mf[ll])

    # Backwards support for the now-depreciated etype argument
    etype <- model.extract(mf, "etype")
    if (!is.null(etype)) {
        if (attr(Y, "type") == "mcounting" ||
            attr(Y, "type") == "mright")
            stop("cannot use both the etype argument and mstate survival type")
        if (length(istate)) 
            stop("cannot use both the etype and istate arguments")
        status <- Y[,ncol(Y)]
        etype <- as.factor(etype)
        temp <- table(etype, status==0)

        if (all(rowSums(temp==0) ==1)) {
            # The user had a unique level of etype for the censors
            newlev <- levels(etype)[order(-temp[,2])] #censors first
        }
        else newlev <- c(" ", levels(etype)[temp[,1] >0])
        status <- factor(ifelse(status==0,0, as.numeric(etype)),
                             labels=newlev)

        if (attr(Y, 'type') == "right")
            Y <- Surv(Y[,1], status, type="mstate")
        else if (attr(Y, "type") == "counting")
            Y <- Surv(Y[,1], Y[,2], status, type="mstate")
        else stop("etype argument incompatable with survival type")
    }
                         
    # Deal with the near-ties problem
    if (!is.logical(timefix) || length(timefix) > 1)
        stop("invalid value for timefix option")
    if (timefix) newY <- aeqSurv(Y) else newY <- Y
    
    if (missing(robust)) robust <- NULL
    # Call the appropriate helper function
    if (attr(Y, 'type') == 'left' || attr(Y, 'type') == 'interval')
        temp <-  survfitTurnbull(X, newY, casewt, cluster= cluster,
                                 robust= robust, ...)
    else if (attr(Y, 'type') == "right" || attr(Y, 'type')== "counting")
        temp <- survfitKM(X, newY, casewt, stype=stype, ctype=ctype, id=id, 
                          cluster=cluster, robust=robust, ...)
    else if (attr(Y, 'type') == "mright" || attr(Y, "type")== "mcounting")
        temp <- survfitCI(X, newY, weights=casewt, stype=stype, ctype=ctype, 
                          id=id, cluster=cluster, robust=robust, 
                          istate=istate, ...)
    else {
        # This should never happen
        stop("unrecognized survival type")
    }

    # If a stratum had no one beyond start.time, the length 0 gives downstream
    #  failure, e.g., there is no sensible printout for summary(fit, time= 100)
    #  for such a curve
    temp$strata <- temp$strata[temp$strata >0]  
    if (is.null(temp$states)) class(temp) <- 'survfit'
    else class(temp) <- c("survfitms", "survfit")

    if (!is.null(attr(mf, 'na.action')))
            temp$na.action <- attr(mf, 'na.action')

    temp$call <- Call
    temp
    }
dim.survfit <- function(x) {
    d1name <- "strata"
    d2name <- "data"
    d3name <- "states"
    if (is.null(x$strata))  {d1 <- d1name <- NULL} else d1 <- length(x$strata)
    if (is.null(x$newdata)) {d2 <- d2name <- NULL} else d2 <- nrow(x$newdata)
    if (is.null(x$states))  {d3 <- d3name <- NULL} else d3 <- length(x$states)
    
    if (inherits(x, "survfitcox") && is.null(d2) && is.null(d3) &&
        is.matrix(x$surv)) {
        # older style survfit.coxph object, before I added newdata to the output
        d2name <- "data"
        d2 <- ncol(x$surv)
    }

    dd <- c(d1, d2, d3)
    names(dd) <- c(d1name, d2name, d3name)
    dd
}

# there is a separate function for survfitms objects
"[.survfit" <- function(x, ... , drop=TRUE) {
    nmatch <- function(indx, target) { 
        # This function lets R worry about character, negative, or 
        #  logical subscripts.
        #  It always returns a set of positive integer indices
        temp <- 1:length(target)
        names(temp) <- target
        temp[indx]
    }
    
    if (!inherits(x, "survfit")) stop("[.survfit called on non-survfit object")
    ndots <- ...length()      # the simplest, but not avail in R 3.4
    # ndots <- length(list(...))# fails if any are missing, e.g. fit[,2]
    # ndots <- if (missing(drop)) nargs()-1 else nargs()-2  # a workaround

    dd <- dim(x)
    # for dd=NULL, an object with only one curve, x[1] is always legal
    if (is.null(dd)) dd <- c(strata=1L) # survfit object with only one curve
    dtype <- match(names(dd), c("strata", "data", "states"))

    if (ndots >0 && !missing(..1)) i <- ..1 else i <- NULL
    if (ndots> 1 && !missing(..2)) j <- ..2 else j <- NULL
    
    if (ndots > length(dd)) 
        stop("incorrect number of dimensions")
    if (length(dtype) > 2) stop("invalid survfit object")  # should never happen
    if (is.null(i) && is.null(j)) {
        # called with no subscripts given -- return x untouched
        return(x)
    }
    
    # Code below is easier if "i" is always the strata
    if (dtype[1] !=1) {
        dtype <- c(1, dtype)
        j <- i; i <- NULL
        dd <- c(1, dd)
        ndots <- ndots +1
    }       

   # We need to make a new one
    newx <- vector("list", length(x))
    names(newx) <- names(x)
    for (k in c("logse", "version", "conf.int", "conf.type", "type", "call"))
        if (!is.null(x[[k]])) newx[[k]] <- x[[k]]
    class(newx) <- class(x)
    
    if (ndots== 1 && length(dd)==2) {
        # one subscript given for a two dimensional object
        # If one of the dimensions is 1, it is easier for me to fill in i and j
        if (dd[1]==1) {j <- i; i<- 1}
        else if (dd[2]==1) j <- 1
        else {
            #  the user has a mix of rows/cols
            index <- 1:prod(dd)
            itemp <- matrix(index, nrow=dd[1])
            keep <- itemp[i]   # illegal subscripts will generate an error
            if (length(keep) == length(index) && all(keep==index)) return(x)

            ii <- row(itemp)[keep]
            jj <- col(itemp)[keep]
            # at this point we have a matrix subscript of (ii, jj)
            # expand into a long pair of rows and cols
            temp <- split(seq(along.with=x$time), 
                          rep(1:length(x$strata), x$strata))
            indx1 <- unlist(temp[ii])   # rows of the surv object
            indx2 <- rep(jj, x$strata[ii])
        
            # return with each curve as a separate strata
            newx$n <- x$n[ii]
            for (k in c("time", "n.risk", "n.event", "n.censor", "n.enter"))
                if (!is.null(x[[k]])) newx[[k]] <- (x[[k]])[indx1]
            k <- cbind(indx1, indx2)
            for (j in c("surv", "std.err", "upper", "lower", "cumhaz",
                        "std.chaz", "influence.surv", "influence.chaz"))
                if (!is.null(x[[j]])) newx[[j]] <- (x[[j]])[k]
            temp <- x$strata[ii]
            names(temp) <- 1:length(ii)
            newx$strata <- temp
            return(newx)
        }
    }
    
    # irow will be the rows that need to be taken
    #  j the columns (of present)
    if (is.null(x$strata)) {
           if (is.null(i) || all(i==1)) irow <- seq(along.with=x$time)
           else stop("subscript out of bounds")
           newx$n <- x$n
    }
    else { 
        if (is.null(i)) indx <- seq(along.with= x$strata)
        else indx <- nmatch(i, names(x$strata)) #strata to keep
        if (any(is.na(indx))) 
            stop(paste("strata", 
                       paste(i[is.na(indx)], collapse=' '),
                       'not matched'))
        # Now, indx may not be in order: some can use curve[3:2] to reorder
        #  The list/unlist construct will reorder the data
        temp <- split(seq(along.with =x$time), 
                      rep(1:length(x$strata), x$strata))
        irow <- unlist(temp[indx])
        
        if (length(indx) <=1 && drop) newx$strata <- NULL
        else               newx$strata  <- x$strata[i]

        newx$n <- x$n[indx]
        if (length(indx) ==1 & drop) x$strata <- NULL
        else    newx$strata <- x$strata[indx]
    }

    if (length(dd)==1) {  # no j dimension
        for (k in c("time", "n.risk", "n.event", "n.censor", "n.enter",
               "surv", "std.err", "cumhaz", "std.chaz", "upper", "lower",
               "influence.surv", "influence.chaz"))
            if (!is.null(x[[k]])) newx[[k]] <- (x[[k]])[irow]
    }
       
    else { # 2 dimensional object
        if (is.null(j)) j <- seq.int(ncol(x$surv))
        # If the curve has been selected by strata and keep has only
        #  one row, we don't want to lose the second subscript too
        if (length(irow)==1)  drop <- FALSE

        for (k in c("time", "n.risk", "n.event", "n.censor", "n.enter"))
                 if (!is.null(x[[k]])) newx[[k]] <- (x[[k]])[irow]
        for (k in c("surv", "std.err", "cumhaz", "std.chaz", "upper", "lower",
               "influence.surv", "influence.chaz"))
            if (!is.null(x[[k]])) newx[[k]] <- (x[[k]])[irow, j, drop=drop]
    }
    newx
}
survfit.Surv <- function(formula, ...)
    stop("the survfit function requires a formula as its first argument")
survfit_confint <- function(p, se, logse=TRUE, conf.type, conf.int,
                            selow, ulimit=TRUE) {
    zval <- qnorm(1- (1-conf.int)/2, 0,1)
    if (missing(selow)) scale <- 1.0
    else scale <- ifelse(selow==0, 1.0, selow/se)  # avoid 0/0 at the origin
    if (!logse) se <- ifelse(se==0, 0, se/p)   # se of log(survival) = log(p)

    if (conf.type=='plain') {
        se2 <- se* p * zval  # matches equation 4.3.1 in Klein & Moeschberger
        if (ulimit) list(lower= pmax(p -se2*scale, 0), upper = pmin(p + se2, 1))
        else  list(lower= pmax(p -se2*scale, 0), upper = p + se2)
    }
    else if (conf.type=='log') {
        #avoid some "log(0)" messages
        xx <- ifelse(p==0, NA, p)  
        se2 <- zval* se 
        temp1 <- exp(log(xx) - se2*scale)
        temp2 <- exp(log(xx) + se2)
        if (ulimit) list(lower= temp1, upper= pmin(temp2, 1))
        else  list(lower= temp1, upper= temp2)
    }
    else if (conf.type=='log-log') {
        xx <- ifelse(p==0 | p==1, NA, p)
        se2 <- zval * se/log(xx)
        temp1 <- exp(-exp(log(-log(xx)) - se2*scale))
        temp2 <- exp(-exp(log(-log(xx)) + se2))
        list(lower = temp1 , upper = temp2)
    }
    else if (conf.type=='logit') {
        xx <- ifelse(p==0, NA, p)  # avoid log(0) messages
        se2 <- zval * se *(1 + xx/(1-xx))
 
        temp1 <- 1- 1/(1+exp(log(p/(1-p)) - se2*scale))
        temp2 <- 1- 1/(1+exp(log(p/(1-p)) + se2))
        list(lower = temp1, upper=temp2)
    }
    else if (conf.type=="arcsin") {
        xx <- ifelse(p==0, NA, p)
        se2 <- .5 *zval*se * sqrt(xx/(1-xx))
        list(lower= (sin(pmax(0, asin(sqrt(xx)) - se2*scale)))^2,
             upper= (sin(pmin(pi/2, asin(sqrt(xx)) + se2)))^2)
    }
    else stop("invalid conf.int type")
}
