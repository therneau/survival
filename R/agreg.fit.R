# Automatically generated from the noweb directory
agreg.fit <- function(x, y, strata, offset, init, control,
                        weights, method, rownames, resid=TRUE, nocenter=NULL)
    {
    nvar <- ncol(x)
    event <- y[,3]
    if (all(event==0)) stop("Can't fit a Cox model with 0 failures")

    if (missing(offset) || is.null(offset)) offset <- rep(0.0, nrow(y))
    if (missing(weights)|| is.null(weights))weights<- rep(1.0, nrow(y))
    else if (any(weights<=0)) stop("Invalid weights, must be >0")
    else weights <- as.vector(weights)

    # Find rows to be ignored.  We have to match within strata: a
    #  value that spans a death in another stratum, but not it its
    #  own, should be removed.  Hence the per stratum delta
    if (length(strata) ==0) {y1 <- y[,1]; y2 <- y[,2]}
    else  {
        if (is.numeric(strata)) strata <- as.integer(strata)
        else strata <- as.integer(as.factor(strata))
        delta  <-  strata* (1+ max(y[,2]) - min(y[,1]))
        y1 <- y[,1] + delta
        y2 <- y[,2] + delta
    }
    event <- y[,3] > 0
    dtime <- sort(unique(y2[event]))
    indx1 <- findInterval(y1, dtime)
    indx2 <- findInterval(y2, dtime) 
    # indx1 != indx2 for any obs that spans an event time
    ignore <- (indx1 == indx2)
    nused  <- sum(!ignore)

    # Sort the data (or rather, get a list of sorted indices)
    #  For both stop and start times, the indices go from last to first
    if (length(strata)==0) {
        sort.end  <- order(ignore, -y[,2]) -1L #indices start at 0 for C code
        sort.start<- order(ignore, -y[,1]) -1L
        strata <- rep(0L, nrow(y))
        }
    else {
        sort.end  <- order(ignore, strata, -y[,2]) -1L
        sort.start<- order(ignore, strata, -y[,1]) -1L
        }

    if (is.null(nvar) || nvar==0) {
        # A special case: Null model.  Just return obvious stuff
        #  To keep the C code to a small set, we call the usual routines, but
        #  with a dummy X matrix and 0 iterations
        nvar <- 1
        x <- matrix(as.double(1:nrow(y)), ncol=1)  #keep the .C call happy
        maxiter <- 0
        nullmodel <- TRUE
        if (length(init) !=0) stop("Wrong length for inital values")
        init <- 0.0  #dummy value to keep a .C call happy (doesn't like 0 length)
        }
    else {
        nullmodel <- FALSE
        maxiter <- control$iter.max
        
        if (is.null(init)) init <- rep(0., nvar)
        if (length(init) != nvar) stop("Wrong length for inital values")
        }

    # 2021 change: pass in per covariate centering.  This gives
    #  us more freedom to experiment.  Default is to leave 0/1 variables alone
    if (is.null(nocenter)) zero.one <- rep(FALSE, ncol(x))
    zero.one <- apply(x, 2, function(z) all(z %in% nocenter)) 

    # the returned value of agfit$coef starts as a copy of init, so make sure
    #  is is a vector and not a matrix; as.double suffices.
    # Solidify the storage mode of other arguments
    storage.mode(y) <- storage.mode(x) <- "double"
    storage.mode(offset) <- storage.mode(weights) <- "double"
    agfit <- .Call(Cagfit4, nused, 
                   y, x, strata, weights, 
                   offset,
                   as.double(init), 
                   sort.start, sort.end, 
                   as.integer(method=="efron"),
                   as.integer(maxiter), 
                   as.double(control$eps),
                   as.double(control$toler.chol),
                   ifelse(zero.one, 0L, 1L))
    # agfit4 centers variables within strata, so does not return a vector
    #  of means.  Use a fill in consistent with other coxph routines
    agmeans <- ifelse(zero.one, 0, colMeans(x))

    vmat <- agfit$imat
    coef <- agfit$coef
    if (agfit$flag[1] < nvar) which.sing <- diag(vmat)==0
    else which.sing <- rep(FALSE,nvar)

    if (maxiter >1) {
        infs <- abs(agfit$u %*% vmat)
        if (any(!is.finite(coef)) || any(!is.finite(vmat)))
            stop("routine failed due to numeric overflow. This should never happen. Please contact the author.")
        if (agfit$flag[4] > 0)
            warning("Ran out of iterations and did not converge")
        else {
            infs <- (!is.finite(agfit$u) |
                     infs > control$toler.inf*(1+ abs(coef)))
            if (any(infs))
                warning(gettextf("Loglik converged before variable %s; beta may be infinite.",
                              paste((1:nvar)[infs],collapse=",")))
        }
    }
    lp  <- as.vector(x %*% coef + offset - sum(coef * agmeans))
    if (resid) {
        if (any(lp > log(.Machine$double.xmax))) {
            # prevent a failure message due to overflow
            #  this occurs with near-infinite coefficients
            temp <- lp + log(.Machine$double.xmax) - (1 + max(lp))
            score <- exp(temp)
        } else score <- exp(lp)

        residuals <- .Call(Cagmart3, nused,
                       y, score, weights,
                       strata,
                       sort.start, sort.end,
                       as.integer(method=='efron'))
        names(residuals) <- rownames
    }

    # The if-then-else below is a real pain in the butt, but the tccox
    #  package's test suite assumes that the ORDER of elements in a coxph
    #  object will never change.
    #
    if (nullmodel) {
        rval <- list(loglik=agfit$loglik[2],
             linear.predictors = offset,
             method= method,
             class = c("coxph.null", 'coxph') )
        if (resid) rval$residuals <- residuals
    }
    else {
        names(coef) <- dimnames(x)[[2]]
        if (maxiter > 0) coef[which.sing] <- NA  # always leave iter=0 alone
        flag <- agfit$flag
        names(flag) <- c("rank", "rescale", "step halving", "convergence")
        
        if (resid) {
            rval <- list(coefficients  = coef,
                         var    = vmat,
                         loglik = agfit$loglik,
                         score  = agfit$sctest,
                         iter   = agfit$iter,
                         linear.predictors = as.vector(lp),
                         residuals = residuals, 
                         means = agmeans,
                         first = agfit$u,
                         info = flag,
                         method= method,
                         class = "coxph")
        } else {
             rval <- list(coefficients  = coef,
                         var    = vmat,
                         loglik = agfit$loglik,
                         score  = agfit$sctest,
                         iter   = agfit$iter,
                         linear.predictors = as.vector(lp),
                         means = agmeans,
                         first = agfit$u,
                         info = flag,
                         method = method,
                         class = "coxph")
        }
        rval
    }
    rval        
}  
