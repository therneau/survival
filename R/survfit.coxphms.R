#
# Curves from a coxphms object
#
survfit.coxphms <-
function(formula, newdata, se.fit=FALSE, conf.int=.95,
         individual= FALSE, stype=2, ctype, 
         conf.type=c("log", "log-log", "plain", "none", "logit", "arcsin"),
         censor=TRUE, start.time, id, influence=FALSE,
         na.action, type, p0=NULL, time0=FALSE, ...) {

    if (!inherits(formula, "coxphms"))
        stop("argument must be a coxphms object")
    Call <- match.call()
    Call[[1]] <- as.name("survfit")  #nicer output for the user
    object <- formula     #'formula' because it has to match survfit
    se.fit <- FALSE   #still to do
    if (missing(newdata))
        stop("multi-state survival requires a newdata argument")
    if (!missing(id) || individual)
        stop("using a covariate path is not supported for multi-state")
    missid <- TRUE
    individual <- FALSE
    if (se.fit) {
        warning("se.fit not yet implemented for multistate coxph models")
        se.fit <- FALSE
    }

    # process options, set up Y and the model frame for the original data
    if (!missing(type)) {  # old style argument
        if (!missing(stype) || !missing(ctype))
            warning("type argument ignored")
        else {
            temp1 <- c("kalbfleisch-prentice", "aalen", "efron",
                       "kaplan-meier", "breslow", "fleming-harrington",
                       "greenwood", "tsiatis", "exact")
            
            survtype <- match(match.arg(type, temp1), temp1)
            stype <- c(1,2,2,1,2,2,2,2,2)[survtype]
            if (stype!=1) ctype <-c(1,1,2,1,1,2,1,1,1)[survtype]
        }
    }
    if (missing(ctype)) {
        # Use the appropriate one from the model
        temp1 <- match(object$method, c("exact", "breslow", "efron"))
        ctype <- c(1,1,2)[temp1]
    }
    else if (!(ctype %in% 1:2)) stop ("ctype must be 1 or 2")
    if (!(stype %in% 1:2)) stop("stype must be 1 or 2")

    if (!se.fit) conf.type <- "none"
    else conf.type <- match.arg(conf.type)

    Terms  <- terms(object)
    tfac <- attr(Terms, 'factors')
    temp <- attr(Terms, 'specials')$strata 
    has.strata <- !is.null(temp)
    if (has.strata) {
        strata <- object$strata
        stangle = untangle.specials(Terms, "strata")  #used later
        # Toss out strata terms in tfac before doing the test 1 line below, as
        #  strata end up in the model with age:strat(grp) terms or *strata() 
        #  terms  (There might be more than one strata term)
        for (i in temp) tfac <- tfac[,tfac[i,] ==0]  # toss out strata terms
        if (any(tfac >1))
            stop("not able to create a curve for models that contain an interaction without the lower order effect")
        temp <- attr(Terms, "specials")$strata
        factors <- attr(Terms, "factors")[temp,]
        strata.interaction <- any(t(factors)*attr(Terms, "order") >1)
    } else strata <- NULL

    n <- object$n[1]

    # We need Y, X, istate, and id, also strata if present
    # If it was a Surv2() call, istate will be in the model frame and not
    #  explicit.  We have to get the model frame
    mf <- model.frame(object) 
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (is.null(offset)) offset.mean <- 0
    else {
        if (is.null(weights)) offset.mean <- mean(offset)
        else offset.mean <- sum(offset * (weights/sum(weights)))
    }
    X <- model.matrix(object, data=mf)
    Y <- object[['y']]
    if (is.null(Y)) {
        Y <- model.response(mf)
        if (is.null(object$timefix) || object$timefix) Y <- aeqSurv(Y)
    }
    if (nrow(Y) != object$n[1]) 
        stop("Failed to reconstruct the original data set")
 
    istate <- model.extract(mf, "istate")
    oldid <- model.extract(mf, "id")
    if (length(oldid) && ncol(Y)==3) position <- survflag(Y, oldid)
    else position <- NULL

    if (has.strata) {
        if (length(strata)==0) {
            if (length(stangle$vars) ==1) strata <- mf[[stangle$vars]]
            else strata <- strata(mf[, stangle$vars], shortlabel=TRUE)
        }
    } else strata <- NULL

    #deal with start time, by throwing out observations that end before then
    if (!missing(start.time)) {
        if (!is.numeric(start.time) || length(start.time) !=1
            || !is.finite(start.time))
            stop("start.time must be a single numeric value")
        toss <- which(Y[,ncol(Y)-1] <= start.time)
        if (length(toss)) {
            n <- nrow(Y)
            if (length(toss)==n) stop("start.time has removed all observations")
            Y <- Y[-toss,,drop=FALSE]
            X <- X[-toss,,drop=FALSE]
            weights <- weights[-toss]
            oldid <- oldid[-toss]
            istate <- istate[-toss]
            if (!is.null(strata)) strata <- strata[-toss]
        }
    }
    
    # Rebuild istate using the survcheck routine, as a double check
    # that the data set hasn't been modified
    mcheck <- survcheck2(Y, oldid, istate)
    transitions <- mcheck$transitions
    if (!identical(object$states, mcheck$states))
        stop("failed to rebuild the data set")
    istate <- mcheck$istate
    states <- object$states
    nstate <- length(states)

    #
    # Build the X matrix for the prototype subject(s)
    #
    if (missing(newdata)) {
        stop("a newdata argument is required for multistate model predictions")
    }
    Terms2 <- delete.response(Terms)

    # A named vector or list is allowed, sometimes used for a "1 row" newdata
    #  in my own code. (But I never should have done that, and make sure to
    #  never mention it's legality to users.)
    if (!inherits(newdata, "data.frame")) {
        if (is.list(newdata)) newdata <- data.frame(newdata)
        else if (is.numeric(newdata)) {
            if (is.null(names(newdata))) {
                stop("Newdata argument must be a data frame")
            }
            newdata <- data.frame(as.list(newdata), stringsAsFactors=FALSE)
        }  
    }

    # Now we have a valid newdata --- almost
    # If there were shared hazards, then phZ covariates will have been omitted
    #  from newdata by the user (the values in newdata won't be used), but
    #  the model.frame and model.matrix routines will be unhappy without it.
    # However, it turns out to be tricky to create a sensible covariate when you
    #  have the terms and model frame, but not the original data: imagine
    #  ns(df=5, x=zed).  The code below works for simple names or factors, 
    #  anything else the user will need to put in their own dummy.
    if (!is.null(object$phZ)) {
        # Zindex = which rows of cmap = cols of final X, are "Z" variables?
        Zindex <- match(rownames(object$phZ), rownames(object$cmap))
        # this code works for simple variables and factors
        vname <- all.vars(Terms2)
        # covert from Splus style assign to R style assign
        asn  <- rep(seq(along=object$assign), sapply(object$assign, length))
        Zname <- unique((names(object$assign))[asn[Zindex]])
        # Zname = list of variables we hope is in newdata
        for (k in 1:length(Zname)) {
            if (is.na(match(Zname[k], names(newdata)))) {
                # the name isn't there
                if (!is.null(object$xlevels) && 
                    !is.na(kk <- match(Zname[k], names(object$xlevels)))){
                    temp <- object$xlevels[[kk]]
                    dummy <- data.frame(x= factor(temp[1], temp))
                    names(dummy) <- Zname[k]
                    newdata <- cbind(newdata, dummy)
                }
                else if (!is.na(kk <- match(Zname[k], names(object$means)))) {
                    dummy <- data.frame(x= unname(object$means[kk]))
                    names(dummy) <- Zname[k]
                    newdata <- cbind(newdata, dummy)
                }
            }
        }
    }

    if (has.strata) {
        found.strata <- TRUE
        tempenv <- new.env(, parent=emptyenv())
        assign("strata", function(..., na.group, shortlabel, sep)
            list(...), envir=tempenv)
        assign("list", list, envir=tempenv)
        for (svar in stangle$vars) {
            temp <- try(eval(parse(text=svar), newdata, tempenv),
                        silent=TRUE)
            if (!is.list(temp) || 
                any(unlist(lapply(temp, class))== "function"))
                found.strata <- FALSE
        }
        
        if (!found.strata) {
            ss <- untangle.specials(Terms2, "strata")
            Terms2 <- Terms2[-ss$terms]
        }
    }
        
    # create a model frame and X matrix for the new data
    tcall <- Call[c(1, match(c("formula", "na.action"), 
                                     names(Call), nomatch=0))]
    tcall$data <- newdata
    tcall$formula <- Terms2
    tcall$xlev <- object$xlevels[match(attr(Terms2,'term.labels'),
                                       names(object$xlevels), nomatch=0)]
    tcall[[1L]] <- quote(stats::model.frame)
 
    mf2 <- eval(tcall)
    if (nrow(mf2) ==0)
        stop("all rows of newdata have missing values")
    x2 <- model.matrix(object, mf2)

    offset2 <- model.offset(mf2)
    if (length(offset2)==0 ) offset2 <- 0

    if (has.strata && found.strata) { #pull them off
        temp <- untangle.specials(Terms2, 'strata')
        strata2 <- strata(mf2[temp$vars], shortlabel=TRUE)
        strata2 <- factor(strata2, levels=levels(strata))
        if (any(is.na(strata2)))
            stop("New data set has strata levels not found in the original")
        # An expression like age:strata(sex) will have temp$vars= "strata(sex)"
        #  and temp$terms = integer(0).  This does not work as a subscript
        if (length(temp$terms) >0) Terms2 <- Terms2[-temp$terms]
    }
    else strata2 <- factor(rep(0, nrow(mf2)))

    offset2 <- model.offset(mf2)
    if (length(offset2)==0 ) offset2 <- 0

    # Let the survfitAJ routine do the work of creating the
    #  overall counts (n.risk, etc).  The rest of this code then
    #  replaces the surv and hazard components.
    if (missing(start.time)) start.time <- min(Y[,2], 0)

    if (is.null(weights)) weights <- rep(1.0, nrow(Y))
    if (is.null(strata))  tempstrat <- rep(1L, nrow(Y))
    else                  tempstrat <- strata


    cifit <- survfitAJ(as.factor(tempstrat), Y, weights, 
                        id= oldid, istate = istate, se.fit=FALSE, 
                        start.time=start.time, p0=p0, time0= time0)
    # if the user did not provide p0, survfitAJ fills one in for us
    p0 <- cifit$p0
    if (is.null(cifit$strata)) cifit.index <- 1:length(cifit$time)
    else cifit.index <- rep(1:length(cifit$strata), cifit$strata)

    # Create matrices for the final cumhaz and pstate
    ntime <- length(cifit$time)
    ndata <- nrow(x2)
    baseline <- object$smap[1,]
    nhaz <- length(unique(baseline))
    cumhaz <- array(0, dim=c(ntime, ndata, nhaz), 
                dimnames= list(NULL, NULL, names(baseline)[unique(baseline)]))
    pstate <- array(0, dim=c(ntime, ndata, nstate),
                    dimnames= list(NULL, NULL, states))

    from.state <- as.numeric(sub(":[0-9]*$", "", names(baseline)))
    to.state   <- as.numeric(sub("^[0-9]*:", "", names(baseline)))
    hindex1 <- cbind(from.state, to.state)
    hindex2 <- match(baseline, unique(baseline))
    Afill <- function(lambda, h1, h2, nstate) {
        # create the A matrix from a set of hazards
        A <- matrix(0, nstate, nstate)
        A[h1] <- lambda[h2]
        diag(A) <- -rowSums(A)
        A
    }
    # this is a trick that sets the defaults for Afill, so that I don't
    #  have to pass these to multihaz
    temp <- formals(Afill)
    temp$h1 <- hindex1
    temp$h2 <- hindex2
    temp$nstate <- nstate
    formals(Afill) <- temp

    # For each row of newdata
    #   1. Update the X matrix, so that each column is recentered at this
    #   row of newdata, i.e., the curves for a subject with covariate x2[i,]
    #   2. Use stacker to create the expanded X and Y
    #   3. For any phZ variables, replace selected columns of the expanded X,
    #       compute the risk weight exp(XB) for each data row
    #      Create the phscale vector
    #   4. Iterate over time (within strata if present)
    #     a. compute dN(jk,t) and the risk sum Y(jk,t) for each hazard jk at
    #         this time, using both case weights and risk weights
    #     b. collapse shared hazards to get lambda(t)
    #     c. incement the cumulative hazard and the p(t)
    #     d. if multiple strata, rezero temps at the start of each stratum
    tempstrat <- as.integer(tempstrat)
    phscale <- model.matrix(~factor(object$smap[1,]) -1)
    for (idata in 1:ndata) {
        tempX <- scale(X, center=x2[idata,], scale=FALSE)
        xstack <- stacker(object$cmap, object$smap, as.integer(istate), tempX,
                          Y, mf=mf, states=object$states, dropzero=FALSE)
        X2 <- xstack$X
        if (!is.null(object$phZ)) {
            phZ <- object$phZ  # less typing below
            cmap <- object$cmap
            ctemp <- cmap[match(dimnames(phZ)[[1]], rownames(cmap)),, drop=FALSE]
            # Go through the individual subshared terms one at a time.
            # There are never very many (commonly < 10). 
            # If there are transitions with no covariates, the xstack$X matrix
            #  will have rows for them, none are relevant to this calculation.
            # Only look at xstack$subshare values that match names in phZ
            for (sub in dimnames(phZ)[[3]]) {
                irow <- which(xstack$subshare == sub) #relevant rows of X
                nr <- length(irow)
                for (j in xstack$submap[sub]) {
                    # j = columns for this shared subterm
                    k <- (ctemp[,j] >0)
                    if (any(k)) {
                        X2[irow, ctemp[k,j]] <- rep(phZ[k,j,sub], each=nr)
                        eta <- sum(phZ[k,j,sub] * object$coef[ctemp[k,j]])
                        browser()
                        phscale[j] <- exp(eta)
                    }
                }
            }
        }
        eta <- X2 %*% object$coefficients

        if (is.null(strata)) {
            kk <- xstack$rindex
            mfit <- multihaz(xstack$Y, X2, position[kk], weights[kk], exp(eta),
                             xstack$shared, ctype=ctype, stype=stype, Afill, 
                             phscale, p0 = p0, utime= cifit$time, nstate=nstate)
            cumhaz[,idata,] <- mfit$cumhaz
            pstate[,idata,] <- mfit$pstate
        }
        else {
            cumhaz <- array(0, dim=c(ntime, ndata, nhaz))
            pstate <- array(0, dim=c(ntime, ndata, nhaz))
            for (tstrat in 1:nstrat) {
                jj <- which(tempstrat[xstack$rindex] == tstrat) # rows of xstacl
                kk <- xstack$rindex[jj] # rows of data
                mfit <- multihaz(xstack$Y[jj,,drop= FALSE], X2[jj,, drop=FALSE],
                             position[kk], weights[kk], exp(eta[jj]), 
                             xstack$stran[jj], ctype=ctype, stype=stype, Afill, 
                             p0 = p0, 
                             utime= cifit$time[cifit.index== tstrat], nstate)
                cumhaz[cifit.index==i, idata,] <- mfit$cumhaz
                pstate[cifit.index==i, idata,] <- mfit$cumhaz
            }
        }
    }
    cifit$cumhaz <- drop(cumhaz)
    cifit$pstate <- drop(pstate)
    cifit$newdata <- newdata

    cifit$call <- Call
    class(cifit) <- c("survfitcoxms", "survfitms", "survfit")
    cifit
}

# Compute the cumulative hazard and probability in state functions 
#  This can only handle one external strata at a time, but we expect 
#  external strata to be rare in multistate models
multihaz <- function(y, x, position, weight, risk, transition, ctype, stype, 
                     Afill, phscale, p0, utime, nstate) {
    ny <- ncol(y)
    sort2 <- order(transition, y[,ny-1L]) -1L
    ntime <- length(utime)
    storage.mode(weight) <- "double"  #failsafe

    # this returns all of the counts we might desire.
    if (ny ==2) {
        fit <- .Call(Ccoxsurv1, utime, y, weight, sort2, transition, x, risk)
        cn <- fit$count  
        dim(cn) <- c(length(utime), fit$ntrans, 10) 
    }
    else {    
        sort1 <- order(transition, y[,1]) -1L
        fit <- .Call(Ccoxsurv2, utime, y, weight, sort1, sort2, position, 
                        transition, x, risk)
        cn <- fit$count  
        dim(cn) <- c(length(utime), fit$ntrans, 12) 
    }

    # cn is returned as a matrix since there is an allocMatrix C macro, but
    #  no allocArray macro.  So we first reset the dimensions.
    # The first dimension is time
    # Second is the transition
    # Third is the count type: 1-3 = at risk (unweighted, with case weights,
    #  with casewt * risk wt), 4-6 = events (unweighted, case, risk), 
    #  7-8 = censored, 9-10 = Efron, 11-12= unused

    # We will use events/(at risk) = cn[,,5]/cn[,,3] a few lines below; but
    #  avoid 0/0. If there is no one at risk there are no events, obviously.
    # cn[,,1] is the safer check since it is an integer, but if there are wts
    #  and a subject with weight=0 were the only one at risk, we need cn[,,2]
    # (Users should never use weights of 0, but someone, somewhere, will do it.)
    none.atrisk <- (cn[,,1]==0 | cn[,,2]==0)
    if (ctype ==1) {
        denom1 <- ifelse(none.atrisk, 1, cn[,,3])   # avoid a later 0/0
        denom2 <- ifelse(none.atrisk, 1, cn[,,3]^2)
    } else {
        denom1 <- ifelse(none.atrisk, 1, cn[,,9])
        denom2 <- ifelse(none.atrisk, 1, cn[,,10])
    }

    if (ctype==1) hazard <- cn[,,5]/ifelse(cn[,,3]<=0, 1, cn[,,3])
    else          hazard <- cn[,,5]/ifelse(cn[,,9] <=0, 1, cn[,,9])

    browser()
    hazard <- hazard * rep(phscale, each=nrow(hazard))
    
    cumhaz <- apply(hazard, 2, cumsum)
    pstate <- matrix(0, length(utime), nstate)
    phat <- p0
    for (i in 1:length(utime)) {
        if (all(hazard[i,] ==0)) pstate[i,] <- phat
        else {
            amat <- Afill(hazard[i,])
            if (stype==1) phat <- phat %*% (diag(nstate) + amat)
            else phat <- phat %*% survexpm(amat)
            pstate[i,] <- phat
        }
    }
    list(cumhaz = cumhaz, pstate=pstate, ntran= cn[,,3])
}
