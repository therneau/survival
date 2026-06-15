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
        stop("survfit.coxph has not been updated for strata")
        stangle = untangle.specials(Terms, "strata")  #used multiple times, later
        # Toss out strata terms in tfac before doing the test 1 line below, as
        #  strata end up in the model with age:strat(grp) terms or *strata() terms
        #  (There might be more than one strata term)
        for (i in temp) tfac <- tfac[,tfac[i,] ==0]  # toss out strata terms
    }
    if (any(tfac >1))
        stop("not able to create a curve for models that contain an interaction without the lower order effect")

    n <- object$n[1]
    if (!has.strata) strata <- NULL
    else strata <- object$strata

    if (has.strata) {
        temp <- attr(Terms, "specials")$strata
        factors <- attr(Terms, "factors")[temp,]
        strata.interaction <- any(t(factors)*attr(Terms, "order") >1)
    }

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
    nstate <- length(object$states)

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

    # Create matrices for the final cumhaz and pstate
    ntime <- length(cifit$time)
    ndata <- nrow(x2)
    baseline <- object$smap[1,]
    nhaz <- length(unique(baseline))
    cumhaz <- array(0, dim=c(ntime, ndata, nhaz))
    pstate <- array(0, dim=c(ntime, ndata, nhaz))
    B <- coef(object, matrix=TRUE)
    # For each row of newdata
    #   1. Update the X matrix, so that each column is recentered at this
    #   row of newdata, i.e., the curves for a subject with covariate x2[i,]
    #   2. Use stacker to create the expanded X and Y
    #   3. For any phZ variables, replace selected columns of the expanded X,
    #       compute the risk weight exp(XB) for each data row
    #   4. Iterate over time (within strata if present)
    #     a. compute dN(jk,t) and the risk sum Y(jk,t) for each hazard jk at
    #         this time, using both case weights and risk weights
    #     b. collapse shared hazards to get lambda(t)
    #     c. incement the cumulative hazard and the p(t)
    #     d. if multiple strata, rezero temps at the start of each stratum
    
    for (idata in 1:ndata) {
        tempX <- scale(X, center=x2[idata,], scale=FALSE)
        xstack <- stacker(object$cmap, object$smap, as.integer(istate), tempX,
                          Y, mf=mf, states=object$states, dropzero=FALSE)
        X2 <- xstack$X
        n.per.strat <- table(xstack$transition)
        if (!is.null(object$phZ)) {
            temp <- t(object$phZ)[as.integer(xstack$transition),]
            X2[, Zindex] <- temp
        }
        browser()
        eta <- X2 %*% B
        # move this to C code
        if (idata==1) { # the sort won't change across iterations, eta will
            # remember that C indices start at 0
            if (!is.null(strata)) {
                sort1 <- order(strata[xstack$rindex], xstack$Y[,1]) -1L
                sort2 <- order(strata[xstace$rindex], xstack$Y[,2]) -1L
            } else {
                sort1 <- order(xstack$Y[,1]) -1L
                sort2 <- order(xstack$Y[,2]) -1L
            }
        }   
        # the eta for each transition
        xeta <- eta[cbind(1:nrow(eta), xstack$transition)]
        hcount  <- .Call("coxcount3", xstack$Y, weights[xstack$rindex], 
                         exp(xeta), 
                         sort1, sort2, tempstrat[xstack$rindex], 
                         as.integer(xstack$stran), nhaz, cifit$time)
        browser()
    }

    if (is.null(strata)) {
        fit <- multihaz(Y, X, position, weights, risk, istrat, ctype, stype,
                        baselinecoef, hfill, x2, risk2, varmat, nstate, se.fit, 
                        cifit$p0, cifit$time)
        cifit$pstate <- fit$pstate
        cifit$cumhaz <- fit$cumhaz
    }
    else {
        if (is.factor(strata)) ustrata <- levels(strata)
        else                   ustrata <- sort(unique(strata))
        nstrata <- length(cifit$strata)
        itemp <- rep(1:nstrata, cifit$strata)
        timelist <- split(cifit$time, itemp)
        ustrata <- names(cifit$strata)
        tfit <- vector("list", nstrata)
        for (i in 1:nstrata) {
            indx <- which(strata== ustrata[i])  # divides the data
            tfit[[i]] <- multihaz(Y[indx,,drop=F], X[indx,,drop=F],
                                  position[indx], weights[indx], risk[indx],
                                  istrat[indx], ctype, stype, baselinecoef, hfill,
                                  x2, risk2, varmat, nstate, se.fit,
                                  cifit$p0[i,], timelist[[i]])
        }

        # do.call(rbind) doesn't work for arrays, it loses a dimension
        ntime <- length(cifit$time)
        cifit$pstate <- array(0., dim=c(ntime, dim(tfit[[1]]$pstate)[2:3]))
        cifit$cumhaz <- array(0., dim=c(ntime, dim(tfit[[1]]$cumhaz)[2:3]))
        rtemp <- split(seq(along=cifit$time), itemp)
        for (i in 1:nstrata) {
            cifit$pstate[rtemp[[i]],,] <- tfit[[i]]$pstate
            cifit$cumhaz[rtemp[[i]],,] <- tfit[[i]]$cumhaz
        }
    }

    cifit$newdata <- newdata

    cifit$call <- Call
    class(cifit) <- c("survfitcoxms", "survfitms", "survfit")
    cifit
}
# Compute the hazard  and survival functions 
multihaz <- function(y, x, position, weight, risk, istrat, ctype, stype, 
                     bcoef, hfill, x2, risk2, vmat, nstate, se.fit, p0, utime) {
    ny <- ncol(y)
    sort2 <- order(istrat, y[,ny-1L]) -1L
    ntime <- length(utime)
    storage.mode(weight) <- "double"  #failsafe

    # this returns all of the counts we might desire.
    if (ny ==2) {
        fit <- .Call(Ccoxsurv1, utime, y, weight, sort2, istrat, x, risk)
        cn <- fit$count  
        dim(cn) <- c(length(utime), fit$ntrans, 10) 
    }
    else {    
        sort1 <- order(istrat, y[,1]) -1L
        fit <- .Call(Ccoxsurv2, utime, y, weight, sort1, sort2, position, 
                        istrat, x, risk)
        cn <- fit$count  
        dim(cn) <- c(length(utime), fit$ntrans, 12) 
    }
    # cn is returned as a matrix since there is an allocMatrix C macro, but
    #  no allocArray macro.  So we first reset the dimensions.
    # The first dimension is time
    # Second is the transition, same order as columns of bcoef
    # Third is the count type: 1-3 = at risk (unweighted, with case weights,
    #  with casewt * risk wt), 4-6 = events (unweighted, case, risk), 
    #  7-8 = censored events, 9-10 = censored, 11-12 = Efron

    # We will use events/(at risk) = cn[,,5]/cn[,,3] a few lines below; avoid 0/0
    # If there is no one at risk there are no events, obviously.
    # cn[,,1] is the safer check since it is an integer, but if there are weights
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

    # We want to avoid 0/0. If there is no one at risk (denominator) then
    # by definition there will be no events (numerator), and that element of
    # the hazard is by definintion also 0.
    if (any(duplicated(bcoef[1,]))) {
        # there are shared hazards: we have to collapse and then expand
        if (all(bcoef[1,] == bcoef[1,1])) design <- matrix(1, nrow=ncol(bcoef))
        else design <- model.matrix(~factor(zed) -1, data.frame(zed=bcoef[1,]))
        colnames(design) <- 1:ncol(design)  # easier to read when debuggin
        events <- cn[,,5] %*% design
        if (ctype==1) atrisk <- cn[,,3]  %*% design
        else          atrisk <- cn[,,9] %*% design
        basehaz <- events/ifelse(atrisk<=0, 1, atrisk)
        hazard <- basehaz[,bcoef[1,]] * rep(bcoef[2,], each=nrow(basehaz))
    }                                  
    else {
        if (ctype==1) hazard <- cn[,,5]/ifelse(cn[,,3]<=0, 1, cn[,,3])
        else          hazard <- cn[,,5]/ifelse(cn[,,9] <=0, 1, cn[,,9])
    }

    # Expand the result, one "hazard set" for each row of x2
    nx2 <- nrow(x2)
    h2 <- array(0, dim=c(nrow(hazard), nx2, ncol(hazard)))
    S <- double(nstate)  # survival at the current time
    S2 <- array(0, dim=c(nrow(hazard), nx2, nstate))
 
    H <- matrix(0, nstate, nstate)
    if (stype==2) {
        H[hfill] <- colMeans(hazard)    # dummy H to drive esetup
        diag(H) <- diag(H) -rowSums(H)
        esetup <- survexpmsetup(H)
    }

    for (i in 1:nx2) {
        h2[,i,] <- apply(hazard * rep(risk2[i,], each=ntime), 2, cumsum)
        if (FALSE) {  # if (se.fit) eventually
            d1 <- fit$xbar - rep(x[i,], each=nrow(fit$xbar))
            d2 <- apply(d1*hazard, 2, cumsum)
            d3 <- rowSums((d2%*% vmat) * d2)
            v2[jj,] <- (apply(varhaz[jj,],2, cumsum) + d3) * (risk2[i])^2
        }

        S <- p0
        for (j in 1:ntime) {
            if (any(hazard[j,] > 0)) { # hazard =0 for censoring times
                H[,] <- 0.0
                H[hfill] <- hazard[j,] *risk2[i,]
                if (stype==1) {
                    diag(H) <- pmax(0, 1.0 - rowSums(H))
                    S <- as.vector(S %*% H)  # don't keep any names
                }
                else {
                    diag(H) <- 0.0 - rowSums(H)
                    #S <- as.vector(S %*% expm(H))  # dgeMatrix issue
                    S <- as.vector(S %*% survexpm(H, 1, esetup))
                }
            }
            S2[j,i,] <- S
        }
    }
    rval <- list(time=utime, xgrp=rep(1:nx2, each=nrow(hazard)),
                 pstate=S2, cumhaz=h2)
    #if (se.fit) rval$varhaz <- v2
    rval
}
