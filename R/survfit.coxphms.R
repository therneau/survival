# Automatically generated from the noweb directory
survfit.coxphms <-
function(formula, newdata, se.fit=FALSE, conf.int=.95, individual=FALSE,
         stype=2, ctype, 
         conf.type=c("log", "log-log", "plain", "none", "logit", "arcsin"),
         censor=TRUE, start.time, id, influence=FALSE,
         na.action=na.pass, type, p0=NULL, ...) {

    Call <- match.call()
    Call[[1]] <- as.name("survfit")  #nicer output for the user
    object <- formula     #'formula' because it has to match survfit
    se.fit <- FALSE   #still to do
    if (missing(newdata))
        stop("multi-state survival requires a newdata argument")
    if (!missing(id)) 
        stop("using a covariate path is not supported for multi-state")
    temp <- object$smap["(Baseline)",] 
    baselinecoef <- rbind(temp, coef= 1.0)
    if (any(duplicated(temp))) {
        # We have shared hazards 
        # If there are k duplicates, then the last k coefficient in beta are
        #  the coefs that match them
        idup <- duplicated(temp)
        ncoef <- length(object$coefficients)
        ndup <- sum(idup)
        # shared baseline coefficients are last in the coefficient vector
        i <- seq(to = ncoef, length=ndup)
        baselinecoef[2, idup] <- exp(object$coefficients[i])

        # which rows of cmap point to scale coefs?
        phbase <- apply(object$cmap, 1, function(i) any(i > ncoef-ndup))
    }
    else phbase <- rep(FALSE, nrow(object$cmap))
      
    # process options, set up Y and the model frame for the original data
    Terms  <- terms(object)
    robust <- !is.null(object$naive.var)   # did the coxph model use robust var?

    if (!is.null(attr(object$terms, "specials")$tt))
        stop("The survfit function can not process coxph models with a tt term")

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

    tfac <- attr(Terms, 'factors')
    temp <- attr(Terms, 'specials')$strata 
    has.strata <- !is.null(temp)
    if (has.strata) {
        stangle = untangle.specials(Terms, "strata")  #used multiple times, later
        # Toss out strata terms in tfac before doing the test 1 line below, as
        #  strata end up in the model with age:strat(grp) terms or *strata() terms
        #  (There might be more than one strata term)
        for (i in temp) tfac <- tfac[,tfac[i,] ==0]  # toss out strata terms
    }
    if (any(tfac >1))
        stop("not able to create a curve for models that contain an interaction without the lower order effect")

    Terms <- object$terms
    n <- object$n[1]
    if (!has.strata) strata <- NULL
    else strata <- object$strata

    if (!missing(individual)) warning("the `id' option supersedes `individual'")
    missid <- missing(id) # I need this later, and setting id below makes
                          # "missing(id)" always false

    if (!missid) individual <- TRUE
    else if (missid && individual) id <- rep(0L,n)  #dummy value
    else id <- NULL

    if (individual & missing(newdata)) {
        stop("the id option only makes sense with new data")
    }
    if (has.strata) {
        temp <- attr(Terms, "specials")$strata
        factors <- attr(Terms, "factors")[temp,]
        strata.interaction <- any(t(factors)*attr(Terms, "order") >1)
    }
    coxms <- inherits(object, "coxphms")
    if (coxms || is.null(object$y) || is.null(object[['x']]) ||
        !is.null(object$call$weights) || !is.null(object$call$id) ||
        (has.strata && is.null(object$strata)) ||
        !is.null(attr(object$terms, 'offset'))) {
        
        mf <- stats::model.frame(object)
        }
    else mf <- NULL  #useful for if statements later
    position <- NULL
    Y <- object[['y']]
    if (is.null(mf)) {
        weights <- object$weights  # let offsets/weights be NULL until needed
        offset <- NULL
        offset.mean <- 0
        X <- object[['x']]
    }
    else {
        weights <- model.weights(mf)
        offset <- model.offset(mf)
        if (is.null(offset)) offset.mean <- 0
        else {
            if (is.null(weights)) offset.mean <- mean(offset)
            else offset.mean <- sum(offset * (weights/sum(weights)))
        }
        X <- model.matrix.coxph(object, data=mf)
        if (is.null(Y) || coxms) {
            Y <- model.response(mf)
            if (is.null(object$timefix) || object$timefix) Y <- aeqSurv(Y)
        }
        oldid <- model.extract(mf, "id")
        if (length(oldid) && ncol(Y)==3) position <- survflag(Y, oldid)
        else position <- NULL
        if (!coxms && (nrow(Y) != object$n[1])) 
            stop("Failed to reconstruct the original data set")
        if (has.strata) {
            if (length(strata)==0) {
                if (length(stangle$vars) ==1) strata <- mf[[stangle$vars]]
                else strata <- strata(mf[, stangle$vars], shortlabel=TRUE)
            }
        }

    }
    istate <- model.extract(mf, "istate")

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
        }
    }
    
    # expansion of the X matrix with stacker, set up shared hazards
    # Rebuild istate using the survcheck routine, as a double check
    # that the data set hasn't been modified
    mcheck <- survcheck2(Y, oldid, istate)
    transitions <- mcheck$transitions
    if (!identical(object$states, mcheck$states))
        stop("failed to rebuild the data set")
    if (is.null(istate)) istate <- mcheck$istate
    else {
        # if istate has unused levels, mcheck$istate won't have them so they
        #  need to be dropped.
        istate <- factor(istate, object$states) 
        # a new level in state should only happen if someone has mucked up the
        #  data set used in the coxph fit
        if (any(is.na(istate))) stop("unrecognized initial state, data changed?")
    }

    # Let the survfitCI routine do the work of creating the
    #  overall counts (n.risk, etc).  The rest of this code then
    #  replaces the surv and hazard components.
    if (missing(start.time)) start.time <- min(Y[,2], 0)

    if (is.null(weights)) weights <- rep(1.0, nrow(Y))
    if (is.null(strata))  tempstrat <- rep(1L, nrow(Y))
    else                  tempstrat <- strata

    cifit <- survfitCI(as.factor(tempstrat), Y, weights, 
                            id= oldid, istate = istate, se.fit=FALSE, 
                            start.time=start.time, p0=p0)

    # For computing the  actual estimates it is easier to work with an
    #  expanded data set.
    # Replicate actions found in the coxph-multi-X chunk
    # Note the dropzero=FALSE argument: if there is a transition with no 
    #  covariates we still need it expanded; this differs from coxph.
    # A second differnence is tstrata: force stacker to think that every
    #  transition is a unique hazard, so that it does proper expansion.
    cluster <- model.extract(mf, "cluster")
    tstrata <- object$smap
    tstrata[1,] <- 1:ncol(tstrata)
    xstack <- stacker(object$cmap, tstrata, as.integer(istate), X, Y,
                      as.integer(strata),
                      states= object$states, dropzero=FALSE)
    if (length(position) >0)
        position <- position[xstack$rindex]   # id was required by coxph
    X <- xstack$X
    Y <- xstack$Y
    strata <- strata[xstack$rindex]  # strat in the model, other than transitions
    transition <- xstack$transition
    istrat <- xstack$strata
    if (length(offset)) offset <- offset[xstack$rindex]
    if (length(weights)) weights <- weights[xstack$rindex]
    if (length(cluster)) cluster <- cluster[xstack$rindex]
    oldid <- oldid[xstack$rindex]
    if (robust & length(cluster)==0) cluster <- oldid

    # risk scores, mf2, and x2
    if (length(object$means) ==0) { # a model with only an offset term
        # Give it a dummy X so the rest of the code goes through
        #  (This case is really rare)
        # se.fit <- FALSE
        X <- matrix(0., nrow=n, ncol=1)
        if (is.null(offset)) offset <- rep(0, n)
        xcenter <- offset.mean
        coef <- 0.0
        varmat <- matrix(0.0, 1, 1)
        risk <- rep(exp(offset- offset.mean), length=n)
    }
    else {
        varmat <- object$var
        beta <- ifelse(is.na(object$coefficients), 0, object$coefficients)
        xcenter <- sum(object$means * beta) + offset.mean
        if (!is.null(object$frail)) {
           keep <- !grepl("frailty(", dimnames(X)[[2]], fixed=TRUE)
           X <- X[,keep, drop=F]
        }
            
        if (is.null(offset)) risk <- c(exp(X%*% beta - xcenter))
        else     risk <- c(exp(X%*% beta + offset - xcenter))
    }
    if (missing(newdata)) {
        # If the model has interactions, print out a long warning message.
        #  People may hate it, but I don't see another way to stamp out these
        #  bad curves without backwards-incompatability.  
        # I probably should complain about factors too (but never in a strata
        #   or cluster term).
        if (any(attr(Terms, "order") > 1) )
            warning("the model contains interactions; the default curve based on columm means of the X matrix is almost certainly not useful. Consider adding a newdata argument.")
        
        if (length(object$means)) {
            mf2 <- as.list(object$means)   #create a dummy newdata
            names(mf2) <- names(object$coefficients)
            mf2 <- as.data.frame(mf2)
            x2 <- matrix(object$means, 1)
        }
        else { # nothing but an offset
            mf2 <- data.frame(X=0)
            x2 <- 0
        }
        offset2 <- 0
        found.strata <- FALSE  
    }
    else {
        if (!is.null(object$frail))
            stop("Newdata cannot be used when a model has frailty terms")

        Terms2 <- Terms 
        if (!individual)  {
            Terms2 <- delete.response(Terms)
            y2 <- NULL  # a dummy to carry along, for the call to coxsurv.fit
        }
        if (is.vector(newdata, "numeric")) {
            if (individual) stop("newdata must be a data frame")
            if (is.null(names(newdata))) {
                stop("Newdata argument must be a data frame")
            }
            newdata <- data.frame(as.list(newdata), stringsAsFactors=FALSE)
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

        tcall <- Call[c(1, match(c('id', "na.action"), 
                                     names(Call), nomatch=0))]
        tcall$data <- newdata
        tcall$formula <- Terms2
        tcall$xlev <- object$xlevels[match(attr(Terms2,'term.labels'),
                                           names(object$xlevels), nomatch=0)]
        tcall[[1L]] <- quote(stats::model.frame)
        mf2 <- eval(tcall)
    }
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

    if (!robust) cluster <- NULL
    if (individual) {
        if (missing(newdata)) 
            stop("The newdata argument must be present when individual=TRUE")
        if (!missid) {  #grab the id variable
            id2 <- model.extract(mf2, "id")
            if (is.null(id2)) stop("id=NULL is an invalid argument")
            }
        else id2 <- rep(1, nrow(mf2))
        
        x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]  #no intercept
        if (length(x2)==0) stop("Individual survival but no variables")

        offset2 <- model.offset(mf2)
        if (length(offset2) ==0) offset2 <- 0
                     
        y2 <- model.extract(mf2, 'response')
        if (attr(y2,'type') != type)
            stop("Survival type of newdata does not match the fitted model")
        if (attr(y2, "type") != "counting")
            stop("Individual=TRUE is only valid for counting process data")
        y2 <- y2[,1:2, drop=F]  #throw away status, it's never used
    }
    else if (missing(newdata)) {
        if (has.strata && strata.interaction)
            stop ("Models with strata by covariate interaction terms require newdata")
        offset2 <- 0
        if (length(object$means)) {
            x2 <- matrix(object$means, nrow=1, ncol=ncol(X))
        } else {
            # model with only an offset and no new data: very rare case 
            x2 <- matrix(0.0, nrow=1, ncol=1)   # make a dummy x2
        }
    } else {
        offset2 <- model.offset(mf2)
        if (length(offset2)==0 ) offset2 <- 0
        # a model with only an offset, but newdata containing a value for it
        if (length(object$means)==0) x2 <- 0
        else x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]  #no intercept
    }

    if (has.strata && !is.null(mf2[[stangle$vars]])){
        mf2 <- mf2[is.na(match(names(mf2), stangle$vars))]
        mf2 <- unique(mf2)
        x2 <- unique(x2)
    }
    temp <- coef(object, matrix=TRUE)[!phbase,,drop=FALSE] # ignore missing coefs
    # temp will be a matrix of coefficients, with ncol = number of transtions
    #  and nrow = the covariate set of a "normal" Cox model.
    # x2 will have one row per desired curve and one col per 'normal' covariate.
    risk2 <- exp(x2 %*% ifelse(is.na(temp), 0, temp) - xcenter)
    # risk2 has a risk score with rows= curve and cols= transition
    # make the expansion map.  
    # The H matrices we will need are nstate by nstate, at each time, with
    # elements that are non-zero only for observed transtions.
    states <- object$states
    nstate <- length(states)
    from <- as.numeric(sub(":.*$", "", colnames(object$smap)))
    to   <- as.numeric(sub("^.*:", "", colnames(object$smap)))
    hfill <- cbind(from, to)

    if (individual) {
        stop("time dependent survival curves are not supported for multistate")
    }
    ny <- ncol(Y)
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
    cifit$newdata <- mf2

    cifit$call <- Call
    class(cifit) <- c("survfitms", "survfit")
    cifit
}
# Compute the hazard  and survival functions 
multihaz <- function(y, x, position, weight, risk, istrat, ctype, stype, 
                     bcoef, hfill, x2, risk2, vmat, nstate, se.fit, p0, utime) {
    sort2 <- order(istrat, y[,2]) -1L
    ntime <- length(utime)
    storage.mode(weight) <- "double"  #failsafe

    # this returns all of the counts we might desire.
    if (ncol(y) ==2) 
    if (ncol(y) ==2) {
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
        denom1 <- ifelse(none.atrisk, 1, cn[,,11])
        denom2 <- ifelse(none.atrisk, 1, cn[,,12])
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
        else          atrisk <- cn[,,11] %*% design
        basehaz <- events/ifelse(atrisk<=0, 1, atrisk)
        hazard <- basehaz[,bcoef[1,]] * rep(bcoef[2,], each=nrow(basehaz))
    }                                  
    else {
        if (ctype==1) hazard <- cn[,,5]/ifelse(cn[,,3]<=0, 1, cn[,,3])
        else          hazard <- cn[,,5]/ifelse(cn[,,11] <=0, 1, cn[,,11])
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
