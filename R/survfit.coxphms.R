# Automatically generated from the noweb directory
survfit.coxphms <-
function(formula, newdata, se.fit=TRUE, conf.int=.95, individual=FALSE,
         stype=2, ctype, 
         conf.type=c("log", "log-log", "plain", "none", "logit", "arcsin"),
         censor=TRUE, start.time, id, influence=FALSE,
         na.action=na.pass, type, ...) {

    Call <- match.call()
    Call[[1]] <- as.name("survfit")  #nicer output for the user
    object <- formula     #'formula' because it has to match survfit
    se.fit <- FALSE   #still to do

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
    if (!has.strata) strata <- rep(0L, n)
    else strata <- object$strata

    missid <- missing(id) # I need this later, and setting id below makes
                          # "missing(id)" always false
    if (!missid & !missing(individual))
        warning("the `id' option supersedes `individual'")

    if (!missid) individual <- TRUE
    else if (missid && individual) id <- rep(0,n)  #dummy value
    else id <- NULL

    if (individual & missing(newdata)) {
        stop("the id and/or individual options only make sense with new data")
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
        weights <- rep(1., n)
        offset <- rep(0., n)
        X <- object[['x']]
    }
    else {
        weights <- model.weights(mf)
        if (is.null(weights)) weights <- rep(1.0, n)
        offset <- model.offset(mf)
        if (is.null(offset)) offset <- rep(0., n)
        X <- model.matrix.coxph(object, data=mf)
        if (is.null(Y)) Y <- aeqSurv(model.response(mf))
        oldid <- model.extract(mf, "id")
        if (length(oldid) && ncol(Y)==3) position <- survflag(Y, oldid)
        else position <- NULL
        if (nrow(Y) != object$n[1]) 
            stop("Failed to reconstruct the original data set")
        if (has.strata) {
            if (length(strata)==0) {
                if (length(stangle$vars) ==1) strata <- mf[[stangle$vars]]
                else strata <- strata(mf[, stangle$vars], shortlabel=TRUE)
            }
        }

    }
    # Rebuild istate using the mcheck routine
    istate <- model.extract(mf, "istate")
    mcheck <- survcheck2(Y, oldid, istate)
    transitions <- mcheck$transitions
    istate <- mcheck$istate
    if (!identical(object$states, mcheck$states))
        stop("failed to rebuild the data set")

    # Let the survfitCI routine do the work of creating the
    #  overall counts (n.risk, etc).  The rest of this code then
    #  replaces the surv and hazard components.
    if (missing(start.time)) start.time <- min(Y[,2], 0)
    # If the data has absorbing states (ones with no transitions out), then
    #  remove those rows first since they won't be in the final output.
    t2 <- transitions[, is.na(match(colnames(transitions), "(censored)"))]
    absorb <- row.names(t2)[rowSums(t2)==0]
    if (length(absorb)) droprow <- istate %in% absorb  else droprow <- FALSE
    if (any(droprow)) {
        j <- which(!droprow)
        cifit <- survfitCI(as.factor(strata[j]), Y[j,], weights[j], oldid[j], 
                           istate[j], stype=stype,
                           ctype=ctype, se.fit=FALSE, start.time=start.time)
        }
    else cifit <- survfitCI(as.factor(strata), Y, weights, oldid, istate, 
                            stype=stype, ctype=ctype, se.fit=FALSE, 
                            start.time=start.time)

    # For computing the  actual estimates it is easier to work with an
    #  expanded data set.
    # Replicate actions found in the coxph-multi-X chunk,
    cluster <- model.extract(mf, "cluster")
    xstack <- stacker(object$cmap, as.integer(istate), X, Y,
                      as.integer(strata),
                      states= object$states)
    position <- position[xstack$rindex]   # id was required by coxph
    X <- xstack$X
    Y <- xstack$Y
    strata <- strata[xstack$rindex]
    transition <- xstack$transition
    if (length(offset)) offset <- offset[xstack$rindex]
    if (length(weights)) weights <- weights[xstack$rindex]
    if (length(cluster)) cluster <- cluster[xstack$rindex]
    oldid <- oldid[xstack$rindex]
    if (robust & length(cluster)==0) cluster <- oldid
    if (length(object$means) ==0) { # a model with only an offset term
        # Give it a dummy X so the rest of the code goes through
        #  (This case is really rare)
        # se.fit <- FALSE
        X <- matrix(0., nrow=n, ncol=1)
        xcenter <- mean(offset)
        coef <- 0.0
        varmat <- matrix(0.0, 1, 1)
        risk <- rep(exp(offset- mean(offset)), length=n)
    }
    else {
        varmat <- object$var
        beta <- ifelse(is.na(object$coefficients), 0, object$coefficients)
        xcenter <- sum(object$means * beta)+ mean(offset)
        if (!is.null(object$frail)) {
           keep <- !grepl("frailty(", dimnames(X)[[2]], fixed=TRUE)
           X <- X[,keep, drop=F]
        }
            
        risk <- c(exp(X%*% beta + offset - xcenter))
    }
    if (missing(newdata)) {
        # If the model has factor predictors or it has interactions, print
        #  out a long warning message.  People are going to hate it, but
        #  I don't see another way to stamp out these bad curves without
        #  backwards-incompatability.  Don't complain about factors in any
        #  strata or cluster terms, however.
        sdrop <- c(attr(Terms, "specials")$strata, attr(Terms, "specials")$cluster)
        if (is.null(sdrop)) dc <- attr(Terms, "dataClasses")
        else dc <- attr(Terms, "dataClasses")[-sdrop]
        if (any(attr(Terms, "order") > 1) || any(dc %in% c("factor", "character")))
            warning("the model contains factor variables and/or interactions; the default curve based on columm means of the X matrix is almost certainly not useful. Consider adding a newdata argument.")
        
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
        if (!individual)  Terms2 <- delete.response(Terms)
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
        if (length(offset2) >0) offset2 <- offset2 
        else offset2 <- 0
        x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]  #no intercept
    }
    risk2 <- exp(x2 %*% coef(object, type="matrix") - xcenter)
    if (individual) {
        stop("time dependent survival curves not yet supported for multistate")
        result <- coxsurv.fit2(ctype, stype, se.fit, varmat, cluster, start.time,
                               object$cmap, object$transitions, object$states,
                               Y, X, weights, risk, position, strata, oldid,
                               transition, y2, x2, risk2, strata2, id2)
                              
    } else {
        if (is.null(cifit$strata)) p0 <- cifit$pstate[1,, drop=FALSE]
        else {
            last <- cumsum(cifit$strata)  # last obs of each strata
            first<- 1 + c(0, last[-length(last)])
            p0 <- cifit$pstate[first,, drop=FALSE]
        }
        cifit <- coxsurv.fit2(ctype, stype, se.fit, varmat, cluster, start.time,
                               object$cmap, object$transitions, object$states,
                               Y, X, weights, risk, position, strata, oldid,
                               transition, y2, x2, risk2, cifit=cifit)

        cifit$newdata <- mf2
    }


    cifit$call <- Call
    class(cifit) <- c("survfitms", "survfit")
    cifit
}
coxsurv.fit2 <- function (ctype, stype, se.fit, varmat, cluster, 
                          start.time, cmap, tmat, states,
                          y, x, weights, risk, position, strata, id,
                          transition, y2, x2, risk2, strata2, id2, cifit) {
    # args are the options (ctype, stype, se.fit), args info from the prior fit
    # (varmat, ..., tmat), original data (Y, ..., transition), and data for the
    # new subjects

    if (length(strata)==0) strata <- rep(0L, nrow(y))

    if (is.factor(strata)) ustrata <- levels(strata)
    else                   ustrata <- sort(unique(strata))
    nstrata <- length(ustrata)

    # make the expansion map.  
    #  cmap[1,] will contain integers 1, 2,... which match the values in
    # the transtion vector, which in turn is the set of hazard functions that
    # come back from the .Call
    #  The H matrices we will need are nstate by nstate, at each time, with
    # elements that are non-zero only for observed transtions.  Some elements
    # may be the same: cmat[1,] can have repeats.
    nstate <- length(states)
    tmat <- tmat[,is.na(match(colnames(tmat), "(censored)")), drop=FALSE]
    from <- row(tmat)[tmat>0]  # tmat contains fit$transitions matrix
    from <- match(rownames(tmat), states)[from]  # actual row of H
    to   <- col(tmat)[tmat>0]
    to   <- match(colnames(tmat), states)[to]    # actual col of H
    hfill <- cbind(from, to)

    storage.mode(position) <- "integer"  # failsafe
    storage.mode(weights) <- "double"
 
    if (nstrata==1) {
        temp <- multihaz(y, x, position, weights, risk, transition,
                                  ctype, stype, hfill, cmap[1,], 
                                  x2, risk2, varmat, nstate, se.fit, 
                                  cifit$pstate[1,], cifit$time)
        cifit$pstate <- temp$pstate
        cifit$cumhaz <- temp$cumhaz
    } 
    else {
        itemp <- rep(1:nstrata, cifit$strata)
        timelist <- split(cifit$time, itemp)
        firstrow <- match(1:nstrata, itemp)
        ustrata <- names(cifit$strata)
        survlist <- vector("list", nstrata)
        for (i in 1:nstrata) {
            indx <- which(strata== ustrata[i])  # divides the data
            survlist[[i]] <- multihaz(y[indx,,drop=F], x[indx,,drop=F],
                                  position[indx], weights[indx], risk[indx],
                                  transition[indx], ctype, stype, hfill,
                                  cmap[1,], x2, risk2, varmat, nstate, se.fit, 
                                  cifit$pstate[firstrow[i],], timelist[[i]])
                                  
            }
        cifit$pstate <- do.call(rbind, lapply(survlist, function(x) x$pstate))
        cifit$cumhaz <- do.call(rbind, lapply(survlist, function(x) x$cumhaz))
    }
    cifit
}
# Compute the hazard  and survival functions 
multihaz <- function(y, x, position, weight, risk, transition, ctype, stype, 
                     hfill, cmap, x2, risk2, vmat, nstate, se.fit, p0, utime) {
    if (ncol(y) ==2) {
       sort1 <- seq.int(nrow(y))   # any order will be the same
       y <- cbind(0.0, y)          # add a start.time column
    }
    else sort1 <- order(transition, y[,1]) -1L
    sort2 <- order(transition, y[,2]) -1L
    ntime <- length(utime)

    # this returns all of the counts we might desire.
    fit <- .Call(Ccoxsurv2, utime, y, weight, sort1, sort2, position, 
                        transition, x, risk)
    cn <- fit$count  # 1-3 = at risk, 4-6 = events, 7-8 = censored events
                     # 9-10 = censored, 11-12 = Efron, 13-15 = entry
    if (ctype ==1) {
        denom1 <- ifelse(cn[,4]==0, 1, cn[,3])
        denom2 <- ifelse(cn[,4]==0, 1, cn[,3]^2)
    } else {
        denom1 <- ifelse(cn[,4]==0, 1, cn[,11])
        denom2 <- ifelse(cn[,4]==0, 1, cn[,12])
    }

    hazard <- matrix(cn[,5] / denom1, ncol = fit$ntrans)
    varhaz <- matrix(cn[,5] / denom2, ncol = fit$ntrans)
    if (any(cmap != seq(along=cmap))) {
        hazard <- hazard[, cmap]
        varhaz <- varhaz[, cmap]
    }

    # Expand the result, one "hazard set" for each row of x2
    nx2 <- nrow(x2)
    h2 <- array(0, dim=c(nrow(hazard), nx2, ncol(hazard)))
    if (se.fit) v2 <- h2
    S <- matrix(0, nrow(hazard), nstate)
    S2 <- array(0, dim=c(nrow(hazard), nx2, nstate))
 
    H <- matrix(0, nstate, nstate)
    for (i in 1:nx2) {
        h2[,i,] <- apply(hazard %*% diag(risk2[i,]), 2, cumsum)
        if (se.fit) {
            d1 <- fit$xbar - rep(x[i,], each=nrow(fit$xbar))
            d2 <- apply(d1*hazard, 2, cumsum)
            d3 <- rowSums((d2%*% vmat) * d2)
#            v2[jj,] <- (apply(varhaz[jj,],2, cumsum) + d3) * (risk2[i])^2
        }

        S[1,] <- p0
        for (j in 2:ntime) {
            H[,] <- 0.0
            H[hfill] <- hazard[j,] *risk2[i,]
            if (stype==1) {
                diag(H) <- pmin(0, 1 + diag(H)- rowSums(H))
                S[j,] <- drop(S[j-1,] %*% H)  
            }
            else {
                diag(H) <- diag(H) - rowSums(H)
                S[j,] <- as.vector(S[j-1,] %*% expm(H))  # dgeMatrix issue
            }
        }

        S2[,i,] <- S
    }

    rval <- list(time=utime, xgrp=rep(1:nx2, each=nrow(hazard)),
                 pstate=S2, cumhaz=h2)
    if (se.fit) rval$varhaz <- v2
    rval
}
