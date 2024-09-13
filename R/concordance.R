# Automatically generated from the noweb directory
concordance <- function(object, ...) 
    UseMethod("concordance")

concordance.formula <- function(object, data,
                                weights, subset, na.action, cluster,
                                ymin, ymax, 
                                timewt=c("n", "S", "S/G", "n/G2", "I"),
                               influence=0, ranks=FALSE, reverse=FALSE,
                                timefix=TRUE, keepstrata=10, ...) {
    Call <- match.call()  # save a copy of of the call, as documentation
    timewt <- match.arg(timewt)
    if (missing(ymin)) ymin <- NULL
    if (missing(ymax)) ymax <- NULL
    
    index <- match(c("data", "weights", "subset", "na.action", 
                     "cluster"),
                   names(Call), nomatch=0)
    temp <- Call[c(1, index)]
    temp[[1L]] <-  quote(stats::model.frame)
    special <- c("strata", "cluster")
    temp$formula <- if(missing(data)) terms(object, special)
                    else              terms(object, special, data=data)
    mf <- eval(temp, parent.frame())  # model frame
    if (nrow(mf) ==0) stop("No (non-missing) observations")
    Terms <- terms(mf)

    Y <- model.response(mf)
    if (inherits(Y, "Surv")) {
        if (timefix) Y <- aeqSurv(Y)
        if (ncol(Y) == 3 && timewt %in% c("S/G", "n/G", "n/G2"))
            stop(gettext("'%s' timewt option not supported for (time1, time2) data", timewt))
    } else {
        if (is.factor(Y) && (is.ordered(Y) || length(levels(Y))==2))
            Y <- Surv(as.numeric(Y))
        else if (is.numeric(Y) && is.vector(Y))  Y <- Surv(Y)
        else if (is.logical(Y)) Y <- Surv(as.numeric(Y))
        else stop("left hand side of the formula must be a numeric vector,
 survival object, or an orderable factor")
        timewt <- "n"
        if (timefix) Y <- aeqSurv(Y)
    }
    n <- nrow(Y)
    
    wt <- model.weights(mf)
    offset<- attr(Terms, "offset")
    if (length(offset)>0) stop("Offset terms not allowed")

    stemp <- untangle.specials(Terms, "strata")
    if (length(stemp$vars)) {
        if (length(stemp$vars)==1) strat <- mf[[stemp$vars]]
        else strat <- strata(mf[,stemp$vars], shortlabel=TRUE)
        Terms <- Terms[-stemp$terms]
    }
    else strat <- NULL
    
    # if "cluster" was an argument, use it, otherwise grab it from the model
    group <- model.extract(mf, "cluster")
    cluster<- attr(Terms, "specials")$cluster
    if (length(cluster)) {
        tempc <- untangle.specials(Terms, 'cluster', 1:10)
        ord <- attr(Terms, 'order')[tempc$terms]
        if (any(ord>1)) stop("Cluster can not be used in an interaction")
        cluster <- strata(mf[,tempc$vars], shortlabel=TRUE)  #allow multiples
        Terms <- Terms[-tempc$terms]  # toss it away
    }
    if (length(group)) cluster <- group
                                            
    x <- model.matrix(Terms, mf)[,-1, drop=FALSE]  #remove the intercept

    if (!is.null(ymin) & (length(ymin)> 1 || !is.numeric(ymin)))
        stop("ymin must be a single number")
    if (!is.null(ymax) & (length(ymax)> 1 || !is.numeric(ymax)))
        stop("ymax must be a single number")
    if (!is.logical(reverse)) 
        stop("the reverse argument must be TRUE/FALSE")
 
    fit <- concordancefit(Y, x, strat, wt, ymin, ymax, timewt, cluster,
                           influence, ranks, reverse, keepstrata=keepstrata)
    na.action <- attr(mf, "na.action")
    if (length(na.action)) fit$na.action <- na.action
    fit$call <- Call

    class(fit) <- 'concordance'
    fit
}

print.concordance <- function(x, digits= max(1L, getOption("digits") - 3L), 
                              ...) {
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
        }
    omit <- x$na.action
    if(length(omit))
        cat("n=", x$n, " (", naprint(omit), ")\n", sep = "")
    else cat("n=", x$n, "\n")
    
    if (is.null(x$var)) {
        # Result of a call with std.err = FALSE
        cat("Concordance= ", format(x$concordance, digits=digits), "\n")
    } else {
        if (length(x$concordance) > 1) {
            # result of a call with multiple fits
            tmat <- cbind(concordance= x$concordance, se=sqrt(diag(x$var)))
            print(round(tmat, digits=digits), ...)
            cat("\n")
        }
        else cat("Concordance= ", format(x$concordance, digits=digits), 
                 " se= ", format(sqrt(x$var), digits=digits), '\n', sep='')
    }
    if (!is.matrix(x$count) || nrow(x$count < 11)) 
        print(round(x$count,2))
    invisible(x)
}

concordancefit <- function(y, x, strata, weights, ymin=NULL, ymax=NULL, 
                            timewt=c("n", "S", "S/G", "n/G2", "I"),
                            cluster, influence=0, ranks=FALSE, reverse=FALSE,
                            timefix=TRUE, keepstrata=10, std.err =TRUE) {
    # The coxph program may occassionally fail, and this will kill the C
    #  routine further below.  So check for it.
    if (any(is.na(x)) || any(is.na(y))) return(NULL)
    timewt <- match.arg(timewt)
    if (!is.null(ymin) && !is.numeric(ymin)) stop("ymin must be numeric")
    if (!is.null(ymax) && !is.numeric(ymax))    stop("ymax must be numeric")
    if (!std.err) {ranks <- FALSE; influence <- 0;}

    n <- length(y)
    X <- as.matrix(x)
    nvar <- ncol(X)
    if (nvar >1) {
        Xname <- colnames(X)
        if (is.null(Xname)) Xname <- paste0("X", 1:nvar)
    }
    if (nrow(X) != n) stop("x and y are not the same length")

    if (missing(strata) || length(strata)==0) strata <- rep(1L, n)
    if (length(strata) != n)
        stop("y and strata are not the same length")
    if (missing(weights) || length(weights)==0) weights <- rep(1.0, n)
    else {
        if (length(weights) != n) stop("y and weights are not the same length")
        storage.mode(weights) <- "double" # for the .Call, in case of integers
    }
    if (is.Surv(y)) {
        ny <- ncol(y)
        if (ny == 3 && timewt %in% c("S/G", "n/G2"))
            stop(gettext("'%s' timewt option not supported for (time1, time2) data", timewt))
        if (!is.null(attr(y, "states"))) 
            stop("concordance not defined for a multi-state outcome")
        if (timefix) y <- aeqSurv(y)
        if (!is.null(ymin)) {
            censored <- (y[,ny] ==0)
            # relaxed this rule, 30 March 2023
            #if (any(y[censored, ny-1] < ymin))
            #    stop("data has a censored value less than ymin")
            #else y[,ny-1] <- pmax(y[,ny-1], ymin)
            y[,ny-1] <- pmax(y[,ny-1], ymin)
        }
    } else {
        # should only occur if another package calls this routine
        if (is.factor(y) && (is.ordered(y) || length(levels(y))==2))
            y <- Surv(as.numeric(y))
        else if (is.numeric(y) && is.vector(y))  y <- Surv(y)
        else stop("left hand side of the formula must be a numeric vector,
 survival object, or an orderable factor")
    }

    type <- attr(y, "type")
    if (type %in% c("left", "interval"))
        stop("left or interval censored data is not supported")
    if (type %in% c("mright", "mcounting"))
        stop("multiple state survival is not supported")
    storage.mode(y) <- "double"   # in case of integers, for the .Call

    if (!is.null(ymin) && any(y[, ny-1] < ymin))
        y[,ny-1] <- pmax(y[,ny-1], ymin)
    # ymax is dealt with in the docount routine, as shifting end of a (t1, t2)
    #  interval could generate invalid data

    nstrat <- length(unique(strata))
    if (!is.logical(keepstrata)) {
        if (!is.numeric(keepstrata))
            stop("keepstrata argument must be logical or numeric")
        else keepstrata <- (nstrat <= keepstrata)
    }
    if (nvar>1 || nstrat ==1) keepstrata <- FALSE  #keeping both is difficult

    if (timewt %in% c("n", "I") && nstrat > 10 && !keepstrata) {
        # Special trickery for matched case-control data, where the
        #  number of strata is huge, n per strata is small, and compute
        #  time becomes excessive.  Make the data all one strata, but over
        #  disjoint time intervals
        stemp <- as.numeric(as.factor(strata)) -1
        if (ncol(y) ==3) {
            delta <- 2+ max(y[,2]) - min(y[,1])
            y[,1] <- y[,1] + stemp*delta
            y[,2] <- y[,2] + stemp*delta
        }
        else {
            delta <- max(y[,1]) +2
            m1 <- rep(-1L, nrow(y))
            y <- Surv(m1 + stemp*delta, y[,1] + stemp*delta, y[,2])
        }
        strata <- rep(1L, n)
        nstrat <- 1
    }

    # This routine is called once per stratum
    docount <- function(y, risk, wts, timeopt= 'n') {
        n <- length(risk)
        ny <- ncol(y)   # 2 or 3

        if (sum(y[,ncol(y)]) ==0) {
            # the special case of a stratum with no events (it happens)
            # No need to do any more work
            return(list(count= rep(0.0, 6), influence=matrix(0.0, n, 5),
                        resid=NULL))
        }
     
        # this next line is mostly invoked in stratified logistic, where
        #  only 1 event per stratum occurs.  All time weightings are the same
        # don't waste time even if the user asked for something different
        if (sum(y[,ny]) <2) timeopt <- 'n'
        
        # order the data: reverse time, censors before deaths
        if (ny ==2) {
            sort.stop <- order(-y[,1], y[,2], risk) -1L 
        } else {
            sort.stop  <- order(-y[,2], y[,3], risk) -1L   #order by endpoint
            sort.start <- order(-y[,1]) -1L       
        }

        if (timeopt == 'n') {
            deaths <- y[,ny] > 0
            etime  <- sort(unique(y[deaths, ny-1]))  # event times
        }
        else {
            if (ny==2) {
                sort.stop <- order(-y[,1], y[,2], risk) -1L 
                gfit <- .Call(Cfastkm1, y, wts, sort.stop)
            } else {
                sort.stop  <- order(-y[,2], y[,3], risk) -1L   #order by endpoint
                sort.start <- order(-y[,1]) -1L       
                gfit <- .Call(Cfastkm2, y, wts, sort.start, sort.stop)
            }
            etime <- gfit$etime
        }
         
        timewt <- switch(timeopt,
                         "S"   = sum(wts)* gfit$S/gfit$nrisk,
                         "S/G" = sum(wts)* gfit$S/ (gfit$G * gfit$nrisk),
                         "n" =   rep(1.0, length(etime)),
                         "n/G2"= 1/gfit$G^2,
                         "I"  =  1/gfit$nrisk)
        if (any(!is.finite(timewt))) stop("program error, notify author")
                    
        if (!is.null(ymax)) timewt[etime > ymax] <- 0
 
        # match each prediction score to the unique set of scores
        # (to deal with ties)
        utemp <- match(risk, sort(unique(risk)))
        bindex <- btree(max(utemp))[utemp]
        
        if (std.err) {
            if (ncol(y) ==2)
                fit <- .Call(Cconcordance3, y, bindex, wts, rev(timewt), 
                             sort.stop, ranks)
            else fit <- .Call(Cconcordance4, y, bindex, wts, rev(timewt), 
                              sort.start, sort.stop, ranks)
            if (ranks) {
                if (ncol(y)==2) dtime <- y[y[,2]==1, 1]
                else dtime <- y[y[,3]==1, 2]
                temp <- cbind(sort(dtime), fit$resid)
                colnames(temp) <- c("time", "rank", "timewt", "casewt")
                fit$resid <- temp[temp[,3] > 0,]  # don't return zeros
            }
        }
        else {
            if (ncol(y) ==2)
                fit <- .Call(Cconcordance5, y, bindex, wts, rev(timewt), 
                             sort.stop)
            else fit <- .Call(Cconcordance6, y, bindex, wts, rev(timewt), 
                              sort.start, sort.stop)
        }
        fit
    }
        
    # iterate over predictors (x) if necessary, calling docount for each
    #  then repack the returned list
    cfun <- function(y, x, weights, timewt) {
        if (!is.matrix(x) || ncol(x) ==1) {
            fit <- docount(y, as.vector(x), weights, timewt)
        } else {
            temp <- lapply(1:ncol(x), function(i) 
                docount(y, x[,i], weights, timewt))
            fit <- list(count = t(sapply(temp, function(x) x$count)))
            nvar <- ncol(X)
            if (std.err) 
                fit$influence <-array(unlist(lapply(temp, function(x) x$influence)),
                         dim=c(dim(temp[[1]]$influence), nvar))
            if (ranks) {
                fit$resid <- array(unlist(lapply(temp, function(x) x$resid)),
                                   dim=c(dim(temp[[1]]$resid), nvar),
                                   dimnames =c(dimnames(temp[[1]]$resid),
                                               list(Xname)))
            }
        }
        fit
    }       
        
    # unpack the strata, if needed
    if (nstrat < 2) fit <- cfun(y, drop(X), weights, timewt)
    else {
        # iterate over strata, calling cfun for each
        strata <- as.factor(strata)
        ustrat <- levels(strata)[table(strata) >0]  #some strata may have 0 obs
        tfit <- lapply(ustrat, function(i) {
            keep <- which(strata== i)
            cfun(y[keep,,drop=FALSE], X[keep,,drop=FALSE], weights[keep], timewt)
        })

        # note that both keepstrata and ncol(X) >1 won't both be true, so
        #  count will be a vector or matrix, never an array
        temp <- unlist(lapply(tfit, function(x) x$count))
        if (!keepstrata & nvar ==1)  # collapse over strata
            fit <- list(count= colSums(matrix(temp, ncol=6, byrow=TRUE)))
        else fit <- list(count = matrix(temp, ncol=6, byrow=TRUE))

        if (std.err || influence >0) {
            # the influence array is per observation, strata splits those rows
            #  up for computation, but not the result.  Each strata will be a
            #  different size.
            if (nvar ==1) {
                temp <- matrix(0, n, 5)
                for (i in 1:nstrat) 
                    temp[strata==ustrat[i],] <- tfit[[i]]$influence
            } else {
                temp <- array(0, dim=c(n, 5, nvar))
                for (i in 1:nstrat)
                    temp[strata==ustrat[i],,] <- tfit[[i]]$influence
            }
            fit$influence <- temp
        }
        if (ranks) {
            # same reassemby task for resid as for influence
            if (nvar ==1) {
                temp <- matrix(0, n, 4,
                               dimnames=list(NULL, c("time", "rank", "timewt",
                                                     "casewt")))
                for (i in 1:nstrat) 
                    temp[strata==ustrat[i],] <- tfit[[i]]$resid
            } else {        
                temp <- array(0, dim=c(n, 4, nvar),
                               dimnames=list(NULL, c("time", "rank", "timewt",
                                                     "casewt"), Xname))
                for (i in 1:nstrat)
                    temp[strata==ustrat[i],,] <- tfit[[i]]$resid
            }
            fit$resid <- temp
        }
    }
    
    # Assemble the result
    cname <- c("concordant", "discordant", "tied.x", "tied.y","tied.xy")
    if (nvar > 1) {  # concordance per outcome (count has one row per X col)
        npair <- rowSums(fit$count[,1:3])
        somer <- (fit$count[,1] - fit$count[,2])/ npair
        names(somer) <- Xname
        coxvar <- fit$count[,6]/(4*npair^2)
        rval <- list(concordance = (somer +1)/2, count=fit$count[,1:5], n=n)
        dimnames(rval$count) <- list(Xname, cname)
        if (std.err || influence ==1 || influence > 3) {
            # dfbeta will have one column per outcome, one row per obs
            #  influence has dimension (n, 5, nvar)
            dfbeta <- ((fit$influence[,1,]- fit$influence[,2,]) -
                       (apply(fit$influence[,1:3,], c(1,3), sum) *
                        rep(somer, each=n)))* weights/ (2* rep(npair, each= n))
            if (!missing(cluster) && length(cluster) > 0)
                dfbeta <- rowsum(dfbeta, cluster)
            cvar <- crossprod(dfbeta)
        }
    } 
    else {
        if (nstrat > 1) ctemp <- colSums(fit$count) else ctemp <- fit$count
        npair <- sum(ctemp[1:3])
        somer <- (ctemp[1] - ctemp[2])/ npair  # no concordance per stratum
        coxvar <- sum(ctemp[6])/(4*npair^2)  # cox variance
        if (keepstrata) {
            rval <- list(concordance = (somer+1)/2, count=fit$count[,1:5], n=n)
            dimnames(rval$count) <- list(levels(strata), cname)
        }
        else {
            rval <- list(concordance= (somer+1)/2, count=ctemp[1:5], n=n)
            names(rval$count) <- cname
        }

        if (std.err || influence ==1 || influence ==3) {
            # deriv of (A/B) = dA/B - dB*A/B^2 = (dA - db*A/B)) /B
            dfbeta <- ((fit$influence[,1]- fit$influence[,2]) - 
                       (rowSums(fit$influence[,1:3])*somer))*weights/(2*npair)
            # dfbeta is a vector of length n
            if (!missing(cluster) && length(cluster)>0) 
                dfbeta <- drop(rowsum(dfbeta, cluster))
            cvar <- sum(dfbeta^2)
        }
    }


    if (std.err) {
        rval$var <- cvar
        rval$cvar <- coxvar
    }
    if (influence == 1 || influence==3) rval$dfbeta <- dfbeta
    if (influence >=2) {
        rval$influence <- fit$influence
        if (nvar > 1) dimnames(rval$influence) <- list(NULL, cname, Xname)
        else dimnames(rval$influence) <- list(NULL, cname)
    }

    if (ranks) rval$ranks <- data.frame(fit$resid)

    if (reverse) {
        # flip concordant/discordant values but not the labels
        rval$concordance <- 1- rval$concordance
        if (!is.null(rval$dfbeta)) rval$dfbeta <- -rval$dfbeta
            
        if (is.matrix(rval$count)) {
            rval$count <- rval$count[, c(2,1,3,4,5)]
            colnames(rval$count) <- colnames(rval$count)[c(2,1,3,4,5)]
        }
        else {
            rval$count <- rval$count[c(2,1,3,4,5)]
            names(rval$count) <- names(rval$count)[c(2,1,3,4,5)]
        }

        if (!is.null(rval$influence)) {
            if (nvar >1) {
                rval$influence <- rval$influence[,c(2,1,3,4,5),]
                dimnames(rval$influence) <- list(NULL, cname, Xname)
            } else {
                rval$influence <- rval$influence[, c(2,1,3,4,5)]
                dimnames(rval$influence) <- list(NULL, cname)
            }
        }
        if (ranks) rval$ranks[,"rank"] <- -rval$ranks[,"rank"]
    }
    rval
}

btree <- function(n) {
   tfun <- function(n, id, power) {
       if (n==1L) id
       else if (n==2L) c(2L *id + 1L, id)
       else if (n==3L) c(2L*id + 1L, id, 2L*id +2L)
       else {
           nleft <- if (n== power*2L) power  else min(power-1L, n-power%/%2L)
           c(tfun(nleft, 2L *id + 1L, power%/%2), id,
             tfun(n-(nleft+1L), 2L*id +2L, power%/%2))
       }
   }
   tfun(as.integer(n), 0L, as.integer(2^(floor(logb(n-1,2)))))
}
cord.getdata <- function(object, newdata=NULL, cluster=NULL, need.wt, timefix=TRUE) {
    # For coxph object, don't reconstruct the model frame unless we must.
    # This will occur if weights, strata, or cluster are needed, or if
    #  there is a newdata argument.  Of course, if the model frame is 
    #  already present, then use it!
    Terms <- terms(object)
    specials <- attr(Terms, "specials")
    if (!is.null(specials$tt)) 
        stop("cannot yet handle models with tt terms")
 
    if (!is.null(newdata)) {
        mf <- model.frame(Terms, data=newdata)
        y <- model.response(mf)
        # why not model.frame(object, newdata)?  Because model.frame.coxph
        #  uses the Call to fill in names for subset, if present, which of
        #  course is not relevant to the new data. model.frame.lm does the
        #  same, btw.
        y <- model.response(mf)
        if (!is.Surv(y)) {
            if (is.numeric(y) && is.vector(y))  y <- Surv(y)
            else stop("left hand side of the formula  must be a numeric vector or a survival object")
        }
        if (timefix) y <- aeqSurv(y)

        # special handling for coxph objects: object to tt()
        specials <- attr(Terms, "specials")
        if (!is.null(specials$tt)) stop("newdata not allowed for tt() terms")

            
        yhat <- predict(object, newdata=newdata, na.action = na.omit) 
        # yes, the above will construct mf again, but all the special processing
        #  we get for lm, glm, coxph  is worth it.
        rval <- list(y= y, x= as.vector(yhat))

        # recreate the strata, if needed
        if (length(attr(Terms, "specials")$strata) > 0) {
            stemp <- untangle.specials(Terms, 'strata', 1)
            if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
            else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
            rval$strata <- as.integer(strata.keep)
        }
    } 
    else {
        mf <- object$model
        y <- object$y
        if (is.null(y)) {
            if (is.null(mf)) mf <- model.frame(object)
            y <- model.response(mf)
        }
        x <- object$linear.predictors    # used by most
        if (is.null(x)) x <- object$fitted.values # used by lm
        if (is.null(x)) {object$na.action <- NULL; x <- predict(object)}
        rval <- list(y = y, x= x)

        if (!is.null(specials$strata)) {
            if (is.null(mf)) mf <- model.frame(object)
            stemp <- untangle.specials(Terms, 'strata', 1)
            if (length(stemp$vars)==1) rval$strata <- mf[[stemp$vars]]
            else rval$strata <- strata(mf[,stemp$vars], shortlabel=TRUE)
        } 
    }
        
    if (need.wt) {
        if (is.null(mf)) mf <- model.frame(object)
        rval$weights <- model.weights(mf)
    }

    if (is.null(cluster)) {
        if (!is.null(specials$cluster)) {
            if (is.null(mf)) mf <- model.frame(object)
            tempc <- untangle.specials(Terms, 'cluster', 1:10)
            ord <- attr(Terms, 'order')[tempc$terms]
            rval$cluster <- strata(mf[,tempc$vars], shortlabel=TRUE) 
        }
        else if (!is.null(object$call$cluster)) {
            if (is.null(mf)) mf <- model.frame(object)
            rval$cluster <- model.extract(mf, "cluster")
        }
    }
    else rval$cluster <- cluster
    rval
}
concordance.lm <- function(object, ..., newdata, cluster, ymin, ymax, 
                           influence=0, ranks=FALSE, timefix=TRUE,
                           keepstrata=10) {
    Call <- match.call()
    fits <- list(object, ...)
    nfit <- length(fits)
    fname <- as.character(Call)  # like deparse(substitute()) but works for ...
    fname <- fname[1 + 1:nfit]
    notok <- sapply(fits, function(x) !inherits(x, "lm"))
    if (any(notok)) {
        # a common error is to mistype an arg, "ramk=TRUE" for instance,
        #  and it ends up in the ... list
        # try for a nice message in this case: the name of the arg if it
        #  has one other than "object", fname otherwise
        indx <- which(notok)
        id2 <- names(Call)[indx+1]
        temp <- ifelse(id2 %in% c("","object"), fname, id2)
        stop(temp, " argument is not an appropriate fit object")
    }
        
    cargs <- c("ymin", "ymax","influence", "ranks", "keepstrata")
    cfun <- Call[c(1, match(cargs, names(Call), nomatch=0))]
    cfun[[1]] <- cord.work   # or quote(survival:::cord.work)
    cfun$fname <- fname
    
    if (missing(newdata)) newdata <- NULL
    if (missing(cluster)) cluster <- NULL
    need.wt <- any(sapply(fits, function(x) !is.null(x$call$weights)))
    
    cfun$data <- lapply(fits, cord.getdata, newdata=newdata, cluster=cluster,
                        need.wt=need.wt, timefix=timefix)
    rval <- eval(cfun, parent.frame())
    rval$call <- Call
    rval
}

concordance.survreg <- function(object, ..., newdata, cluster, ymin, ymax,
                                timewt=c("n", "S", "S/G", "n/G2", "I"),
                                influence=0, ranks=FALSE, timefix= TRUE,
                                keepstrata=10) {
    Call <- match.call()
    fits <- list(object, ...)
    nfit <- length(fits)
    fname <- as.character(Call)  # like deparse(substitute()) but works for ...
    fname <- fname[1 + 1:nfit]
    notok <- sapply(fits, function(x) !inherits(x, "survreg"))
    if (any(notok)) {
        # a common error is to mistype an arg, "ramk=TRUE" for instance,
        #  and it ends up in the ... list
        # try for a nice message in this case: the name of the arg if it
        #  has one other than "object", fname otherwise
        indx <- which(notok)
        id2 <- names(Call)[indx+1]
        temp <- ifelse(id2 %in% c("","object"), fname, id2)
        stop(temp, " argument is not an appropriate fit object")
    }
        
    cargs <- c("ymin", "ymax","influence", "ranks", "timewt", "keepstrata")
    cfun <- Call[c(1, match(cargs, names(Call), nomatch=0))]
    cfun[[1]] <- cord.work
    cfun$fname <- fname
    
    if (missing(newdata)) newdata <- NULL
    if (missing(cluster)) cluster <- NULL
    need.wt <- any(sapply(fits, function(x) !is.null(x$call$weights)))
    
    cfun$data <- lapply(fits, cord.getdata, newdata=newdata, cluster=cluster,
                        need.wt=need.wt, timefix=timefix)
    rval <- eval(cfun, parent.frame())
    rval$call <- Call
    rval
}
    
concordance.coxph <- function(object, ..., newdata, cluster, ymin, ymax, 
                               timewt=c("n", "S", "S/G", "n/G2", "I"),
                               influence=0, ranks=FALSE, timefix=TRUE,
                               keepstrata=10) {
    Call <- match.call()
    fits <- list(object, ...)
    nfit <- length(fits)
    fname <- as.character(Call)  # like deparse(substitute()) but works for ...
    fname <- fname[1 + 1:nfit]
    notok <- sapply(fits, function(x) !inherits(x, "coxph"))
    if (any(notok)) {
        # a common error is to mistype an arg, "ramk=TRUE" for instance,
        #  and it ends up in the ... list
        # try for a nice message in this case: the name of the arg if it
        #  has one other than "object", fname otherwise
        indx <- which(notok)
        id2 <- names(Call)[indx+1]
        temp <- ifelse(id2 %in% c("","object"), fname, id2)
        stop(temp, " argument is not an appropriate fit object")
    }
        
    # the cargs trick is a nice one, but it only copies over arguments that
    #  are present.  If 'ranks' was not specified, the default of FALSE is
    #  not set.  We keep it in the arg list only to match the documentation.
    cargs <- c("ymin", "ymax","influence", "ranks", "timewt", "keepstrata")
    cfun <- Call[c(1, match(cargs, names(Call), nomatch=0))]
    cfun[[1]] <- cord.work   # a copy of the function
    cfun$fname <- fname
    cfun$reverse <- TRUE

    if (missing(newdata)) newdata <- NULL
    if (missing(cluster)) cluster <- NULL
    need.wt <- any(sapply(fits, function(x) !is.null(x$call$weights)))
    
    cfun$data <- lapply(fits, cord.getdata, newdata=newdata, cluster=cluster,
                        need.wt=need.wt, timefix=timefix)
    rval <- eval(cfun, parent.frame())
    rval$call <- Call
    rval
}
cord.work <- function(data, timewt, ymin, ymax, influence=0, ranks=FALSE, 
                      reverse, fname, keepstrata) {
    Call <- match.call()
    fargs <- c("timewt", "ymin", "ymax", "influence", "ranks", "reverse",
               "keepstrata")
    fcall <- Call[c(1, match(fargs, names(Call), nomatch=0))]
    fcall[[1L]] <- concordancefit

    nfit <- length(data)
    if (nfit==1) {
        dd <- data[[1]] 
        fcall$y <- dd$y
        fcall$x <- dd$x
        fcall$strata <- dd$strata
        fcall$weights <- dd$weights
        fcall$cluster  <- dd$cluster
        rval <- eval(fcall, parent.frame())
    }
    else {
        # Check that all of the models used the same data set, in the same
        #  order, to the best of our abilities
        n <- length(data[[1]]$x)
        for (i in 2:nfit) {
            if (length(data[[i]]$x) != n)
                stop("all models must have the same sample size")
            
            if (!identical(data[[1]]$y, data[[i]]$y))
                warning("models do not have the same response vector")
            
            if (!identical(data[[1]]$weights, data[[i]]$weights))
                stop("all models must have the same weight vector")
        }
        
        if (influence==2) fcall$influence <-3 else fcall$influence <- 1
        flist <- lapply(data, function(d) {
                         temp <- fcall
                         temp$y <- d$y
                         temp$x <- d$x
                         temp$strata <- d$strata
                         temp$weights <- d$weights
                         temp$cluster <- d$cluster
                         eval(temp, parent.frame())
                     })
            
        for (i in 2:nfit) {
            if (length(flist[[1]]$dfbeta) != length(flist[[i]]$dfbeta))
                stop("models must have identical clustering")
        }
        count = do.call(rbind, lapply(flist, function(x) {
            if (is.matrix(x$count)) colSums(x$count) else x$count}))

        concordance <- sapply(flist, function(x) x$concordance)
        dfbeta <- sapply(flist, function(x) x$dfbeta)

        names(concordance) <- fname
        rownames(count) <- fname

        wt <- data[[1]]$weights
        if (is.null(wt)) vmat <- crossprod(dfbeta)
        else vmat <- t(wt * dfbeta) %*% dfbeta
        rval <- list(concordance=concordance, count=count, 
                     n=flist[[1]]$n, var=vmat,
                     cvar= sapply(flist, function(x) x$cvar))

        if (influence==1) rval$dfbeta <- dfbeta
        else if (influence ==2) {
            temp <- unlist(lapply(flist, function(x) x$influence))
            rval$influence <- array(temp, 
                                    dim=c(dim(flist[[1]]$influence), nfit))
        }
        
        if (ranks) {
            temp <- lapply(flist, function(x) x$ranks)
            rdat <- data.frame(fit= rep(fname, sapply(temp, nrow)),
                               do.call(rbind, temp))
            row.names(rdat) <- NULL
            rval$ranks <- rdat
        }
     }
    
    class(rval) <- "concordance"
    rval
}
coef.concordance <- function(object, ...) object$concordance
vcov.concordance <- function(object, ...) object$var
