# Automatically generated from the noweb directory
concordance <- function(object, ...) 
    UseMethod("concordance")

concordance.formula <- function(object, data,
                                weights, subset, na.action, cluster,
                                ymin, ymax, 
                                timewt=c("n", "S", "S/G", "n/G", "n/G2", "I"),
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
    } else {
        if (is.factor(Y) && (is.ordered(Y) || length(levels(Y))==2))
            Y <- Surv(as.numeric(Y))
        else if (is.numeric(Y) && is.vector(Y))  Y <- Surv(Y)
        else stop("left hand side of the formula must be a numeric vector,
 survival object, or an orderable factor")
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
        if (any(ord>1)) stop ("Cluster can not be used in an interaction")
        cluster <- strata(mf[,tempc$vars], shortlabel=TRUE)  #allow multiples
        Terms <- Terms[-tempc$terms]  # toss it away
    }
    if (length(group)) cluster <- group
                                            
    x <- model.matrix(Terms, mf)[,-1, drop=FALSE]  #remove the intercept
    if (ncol(x) > 1) stop("Only one predictor variable allowed")

    if (!is.null(ymin) & (length(ymin)> 1 || !is.numeric(ymin)))
        stop("ymin must be a single number")
    if (!is.null(ymax) & (length(ymax)> 1 || !is.numeric(ymax)))
        stop("ymax must be a single number")
    if (!is.logical(reverse)) 
        stop ("the reverse argument must be TRUE/FALSE")
 
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
    
    if (length(x$concordance) > 1) {
        # result of a call with multiple fits
        tmat <- cbind(concordance= x$concordance, se=sqrt(diag(x$var)))
        print(round(tmat, digits=digits), ...)
        cat("\n")
    }
    else cat("Concordance= ", format(x$concordance, digits=digits), " se= ", 
             format(sqrt(x$var), digits=digits), '\n', sep='')

    if (!is.matrix(x$count) || nrow(x$count < 11)) 
        print(round(x$count,2))
    invisible(x)
    }

concordancefit <- function(y, x, strata, weights, ymin=NULL, ymax=NULL, 
                            timewt=c("n", "S", "S/G", "n/G", "n/G2", "I"),
                            cluster, influence=0, ranks=FALSE, reverse=FALSE,
                            timefix=TRUE, keepstrata=10, robustse =TRUE) {
    # The coxph program may occassionally fail, and this will kill the C
    #  routine further below.  So check for it.
    if (any(is.na(x)) || any(is.na(y))) return(NULL)
    timewt <- match.arg(timewt)

    if (!robustse) {ranks <- FALSE; influence =0;}

    # these should only occur if something other package calls this routine
    if (!is.Surv(y)) {
        if (is.factor(y) && (is.ordered(y) || length(levels(y))==2))
            y <- Surv(as.numeric(y))
        else if (is.numeric(y) && is.vector(y))  y <- Surv(y)
        else stop("left hand side of the formula must be a numeric vector,
 survival object, or an orderable factor")
    }
     
    # When concordance is called with new data, this is the first point where
    #  aeqSurv will be encountered.  Survival data or not, near ties will be
    #  fatal.  
    if (timefix) y <- aeqSurv(y) 
    n <- length(y)
    if (length(x) != n) stop("x and y are not the same length")
    if (missing(strata) || length(strata)==0) strata <- rep(1L, n)
    if (length(strata) != n)
        stop("y and strata are not the same length")
    if (missing(weights) || length(weights)==0) weights <- rep(1.0, n)
    else if (length(weights) != n) stop("y and weights are not the same length")

    type <- attr(y, "type")
    if (type %in% c("left", "interval"))
        stop("left or interval censored data is not supported")
    if (type %in% c("mright", "mcounting"))
        stop("multiple state survival is not supported")

    nstrat <- length(unique(strata))
    if (!is.logical(keepstrata)) {
        if (!is.numeric(keepstrata))
            stop("keepstrat argument must be logical or numeric")
        else keepstrata <- (nstrat <= keepstrata)
    }

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
    docount <- function(y, risk, wts, timeopt= 'n', timefix) {
        n <- length(risk)
        # this next line is mostly invoked in stratified logistic, where
        #  only 1 event per stratum occurs.  All time weightings are the same
        # don't waste time even if the user asked for something different
        if (sum(y[,ncol(y)]) <2) timeopt <- 'n'
        
        sfit <- survfit(y~1, weights=wts, se.fit=FALSE, timefix=timefix)
        etime <- sfit$time[sfit$n.event > 0]
        esurv <- sfit$surv[sfit$n.event > 0]
        
        if (length(etime)==0) {
            # the special case of a stratum with no events (it happens)
            # No need to do any more work
            return(list(count= rep(0.0, 6), influence=matrix(0.0, n, 5),
                        resid=NULL))
        }

       if (timeopt %in% c("S/G", "n/G", "n/G2")) {
            temp <- y
            temp[,ncol(temp)] <- 1- temp[,ncol(temp)] # switch event/censor
            gfit <- survfit(temp~1, weights=wts, se.fit=FALSE, timefix=timefix)
            # G has the exact same time values as S
            gsurv <- c(1, gfit$surv)  # We want G(t-)
            gsurv <- gsurv[which(sfit$n.event > 0)]
        }

        npair <- (sfit$n.risk- sfit$n.event)[sfit$n.event>0]
        temp  <- ifelse(esurv==0, 0, esurv/npair)  # avoid 0/0
        timewt <- switch(timeopt,
                         "S" =  sum(wts)*temp,
                         "S/G" = sum(wts)* temp/ gsurv,
                         "n" =   rep(1.0, length(npair)),
                         "n/G" = 1/gsurv,
                         "n/G2"= 1/gsurv^2,
                         "I"  =  rep(1.0, length(esurv))
                     )
        if (!is.null(ymin)) timewt[etime < ymin] <- 0
        if (!is.null(ymax)) timewt[etime > ymax] <- 0
        timewt <- ifelse(is.finite(timewt), timewt, 0)  # 0 at risk case

        # order the data: reverse time, censors before deaths
        if (ncol(y)==2) { 
            sort.stop <- order(-y[,1], y[,2], risk) -1L 
        } else {
            sort.stop  <- order(-y[,2], y[,3], risk) -1L   #order by endpoint
            sort.start <- order(-y[,1]) -1L       
        }
 
        # match each prediction score to the unique set of scores
        # (to deal with ties)
        utemp <- match(risk, sort(unique(risk)))
        bindex <- btree(max(utemp))[utemp]
        
        storage.mode(y) <- "double"  # just in case y is integer
        storage.mode(wts) <- "double"
        if (robustse) {
            if (ncol(y) ==2)
                fit <- .Call(Cconcordance3, y, bindex, wts, rev(timewt), 
                             sort.stop, ranks)
            else fit <- .Call(Cconcordance4, y, bindex, wts, rev(timewt), 
                              sort.start, sort.stop, ranks)
            
            # The C routine gives back an influence matrix which has columns for
            #  concordant, discordant, tied on x but not y, tied on y, and tied
            #  on both x and y. 
            dimnames(fit$influence) <- list(NULL, 
                   c("concordant", "discordant", "tied.x", "tied.y", "tied.xy"))
            if (ranks) {
                if (ncol(y)==2) dtime <- y[y[,2]==1, 1]
                else dtime <- y[y[,3]==1, 2]
                temp <- data.frame(time= sort(dtime), fit$resid)
                names(temp) <- c("time", "rank", "timewt", "casewt", "variance")
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
        
    if (nstrat < 2) {
        fit <- docount(y, x, weights, timewt, timefix=timefix)
        count2 <- fit$count[1:5]
        vcox <- fit$count[6]
        fit$count <- fit$count[1:5]
        if (robustse) imat <- fit$influence
        if (ranks) resid <- fit$resid
    } else {
        strata <- as.factor(strata)
        ustrat <- levels(strata)[table(strata) >0]  #some strata may have 0 obs
        tfit <- lapply(ustrat, function(i) {
            keep <- which(strata== i)
            docount(y[keep,,drop=F], x[keep], weights[keep], timewt,
                    timefix=timefix)
        })
        temp <-  t(sapply(tfit, function(x) x$count))
        fit <- list(count = temp[,1:5])
        count2 <- colSums(fit$count)
        if (!keepstrata) fit$count <- count2
        vcox <- sum(temp[,6])
        if (robustse) {
            imat <- do.call("rbind", lapply(tfit, function(x) x$influence))
            # put it back into data order
            index <- match(1:n, (1:n)[order(strata)])
            imat <- imat[index,]
            if (ranks) {
                nr <- lapply(tfit, function(x) nrow(x$resid))
                resid <- do.call("rbind", lapply(tfit, function(x) x$resid))
                resid$strata <- rep(ustrat, nr)
            }
        }
    }
        
    npair <- sum(count2[1:3])
    if (!keepstrata && is.matrix(fit$count)) fit$count <- colSums(fit$count)
    somer <- (count2[1] - count2[2])/npair
    if (robustse) {
        dfbeta <- weights*((imat[,1]- imat[,2])/npair -
                           (somer/npair)* rowSums(imat[,1:3]))
        if (!missing(cluster) && length(cluster)>0) {
            dfbeta <- tapply(dfbeta, cluster, sum)
            dfbeta <- ifelse(is.na(dfbeta),0, dfbeta)  # if cluster is a factor
        }
        var.somer <- sum(dfbeta^2)
        rval <- list(concordance = (somer+1)/2, count=fit$count, n=n,
                     var = var.somer/4, cvar=vcox/(4*npair^2))
        }
    else  rval <- list(concordance = (somer+1)/2, count=fit$count, n=n,
                     cvar=vcox/(4*npair^2))
    if (is.matrix(rval$count))
        colnames(rval$count) <- c("concordant", "discordant", "tied.x", 
                                   "tied.y", "tied.xy")
    else names(rval$count) <- c("concordant", "discordant", "tied.x", "tied.y",
                           "tied.xy")

    if (influence == 1 || influence==3) rval$dfbeta <- dfbeta/2
    if (influence >=2) rval$influence <- imat
         
    if (ranks) rval$ranks <- resid
    if (reverse) {
        # flip concordant/discordant values but not the labels
        rval$concordance <- 1- rval$concordance
        if (!is.null(rval$dfbeta)) rval$dfbeta <- -rval$dfbeta
        if (!is.null(rval$influence)) {
            rval$influence <- rval$influence[,c(2,1,3,4,5)]
            colnames(rval$influence) <- colnames(rval$influence)[c(2,1,3,4,5)]
        }
        if (is.matrix(rval$count)) {
            rval$count <- rval$count[, c(2,1,3,4,5)]
            colnames(rval$count) <- colnames(rval$count)[c(2,1,3,4,5)]
        }
        else {
            rval$count <- rval$count[c(2,1,3,4,5)]
            names(rval$count) <- names(rval$count)[c(2,1,3,4,5)]
        }
        if (ranks) rval$ranks$rank <- -rval$ranks$rank
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
<<<<<<< working copy
        # do not use model.frame(object, newdata) --- if there are other
        # terms such as offset, weights, id,... it doesn't do what you
        # would expect.
=======
        y <- model.response(mf)
        # why not model.frame(object, newdata)?  Because model.frame.coxph
        #  uses the Call to fill in names for subset, if present, which of
        #  course is not relevant to the new data. model.frame.lm does the
        #  same, btw.
>>>>>>> merge rev
        y <- model.response(mf)
        if (!is.Surv(y)) {
            if (is.numeric(y) && is.vector(y))  y <- Surv(y)
            else stop("left hand side of the formula  must be a numeric vector or a survival object")
        }
        if (timefix) y <- aeqSurv(y)
<<<<<<< working copy
        xhat <- model.matrix(object,newdata)%*% coef(object)
        rval <- list(y= y, x= xhat)
        # the type of prediction does not matter, as long as it is a 
        #  monotone transform of the linear predictor
=======

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
>>>>>>> merge rev
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
                                timewt=c("n", "S", "S/G", "n/G", "n/G2", "I"),
<<<<<<< working copy
                                influence=0, ranks=FALSE, timefix=TRUE,
=======
                                influence=0, ranks=FALSE, timefix= TRUE,
>>>>>>> merge rev
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
                               timewt=c("n", "S", "S/G", "n/G", "n/G2", "I"),
<<<<<<< working copy
                               influence=0, ranks=FALSE, timefix=TRUE,
=======
                               influence=0, ranks=FALSE, 
                               timefix=!missing(newdata),
>>>>>>> merge rev
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
