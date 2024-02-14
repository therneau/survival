# Automatically generated from the noweb directory
# residuals for a survfit object
residuals.survfit <- function(object, times, 
                              type= "pstate",
                              collapse, weighted=FALSE, method=1, ...){

    if (!inherits(object, "survfit"))
        stop("argument must be a survfit object")
    if (object$type=="interval") {
        # trial code to support it
        # reconstruct the data set
        # create dummy time/status for all interval or left censored
        #   over the span of jump points in S, non-censored obs with
        #   weights proportional to the jumps
        # combine dummy + (exact, right) from original, compute KM
        # get pseudo for this new KM
        # collapse dummy obs back to a single
        stop("residuals for interval-censored data are not available")
        }

    if (missing(times)) stop("the times argument is required")
    # allow a set of alias
    temp <- c("pstate", "cumhaz", "sojourn", "survival",
                              "chaz", "rmst", "rmts", "auc")
    type <- match.arg(casefold(type), temp)
    itemp <-  c(1,2,3,1,2,3,3,3)[match(type, temp)]
    type <- c("pstate", "cumhaz", "auc")[itemp]

    if (missing(collapse)) 
         fit <- survresid.fit(object, times, type, weighted=weighted, 
                              method= method)
    else fit <- survresid.fit(object, times, type, collapse= collapse, 
                              weighted= weighted, method= method)

    fit$residuals
}

survresid.fit <- function(object, times, 
                              type= "pstate",
                              collapse, weighted=FALSE, method=1) {
    if (object$type=="interval") stop("interval censored not yet supported")
    survfitms <- inherits(object, "survfitms")
    coxsurv <- inherits(object, "survfitcox")  # should never be true, as there
                                               #  is a residuals.survfitcox
    timefix <- (is.null(object$timefix) || object$timefix)
    
    start.time <- object$start.time
    if (is.null(start.time)) start.time <- min(c(0, object$time))

    # check input arguments
    if (missing(times)) 
        stop ("the times argument is required")
    else {
        if (!is.numeric(times)) stop("times must be a numeric vector")
        times <- sort(unique(times))
        if (timefix) times <- aeqSurv(Surv(times))[,1]
    }

    # get the data
    Call <- object$call
    Terms <- object$terms

    # remember the name of the id variable, if present.
    #  but we don't try to parse it:  id= mydata$clinic becomes NULL
    idname <- Call$id
    if (is.name(idname)) idname <- as.character(idname)
    else idname <- NULL   
    # I always need the model frame
    mf <- model.frame(object)
    if (is.null(object$y)) Y <- model.response(mf)
    else Y <- object$y

    formula <- formula(object)
    # the chunk below is shared with survfit.formula 
    na.action <- getOption("na.action")
    if (is.character(na.action))
        na.action <- get(na.action)  # a hack to allow the shared code
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
    # end of shared code 

    xlev <- levels(X)

    # Deal with ties
    if (is.null(Call$timefix) || Call$timefix) newY <- aeqSurv(Y) else newY <- Y

    if (missing(collapse)) collapse <- (!(is.null(id)) && any(duplicated(id)))
    if (collapse && is.null(id)) stop("collapse argument requires an id or cluster argument in the survfit call")

    ny <- ncol(newY)
    if (collapse && any(X != X[1])) {
        # If the same id shows up in multiple curves, we just can't deal
        #  with it.
        temp <- unlist(lapply(split(id, X), unique))
        if (any(duplicated(temp)))
            stop("same id appears in multiple curves, cannot collapse")
    }
    
    timelab <- signif(times, 3)  # used for dimnames
    # What type of survival curve?
    stype <- Call$stype
    if (is.null(stype)) stype <- 1
    ctype <- Call$ctype
    if (is.null(ctype)) ctype <- 1
    if (!survfitms) {
        resid <- rsurvpart1(newY, X, casewt, times,
                            type, stype, ctype, object)
        if (collapse) {
            resid <- rowsum(resid, id, reorder=FALSE)
            dimnames(resid) <- list(id= unique(id), times=timelab)
            curve <- (as.integer(X))[!duplicated(id)] #which curve for each
        } 
        else {
            if (length(id) >0) dimnames(resid) <- list(id=id, times=timelab)
            curve <- as.integer(X)
        }
    }
    else {  # multi-state
        if (!collapse) {
            if (length(id >0)) d1name <- id else d1name <- NULL
            cluster <- d1name
            curve <- as.integer(X)
        }       
        else {
            d1name <- unique(id)
            cluster <- match(id, d1name)
            curve <- (as.integer(X))[!duplicated(id)]
        }
        resid <- rsurvpart2(newY, X, casewt, istate, times, cluster,
                            type, object, method=method, collapse=collapse)

        if (type == "cumhaz") {
            ntemp <- colnames(object$cumhaz)
            if (length(dim(resid)) ==3)
                 dimnames(resid) <- list(id=d1name, times=timelab, 
                                         cumhaz= ntemp)
            else dimnames(resid) <- list(id=d1name, cumhaz=ntemp)
        }
        else {
            ntemp <- object$states
            if (length(dim(resid)) ==3) 
                dimnames(resid) <- list(id=d1name, times=timelab, 
                                        state= ntemp)
            else dimnames(resid) <- list(id=d1name, state= ntemp)
        }
    }

    if (weighted && any(casewt !=1)) resid <- resid*casewt

    list(residuals= resid, curve= curve, id= id, idname=idname)
}
rsurvpart1 <- function(Y, X, casewt, times,
         type, stype, ctype, fit) {
     
    ntime <- length(times)
    etime <- (fit$n.event >0)
    ny <- ncol(Y)
    event <- (Y[,ny] >0)
    status <- Y[,ny]

    # 
    #  Create a list whose first element contains the location of
    #   the death times in curve 1, second element the death times for curve 2,
    #  
    if (is.null(fit$strata)) {
        fitrow <- list(which(etime))
    }
    else {
        temp1 <- cumsum(fit$strata)
        temp2 <- c(1, temp1+1)
        fitrow <- lapply(1:length(fit$strata), function(i) {
            indx <- seq(temp2[i], temp1[i])
            indx[etime[indx]] # keep the death times
        })
    }
    ff <- unlist(fitrow) 
 
    # for each time x, the index of the last death time which is <=x.
    #  0 if x is before the first death time in the fit object.
    #  The result is an index to the survival curve
    matchfun <- function(x, fit, index) {
        dtime <- fit$time[index]  # subset to this curve
        i2 <- findInterval(x, dtime, left.open=FALSE)
        c(0, index)[i2 +1]
    }
     
    # output matrix D will have one row per observation, one col for each
    #  reporting time. tindex and yindex have the same dimension as D.
    # tindex points to the last death time in fit which
    #  is <= the reporting time.  (If there is only 1 curve, each col of
    #  tindex will be a repeat of the same value.)
    tindex <- matrix(0L, nrow(Y), length(times))
    for (i in 1:length(fitrow)) {
        yrow <- which(as.integer(X) ==i)
        temp <- matchfun(times, fit, fitrow[[i]])
        tindex[yrow, ] <- rep(temp, each= length(yrow))
    }
    tindex[,] <- match(tindex, c(0,ff)) -1L  # the [,] preserves dimensions

    # repeat the indexing for Y onto fit$time.  Each row of yindex points
    #  to the last row of fit with death time <= Y[,ny]
    ny <- ncol(Y)
    yindex <- matrix(0L, nrow(Y), length(times))
    event <- (Y[,ny] >0)
    if (ny==3) startindex <- yindex
    for (i in 1:length(fitrow)) {
        yrow <- (as.integer(X) ==i)  # rows of Y for this curve
        temp <- matchfun(Y[yrow,ny-1], fit, fitrow[[i]])
        yindex[yrow,] <- rep(temp, ncol(yindex))
        if (ny==3) {
            temp <- matchfun(Y[yrow,1], fit, fitrow[[i]])
            startindex[yrow,] <- rep(temp, ncol(yindex))
        }
    }                    
    yindex[,] <- match(yindex, c(0,ff)) -1L
    if (ny==3) {
        startindex[,] <- match(startindex, c(0,ff)) -1L
        # no subtractions for report times before subject's entry
        startindex <- pmin(startindex, tindex) 
    }
    
    # Now do the work
    if (type=="cumhaz" || stype==2) {  # result based on hazards
        if (ctype==1) {
            death <- (yindex <= tindex & rep(event, ntime)) # an event occured at <= t

            term1 <- 1/fit$n.risk[ff]
            term2 <- lapply(fitrow, function(i) fit$n.event[i]/fit$n.risk[i]^2)
            term3 <- unlist(lapply(term2, cumsum))

            sum1 <- c(0, term1)[ifelse(death, 1+yindex, 1)]
            sum2 <- c(0, term3)[1 + pmin(yindex, tindex)]
            if (ny==3) sum3 <- c(0, term3)[1 + pmin(startindex, tindex)]

            if (ny==2) D <- matrix(sum1 -  sum2, ncol=ntime)
            else       D <- matrix(sum1 + sum3 - sum2, ncol=ntime)

            # survival is exp(-H) so the derivative is a simple transform of D
            if (type== "pstate") D <- -D* c(1,fit$surv[ff])[1+ tindex]
            else if (type == "auc") {
                auc1 <- lapply(fitrow, function(i) {
                             if (length(i) <=1) 0
                             else c(0, cumsum(diff(fit$time[i]) * (fit$surv[i])[-length(i)]))
                                 })  # AUC at each event time
                auc2 <- lapply(fitrow, function(i) {
                             if (length(i) <=1) 0
                             else {
                                 xx <- sort(unique(c(fit$time[i], times))) # all the times
                                 yy <- (fit$surv[i])[findInterval(xx, fit$time[i])]
                                 auc <- cumsum(c(diff(xx),0) * yy)
                                 c(0, auc)[match(times, xx)]
                                 }})  # AUC at the output times

                # Most often this function is called with a single curve, so make that case
                #  faster.  (Or I presume so: mapply and do.call may be more efficient than 
                #  I think for lists of length 1).
                if (length(fitrow)==1) { # simple case, most common to ask for auc 
                    wtmat <- pmin(outer(auc1[[1]], -auc2[[1]], '+'),0)
                    term1 <- term1 * wtmat
                    term2 <- unlist(term2) * wtmat
                    term3 <- apply(term2, 2, cumsum)
                }
                else { #more than one curve, compute weighted cumsum per curve
                    wtmat <- mapply(function(x, y) pmin(outer(x, -y, "+"), 0), auc1, auc2)
                    term1 <- term1 * do.call(rbind, wtmat)
                    temp <- mapply(function(x, y) apply(x*y, 2, cumsum), term2, wtmat)
                    term3 <- do.call(rbind, temp)
                }

                sum1 <- sum2 <- matrix(0, nrow(yindex), ntime)
                if (ny ==3) sum3 <- sum1
                for (i in 1:ntime) {
                    sum1[,i] <- c(0, term1[,i])[ifelse(death[,i], 1 + yindex[,i], 1)]
                    sum2[,i] <- c(0, term3[,i])[1 + pmin(yindex[,i], tindex[,i])]
                    if (ny==3) sum3[,i] <- c(0, term3[,i])[1 + pmin(startindex[,i], tindex[,i])]
                }
                # Perhaps a bit faster(?), but harder to read. And for AUC people usually only
                #  ask for one time point
                #sum1 <- rbind(0, term1)[cbind(c(ifelse(death, 1+yindex, 1)), c(col(yindex)))]
                #sum2 <- rbind(0, term3)[cbind(c(1 + pmin(yindex, tindex)), c(col(yindex)))]
                #if (ny==3) sum3 <- 
                #             rbind(0, term3)[c(cbind(1 + pmin(startindex, tindex)), 
                #                               c(col(yindex)))]
                if (ny==2) D <- matrix(sum1 -  sum2, ncol=ntime)
                else       D <- matrix(sum1 + sum3 - sum2, ncol=ntime)
            }
        } else {
            stop("residuals function still imcomplete, for FH estimate")
            if (any(casewt != casewt[1])) {
                # Have to reconstruct the number of obs with an event, the curve only
                # contains the weighted sum
                nevent <- unlist(lapply(seq(along.with=levels(X)), function(i) {
                    keep <- which(as.numeric(X) ==i)
                    counts <- table(Y[keep, ny-1], status)
                    as.vector(counts[, ncol(counts)])
                    }))
            } else nevent <- fit$n.event

            n2 <- fit$n.risk
            risk2 <- 1/fit$n.risk
            ltemp <- risk2^2
            for (i in which(nevent>1)) {  # assume not too many ties
                denom <- fit$n.risk[i] - fit$n.event[i]*(0:(nevent[i]-1))/nevent[i] 
                risk2[i] <- mean(1/denom) # multiplier for the event
                ltemp[i] <- mean(1/denom^2)
                n2[i] <- mean(denom)
            }

            death <- (yindex <= tindex & rep(event, ntime))
            term1 <- risk2[ff]
            term2 <- lapply(fitrow, function(i) event[i]*ltemp[i])
            term3 <- unlist(lapply(term2, cumsum))

            sum1 <- c(0, term1)[ifelse(death, 1+yindex, 1)]
            sum2 <- c(0, term3)[1 + pmin(yindex, tindex)]
            if (ny==3) sum3 <- c(0, term3)[1 + pmin(startindex, tindex)]

            if (ny==2) D <- matrix(sum1 -  sum2, ncol=ntime)
            else       D <- matrix(sum1 + sum3 - sum2, ncol=ntime)

            if (type=="pstate") D <- -D* c(0,fit$surv[ff])[1+ tindex]
            else if (type=="auc") {
                auc1 <- lapply(fitrow, function(i) {
                             if (length(i) <=1) 0
                             else c(0, cumsum(diff(fit$time[i]) * (fit$surv[i])[-length(i)]))
                                 })  # AUC at each event time
                auc2 <- lapply(fitrow, function(i) {
                             if (length(i) <=1) 0
                             else {
                                 xx <- sort(unique(c(fit$time[i], times))) # all the times
                                 yy <- (fit$surv[i])[findInterval(xx, fit$time[i])]
                                 auc <- cumsum(c(diff(xx),0) * yy)
                                 c(0, auc)[match(times, xx)]
                                 }})  # AUC at the output times

                # Most often this function is called with a single curve, so make that case
                #  faster.  (Or I presume so: mapply and do.call may be more efficient than 
                #  I think for lists of length 1).
                if (length(fitrow)==1) { # simple case, most common to ask for auc 
                    wtmat <- pmin(outer(auc1[[1]], -auc2[[1]], '+'),0)
                    term1 <- term1 * wtmat
                    term2 <- unlist(term2) * wtmat
                    term3 <- apply(term2, 2, cumsum)
                }
                else { #more than one curve, compute weighted cumsum per curve
                    wtmat <- mapply(function(x, y) pmin(outer(x, -y, "+"), 0), auc1, auc2)
                    term1 <- term1 * do.call(rbind, wtmat)
                    temp <- mapply(function(x, y) apply(x*y, 2, cumsum), term2, wtmat)
                    term3 <- do.call(rbind, temp)
                }

                sum1 <- sum2 <- matrix(0, nrow(yindex), ntime)
                if (ny ==3) sum3 <- sum1
                for (i in 1:ntime) {
                    sum1[,i] <- c(0, term1[,i])[ifelse(death[,i], 1 + yindex[,i], 1)]
                    sum2[,i] <- c(0, term3[,i])[1 + pmin(yindex[,i], tindex[,i])]
                    if (ny==3) sum3[,i] <- c(0, term3[,i])[1 + pmin(startindex[,i], tindex[,i])]
                }
                # Perhaps a bit faster(?), but harder to read. And for AUC people usually only
                #  ask for one time point
                #sum1 <- rbind(0, term1)[cbind(c(ifelse(death, 1+yindex, 1)), c(col(yindex)))]
                #sum2 <- rbind(0, term3)[cbind(c(1 + pmin(yindex, tindex)), c(col(yindex)))]
                #if (ny==3) sum3 <- 
                #             rbind(0, term3)[c(cbind(1 + pmin(startindex, tindex)), 
                #                               c(col(yindex)))]
                if (ny==2) D <- matrix(sum1 -  sum2, ncol=ntime)
                else       D <- matrix(sum1 + sum3 - sum2, ncol=ntime)
            }
        }
    } else { # not hazard based
        death <- (yindex <= tindex & rep(event, ntime))
        # dtemp avoids 1/0.  (When this occurs the influence is 0, since
        #  the curve has dropped to zero; and this avoids Inf in term1 and term2).
        dtemp <- ifelse(fit$n.risk==fit$n.event, 0, 1/(fit$n.risk- fit$n.event))
        term1 <- dtemp[ff]
        term2 <- lapply(fitrow, function(i) dtemp[i]*fit$n.event[i]/fit$n.risk[i])
        term3 <- unlist(lapply(term2, cumsum))

        add1 <- c(0, term1)[ifelse(death, 1+yindex, 1)]
        add2 <- c(0, term3)[1 + pmin(yindex, tindex)]
        if (ny==3) add3 <- c(0, term3)[1 + pmin(startindex, tindex)]

        if (ny==2) D <- matrix(add1 -  add2, ncol=ntime)
        else       D <- matrix(add1 + add3 - add2, ncol=ntime)

        # survival is exp(-H) so the derivative is a simple transform of D
        if (type== "pstate") D <- -D* c(1,fit$surv[ff])[1+ tindex]
        else if (type == "auc") {
            auc1 <- lapply(fitrow, function(i) {
                         if (length(i) <=1) 0
                         else c(0, cumsum(diff(fit$time[i]) * (fit$surv[i])[-length(i)]))
                             })  # AUC at each event time
            auc2 <- lapply(fitrow, function(i) {
                         if (length(i) <=1) 0
                         else {
                             xx <- sort(unique(c(fit$time[i], times))) # all the times
                             yy <- (fit$surv[i])[findInterval(xx, fit$time[i])]
                             auc <- cumsum(c(diff(xx),0) * yy)
                             c(0, auc)[match(times, xx)]
                             }})  # AUC at the output times

            # Most often this function is called with a single curve, so make that case
            #  faster.  (Or I presume so: mapply and do.call may be more efficient than 
            #  I think for lists of length 1).
            if (length(fitrow)==1) { # simple case, most common to ask for auc 
                wtmat <- pmin(outer(auc1[[1]], -auc2[[1]], '+'),0)
                term1 <- term1 * wtmat
                term2 <- unlist(term2) * wtmat
                term3 <- apply(term2, 2, cumsum)
            }
            else { #more than one curve, compute weighted cumsum per curve
                wtmat <- mapply(function(x, y) pmin(outer(x, -y, "+"), 0), auc1, auc2)
                term1 <- term1 * do.call(rbind, wtmat)
                temp <- mapply(function(x, y) apply(x*y, 2, cumsum), term2, wtmat)
                term3 <- do.call(rbind, temp)
            }

            sum1 <- sum2 <- matrix(0, nrow(yindex), ntime)
            if (ny ==3) sum3 <- sum1
            for (i in 1:ntime) {
                sum1[,i] <- c(0, term1[,i])[ifelse(death[,i], 1 + yindex[,i], 1)]
                sum2[,i] <- c(0, term3[,i])[1 + pmin(yindex[,i], tindex[,i])]
                if (ny==3) sum3[,i] <- c(0, term3[,i])[1 + pmin(startindex[,i], tindex[,i])]
            }
            # Perhaps a bit faster(?), but harder to read. And for AUC people usually only
            #  ask for one time point
            #sum1 <- rbind(0, term1)[cbind(c(ifelse(death, 1+yindex, 1)), c(col(yindex)))]
            #sum2 <- rbind(0, term3)[cbind(c(1 + pmin(yindex, tindex)), c(col(yindex)))]
            #if (ny==3) sum3 <- 
            #             rbind(0, term3)[c(cbind(1 + pmin(startindex, tindex)), 
            #                               c(col(yindex)))]
            if (ny==2) D <- matrix(sum1 -  sum2, ncol=ntime)
            else       D <- matrix(sum1 + sum3 - sum2, ncol=ntime)
        }
    }
    D
}
rsurvpart2 <- function(Y, X, casewt, istate, times, cluster, type, fit,
                       method, collapse) {
    ny <- ncol(Y)
    ntime <- length(times)
    nstate <- length(fit$states)
    
    # ensure that Y, istate, and fit all use the same set of states
    states <- fit$states
    if (!identical(attr(Y, "states"), fit$states)) {
        map <- match(attr(Y, "states"), fit$states)
        Y[,ny] <- c(0, map)[1+ Y[,ny]]    # 0 = censored
        attr(Y, "states") <- fit$states
    }
    if (is.null(istate)) istate <- rep(1L, nrow(Y)) #everyone starts in s0
    else {
        if (is.character(istate)) istate <- factor(istate)
        if (is.factor(istate)) {
            if (!identical(levels(istate), fit$states)) {
                map <- match(levels(istate), fit$states)
                if (any(is.na(map))) stop ("invalid levels in istate")
                istate <- map[istate]
            }       
        } # istate is numeric, we take what we get and hope it is right
    }

    # collapse redundant rows in Y, for efficiency
    #  a redundant row is a censored obs in the middle of a chain of times
    #  If the user wants individial obs, however, we would just have to
    #  expand it again
    if (ny==3 && collapse & any(duplicated(cluster))) {
        ord <- order(cluster, Y[,1])  # time within subject
        cfit <- .Call(Ccollapse, Y, X, istate, cluster, casewt, ord -1L) 
        if (nrow(cfit) < .8*length(X))  {
            # shrinking the data by 20 percent is worth it
            temp <- Y[ord,]
            Y <- cbind(temp[cfit[,1], 1], temp[cfit[2], 2:3])
            X <- X[cfit[,1]]
            istate <- istate[cfit[1,]]
            cluster <- cluster[cfit[1,]]
        }       
    }

    # Compute the initial leverage
    inf0 <- NULL
    if (is.null(fit$call$p0) && any(istate != istate[1])) { 
        #p0 was not supplied by the user, and the intitial states vary
        inf0 <- matrix(0., nrow=nrow(Y), ncol=nstate)
        i0fun <- function(i, fit, inf0) {
            # reprise algorithm in survfitCI
            p0 <- fit$p0
            t0 <- fit$time[1]
            if (ny==2) at.zero <- which(as.numeric(X) ==i)
            else       
                at.zero <- which(as.numeric(X) ==i &
                          (Y[,1] < t0 & Y[,2] >= t0))
            for (j in 1:nstate) {
                inf0[at.zero, j] <- (ifelse(istate[at.zero]==states[j], 1, 0) -
                                     p0[j])/sum(casewt[at.zero])
            }
            inf0
        }

        if (is.null(fit$strata)) inf0 <- i0fun(1, fit, inf0)
        else for (i in 1:length(levels(X)))
            inf0 <- i0fun(i, fit[i], inf0)  # each iteration fills in some rows
    }

    p0 <- fit$p0          # needed for method==1, type != cumhaz
    fit <- survfit0(fit)  # package the initial state into the picture
    start.time <- fit$time[1]

    # This next block is identical to the one in rsurvpart1, more comments are
    #  there
    etime <- (rowSums(fit$n.event) >0)
    event <- (Y[,ny] >0)
    # 
    #  Create a list whose first element contains the location of
    #   the death times in curve 1, second element for curve 2, etc.
    #  
    if (is.null(fit$strata)) fitrow <- list(which(etime))
    else {
        temp1 <- cumsum(fit$strata)
        temp2 <- c(1, temp1+1)
        fitrow <- lapply(1:length(fit$strata), function(i) {
            indx <- seq(temp2[i], temp1[i])
            indx[etime[indx]] # keep the death times
        }) 
    }
    ff <- unlist(fitrow)

    # for each time x, the index of the last death time which is <=x.
    #  0 if x is before the first death time
    matchfun <- function(x, fit, index) {
        dtime <- fit$time[index]  # subset to this curve
        i2 <- findInterval(x, dtime, left.open=FALSE)
        c(0, index)[i2 +1]
    }
     

    if (type== "cumhaz") {
        # output matrix D will have one row per observation, one col for each
        #  reporting time. tindex and yindex have the same dimension as D.
        # tindex points to the last death time in fit which
        #  is <= the reporting time.  (If there is only 1 curve, each col of
        #  tindex will be a repeat of the same value.)
        tindex <- matrix(0L, nrow(Y), length(times))
        for (i in 1:length(fitrow)) {
            yrow <- which(as.integer(X) ==i)
            temp <- matchfun(times, fit, fitrow[[i]])
            tindex[yrow, ] <- rep(temp, each= length(yrow))
        }
        tindex[,] <- match(tindex, c(0,ff)) -1L  # the [,] preserves dimensions

        # repeat the indexing for Y onto fit$time.  Each row of yindex points
        #  to the last row of fit with death time <= Y[,ny]
        ny <- ncol(Y)
        yindex <- matrix(0L, nrow(Y), length(times))
        event <- (Y[,ny] >0)
        if (ny==3) startindex <- yindex
        for (i in 1:length(fitrow)) {
            yrow <- (as.integer(X) ==i)  # rows of Y for this curve
            temp <- matchfun(Y[yrow,ny-1], fit, fitrow[[i]])
            yindex[yrow,] <- rep(temp, ncol(yindex))
            if (ny==3) {
                temp <- matchfun(Y[yrow,1], fit, fitrow[[i]])
                startindex[yrow,] <- rep(temp, ncol(yindex))
            }
        }                    
        yindex[,] <- match(yindex, c(0,ff)) -1L
        if (ny==3) {
            startindex[,] <- match(startindex, c(0, ff)) -1L
            # no subtractions for report times before subject's entry
            startindex <- pmin(startindex, tindex) 
        }

        dstate <- Y[,ncol(Y)]
        istate <- as.integer(istate)
        ntrans <- ncol(fit$cumhaz)  # the number of possible transitions
        D <- array(0, dim=c(nrow(Y), ntime, ntrans))

        scount <- table(istate[dstate!=0], dstate[dstate!=0]) # observed transitions
        state1 <- row(scount)[scount>0]
        state2 <- col(scount)[scount>0]
        temp <- paste(rownames(scount)[state1], 
                      colnames(scount)[state2], sep='.')
        if (!identical(temp, colnames(fit$cumhaz))) stop("setup error")

        for (k in length(state1)) {
            e2 <- Y[,ny] == state2[k]
            add1 <- (yindex <= tindex & rep(e2, ntime))
            lsum <- unlist(lapply(fitrow, function(i) 
                     cumsum(fit$n.event[i,k]/fit$n.risk[i,k]^2)))
            
            term1 <- c(0, 1/fit$n.risk[ff,k])[ifelse(add1, 1+yindex, 1)]
            term2 <- c(0, lsum)[1+pmin(yindex, tindex)]
            if (ny==3) term3 <- c(0, lsum)[1 + startindex]

            if (ny==2) D[,,k] <- matrix(term1 -  term2, ncol=ntime)
            else       D[,,k] <- matrix(term1 + term3 - term2, ncol=ntime)
        }
    } else {
        if (method==1) {
            # Compute the result using the direct method, in C code
            # the routine is called separately for each curve, data in sorted order
            #
            is1 <- as.integer(istate) -1L  # 0 based subscripts for C
            if (is.null(inf0)) inf0 <- matrix(0, nrow=nrow(Y), ncol=nstate)
            if (all(as.integer(X) ==1)) { # only one curve
                if (ny==2) asort1 <- 0L else asort1 <- order(Y[,1], Y[,2]) -1L
                asort2 <- order(Y[,ny-1]) -1L
                tfit <- .Call(Csurvfitresid, Y, asort1, asort2, is1, 
                              casewt, p0, inf0, times, start.time, 
                              type== "auc")

                if (ntime==1) {
                    if (type=="auc") D <- tfit[[2]] else D <- tfit[[1]]
                }
                else {
                    if (type=="auc") D <- array(tfit[[2]], dim=c(nrow(Y), nstate, ntime))
                    else         D <- array(tfit[[1]], dim=c(nrow(Y), nstate, ntime))
                }
            }
            else { # one curve at a time
                ix <- as.numeric(X)  # 1, 2, etc
                if (ntime==1) D <- matrix(0, nrow(Y), nstate)
                else D <- array(0, dim=c(nrow(Y), nstate, ntime))
                for (curve in 1:max(ix)) {
                    j <- which(ix==curve)
                    ytemp <- Y[j,,drop=FALSE]
                    if (ny==2) asort1 <- 0L 
                    else asort1 <- order(ytemp[,1], ytemp[,2]) -1L
                    asort2 <- order(ytemp[,ny-1]) -1L

                    # call with a subset of the data
                    j <- which(ix== curve)
                    tfit <- .Call(Csurvfitresid, ytemp, asort1, asort2, is1[j],
                                  casewt[j], p0[curve,], inf0[j,], times, 
                                  start.time, type=="auc")
                    if (ntime==1) {
                        if (type=="auc") D[j,] <- tfit[[2]] else D[j,] <- tfit[[1]]
                    } else {
                        if (type=="auc") D[j,,] <- tfit[[2]] else D[j,,] <- tfit[[1]]
                    }
                }
            } 
            # the C code makes time the last dimension, we want it to be second
            if (ntime > 1) D <- aperm(D, c(1,3,2))
        }
        else {
            # method 2
            Yold <- Y
            utime  <- fit$time[fit$time <= max(times) & etime] # unique death times
            ndeath <- length(utime)    # number of unique event times
            delta <- diff(c(start.time, utime))

            # Expand Y
            if (ny==2) split <- .Call(Csurvsplit, rep(0., nrow(Y)), Y[,1], times)
            else split <- .Call(Csurvsplit, Y[,1], Y[,2], times)
            X <- X[split$row]
            casewt <- casewt[split$row]
            istate <- istate[split$row]
            Y <- cbind(split$start, split$end, 
                        ifelse(split$censor, 0, Y[split$row,ny]))
            ny <- 3

            # Create a vector containing the index of each end time into the fit object
            yindex <- ystart <- double(nrow(Y))
            for (i in 1:length(fitrow)) {
                yrow <- (as.integer(X) ==i)  # rows of Y for this curve
                yindex[yrow] <- matchfun(Y[yrow, 2], fit, fitrow[[i]])
                ystart[yrow] <- matchfun(Y[yrow, 1], fit, fitrow[[i]])
            }
            # And one indexing the reporting times into fit
            tindex <- matrix(0L, nrow=length(fitrow), ncol=ntime)
            for (i in 1:length(fitrow)) {
                tindex[i,] <- matchfun(times, fit, fitrow[[i]])
            }
            yindex[,] <- match(yindex, c(0,ff)) -1L
            tindex[,] <- match(tindex, c(0,ff)) -1L
            ystart[,] <- pmin(match(ystart, c(0,ff)) -1L, tindex)

            # Create the array of C matrices
            cmat <- array(0, dim=c(nstate, nstate, ndeath)) # max(i2) = ndeath, by design
            Hmat <- cmat

            # We only care about observations that had a transition; any transitions
            #  after the last reporting time are not relevant
            transition <- (Y[,ny] !=0 & Y[,ny] != istate &
                           Y[,ny-1] <= max(times)) # obs that had a transition
            i2 <- match(yindex, sort(unique(yindex)))  # which C matrix this obs goes to
            i2 <- i2[transition]
            from <- as.numeric(istate[transition])  # from this state
            to   <- Y[transition, ny]   # to this state
            nrisk <- fit$n.risk[cbind(yindex[transition], from)]  # number at risk
            wt <- casewt[transition]
            for (i in seq(along.with =from)) {
                j <- c(from[i], to[i])
                haz <- wt[i]/nrisk[i]
                cmat[from[i], j, i2[i]] <- cmat[from[i], j, i2[i]] + c(-haz, haz)
            }
            for (i in 1:ndeath) Hmat[,,i] <- cmat[,,i] + diag(nstate)

            # The transformation matrix H(t) at time t  is cmat[,,t] + I
            # Create the set of W and V matrices.
            # 
            dindex <- which(etime & fit$time <= max(times))
            Wmat <- Vmat <- array(0, dim=c(nstate, nstate, ndeath))
            for (i in ndeath:1) {
                j <- match(dindex[i], tindex, nomatch=0) 
                if (j > 0) {
                    # this death matches one of the reporting times
                    Wmat[,,i] <- diag(nstate)
                    Vmat[,,i] <- matrix(0, nstate, nstate)
                } 
                else {
                    Wmat[,,i] <- Hmat[,,i+1] %*% Wmat[,,i+1]
                    Vmat[,,i] <- delta[i] +  Hmat[,,i+1] %*% Wmat[,,i+1]
                }
            }
            iterm <- array(0, dim=c(nstate, nstate, ndeath)) # term in equation
            itemp <- vtemp <- matrix(0, nstate, nstate)  # cumulative sum, temporary
            isum  <- isum2 <- iterm  # cumulative sum
            vsum  <- vsum2 <- vterm <- iterm
            for (i in 1:ndeath) {
                j <- dindex[i]
                n0 <- ifelse(fit$n.risk[j,] ==0, 1, fit$n.risk[j,]) # avoid 0/0
                iterm[,,i] <- ((fit$pstate[j-1,]/n0) * cmat[,,i]) %*% Wmat[,,i]
                vterm[,,i] <- ((fit$pstate[j-1,]/n0) * cmat[,,i]) %*% Vmat[,,i]
                itemp <- itemp + iterm[,,i]
                vtemp <- vtemp + vterm[,,i]
                isum[,,i] <- itemp
                vsum[,,i] <- vtemp
                j <- match(dindex[i], tindex, nomatch=0)
                if (j>0) itemp <- vtemp <- matrix(0, nstate, nstate)  # reset
                isum2[,,i] <- itemp
                vsum2[,,i] <- vtemp
            }

            # We want to add isum[state,, entry time] - isum[state,, exit time] for
            #  each subject, and for those with an a:b transition there will be an 
            #  additional vector with -1, 1 in the a and b position.
            i1 <- match(ystart, sort(unique(yindex)), nomatch=0) # start at 0 gives 0
            i2 <- match(yindex, sort(unique(yindex)))
            D <- matrix(0., nrow(Y), nstate)
            keep <- (Y[,2] <= max(times))  # any intervals after the last reporting time
                                            # will have 0 influence
            for (i in which(keep)) {
                if (Y[i,3] !=0 && istate[i] != Y[i,3]) {
                    z <- fit$pstate[yindex[i]-1, istate[i]]/fit$n.risk[yindex[i], istate[i]]
                    temp <- double(nstate)
                    temp[istate[i]] = -z
                    temp[Y[i,3]]    =  z
                    temp <- temp %*% Wmat[,,i2[i]] - isum[istate[i],,i2[i]]
                    if (i1[i] >0) temp <- temp + isum2[istate[i],, i1[i]]
                    D[i,] <- temp
                }
                else {
                    if (i1[i] >0) D[i,] = isum2[istate[i],,i1[i]] - isum[istate[i],, i2[i]]
                    else  D[i,] =  -isum[istate[i],, i2[i]]
                }
            }
            Dsave <- D
            if (!is.null(inf0)) {
                # add in the initial influence, to the first row of each obs
                #   (inf0 was created on unsplit data)
                j <- which(!duplicated(split$row))
                D[j,] <- D[j,] + (inf0%*% Hmat[,,1] %*% Wmat[,,1])
            }
            if (ntime > 1) {
                interval <- findInterval(yindex, tindex, left.open=TRUE)
                D2 <- array(0., dim=c(dim(D), ntime))
                D2[interval==0,,1] <- D[interval==0,]
                for (i in 1:(ntime-1)) {
                    D2[interval==i,,i+1] = D[interval==i,]
                    j <- tindex[i]
                    D2[,,i+1] = D2[,,i+1] + D2[,,i] %*% (Hmat[,,j] %*% Wmat[,,j])
                } 
                D <- D2
            }

            # undo any artificial split
            if (any(duplicated(split$row))) {
                if (ntime==1) D <- rowsum(D, split$row)
                else {
                    # rowsums has to be fooled
                    temp <- rowsum(matrix(D, ncol=(nstate*ntime)), split$row)
                    # then undo it
                    D <- array(temp, dim=c(nrow(temp), nstate, ntime))
                }
            }
        }
    }   

    # since we may have done a partial collapse (removing redundant rows), the
    # parent routine can't collapse the data
    if (collapse & any(duplicated(cluster))) {
        if (length(dim(D)) ==2)
            D <- rowsum(D, cluster, reorder=FALSE)
        else { #rowsums has to be fooled
            dd <- dim(D)
            temp <- rowsum(matrix(D, nrow=dd[1]), cluster)
            D <- array(temp, dim=c(nrow(temp), dd[2:3]))
        }       
    }
    D
}
